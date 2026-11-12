#!/usr/bin/env python3
"""Plot measured, accepted solver bounds for a fixed cohort of paired TAL runs.

Uses the original cross-TAL normalization and TikZ renderer. The x axis is
shared elapsed time, including SCOPF and main-model setup, rather than a sum of
solver times from concurrent processes. No cross-run minimum or backfill is used.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd


def step(times, values, grid, initial=np.nan, minimum=True):
    good = np.isfinite(values) & (np.abs(values) < 1e90)
    times, values = times[good], values[good]
    if not len(times):
        return np.full_like(grid, initial)
    order = np.argsort(times, kind="stable")
    times, values = times[order], values[order]
    values = np.minimum.accumulate(values) if minimum else np.maximum.accumulate(values)
    indices = np.searchsorted(times, grid, side="right") - 1
    return np.where(indices >= 0, values[np.maximum(indices, 0)], initial)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--horizon", type=float, default=120)
    parser.add_argument("--style-script", type=Path, required=True)
    args = parser.parse_args()
    if not np.isfinite(args.horizon) or args.horizon <= 1:
        parser.error("The plot horizon must be finite and greater than one second")
    args.output.mkdir(parents=True, exist_ok=True)
    grid = np.geomspace(1, args.horizon, 701)
    curves = {1: [[], []], 2: [[], []]}
    rows, provenance, events = [], [], []
    seen_instances = set()
    for run in args.runs:
        manifest = json.loads((run / "manifest.json").read_text())
        if manifest.get("helper", "tal1") != "tal1" or manifest.get("main_tal", 2) != 2:
            raise ValueError("This renderer requires TAL=1 -> TAL=2; SCOPF/RSCOTS helpers need their own labels and reference")
        if not manifest.get("complete"):
            raise ValueError(f"Incomplete pair: {run}")
        identity = manifest.get("case_sha256", manifest["case"])
        if identity in seen_instances:
            raise ValueError(f"Duplicate instance in fixed cohort: {run}")
        seen_instances.add(identity)
        donor = pd.read_csv(run / "donor.events.csv")
        receiver_role = manifest["roles"][1]
        receiver = pd.read_csv(run / f"{receiver_role}.events.csv")
        for role, df in (("donor", donor), (receiver_role, receiver)):
            result = json.loads((run / role / "result.json").read_text())
            info = result["info"]
            failures = {key: result.get("callback", {}).get(key, 0)
                        for key in ("n_benders_cut_value_fail",)}
            if any(failures.values()):
                raise ValueError(f"Unresolved solver validation failures: {run}, {role}: {failures}")
            if info["status"] not in ("OPTIMAL", "TIME_LIMIT", "MEMORY_LIMIT"):
                raise ValueError(f"Invalid solver result: {run}, {role}: {info['status']}")
            reported = float(df[df.event == "final"].iloc[-1].primal)
            if info["objective"] is None or not np.isclose(reported, info["objective"], rtol=1e-6, atol=1e-6):
                raise ValueError(f"Final bound disagrees with result: {run}, {role}")
        # Final observed donor incumbent is the common reference, even if the
        # short run did not prove optimality. Record that status explicitly.
        final = donor[donor.event == "final"].iloc[-1]
        reference = float(final.primal)
        if not np.isfinite(reference):
            raise ValueError(f"No finite donor incumbent: {run}")
        provenance.append({"run": str(run.resolve()), "reference": reference,
                           "donor_status": final.detail, "manifest": manifest})
        for z, df in ((1, donor), (2, receiver)):
            accepted = df[df.event.isin(["scopf", "bounds", "final"])]
            p = step(accepted.time_s.to_numpy(), accepted.primal.to_numpy(), grid)
            main = accepted[accepted.event != "scopf"]
            d = step(main.time_s.to_numpy(), main.dual.to_numpy(), grid, minimum=False)
            if np.any(d > p + 1e-5 * max(abs(reference), 1.0)):
                raise ValueError(f"Dual exceeds primal: {run}, TAL={z}")
            if z == 2 and np.any(d > reference + 1e-5 * max(abs(reference), 1.0)):
                raise ValueError(f"TAL=2 dual contradicts feasible TAL=1 reference: {run}")
            curves[z][0].append(1 + (p - reference) / max(abs(reference), 1.0))
            curves[z][1].append(1 + (d - reference) / max(abs(reference), 1.0))
            rows.append({"instance": manifest["case"], "tal": z, "reference": reference,
                         "primal": p[-1], "dual": d[-1], "relative_primal": curves[z][0][-1][-1],
                         "final_status": df[df.event == "final"].iloc[-1].detail})
        for _, e in receiver[receiver.event.isin(["submitted", "confirmed", "not_improving"])].iterrows():
            events.append({"run": str(run), **e.to_dict()})
    # Every point represents the same fixed cohort. Wait until every pair has
    # an observed feasible incumbent; do not silently average fewer instances.
    aggregate = {z: (np.mean(c[0], axis=0), np.mean(c[1], axis=0)) for z, c in curves.items()}
    visible = np.logical_and.reduce([np.isfinite(c[0]) for c in aggregate.values()])
    if not visible.any():
        raise ValueError("No common feasible interval within horizon")
    grid = grid[visible]
    aggregate = {z: (p[visible], d[visible]) for z, (p, d) in aggregate.items()}
    pd.DataFrame({"time_s": grid, **{f"z{z}_{kind}": bounds[i] for z, bounds in aggregate.items()
                                   for i, kind in enumerate(("primal", "dual"))}}).to_csv(args.output / "mean_bounds.csv", index=False)
    pd.DataFrame(rows).to_csv(args.output / "per_instance_endpoints.csv", index=False)
    pd.DataFrame(events).to_csv(args.output / "injections.csv", index=False)
    (args.output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    (args.output / "plot_metadata.json").write_text(json.dumps({
        "horizon": args.horizon,
        "plot_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "style_script": str(args.style_script.resolve()),
        "style_script_sha256": hashlib.sha256(args.style_script.read_bytes()).hexdigest(),
        "missing_bounds": "Not drawn until observed; fixed cohort throughout.",
    }, indent=2) + "\n")
    sys.path.insert(0, str(args.style_script.parent))
    spec = importlib.util.spec_from_file_location("original_cross_tal", args.style_script)
    style = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(style)
    style.write_plot(grid, aggregate, len(args.runs), args.output)
    tikz = args.output / "scopf_cross_tal.tikz.tex"
    tex = tikz.read_text().replace("Cumulative solver time (s)", "Elapsed solve time (s)")
    tex = tex.replace("cohort_provenance.json", "provenance.json")
    tex = tex.replace("no marks,color=", "const plot,no marks,color=")
    tex = tex.replace(r"width=0.72\linewidth, height=0.42\linewidth,",
                      r"width=0.92\linewidth, height=0.48\linewidth,ytick distance=0.2,")
    tex = tex.replace("legend columns=1", "legend columns=2")
    tex = tex.replace("at={(1.04,0.5)},anchor=west", "at={(0.5,1.025)},anchor=south")
    # Missing dual bounds form only a leading interval. Omit those coordinates
    # so the original renderer draws no invented leading bound.
    for z, (_, d) in aggregate.items():
        observed = grid[np.isfinite(d)]
        if len(observed):
            tex = tex.replace(f"fill between[of=crossP{z} and crossD{z}]",
                              f"fill between[of=crossP{z} and crossD{z},soft clip={{domain={observed[0]:.10g}:{grid[-1]:.10g}}}]")
        else:
            tex = "\n".join(line for line in tex.splitlines() if f"of=crossP{z} and crossD{z}" not in line)
    tikz.write_text("\n".join(line for line in tex.splitlines() if ",nan)" not in line) + "\n")
    # A directly viewable artifact as well as the original TikZ/PDF.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 11})
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for z, (p, d) in aggregate.items():
        color = "#" + style.COLORS[z]
        ax.step(grid, p, where="post", color=color, linewidth=1.4, label=rf"$Z={z}$")
        ax.step(grid, d, where="post", color=color, linewidth=1, linestyle=":")
        ax.fill_between(grid, d, p, step="post", color=color, alpha=.045)
    ax.axhline(1, color="gray", linewidth=.8, linestyle=(0, (4, 4)))
    ax.set(xscale="log", xlim=(1, args.horizon), ylim=(0, max(1.05, np.ceil(max(c[0].max() for c in aggregate.values())*10)/10)),
           xlabel="Elapsed solve time (s)", ylabel=r"Objective bound relative to $Z=1$")
    ticks = [t for t in (1, 10, 100, 1000, 10000) if t <= args.horizon / 2] + [args.horizon]
    ax.set_xticks(ticks, [f"{t:g}" for t in ticks])
    ax.grid(True, which="major", alpha=.3)
    ax.legend(loc="lower center", bbox_to_anchor=(.5, 1.01), ncol=2, frameon=False)
    ax.text(.03, .045, rf"$n={len(args.runs)}$; solid: primal, dotted: dual", transform=ax.transAxes, color="gray")
    fig.tight_layout()
    fig.savefig(args.output / "parallel_tal.png", dpi=180)
    fig.savefig(args.output / "parallel_tal.pdf")
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
