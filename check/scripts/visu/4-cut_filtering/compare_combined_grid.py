#!/usr/bin/env python3
"""Stage 4 — re-run the OLD combined CF×SP grid search on **random** SP selection
through the *current* (fixed) parsers, to confirm the best `cfNN_spNN` combo
still holds after the 1-baseline cut-error / time-limit fixes.

Context (EXPERIMENTS.md §12.3b): the combined cut-filtering × SP-selection grid
was already searched on random, and is **cold** (no `[WARMSTARTING]`) on the
31-inst grid. If random beats diversity (§11.0) we reuse it as-is — but the
selection of the best combo was made with the *old* reader, before the cut-error
majority-vote and time-limit-capping fixes landed here. This re-derives the best
combo under the fixed reader, with proper per-experiment grid cut-err scoping
(the old auto scripts applied none).

Zero-touch: reads the read-only ``4-cut_filtering`` results, writes only into
``1-baseline/check/results/paper/4-cut_filtering/<exp>/``.

Random combined grids (un-suffixed = ``SORTING_MODE = :RND``; the ``-spdiv`` /
``-cfspdiv`` siblings are the warm diversity versions, excluded here):

    3-cutfiltering-abs-spsel-abs   abs × abs
    4-cutfiltering-rel-spsel-rel   rel × rel
    5-cutfiltering-rel-spsel-abs   rel × abs
    6-cutfiltering-abs-spsel-rel   abs × rel
    6-cutfiltering-both-both       both × both (full merge)

Baseline for each = the loosest filter (largest limits ≈ "keep everything"), so
the table reads as "how much does filtering to cfX/spY help vs no filtering".

Run:  python3 compare_combined_grid.py
"""

import re
import sys
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]          # visu/4-cut_filtering -> 1-baseline
SHARED = PROJECT_ROOT.parent
sys.path.insert(0, str(SCRIPT_DIR.parent))    # visu/ for reader, tables
sys.path.insert(0, str(SCRIPT_DIR))           # for _cutfiltering_common

from reader import CUT_ERR_OBJ_TOLERANCE        # noqa: E402
from selection import TIME_LIMIT                 # noqa: E402
from tables import run_full_comparison           # noqa: E402
from _cutfiltering_common import (               # noqa: E402
    build_label_map, build_style_map, discover_settings, parse_setting,
)

STAGE = "4-cut_filtering"
SRC_ROOT = SHARED / STAGE / "check/results"
PAPER_ROOT = PROJECT_ROOT / "check/results/paper" / STAGE / "grid"

# exp -> baseline (loosest filter = largest cf & sp limit, the "no-filter" ref)
GRIDS = {
    "3-cutfiltering-abs-spsel-abs": "cf10_sp10",
    "4-cutfiltering-rel-spsel-rel": "cf100_sp100",
    "5-cutfiltering-rel-spsel-abs": "cf100_sp10",
    "6-cutfiltering-abs-spsel-rel": "cf10_sp100",
    "6-cutfiltering-both-both":     "cf100_sp100",
}
MIN_BUSES = 73        # match the original stage-4 grid analysis (dispatch default)


def _penalize_cut_errors(target: str):
    """df_transform: charge every cut-erroring *cell* a TIME_LIMIT penalty
    instead of evicting the whole instance.

    Two kinds of cut error, both penalized on the offending (instance, setting)
    cell only (§12.3d):
      1. explicit ``status == CUT_ERROR`` — already clipped to TIME_LIMIT by the
         inclusive `_all` view, penalized here too for the strict view;
      2. objective deviation — a run that returned ``OPTIMAL`` but with an
         objective off the per-instance consensus (median of the OPTIMAL runs) by
         > CUT_ERR_OBJ_TOLERANCE. These are the dangerous ones: a setting that
         cuts off the true optimum finishes fast and looks *best*, so it must be
         penalized, not rewarded.

    Keeps every instance in the comparison — the whole point, since dropping all
    cut-err instances collapses the common set of an 80-wide grid to ~12 (§12.3c).
    Consensus is the median (baseline is not privileged), matching
    reader.collect_cut_err_instances' philosophy.
    """
    def _tf(df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()
        n_expl = n_dev = 0
        if "status" in df.columns:
            # Explicit cut errors need no touch: CUT_ERROR is not a BROKEN status,
            # so it stays in `ran` and the inclusive `_all` view already clips it to
            # TIME_LIMIT in the values while the CutErr column still counts it.
            n_expl = int((df["status"].astype(str) == "CUT_ERROR").sum())
        if {"status", "instance", "objective", "settings"}.issubset(df.columns):
            # Objective-deviation cut errors: status is OPTIMAL (fast but wrong).
            # Relabel to CUT_ERROR (NOT TIME_LIMIT) so they leave the solved set,
            # show in CutErr, and get the same TIME_LIMIT value-penalty in `_all`.
            opt = df["status"].astype(str) == "OPTIMAL"
            obj = pd.to_numeric(df["objective"], errors="coerce")
            cons = obj.where(opt).groupby(df["instance"]).transform("median")
            rel = (obj - cons).abs() / cons.abs().clip(lower=1e-10)
            bad = opt & rel.notna() & (rel > CUT_ERR_OBJ_TOLERANCE)
            n_dev = int(bad.sum())
            df.loc[bad, "status"] = "CUT_ERROR"
        print(f"    cut errors -> CutErr (time-penalized in _all): "
              f"{n_expl} explicit + {n_dev} obj-deviation cell(s)")
        return df
    return _tf


def _uniq_labels(settings: list[str]) -> dict[str, str]:
    """Label map guaranteed injective — the both-both merge has settings that
    format to the same label (cf100_sp100 / _abs / _rel), which would create
    duplicate table columns. Disambiguate collisions with the raw name."""
    base = build_label_map(settings)
    used: set[str] = set()
    out: dict[str, str] = {}
    for s, lbl in base.items():
        cand, n = lbl, 0
        while cand in used:                       # LaTeX-safe suffix (letter) inside braces
            suf = chr(ord("b") + n)
            cand = (lbl[:-1] + rf"~{suf}}}") if lbl.endswith("}") else f"{lbl} {suf}"
            n += 1
        used.add(cand)
        out[s] = cand
    return out


MAX_CUT_ERR = 2   # a run with >2 cut errors is too unreliable to be "best"


def _num(cell):
    m = re.search(r"[-+]?\d+\.?\d*", str(cell).replace("**", ""))
    return float(m.group()) if m else None


def _best_from_short(md_path: Path, label_map: dict) -> None:
    """Rank settings from the short `_all` table (rows=settings), gating out any
    run with more than MAX_CUT_ERR cut errors before picking min Sum / Geo."""
    if not md_path.is_file():
        print(f"    (no {md_path.name} to rank)")
        return
    lines = [l for l in md_path.read_text().splitlines()
             if l.lstrip().startswith("|")]
    if len(lines) < 3:
        return
    cols = [c.strip() for c in lines[0].strip("|").split("|")]
    idx = {c: i for i, c in enumerate(cols)}
    inv = {v: k for k, v in label_map.items()}
    ce_i = idx.get("CutErr")

    rows = []  # (setting, cut_err, sum, geo) — one per setting row
    for l in lines[2:]:
        c = [x.strip() for x in l.strip("|").split("|")]
        if len(c) < len(cols):
            continue
        lbl = c[0].replace("**", "").strip()
        name = inv.get(lbl)          # only real settings map back; skip anything else
        if name is None:
            continue
        ce = int(_num(c[ce_i]) or 0) if ce_i is not None else 0
        rows.append((name, ce, _num(c[idx.get("Sum", -1)]),
                     _num(c[idx.get("Geo. Mean", -1)])))
    elig = [r for r in rows if r[1] <= MAX_CUT_ERR]
    dropped = len(rows) - len(elig)
    print(f"    eligible (≤{MAX_CUT_ERR} cut-err): {len(elig)}/{len(rows)} settings"
          + (f"  [{dropped} excluded for >{MAX_CUT_ERR} cut errors]" if dropped else ""))
    for key, name in ((2, "Sum"), (3, "Geo. Mean")):
        cand = [r for r in elig if r[key] is not None]
        if cand:
            b = min(cand, key=lambda r: r[key])
            print(f"    best by {name:<9}: {b[0]}  ({b[key]:g}, {b[1]} cut-err)")


def main() -> None:
    for exp, baseline in GRIDS.items():
        src = SRC_ROOT / exp
        if not src.is_dir():
            print(f"[{exp}] missing source; skip")
            continue
        settings = [s for s in discover_settings(src)
                    if parse_setting(s).get("cf") or parse_setting(s).get("sp")]
        if baseline not in settings:
            # fall back to the largest-limit setting actually present
            baseline = settings[-1]
        out = PAPER_ROOT / exp
        print(f"[{exp}] {len(settings)} settings, baseline {baseline} "
              f"(cut errors -> CutErr, time-penalized, instances kept)")
        run_full_comparison(
            src, baseline=baseline, output_dir=out, target="time_solve",
            df_transform=_penalize_cut_errors("time_solve"),
            pd_convergence=False, min_buses=MIN_BUSES,
            label_map=_uniq_labels(settings), style_map=build_style_map(settings),
            include_settings=settings,
        )
        # Rank from the inclusive short `_all` table: cut errors are time-penalized
        # (kept, not dropped) and counted in CutErr; we then pick the best run
        # among those with ≤ MAX_CUT_ERR cut errors.
        _best_from_short(out / "comparison_short_all.md", _uniq_labels(settings))
        _heatmap(out, baseline)


def _heatmap(out, baseline) -> None:
    """2-D CF×SP heatmap built from the table's own ``Geo. Mean`` values, so the
    cells match the table exactly (same instance selection + cut-error policy).
    Labels read ``CF=<X> / SP=<Y>[%]``; a dropped setting is a blank cell."""
    from heatmap import grid_heatmap_from_values, geomeans_from_short_table

    def coord(label):
        m = re.search(r"CF=(\d+)\s*/\s*SP=(\d+)", label)
        return (int(m.group(1)), int(m.group(2))) if m else None

    md = out / "comparison_short_all.md"
    if not md.is_file():
        return
    gm, nn, note = geomeans_from_short_table(md, coord)
    if not gm:
        print(f"    (no 2-D heatmap for {out.name})")
        return
    xs = sorted({x for x, _ in gm})
    ys = sorted({y for _, y in gm})
    bm = re.search(r"cf(\d+)_sp(\d+)", baseline or "")
    bxy = (int(bm.group(1)), int(bm.group(2))) if bm else None
    grid_heatmap_from_values(
        gm, nn, out / "heatmap_time_solve_all.pdf", x_levels=xs, y_levels=ys,
        baseline_xy=bxy, xlabel=r"CF threshold (\%)", ylabel=r"SP threshold (\%)",
        n_note=note,
    )
    print(f"    heatmap -> {out / 'heatmap_time_solve_all.pdf'} ({len(gm)} cells from table)")


if __name__ == "__main__":
    main()
