#!/usr/bin/env python3
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                           #
#  Copyright 2026 Operations Research, RWTH Aachen University               #
#                                                                           #
#  Licensed under the Apache License, Version 2.0 (the "License");          #
#  you may not use this file except in compliance with the License.         #
#  You may obtain a copy of the License at                                  #
#                                                                           #
#      http://www.apache.org/licenses/LICENSE-2.0                           #
#                                                                           #
#  Unless required by applicable law or agreed to in writing, software      #
#  distributed under the License is distributed on an "AS IS" BASIS,        #
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. #
#  See the License for the specific language governing permissions and      #
#  limitations under the License.                                           #
#                                                                           #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

""" Build all tables and figures.

Usage: python check/scripts/visu/compare.py <experiment>

Output:
  comparison_short.md, comparison_extra.tex, time_decomp/table.tex,
  performance_profiles/{total,SP,solve_time,solve_time_scaled}/factor.{pdf,tikz},
  performance_profiles/solve_time_scaled/scatter_sorted_by_sle_runtime.{pdf,tikz},
  time_decomp_boxplot/{detailed,detailed_relative}.{pdf,tikz}
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plots import TIME_COMPONENTS, box_height_cm, performance_profile, runtime_scatter, time_boxplot

PROJECT_ROOT = Path(__file__).resolve().parents[3]
EXPERIMENT_DIR = PROJECT_ROOT / "check/results" / sys.argv[1]
INSTANCES_DIR = PROJECT_ROOT / "check/instances/main"
MIN_BUSES = 73  # drop instances with fewer buses
MIN_TIME = 50  # drop instances where all settings finish below this (seconds)
TIME_LIMIT = 10800
BAD_STATUSES = {"CUT_ERROR", "INFEASIBLE", "INFEASIBLE_OR_UNBOUNDED", "DID_NOT_RUN"}

FAMILIES = {
    "FLLP": (r"\texttt{GRB}", {"color": "#000000", "linestyle": "-"}),
    "LEQ": (r"\texttt{SLE}", {"color": "#E69F00", "linestyle": "--"}),
    "LEQRF": (r"\texttt{SLE-U}", {"color": "#0072B2", "linestyle": ":"}),
}


def family(name: str) -> str | None:
    return next((t for t in FAMILIES if f"_{t}_" in f"_{name.upper()}_"), None)


def fmt(v) -> str:
    return "--" if pd.isna(v) else f"{float(v):.4f}"


def geo_mean(x: pd.Series) -> float:
    return np.exp(np.log(x[x > 0]).mean()) if (x > 0).any() else np.nan


def keep_groups(df: pd.DataFrame, keys: list[str], flags: pd.Series) -> pd.DataFrame:
    return df.merge(flags[flags].reset_index()[keys], on=keys)


def add_time_components(df: pd.DataFrame) -> pd.DataFrame:
    """Add the "_sp_*" columns of plots.TIME_COMPONENTS."""
    df = df.copy()
    col = lambda c: df[c].fillna(0.0) if c in df.columns else pd.Series(0.0, index=df.index)
    fallback = col("time_SP_N1_fallback")
    cut = col("time_SP_N1_cut_build") + col("time_SP_N1_cut_check")
    df["_sp_dual"] = col("time_SP_N1_dual") + 0.5 * fallback
    df["_sp_primal"] = col("time_SP_N1_primal") + 0.5 * fallback
    accounted = (col("time_SP_conn") + col("time_SP_N1_dual") + col("time_SP_N1_primal") + fallback + cut
                 + col("time_SP_N1_refactorization"))
    df["_sp_remaining"] = (col("time_SP") - accounted).clip(lower=0.0) + col("time_SP_conn") + cut
    return df


def short_table_md(df, keys, settings, labels) -> str:
    """Solved count, sum, mean, geo. mean."""
    time = df["time_solve"].where(df["status"] != "TIME_LIMIT", TIME_LIMIT)
    vals = (df.assign(_t=time).pivot_table(index=keys, columns="settings", values="_t", aggfunc="first")
            .reindex(columns=settings).apply(pd.to_numeric, errors="coerce"))
    solved = (vals < TIME_LIMIT) & vals.notna()
    common = vals[solved.all(axis=1)]
    table = pd.DataFrame(
        [[fmt(v) for v in (solved[s].sum(), common[s].sum(), common[s].mean(), geo_mean(common[s]))] for s in settings],
        index=[labels[s] for s in settings], columns=["Solved", "Sum", "Mean", "Geo. Mean"],
    )
    return table.to_markdown(disable_numparse=True) + f"\n\n*n = {len(common)} instances (solved by all settings)*"


def _matpower_rows(lines: list[str], table: str) -> int | None:
    it = iter(lines)
    for line in it:
        if line.startswith(f"mpc.{table}") and "=" in line and "[" in line:
            n = 0
            for line in it:
                if line.startswith("];"):
                    return n
                n += bool(line) and not line.startswith("%")
            return None
    return None


def graph_stats(instances: set[str]) -> pd.DataFrame:
    """Instance file stats."""
    rows, missing = {}, []
    for name in sorted(instances):
        path = next((p for p in (INSTANCES_DIR / name, INSTANCES_DIR / f"{name}.m") if p.is_file()), None)
        lines = [l.strip() for l in path.read_text(encoding="utf-8", errors="ignore").splitlines()] if path else []
        v, a = _matpower_rows(lines, "bus"), _matpower_rows(lines, "branch")
        if not v or a is None:
            missing.append(name)
            continue
        v, a = float(v), float(a)
        rows[name] = [v, a, a / v, (v + 2.0 * a) / (v ** 2)]
    if missing:
        print(f"WARNING: no graph stats for {len(missing)} instances in {INSTANCES_DIR}: {', '.join(missing[:8])}")
    return pd.DataFrame.from_dict(rows, orient="index", columns=["nodes", "branches", "ratio", "fill"])


def extra_table_tex(df, settings, labels) -> str:
    """Instance file stats -> tex."""
    stem = lambda name: Path(name).name.removesuffix(".m")
    stats = graph_stats({stem(i) for i in df["instance"].dropna().astype(str)})
    cols = ["Solved Subs Sum", "Solved Subs Mean", "Solved Subs Geo. Mean",
            "Avg Nodes", "Avg Branches", "Avg A/V Ratio", "Avg Laplacian Fill"]
    lines = [r"\begin{tabular}{l" + "r" * len(cols) + "}", r"\toprule", " & ".join([""] + cols) + r" \\", r"\midrule"]
    for s in settings:
        sdf = df[df["settings"] == s]
        subs = pd.to_numeric(sdf["n_solved_subs_total"], errors="coerce").dropna()
        graph = stats.reindex(sorted({stem(i) for i in sdf["instance"].dropna().astype(str)})).dropna(how="all")
        vals = [subs.sum() if len(subs) else np.nan, subs.mean(), geo_mean(subs)] + [graph[c].mean() for c in graph]
        lines.append(" & ".join([labels[s]] + [fmt(v) for v in vals]) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}", f"% graph stats parsed for {len(stats)} unique instances"]
    return "\n".join(lines)


def time_table_tex(df, settings, labels) -> str:
    cols = [c for c, _ in TIME_COMPONENTS]
    mean = df.groupby("settings")[cols].mean().reindex(settings)
    total = mean.sum(axis=1).replace(0, 1.0)
    pct = mean.div(total, axis=0) * 100.0
    lines = [r"\begin{tabular}{l" + "r" * (len(cols) + 1) + "}", r"\toprule",
             " & ".join(["Setting"] + [label for _, label in TIME_COMPONENTS] + ["Total"]) + r" \\", r"\midrule"]
    for s in settings:
        cells = [f"{mean.loc[s, c]:.4f} ({pct.loc[s, c]:.4f}\\%)" for c in cols]
        lines.append(" & ".join([labels[s]] + cells + [f"{total[s]:.4f}"]) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines)


def main() -> None:
    exp = EXPERIMENT_DIR
    settings = sorted(
        (d.name for d in exp.iterdir()
         if d.is_dir() and family(d.name) and not d.name.endswith("_infeasible") and any(d.glob("*.csv"))),
        key=lambda s: (list(FAMILIES).index(family(s)), s),
    )
    labels = {s: FAMILIES[family(s)][0] for s in settings}
    styles = {s: FAMILIES[family(s)][1] for s in settings}
    height_cm = box_height_cm(len(settings))

    # Load & filter.
    df = pd.concat([pd.read_csv(f) for s in sorted(settings) for f in sorted((exp / s).glob("*.csv"))],
                   ignore_index=True)
    buses = df["instance"].str.extract(r"case(\d+)", expand=False).fillna(0).astype(int)
    df = df[(buses >= MIN_BUSES) & df["settings"].isin(settings) & ~df["status"].isin(BAD_STATUSES)]
    df = df.reset_index(drop=True)
    settings = [s for s in settings if s in set(df["settings"])]
    keys = ["instance", "tal"] if "tal" in df.columns else ["instance"]

    df = keep_groups(df, keys, df.groupby(keys)["settings"].nunique() >= min(2, len(settings)))
    g = df.assign(_t=df["time_MP"].fillna(0) + df["time_SP_N1"].fillna(0)).groupby(keys)
    df = keep_groups(df, keys, np.fmax(g["_t"].max(), g["time_solve"].max()) >= MIN_TIME)
    print(f"tables: {df.groupby(keys).ngroups} instances")

    (exp / "comparison_short.md").write_text(short_table_md(df, keys, settings, labels))
    (exp / "comparison_extra.tex").write_text(extra_table_tex(df, settings, labels))

    # Plots use instances with at least one OPTIMAL setting and results for all settings.
    df = keep_groups(df, keys, df.assign(_o=df["status"] == "OPTIMAL").groupby(keys)["_o"].any())
    df = keep_groups(df, keys, df.groupby(keys)["settings"].nunique() >= len(settings))
    print(f"plots: {df.groupby(keys).ngroups} instances")

    col = lambda c: df[c].fillna(0) if c in df.columns else pd.Series(0, index=df.index)
    df["_sp_solve_time"] = (col("time_SP_N1_primal") + col("time_SP_N1_dual") + col("time_SP_N1_fallback")
                            + col("time_SP_N1_refactorization"))
    df["_sp_dual_time"] = col("time_SP_N1_dual") + 0.5 * col("time_SP_N1_fallback")

    pp = exp / "performance_profiles"
    for target, sub in [("time_solve", "total"), ("time_SP", "SP"), ("_sp_solve_time", "solve_time")]:
        performance_profile(df, target, pp / sub / "factor.pdf", settings, labels, styles, height_cm,
                            legend_left=sub == "solve_time")

    # SP solve time of the settings with dual solving, SLE^{1/2, 2/3}.
    sle = next((s for s in settings if family(s) == "LEQ"), None)
    if sle is not None:
        scaled = [s for s in settings if df.loc[df["settings"] == s, "_sp_dual_time"].sum()
                  / max(df.loc[df["settings"] == s, "time_SP"].sum(), 1e-12) > 0.01]
        parts = [df[df["settings"].isin(scaled)]]
        dash_1dot, dash_2dot = (0, (5, 1.5, 1, 1.5)), (0, (5, 1.5, 1, 1.5, 1, 1.5))
        for tag, exponent, color, dash in [("23", 2.0 / 3.0, "#666666", dash_1dot), ("12", 1.0 / 2.0, "#C8C8C8", dash_2dot)]:
            virtual = f"{sle}__scaled_{tag}"
            rows = df[df["settings"] == sle].copy()
            rows["settings"] = virtual
            rows["_sp_solve_time"] = rows["_sp_solve_time"].pow(exponent)
            parts.append(rows)
            scaled.append(virtual)
            labels[virtual] = rf"{labels[sle]}$^{{{tag[0]}/{tag[1]}}}$"
            styles[virtual] = {"color": color, "linestyle": dash, "marker": "x"}
        df_scaled = pd.concat(parts, ignore_index=True)
        performance_profile(df_scaled, "_sp_solve_time", pp / "solve_time_scaled" / "factor.pdf",
                            scaled, labels, styles, height_cm)
        runtime_scatter(df_scaled, "_sp_solve_time", pp / "solve_time_scaled" / "scatter_sorted_by_sle_runtime.pdf",
                        sle, [s for s in scaled if s != sle], labels, styles, height_cm)

    df = add_time_components(df)
    (exp / "time_decomp").mkdir(exist_ok=True)
    (exp / "time_decomp" / "table.tex").write_text(time_table_tex(df, settings, labels))
    time_boxplot(df, exp / "time_decomp_boxplot" / "detailed.pdf", settings, labels)
    time_boxplot(df, exp / "time_decomp_boxplot" / "detailed_relative.pdf", settings, labels, relative=True)


if __name__ == "__main__":
    main()
