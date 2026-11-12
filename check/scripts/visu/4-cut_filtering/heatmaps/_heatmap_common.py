"""Shared helpers for the cut-filtering grid-search heatmaps.

Each per-experiment call script in this folder is intentionally tiny: it just
imports :func:`render_heatmap` and calls it with the experiment subdir,
baseline setting name, and a flag for the single-axis pb experiments.

All data loading, cleanup, axis synthesis and plotting lives here so the call
scripts stay free of duplicated boilerplate.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

# --- Path wiring -----------------------------------------------------------
# This file lives at: <repo>/check/scripts/visu/4-cut_filtering/heatmaps/
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[4]          # -> repo root (4-cut_filtering)
VISU_DIR = SCRIPT_DIR.parents[1]              # -> check/scripts/visu/
if str(VISU_DIR) not in sys.path:
    sys.path.insert(0, str(VISU_DIR))

from heatmap import grid_search_heatmap, grid_search_heatmap_stacked   # noqa: E402
from reader import read_experiment, cleanup_data                        # noqa: E402

# --- Configuration ---------------------------------------------------------
RESULTS_ROOT = PROJECT_ROOT / "check/results"
MIN_BUSES = 73
MIN_TIME = 30.0
TARGET = "time_solve"

# Setting-name pattern: group 1 -> X (CF), group 2 -> Y (SP).
_PATTERN_2D = r"cf(\d+)_sp(\d+)"

# Single-axis (pb) experiments have setting names like ``sp005`` / ``sp10``.
# We synthesise ``cf<X>_sp00`` so the existing 2-D helper renders them as a
# single row (y = 0) keyed by SP on the X axis.
_PATTERN_PB = re.compile(r"^sp(\d+)$")


def _filter_min_time(df, min_time: float, target: str):
    """Drop instance groups whose max ``target`` across all settings is below ``min_time``."""
    if min_time <= 0:
        return df
    group_cols = ["instance", "tal"] if "tal" in df.columns else ["instance"]
    keep = df.groupby(group_cols)[target].max() >= min_time
    keep_keys = keep[keep].reset_index()[group_cols]
    return df.merge(keep_keys, on=group_cols)


def _rename_pb(df, baseline: str | None):
    """Map ``sp<N>`` -> ``cf<N>_sp00`` in place; return new df and remapped baseline."""
    mapping = {}
    for s in df["settings"].unique():
        m = _PATTERN_PB.fullmatch(str(s))
        if m:
            mapping[s] = f"cf{int(m.group(1)):03d}_sp00"
    if not mapping:
        return df, baseline
    df = df.copy()
    df["settings"] = df["settings"].map(lambda s: mapping.get(s, s))
    if baseline is not None and baseline in mapping:
        baseline = mapping[baseline]
    return df, baseline


def render_heatmap(
    experiment_subdir: str,
    *,
    baseline: str | None = None,
    is_pb: bool = False,
    xlabel: str | None = None,
    ylabel: str | None = None,
    output_name: str = "heatmap_time_solve.pdf",
    include_settings: list[str] | None = None,
    pattern: str | None = None,
    common_only: bool = True,
    penalize_nonoptimal: bool = False,
    min_time: float | None = None,
) -> Path | None:
    """Render the grid-search heatmap for one cut-filtering experiment.

    Returns the output path on success, or ``None`` if the experiment has no
    results yet (so the script can be invoked in a loop without crashing).
    """
    exp_dir = RESULTS_ROOT / experiment_subdir
    if not exp_dir.is_dir() or not any(exp_dir.iterdir()):
        print(f"warning: no results in {exp_dir}; skipping heatmap")
        return None

    try:
        df = read_experiment(exp_dir)
    except FileNotFoundError:
        print(f"warning: no CSV data under {exp_dir}; skipping heatmap")
        return None
    if df.empty:
        print(f"warning: empty dataframe for {experiment_subdir}; skipping heatmap")
        return None

    if include_settings is not None:
        df = df[df["settings"].isin(include_settings)].reset_index(drop=True)
        if df.empty:
            print(f"warning: no rows after include_settings filter for {experiment_subdir}")
            return None

    if is_pb:
        df, baseline = _rename_pb(df, baseline)
        # pb experiments only vary SP; relabel axes accordingly.
        xlabel = xlabel if xlabel is not None else r"SP threshold (\%)"
        ylabel = ylabel if ylabel is not None else ""
    else:
        xlabel = xlabel if xlabel is not None else r"CF threshold (\%)"
        ylabel = ylabel if ylabel is not None else r"SP threshold (\%)"

    if not penalize_nonoptimal:
        df, _ = cleanup_data(df, exclude_all_timelimit=False, baseline=baseline)

    min_time_eff = MIN_TIME if min_time is None else min_time
    df = _filter_min_time(df, min_time_eff, target=TARGET)

    out = exp_dir / output_name
    grid_search_heatmap(
        df,
        output=out,
        pattern=pattern if pattern is not None else _PATTERN_2D,
        target=TARGET,
        baseline=baseline,
        common_only=common_only,
        penalize_nonoptimal=penalize_nonoptimal,
        xlabel=xlabel,
        ylabel=ylabel,
    )
    print(f"Wrote {out}")
    return out


def render_heatmap_stacked(
    experiment_subdir: str,
    *,
    groups: list[tuple[str, list[str]]],
    baseline: str | None = None,
    output_name: str,
    pattern: str | None = None,
    common_only: bool = True,
    penalize_nonoptimal: bool = False,
    min_time: float | None = None,
    suptitle: str | None = None,
) -> Path | None:
    """Render a stacked (per-group) 2-column heatmap grid."""
    exp_dir = RESULTS_ROOT / experiment_subdir
    if not exp_dir.is_dir() or not any(exp_dir.iterdir()):
        print(f"warning: no results in {exp_dir}; skipping stacked heatmap")
        return None

    try:
        df = read_experiment(exp_dir)
    except FileNotFoundError:
        print(f"warning: no CSV data under {exp_dir}; skipping stacked heatmap")
        return None
    if df.empty:
        print(f"warning: empty dataframe for {experiment_subdir}; skipping stacked heatmap")
        return None

    all_include = sorted({s for _, inc in groups for s in inc})
    df = df[df["settings"].isin(all_include)].reset_index(drop=True)
    if df.empty:
        print(f"warning: no rows after include filter for {experiment_subdir}")
        return None

    if not penalize_nonoptimal:
        df, _ = cleanup_data(df, exclude_all_timelimit=False, baseline=baseline)

    min_time_eff = MIN_TIME if min_time is None else min_time
    df = _filter_min_time(df, min_time_eff, target=TARGET)

    out = exp_dir / output_name
    grid_search_heatmap_stacked(
        df,
        output=out,
        groups=groups,
        pattern=pattern if pattern is not None else _PATTERN_2D,
        target=TARGET,
        baseline=baseline,
        common_only=common_only,
        penalize_nonoptimal=penalize_nonoptimal,
        suptitle=suptitle,
    )
    print(f"Wrote {out}")
    return out
