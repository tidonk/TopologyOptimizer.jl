"""Shared helpers for 6-seq_warmstarting grid-search heatmaps.

The ``1-tal1`` sweep is a 2-D grid over the pool size (``n``) and the
weight scheme (``w``). The weight scheme is categorical, so we synthesise
``n{N:02d}_w{IDX:02d}`` from each ``rscots_n{N}_w_{W}`` setting and feed it
to the existing two-integer-group heatmap helper. The Y-axis tick labels
are restored to the weight names afterwards via ``matplotlib`` post-processing.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[4]
VISU_DIR = SCRIPT_DIR.parents[1]
if str(VISU_DIR) not in sys.path:
    sys.path.insert(0, str(VISU_DIR))
if str(SCRIPT_DIR.parent) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR.parent))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from heatmap import grid_search_heatmap, _shifted_geomean  # noqa: E402
from reader import read_experiment, cleanup_data  # noqa: E402
from tables import TIME_LIMIT  # noqa: E402
from _seqwarmstart_common import W_ORDER  # noqa: E402

RESULTS_ROOT = PROJECT_ROOT / "check/results"
MIN_BUSES = 73
MIN_TIME = 30.0
TARGET = "time_solve"

_PATTERN_SRC = re.compile(r"^rscots_n(\d+)_w_([a-z]+)$")
_PATTERN_RENAMED = r"n(\d+)_w(\d+)"


def _filter_min_time(df, min_time: float, target: str):
    if min_time <= 0:
        return df
    group_cols = ["instance", "tal"] if "tal" in df.columns else ["instance"]
    keep = df.groupby(group_cols)[target].max() >= min_time
    keep_keys = keep[keep].reset_index()[group_cols]
    return df.merge(keep_keys, on=group_cols)


def _rename_grid(df, baseline: str | None):
    """Map ``rscots_n{N}_w_{W}`` -> ``n{N:02d}_w{IDX:02d}`` in place."""
    mapping: dict[str, str] = {}
    for s in df["settings"].unique():
        m = _PATTERN_SRC.match(str(s))
        if not m:
            continue
        n_val = int(m.group(1))
        w_name = m.group(2)
        if w_name not in W_ORDER:
            continue
        w_idx = W_ORDER.index(w_name) + 1   # 1..len(W_ORDER), keeps zero free
        mapping[s] = f"n{n_val:02d}_w{w_idx:02d}"
    if not mapping:
        return df, baseline
    df = df.copy()
    df["settings"] = df["settings"].map(lambda s: mapping.get(s, s))
    if baseline is not None and baseline in mapping:
        baseline = mapping[baseline]
    return df, baseline


_TARGET_TITLE: dict[str, str] = {
    "time_solve": r"Geom.\ mean of $t_{\mathrm{solve}}$ [s]",
    "time_with_warmstart": r"Geom.\ mean of $t_{\mathrm{solve}} + t_{\mathrm{ws}}$ [s]",
}


def _relabel_y_axis(fig: plt.Figure, out: Path):
    """Replace numeric y tick labels (1..|W_ORDER|) with the weight names and
    save the final figure (closes the figure afterwards)."""
    for ax in fig.axes:
        labels = [t.get_text() for t in ax.get_yticklabels()]
        new_labels = []
        for lbl in labels:
            try:
                idx = int(lbl)
                new_labels.append(W_ORDER[idx - 1] if 1 <= idx <= len(W_ORDER) else lbl)
            except (ValueError, TypeError):
                new_labels.append(lbl)
        ax.set_yticklabels(new_labels)
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), bbox_inches="tight")
    plt.close(fig)


def render_heatmap(
    experiment_subdir: str,
    *,
    baseline: str | None = None,
    output_name: str = "heatmap_time_solve.pdf",
    common_only: bool = True,
    penalize_nonoptimal: bool = False,
    min_time: float | None = None,
    target: str = TARGET,
    df_transform=None,
) -> Path | None:
    """Render the n x w grid-search heatmap for a 6-seq_warmstart experiment."""
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

    if df_transform is not None:
        df = df_transform(df)

    df, baseline = _rename_grid(df, baseline)

    # Compute baseline geo mean BEFORE the pattern filter (while baseline setting is still present)
    baseline_gm: float | None = None
    if baseline is not None:
        bl_vals = np.asarray(df.loc[df["settings"] == baseline, target], dtype=float)
        bl_vals = bl_vals[np.isfinite(bl_vals) & (bl_vals < TIME_LIMIT)]
        bl_vals = np.clip(bl_vals, 0, None)
        if len(bl_vals) > 0:
            baseline_gm = _shifted_geomean(bl_vals)
    title_rel = rf"Relative to \texttt{{{baseline or 'baseline'}}} (\%)"

    df = df[df["settings"].str.match(_PATTERN_RENAMED)].reset_index(drop=True)
    if df.empty:
        print(f"warning: no grid settings after rename for {experiment_subdir}")
        return None

    if not penalize_nonoptimal:
        df, _ = cleanup_data(df, exclude_all_timelimit=False, baseline=baseline)

    min_time_eff = MIN_TIME if min_time is None else min_time
    df = _filter_min_time(df, min_time_eff, target=target)

    out = exp_dir / output_name
    title_abs = _TARGET_TITLE.get(target, rf"Geom.\ mean of $t_{{{target}}}$ [s]")
    fig = grid_search_heatmap(
        df,
        output=out,
        pattern=_PATTERN_RENAMED,
        target=target,
        baseline=baseline,
        common_only=common_only,
        penalize_nonoptimal=penalize_nonoptimal,
        xlabel=r"pool size $n$",
        ylabel=r"weight scheme",
        title_abs=title_abs,
        title_rel=title_rel,
        baseline_value=baseline_gm,
    )
    _relabel_y_axis(fig, out)
    print(f"Wrote {out}")
    return out
