"""Publication-ready heatmaps.

Two layouts are supported:

* :func:`grid_search_heatmap` / :func:`grid_search_heatmap_stacked` --- the
  legacy two-parameter sweep (e.g.\ ``cf001_sp025``).  Setting names are
  parsed with a regex capturing two integer groups (X, Y).
* :func:`method_level_heatmap` / :func:`method_level_heatmap_2x2` --- the
  partial-decomposition family, with a categorical y-axis (scoring method)
  and an integer x-axis (subproblem count / percentage).  Each cell is
  built from the shared :class:`selection.InstanceSelection`, guaranteeing
  numerical agreement with the comparison tables produced by
  :func:`tables.run_full_comparison`.
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import patches as mpatches

from reader import read_experiment
from selection import InstanceSelection, build_selection
from style import setup_latex_style
from tables import TIME_LIMIT


GEO_SHIFT = 10.0  # standard MIPLIB shift (seconds); matches tables.GEO_SHIFT


def _shifted_geomean(values: np.ndarray, shift: float = GEO_SHIFT) -> float:
    """Shifted geometric mean: exp(mean(log(v + s))) - s."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return np.nan
    return float(np.exp(np.mean(np.log(v + shift))) - shift)


def _parse_grid(settings: list[str], pattern: str) -> dict[str, tuple[int, int]]:
    """Return {setting_name: (x, y)} parsed from ``pattern`` (two int groups)."""
    rx = re.compile(pattern)
    out: dict[str, tuple[int, int]] = {}
    for s in settings:
        m = rx.fullmatch(s)
        if m:
            out[s] = (int(m.group(1)), int(m.group(2)))
    return out


def _build_grid(
    df: pd.DataFrame,
    pattern: str,
    target: str,
    baseline: str | None,
    common_only: bool,
    penalize_nonoptimal: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, list[int], list[int], int, int]:
    """Compute geomean and N-solved per (x, y) cell using the same
    preprocessing as ``comparison_short_table`` (tables.py):

    * Build a (instance[, tal]) x settings pivot of ``target``.
    * Rows with ``status == TIME_LIMIT`` are set to ``TIME_LIMIT``.
    * If ``common_only`` is True (the default), only instances finite and
      below ``TIME_LIMIT`` for *every* setting in the grid contribute to
      the geometric mean.
    * Geometric mean is the shifted geometric mean
      ``exp(mean(log(v + s))) - s`` with ``s = 10`` — matching the
      markdown/LaTeX comparison table (see ``tables.GEO_SHIFT``).

    Returns
    -------
    (gm, nsolved, x_levels, y_levels, n_common, n_total)
        ``gm``/``nsolved`` are DataFrames indexed by y_levels, columns x_levels.
        ``n_common`` is the number of instances solved by all settings.
        ``n_total`` is the total number of instance groups in the grid.
    """
    settings_all = sorted(df["settings"].unique())
    coords = _parse_grid(settings_all, pattern)
    if not coords:
        raise ValueError(f"No settings matched pattern {pattern!r}.")
    settings = sorted(coords)
    df = df[df["settings"].isin(settings)].copy()

    has_tal = "tal" in df.columns
    group_cols = ["instance", "tal"] if has_tal else ["instance"]

    pivoted = df.pivot_table(
        index=group_cols, columns="settings", values=target, aggfunc="first"
    ).reindex(columns=settings)

    status = df.pivot_table(
        index=group_cols, columns="settings", values="status", aggfunc="first"
    ).reindex(columns=settings)

    vals = pivoted.apply(pd.to_numeric, errors="coerce")
    if penalize_nonoptimal:
        # All non-OPTIMAL and missing rows are treated as hard time-limit hits.
        vals = vals.where(status.eq("OPTIMAL"), TIME_LIMIT).fillna(TIME_LIMIT)
    else:
        # Preserve historical behavior: only explicit TIME_LIMIT rows are clipped.
        for s in settings:
            if s not in status.columns:
                continue
            mask = status[s].isin({"TIME_LIMIT", "MEMORY_LIMIT"})
            if mask.any():
                vals.loc[mask[mask].index, s] = TIME_LIMIT

    n_total = len(vals)
    all_solved = (vals < TIME_LIMIT).all(axis=1) & vals.notna().all(axis=1)
    n_common = int(all_solved.sum())
    solved = vals[all_solved]

    x_levels = sorted({xy[0] for xy in coords.values()})
    y_levels = sorted({xy[1] for xy in coords.values()})

    gm = pd.DataFrame(np.nan, index=y_levels, columns=x_levels, dtype=float)
    nsolved = pd.DataFrame(0, index=y_levels, columns=x_levels, dtype=int)

    for s, (x, y) in coords.items():
        n_s = int(((vals[s] < TIME_LIMIT) & vals[s].notna()).sum())
        nsolved.loc[y, x] = n_s
        if common_only:
            col = solved[s].dropna()
        else:
            col = vals[s].dropna()
        if not col.empty:
            gm.loc[y, x] = _shifted_geomean(col.clip(lower=0).to_numpy())

    return gm, nsolved, x_levels, y_levels, n_common, n_total


def grid_search_heatmap(
    df: pd.DataFrame,
    output: str | Path,
    pattern: str = r"cf(\d+)_sp(\d+)",
    target: str = "time_solve",
    baseline: str | None = None,
    common_only: bool = True,
    penalize_nonoptimal: bool = False,
    xlabel: str = r"CF threshold (\%)",
    ylabel: str = r"SP threshold (\%)",
    title_abs: str = r"Geom.\ mean of $t_{\mathrm{solve}}$ [s]",
    title_rel: str = r"Relative to baseline (\%)",
    cmap_abs: str = "viridis_r",
    cmap_rel: str = "RdBu_r",
    baseline_value: float | None = None,
    xtick_labels: list[str] | None = None,
    ytick_labels: list[str] | None = None,
) -> plt.Figure:
    """Render a side-by-side heatmap (absolute geomean and relative-to-baseline %).

    Preprocessing mirrors ``comparison_short_table`` by default (common_only).
    If ``penalize_nonoptimal`` is true, non-OPTIMAL and missing rows are set
    to ``TIME_LIMIT`` before aggregation (comparison_short_all behavior).
    The setting-name pattern must capture two integer groups; group 1 -> X,
    group 2 -> Y.
    """
    setup_latex_style()
    output = Path(output)

    gm, nsolved, x_levels, y_levels, n_common, n_total = _build_grid(
        df, pattern=pattern, target=target, baseline=baseline,
        common_only=common_only,
        penalize_nonoptimal=penalize_nonoptimal,
    )

    # Baseline coordinates (if any)
    base_xy: tuple[int, int] | None = None
    if baseline is not None:
        coords = _parse_grid([baseline], pattern)
        if baseline in coords:
            base_xy = coords[baseline]

    # --- Relative grid (% change vs baseline) ----------------------------
    rel: pd.DataFrame | None = None
    if base_xy is not None and not np.isnan(gm.loc[base_xy[1], base_xy[0]]):
        b = gm.loc[base_xy[1], base_xy[0]]
        rel = (gm - b) / b * 100.0
    elif base_xy is None and baseline_value is not None and np.isfinite(baseline_value) and baseline_value > 0:
        rel = (gm - baseline_value) / baseline_value * 100.0

    # --- Figure layout ---------------------------------------------------
    ncols = 2 if rel is not None else 1
    fig, axes = plt.subplots(
        1, ncols,
        figsize=(5.5 * ncols, 4.5),
        squeeze=False,
        gridspec_kw={"wspace": 0.45},
    )
    axes = axes[0]

    # When common_only, every cell's geomean uses the same n_common instances;
    # replace per-setting counts with the true included count.
    cell_n = nsolved.copy()
    if common_only:
        cell_n[:] = n_common

    _draw_heatmap(
        axes[0], gm, cell_n, x_levels, y_levels,
        cmap=cmap_abs, title=title_abs, fmt="{:.1f}",
        xlabel=xlabel, ylabel=ylabel, baseline_xy=base_xy,
        xtick_labels=xtick_labels, ytick_labels=ytick_labels,
    )
    if rel is not None:
        vmax = float(np.nanmax(np.abs(rel.to_numpy())))
        _draw_heatmap(
            axes[1], rel, cell_n, x_levels, y_levels,
            cmap=cmap_rel, title=title_rel, fmt="{:+.1f}\\%",
            xlabel=xlabel, ylabel=ylabel, baseline_xy=base_xy,
            vmin=-vmax, vmax=vmax,
            xtick_labels=xtick_labels, ytick_labels=ytick_labels,
        )

    if common_only:
        subtitle = rf"$n = {n_common}$ instances (solved by all settings)"
    elif penalize_nonoptimal:
        subtitle = rf"$n = {n_total}$ instances (all; non-OPTIMAL/missing set to {TIME_LIMIT}s)"
    else:
        subtitle = rf"$n = {n_total}$ instances"

    fig.suptitle(subtitle, fontsize=10, y=1.02)
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    fig.savefig(output.with_suffix(".png"), bbox_inches="tight")
    return fig


def _draw_heatmap(
    ax,
    data: pd.DataFrame,
    nsolved: pd.DataFrame,
    x_levels: list[int],
    y_levels: list[int],
    cmap: str,
    title: str,
    fmt: str,
    xlabel: str,
    ylabel: str,
    baseline_xy: tuple[int, int] | None,
    vmin: float | None = None,
    vmax: float | None = None,
    xtick_labels: list[str] | None = None,
    ytick_labels: list[str] | None = None,
) -> None:
    """Render a single heatmap panel with annotations."""
    arr = data.to_numpy(dtype=float)
    # imshow: origin='lower' so y_levels increase upward
    im = ax.imshow(
        arr, origin="lower", aspect="auto",
        cmap=cmap, vmin=vmin, vmax=vmax,
        interpolation="nearest",
    )

    if xtick_labels is None:
        xtick_labels = [str(x) for x in x_levels]
    ax.set_xticks(range(len(x_levels)), xtick_labels)
    if ytick_labels is None:
        ytick_labels = [str(y) for y in y_levels]
    ax.set_yticks(range(len(y_levels)), ytick_labels)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(False)

    # Annotate each cell with value (top) and n=... (smaller, faded, below)
    norm = im.norm
    cmap_obj = im.cmap
    for iy, y in enumerate(y_levels):
        for ix, x in enumerate(x_levels):
            val = arr[iy, ix]
            if not np.isfinite(val):
                continue
            n = int(nsolved.loc[y, x])
            rgba = cmap_obj(norm(val))
            luminance = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
            text_color = "white" if luminance < 0.5 else "black"
            ax.text(
                ix, iy + 0.08,
                fmt.format(val),
                ha="center", va="center",
                fontsize=9, color=text_color,
            )
            ax.text(
                ix, iy - 0.28,
                rf"$n{{=}}{n}$",
                ha="center", va="center",
                fontsize=6, color=text_color, alpha=0.55,
            )

    # Baseline outline
    if baseline_xy is not None:
        bx, by = baseline_xy
        if bx in x_levels and by in y_levels:
            ix = x_levels.index(bx)
            iy = y_levels.index(by)
            ax.add_patch(mpatches.Rectangle(
                (ix - 0.5, iy - 0.5), 1, 1,
                fill=False, edgecolor="black", linewidth=1.8,
            ))

    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


def geomeans_from_short_table(md_path, coord_fn):
    """Parse a ``comparison_short_all.md`` into ``(gm_map, n_map, n_note)``.

    ``coord_fn(plain_label) -> (x, y) | None`` maps each row's (de-LaTeX'd) label
    to a grid coordinate; rows it returns ``None`` for are skipped (baseline-only
    refs, non-grid settings).  Because the values come straight from the table,
    a heatmap built from them is guaranteed to match it exactly — and a setting
    the table dropped (too many cut errors) is simply absent, i.e. a blank cell.
    """
    from pathlib import Path as _P
    text = _P(md_path).read_text()
    lines = [l for l in text.splitlines() if l.startswith("|")]
    if len(lines) < 3:
        return {}, {}, ""
    cols = [c.strip() for c in lines[0].strip("|").split("|")]
    gi = cols.index("Geo. Mean") if "Geo. Mean" in cols else -2
    ni = cols.index("Solved") if "Solved" in cols else 1
    gm, nn = {}, {}
    for l in lines[2:]:
        c = [x.strip() for x in l.strip("|").split("|")]
        if len(c) <= gi:
            continue
        plain = re.sub(r"\\texttt\{([^{}]*)\}", r"\1", c[0]).replace("\\%", "%").replace("\\_", "_").strip()
        xy = coord_fn(plain)
        if xy is None:
            continue
        g = re.search(r"[-+]?\d+\.?\d*", c[gi].replace("**", ""))
        n = re.search(r"\d+", c[ni]) if ni < len(c) else None
        if g:
            gm[xy] = float(g.group())
            nn[xy] = int(n.group()) if n else 0
    m = re.search(r"\*n = (\d+)", text)
    note = rf"$n = {m.group(1)}$ instances (all; non-OPTIMAL/missing $\to$ TL)" if m else ""
    return gm, nn, note


def grid_heatmap_from_values(
    gm_map: dict,
    n_map: dict,
    output: str | Path,
    *,
    x_levels: list[int],
    y_levels: list[int],
    baseline_xy: tuple[int, int] | None = None,
    baseline_value: float | None = None,
    xlabel: str = "",
    ylabel: str = "",
    xtick_labels: list[str] | None = None,
    ytick_labels: list[str] | None = None,
    title_abs: str = r"Geom.\ mean of $t_{\mathrm{solve}}$ [s]",
    title_rel: str = r"Relative to baseline (\%)",
    cmap_abs: str = "viridis_r",
    cmap_rel: str = "RdBu_r",
    n_note: str = "",
) -> "plt.Figure":
    """Render a grid heatmap from **pre-computed** per-cell geomeans.

    Unlike :func:`grid_search_heatmap` (which recomputes the geomean from raw
    runs), this takes ``gm_map``/``n_map`` = ``{(x, y): value}`` directly — use
    it to guarantee the cells equal a table's already-published ``Geo. Mean``
    (same instance selection, same cut-error dropping/penalising).  Missing
    cells are left blank.
    """
    setup_latex_style()
    output = Path(output)
    gm = pd.DataFrame(index=y_levels, columns=x_levels, dtype=float)
    nn = pd.DataFrame(0, index=y_levels, columns=x_levels, dtype=int)
    for (x, y), v in gm_map.items():
        if x in x_levels and y in y_levels and v is not None:
            gm.loc[y, x] = float(v)
    for (x, y), v in n_map.items():
        if x in x_levels and y in y_levels and v is not None:
            nn.loc[y, x] = int(v)

    rel = None
    b = None
    if baseline_xy is not None and baseline_xy[0] in x_levels and baseline_xy[1] in y_levels:
        b = gm.loc[baseline_xy[1], baseline_xy[0]]
    elif baseline_value is not None:
        b = baseline_value                       # off-grid reference (e.g. "No PD")
    if b is not None and np.isfinite(b) and b > 0:
        rel = (gm - b) / b * 100.0

    ncols = 2 if rel is not None else 1
    fig, axes = plt.subplots(1, ncols, figsize=(5.5 * ncols, 4.5),
                             squeeze=False, gridspec_kw={"wspace": 0.45})
    axes = axes[0]
    _draw_heatmap(axes[0], gm, nn, x_levels, y_levels, cmap=cmap_abs,
                  title=title_abs, fmt="{:.1f}", xlabel=xlabel, ylabel=ylabel,
                  baseline_xy=baseline_xy, xtick_labels=xtick_labels,
                  ytick_labels=ytick_labels)
    if rel is not None:
        vmax = float(np.nanmax(np.abs(rel.to_numpy())))
        _draw_heatmap(axes[1], rel, nn, x_levels, y_levels, cmap=cmap_rel,
                      title=title_rel, fmt="{:+.1f}\\%", xlabel=xlabel,
                      ylabel=ylabel, baseline_xy=baseline_xy,
                      vmin=-vmax, vmax=vmax, xtick_labels=xtick_labels,
                      ytick_labels=ytick_labels)
    if n_note:
        fig.suptitle(n_note, fontsize=10, y=1.02)
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    fig.savefig(output.with_suffix(".png"), bbox_inches="tight")
    return fig


def grid_search_heatmap_stacked(
    df: pd.DataFrame,
    output: str | Path,
    groups: list[tuple[str, list[str]]],
    pattern: str = r"cf(\d+)_sp(\d+)",
    target: str = "time_solve",
    baseline: str | None = None,
    common_only: bool = True,
    penalize_nonoptimal: bool = False,
    xlabel: str = r"CF threshold (\%)",
    ylabel: str = r"Number of cuts per callback",
    title_abs: str = r"Geom.\ mean of $t_{\mathrm{solve}}$ [s]",
    title_rel: str = r"Relative to baseline (\%)",
    cmap_abs: str = "viridis_r",
    cmap_rel: str = "RdBu_r",
    y_inf_at: int | None = 100,
    suptitle: str | None = None,
) -> None:
    """Render a (len(groups) x 2) grid of heatmaps stacked vertically.

    Each group is rendered as one row (absolute geomean + relative-to-baseline)
    using only the settings in its include list. The y-axis label `100` is
    replaced with `\\infty` if ``y_inf_at`` matches the actual level.

    Parameters
    ----------
    groups : list of (label, include_settings)
        Each tuple drives one row; ``label`` becomes the row's left-side title.
    """
    setup_latex_style()
    output = Path(output)

    nrows = len(groups)
    fig, axes = plt.subplots(
        nrows, 2, figsize=(11.0, 4.5 * nrows), squeeze=False,
    )

    for row, (label, include) in enumerate(groups):
        sub = df[df["settings"].isin(include)].copy()
        if sub.empty:
            for c in range(2):
                axes[row, c].axis("off")
                axes[row, c].set_title(f"{label}: no data")
            continue

        gm, nsolved, x_levels, y_levels, n_common, n_total = _build_grid(
            sub, pattern=pattern, target=target, baseline=baseline,
            common_only=common_only, penalize_nonoptimal=penalize_nonoptimal,
        )

        base_xy: tuple[int, int] | None = None
        if baseline is not None and baseline in include:
            coords = _parse_grid([baseline], pattern)
            base_xy = coords.get(baseline)

        rel: pd.DataFrame | None = None
        if base_xy is not None and base_xy[1] in gm.index and base_xy[0] in gm.columns \
                and not np.isnan(gm.loc[base_xy[1], base_xy[0]]):
            b = gm.loc[base_xy[1], base_xy[0]]
            rel = (gm - b) / b * 100.0

        cell_n = nsolved.copy()
        if common_only:
            cell_n[:] = n_common

        ytick_labels = [
            (r"$\infty$" if (y_inf_at is not None and y == y_inf_at) else str(y))
            for y in y_levels
        ]

        if common_only:
            n_note = rf"$n = {n_common}$ instances (solved by all)"
        elif penalize_nonoptimal:
            n_note = rf"$n = {n_total}$ instances (all; non-OPTIMAL$\to${TIME_LIMIT}s)"
        else:
            n_note = rf"$n = {n_total}$ instances"

        _draw_heatmap(
            axes[row, 0], gm, cell_n, x_levels, y_levels,
            cmap=cmap_abs, title=f"{label}: {title_abs}", fmt="{:.1f}",
            xlabel=xlabel, ylabel=ylabel, baseline_xy=base_xy,
            ytick_labels=ytick_labels,
        )
        if rel is not None:
            vmax = float(np.nanmax(np.abs(rel.to_numpy())))
            _draw_heatmap(
                axes[row, 1], rel, cell_n, x_levels, y_levels,
                cmap=cmap_rel, title=f"{label}: {title_rel}",
                fmt="{:+.1f}\\%", xlabel=xlabel, ylabel=ylabel,
                baseline_xy=base_xy, vmin=-vmax, vmax=vmax,
                ytick_labels=ytick_labels,
            )
        else:
            axes[row, 1].axis("off")
            axes[row, 1].set_title(f"{label}: (no baseline)")

        # Per-row n note below the row, on the left panel.
        axes[row, 0].text(
            0.0, -0.18, n_note, transform=axes[row, 0].transAxes,
            ha="left", va="top", fontsize=8,
        )

    if suptitle:
        fig.suptitle(suptitle, fontsize=11, y=1.0)
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    fig.savefig(output.with_suffix(".png"), bbox_inches="tight")
    plt.close(fig)


# =====================================================================
#  Method × level heatmaps (partial-decomposition family of experiments)
# =====================================================================
#
# y-axis: categorical (scoring method) ; x-axis: integer level.
# All data flows through :class:`selection.InstanceSelection` so every cell
# uses exactly the same instance partition as the corresponding table row.

# Default configuration shared with the compare scripts.
DEFAULT_RESULTS_ROOT = Path(__file__).resolve().parents[2] / "results"
DEFAULT_MIN_BUSES = 73
DEFAULT_MIN_TIME = 0.0
DEFAULT_TARGET = "time_solve"

METHOD_ORDER = ["cprox", "div_cos", "div_euc", "eta", "rnd", "sens"]
METHOD_LABELS = {
    "cprox":   r"\textsc{CProx}",
    "div_cos": r"\textsc{div-cos}",
    "div_euc": r"\textsc{div-euc}",
    "eta":     r"\textsc{eta}",
    "rnd":     r"\textsc{rnd}",
    "sens":    r"\textsc{sens}",
}
DEFAULT_LEVELS = [1, 2, 3, 4, 5]


def _build_method_level_grid(
    sel: InstanceSelection,
    *,
    axis_tag: str,
    methods: list[str],
    levels: list[int],
    baseline_setting: str,
    mode: str,
) -> tuple[pd.DataFrame, pd.DataFrame, float, int, int, int]:
    """Build the method × level geomean / N-solved grid from a selection.

    ``mode='optimal'`` uses :attr:`InstanceSelection.optimal` (common-solved).
    ``mode='ran'``     uses :attr:`InstanceSelection.ran` with non-OPTIMAL
    clipped to ``TIME_LIMIT`` for the geomean.

    Returns ``(gm, nsolved, base_gm, base_n, n_stats, n_listing)``.
    """
    rx = re.compile(
        rf"^({'|'.join(re.escape(m) for m in methods)})_{re.escape(axis_tag)}_(\d+)$"
    )
    coords: dict[str, tuple[str, int]] = {}
    for s in sel.settings:
        m = rx.fullmatch(str(s))
        if m and int(m.group(2)) in levels:
            coords[s] = (m.group(1), int(m.group(2)))
    if not coords:
        raise ValueError(
            f"No settings matched <method>_{axis_tag}_<level> for methods={methods}"
        )

    stats_vals = sel.vals_for(mode, penalize_nonoptimal=(mode == "ran"))
    n_stats = len(stats_vals)
    n_listing = len(sel.listing)

    gm = pd.DataFrame(np.nan, index=methods, columns=levels, dtype=float)
    nsolved = pd.DataFrame(0, index=methods, columns=levels, dtype=int)

    for s, (method, level) in coords.items():
        if s in sel.status.columns:
            nsolved.loc[method, level] = int(sel.status[s].eq("OPTIMAL").sum())
        if s in stats_vals.columns and not stats_vals[s].empty:
            col = stats_vals[s].dropna()
            if not col.empty:
                gm.loc[method, level] = _shifted_geomean(
                    col.clip(lower=0).to_numpy()
                )

    base_gm = float("nan")
    base_n = 0
    if baseline_setting in sel.status.columns:
        base_n = int(sel.status[baseline_setting].eq("OPTIMAL").sum())
    if baseline_setting in stats_vals.columns:
        b_col = stats_vals[baseline_setting].dropna()
        if not b_col.empty:
            base_gm = _shifted_geomean(b_col.clip(lower=0).to_numpy())

    return gm, nsolved, base_gm, base_n, n_stats, n_listing


def _load_method_level(
    experiment_subdir: str,
    *,
    axis_tag: str,
    methods: list[str],
    levels: list[int],
    baseline_setting: str,
    common_only: bool,
    penalize_nonoptimal: bool,
    min_time: float,
    min_buses: int = DEFAULT_MIN_BUSES,
    target: str = DEFAULT_TARGET,
    results_root: Path = DEFAULT_RESULTS_ROOT,
    df_transform=None,
) -> "tuple | None":
    """Load results and build the method×level grid for one experiment.

    Returns ``(gm, nsolved, base_gm, base_n, n_stats, n_listing, exp_dir, sel)``
    or ``None`` when the experiment directory is missing / empty.
    """
    exp_dir = Path(results_root) / experiment_subdir
    if not exp_dir.is_dir() or not any(exp_dir.iterdir()):
        print(f"warning: no results in {exp_dir}; skipping")
        return None
    try:
        df = read_experiment(exp_dir)
    except FileNotFoundError:
        print(f"warning: no CSV data under {exp_dir}; skipping")
        return None
    if df.empty:
        print(f"warning: empty dataframe for {experiment_subdir}; skipping")
        return None

    if df_transform is not None:
        df = df_transform(df)

    keep_settings = [f"{m}_{axis_tag}_{lv}" for m in methods for lv in levels]
    if baseline_setting in df["settings"].unique():
        keep_settings.append(baseline_setting)
    keep_settings = [s for s in keep_settings if s in set(df["settings"].unique())]
    if not keep_settings:
        print(f"warning: no matching settings in {experiment_subdir}")
        return None

    sel = build_selection(
        df,
        settings=keep_settings,
        baseline=baseline_setting if baseline_setting in keep_settings else None,
        target=target,
        min_buses=min_buses,
        min_time=min_time,
    )

    mode = "ran" if penalize_nonoptimal else "optimal" if common_only else "ran"
    try:
        gm, nsolved, base_gm, base_n, n_stats, n_listing = _build_method_level_grid(
            sel,
            axis_tag=axis_tag, methods=methods, levels=levels,
            baseline_setting=baseline_setting, mode=mode,
        )
    except ValueError as exc:
        print(f"warning: {exc}; skipping {experiment_subdir}")
        return None
    return gm, nsolved, base_gm, base_n, n_stats, n_listing, exp_dir, sel


def _method_level_subtitle(
    common_only: bool,
    penalize_nonoptimal: bool,
    n_stats: int,
    n_listing: int,
    base_gm: float,
    base_n: int,
    baseline_setting: str,
) -> str:
    """Subtitle string for a method × level heatmap row."""
    parts: list[str] = []
    if common_only:
        parts.append(
            rf"$n = {n_stats}$ of ${n_listing}$ instances "
            r"(solved to OPTIMAL by every setting incl.\ baseline)"
        )
    elif penalize_nonoptimal:
        parts.append(
            rf"$n = {n_stats}$ of ${n_listing}$ instances "
            rf"(all settings ran; non-OPTIMAL/missing $\to$ {TIME_LIMIT}\,s)"
        )
    else:
        parts.append(rf"$n = {n_stats}$ of ${n_listing}$ instances")
    if np.isfinite(base_gm):
        baseline_tex = baseline_setting.replace("_", r"\_")
        parts.append(
            rf"\texttt{{{baseline_tex}}}: "
            rf"geo.\ mean $= {base_gm:.1f}$\,s ($n_{{\mathrm{{solved}}}} = {base_n}$)"
        )
    return " --- ".join(parts)


def _draw_method_level_panels(
    ax_abs,
    ax_rel,
    gm: pd.DataFrame,
    nsolved: pd.DataFrame,
    base_gm: float,
    n_stats: int,
    *,
    levels: list[int],
    methods: list[str],
    common_only: bool,
    xlabel: str,
    ylabel: str,
    title_abs: str,
    title_rel: str,
    cmap_abs: str,
    cmap_rel: str,
    ytick_labels: list[str],
) -> None:
    """Draw one (absolute, relative-to-baseline) panel pair onto two axes.

    Hides ``ax_rel`` when no finite baseline is available.
    """
    cell_n = nsolved.copy()
    if common_only:
        cell_n[:] = n_stats

    _draw_heatmap(
        ax_abs, gm, cell_n, levels, methods,
        cmap=cmap_abs, title=title_abs, fmt="{:.1f}",
        xlabel=xlabel, ylabel=ylabel, baseline_xy=None,
        ytick_labels=ytick_labels,
    )

    if not (np.isfinite(base_gm) and base_gm > 0):
        if ax_rel is not None:
            ax_rel.set_visible(False)
        return

    rel = (gm - base_gm) / base_gm * 100.0
    vmax = float(np.nanmax(np.abs(rel.to_numpy())))
    if not np.isfinite(vmax) or vmax == 0.0:
        vmax = 1.0
    _draw_heatmap(
        ax_rel, rel, cell_n, levels, methods,
        cmap=cmap_rel, title=title_rel, fmt=r"{:+.1f}\%",
        xlabel=xlabel, ylabel=ylabel, baseline_xy=None,
        vmin=-vmax, vmax=vmax, ytick_labels=ytick_labels,
    )


def _default_xlabel(axis_tag: str) -> str:
    return (r"Subproblem count $N$" if axis_tag == "num"
            else r"Subproblem percentage (\%)")


def method_level_heatmap(
    experiment_subdir: str,
    *,
    axis_tag: str,
    output_name: str,
    baseline_setting: str = "cf05_sp005",
    methods: list[str] = METHOD_ORDER,
    levels: list[int] = DEFAULT_LEVELS,
    common_only: bool = True,
    penalize_nonoptimal: bool = False,
    min_time: float | None = None,
    min_buses: int = DEFAULT_MIN_BUSES,
    target: str = DEFAULT_TARGET,
    results_root: Path = DEFAULT_RESULTS_ROOT,
    xlabel: str | None = None,
    ylabel: str = r"Subproblem-selection score",
    title_abs: str = r"Geom.\ mean of $t_{\mathrm{solve}}$ [s]",
    title_rel: str = r"Relative to \texttt{cf05\_sp005} (\%)",
    cmap_abs: str = "viridis_r",
    cmap_rel: str = "RdBu_r",
    suptitle: str | None = None,
    df_transform=None,
) -> Path | None:
    """Render a single method × level heatmap (abs + rel-to-baseline panels)."""
    result = _load_method_level(
        experiment_subdir,
        axis_tag=axis_tag, methods=methods, levels=levels,
        baseline_setting=baseline_setting,
        common_only=common_only, penalize_nonoptimal=penalize_nonoptimal,
        min_time=DEFAULT_MIN_TIME if min_time is None else min_time,
        min_buses=min_buses, target=target, results_root=results_root,
        df_transform=df_transform,
    )
    if result is None:
        return None
    gm, nsolved, base_gm, base_n, n_stats, n_listing, exp_dir, _sel = result

    setup_latex_style()
    has_rel = np.isfinite(base_gm) and base_gm > 0
    ncols = 2 if has_rel else 1
    fig, axes = plt.subplots(
        1, ncols, figsize=(5.5 * ncols, 0.55 * len(methods) + 2.0), squeeze=False,
    )
    ax_abs = axes[0, 0]
    ax_rel = axes[0, 1] if has_rel else None

    _draw_method_level_panels(
        ax_abs, ax_rel, gm, nsolved, base_gm, n_stats,
        levels=levels, methods=methods, common_only=common_only,
        xlabel=xlabel or _default_xlabel(axis_tag), ylabel=ylabel,
        title_abs=title_abs, title_rel=title_rel,
        cmap_abs=cmap_abs, cmap_rel=cmap_rel,
        ytick_labels=[METHOD_LABELS.get(m, m) for m in methods],
    )

    subtitle = _method_level_subtitle(
        common_only, penalize_nonoptimal,
        n_stats, n_listing, base_gm, base_n, baseline_setting,
    )
    if suptitle is not None:
        subtitle = f"{suptitle}\n{subtitle}"
    fig.suptitle(subtitle, fontsize=10, y=1.02)
    fig.tight_layout()

    out = exp_dir / output_name
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")
    return out


def method_level_heatmap_2x2(
    experiment_top: str,
    experiment_bottom: str,
    *,
    axis_tag: str,
    axis_tag_bottom: str | None = None,
    output_name: str,
    baseline_setting: str = "cf05_sp005",
    methods: list[str] = METHOD_ORDER,
    levels_top: list[int] = DEFAULT_LEVELS,
    levels_bottom: list[int] | None = None,
    common_only: bool = True,
    penalize_nonoptimal: bool = False,
    min_time: float | None = None,
    min_buses: int = DEFAULT_MIN_BUSES,
    target: str = DEFAULT_TARGET,
    results_root: Path = DEFAULT_RESULTS_ROOT,
    row_label_top: str = "",
    row_label_bottom: str = "",
    xlabel: str | None = None,
    ylabel: str = r"Subproblem-selection score",
    title_abs: str = r"Geom.\ mean of $t_{\mathrm{solve}}$ [s]",
    title_rel: str = r"Relative to \texttt{cf05\_sp005} (\%)",
    cmap_abs: str = "viridis_r",
    cmap_rel: str = "RdBu_r",
    suptitle: str | None = None,
    output_dir: str | None = None,
    df_transform=None,
) -> Path | None:
    """Render a 2×2 method × level heatmap: two experiments, abs + rel panels."""
    _min_time = DEFAULT_MIN_TIME if min_time is None else min_time
    levels_bot = levels_bottom if levels_bottom is not None else levels_top
    axis_tag_bot = axis_tag_bottom if axis_tag_bottom is not None else axis_tag

    common = dict(
        methods=methods, baseline_setting=baseline_setting,
        common_only=common_only, penalize_nonoptimal=penalize_nonoptimal,
        min_time=_min_time, min_buses=min_buses,
        target=target, results_root=results_root, df_transform=df_transform,
    )
    top = _load_method_level(experiment_top,    axis_tag=axis_tag,     levels=levels_top, **common)
    bot = _load_method_level(experiment_bottom, axis_tag=axis_tag_bot, levels=levels_bot, **common)

    if top is None and bot is None:
        print("warning: no data for either experiment; skipping 2x2 heatmap")
        return None

    setup_latex_style()
    xlabel = xlabel or _default_xlabel(axis_tag)
    ytick_labels = [METHOD_LABELS.get(m, m) for m in methods]
    row_height = 0.55 * len(methods) + 2.0
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 2 * row_height), squeeze=False)

    for row_idx, (row_data, row_levels, row_label) in enumerate([
        (top, levels_top, row_label_top),
        (bot, levels_bot, row_label_bottom),
    ]):
        ax_abs, ax_rel = axes[row_idx]
        if row_data is None:
            ax_abs.set_visible(False)
            ax_rel.set_visible(False)
            continue

        gm, nsolved, base_gm, base_n, n_stats, n_listing, _exp_dir, _sel = row_data
        row_title_abs = (rf"\textbf{{{row_label}}}~" if row_label else "") + title_abs

        _draw_method_level_panels(
            ax_abs, ax_rel, gm, nsolved, base_gm, n_stats,
            levels=row_levels, methods=methods, common_only=common_only,
            xlabel=xlabel, ylabel=ylabel,
            title_abs=row_title_abs, title_rel=title_rel,
            cmap_abs=cmap_abs, cmap_rel=cmap_rel,
            ytick_labels=ytick_labels,
        )

        ax_abs.annotate(
            _method_level_subtitle(
                common_only, penalize_nonoptimal,
                n_stats, n_listing, base_gm, base_n, baseline_setting,
            ),
            xy=(0.5, -0.22), xycoords="axes fraction",
            ha="center", va="top", fontsize=7, style="italic",
        )

    if suptitle is not None:
        fig.suptitle(suptitle, fontsize=10, y=1.01)
    fig.tight_layout()

    # Output dir: explicit > top exp dir > bottom exp dir.
    out_dir = Path(output_dir) if output_dir is not None else (top or bot)[6]
    out = out_dir / output_name
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")
    return out
