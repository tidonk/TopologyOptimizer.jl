"""Stacked bar charts showing time decomposition.

- ``per_instance_time_bar``: one row per instance, two horizontal bars per
  setting, with a visual divider between settings.  Designed for portrait A4.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from reader import clean_instance_name, instance_sort_key
from style import COLORS, setup_latex_style, fmt_setting

# Top-level, non-overlapping time components
DEFAULT_TIME_COMPONENTS = [
    ("time_MP", "Master Problem"),
    ("time_SP_N1", "Subproblems (N-1)"),
]

# Greyscale palette for time components (indexed by component order)
_COMPONENT_COLORS = ["#000000", "#808080", "#B0B0B0", "#D0D0D0"]

# Hatch patterns per setting index (first setting is solid)
_SETTING_HATCHES = [None, "...", "//", "xx", "++", "\\\\", "OO"]


def _order_settings(settings: list[str], baseline: str | None,
                     include_settings: list[str] | None = None) -> list[str]:
    """Return *settings* ordered by *include_settings* (if given), else alphabetically.

    *baseline* is always placed first when present.
    """
    if include_settings is not None:
        # Use include_settings order, restricted to what's actually in settings
        present = set(settings)
        ordered = [s for s in include_settings if s in present]
        # Append any settings not in include_settings (alphabetically)
        extra = sorted(s for s in settings if s not in set(include_settings))
        ordered.extend(extra)
        return ordered
    rest = sorted(s for s in settings if s != baseline)
    if baseline and baseline in settings:
        return [baseline] + rest
    return rest


def per_instance_time_bar(
    df: pd.DataFrame,
    output: str | Path | None = None,
    label_map: dict[str, str] | None = None,
    time_components: list[tuple[str, str]] | None = None,
    bar_height: float = 0.35,
    axes_limits: list[float] | None = None,
    baseline: str | None = None,
    min_time: float = 0.0,
    include_settings: list[str] | None = None,
) -> plt.Figure:
    """Plot a per-instance stacked horizontal bar chart (portrait, A4-ready).

    Each instance gets one row with two horizontal bars (one per setting),
    separated by a thin gap.  A horizontal divider line separates the two
    settings visually.

    Parameters
    ----------
    df : experiment DataFrame (cleaned, all settings present for every instance).
    output : save path.
    label_map : setting name → display label.
    time_components : list of (column, display_label) pairs.
    bar_height : height of each individual bar.
    axes_limits : if given, list of upper x-axis cutoffs (e.g. [600, 3600, inf])
        that splits instances into panels.  Instances are assigned to a panel
        based on their max total time across settings.
    min_time : Drop instances where all settings finish below this cutoff.
    baseline : setting name used for sorting (ascending by total time).
    """
    setup_latex_style()

    if time_components is None:
        time_components = DEFAULT_TIME_COMPONENTS
    if label_map is None:
        label_map = {}

    col_names = [tc[0] for tc in time_components]
    col_labels = [tc[1] for tc in time_components]
    settings = _order_settings(list(df["settings"].unique()), baseline, include_settings)
    n_settings = len(settings)

    # Build display labels
    has_tal = "tal" in df.columns
    group_cols = ["instance", "tal"] if has_tal else ["instance"]

    instances = df.drop_duplicates(subset=group_cols)[group_cols].copy()

    def _total_time(key, setting):
        """Total time of *setting* for a given instance key."""
        if has_tal and len(key) > 1:
            mask = (df["instance"] == key[0]) & (df["tal"] == key[1]) & (df["settings"] == setting)
        else:
            mask = (df["instance"] == key[0]) & (df["settings"] == setting)
        row = df[mask]
        if row.empty:
            return float("inf")
        row = row.iloc[0]
        return sum(float(row.get(c, 0) or 0) for c in col_names)

    # Sort by baseline total time (ascending); fall back to instance number
    sort_setting = baseline if baseline and baseline in settings else settings[0]
    instances["_sort"] = instances.apply(
        lambda r: _total_time(tuple(r[group_cols]), sort_setting), axis=1)
    instances = instances.sort_values("_sort").drop(columns="_sort")
    instance_keys = list(instances.itertuples(index=False, name=None))

    def _display_name(key):
        if has_tal and len(key) > 1 and key[1] not in (None, 0):
            return f"{clean_instance_name(key[0])}\\_tal{int(key[1])}"
        return clean_instance_name(key[0])

    # Compute max total time per instance (across settings) for bucket assignment
    def _max_total(key):
        return max(_total_time(key, s) for s in settings)

    # Drop trivially fast instances
    if min_time > 0:
        instance_keys = [k for k in instance_keys if _max_total(k) >= min_time]

    if axes_limits is None:
        buckets = [(float("inf"), instance_keys)]
    else:
        limits = sorted(axes_limits)
        if limits[-1] != float("inf"):
            limits.append(float("inf"))
        bucket_map: dict[float, list] = {lim: [] for lim in limits}
        for key in instance_keys:
            mt = _max_total(key)
            for lim in limits:
                if mt <= lim:
                    bucket_map[lim].append(key)
                    break
            else:
                bucket_map[limits[-1]].append(key)
        buckets = [(lim, keys) for lim, keys in bucket_map.items() if keys]

    n_panels = len(buckets)

    instance_gap = 0.4   # gap between instance groups
    bar_gap = 0.02       # tiny gap between bars within an instance
    group_height = n_settings * bar_height + (n_settings - 1) * bar_gap
    row_pitch = group_height + instance_gap

    # Fixed physical height per row ensures identical bar sizes across panels
    INCH_PER_ROW = 0.42
    MARGIN_INCHES = 1.5  # top/bottom margins + legend
    panel_counts = [len(keys) for _, keys in buckets]
    total_count = sum(panel_counts)
    fig_height = max(4, total_count * INCH_PER_ROW + MARGIN_INCHES)

    fig, axes = plt.subplots(
        n_panels, 1,
        figsize=(7, fig_height),
        gridspec_kw={"height_ratios": panel_counts} if n_panels > 1 else {},
        squeeze=False,
    )
    axes = [ax for ax in axes[:, 0]]

    for panel_idx, (xlim, panel_keys) in enumerate(buckets):
        ax = axes[panel_idx]
        ax.set_axisbelow(True)
        n_inst = len(panel_keys)

        y_positions = {}
        y_ticks, y_labels = [], []

        for i, key in enumerate(panel_keys):
            base_y = i * row_pitch
            for j in range(n_settings):
                y_positions[(i, j)] = base_y + j * (bar_height + bar_gap)
            row_center = base_y + (group_height - bar_height) / 2
            y_ticks.append(row_center)
            y_labels.append(_display_name(key))

        # Draw bars
        for j, s in enumerate(settings):
            hatch = _SETTING_HATCHES[j % len(_SETTING_HATCHES)]
            for i, key in enumerate(panel_keys):
                if has_tal and len(key) > 1:
                    mask = (df["instance"] == key[0]) & (df["tal"] == key[1]) & (df["settings"] == s)
                else:
                    mask = (df["instance"] == key[0]) & (df["settings"] == s)
                row = df[mask]
                if row.empty:
                    continue
                row = row.iloc[0]
                y = y_positions[(i, j)]
                left = 0.0
                for k, col in enumerate(col_names):
                    val = float(row.get(col, 0) or 0)
                    lbl = None
                    if panel_idx == 0 and i == 0 and j == 0:
                        lbl = col_labels[k]
                    ax.barh(y, val, left=left, height=bar_height,
                            color=_COMPONENT_COLORS[k % len(_COMPONENT_COLORS)],
                            edgecolor="white", linewidth=0.3,
                            hatch=hatch,
                            label=lbl)
                    left += val

        # Symmetric padding: half-gap above first bar group, half-gap below last
        total_y = n_inst * row_pitch - instance_gap
        pad = instance_gap / 2
        ax.set_ylim(-pad, total_y + pad)

        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_labels, fontsize=7)
        ax.invert_yaxis()

        # X-axis: ticks divisible by 300, last tick pinned to panel limit
        if np.isinf(xlim):
            panel_max = 12000
        else:
            panel_max = int(xlim)
        ax.set_xlim(left=0, right=panel_max)
        # Tick step depends on panel range
        if panel_max <= 300:
            step = 60
        elif panel_max <= 3600:
            step = 300
        else:
            step = 1200
        xticks = list(range(0, panel_max + 1, step))
        if xticks[-1] != panel_max:
            xticks.append(panel_max)
        ax.set_xticks(xticks)

        # Dashed vertical line at TIME_LIMIT (10800 s) when visible
        if panel_max >= 10800:
            ax.axvline(10800, color="gray", linewidth=0.8, linestyle="--")
            ax.text(10800, -pad * 0.4, "TL", ha="center", va="top",
                    fontsize=7, color="gray")

        if panel_idx == n_panels - 1:
            ax.set_xlabel(r"Time (s)")

        ax.text(0.99, 0.01, rf"$n = {n_inst}$",
                transform=ax.transAxes, va="bottom", ha="right", fontsize=8, color="gray")

    # Build combined legend from first panel + setting patches
    from matplotlib.patches import Patch
    setting_patches = []
    for j, s in enumerate(settings):
        hatch = _SETTING_HATCHES[j % len(_SETTING_HATCHES)]
        setting_patches.append(
            Patch(facecolor="lightgray", edgecolor="black", linewidth=0.5,
                  hatch=hatch,
                  label=fmt_setting(s, label_map))
        )
    handles, labels = axes[0].get_legend_handles_labels()
    all_handles = handles + setting_patches
    all_labels = labels + [p.get_label() for p in setting_patches]
    axes[-1].legend(all_handles, all_labels, frameon=False, loc="upper center",
                    bbox_to_anchor=(0.5, -0.15), ncol=min(len(all_handles), 4), fontsize=7)

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.08)
    if output:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output)
        if output.suffix == ".pdf":
            fig.savefig(output.with_suffix(".png"))
        print(f"Saved per-instance time bar to {output}")
    return fig


# Greyscale shades for status categories
_STATUS_COLORS = {
    "OPTIMAL": "#333333",
    "TIME_LIMIT": "#888888",
    "DID_NOT_RUN": "#FFFFFF",
    "INFEASIBLE": "#AAAAAA",
    "CUT_ERROR": "#E0E0E0",
}

_STATUS_HATCHES = {
    "OPTIMAL": None,
    "TIME_LIMIT": None,
    "DID_NOT_RUN": None,
    "INFEASIBLE": None,
    "CUT_ERROR": "xx",
}


def status_bar_chart(
    df: pd.DataFrame,
    output: str | Path | None = None,
    label_map: dict[str, str] | None = None,
    baseline: str | None = None,
    include_settings: list[str] | None = None,
) -> plt.Figure:
    """Stacked bar chart showing instance status counts per setting.

    Each bar has the same total height (all instances), with segments for
    OPTIMAL, TIME_LIMIT, CUT_ERROR, INFEASIBLE, and DID_NOT_RUN stacked.
    """
    setup_latex_style()
    if label_map is None:
        label_map = {}

    settings = _order_settings(list(df["settings"].unique()), baseline, include_settings)
    statuses = sorted(
        df["status"].unique(),
        key=lambda s: list(_STATUS_COLORS.keys()).index(s)
        if s in _STATUS_COLORS else 999,
    )

    # Compute DID_NOT_RUN: instances present in at least one setting but missing from another
    has_tal = "tal" in df.columns
    group_cols = ["instance", "tal"] if has_tal else ["instance"]
    all_instances = df[group_cols].drop_duplicates()
    n_total = len(all_instances)
    did_not_run_counts = {}
    for s in settings:
        s_instances = df[df["settings"] == s][group_cols].drop_duplicates()
        did_not_run_counts[s] = n_total - len(s_instances)

    if any(v > 0 for v in did_not_run_counts.values()):
        if "DID_NOT_RUN" not in statuses:
            statuses.append("DID_NOT_RUN")

    counts = df.groupby(["settings", "status"]).size().unstack(fill_value=0).reindex(
        index=settings, columns=statuses, fill_value=0)
    if "DID_NOT_RUN" in statuses:
        for s in settings:
            counts.loc[s, "DID_NOT_RUN"] = did_not_run_counts.get(s, 0)

    n_settings = len(settings)
    x = np.arange(n_settings)
    width = 0.6

    fig, ax = plt.subplots(figsize=(max(4, n_settings * 1.2), 3.0))
    ax.set_axisbelow(True)

    bottom = np.zeros(n_settings)
    for status in statuses:
        vals = counts[status].values.astype(float)
        color = _STATUS_COLORS.get(status, "#CCCCCC")
        hatch = _STATUS_HATCHES.get(status)
        ax.bar(x, vals, width=width, bottom=bottom,
               label=status.replace("_", r"\_"),
               color=color, edgecolor="black", linewidth=0.5,
               hatch=hatch)
        # Value labels inside large segments
        for i, (v, b) in enumerate(zip(vals, bottom)):
            if v >= 10:
                cy = b + v / 2
                dark = color <= "#666666"
                fontcolor = "white" if dark else "black"
                text_kw = dict(ha="center", va="center", fontsize=7,
                               color=fontcolor)
                if hatch:
                    text_kw["bbox"] = dict(facecolor="white", alpha=0.7,
                                           edgecolor="none", pad=1.0)
                ax.text(x[i], cy, str(int(v)), **text_kw)
        bottom += vals

    # Total label on top of each bar
    for i in range(n_settings):
        ax.text(x[i], bottom[i] + 1, str(int(bottom[i])),
                ha="center", va="bottom", fontsize=7, color="gray")

    ax.set_xticks(x)
    ax.set_xticklabels([fmt_setting(s, label_map) for s in settings],
                       rotation=30, ha="right")
    ax.set_ylabel("Number of instances")
    ax.legend(frameon=False, fontsize=7, loc="center left",
             bbox_to_anchor=(1.01, 0.5))
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

    fig.tight_layout()
    if output:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output)
        if output.suffix == ".pdf":
            fig.savefig(output.with_suffix(".png"))
        print(f"Saved status bar chart to {output}")
    return fig


def time_boxplot(
    df: pd.DataFrame,
    output: str | Path | None = None,
    label_map: dict[str, str] | None = None,
    time_components: list[tuple[str, str]] | None = None,
    log_scale: bool = False,
    baseline: str | None = None,
    include_settings: list[str] | None = None,
) -> plt.Figure:
    """Grouped boxplot of per-instance times: one pair of boxes per setting.

    For each setting, two side-by-side boxes are drawn (e.g. Master and
    Subproblem time), allowing direct visual comparison of time spent in
    each phase.
    """
    setup_latex_style()
    if label_map is None:
        label_map = {}
    if time_components is None:
        time_components = DEFAULT_TIME_COMPONENTS

    settings = _order_settings(list(df["settings"].unique()), baseline, include_settings)
    col_names = [tc[0] for tc in time_components]
    col_labels = [tc[1] for tc in time_components]
    n_comp = len(col_names)

    # Collect data: list of arrays per (setting, component)
    all_data = []
    positions = []
    group_centers = []
    group_width = n_comp + 1  # spacing between setting groups
    box_width = 0.6

    for i, s in enumerate(settings):
        s_df = df[df["settings"] == s]
        center = i * group_width
        group_centers.append(center)
        for j, col in enumerate(col_names):
            vals = s_df[col].dropna().values
            all_data.append(vals)
            positions.append(center + (j - (n_comp - 1) / 2) * 0.8)

    fig, ax = plt.subplots(figsize=(max(4, len(settings) * 1.2), 3.0))
    ax.set_axisbelow(True)

    bp = ax.boxplot(
        all_data,
        positions=positions,
        widths=box_width,
        patch_artist=True,
        showfliers=True,
        flierprops=dict(marker=".", markersize=3, markerfacecolor="gray",
                        markeredgecolor="gray"),
        medianprops=dict(color="black", linewidth=1.2),
        whiskerprops=dict(color="black", linewidth=0.8),
        capprops=dict(color="black", linewidth=0.8),
    )

    # Alternate fill for the two components
    comp_fills = ["#404040", "#B0B0B0"]
    for idx, patch in enumerate(bp["boxes"]):
        comp_idx = idx % n_comp
        patch.set_facecolor(comp_fills[comp_idx % len(comp_fills)])
        patch.set_edgecolor("black")
        patch.set_linewidth(0.8)

    # Sample-size annotation per setting
    for i, s in enumerate(settings):
        n = len(df[df["settings"] == s])
        ax.text(group_centers[i], 1.02, f"$n\\!=\\!{n}$",
                transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=6, color="gray")

    ax.set_xticks(group_centers)
    ax.set_xticklabels([fmt_setting(s, label_map) for s in settings],
                       rotation=30, ha="right")
    ax.set_ylabel(r"Time (s)")

    # Same fragility as cuts_boxplot: log scale with all-zero data lets
    # matplotlib stretch the canvas to absurd dimensions during savefig.
    has_positive = any((d > 0).any() for d in all_data if len(d) > 0)
    if log_scale and has_positive:
        ax.set_yscale("log")

    # Legend for components
    from matplotlib.patches import Patch
    legend_handles = [Patch(facecolor=comp_fills[j % len(comp_fills)],
                            edgecolor="black", linewidth=0.8, label=col_labels[j])
                      for j in range(n_comp)]
    ax.legend(handles=legend_handles, frameon=False, loc="upper left")

    fig.tight_layout()
    if output:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output)
        if output.suffix == ".pdf":
            fig.savefig(output.with_suffix(".png"))
        elif output.suffix == ".png":
            fig.savefig(output.with_suffix(".pdf"))
        print(f"Saved time boxplot to {output}")
    return fig


def cuts_boxplot(
    df: pd.DataFrame,
    column: str = "n_benders_cuts_added",
    output: str | Path | None = None,
    label_map: dict[str, str] | None = None,
    log_scale: bool = True,
    baseline: str | None = None,
    include_settings: list[str] | None = None,
) -> plt.Figure:
    """Boxplot of per-instance cut counts, one box per setting."""
    setup_latex_style()
    if label_map is None:
        label_map = {}

    settings = _order_settings(list(df["settings"].unique()), baseline, include_settings)
    data = [df.loc[df["settings"] == s, column].dropna().values for s in settings]

    fig, ax = plt.subplots(figsize=(max(4, len(settings) * 1.2), 4))
    ax.set_axisbelow(True)

    bp = ax.boxplot(
        data,
        labels=[fmt_setting(s, label_map) for s in settings],
        patch_artist=True,
        showfliers=True,
        flierprops=dict(marker=".", markersize=3, markerfacecolor="gray",
                        markeredgecolor="gray"),
        medianprops=dict(color="black", linewidth=1.2),
        whiskerprops=dict(color="black", linewidth=0.8),
        capprops=dict(color="black", linewidth=0.8),
    )

    # Greyscale fill for boxes
    greys = ["#D0D0D0", "#A0A0A0", "#808080", "#606060", "#B0B0B0",
             "#909090", "#C0C0C0", "#707070"]
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(greys[i % len(greys)])
        patch.set_edgecolor("black")
        patch.set_linewidth(0.8)

    # Annotate median values and sample sizes
    # Guard against placing text at y=0 on a log axis (causes degenerate
    # figure dimensions on savefig).
    has_positive = any((d > 0).any() for d in data if len(d) > 0)
    use_log = log_scale and has_positive
    for i, d in enumerate(data):
        if len(d) > 0:
            med = float(np.median(d))
            if use_log and med <= 0:
                med = max(float(d[d > 0].min()) if (d > 0).any() else 1.0, 1e-3)
            ax.text(i + 1, med, f" {med:,.0f}", va="center", ha="left",
                    fontsize=6, color="black")
        ax.text(i + 1, 1.02, f"$n\!=\!{len(d)}$", transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=6, color="gray")

    ax.set_xticklabels([fmt_setting(s, label_map) for s in settings],
                       rotation=30, ha="right")
    ax.set_ylabel(r"Number of Benders cuts added")

    if use_log:
        ax.set_yscale("log")

    fig.tight_layout()
    if output:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output)
        if output.suffix == ".pdf":
            fig.savefig(output.with_suffix(".png"))
        print(f"Saved cuts boxplot to {output}")
    return fig
