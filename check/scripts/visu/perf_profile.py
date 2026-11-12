"""Performance profile plots (Dolan-Moré and cumulative-over-time).

Default: Dolan-Moré — fraction of instances solved vs. slowdown factor (log-x).
Alternative: cumulative percentage of instances solved vs. absolute time.

Publication-ready: LaTeX fonts, no figure title, clean axis labels.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from reader import read_experiment
from style import COLORS, setup_latex_style, fmt_setting

TIME_LIMIT = 10800
MAX_FACTOR = 100  # cap for unsolved/failed


def _compute_slowdown(df: pd.DataFrame, target: str, settings: list[str]) -> pd.DataFrame:
    """Add slowdown factor column: value / best across settings per instance(+tal)."""
    group_cols = ["instance", "tal"] if "tal" in df.columns else ["instance"]

    # Best per instance across settings (only OPTIMAL)
    solved = df[df["status"].isin(["OPTIMAL"])].copy()
    best = solved.groupby(group_cols)[target].min().rename("_best")
    df = df.merge(best, on=group_cols, how="left")

    df["slowdown"] = np.where(
        df["status"].isin(["OPTIMAL"]) & df["_best"].notna(),
        df[target] / df["_best"],
        MAX_FACTOR + 1,  # beyond x_grid range so unsolved lines stay < 1.0
    )
    df.drop(columns="_best", inplace=True)
    return df


def _ecdf(values: np.ndarray, x_grid: np.ndarray, n_total: int | None = None) -> np.ndarray:
    """Compute empirical CDF of *values* evaluated at *x_grid*.

    Parameters
    ----------
    n_total : if given, use as denominator instead of len(values).
              This accounts for instances missing from this setting (DID_NOT_RUN).
    """
    sorted_vals = np.sort(values)
    denom = n_total if n_total is not None else len(sorted_vals)
    return np.searchsorted(sorted_vals, x_grid, side="right") / denom


def performance_profile(
    df: pd.DataFrame,
    target: str = "time_solve",
    mode: str = "factor",
    output: str | Path | None = None,
    label_map: dict[str, str] | None = None,
    legend_title: str = "Setting",
    style_map: dict[str, dict] | None = None,
    baseline: str | None = None,
    include_settings: list[str] | None = None,
) -> plt.Figure:
    """Plot a performance profile.

    Parameters
    ----------
    df : experiment DataFrame (from reader.read_experiment).
    target : metric column.
    mode : 'factor' (Dolan-Moré, default) or 'time' (cumulative over absolute time).
    output : if given, save to this path (extension determines format).
    label_map : optional mapping from setting name → display label.
    legend_title : title for the legend box.
    style_map : optional mapping from setting name → dict with 'color', 'linestyle' keys.
    include_settings : if given, use this order for settings in the plot.
    """
    setup_latex_style()

    all_settings = list(df["settings"].unique())
    if include_settings is not None:
        # Hard filter: plot ONLY the listed settings, in the listed order.
        # This is what the top-N selection in run_full_comparison relies on.
        present = set(all_settings)
        settings = [s for s in include_settings if s in present]
    else:
        rest = sorted(s for s in all_settings if s != baseline)
        settings = ([baseline] + rest) if baseline and baseline in all_settings else rest
    if label_map is None:
        label_map = {}

    def _line_style(i: int, s: str) -> dict:
        if style_map and s in style_map:
            return {k: v for k, v in style_map[s].items() if k in ("color", "linestyle")}
        return {"color": COLORS[i % len(COLORS)]}

    group_cols = ["instance", "tal"] if "tal" in df.columns else ["instance"]
    n_total = len(df.groupby(group_cols))

    # Bigger canvas; the legend reservation below trims it down to a usable
    # plotting width even when there are many settings.
    fig, ax = plt.subplots(figsize=(9, 5.5))

    if mode == "factor":
        df = _compute_slowdown(df, target, settings)
        x_grid = np.logspace(0, np.log10(MAX_FACTOR), 500)

        for i, s in enumerate(settings):
            vals = df[df["settings"] == s]["slowdown"].values
            ax.plot(x_grid, _ecdf(vals, x_grid, n_total=n_total),
                    label=fmt_setting(s, label_map), linewidth=1.5, **_line_style(i, s))

        ax.set_xscale("log")
        ax.set_xlabel(r"Performance ratio $\tau$")
        ax.set_ylabel(r"Fraction of instances $P(r \leq \tau)$")

    elif mode == "time":
        x_max = TIME_LIMIT
        x_grid = np.linspace(0, x_max, 1000)

        for i, s in enumerate(settings):
            sub = df[df["settings"] == s]
            # Unsolved get TIME_LIMIT + 1 (beyond grid) so line stays < 1.0
            times = np.where(
                sub["status"].isin(["OPTIMAL"]),
                sub[target].values,
                TIME_LIMIT + 1,
            )
            ax.plot(x_grid, _ecdf(times, x_grid, n_total=n_total),
                    label=fmt_setting(s, label_map), linewidth=1.5, **_line_style(i, s))

        ax.set_xlabel(r"Time (s)")
        ax.set_ylabel(r"Fraction of instances solved")
        ax.set_xticks(range(0, int(x_max) + 1, 1200))

    else:
        raise ValueError(f"Unknown mode '{mode}', use 'factor' or 'time'.")

    n_pairs = n_total
    ax.set_ylim(-0.02, 1.02)
    # Legend placement: side legend for short lists (preserves the original
    # paper-style layout); bottom multi-column legend when the side legend
    # would push the axes too narrow.
    n_entries = len(settings)
    if n_entries <= 8:
        ax.legend(title=legend_title, frameon=False,
                  loc="center left", bbox_to_anchor=(1.02, 0.5))
        layout_rect = [0, 0, 0.82, 1]
    else:
        # 4 entries per column keeps the legend block roughly 2 lines tall
        # for up to 16 entries, 3 lines for up to 24 (we never go higher).
        ncol = min(4, (n_entries + 3) // 4)
        # Reserve space below the axes proportional to the row count.
        n_rows_leg = -(-n_entries // ncol)  # ceil
        bottom = 0.02 + 0.045 * n_rows_leg
        ax.legend(title=legend_title, frameon=False,
                  loc="upper center",
                  bbox_to_anchor=(0.5, -0.13),
                  ncol=ncol, columnspacing=1.4, handletextpad=0.5)
        layout_rect = [0, bottom, 1, 1]
    ax.text(0.01, 0.99, f"$n = {n_pairs}$",
            transform=ax.transAxes, va="top", ha="left", color="gray")

    fig.tight_layout(rect=layout_rect)
    if output:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output)
        # Also save sibling format
        if output.suffix == ".pdf":
            fig.savefig(output.with_suffix(".png"))
        elif output.suffix == ".png":
            fig.savefig(output.with_suffix(".pdf"))
        print(f"Saved performance profile to {output}")
    return fig


def plot_from_experiment(
    experiment_dir: str | Path,
    target: str = "time_solve",
    mode: str = "factor",
    output: str | Path | None = None,
    label_map: dict[str, str] | None = None,
    legend_title: str = "Setting",
    style_map: dict[str, dict] | None = None,
) -> plt.Figure:
    """Convenience: read experiment and plot."""
    df = read_experiment(experiment_dir)
    if output is None:
        experiment_dir = Path(experiment_dir)
        output = experiment_dir / f"performance_profile_{target}.pdf"
    return performance_profile(df, target=target, mode=mode, output=output,
                               label_map=label_map, legend_title=legend_title,
                               style_map=style_map)


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python perf_profile.py <experiment_dir> [target] [mode]")
        sys.exit(1)
    exp = sys.argv[1]
    target = sys.argv[2] if len(sys.argv) > 2 else "time_solve"
    mode = sys.argv[3] if len(sys.argv) > 3 else "factor"
    plot_from_experiment(exp, target=target, mode=mode)
