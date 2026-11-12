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

### Figures

import re

import matplotlib.backends.backend_pgf as backend_pgf
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

# tikzplotlib 0.10 expects names that newer numpy / matplotlib removed.
if not hasattr(np, "float_"):
    np.float_ = np.float64
if not hasattr(backend_pgf, "common_texification"):
    backend_pgf.common_texification = backend_pgf._tex_escape
try:
    import tikzplotlib
except ImportError as exc:
    raise ImportError(
        f"tikzplotlib is required for TikZ exports ({exc}), install check/scripts/visu/requirements.txt"
    ) from exc

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.grid": True,
    "grid.alpha": 0.3,
})

COLORS = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9"]
MAX_WIDTH_IN = 7 / 2.54 
PROFILE_WIDTH_IN = MAX_WIDTH_IN + 0.45
MAX_FACTOR = 100  # performance-ratio axis
AXIS_FONTSIZE = "small"

# Detailed time decomposition (column, label)
TIME_COMPONENTS = [
    ("_sp_primal", "Primal"),
    ("_sp_dual", "Dual"),
    ("time_SP_N1_refactorization", "Update"),
    ("_sp_remaining", "Misc. SP"),
    ("time_MP", "MP"),
]


def box_height_cm(n_settings: int) -> float:
    return 7 * 5.1 / max(5.0, n_settings * 1.35)


def _save(fig, path, boxplot: bool = False) -> None:
    """Save fig as PDF and TikZ file"""
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches=None)

    for ax in fig.axes:
        if (leg := ax.get_legend()) is not None:
            leg._ncol, leg.legendHandles = leg._ncols, leg.legend_handles
    for line in fig.findobj(match=Line2D):
        line._us_dashOffset, line._us_dashSeq = line._dash_pattern

    w, h = fig.get_size_inches()
    tikz = path.with_suffix(".tikz")
    tikzplotlib.save(str(tikz), figure=fig, axis_width=f"{w * 2.54:.5f}cm", axis_height=f"{h * 2.54:.5f}cm")
    t = tikz.read_text(encoding="utf-8", errors="ignore")
    t = t.replace("legend style={", "legend style={\n  font=\\footnotesize,\n  fill=none,")
    t = re.sub(r"fill opacity=\s*[-+0-9.]+", "fill opacity=1", t)
    t = t.replace("\\begin{tikzpicture}", "\\begin{tikzpicture}\n\\providecommand{\\mathdefault}[1]{#1}", 1)
    if boxplot:
        # add background rectangles and legend entries
        ticks = [float(v) for v in re.search(r"xtick=\{([^}]*)\}", t).group(1).split(",") if v.strip()]
        y0 = float(re.search(r"ymin=([-+0-9.eE]+)", t).group(1))
        y1 = float(re.search(r"ymax=([-+0-9.eE]+)", t).group(1))
        trim = (y1 - y0) * 2.0 / max(h * 300.0, 1.0)
        half = 0.5 * (ticks[1] - ticks[0])
        strips = "".join(
            f"\\path [draw=none, fill=whitesmoke242] (axis cs:{c - half},{y0 + trim}) "
            f"rectangle (axis cs:{c + half},{y1 - trim});\n"
            for c in ticks[1::2]
        )
        header = re.search(r"(\\begin\{axis\}\[[\s\S]*?\]\n)", t).group(1)
        t = t.replace(header, header + "\n" + strips, 1)
        leg = fig.axes[0].get_legend()
        entries = ""
        for text, handle in zip(leg.get_texts(), leg.legend_handles):
            r, g, b = (round(c * 255) for c in mcolors.to_rgb(handle.get_markerfacecolor()))
            entries += (f"\\addlegendimage{{only marks, mark=square*, mark size=2.5pt, draw=black, "
                        f"fill={{rgb,255:red,{r};green,{g};blue,{b}}}}}\n\\addlegendentry{{{text.get_text()}}}\n")
        t = t.replace("\\end{axis}", "\n" + entries + "\\end{axis}", 1)
    else:
        t = re.sub(r"x grid style=\{[^}]+\}", "x grid style={gray,opacity=0.3}", t)
        t = re.sub(r"y grid style=\{[^}]+\}", "y grid style={gray,opacity=0.3}", t)
        t = t.replace("[semithick,", "[line width=1.5pt,")
    tikz.write_text(t)
    plt.close(fig)
    print(f"Saved {path}")


def performance_profile(df, target, path, settings, labels, styles, height_cm, legend_left=False) -> None:
    keys = ["instance", "tal"] if "tal" in df.columns else ["instance"]
    n_total = df.groupby(keys).ngroups
    best = df[df["status"] == "OPTIMAL"].groupby(keys)[target].min().rename("_best")
    df = df.merge(best, on=keys, how="left")
    slowdown = np.where((df["status"] == "OPTIMAL") & df["_best"].notna(), df[target] / df["_best"], MAX_FACTOR + 1)
    x = np.logspace(0, np.log10(MAX_FACTOR), 500)

    fig, ax = plt.subplots(figsize=(PROFILE_WIDTH_IN, height_cm / 2.54))
    for s in settings:
        y = np.searchsorted(np.sort(slowdown[df["settings"] == s]), x, side="right") / n_total
        ax.plot(x, y, label=labels[s], linewidth=1.5, color=styles[s]["color"], linestyle=styles[s]["linestyle"])
    ax.set_xscale("log")
    ax.set_xlabel(r"Performance ratio $\tau$", fontsize=AXIS_FONTSIZE)
    ax.set_ylabel(r"Fraction of instances $P(r \leq \tau)$", fontsize=AXIS_FONTSIZE)
    ax.set_ylim(-0.02, 1.02)
    ax.set_box_aspect(1.0)
    if legend_left:  # move legend if overlapping
        ax.legend(title="Setting", frameon=False, loc="lower left", bbox_to_anchor=(0.0, 0.42))
    else:
        ax.legend(title="Setting", frameon=False, loc="lower right")
    fig.tight_layout()
    _save(fig, path)


def runtime_scatter(df, target, path, sort_setting, settings, labels, styles, height_cm) -> None:
    keys = ["instance", "tal"] if "tal" in df.columns else ["instance"]
    pivot = df.pivot_table(index=keys, columns="settings", values=target, aggfunc="first")
    order = pd.to_numeric(pivot[sort_setting], errors="coerce").dropna().sort_values().index
    pivot = pivot.reindex(columns=settings).loc[order]
    n = len(pivot)
    step = min(0.8 / len(settings), 0.1)

    fig, ax = plt.subplots(figsize=(PROFILE_WIDTH_IN, height_cm / 2.54))
    ys = []
    for i, s in enumerate(settings):
        y = pd.to_numeric(pivot[s], errors="coerce").values
        ok = np.isfinite(y) & (y > 0)
        x = np.arange(n) + (i - (len(settings) - 1) / 2.0) * step
        ax.scatter(x[ok], y[ok], s=10.08, alpha=0.9, linewidths=0, label=labels[s],
                   color=styles[s]["color"], marker=styles[s].get("marker"))
        ys.append(y[ok])
    ys = np.concatenate(ys)
    ax.set_xlim(-4.0, n - 1 + 4.0)
    ax.set_ylim(float(ys.min()) / 4.05, float(ys.max()) * 4.05)
    ax.set_xticks([])
    ax.set_xlabel("Instance (sorted by running time)", fontsize=AXIS_FONTSIZE)
    ax.set_ylabel(r"Runtime (s)", fontsize=AXIS_FONTSIZE)
    ax.set_yscale("log")
    ax.set_box_aspect(1.0)
    ax.legend(title="Setting", frameon=False, loc="lower right")
    fig.tight_layout()
    _save(fig, path)


def time_boxplot(df, path, settings, labels, relative=False) -> None:
    cols = [c for c, _ in TIME_COMPONENTS]
    spacing = 0.75
    group_width = max(1.8, (len(cols) - 1) * spacing + 1.2)
    data, positions, centers = [], [], []
    for i, s in enumerate(settings):
        centers.append(i * group_width)
        comp = df[df["settings"] == s][cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        total = comp.sum(axis=1).replace(0, np.nan)
        for j, c in enumerate(cols):
            data.append((comp[c] / total * 100.0).dropna().values if relative else comp[c].dropna().values / 3600.0)
            positions.append(centers[-1] + (j - (len(cols) - 1) / 2) * spacing)

    fig, ax = plt.subplots(figsize=(max(5, len(settings) * 1.35), 5.1))
    ax.set_axisbelow(True)
    ax.set_facecolor("white")
    ax.grid(False)
    for c in centers[1::2]:
        ax.add_patch(Rectangle((c - group_width / 2.0, 0.0), group_width, 1.0, facecolor="#F2F2F2", edgecolor="none",
                               zorder=0, transform=ax.get_xaxis_transform(), clip_on=True))
    bp = ax.boxplot(
        data, positions=positions, widths=min(0.55, spacing * 0.75), patch_artist=True, showfliers=True,
        flierprops=dict(marker=".", markersize=3, markerfacecolor="gray", markeredgecolor="gray"),
        medianprops=dict(color="black", linewidth=1.2),
        whiskerprops=dict(color="black", linewidth=0.8),
        capprops=dict(color="black", linewidth=0.8),
    )
    for k, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(COLORS[k % len(cols)])
        patch.set_edgecolor("black")
        patch.set_linewidth(0.8)
    ax.set_xticks(centers)
    ax.set_xticklabels([labels[s] for s in settings], rotation=0, ha="center")
    ax.set_ylabel(r"Share of solve time (\%)" if relative else r"Time (h)")
    upper = 100.0 if relative else 3.0
    gap = upper * 0.03
    ax.set_ylim(-gap, upper + gap)
    handles = [Line2D([0], [0], marker="s", linestyle="None", markerfacecolor=COLORS[j], markeredgecolor="black",
                      markeredgewidth=0.8, markersize=6, label=label) for j, (_, label) in enumerate(TIME_COMPONENTS)]
    x_mid = 0.5 * (centers[0] + centers[1]) if len(centers) > 1 else centers[0]
    x_anchor = ax.transAxes.inverted().transform(ax.transData.transform((x_mid, 0.0)))[0]
    ax.legend(handles=handles, frameon=False, loc="upper center", bbox_to_anchor=(x_anchor, 0.995), fontsize=7)
    fig.tight_layout()
    w, h = fig.get_size_inches()
    fig.set_size_inches(MAX_WIDTH_IN, h * (MAX_WIDTH_IN / w), forward=True)
    _save(fig, path, boxplot=True)
