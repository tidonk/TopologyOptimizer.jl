"""Plot calculations and styles used for the published figures."""
import numpy as np
import pandas as pd
HORIZON = 21600
MAX_FACTOR = 100
METHODS = ("direct", "sequential", "full", "warmstart", "parallel")
SEQUENTIAL_COLOR = "#DC267F"

STYLE = r'''grid=major,major grid style={draw=black!12,line width=0.3pt},
xmajorgrids=false,axis line style={draw=black!65,line width=0.4pt},
tick align=outside,tick style={draw=black!65},
label style={font=\small},tick label style={font=\footnotesize},
title style={font=\small},clip mode=individual,
'''.strip()


def coordinates(x, y):
    return ' '.join(f'({a:.10g},{b:.10g})' for a, b in zip(x, y))


def render_time_selection(curves, selection, save):
    """Render matched time curves, leaving panels without results empty."""
    lines = [
        "% Clock: solver time excluding each run's SCOPF phase.",
        "% Events at or before 1 s are included at the left edge of the logarithmic axis.",
        "% Horizontal tails retain attained solutions after a run stops.",
        "% Common instances and best objective across all five configurations.",
        "% Parallel heuristic uses the receiver clock; concurrent helper time is not added.",
        r"\definecolor{studyDirect}{HTML}{FE6100}",
        r"\definecolor{studyHeuristic}{HTML}{FFB000}",
        r"\definecolor{studyExact}{HTML}{648FFF}",
        r"\definecolor{studyWarmstart}{HTML}{DC267F}",
        r"\definecolor{studyParallel}{HTML}{785EF0}",
        r"\pgfplotslegendfromname{heuristicSelectionLegend}\par\smallskip",
        r"\noindent\hspace*{0.08\linewidth}%",
    ]
    for z in range(1, 6):
        lines.extend([
            r"\begin{subfigure}[t]{0.18\linewidth}",
            r"\centering",
            r"\captionsetup{justification=centering}",
            r"\begin{tikzpicture}[baseline=(current axis.north),trim axis left,trim axis right]",
            r"\begin{axis}[scale only axis,width=0.16\practicalFigureWidth,height=0.18\practicalFigureWidth,",
            r"xmode=log,xmin=1,xmax=21600,ymin=0,ymax=105,ytick={0,25,50,75,100},",
            r"xtick={1,10,100,1000,10000},xticklabels={$10^0$,$10^1$,$10^2$,$10^3$,$10^4$},xminorticks=false,",
            r"xlabel={Solver time (s)},",
            STYLE,
            r"xmajorgrids=true,",
            r"xticklabel style={font=\scriptsize,inner xsep=0pt},",
        ])
        if z == 1:
            lines.extend([
                r"ylabel={Instances (\%)},legend to name=heuristicSelectionLegend,",
                r"legend columns=5,legend cell align={left},",
                r"legend style={draw=none,fill=none,inner sep=0pt,font=\footnotesize,",
                r"/tikz/every even column/.append style={column sep=0.6em}},",
            ])
        else:
            lines.append(r"yticklabels={},")
        lines.append("]")
        n = next(r["n"] for r in selection if r["z"] == z)
        if not n:
            lines.append(r"\node[font=\scriptsize,align=center,text=gray] at (rel axis cs:0.5,0.5) {Results\\pending};")
        for method, color, style, label in [
                ("direct", "studyDirect", "dashed", "Single-stage"),
                ("sequential", "studyHeuristic", "solid", "Sequential"),
                ("full", "studyExact", "solid", "Exact"),
                ("warmstart", "studyWarmstart", "solid", "Heuristic warmstart"),
                ("parallel", "studyParallel", "solid", "Parallel heuristic")]:
            if not n or (z, method) not in curves:
                continue
            coords = " ".join(f"({t:.12g},{p:.12g})" for t, p in curves[z, method])
            lines.append(rf"\addplot[const plot,no marks,color={color},{style},line width=0.9pt] coordinates {{{coords}}};")
            if z == 1:
                lines.append(rf"\addlegendentry{{{label}}}")
        lines.extend([
            r"\end{axis}",
            r"\end{tikzpicture}",
            rf"\caption{{$Z={z}$ ($n={next(r['n'] for r in selection if r['z'] == z)}$).}}",
            rf"\label{{fig:heuristic-selection-z{z}}}",
            r"\end{subfigure}" + (r"\hfill" if z < 5 else ""),
        ])
    save("heuristic-selection-within-one-percent.tikz.tex", lines)


def tikz_figure(summaries, n):
    maximum = len(summaries['mean']['achieved']) - 1
    limits = range(maximum + 1)
    def coords(x, y):
        return ' '.join(f'({a:.10g},{b:.10g})' for a, b in zip(x, y))

    upper = max(5, 5*np.ceil(max(max(row['upper_bound']) for row in summaries.values())/5))
    # Use the same axis dimensions as one panel of the attainment figure.
    lines = [rf'\def\cumulativeSavingsN{{{n}}}',
             r'\ifdefined\practicalFigureWidth\else',
             r'\newlength{\practicalFigureWidth}\setlength{\practicalFigureWidth}{\linewidth}\fi',
             r'\begin{tikzpicture}[baseline=(current axis.north)]',
             rf'\definecolor{{practicalSequential}}{{HTML}}{{{SEQUENTIAL_COLOR[1:]}}}',
             r'\begin{axis}[scale only axis,width=0.255\practicalFigureWidth,height=0.18\practicalFigureWidth,',
             rf'xmin=0,xmax={maximum},ymin=0,ymax={upper:g},xtick={{{",".join(map(str, limits))}}},ytick distance=5,',
             r'xlabel={Switching limit $Z$},ylabel={Objective reduction\\relative to SCOPF (\%)},',
             r'ylabel style={align=center},',
             r'grid=major,major grid style={draw=black!12,line width=0.3pt},xmajorgrids=false,',
             r'axis line style={draw=black!65,line width=0.4pt},tick align=outside,',
             r'tick style={draw=black!65},label style={font=\small},tick label style={font=\footnotesize},',
             r'legend columns=2,legend cell align={left},',
             r'legend style={at={(0.5,1.08)},anchor=south,draw=none,fill=none,font=\footnotesize,',
             r'/tikz/every even column/.append style={column sep=0.9em}},]']
    for key, shade in [('mean', 14), ('median', 7)]:
        row = summaries[key]
        polygon = coords(limits, row['achieved']) + ' ' + coords(range(maximum, -1, -1), row['upper_bound'][::-1])
        lines.append(rf'\addplot[draw=none,fill=practicalSequential!{shade},forget plot] coordinates {{{polygon}}} \closedcycle;')
    for key, label, mark, options in [('mean', 'Mean', '*', 'fill=practicalSequential'),
                                       ('median', 'Median', 'square*', 'fill=white')]:
        row = summaries[key]
        lines.append(rf'\addplot[color=practicalSequential,dashed,no marks,line width=0.7pt,forget plot] coordinates {{{coords(limits, row["upper_bound"])}}};')
        lines.append(rf'\addplot[color=practicalSequential,solid,mark={mark},mark options={{{options}}},line width=0.9pt,mark size=1.8pt] coordinates {{{coords(limits, row["achieved"])}}};')
        lines.append(rf'\addlegendentry{{{label}}}')
    lines.extend([r'\end{axis}', r'\end{tikzpicture}'])
    return '\n'.join(lines)+'\n'


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
