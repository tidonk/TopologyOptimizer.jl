#!/usr/bin/env python3
"""Heatmaps for stage 5 partial-decomposition num grids (method × N).

Table-driven so the cells **match the tables exactly**: each heatmap reads the
per-setting ``Geo. Mean`` (on ``time_with_warmstart``) out of the leaf's own
``comparison_short_all.md`` and lays the scorers (rows) against the subproblem
count N (cols).  The off-grid ``No PD`` (cf05_sp005) row is the relative-panel
baseline.  A setting the table dropped for cut errors is a blank cell.

Two leaves, one heatmap each (into the same paper leaf as the table):
  num          -> paper/5-partial_decomp/grid/num/heatmap_time_solve_all.pdf
  num_no_base  -> paper/5-partial_decomp/grid/num_no_base/heatmap_time_solve_all.pdf

Run *after* ``5-partial_decomp/compare_grid.py``.
"""
import re
import sys
from pathlib import Path

VISU_DIR = Path(__file__).resolve().parents[2]   # check/scripts/visu/
sys.path.insert(0, str(VISU_DIR))
from heatmap import grid_heatmap_from_values, geomeans_from_short_table   # noqa: E402

PAPER_GRID = VISU_DIR.parents[2] / "check/results/paper/5-partial_decomp/grid"
METHODS = ["Viol", "Rnd", "CProx", "Div-Cos", "Div-Euc", "Sens"]   # row order
NUMS = [1, 2, 3, 4, 5]                                             # col order (N)


def _coord(label: str):
    """`Viol3` -> (N=3, method='Viol'); `No PD` / `Base only` -> None (off-grid)."""
    m = re.fullmatch(r"(.+?)(\d+)", label)
    if m and m.group(1) in METHODS:
        return int(m.group(2)), m.group(1)
    return None


for leaf in ("num", "num_no_base"):
    md = PAPER_GRID / leaf / "comparison_short_all.md"
    if not md.is_file():
        print(f"warning: {md} missing (run compare_grid.py first); skipping")
        continue
    gm, nn, note = geomeans_from_short_table(md, _coord)
    if not gm:
        print(f"warning: no method×N cells parsed from {leaf}; skipping")
        continue
    # off-grid baseline: the "No PD" (cf05_sp005) geomean, parsed straight from the table.
    base_gm, _, _ = geomeans_from_short_table(md, lambda l: (0, "b") if l == "No PD" else None)
    baseline_value = base_gm.get((0, "b"))
    grid_heatmap_from_values(
        gm, nn, PAPER_GRID / leaf / "heatmap_time_solve_all.pdf",
        x_levels=NUMS, y_levels=METHODS, baseline_value=baseline_value,
        xlabel=r"Subproblems left in master $N$", ylabel="Scoring method",
        title_abs=r"Geom.\ mean of $t_{\mathrm{solve}}+t_{\mathrm{ws}}$ [s]",
        title_rel=r"Relative to \texttt{No PD} (\%)", n_note=note,
    )
    print(f"Wrote {PAPER_GRID / leaf / 'heatmap_time_solve_all.pdf'} ({len(gm)} cells from table)")
