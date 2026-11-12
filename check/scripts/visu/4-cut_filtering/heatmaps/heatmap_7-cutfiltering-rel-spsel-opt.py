#!/usr/bin/env python3
"""Heatmap for Exp. 4 Stage 2 - 7-cutfiltering-rel-spsel-opt.

Two output files (each PDF + PNG), each a 2x2 grid:
  rows = STOP_AFTER_N_HITS_REL fixed at 0.05 (top) and 0.10 (bottom)
  cols = absolute geomean of t_solve (left) / relative-to-baseline % (right)

  heatmap_time_solve.{pdf,png}     - common-solved subset
  heatmap_time_solve_all.{pdf,png} - all instances, non-OPTIMAL -> TIME_LIMIT

Y axis: "Number of cuts per callback" with ticks 2, 5, 10, inf (cf=100% / sp=100
are the no-CF-filtering / no-SP-selection baselines).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _heatmap_common import render_heatmap_stacked  # noqa: E402

# SP=100 row (CF only, no SP selection) - shared by both rel-split groups.
SP100_BASELINES = [
    "cf001_sp100", "cf005_sp100", "cf010_sp100", "cf100_sp100",
]

# rel=0.05 grid (cf x sp). cf010_sp02_005 / cf010_sp05_005 are queued reruns
# (originals were mislabeled rel=0.10 and renamed to cf010_sp02_010 /
# cf010_sp05_010_dup; new correct rel=0.05 tomls await cluster execution).
INCLUDE_REL005 = [
    "cf001_sp02_005", "cf005_sp02_005", "cf010_sp02_005", "cf100_sp02_005",
    "cf001_sp05_005", "cf005_sp05_005", "cf010_sp05_005", "cf100_sp05_005",
    "cf001_sp10_005", "cf005_sp10_005", "cf010_sp10_005", "cf100_sp10_005",
    *SP100_BASELINES,
]

# rel=0.10 grid (cf x sp). cf100_sp{02,05,10}_010 are symlinks to the _005
# runs (cf=100% means no CF filtering and the abs N_HITS dominates over REL on
# the typical instance sizes; see check/results/.../cf100_sp*_010).
INCLUDE_REL010 = [
    "cf001_sp02_010", "cf005_sp02_010", "cf010_sp02_010", "cf100_sp02_010",
    "cf001_sp05_010", "cf005_sp05_010", "cf010_sp05_010", "cf100_sp05_010",
    "cf001_sp10_010", "cf005_sp10_010", "cf010_sp10_010", "cf100_sp10_010",
    *SP100_BASELINES,
]

GROUPS = [
    (r"REL$=0.05$", INCLUDE_REL005),
    (r"REL$=0.10$", INCLUDE_REL010),
]

PATTERN = r"cf(\d+)_sp(\d+)(?:_\d+)?"

# Common-solved subset.
render_heatmap_stacked(
    "7-cutfiltering-rel-spsel-opt",
    groups=GROUPS,
    baseline="cf100_sp100",
    output_name="heatmap_time_solve.pdf",
    pattern=PATTERN,
    common_only=True,
    suptitle=r"Exp.~7: CF $\times$ \#cuts grid, common-solved subset",
)

# All-instance, TL-penalised view.
render_heatmap_stacked(
    "7-cutfiltering-rel-spsel-opt",
    groups=GROUPS,
    baseline="cf100_sp100",
    output_name="heatmap_time_solve_all.pdf",
    pattern=PATTERN,
    common_only=False,
    penalize_nonoptimal=True,
    min_time=0.0,
    suptitle=r"Exp.~7: CF $\times$ \#cuts grid, all instances (TL penalised)",
)
