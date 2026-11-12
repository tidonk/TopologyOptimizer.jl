#!/usr/bin/env python3
"""Heatmap for Exp. 4 Stage 3 - 8-cutfiltering-hyb-spsel-opt.

Mirror of the exp-7 layout: two 2x2 grids (rows = REL=0.05 / REL=0.10;
cols = absolute geomean / relative-to-baseline %). cf100 columns and the sp=100
row are sourced via symlinks (see check/results/8-cutfiltering-hyb-spsel-opt).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _heatmap_common import render_heatmap_stacked  # noqa: E402

SP100_BASELINES = [
    "cf001_sp100", "cf005_sp100", "cf010_sp100", "cf100_sp100",
]

# rel=0.05 grid (9 hybrid tomls + cf100 symlinks).
INCLUDE_REL005 = [
    "cf001_sp02_005", "cf005_sp02_005", "cf010_sp02_005", "cf100_sp02_005",
    "cf001_sp05_005", "cf005_sp05_005", "cf010_sp05_005", "cf100_sp05_005",
    "cf001_sp10_005", "cf005_sp10_005", "cf010_sp10_005", "cf100_sp10_005",
    *SP100_BASELINES,
]

# rel=0.10 grid (6 hybrid tomls + cf100 symlinks; no sp02_010 hybrid tomls).
INCLUDE_REL010 = [
    "cf100_sp02_010",
    "cf001_sp05_010", "cf005_sp05_010", "cf010_sp05_010", "cf100_sp05_010",
    "cf001_sp10_010", "cf005_sp10_010", "cf010_sp10_010", "cf100_sp10_010",
    *SP100_BASELINES,
]

GROUPS = [
    (r"REL$=0.05$", INCLUDE_REL005),
    (r"REL$=0.10$", INCLUDE_REL010),
]

PATTERN = r"cf(\d+)_sp(\d+)(?:_\d+)?"

render_heatmap_stacked(
    "8-cutfiltering-hyb-spsel-opt",
    groups=GROUPS,
    baseline="cf100_sp100",
    output_name="heatmap_time_solve.pdf",
    pattern=PATTERN,
    common_only=True,
    suptitle=r"Exp.~8: CF $\times$ \#cuts grid (hybrid), common-solved subset",
)

render_heatmap_stacked(
    "8-cutfiltering-hyb-spsel-opt",
    groups=GROUPS,
    baseline="cf100_sp100",
    output_name="heatmap_time_solve_all.pdf",
    pattern=PATTERN,
    common_only=False,
    penalize_nonoptimal=True,
    min_time=0.0,
    suptitle=r"Exp.~8: CF $\times$ \#cuts grid (hybrid), all instances (TL penalised)",
)
