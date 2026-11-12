#!/usr/bin/env python3
"""Grid-search heatmap for the 1-cf-spsel-hybrid experiment (CF x SP thresholds).

Uses the shared preprocessing helpers from ``reader`` so the displayed
geometric means and instance counts match ``comparison_short.md``.
"""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[4]  # check/scripts/visu/4-cut_filtering -> project root
sys.path.insert(0, str(SCRIPT_DIR.parents[1]))  # visu/ for heatmap, reader, etc.

from heatmap import grid_search_heatmap
from reader import read_experiment, cleanup_data

# --- Configuration (mirrors compare_1-cf-spsel-hybrid.py) ---
EXPERIMENT_DIR = PROJECT_ROOT / "check/results/1-cf-spsel-hybrid"
BASELINE = "cf100_sp100"          # cell to highlight in the heatmap
PATTERN = r"cf(\d+)_sp(\d+)"   # group 1 -> X (CF%), group 2 -> Y (SP%)
TARGET = "time_solve"
MIN_BUSES = 73
MIN_TIME = 30.0

INCLUDE_SETTINGS = [
    "cf001_sp001", "cf001_sp005",  "cf001_sp010",  "cf001_sp025",  "cf001_sp050", "cf001_sp075",  "cf001_sp100",
    "cf005_sp005", "cf005_sp010",  "cf005_sp025",  "cf005_sp050",  "cf005_sp075", "cf005_sp100",
    "cf010_sp010", "cf010_sp025",  "cf010_sp050",  "cf010_sp075",  "cf010_sp100",
    "cf025_sp025", "cf025_sp050",  "cf025_sp075",  "cf025_sp100",
    "cf050_sp050", "cf050_sp075",  "cf050_sp100",
    "cf075_sp075", "cf075_sp100",
    "cf100_sp100",
]

# ------------------------------------------------------------

df = read_experiment(EXPERIMENT_DIR)
df = df[df["settings"].isin(INCLUDE_SETTINGS)].reset_index(drop=True)
df, _ = cleanup_data(df, exclude_all_timelimit=False, baseline=BASELINE)

grid_search_heatmap(
    df,
    output=EXPERIMENT_DIR / "heatmap_time_solve.pdf",
    pattern=PATTERN,
    target=TARGET,
    baseline=BASELINE,
    common_only=True,
    xlabel=r"CF threshold (\%)",
    ylabel=r"SP threshold (\%)",
)
print(f"Wrote {EXPERIMENT_DIR / 'heatmap_time_solve.pdf'}")
