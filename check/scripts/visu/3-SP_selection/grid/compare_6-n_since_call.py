#!/usr/bin/env python3
"""Run full comparison analysis for the 5-score-minprio experiment."""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[4]  # check/scripts/visu/3-SP_selection -> project root
sys.path.insert(0, str(SCRIPT_DIR.parents[1]))  # visu/ for tables, style, etc.

from tables import run_full_comparison

# --- Configuration ---
EXPERIMENT_DIR = PROJECT_ROOT / "check/results/5-score-minprio"
BASELINE = "score_minprio_050"
PD_CONVERGENCE = False

STYLE_MAP = {
    "score_minprio_025":  {"color": "#E69F00"},               # orange
    "score_minprio_050":  {"color": "#000000"},               # black
    "score_minprio_075":  {"color": "#009E73"},               # green
    "score_minprio_090":  {"color": "#0072B2"},               # blue
}

LABEL_MAP = {
    "score_minprio_025":  r"\texttt{Minprio 0.25}",
    "score_minprio_050":  r"\texttt{Minprio 0.50}",
    "score_minprio_075":  r"\texttt{Minprio 0.75}",
    "score_minprio_090":  r"\texttt{Minprio 0.90}",
}

INCLUDE_SETTINGS = [
    "score_minprio_025",
    "score_minprio_050",
    "score_minprio_075",
    "score_minprio_090",
]
# ---------------------

MIN_TIME = 0  # drop instances where all settings finish below this (seconds)

run_full_comparison(EXPERIMENT_DIR, baseline=BASELINE,
                    pd_convergence=PD_CONVERGENCE, style_map=STYLE_MAP,
                    label_map=LABEL_MAP, min_buses=73,
                    include_settings=INCLUDE_SETTINGS,
                    min_time=MIN_TIME)
