#!/usr/bin/env python3
"""Run full comparison analysis for an experiment."""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]  # check/scripts/visu/1-baseline -> project root
sys.path.insert(0, str(SCRIPT_DIR.parent))  # visu/ for tables, style, etc.

from tables import run_full_comparison

# --- Configuration ---
EXPERIMENT_DIR = PROJECT_ROOT / "check/results/1-baseline-compact1"
BASELINE = "benders_FLLP"
PD_CONVERGENCE = False

STYLE_MAP = {
    "benders_FLLP": {"color": "#000000"},               # black
    "compact":      {"color": "#E69F00"},               # orange
}

LABEL_MAP = {
    "benders_FLLP":     r"\texttt{Benders}",
    "compact":          r"\texttt{Compact}",
}

INCLUDE_SETTINGS = [
    "benders_FLLP",
    "compact",
]
# ---------------------

MIN_TIME = 0  # drop instances where all settings finish below this (seconds)

run_full_comparison(EXPERIMENT_DIR, baseline=BASELINE,
                    pd_convergence=PD_CONVERGENCE, style_map=STYLE_MAP,
                    label_map=LABEL_MAP, min_buses=73,
                    include_settings=INCLUDE_SETTINGS,
                    min_time=MIN_TIME)
