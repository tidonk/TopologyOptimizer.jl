#!/usr/bin/env python3
"""Compare all settings from the 1-cutfiltering experiment (Exp. 4 Stage 1)."""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[4]  # check/scripts/visu/4-cut_filtering -> project root
sys.path.insert(0, str(SCRIPT_DIR.parents[1]))  # visu/ for tables, style, etc.

from tables import run_full_comparison

# --- Configuration ---
EXPERIMENT_DIR = PROJECT_ROOT / "check/results/1-cutfiltering"
BASELINE = "viol"
PD_CONVERGENCE = False

STYLE_MAP = {
    "viol":        {"color": "#000000"},
    "rnd":         {"color": "#E69F00"},
    "agg_hyb_k1":  {"color": "#0072B2"},
}

LABEL_MAP = {
    "viol":        r"\texttt{violation}",
    "rnd":         r"\texttt{random}",
    "agg_hyb_k1":  r"\texttt{hybrid+}",
}

INCLUDE_SETTINGS = [
    "viol",
    "rnd",
    "agg_hyb_k1",
]
# ---------------------

MIN_TIME = 0  # keep all currently available instances for this small experiment

run_full_comparison(EXPERIMENT_DIR, baseline=BASELINE,
                    pd_convergence=PD_CONVERGENCE, style_map=STYLE_MAP,
                    label_map=LABEL_MAP, min_buses=73,
                    include_settings=INCLUDE_SETTINGS,
                    min_time=MIN_TIME)
