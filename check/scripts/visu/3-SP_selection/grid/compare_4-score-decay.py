#!/usr/bin/env python3
"""Run full comparison analysis for the 4-score-decay experiment."""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[4]  # check/scripts/visu/3-SP_selection -> project root
sys.path.insert(0, str(SCRIPT_DIR.parents[1]))  # visu/ for tables, style, etc.

from tables import run_full_comparison

# --- Configuration ---
EXPERIMENT_DIR = PROJECT_ROOT / "check/results/4-score-decay"
BASELINE = "rel_nodecay"
PD_CONVERGENCE = False

STYLE_MAP = {
    "decay_0001":  {"color": "#000000"},               # black
    "decay_010":   {"color": "#E69F00"},               # orange
    "decay_025":   {"color": "#56B4E9"},               # light blue
    "decay_050":   {"color": "#009E73"},               # green
    "decay_075":   {"color": "#0072B2"},               # blue
    "rel_decay":   {"color": "#D55E00"},               # red
    "rel_nodecay": {"color": "#CC79A7", "linestyle": "dashed"}  # pink, dashed
}

LABEL_MAP = {
    "decay_0001":  r"\texttt{Decay 0.001}",
    "decay_010":   r"\texttt{Decay 0.1}",
    "decay_025":   r"\texttt{Decay 0.25}",
    "decay_050":   r"\texttt{Decay 0.5}",
    "decay_075":   r"\texttt{Decay 0.75}",
    "rel_decay":   r"\texttt{Decay 0.90}",
    "rel_nodecay": r"\texttt{No Decay}"
}

INCLUDE_SETTINGS = [
    "decay_0001",
    "decay_010",
    "decay_025",
    "decay_050",
    "decay_075",
    "rel_decay",
    "rel_nodecay"
]
# ---------------------

MIN_TIME = 0  # drop instances where all settings finish below this (seconds)

run_full_comparison(EXPERIMENT_DIR, baseline=BASELINE,
                    pd_convergence=PD_CONVERGENCE, style_map=STYLE_MAP,
                    label_map=LABEL_MAP, min_buses=73,
                    include_settings=INCLUDE_SETTINGS,
                    min_time=MIN_TIME)
