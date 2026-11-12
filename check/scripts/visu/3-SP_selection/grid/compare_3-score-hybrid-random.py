#!/usr/bin/env python3
"""Run full comparison analysis for the 3-score-hybrid-random experiment (Rel vs Random)."""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[4]  # check/scripts/visu/3-SP_selection -> project root
sys.path.insert(0, str(SCRIPT_DIR.parents[1]))  # visu/ for tables, style, etc.

from tables import run_full_comparison

# --- Configuration ---
EXPERIMENT_DIR = PROJECT_ROOT / "check/results/3-score-hybrid-random"
BASELINE = "hybrid_rel_rnd_50_50"
PD_CONVERGENCE = False

STYLE_MAP = {
    "hybrid_rel_rnd_90_10":  {"color": "#005a32"},               # dark green (rel-heavy)
    "hybrid_rel_rnd_75_25":  {"color": "#238b45"},
    "hybrid_rel_rnd_50_50":  {"color": "#000000"},               # black (equal)
    "hybrid_rel_rnd_25_75":  {"color": "#E69F00"},               # orange
    "hybrid_rel_rnd_10_90":  {"color": "#D55E00"},               # red   (rnd-heavy)
    "rnd_rel_010":          {"color": "#000000", "linestyle": "dotted"}, # dotted black
    "rel_decay":            {"color": "#009E73"},               # green
}

LABEL_MAP = {
    "hybrid_rel_rnd_90_10":  r"\texttt{Rel 90 / Rnd 10}",
    "hybrid_rel_rnd_75_25":  r"\texttt{Rel 75 / Rnd 25}",
    "hybrid_rel_rnd_50_50":  r"\texttt{Rel 50 / Rnd 50}",
    "hybrid_rel_rnd_25_75":  r"\texttt{Rel 25 / Rnd 75}",
    "hybrid_rel_rnd_10_90":  r"\texttt{Rel 10 / Rnd 90}",
    "rnd_rel_010":           r"\texttt{Rnd 100}",
    "rel_decay":             r"\texttt{Rel 100}",
}

INCLUDE_SETTINGS = [
    "hybrid_rel_rnd_90_10",
    "hybrid_rel_rnd_75_25",
    "hybrid_rel_rnd_50_50",
    "hybrid_rel_rnd_25_75",
    "hybrid_rel_rnd_10_90",
    "rnd_rel_010",
    "rel_decay"
]
# ---------------------

MIN_TIME = 0  # drop instances where all settings finish below this (seconds)

run_full_comparison(EXPERIMENT_DIR, baseline=BASELINE,
                    pd_convergence=PD_CONVERGENCE, style_map=STYLE_MAP,
                    label_map=LABEL_MAP, min_buses=73,
                    include_settings=INCLUDE_SETTINGS,
                    min_time=MIN_TIME)
