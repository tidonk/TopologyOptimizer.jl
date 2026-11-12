#!/usr/bin/env python3
"""Run full comparison analysis for the 1-n_hits experiment."""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[4]  # check/scripts/visu/3-SP_selection -> project root
sys.path.insert(0, str(SCRIPT_DIR.parents[1]))  # visu/ for tables, style, etc.

from tables import run_full_comparison

# --- Configuration ---
EXPERIMENT_DIR = PROJECT_ROOT / "check/results/1-n_hits_pb"
BASELINE = "rnd_rel_100"
PD_CONVERGENCE = False

STYLE_MAP = {
    "rnd_rel_001":  {"color": "#000000"},               # black
    "rnd_rel_005":  {"color": "#E69F00"},               # orange
    "rnd_rel_010":  {"color": "#56B4E9"},               # light blue
    "rnd_rel_025":  {"color": "#009E73"},               # green
    "rnd_rel_050":  {"color": "#0072B2"},               # blue
    "rnd_rel_075":  {"color": "#CC79A7"},               # purple
    "rnd_rel_100":  {"color": "#D55E00"},               # red
    "rnd_abs_1":    {"color": "#000000", "linestyle": "dashed"}, # dashed black
    "rnd_abs_2":    {"color": "#E69F00", "linestyle": "dashed"}, # dashed orange
    "rnd_abs_5":    {"color": "#56B4E9", "linestyle": "dashed"}, # dashed light blue
    "rnd_abs_10":   {"color": "#009E73", "linestyle": "dashed"}, # dashed green
    "rnd_pb":       {"color": "#0072B2", "linestyle": "dashed"}, # dashed blue
}

LABEL_MAP = {
    "rnd_rel_001":  r"\texttt{1\%}",
    "rnd_rel_005":  r"\texttt{5\%}",
    "rnd_rel_010":  r"\texttt{10\%}",
    "rnd_rel_025":  r"\texttt{25\%}",
    "rnd_rel_050":  r"\texttt{50\%}",
    "rnd_rel_075":  r"\texttt{75\%}",
    "rnd_rel_100":  r"\texttt{100\%}",
    "rnd_abs_1":    r"\texttt{1}",
    "rnd_abs_2":    r"\texttt{2}",
    "rnd_abs_5":    r"\texttt{5}",
    "rnd_abs_10":   r"\texttt{10}",
    "rnd_pb":       r"\texttt{PB}"
}

INCLUDE_SETTINGS = [
    "rnd_rel_001",
    #"rnd_rel_005",
    #"rnd_rel_010",
    #"rnd_rel_025",
    #"rnd_rel_050",
    #"rnd_rel_075",
    "rnd_rel_100",  # baseline
    #"rnd_abs_1",
    "rnd_abs_2",
    "rnd_abs_5",
    "rnd_abs_10",
    "rnd_pb"
]
# ---------------------

MIN_TIME = 0  # drop instances where all settings finish below this (seconds)

run_full_comparison(EXPERIMENT_DIR, baseline=BASELINE,
                    pd_convergence=PD_CONVERGENCE, style_map=STYLE_MAP,
                    label_map=LABEL_MAP, min_buses=73,
                    include_settings=INCLUDE_SETTINGS,
                    min_time=MIN_TIME)
