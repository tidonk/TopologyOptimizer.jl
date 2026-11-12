#!/usr/bin/env python3
"""Compare all settings from the 1-cf-spsel-hybrid grid search (Exp. 4 Stage 1)."""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[4]  # check/scripts/visu/4-cut_filtering -> project root
sys.path.insert(0, str(SCRIPT_DIR.parents[1]))  # visu/ for tables, style, etc.

from tables import run_full_comparison

# --- Configuration ---
EXPERIMENT_DIR = PROJECT_ROOT / "check/results/1-cf-spsel-hybrid"
BASELINE = "cf001_sp100"
PD_CONVERGENCE = False

# Colors: darker = more CF, cooler hue = more SP
STYLE_MAP = {
    "cf001_sp005":  {"color": "#000000", "linestyle": "--"},    # cf=1%
    "cf001_sp010":  {"color": "#252525", "linestyle": "--"},
    "cf001_sp025":  {"color": "#525252", "linestyle": "--"},
    "cf001_sp050":  {"color": "#737373", "linestyle": "--"},
    "cf001_sp075":  {"color": "#969696", "linestyle": "--"},
    "cf001_sp100":  {"color": "#bdbdbd", "linestyle": "--"},
    "cf005_sp010":  {"color": "#005a32", "linestyle": "-."},    # cf=5%
    "cf005_sp025":  {"color": "#238b45"},
    "cf005_sp050":  {"color": "#41ab5d"},
    "cf005_sp075":  {"color": "#74c476"},
    "cf005_sp100":  {"color": "#a1d99b"},
    "cf010_sp025":  {"color": "#084594"},                       # cf=10%
    "cf010_sp050":  {"color": "#2171b5"},
    "cf010_sp075":  {"color": "#4292c6"},
    "cf010_sp100":  {"color": "#6baed6"},
    "cf025_sp050":  {"color": "#7f0000"},                       # cf=25%
    "cf025_sp075":  {"color": "#cb181d"},
    "cf025_sp100":  {"color": "#fb6a4a"},
    "cf050_sp075":  {"color": "#7a0177"},                       # cf=50%
    "cf050_sp100":  {"color": "#c51b8a"},
    "cf075_sp100":  {"color": "#E69F00"},                       # cf=75%
}

LABEL_MAP = {
    "cf001_sp005":  r"\texttt{CF=1\% / SP=5\%}",
    "cf001_sp010":  r"\texttt{CF=1\% / SP=10\%}",
    "cf001_sp025":  r"\texttt{CF=1\% / SP=25\%}",
    "cf001_sp050":  r"\texttt{CF=1\% / SP=50\%}",
    "cf001_sp075":  r"\texttt{CF=1\% / SP=75\%}",
    "cf001_sp100":  r"\texttt{CF=1\% / SP=100\%}",
    "cf005_sp010":  r"\texttt{CF=5\% / SP=10\%}",
    "cf005_sp025":  r"\texttt{CF=5\% / SP=25\%}",
    "cf005_sp050":  r"\texttt{CF=5\% / SP=50\%}",
    "cf005_sp075":  r"\texttt{CF=5\% / SP=75\%}",
    "cf005_sp100":  r"\texttt{CF=5\% / SP=100\%}",
    "cf010_sp025":  r"\texttt{CF=10\% / SP=25\%}",
    "cf010_sp050":  r"\texttt{CF=10\% / SP=50\%}",
    "cf010_sp075":  r"\texttt{CF=10\% / SP=75\%}",
    "cf010_sp100":  r"\texttt{CF=10\% / SP=100\%}",
    "cf025_sp050":  r"\texttt{CF=25\% / SP=50\%}",
    "cf025_sp075":  r"\texttt{CF=25\% / SP=75\%}",
    "cf025_sp100":  r"\texttt{CF=25\% / SP=100\%}",
    "cf050_sp075":  r"\texttt{CF=50\% / SP=75\%}",
    "cf050_sp100":  r"\texttt{CF=50\% / SP=100\%}",
    "cf075_sp100":  r"\texttt{CF=75\% / SP=100\%}",
}

INCLUDE_SETTINGS = [
    "cf001_sp005",  "cf001_sp010",  "cf001_sp025",  "cf001_sp050",
    "cf001_sp075",  "cf001_sp100",
    "cf005_sp010",  "cf005_sp025",  "cf005_sp050",  "cf005_sp075",
    "cf005_sp100",
    "cf010_sp025",  "cf010_sp050",  "cf010_sp075",  "cf010_sp100",
    "cf025_sp050",  "cf025_sp075",  "cf025_sp100",
    "cf050_sp075",  "cf050_sp100",
    "cf075_sp100",
]
# ---------------------

MIN_TIME = 0  # drop instances where all settings finish below this (seconds)

run_full_comparison(EXPERIMENT_DIR, baseline=BASELINE,
                    pd_convergence=PD_CONVERGENCE, style_map=STYLE_MAP,
                    label_map=LABEL_MAP, min_buses=73,
                    include_settings=INCLUDE_SETTINGS,
                    min_time=MIN_TIME)
