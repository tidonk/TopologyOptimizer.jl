#!/usr/bin/env python3
"""Run full comparison analysis for the 3-score-hybrid experiment (all pairs)."""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[4]  # check/scripts/visu/3-SP_selection -> project root
sys.path.insert(0, str(SCRIPT_DIR.parents[1]))  # visu/ for tables, style, etc.

from tables import run_full_comparison

# --- Configuration ---
EXPERIMENT_DIR = PROJECT_ROOT / "check/results/3-score-hybrid"
BASELINE = "rnd_abs_2"
PD_CONVERGENCE = False

STYLE_MAP = {
    # viol vs rel (reds/oranges)
    "hybrid_viol_rel_90_10":   {"color": "#D55E00"},               # red    (viol-heavy)
    "hybrid_viol_rel_75_25":   {"color": "#E69F00"},               # orange
    "hybrid_viol_rel_50_50":   {"color": "#000000"},               # black  (equal)
    "hybrid_viol_rel_25_75":   {"color": "#56B4E9"},               # light blue
    "hybrid_viol_rel_10_90":   {"color": "#009E73"},               # green  (rel-heavy)
    # viol vs share (purples/mauves)
    "hybrid_viol_share_90_10": {"color": "#7b2d8b"},               # purple (viol-heavy)
    "hybrid_viol_share_75_25": {"color": "#9e4fb5"},
    "hybrid_viol_share_50_50": {"color": "#b07cc6"},               # mid purple
    "hybrid_viol_share_25_75": {"color": "#c9a3d8"},
    "hybrid_viol_share_10_90": {"color": "#0072B2"},               # blue   (share-heavy)
    # rel vs share (greens/teals)
    "hybrid_rel_share_90_10":  {"color": "#005a32"},               # dark green (rel-heavy)
    "hybrid_rel_share_75_25":  {"color": "#238b45"},
    "hybrid_rel_share_50_50":  {"color": "#41ab5d"},               # mid green
    "hybrid_rel_share_25_75":  {"color": "#CC79A7"},               # pink
    "hybrid_rel_share_10_90":  {"color": "#2c7fb8"},               # teal  (share-heavy)
    # pure scores
    "viol_decay":              {"color": "#D55E00"},               # red
    "rel_decay":               {"color": "#009E73"},               # green
    "share_decay":             {"color": "#0072B2"},               # blue
    "rnd_abs_2":             {"color": "#000000", "linestyle": "dotted"}, # dotted black
}

LABEL_MAP = {
    "hybrid_viol_rel_90_10":   r"\texttt{Viol 90 / Rel 10}",
    "hybrid_viol_rel_75_25":   r"\texttt{Viol 75 / Rel 25}",
    "hybrid_viol_rel_50_50":   r"\texttt{Viol 50 / Rel 50}",
    "hybrid_viol_rel_25_75":   r"\texttt{Viol 25 / Rel 75}",
    "hybrid_viol_rel_10_90":   r"\texttt{Viol 10 / Rel 90}",
    "hybrid_viol_share_90_10": r"\texttt{Viol 90 / Share 10}",
    "hybrid_viol_share_75_25": r"\texttt{Viol 75 / Share 25}",
    "hybrid_viol_share_50_50": r"\texttt{Viol 50 / Share 50}",
    "hybrid_viol_share_25_75": r"\texttt{Viol 25 / Share 75}",
    "hybrid_viol_share_10_90": r"\texttt{Viol 10 / Share 90}",
    "hybrid_rel_share_90_10":  r"\texttt{Rel 90 / Share 10}",
    "hybrid_rel_share_75_25":  r"\texttt{Rel 75 / Share 25}",
    "hybrid_rel_share_50_50":  r"\texttt{Rel 50 / Share 50}",
    "hybrid_rel_share_25_75":  r"\texttt{Rel 25 / Share 75}",
    "hybrid_rel_share_10_90":  r"\texttt{Rel 10 / Share 90}",
    "hybrid_rel_rnd_90_10":    r"\texttt{Rel 90 / Rnd 10}",
    "hybrid_rel_rnd_75_25":    r"\texttt{Rel 75 / Rnd 25}",
    "hybrid_rel_rnd_50_50":    r"\texttt{Rel 50 / Rnd 50}",
    "hybrid_rel_rnd_25_75":    r"\texttt{Rel 25 / Rnd 75}",
    "hybrid_rel_rnd_10_90":    r"\texttt{Rel 10 / Rnd 90}",
    "hybrid_viol_rnd_90_10":   r"\texttt{Viol 90 / Rnd 10}",
    "hybrid_viol_rnd_75_25":   r"\texttt{Viol 75 / Rnd 25}",
    "hybrid_viol_rnd_50_50":   r"\texttt{Viol 50 / Rnd 50}",
    "hybrid_viol_rnd_25_75":   r"\texttt{Viol 25 / Rnd 75}",
    "hybrid_viol_rnd_10_90":   r"\texttt{Viol 10 / Rnd 90}",
    "hybrid_share_rnd_90_10":  r"\texttt{Share 90 / Rnd 10}",
    "hybrid_share_rnd_75_25":  r"\texttt{Share 75 / Rnd 25}",
    "hybrid_share_rnd_50_50":  r"\texttt{Share 50 / Rnd 50}",
    "hybrid_share_rnd_25_75":  r"\texttt{Share 25 / Rnd 75}",
    "hybrid_share_rnd_10_90":  r"\texttt{Share 10 / Rnd 90}",
    "viol_decay":              r"\texttt{Viol 100}",
    "rel_decay":               r"\texttt{Rel 100}",
    "share_decay":             r"\texttt{Share 100}",
    "rnd_abs_2":             r"\texttt{Random}"
}

INCLUDE_SETTINGS = [
    "rnd_abs_2", "viol_decay", "rel_decay", "share_decay"
    "hybrid_viol_rel_90_10",   "hybrid_viol_rel_75_25",   "hybrid_viol_rel_50_50",    "hybrid_viol_rel_25_75",   "hybrid_viol_rel_10_90",
    "hybrid_viol_share_90_10", "hybrid_viol_share_75_25", "hybrid_viol_share_50_50",  "hybrid_viol_share_25_75", "hybrid_viol_share_10_90",
    "hybrid_rel_share_90_10",  "hybrid_rel_share_75_25",  "hybrid_rel_share_50_50",   "hybrid_rel_share_25_75",  "hybrid_rel_share_10_90",
    "hybrid_viol_rnd_90_10",   "hybrid_viol_rnd_75_25",   "hybrid_viol_rnd_50_50",    "hybrid_viol_rnd_25_75",   "hybrid_viol_rnd_10_90",
    "hybrid_share_rnd_90_10",  "hybrid_share_rnd_75_25",  "hybrid_share_rnd_50_50",   "hybrid_share_rnd_25_75",  "hybrid_share_rnd_10_90",
    "hybrid_rel_rnd_90_10",    "hybrid_rel_rnd_75_25",    "hybrid_rel_rnd_50_50",     "hybrid_rel_rnd_25_75",    "hybrid_rel_rnd_10_90",
]
# ---------------------

MIN_TIME = 0  # drop instances where all settings finish below this (seconds)

run_full_comparison(EXPERIMENT_DIR, baseline=BASELINE,
                    pd_convergence=PD_CONVERGENCE, style_map=STYLE_MAP,
                    label_map=LABEL_MAP, min_buses=73,
                    include_settings=INCLUDE_SETTINGS,
                    min_time=MIN_TIME)
