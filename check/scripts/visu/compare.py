#!/usr/bin/env python3
"""Run full comparison analysis for an experiment."""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[2]  # check/scripts/visu -> project root
sys.path.insert(0, str(SCRIPT_DIR))

from tables import run_full_comparison

# --- Configuration ---
EXPERIMENT_DIR = PROJECT_ROOT / "check/results/experiment"
BASELINE = "no_cutsel"
PD_CONVERGENCE = False

STYLE_MAP = {
    "no_cutsel":                        {"color": "#000000"},               # black
    "cutsel_random":                    {"color": "#E69F00"},               # orange
    "cutsel_viol_abs_num":              {"color": "#0072B2", "linestyle": "--"},
    "cutsel_viol_abs_num_feasfirst":    {"color": "#0072B2"},
    "cutsel_viol_abs_num_feasfirst2":   {"color": "#0072B2", "linestyle": ":"},
    "cutsel_viol_abs_num_feasfirst3":   {"color": "#0072B2"},               # blue
    "cutsel_viol_pb":                   {"color": "#56B4E9"},
    "cutsel_diversity_kmeans":          {"color": "#009E73"},
    "cutsel_diversity_dbscan":          {"color": "#009E73", "linestyle": "--"},
    "cutsel_diversity_dbscan_eqfreq":   {"color": "#009E73", "linestyle": ":"},
    "cutsel_diversity_dbscan_jaccard":  {"color": "#009E73", "linestyle": "-."},
    "cutsel_diversity_kmeans_eqfreq":   {"color": "#009E73"},
    "cutsel_diversity_kmeans_jaccard":  {"color": "#009E73", "linestyle": "--"},  # green dashed
    "cutsel_diversity_kmeans_cosine":   {"color": "#009E73", "linestyle": ":"},   # green dotted
    "cutsel_hybrid_dbscan_eqfreq":     {"color": "#D55E00", "linestyle": "--"},
    "cutsel_hybrid_dbscan_jaccard":     {"color": "#D55E00", "linestyle": ":"},
    "cutsel_hybrid_kmeans_eqfreq":     {"color": "#D55E00"},                # vermillion
    "cutsel_hybrid_kmeans_jaccard":    {"color": "#D55E00", "linestyle": "--"},  # vermillion dashed
    "cutsel_hybrid_kmeans_cosine":     {"color": "#CC79A7"},                # pink
}

LABEL_MAP = {
    "no_cutsel":                r"\texttt{default}",
    "cutsel_viol_rel_005":      r"\texttt{violation}",
    "cutsel_viol_agg_rel_005":      r"\texttt{violation+}",
    "cutsel_viol_agg_rel_005_w":      r"\texttt{violation+w}",
    "cutsel_viol_agg_rel_005_wn":      r"\texttt{violation+}",
    "cutsel_viol_pb_101":       r"\texttt{violation\_adaptive}",
    "cutsel_random":            r"\texttt{random}",
    "cutsel_random_rel_005":            r"\texttt{random}",
    "cutsel_random_agg_rel_005":            r"\texttt{random+}",
    "cutsel_random_agg_rel_005_w":            r"\texttt{random+w}",
    "cutsel_random_agg_rel_005_wn":            r"\texttt{random+}",
    "cutsel_diversity_kmeans_cosine_rel_005_centroid": r"\texttt{diversity}",
    "cutsel_diversity_kmeans_cosine_agg_hyb_k1_centroid": r"\texttt{diversity+}",
    "cutsel_diversity_kmeans_cosine_agg_hyb_k1_centroid_w": r"\texttt{diversity+w}",
    "cutsel_diversity_kmeans_cosine_agg_hyb_k1_centroid_wn": r"\texttt{diversity+}",
    "cutsel_diversity_kmeans_cosine_rel_005": r"\texttt{hybrid}",
    "cutsel_diversity_kmeans_cosine_agg_hyb_k1": r"\texttt{hybrid+}",
    "cutsel_diversity_kmeans_cosine_agg_hyb_k1_w": r"\texttt{hybrid+w}",
    "cutsel_diversity_kmeans_cosine_agg_hyb_k1_wn": r"\texttt{hybrid+}",
}

INCLUDE_SETTINGS = [
    "no_cutsel",
    "cutsel_random_rel_005",
    #"cutsel_random_agg_rel_005",
    #"cutsel_random_agg_rel_005_w",
    "cutsel_random_agg_rel_005_wn",
    "cutsel_viol_rel_005",
    #"cutsel_viol_agg_rel_005",
    #"cutsel_viol_agg_rel_005_w",
    "cutsel_viol_agg_rel_005_wn",
    "cutsel_diversity_kmeans_cosine_rel_005_centroid",
    #"cutsel_diversity_kmeans_cosine_agg_hyb_k1_centroid",
    #"cutsel_diversity_kmeans_cosine_agg_hyb_k1_centroid_w",
    "cutsel_diversity_kmeans_cosine_agg_hyb_k1_centroid_wn",
    "cutsel_diversity_kmeans_cosine_rel_005",
    #"cutsel_diversity_kmeans_cosine_agg_hyb_k1",
    #"cutsel_diversity_kmeans_cosine_agg_hyb_k1_w",
    "cutsel_diversity_kmeans_cosine_agg_hyb_k1_wn",
]
# ---------------------

MIN_TIME = 0  # drop instances where all settings finish below this (seconds)

run_full_comparison(EXPERIMENT_DIR, baseline=BASELINE,
                    pd_convergence=PD_CONVERGENCE,
                    label_map=LABEL_MAP, min_buses=73,
                    include_settings=INCLUDE_SETTINGS,
                    min_time=MIN_TIME)
