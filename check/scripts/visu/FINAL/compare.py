#!/usr/bin/env python3
"""Compare settings from Exp. 5 - 5-partial_decomp.

Auto-generated comparison for Exp. 5 - 5-partial_decomp.

Outputs include both compact summaries:
- comparison_short.md      (instances solved by all settings)
- comparison_short_all.md  (all instances; non-OPTIMAL/missing penalized)
"""

import sys
import pandas as pd
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]  # check/scripts/visu/5-partial_decomp -> project root
sys.path.insert(0, str(SCRIPT_DIR.parent))  # visu/ for tables, style, etc.

from tables import run_full_comparison

# --- Configuration ---
EXP_NAME = "X-final-full"
EXPERIMENT_DIR = PROJECT_ROOT / "check/results" / EXP_NAME
BASELINE = ""
PD_CONVERGENCE = False
MIN_TIME = 0

STYLE_MAP = {
    "vad_num_1": {"color": "#000000"},
    "overload_num_1": {"color": "#0072B2"},
    "overload_pct_5": {"color": "#56B4E9"},
    "overload_pct_10": {"color": "#009E73"},
    "overload_pct_25": {"color": "#CC79A7"},
    "overload_pct_50": {"color": "#D55E00"},
    "overload_pct_75": {"color": "#E69F00"},
    "overload_pct_90": {"color": "#7F3C8D"},
}

LABEL_MAP = {
    "compact":         r"\texttt{Compact}",
    "benders_FLLP":    r"\texttt{Benders}",
    "benders_LEQRF":   r"\texttt{Benders SLE-U}",
    "div_decay_2":     r"\texttt{SP Selection}",
    "cutfiltering":    r"\texttt{Cut Filtering}",
    "div_decay_5_005": r"\texttt{SP hybrid, no CF}",
    "base_only_div":   r"\texttt{Partial Decomp.}"
}

INCLUDE_SETTINGS = [
    "compact",
    "benders_FLLP",
    "benders_LEQRF",
    "div_decay_2",
    "cutfiltering",
    "div_decay_5_005",
    "base_only_div"
]
# ---------------------

def _add_time_with_warmstart(df: pd.DataFrame) -> pd.DataFrame:
    """Derive time_with_warmstart = time_solve + time_warmstart.

    Only cprox, eta and sens settings use warmstarting; for all others
    warmstart overhead is zero and time_with_warmstart equals time_solve.
    """
    df = df.copy()
    t_solve = pd.to_numeric(df.get("time_solve", 0), errors="coerce").fillna(0)
    t_ws = pd.to_numeric(df.get("time_warmstart", 0), errors="coerce").fillna(0)
    # include warmstart time for EVERY run (time_warmstart is 0 for
    # non-warmstarted settings, so this is a no-op there)
    df["time_with_warmstart"] = t_solve + t_ws
    return df


if not EXPERIMENT_DIR.is_dir() or not any(EXPERIMENT_DIR.iterdir()):
    print(f"warning: no results in {EXPERIMENT_DIR}; skipping comparison")
else:
    try:
        run_full_comparison(EXPERIMENT_DIR, baseline=BASELINE,
                            pd_convergence=PD_CONVERGENCE,# style_map=STYLE_MAP,
                            label_map=LABEL_MAP, min_buses=73,
                            include_settings=INCLUDE_SETTINGS,
                            min_time=MIN_TIME,
                            target="time_with_warmstart",
                            perf_profile_top_n=10,
                            df_transform=_add_time_with_warmstart)
    except FileNotFoundError as e:
        print(f"warning: {e}")