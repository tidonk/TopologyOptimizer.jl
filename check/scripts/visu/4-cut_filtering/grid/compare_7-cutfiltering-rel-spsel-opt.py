#!/usr/bin/env python3
"""Compare settings from Exp. 4 Stage 2 - 7-cutfiltering-rel-spsel-opt.

Restricted to the ``_005`` subset (``STOP_AFTER_N_HITS_REL = 0.05``); baseline
is the SP-only ``cf100_sp10_005`` import from 3-SP_selection.
"""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[4]  # check/scripts/visu/4-cut_filtering -> project root
sys.path.insert(0, str(SCRIPT_DIR.parents[1]))  # visu/ for tables, style, etc.
sys.path.insert(0, str(SCRIPT_DIR.parent))         # for _cutfiltering_common

from tables import run_full_comparison  # noqa: E402
from _cutfiltering_common import build_label_map, build_style_map  # noqa: E402

# --- Configuration ---
EXP_NAME = "7-cutfiltering-rel-spsel-opt"
EXPERIMENT_DIR = PROJECT_ROOT / "check/results" / EXP_NAME
BASELINE = "cf100_sp100"  # raw reference: no CF + no SP selection (imported from 4-cutfiltering-rel-spsel-rel)
PD_CONVERGENCE = False
MIN_TIME = 0

INCLUDE_SETTINGS = [
    "cf001_sp02_005", "cf005_sp02_005", "cf010_sp02_005", "cf010_sp02_010", "cf100_sp02_005",
    "cf001_sp05_005", "cf001_sp05_010", "cf001_sp10_005", "cf001_sp10_010",
    "cf005_sp05_005", "cf005_sp05_010", "cf005_sp10_005", "cf005_sp10_010",
    "cf010_sp05_005", "cf010_sp05_010", "cf010_sp10_005", "cf010_sp10_010",
    "cf100_sp05_005", "cf100_sp10_005",
    # CF-only baselines (no SP selection: STOP_AFTER_N_HITS_REL = 1.0,
    # STOP_ON_PB = false, STOP_AFTER_N_FEAS_CUTS = Inf) imported from
    # 4-cutfiltering-rel-spsel-rel; cf100_sp100 = no CF + no SP selection.
    "cf001_sp100", "cf005_sp100", "cf010_sp100", "cf100_sp100",
]
STYLE_MAP = build_style_map(INCLUDE_SETTINGS)
LABEL_MAP = build_label_map(INCLUDE_SETTINGS)
# ---------------------

if not EXPERIMENT_DIR.is_dir() or not any(EXPERIMENT_DIR.iterdir()):
    print(f"warning: no results in {EXPERIMENT_DIR}; skipping comparison")
else:
    try:
        run_full_comparison(EXPERIMENT_DIR, baseline=BASELINE,
                            pd_convergence=PD_CONVERGENCE, style_map=STYLE_MAP,
                            label_map=LABEL_MAP, min_buses=73,
                            include_settings=INCLUDE_SETTINGS,
                            min_time=MIN_TIME)
    except FileNotFoundError as e:
        print(f'warning: {e}')
