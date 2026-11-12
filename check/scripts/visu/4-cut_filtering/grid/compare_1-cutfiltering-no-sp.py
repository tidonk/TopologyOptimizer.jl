#!/usr/bin/env python3
"""Compare settings from Exp. 4 Stage 1, no subproblem selection.

Edited auto-generated comparison for Exp. 4 Stage 1 - 3-cutfiltering-abs-spsel-abs.
"""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[4]  # check/scripts/visu/4-cut_filtering -> project root
sys.path.insert(0, str(SCRIPT_DIR.parents[1]))  # visu/ for tables, style, etc.
sys.path.insert(0, str(SCRIPT_DIR.parent))         # for _cutfiltering_common

from tables import run_full_comparison  # noqa: E402
from _cutfiltering_common import (       # noqa: E402
    build_label_map, build_style_map, discover_settings,
)

# --- Configuration ---
EXP_NAME = "1-cutfiltering-no-sp"
EXPERIMENT_DIR = PROJECT_ROOT / "check/results" / EXP_NAME
SETTINGS_DIR = PROJECT_ROOT / "settings/4-cut_filtering" / EXP_NAME
BASELINE = "cf100_sp100_abs"
PD_CONVERGENCE = False
MIN_TIME = 0

INCLUDE_SETTINGS = [setting for setting in discover_settings(EXPERIMENT_DIR, SETTINGS_DIR) if not setting in ["cf01_sp100", "cf02_sp100", "cf100_sp100_rel", "cf075_sp100"]]
print(INCLUDE_SETTINGS)
STYLE_MAP = build_style_map(INCLUDE_SETTINGS)
LABEL_MAP = build_label_map(INCLUDE_SETTINGS)
# ---------------------

if not INCLUDE_SETTINGS:
    print(f"warning: no settings discovered for {EXP_NAME}")
elif not EXPERIMENT_DIR.is_dir() or not any(EXPERIMENT_DIR.iterdir()):
    print(f"warning: no results in {EXPERIMENT_DIR}; skipping comparison")
else:
    try:
        print("4")
        run_full_comparison(EXPERIMENT_DIR, baseline=BASELINE,
                            pd_convergence=PD_CONVERGENCE, style_map=STYLE_MAP,
                            label_map=LABEL_MAP, min_buses=73,
                            include_settings=INCLUDE_SETTINGS,
                            min_time=MIN_TIME)
    except FileNotFoundError as e:
        print(f'warning: {e}')
