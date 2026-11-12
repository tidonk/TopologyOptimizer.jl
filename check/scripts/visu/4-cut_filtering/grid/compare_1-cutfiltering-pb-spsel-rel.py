#!/usr/bin/env python3
"""Compare settings from Exp. 4 Stage 1 - 1-cutfiltering-pb-spsel-rel.

Cut filter: PB-aware (fixed). SP selection: relative N_HITS_REL sweep.
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
EXP_NAME = "1-cutfiltering-pb-spsel-rel"
EXPERIMENT_DIR = PROJECT_ROOT / "check/results" / EXP_NAME
SETTINGS_DIR = PROJECT_ROOT / "settings/4-cut_filtering" / EXP_NAME
BASELINE = "sp100"
PD_CONVERGENCE = False
MIN_TIME = 0

INCLUDE_SETTINGS = discover_settings(EXPERIMENT_DIR, SETTINGS_DIR)
STYLE_MAP = build_style_map(INCLUDE_SETTINGS)
LABEL_MAP = build_label_map(INCLUDE_SETTINGS)
# ---------------------

if not INCLUDE_SETTINGS:
    print(f"warning: no settings discovered for {EXP_NAME}")
elif not EXPERIMENT_DIR.is_dir() or not any(EXPERIMENT_DIR.iterdir()):
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
