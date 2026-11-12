#!/usr/bin/env python3
"""Compare settings from Exp. 6 Stage 1 - 1-tal1 (pool size x weight scheme)."""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]            # 1-baseline (this repo)
SHARED = PROJECT_ROOT.parent
STAGE_REPO = SHARED / "6-seq_warmstarting"      # read-only data source
sys.path.insert(0, str(SCRIPT_DIR.parent))  # visu/
sys.path.insert(0, str(SCRIPT_DIR))         # for _seqwarmstart_common

from _seqwarmstart_common import run_experiment  # noqa: E402

EXP_NAME = "1-tal1"
run_experiment(
    experiment_dir=STAGE_REPO / "check/results" / f"6-seq_warmstart--{EXP_NAME}",
    settings_dir=STAGE_REPO / "settings/6-seq_warmstart" / EXP_NAME,
    output_dir=PROJECT_ROOT / "check/results/paper/6-seq_warmstarting/full" / EXP_NAME,
    pd_convergence=True,   # produces aggregated relative-primal plot for the paper
    pd_top_n=5,
)
