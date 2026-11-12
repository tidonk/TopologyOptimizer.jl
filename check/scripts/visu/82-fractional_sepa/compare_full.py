#!/usr/bin/env python3
"""Stage 82 (fractional separation) — FULL-run comparison.

Naming decoded from the setting TOMLs:

* ``_fulllp``          -> Oli's ``_GRB`` label (full-LP subproblem solve)
* ``_1``/``_2``/``_5`` -> ``FRACTIONAL_SEPA_ITERLIM``; **un-suffixed = 10**
  (verified: ``root_node_fulllp.toml`` has ITERLIM=10, ``..._5.toml`` has 5)

``2-fractional_sepa-full{1,2,3}`` holds exactly the three settings of Oli's
table and has 3 repetitions, so it is averaged.

Run:  python3 compare_full.py
"""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]        # visu/82-fractional_sepa -> 1-baseline
SHARED = PROJECT_ROOT.parent
sys.path.insert(0, str(SCRIPT_DIR.parent))

from reader import discover_reps, load_excluded_instances, CUT_ERR_EXCLUSIONS_FILE
from tables import run_full_comparison

STAGE = "82-fractional_sepa"
FULL_BASE = SHARED / STAGE / "check/results/2-fractional_sepa-full"   # -> 1/2/3
OUTPUT_DIR = PROJECT_ROOT / "check/results/paper" / STAGE / "full"
BASELINE = ""          # Oli's table has no (baseline) row / Wilcoxon columns

LABEL_MAP = {
    "base_only_fulllp":   r"\texttt{baseline\_GRB}",
    "root_node_fulllp_1": r"\texttt{root\_1\_GRB}",
    "all_nodes_fulllp_1": r"\texttt{all\_1\_GRB}",
}
INCLUDE_SETTINGS = ["base_only_fulllp", "root_node_fulllp_1", "all_nodes_fulllp_1"]

REPS = discover_reps(FULL_BASE)
# (dir, settings reported from it) — consumed by build_cut_err_exclusions.py
SOURCES = [(d, INCLUDE_SETTINGS) for d in REPS]


def main() -> None:
    if not REPS:
        print(f"warning: no repetitions found for {FULL_BASE}; skipping")
        return
    print(f"[compare_full] averaging {len(REPS)} reps: {[p.name for p in REPS]}")
    run_full_comparison(
        REPS, baseline=BASELINE, output_dir=OUTPUT_DIR,
        exclude_instances=load_excluded_instances(CUT_ERR_EXCLUSIONS_FILE),
        pd_convergence=False, min_buses=73,
        label_map=LABEL_MAP, include_settings=INCLUDE_SETTINGS,
    )


if __name__ == "__main__":
    main()
