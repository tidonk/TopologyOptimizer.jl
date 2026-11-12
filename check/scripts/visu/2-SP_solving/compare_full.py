#!/usr/bin/env python3
"""Stage 2 (SP solving) — FULL-run comparison (GRB vs SLE vs SLE-U).

Reproducible replacement for the on-the-fly comparison whose output Oli moved
into ``across_all/2-SP_solving/full/``.  Averages the three full-run repetitions
(``2-SP_solving1``/``2``/``3``, auto-discovered) per instance.

Zero-touch: lives in and runs from **1-baseline**'s visu tree (which carries the
averaging-capable toolchain), *reads* the 2-SP_solving repo, and writes into the
1-baseline staging dir ``check/results/paper/2-SP_solving/full/``.
``build_paper_resources.sh`` then harvests it into ``across_all/2-SP_solving/full/``.

Importing this module has no side effects: the config below (notably ``SOURCES``)
is read by ``build_cut_err_exclusions.py`` to scope the CUT_ERROR scan to exactly
the settings we actually report.

Run:  python3 compare_full.py
"""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]      # visu/2-SP_solving -> 1-baseline
SHARED = PROJECT_ROOT.parent              # experiment root
sys.path.insert(0, str(SCRIPT_DIR.parent))

from reader import discover_reps, load_excluded_instances, CUT_ERR_EXCLUSIONS_FILE
from tables import run_full_comparison

STAGE = "2-SP_solving"
FULL_BASE = SHARED / STAGE / "check/results/2-SP_solving"   # -> 2-SP_solving1/2/3
OUTPUT_DIR = PROJECT_ROOT / "check/results/paper" / STAGE / "full"
BASELINE = "benders_FLLP"

LABEL_MAP = {
    "benders_FLLP":  r"\texttt{GRB}",
    "benders_LEQ":   r"\texttt{SLE}",
    "benders_LEQRF": r"\texttt{SLE-U}",
}
INCLUDE_SETTINGS = ["benders_FLLP", "benders_LEQ", "benders_LEQRF"]

REPS = discover_reps(FULL_BASE)
# (source dir, settings reported from it) — consumed by build_cut_err_exclusions.py
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
