#!/usr/bin/env python3
"""Stage 85 (Benders primal heuristic) — GRID comparison.

Single two-setting experiment ``benders-primal-heur``
(``SUBMIT_HEUR_SOLUTION`` on vs off); config verbatim from
``85-benders-primal-heur/.../compare_85-benders-primal-heur.py``.

Carry-forward: the ``base_only`` baseline is a **symlink** into
``5-partial_decomp/.../5-partial_decomp_num_grid/base_only``, so it resolves
without any extra source mapping — but note that it therefore shares the stage-5
grid baseline, not a stage-85 run of its own.

Grid runs are single dirs and use **per-experiment** cut-error scoping.

Run:  python3 compare_grid.py
"""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]      # visu/85-benders-primal-heur -> 1-baseline
SHARED = PROJECT_ROOT.parent
sys.path.insert(0, str(SCRIPT_DIR.parent))

from reader import grid_cut_err_instances
from tables import run_full_comparison

STAGE = "85-benders-primal-heur"
SRC = SHARED / STAGE / "check/results/benders-primal-heur"
PAPER = PROJECT_ROOT / "check/results/paper" / STAGE / "grid" / "1-benders-primal-heur"

LABEL_MAP = {
    "base_only":           r"\texttt{baseline}",
    "benders-primal-heur": r"\texttt{heur submit}",
}
SETTINGS = ["base_only", "benders-primal-heur"]
BASELINE = "base_only"

SOURCES = [(SRC, SETTINGS)]


def main() -> None:
    if not SRC.is_dir():
        print(f"warning: {SRC} missing; skipping")
        return
    excl = grid_cut_err_instances(SRC, SETTINGS)
    print(f"[compare_grid] 1-benders-primal-heur: {len(SETTINGS)} settings, "
          f"baseline {BASELINE}, {len(excl)} own cut-err instance(s) excluded")
    run_full_comparison(
        SRC, baseline=BASELINE, output_dir=PAPER,
        exclude_instances=excl, pd_convergence=False, min_buses=0,
        label_map=LABEL_MAP, include_settings=SETTINGS,
    )


if __name__ == "__main__":
    main()
