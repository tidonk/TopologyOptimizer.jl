#!/usr/bin/env python3
"""Stage 82 (fractional separation) — GRID comparisons: all / grb / les.

Source ``82-fractional_sepa`` holds exactly the 14 settings of Oli's
``grid/all`` table.  The three leaves are subsets of it:

* ``all`` — all 14 settings (no baseline / Wilcoxon in Oli's table)
* ``grb`` — the 9 full-LP (``_fulllp`` = ``_GRB``) settings, baseline
  ``base_only_fulllp``
* ``les`` — the 5 non-full-LP settings, baseline ``base_only``

Naming: ``_fulllp`` -> ``_GRB``; ``_1``/``_2``/``_5`` = ``FRACTIONAL_SEPA_ITERLIM``,
**un-suffixed = 10** (verified against the setting TOMLs).

Grid runs are single dirs and use **per-experiment** cut-error scoping (this
grid's own proven CUT_ERRORs only) — never the full-set list.

Run:  python3 compare_grid.py
"""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]        # visu/82-fractional_sepa -> 1-baseline
SHARED = PROJECT_ROOT.parent
sys.path.insert(0, str(SCRIPT_DIR.parent))

from reader import grid_cut_err_instances
from tables import run_full_comparison

STAGE = "82-fractional_sepa"
SRC = SHARED / STAGE / "check/results/82-fractional_sepa"
PAPER = PROJECT_ROOT / "check/results/paper" / STAGE / "grid"

# setting -> Oli's label.  Un-suffixed = ITERLIM 10.
LABEL_MAP = {
    "base_only":            r"\texttt{baseline}",
    "root_node_1":          r"\texttt{root\_1}",
    "root_node":            r"\texttt{root\_10}",
    "all_nodes_1":          r"\texttt{all\_1}",
    "all_nodes":            r"\texttt{all\_10}",
    "base_only_fulllp":     r"\texttt{baseline\_GRB}",
    "root_node_fulllp_1":   r"\texttt{root\_1\_GRB}",
    "root_node_fulllp_2":   r"\texttt{root\_2\_GRB}",
    "root_node_fulllp_5":   r"\texttt{root\_5\_GRB}",
    "root_node_fulllp":     r"\texttt{root\_10\_GRB}",
    "all_nodes_fulllp_1":   r"\texttt{all\_1\_GRB}",
    "all_nodes_fulllp_2":   r"\texttt{all\_2\_GRB}",
    "all_nodes_fulllp_5":   r"\texttt{all\_5\_GRB}",
    "all_nodes_fulllp":     r"\texttt{all\_10\_GRB}",
}

LES = ["base_only", "root_node_1", "root_node", "all_nodes_1", "all_nodes"]
GRB = ["base_only_fulllp",
       "root_node_fulllp_1", "root_node_fulllp_2", "root_node_fulllp_5", "root_node_fulllp",
       "all_nodes_fulllp_1", "all_nodes_fulllp_2", "all_nodes_fulllp_5", "all_nodes_fulllp"]

# (leaf, baseline, settings in Oli's column order)
LEAVES = [
    ("all", "", LES + GRB),
    ("grb", "base_only_fulllp", GRB),
    ("les", "base_only", LES),
]

# (dir, settings reported from it) — consumed by build_cut_err_exclusions.py
SOURCES = [(SRC, LES + GRB)]


def main() -> None:
    if not SRC.is_dir():
        print(f"warning: {SRC} missing; skipping")
        return
    for name, baseline, settings in LEAVES:
        excl = grid_cut_err_instances(SRC, settings)
        print(f"[compare_grid] {name}: {len(settings)} settings, baseline "
              f"{baseline or '(none)'}, {len(excl)} own cut-err instance(s) excluded")
        run_full_comparison(
            SRC, baseline=baseline, output_dir=PAPER / name,
            exclude_instances=excl, pd_convergence=False, min_buses=73,
            label_map=LABEL_MAP, include_settings=settings,
        )


if __name__ == "__main__":
    main()
