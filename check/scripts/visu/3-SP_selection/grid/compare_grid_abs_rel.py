#!/usr/bin/env python3
"""Stage 3 grid — combined abs × rel SP-selection (`1-n_hits_abs_rel`).

STAB=Inf (nostab), random SP selection. Reads the across-folder of symlinks:
20 combined cells (`1-n_hits_abs_rel`) + 9 pure-axis baselines (`1-n_hits_nostab`,
NOT `1-n_hits` which is STAB=100) + the `CF5+SP5%` reference (`4-cut_filtering`).
Cut errors penalized-not-dropped, best among ≤2-cut-err cells (§12.3d).
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))     # visu/
from _gridcompare import run_compare, SHARED, PAPER               # noqa: E402

NAME     = "grid_abs_rel"
SOURCE   = SHARED / "3-SP_selection/check/results/1-n_hits_abs_rel-across"
OUTPUT   = PAPER / "3-SP_selection/grid/1-n_hits_abs_rel"
BASELINE = "cf05_sp005"

# Label convention: SP{abs} (absolute cap), SP{rel}% (relative cap),
# SP{abs}SP{rel}% (both) — e.g. abs 5 + rel 5% -> SP5SP5%.
_ABS, _REL = (1, 2, 5, 10), ("001", "005", "010", "025", "050")
LABELS = {
    **{f"rnd_abs_{a}_rel_{r}": rf"\texttt{{SP{a}SP{int(r)}\%}}"
       for a in _ABS for r in _REL},
    **{f"rnd_abs_{a}": rf"\texttt{{SP{a}}}" for a in _ABS},
    **{f"rnd_rel_{r}": rf"\texttt{{SP{int(r)}\%}}" for r in _REL},
    "rnd_rel_100": r"\texttt{SP100\%}",     # abs off, rel 100% = solve all (the no-cap corner)
    "cf05_sp005": r"\texttt{CF5+SP5\%}",
}

run_compare(name=NAME, source=SOURCE, output=OUTPUT, baseline=BASELINE, labels=LABELS)
