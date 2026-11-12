#!/usr/bin/env python3
"""Stage 3 grid — single-feature scorers vs random (`2-score-individual`).

FINDINGS §3a: each hand-crafted scorer (Viol / Rel / Share, decay & no-decay) is
compared against the closest random setting (`rnd_abs_2`) to show that scoring
does not beat random.  Grid test set, single dir.  Cut errors are *reported*
(CutErr column) but not dropped — matching the published table.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))     # visu/
from _gridcompare import run_compare, SHARED, PAPER               # noqa: E402

NAME     = "grid_score_individual"
SOURCE   = SHARED / "3-SP_selection/check/results/2-score-individual"
OUTPUT   = PAPER / "3-SP_selection/grid/2-score-individual"
BASELINE = "rnd_abs_2"

LABELS = {
    "rnd_abs_2":     r"\texttt{Random}",
    "viol_nodecay":  r"\texttt{Viol (no decay)}",
    "viol_decay":    r"\texttt{Viol (decay)}",
    "rel_nodecay":   r"\texttt{Rel (no decay)}",
    "rel_decay":     r"\texttt{Rel (decay)}",
    "share_nodecay": r"\texttt{Share (no decay)}",
    "share_decay":   r"\texttt{Share (decay)}",
}

run_compare(name=NAME, source=SOURCE, output=OUTPUT, baseline=BASELINE,
            labels=LABELS, cut="none", min_buses=73, rank=False)
