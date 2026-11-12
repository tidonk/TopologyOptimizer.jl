#!/usr/bin/env python3
"""Stage 3 grid — cold diversity decay sweep (`2-score-diversity-cold-grid`).

Settles §10.7/§11.0: diversity run **cold** (warmstart off) vs cold random, on the
31-inst grid. Verdict: random wins — diversity is not competitive cold. Labels use
SP-notation with the method: `SP{n} (Div)` for the decay-diversity settings,
`SP{n} (Rnd)` for the random baselines. Cut errors penalized, ≤2-cut-err gate.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))     # visu/
from _gridcompare import run_compare, SHARED, PAPER               # noqa: E402

NAME     = "grid_cold_diversity"
SOURCE   = SHARED / "3-SP_selection/check/results/2-score-diversity-cold-grid"
OUTPUT   = PAPER / "3-SP_selection/grid/2-score-diversity-cold-grid"
BASELINE = "rnd_abs_2"

LABELS = {
    "div_decay_1":  r"\texttt{SP1 (Div)}",
    "div_decay_2":  r"\texttt{SP2 (Div)}",
    "div_decay_5":  r"\texttt{SP5 (Div)}",
    "div_decay_10": r"\texttt{SP10 (Div)}",
    "rnd_abs_2":    r"\texttt{SP2 (Rnd)}",
    "rnd_abs_5":    r"\texttt{SP5 (Rnd)}",
}

run_compare(name=NAME, source=SOURCE, output=OUTPUT, baseline=BASELINE, labels=LABELS)
