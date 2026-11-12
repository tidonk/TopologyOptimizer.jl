#!/usr/bin/env python3
"""``FINAL/full`` — the paper's headline progression table.

Oli built this from ``5-partial_decomp/check/results/X-final-full``, which is a
**compiled-mean dir** (1 CSV per setting, §3.1 of EXPERIMENTS.md) — so it froze
the aggregation.  Its seven settings were identified by fingerprint-matching the
``Solved``/``TL``/``Sum`` triples of his ``comparison_short_all.tex`` against
every full-run stack on disk: all seven matched X-final-full at 0.0000 relative
error.  Here each is instead re-aggregated from its own real 3-repetition
source, exactly as every other leaf does.

The table carries two extra hand-authored columns encoding the paper's
progression DAG: ``\\#`` (the step number) and ``Builds on``.  They are expressed
as label parameters (``n=`` / ``on=``) and rendered by the standard
``param_columns`` machinery with ``param_keep_label=True`` (§7).

Metric is ``time_solve`` (matching Oli's ``performance_profile_time_solve.pdf``),
even though the X-final-full dir also holds ``time_with_warmstart`` artifacts.

Run:  python3 compare_full.py
"""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]          # visu/FINAL -> 1-baseline
SHARED = PROJECT_ROOT.parent
sys.path.insert(0, str(SCRIPT_DIR.parent))

from reader import discover_reps, load_excluded_instances, CUT_ERR_EXCLUSIONS_FILE
from tables import run_full_comparison

PAPER = PROJECT_ROOT / "check/results/paper/FINAL/full"
EXCL = load_excluded_instances(CUT_ERR_EXCLUSIONS_FILE)   # full test set

# (setting, label, step number, builds on)
# NB (2026-07-29): the "SP Selection" step now uses the RANDOM scorer rnd_abs_2
# (the full-run winner, Geo 68.2 — see paper/3-SP_selection/full/all), replacing
# the old diversity setting div_decay_2 (Geo 88.0).  Steps 6/7 still reference
# diversity variants (div_decay_5_005, base_only_div) — see EXPERIMENTS_INDEX.md
# "Pipeline settings audit" for the remaining repoints to non-div sources.
STEPS = [
    ("compact",         r"\texttt{Compact}",          1, "-"),
    ("benders_FLLP",    r"\texttt{Benders}",          2, "1"),
    ("benders_LEQRF",   r"\texttt{Benders SLE-U}",    3, "2"),
    ("rnd_abs_2",       r"\texttt{SP Selection}",     4, "3"),
    ("cutfiltering",    r"\texttt{Cut Filtering}",    5, "3"),
    ("cf05_sp005",      r"\texttt{SP+CF combined}",   6, "4+5"),
    ("base_only",       r"\texttt{Partial Decomp.}",  7, "6"),
]

LABEL_MAP = {s: f"{lbl}, n={num}, on={on}" for s, lbl, num, on in STEPS}
PARAM_COLS = [(r"\#", "n"), ("Builds on", "on")]
SETTINGS = [s for s, _l, _n, _o in STEPS]

R = SHARED
_SRC = {
    "compact":         discover_reps(R / "1-baseline/check/results/1-baseline-compact"),
    "benders_FLLP":    discover_reps(R / "2-SP_solving/check/results/2-SP_solving"),
    "benders_LEQRF":   discover_reps(R / "2-SP_solving/check/results/2-SP_solving"),
    # SP Selection = RANDOM (rnd_abs_2), the full-run winner — same SP-eval family
    # as the other scorers (7-score-full-eval, cut filtering off at this step),
    # replacing the diversity setting div_decay_2.
    "rnd_abs_2":       discover_reps(R / "3-SP_selection/check/results/7-score-full-eval"),
    "cutfiltering":    discover_reps(R / "4-cut_filtering/check/results/4-cut_filtering"),
    # Step 6 = SP+CF combined (cf05_sp005): random SP (rel-5%) + cut filtering
    # (:VIOLATION, keep 5), no partial decomp — the full run of the grid winner.
    "cf05_sp005":      discover_reps(R / "4-cut_filtering/check/results/9-cutfiltering-full"),
    # Partial Decomp = the NON-diversity partial-decomp run: SP selection :RND +
    # cut filtering ON (base_only.toml), replacing base_only_div (SP diversity +
    # CF off).  Partial-decomp contingency selection may be non-random; SP must not.
    "base_only":       discover_reps(R / "5-partial_decomp/check/results/2-full-eval_num_1"),
}


def _dirs(v):
    return v[0] if isinstance(v, tuple) else v


SOURCES = [(d, [s]) for s in SETTINGS for d in _dirs(_SRC[s])]


def main() -> None:
    spec = {s: _SRC[s] for s in SETTINGS}
    print(f"[compare_full] FINAL: re-aggregating {len(spec)} settings from their own runs")
    run_full_comparison(
        spec, baseline="benders_LEQRF", output_dir=PAPER,
        exclude_instances=EXCL, pd_convergence=False, min_buses=73,
        label_map=LABEL_MAP, include_settings=SETTINGS,
        param_columns=PARAM_COLS, param_keep_label=True,
    )


if __name__ == "__main__":
    main()
