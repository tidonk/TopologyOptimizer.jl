#!/usr/bin/env python3
"""Stage 3 (SP selection) — FULL-run comparisons: single / multi / all / pb.

Reproducible replacement for Oli's on-the-fly runs (whose config lives in
``3-SP_selection/.../compare.py``: the LABEL_MAP + the commented EXPERIMENT_DIRs).

These are inherently **single-run**: the numbered rep dirs
(``7-score-full-eval1/2/3``, ``8-pb-full1/2/3``) are missing the baseline
settings (``benders_LEQRF`` / ``div_decay_2`` / ``rnd_abs_2``), so the complete
comparison only exists in the un-suffixed dirs.  Zero-touch: runs from
1-baseline, reads the 3-SP_selection repo, writes staging under 1-baseline.

Run:  python3 compare_full.py
"""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[4]          # visu/3-SP_selection -> 1-baseline
SHARED = PROJECT_ROOT.parent
sys.path.insert(0, str(SCRIPT_DIR.parents[1]))

from reader import discover_reps, load_excluded_instances, CUT_ERR_EXCLUSIONS_FILE
from tables import run_full_comparison

STAGE = "3-SP_selection"
RES = SHARED / STAGE / "check/results"
PAPER = PROJECT_ROOT / "check/results/paper" / STAGE / "full"
EXCL = load_excluded_instances(CUT_ERR_EXCLUSIONS_FILE)   # full test set

# LABEL_MAP verbatim from 3-SP_selection/compare.py
LABEL_MAP = {
    "benders_LEQRF":          r"\texttt{All}",
    "rnd_abs_2":              r"\texttt{Rnd100}",
    "rel_decay":              r"\texttt{Rel100}",
    "share_decay":            r"\texttt{Shr100}",
    "viol_decay":             r"\texttt{Vio100}",
    "div_decay_2":            r"\texttt{Div100}",
    "hybrid_rel_rnd_50_50":   r"\texttt{Rel50/Rnd50}",
    "hybrid_share_rnd_25_75": r"\texttt{Shr25/Rnd75}",
    "hybrid_share_rnd_50_50": r"\texttt{Shr50/Rnd50}",
    "rnd_pb":                 r"\texttt{pb}",
    "rnd_abs_5":              r"\texttt{abs5}",
    "rnd_abs_10":             r"\texttt{abs10}",
    "rnd_abs_2_pb":           r"\texttt{abs2pb}",
    "rnd_abs_5_pb":           r"\texttt{abs5pb}",
    "rnd_abs_10_pb":          r"\texttt{abs10pb}",
}

# (leaf, source dir, baseline, settings in Oli's column order)
# Per-setting sources.  The un-suffixed dirs (``7-score-full-eval[-single|-multi]``,
# ``8-pb-full``) are NOT runs: each holds a single *pre-computed mean CSV* per
# setting.  The real repetitions live in the numbered dirs — and the two
# "baseline" settings come from other experiments entirely (``benders_LEQRF``
# from stage 2, ``div_decay_2`` from the diversity full run), which is why the
# numbered dirs appear to be "missing" them.
_SCORE = discover_reps(RES / "7-score-full-eval")                  # 3 reps
_PB    = discover_reps(RES / "8-pb-full")                          # 3 reps
# Diversity is now evaluated COLD (warmstart off): the cold full run
# `2-score-diversity-full-cold{1,2,3}` replaces the warm `2-score-diversity-full`.
# discover_reps picks up the three reps once the cluster jobs land.
_DIV   = discover_reps(RES / "2-score-diversity-full-cold")        # 3 reps (cold)
_LEQRF = discover_reps(SHARED / "2-SP_solving/check/results/2-SP_solving")

_SRC = {
    "benders_LEQRF":          _LEQRF,   # carry-forward from stage 2
    "div_decay_2":            _DIV,     # carry-forward from the diversity run
    "rnd_abs_2":              _SCORE,
    "rel_decay":              _SCORE,
    "share_decay":            _SCORE,
    "viol_decay":             _SCORE,
    "hybrid_rel_rnd_50_50":   _SCORE,
    "hybrid_share_rnd_25_75": _SCORE,
    "hybrid_share_rnd_50_50": _SCORE,
    "rnd_pb":                 _PB,
    "rnd_abs_5":              _PB,
    "rnd_abs_10":             _PB,
    "rnd_abs_2_pb":           _PB,
    "rnd_abs_5_pb":           _PB,
    "rnd_abs_10_pb":          _PB,
}

# (leaf, baseline, settings in Oli's column order)
LEAVES = [
    ("single", "benders_LEQRF",
     ["benders_LEQRF", "rnd_abs_2", "rel_decay", "share_decay", "viol_decay", "div_decay_2"]),
    ("multi", "benders_LEQRF",
     ["benders_LEQRF", "rnd_abs_2", "hybrid_rel_rnd_50_50",
      "hybrid_share_rnd_25_75", "hybrid_share_rnd_50_50", "div_decay_2"]),
    # Merged all-settings table. Fair again now that div_decay_2 (Div100) is run
    # COLD (warmstart off) like the rest — no setting gets a hidden warmstart edge.
    ("all", "benders_LEQRF",
     ["benders_LEQRF", "rel_decay", "share_decay", "viol_decay", "rnd_abs_2",
      "div_decay_2", "hybrid_rel_rnd_50_50", "hybrid_share_rnd_25_75", "hybrid_share_rnd_50_50"]),
    ("pb", "rnd_abs_2",
     ["rnd_pb", "rnd_abs_2", "rnd_abs_5", "rnd_abs_10",
      "rnd_abs_2_pb", "rnd_abs_5_pb", "rnd_abs_10_pb"]),
]

# (dir, settings reported from it) — consumed by build_cut_err_exclusions.py
SOURCES = [(d, [s]) for _n, _b, settings in LEAVES for s in settings for d in _SRC[s]]


# Diversity is now run COLD (warmstart off) and no other setting warmstarts, so
# ``time_with_warmstart`` would equal ``time_solve`` everywhere — we emit only the
# plain ``time_solve`` table per leaf.
def main() -> None:
    for name, baseline, settings in LEAVES:
        spec = {s: _SRC[s] for s in settings}
        print(f"[compare_full] {name}: {len(spec)} settings, target time_solve")
        run_full_comparison(
            spec, baseline=baseline, output_dir=PAPER / name,
            target="time_solve", exclude_instances=EXCL,
            pd_convergence=False, min_buses=73,
            label_map=LABEL_MAP, include_settings=settings,
        )


if __name__ == "__main__":
    main()
