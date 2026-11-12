#!/usr/bin/env python3
"""Stage 5 (partial decomposition) — FULL-run comparisons: num_1 / scores.

Config derived from ``5-partial_decomp/.../compare.py`` (EXP_NAME options,
BASELINE ``div_decay_5_005``, ``target="time_with_warmstart"`` + the
``_add_time_with_warmstart`` transform, LABEL_MAP).

Both leaves use the ``_div_full`` variants — i.e. the partial-decomposition
settings re-run on top of the *diversity* SP-selection stack that stages 3/4
settled on (``div_decay_5_005`` rather than the older random ``cf05_sp005``).
Confirmed empirically: ``2-full-eval-scores_div_full`` reproduces Oli's
full/scores numbers exactly, the plain variant does not.

Unlike stages 3/4 these sources DO have 3 repetitions, so they are averaged.

Zero-touch: runs from 1-baseline, reads the 5-partial_decomp repo, writes staging
under 1-baseline.

Run:  python3 compare_full.py
"""

import sys
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]        # visu/5-partial_decomp -> 1-baseline
SHARED = PROJECT_ROOT.parent
sys.path.insert(0, str(SCRIPT_DIR.parent))

from reader import discover_reps, load_excluded_instances, CUT_ERR_EXCLUSIONS_FILE
from tables import run_full_comparison

STAGE = "5-partial_decomp"
RES = SHARED / STAGE / "check/results"
PAPER = PROJECT_ROOT / "check/results/paper" / STAGE / "full"
EXCL = load_excluded_instances(CUT_ERR_EXCLUSIONS_FILE)   # FULL test set

LABEL_MAP = {
    # num_1 leaf — plain label = the *_no_base* variant, "+ Base" = with base
    "div_decay_5_005":      r"\texttt{No PD}",
    "base_only_div":        r"\texttt{Base only}",
    "eta_num_no_base_1":    r"\texttt{Viol1}",
    "rnd_num_no_base_1":    r"\texttt{Rnd1}",
    "eta_num_1_div":        r"\texttt{Viol1 + Base}",
    "rnd_num_1_div":        r"\texttt{Rnd1 + Base}",
    # scores leaf — now the NON-div settings (correct SP/CF = :RND rel-5% + CF).
    # cprox/sens/eta keep warmstart (required); div_cos/div_euc/rnd are cold.
    "eta_num_5":            r"\texttt{Viol5}",
    "rnd_num_5":            r"\texttt{Rnd5}",
    "cprox_num_5":          r"\texttt{CProx5}",
    "div_cos_num_5":        r"\texttt{Div-Cos5}",
    "div_euc_num_5":        r"\texttt{Div-Euc5}",
    "sens_num_5":           r"\texttt{Sens5}",
    # all leaf — base_only (N_MAX_CONTINGENCIES=0) as the reference row.
    "base_only":            r"\texttt{Base only}",
}

# Per-setting sources.  The un-suffixed ``2-full-eval_num_1_div_full`` dir was
# hand-assembled by copying pre-computed *mean* CSVs in from elsewhere: its rep
# dirs hold only base_only_div / eta_num_1_div / rnd_num_1_div, so reading it (or
# averaging its reps) silently loses the other three settings — including the
# ``div_decay_5_005`` baseline.  We therefore map every setting to its **own real
# run dirs** and re-aggregate here; nothing is taken from a copied mean CSV.
_DIV_FULL = discover_reps(RES / "2-full-eval_num_1_div_full")     # 3 reps (OBSOLETE — see num_1 note)
_NO_BASE = discover_reps(RES / "2-full-eval_num_1_no_base_div")   # 3 reps (OBSOLETE — no non-div exists)
# scores leaf, corrected SP/CF base (:RND rel-5% + CF), split by warmstart regime:
_SCORES_WARM = discover_reps(RES / "2-full-eval-scores")          # cprox/sens/eta (warmstart REQUIRED)
_SCORES_COLD = discover_reps(RES / "2-full-eval-scores-cold")     # div_cos/div_euc/rnd (warmstart OFF; lands after resubmit)
# div_decay_5_005 (No PD) has a proper 3-run source in 3-SP_selection, where the
# runs are filed under the folder name ``div_decay_5`` (verified: the mean Oli
# copied in as div_decay_5_005 reproduces these reps exactly, and the CSVs'
# internal settings column reads div_decay_5).
_NO_PD = discover_reps(
    SHARED / "3-SP_selection/check/results/2-score-diversity-full-hybrid")

NUM_1_SOURCES = {
    "eta_num_no_base_1": _NO_BASE,
    "rnd_num_no_base_1": _NO_BASE,
    "div_decay_5_005":   (_NO_PD, "div_decay_5"),   # renamed folder
    "base_only_div":     _DIV_FULL,
    "eta_num_1_div":     _DIV_FULL,
    "rnd_num_1_div":     _DIV_FULL,
}
SCORES_SOURCES = {
    "cprox_num_5":   _SCORES_WARM, "sens_num_5": _SCORES_WARM, "eta_num_5": _SCORES_WARM,
    "div_cos_num_5": _SCORES_COLD, "div_euc_num_5": _SCORES_COLD, "rnd_num_5": _SCORES_COLD,
}
# base_only (N_MAX_CONTINGENCIES=0, cold, correct :RND rel-5% + CF).  Its own
# real 3-rep dirs are the numbered 2-full-eval_num_1{1,2,3}; the un-suffixed
# 2-full-eval_num_1 holds only copied mean CSVs, so we read the reps directly.
_BASE = discover_reps(RES / "2-full-eval_num_1")                  # base_only reps 11/12/13
# The ``all`` leaf places base_only next to the six num_5 scores so the paper
# can show, in one table, that partial decomposition never beats base_only:
#   • cold div/rnd land byte-identical to base_only (n_solved_subs_total,
#     n_callbacks, n_benders_cuts all equal) — PD needs the warmstart to select
#     the contingencies to embed, so cold ⇒ PD is silently inactive (N=0);
#   • warm cprox/sens/eta are the only runs where PD truly activates, and they
#     are ~2x slower than base_only.
ALL_SOURCES = {"base_only": _BASE, **SCORES_SOURCES}

# (leaf, per-setting source map, baseline, settings in Oli's column order)
LEAVES = [
    # ⚠️ OBSOLETE: the num_1 leaf still reads the diversity `_div` stack
    # (SP=:SCORE, CF=off) because the non-div `no_base` full runs do not exist.
    # It must be re-run (non-div base_only/eta/rnd + no_base) before use — see
    # EXPERIMENTS_INDEX.md.  Kept only so the script still produces the table.
    ("num_1", NUM_1_SOURCES, "div_decay_5_005",
     ["div_decay_5_005", "eta_num_no_base_1", "rnd_num_no_base_1",
      "base_only_div", "eta_num_1_div", "rnd_num_1_div"]),
    # scores leaf — corrected: non-div settings, warm cprox/sens/eta + cold div/rnd.
    ("scores", SCORES_SOURCES, "",                               # no baseline
     ["eta_num_5", "rnd_num_5", "cprox_num_5",
      "div_cos_num_5", "div_euc_num_5", "sens_num_5"]),
    # all leaf — base_only (reference) vs every num_5 score in one table.
    # Wilcoxon vs base_only makes the finding explicit: cold div/rnd are not
    # significantly different (same computation, PD inactive), warm scorers lose.
    ("all", ALL_SOURCES, "base_only",
     ["base_only", "rnd_num_5", "div_cos_num_5", "div_euc_num_5",
      "eta_num_5", "cprox_num_5", "sens_num_5"]),
]

def _dirs_folder(value, setting):
    """Normalise a source-map value to ``(dirs, source_folder)``."""
    if isinstance(value, tuple) and len(value) == 2 and isinstance(value[1], str):
        return value[0], value[1]
    return value, setting


# (dir, settings reported from it) — consumed by build_cut_err_exclusions.py.
# Filter on the *folder* name, which is what read_experiment labels rows with.
SOURCES = [(d, [folder])
           for _n, smap, _b, _order in LEAVES
           for s, v in smap.items()
           for dirs, folder in [_dirs_folder(v, s)]
           for d in dirs]


def main() -> None:
    for name, smap, baseline, settings in LEAVES:
        print(f"[compare_full] {name}: re-aggregating {len(smap)} settings from their own runs")
        run_full_comparison(
            smap, baseline=baseline, output_dir=PAPER / name,
            exclude_instances=EXCL, pd_convergence=False, min_buses=73,
            label_map=LABEL_MAP, include_settings=settings,
            target="time_with_warmstart",   # derived centrally in run_full_comparison
        )


if __name__ == "__main__":
    main()
