#!/usr/bin/env python3
"""Stage 9 (nested Benders) — FULL-run comparison.

Config from ``9-nested/.../9-nested/compare.py``, but with its source replaced:
that script points at ``3-nested-bpfi-full``, which is a **compiled-mean dir**
(1 CSV per setting).  The real repetitions are ``3-nested-bpfi-full{1,2,3}``.

The ``base_only_div`` baseline is a **carry-forward** and has no repetition in
this repo.  Its origin is named by the ``instance_file`` column of the copied
mean CSV —
``5-partial_decomp/settings/5-partial_decomp/2-full-eval_num_1_div/base_only_div.toml``
— whose three reps are ``2-full-eval_num_1_div_full{1,2,3}`` (the same source
stage 5 uses for its own ``base_only_div`` row).

Note that ``compare.py`` also lists six settings that were never run at full
scale (``bpfi_cp0``, ``bpfi_gmp``, ``bpfi_ncs``, ``bpfi_no_pbs``, ``bpfi_no_ph``,
``bpfi_spc``); they exist only in the ``2-nested-bpfi-grid`` grid search.

**Cut-error scoping.** Stage 9 is deliberately *opted out* of the global
full-set exclusion list (``build_cut_err_exclusions.EXCLUDED_STAGES``): its BPFI
settings have real, surviving cut errors, and folding them into the shared list
grew it 31 -> 42, shrinking the full test set to 111 instances for *every* other
stage.  Stage 9 is not part of the paper's across_all, so — exactly as for the
grid searches — its defects stay local: exclusions are computed here from its
own ``SOURCES``.

Run:  python3 compare_full.py
"""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]          # visu/9-nested -> 1-baseline
SHARED = PROJECT_ROOT.parent
sys.path.insert(0, str(SCRIPT_DIR.parent))

from reader import collect_cut_err_instances, discover_reps, merge_specs
from tables import run_full_comparison

STAGE = "9-nested"
RES = SHARED / STAGE / "check/results"
PAPER = PROJECT_ROOT / "check/results/paper" / STAGE / "full"

LABEL_MAP = {
    "base_only_div":  r"\texttt{baseline}",
    "bpfi":           r"\texttt{BPFI}",
    "bpfi_gmp_no_ph": r"\texttt{BPFI\_gmp\_no\_ph}",
    "bpfi_no_pd":     r"\texttt{BPFI\_no\_pd}",
    "bpfi_no_sps":    r"\texttt{BPFI\_no\_sps}",
}

_BPFI = discover_reps(RES / "3-nested-bpfi-full")                                    # 3 reps
_BASE = discover_reps(SHARED / "5-partial_decomp/check/results/2-full-eval_num_1_div_full")

SETTINGS = ["base_only_div", "bpfi", "bpfi_gmp_no_ph", "bpfi_no_pd", "bpfi_no_sps"]
_SRC = {s: (_BASE if s == "base_only_div" else _BPFI) for s in SETTINGS}

SOURCES = [(d, [s]) for s in SETTINGS for d in _SRC[s]]


def main() -> None:
    spec = {s: _SRC[s] for s in SETTINGS}
    # Self-scoped exclusions (see module docstring): stage 9 is opted out of the
    # global list, so it computes its own from exactly the runs it reports.
    excl = collect_cut_err_instances(merge_specs(SOURCES), use_obj_deviation=True)
    print(f"[compare_full] 9-nested: re-aggregating {len(spec)} settings "
          f"({len(_BPFI)} bpfi reps, {len(_BASE)} baseline reps), "
          f"{len(excl)} own cut-err instance(s) excluded")
    run_full_comparison(
        spec, baseline="base_only_div", output_dir=PAPER,
        exclude_instances=excl, pd_convergence=False, min_buses=73,
        label_map=LABEL_MAP, include_settings=SETTINGS,
    )


if __name__ == "__main__":
    main()
