#!/usr/bin/env python3
"""Stage 83 (branching priority) — FULL-run comparison.

Reproducible replacement for Oli's on-the-fly run.  Two source stacks:

* the six method settings (``{cprox,rnd,sens}_num_{5,10}_div``) come from the
  three repetitions ``2-bprio-full_div{1,2,3}``.  The *non*-div stacks are not
  usable: ``2-bprio-full{1,2,3}`` only ran ``_num_10`` and
  ``2-bprio-full_5{1,2,3}`` is broken (``sens_num_5`` has 5 of 152 instances).
* the ``base_only`` baseline is a **carry-forward** — it was never re-run here.
  The un-suffixed ``2-bprio-full_div`` holds a single pre-computed mean CSV per
  setting; its ``instance_file`` column names the true origin,
  ``5-partial_decomp/settings/5-partial_decomp/2-full-eval_num_1/base_only.toml``,
  whose three repetitions are ``2-full-eval_num_1{1,2,3}``.  We re-aggregate
  from those instead of reading the copied mean.

Zero-touch: runs from 1-baseline, reads 83-branching + 5-partial_decomp,
writes staging under 1-baseline.

Run:  python3 compare_full.py
"""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]          # visu/83-branching -> 1-baseline
SHARED = PROJECT_ROOT.parent
sys.path.insert(0, str(SCRIPT_DIR.parent))

from reader import discover_reps, load_excluded_instances, CUT_ERR_EXCLUSIONS_FILE
from tables import run_full_comparison

STAGE = "83-branching"
RES = SHARED / STAGE / "check/results"
PAPER = PROJECT_ROOT / "check/results/paper" / STAGE / "full"
EXCL = load_excluded_instances(CUT_ERR_EXCLUSIONS_FILE)   # full test set

# Method name kept as its own column, the priority count as $\bpmaxn$.
LABEL_MAP = {
    "base_only":         r"\texttt{baseline}, n=0",
    "cprox_num_5_div":   r"\texttt{CProx}, n=5",
    "cprox_num_10_div":  r"\texttt{CProx}, n=10",
    "rnd_num_5_div":     r"\texttt{Rnd}, n=5",
    "rnd_num_10_div":    r"\texttt{Rnd}, n=10",
    "sens_num_5_div":    r"\texttt{Sens}, n=5",
    "sens_num_10_div":   r"\texttt{Sens}, n=10",
}
PARAM_COLS = [(r"$\bpmaxn$", "n")]

_DIV  = discover_reps(RES / "2-bprio-full_div")                                # 3 reps
_BASE = discover_reps(SHARED / "5-partial_decomp/check/results/2-full-eval_num_1")  # 3 reps

SETTINGS = ["base_only",
            "cprox_num_5_div", "cprox_num_10_div",
            "rnd_num_5_div", "rnd_num_10_div",
            "sens_num_5_div", "sens_num_10_div"]

_SRC = {s: (_BASE if s == "base_only" else _DIV) for s in SETTINGS}

# (dir, settings reported from it) — consumed by build_cut_err_exclusions.py
SOURCES = [(d, [s]) for s in SETTINGS for d in _SRC[s]]


def main() -> None:
    spec = {s: _SRC[s] for s in SETTINGS}
    print(f"[compare_full] 83-branching: re-aggregating {len(spec)} settings "
          f"({len(_DIV)} div reps, {len(_BASE)} baseline reps)")
    run_full_comparison(
        spec, baseline="base_only", output_dir=PAPER,
        exclude_instances=EXCL, pd_convergence=False, min_buses=73,
        label_map=LABEL_MAP, include_settings=SETTINGS,
        param_columns=PARAM_COLS, param_keep_label=True,
    )


if __name__ == "__main__":
    main()
