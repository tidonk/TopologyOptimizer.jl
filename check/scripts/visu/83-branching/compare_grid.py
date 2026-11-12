#!/usr/bin/env python3
"""Stage 83 (branching priority) — GRID comparisons.

Two leaves, mirroring what Oli kept in ``across_all/83-branching``:

* ``grid`` — the **current** grid, ``1-bprio_grid_div`` (config verbatim from
  ``83-branching/.../83-branching/compare.py``): 21 settings, 4 score types
  (``cprox``/``rnd``/``sens``/``eta``, the latter labelled ``Viol``) x
  ``MAX_NUM`` in {1,2,5,10,20}, plus ``base_only_div``.  Its declared baseline
  (``base_only``) is absent from the dir, so no Wilcoxon columns are produced —
  we state that explicitly here with ``baseline=""`` instead of relying on a
  silent lookup miss.
* ``1-bprio_grid`` — the **superseded** non-div grid
  (``compare_83-bprio_grid.py``), on ``time_with_warmstart``.  Oli dropped the
  ``scopf`` ("Warmstart SCOPF") row and made ``base_only`` the reference; that
  is reproduced here.  ``eta`` and ``*_num_30`` were never part of it, and
  ``rnd`` only ran ``num_{5,10}``.

Grid runs are single dirs and use **per-experiment** cut-error scoping (this
grid's own proven CUT_ERRORs only) — never the full-set list.

Run:  python3 compare_grid.py
"""

import sys
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]        # visu/83-branching -> 1-baseline
SHARED = PROJECT_ROOT.parent
sys.path.insert(0, str(SCRIPT_DIR.parent))

from reader import grid_cut_err_instances
from tables import run_full_comparison

STAGE = "83-branching"
RES = SHARED / STAGE / "check/results"
PAPER = PROJECT_ROOT / "check/results/paper" / STAGE / "grid"

# Method name stays in its own column, MAX_NUM becomes $\bpmaxn$.
METHOD_LABEL = {"cprox": "CProx", "rnd": "Rnd", "sens": "Sens", "eta": "Viol"}
NUMS = [1, 2, 5, 10, 20]


def _label_map(suffix: str, methods, nums=NUMS) -> dict[str, str]:
    lm = {f"base_only{suffix}": r"\texttt{baseline}, n=0"}
    for m in methods:
        for n in nums:
            lm[f"{m}_num_{n}{suffix}"] = rf"\texttt{{{METHOD_LABEL[m]}}}, n={n}"
    return lm


PARAM_COLS = [(r"$\bpmaxn$", "n")]


def _add_time_with_warmstart(df: pd.DataFrame) -> pd.DataFrame:
    """time_with_warmstart = time_solve + time_warmstart (verbatim from the
    original ``compare_83-bprio_grid.py``); only the priority settings warmstart."""
    df = df.copy()
    t_solve = pd.to_numeric(df.get("time_solve", 0), errors="coerce").fillna(0)
    t_ws = pd.to_numeric(df.get("time_warmstart", 0), errors="coerce").fillna(0)
    # include warmstart time for EVERY run (time_warmstart is 0 for
    # non-warmstarted settings, so this is a no-op there)
    df["time_with_warmstart"] = t_solve + t_ws
    return df


# --- leaf definitions -------------------------------------------------------
_DIV_METHODS = ["cprox", "rnd", "sens", "eta"]
_DIV_SETTINGS = ["base_only_div"] + [f"{m}_num_{n}_div"
                                     for m in _DIV_METHODS for n in NUMS]

# The superseded grid: no eta, rnd only ran num_{5,10}, num_30 excluded.
_OLD_SETTINGS = (["base_only"]
                 + [f"sens_num_{n}" for n in NUMS]
                 + [f"cprox_num_{n}" for n in NUMS]
                 + [f"rnd_num_{n}" for n in (5, 10)])
_OLD_LABELS = _label_map("", ["cprox", "sens"])
_OLD_LABELS.update({f"rnd_num_{n}": rf"\texttt{{Rnd}}, n={n}" for n in (5, 10)})

# (leaf, experiment dir, settings, label map, target, df_transform)
LEAVES = [
    ("1-bprio_grid_div", "1-bprio_grid_div", _DIV_SETTINGS,
     _label_map("_div", _DIV_METHODS), "time_solve", None),
    ("1-bprio_grid", "1-bprio_grid", _OLD_SETTINGS,
     _OLD_LABELS, "time_with_warmstart", _add_time_with_warmstart),
]

SOURCES = [(RES / exp, settings) for _n, exp, settings, _l, _t, _f in LEAVES]


def main() -> None:
    for name, exp, settings, label_map, target, transform in LEAVES:
        src = RES / exp
        if not src.is_dir():
            print(f"warning: {src} missing; skipping {name}")
            continue
        excl = grid_cut_err_instances(src, settings)
        print(f"[compare_grid] {name}: {len(settings)} settings, target {target}, "
              f"{len(excl)} own cut-err instance(s) excluded")
        run_full_comparison(
            src, baseline="", output_dir=PAPER / name, target=target,
            exclude_instances=excl, pd_convergence=False, min_buses=0,
            label_map=label_map, include_settings=settings,
            param_columns=PARAM_COLS, param_keep_label=True,
            df_transform=transform,
        )


if __name__ == "__main__":
    main()
