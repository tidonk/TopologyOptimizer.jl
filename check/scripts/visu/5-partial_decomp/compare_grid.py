#!/usr/bin/env python3
"""Stage 5 (partial decomposition) — GRID comparisons: num / num_no_base.

Config from ``compare_5-partial_decomp_num_grid.py`` /
``..._num_no_base_grid.py``: baseline ``cf05_sp005`` (No PD),
``target="time_with_warmstart"`` with the ``_add_time_with_warmstart`` transform.

Grid runs are single dirs (not replicated).  They use the **grid** CUT_ERROR
exclusion list — never the full one: the two test sets reuse instance names for
different problems.

Note: Oli's paper versions of these two tables carry a hand-added ``$\\pdmaxn$``
column (= setting num + 1) and drop CutErr; we emit the standard layout, so the
scorer count stays inside the label (``Viol1`` … ``Viol5``).

Run:  python3 compare_grid.py
"""

import sys
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]        # visu/5-partial_decomp -> 1-baseline
SHARED = PROJECT_ROOT.parent
sys.path.insert(0, str(SCRIPT_DIR.parent))

from reader import grid_cut_err_instances
from tables import run_full_comparison

STAGE = "5-partial_decomp"
RES = SHARED / STAGE / "check/results"
PAPER = PROJECT_ROOT / "check/results/paper" / STAGE / "grid"
BASELINE = "cf05_sp005"

_SCORERS = [("eta", "Viol"), ("rnd", "Rnd"), ("cprox", "CProx"),
            ("div_cos", "Div-Cos"), ("div_euc", "Div-Euc"), ("sens", "Sens")]


def _labels(suffix: str) -> dict[str, str]:
    m = {"cf05_sp005": r"\texttt{No PD}", "base_only": r"\texttt{Base only}"}
    for key, disp in _SCORERS:
        for n in range(1, 6):
            m[f"{key}_num{suffix}_{n}"] = rf"\texttt{{{disp}{n}}}"
    return m


def _settings(suffix: str, with_base: bool) -> list[str]:
    out = ["cf05_sp005"] + (["base_only"] if with_base else [])
    for key, _disp in _SCORERS:
        out += [f"{key}_num{suffix}_{n}" for n in range(1, 6)]
    return out


# (leaf, grid dir, settings, label_map)
LEAVES = [
    ("num", RES / "5-partial_decomp_num_grid",
     _settings("", with_base=True), _labels("")),
    ("num_no_base", RES / "5-partial_decomp_num_no_base_grid",
     _settings("_no_base", with_base=False), _labels("_no_base")),
]

SOURCES = [(d, settings) for _n, d, settings, _lm in LEAVES]


def _add_time_with_warmstart(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    t_solve = pd.to_numeric(df.get("time_solve", 0), errors="coerce").fillna(0)
    t_ws = pd.to_numeric(df.get("time_warmstart", 0), errors="coerce").fillna(0)
    # include warmstart time for EVERY run (time_warmstart is 0 for
    # non-warmstarted settings, so this is a no-op there)
    df["time_with_warmstart"] = t_solve + t_ws
    return df


def main() -> None:
    for name, src, settings, label_map in LEAVES:
        if not src.is_dir():
            print(f"warning: {src} missing; skipping {name}")
            continue
        # Per-experiment scoping: drop only THIS grid search's own proven cut
        # errors, so an unrelated grid's defects don't shrink this table.
        excl = grid_cut_err_instances(src, settings)
        print(f"[compare_grid] {name}: {src.name} (baseline {BASELINE}, "
              f"{len(excl)} own cut-err instance(s) excluded)")
        run_full_comparison(
            src, baseline=BASELINE, output_dir=PAPER / name,
            exclude_instances=excl, pd_convergence=False, min_buses=73,
            label_map=label_map, include_settings=settings,
            target="time_with_warmstart", df_transform=_add_time_with_warmstart,
        )


if __name__ == "__main__":
    main()
