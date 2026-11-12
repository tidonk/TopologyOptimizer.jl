#!/usr/bin/env python3
"""Stage 84 (in-out strengthening) — GRID comparison.

Single grid experiment ``in_out``; config verbatim from
``84-in_out_strengthening/.../84-in_out_strengthening/compare_84-in-out.py``.
Zero-touch: runs from 1-baseline, reads the 84 repo, writes staging here.

Metric is **``time_with_in_out`` = ``time_solve + time_in_out``**: ``time_solve``
covers only the Benders iterations, so the strengthening preprocessing would
otherwise be free for the in-out settings.

Labels stay plain (no parameter columns): the rows sweep three *different*
parameters (lambda / delta / stall) plus two toggles, so there is no single
parameter axis to expand into.

Grid runs are single dirs and use **per-experiment** cut-error scoping (this
grid's own proven CUT_ERRORs only) — never the full-set list.

Run:  python3 compare_grid.py
"""

import sys
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]      # visu/84-in_out_strengthening -> 1-baseline
SHARED = PROJECT_ROOT.parent
sys.path.insert(0, str(SCRIPT_DIR.parent))

from reader import grid_cut_err_instances
from tables import run_full_comparison

STAGE = "84-in_out_strengthening"
SRC = SHARED / STAGE / "check/results/in_out"
PAPER = PROJECT_ROOT / "check/results/paper" / STAGE / "grid" / "1-in_out"

LABEL_MAP = {
    "in_out_off":     r"\texttt{in-out off}",
    "in_out_default": r"\texttt{in-out default}",
    "l000":   r"\texttt{lambda=0.0}",
    "l010":   r"\texttt{lambda=0.1}",
    "l025":   r"\texttt{lambda=0.25}",
    "l050":   r"\texttt{lambda=0.5}",
    "l075":   r"\texttt{lambda=0.75}",
    "l100":   r"\texttt{lambda=1.0}",
    "d00000": r"\texttt{delta=0.0}",
    "d00002": r"\texttt{delta=0.0002}",
    "d00010": r"\texttt{delta=0.001}",
    "d00050": r"\texttt{delta=0.005}",
    "s005":   r"\texttt{stall=5}",
    "s010":   r"\texttt{stall=10}",
    "s030":   r"\texttt{stall=30}",
    "s060":   r"\texttt{stall=60}",
    "esls_off": r"\texttt{ESLS off}",
    "esls_on":  r"\texttt{ESLS on}",
}

# Oli's column order.
SETTINGS = ["in_out_off", "in_out_default",
            "l000", "l010", "l025", "l050", "l075", "l100",
            "d00000", "d00002", "d00010", "d00050",
            "s005", "s010", "s030", "s060",
            "esls_off", "esls_on"]

BASELINE = "in_out_off"

SOURCES = [(SRC, SETTINGS)]


def _add_time_with_in_out(df: pd.DataFrame) -> pd.DataFrame:
    """time_with_in_out = time_solve + time_in_out (verbatim from the original)."""
    df = df.copy()
    t_solve = pd.to_numeric(df.get("time_solve", 0), errors="coerce").fillna(0)
    t_in_out = pd.to_numeric(df.get("time_in_out", 0), errors="coerce").fillna(0)
    df["time_with_in_out"] = t_solve + t_in_out
    return df


def main() -> None:
    if not SRC.is_dir():
        print(f"warning: {SRC} missing; skipping")
        return
    excl = grid_cut_err_instances(SRC, SETTINGS)
    print(f"[compare_grid] 1-in_out: {len(SETTINGS)} settings, baseline {BASELINE}, "
          f"{len(excl)} own cut-err instance(s) excluded")
    run_full_comparison(
        SRC, baseline=BASELINE, output_dir=PAPER, target="time_with_in_out",
        exclude_instances=excl, pd_convergence=False, min_buses=73,
        label_map=LABEL_MAP, include_settings=SETTINGS,
        df_transform=_add_time_with_in_out,
    )


if __name__ == "__main__":
    main()
