#!/usr/bin/env python3
"""Stage 9 (nested Benders) — GRID comparisons.

Two leaves:

* ``1-nested_grid`` — the nested-parameter sweep; config verbatim from
  ``9-nested/.../9-nested/compare_1-nested_grid.py`` (baseline
  ``benders_LEQRF``, i.e. "no nested at all").
* ``2-nested-bpfi-grid`` — the BPFI ablation sweep.  It has **no original
  compare script**; the baseline here is ``bpfi`` (the full variant), against
  which the ``bpfi_*`` ablations are measured.  ``bpfi_spc`` is dropped: the
  setting directory exists but contains **0** result CSVs.

Grid runs are single dirs and use **per-experiment** cut-error scoping (this
grid's own proven CUT_ERRORs only) — never the full-set list.

Run:  python3 compare_grid.py
"""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]      # visu/9-nested -> 1-baseline
SHARED = PROJECT_ROOT.parent
sys.path.insert(0, str(SCRIPT_DIR.parent))

from reader import grid_cut_err_instances
from tables import run_full_comparison

STAGE = "9-nested"
RES = SHARED / STAGE / "check/results"
PAPER = PROJECT_ROOT / "check/results/paper" / STAGE / "grid"

NESTED_LABELS = {
    "benders_LEQRF":          r"\texttt{Benders (no nested)}",
    "base":                   r"\texttt{base (defaults)}",
    "mp_heur_off":            r"\texttt{MP heur off}",
    "base_pf_intermediate":   r"\texttt{base PF in interm.}",
    "no_change_stop":         r"\texttt{no-change stop}",
    "pb_stop_off":            r"\texttt{PB stop off}",
    "cutpool_500":            r"\texttt{cut pool 500}",
    "cutpool_5000":           r"\texttt{cut pool 5000}",
    "same_sol_thresh_1":      r"\texttt{same-sol thresh 1}",
    "same_sol_thresh_10":     r"\texttt{same-sol thresh 10}",
    "gen_in_mp":              r"\texttt{gen in MP}",
    "gen_in_mp_with_sp_cuts": r"\texttt{gen in MP + SP cuts}",
    "sp_cuts_only":           r"\texttt{SP cuts only}",
}

BPFI_LABELS = {
    "bpfi":           r"\texttt{BPFI}",
    "bpfi_cp0":       r"\texttt{BPFI\_cp0}",
    "bpfi_gmp":       r"\texttt{BPFI\_gmp}",
    "bpfi_gmp_no_ph": r"\texttt{BPFI\_gmp\_no\_ph}",
    "bpfi_ncs":       r"\texttt{BPFI\_ncs}",
    "bpfi_no_pbs":    r"\texttt{BPFI\_no\_pbs}",
    "bpfi_no_pd":     r"\texttt{BPFI\_no\_pd}",
    "bpfi_no_ph":     r"\texttt{BPFI\_no\_ph}",
    "bpfi_no_sps":    r"\texttt{BPFI\_no\_sps}",
}

# (leaf, experiment dir, baseline, settings, label map, min_buses)
LEAVES = [
    ("1-nested_grid", "1-nested_grid", "benders_LEQRF",
     list(NESTED_LABELS), NESTED_LABELS, 73),
    ("2-nested-bpfi-grid", "2-nested-bpfi-grid", "bpfi",
     list(BPFI_LABELS), BPFI_LABELS, 73),
]

SOURCES = [(RES / exp, settings) for _n, exp, _b, settings, _l, _m in LEAVES]


def main() -> None:
    for name, exp, baseline, settings, label_map, min_buses in LEAVES:
        src = RES / exp
        if not src.is_dir():
            print(f"warning: {src} missing; skipping {name}")
            continue
        excl = grid_cut_err_instances(src, settings)
        print(f"[compare_grid] {name}: {len(settings)} settings, baseline {baseline}, "
              f"{len(excl)} own cut-err instance(s) excluded")
        run_full_comparison(
            src, baseline=baseline, output_dir=PAPER / name,
            exclude_instances=excl, pd_convergence=False, min_buses=min_buses,
            label_map=label_map, include_settings=settings,
        )


if __name__ == "__main__":
    main()
