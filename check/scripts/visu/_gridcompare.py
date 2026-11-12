"""Shared helpers for the per-experiment grid/full compare scripts.

Keeps each ``visu/<stage>/{grid,full}/compare_*.py`` thin: the cut-error
penalty transform, the "best cell" ranking (with a cut-error gate) and small
markdown helpers live here so no logic is copy-pasted across experiments.

Cut-error policy (see EXPERIMENTS.md §12.3d): cut errors are **penalized to
TIME_LIMIT, not dropped** — every instance stays in the comparison — and routed
to the ``CutErr`` column so timeouts (TL) and correctness failures stay
distinguishable.  ``best_from_short`` then ranks only among runs with at most
``max_cut_err`` cut errors.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from reader import (CUT_ERR_OBJ_TOLERANCE, CUT_ERR_EXCLUSIONS_FILE,
                    discover_reps, grid_cut_err_instances,
                    load_excluded_instances)
from selection import TIME_LIMIT  # noqa: F401  (re-exported for callers)
from tables import run_full_comparison

DEFAULT_MAX_CUT_ERR = 2

# Common paths so compare scripts never compute their own.
_VISU = Path(__file__).resolve().parent               # …/1-baseline/check/scripts/visu
PROJECT_ROOT = _VISU.parents[2]                        # …/1-baseline
SHARED = PROJECT_ROOT.parent                           # experiment root
PAPER = PROJECT_ROOT / "check/results/paper"           # zero-touch output root


def penalize_cut_errors(target: str = "time_solve"):
    """Return a ``df_transform`` that keeps every instance but charges each
    cut-erroring *cell* a TIME_LIMIT penalty, routed to the CutErr column.

    * explicit ``status == CUT_ERROR`` needs no touch — ``CUT_ERROR`` is not a
      BROKEN status, so it stays in ``ran`` and the inclusive ``_all`` view
      already clips it to TIME_LIMIT while ``CutErr`` still counts it;
    * objective-deviation cut errors (``OPTIMAL`` but off the per-instance
      consensus — the median of the OPTIMAL runs — by more than
      :data:`reader.CUT_ERR_OBJ_TOLERANCE`) are relabelled ``CUT_ERROR`` so they
      leave the solved set, show in CutErr and get the same value penalty.  A
      setting that cuts off the true optimum finishes *fast*, so it must be
      penalized, not rewarded.
    """
    def _tf(df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()
        n_expl = n_dev = 0
        if "status" in df.columns:
            n_expl = int((df["status"].astype(str) == "CUT_ERROR").sum())
        if {"status", "instance", "objective", "settings"}.issubset(df.columns):
            opt = df["status"].astype(str) == "OPTIMAL"
            obj = pd.to_numeric(df["objective"], errors="coerce")
            cons = obj.where(opt).groupby(df["instance"]).transform("median")
            rel = (obj - cons).abs() / cons.abs().clip(lower=1e-10)
            bad = opt & rel.notna() & (rel > CUT_ERR_OBJ_TOLERANCE)
            n_dev = int(bad.sum())
            df.loc[bad, "status"] = "CUT_ERROR"
        print(f"    cut errors -> CutErr (time-penalized in _all): "
              f"{n_expl} explicit + {n_dev} obj-deviation cell(s)")
        return df
    return _tf


def _num(cell):
    m = re.search(r"[-+]?\d+\.?\d*", str(cell).replace("**", ""))
    return float(m.group()) if m else None


def best_from_short(md_path: str | Path, max_cut_err: int = DEFAULT_MAX_CUT_ERR) -> None:
    """Print the fastest setting (by Sum and Geo. Mean) from a short ``_all``
    table, considering only rows with at most ``max_cut_err`` cut errors.

    Works with or without ``param_columns``: the leading cells (everything
    before the ``Solved`` column) identify the row and are joined for display.
    """
    md_path = Path(md_path)
    if not md_path.is_file():
        print(f"    (no {md_path.name} to rank)")
        return
    lines = [l for l in md_path.read_text().splitlines() if l.lstrip().startswith("|")]
    if len(lines) < 3:
        return
    cols = [c.strip() for c in lines[0].strip("|").split("|")]
    idx = {c: i for i, c in enumerate(cols)}
    ce_i = idx.get("CutErr")
    id_end = idx.get("Solved", 1)          # leading cells before the metric block
    skip = ("solved", "sum", "mean", "geo", "wilcoxon", "significance", "tl", "cuterr")
    rows = []
    for l in lines[2:]:
        c = [x.strip() for x in l.strip("|").split("|")]
        if len(c) < len(cols):
            continue
        ident = " ".join(x.replace("**", "").strip() for x in c[:id_end]).strip()
        if not ident or ident.lower().startswith(skip):
            continue
        ce = int(_num(c[ce_i]) or 0) if ce_i is not None else 0
        rows.append((ident, ce, _num(c[idx.get("Sum", -1)]), _num(c[idx.get("Geo. Mean", -1)])))
    elig = [r for r in rows if r[1] <= max_cut_err]
    dropped = len(rows) - len(elig)
    print(f"    eligible (≤{max_cut_err} cut-err): {len(elig)}/{len(rows)}"
          + (f"  [{dropped} excluded]" if dropped else ""))
    for key, name in ((2, "Sum"), (3, "Geo. Mean")):
        cand = [r for r in elig if r[key] is not None]
        if cand:
            b = min(cand, key=lambda r: r[key])
            print(f"    best by {name:<9}: {b[0]}  ({b[key]:g}, {b[1]} cut-err)")


def _compose(*fns):
    """Left-to-right function composition, skipping ``None`` (``f(g(x))`` order:
    the first non-None runs first).  Returns ``None`` if nothing to do."""
    fns = [f for f in fns if f is not None]
    if not fns:
        return None
    def _run(df):
        for f in fns:
            df = f(df)
        return df
    return _run


def _present_settings(source) -> list[str]:
    """Setting names that actually carry result CSVs.  Works for a single dir,
    a list of repetition dirs (union), or a dict source (its keys)."""
    if isinstance(source, dict):
        return list(source)
    dirs = [Path(d) for d in source] if isinstance(source, (list, tuple)) else [Path(source)]
    present: set[str] = set()
    for d in dirs:
        if d.is_dir():
            present |= {s.name for s in d.iterdir() if s.is_dir() and any(s.glob("*.csv"))}
    return sorted(present)


def run_compare(*, source, output, baseline, labels=None, param_columns=None,
                target: str = "time_solve", max_cut_err: int = DEFAULT_MAX_CUT_ERR,
                min_buses: int = 0, name: str = "", cut: str = "penalize",
                pre_transform=None, reps: bool = False, min_time: float = 0.0,
                style_map=None, param_keep_label=None, rank: bool = True):
    """One-call grid/full comparison for a thin compare script.

    ``source`` is an experiment dir (often an ``*-across/`` folder of symlinks),
    a list of repetition dirs, or a ``{setting: [dirs]}`` map.  With
    ``reps=True`` a single base dir is expanded to its numbered repetitions via
    :func:`reader.discover_reps` and averaged.

    ``labels`` maps setting name -> LaTeX label (``None`` = use folder names).
    Only settings that both appear in ``labels`` and have result CSVs are used.

    ``cut`` selects the cut-error policy (see EXPERIMENTS.md §12.3d):

    * ``"penalize"`` (default, for *new* grids) — cut errors are kept but
      time-penalized to TIME_LIMIT and routed to CutErr; settings with more than
      ``max_cut_err`` of them are dropped from the short table and the best
      eligible cell is printed.
    * ``"exclude_grid"`` — drop this grid's own proven CUT_ERROR instances
      (per-experiment scoping); reproduces the established paper-grid numbers.
    * ``"exclude_full"`` — drop the shared full-set CUT_ERROR list
      (``CUT_ERR_INSTANCES.txt``); for full (replicated) runs.
    * ``"none"`` — no cut-error handling.

    ``pre_transform`` (e.g. adding a ``time_with_warmstart`` column) runs before
    the cut policy.  Writes tables/plots under ``output`` only.
    """
    output = Path(output)
    if reps:
        source = discover_reps(source)
        if not source:
            print(f"warning: no repetitions for {name or output.name}; skipping"); return
        print(f"[{name or output.name}] averaging {len(source)} reps: {[p.name for p in source]}")
    elif not isinstance(source, (list, tuple, dict)):
        source = Path(source)
        if not source.is_dir():
            print(f"warning: {source} missing; skipping"); return

    present = _present_settings(source)
    if labels is None:
        labels = {s: s for s in present}
    settings = [s for s in labels if s in present]
    if not settings:
        print(f"warning: no labelled settings present for {name or output.name}; skipping"); return
    if baseline and baseline not in settings:
        print(f"warning: baseline {baseline!r} has no data yet; skipping"); return

    # --- cut-error policy -------------------------------------------------
    exclude = None
    penalty = None
    drop = None
    if cut == "penalize":
        penalty, drop = penalize_cut_errors(target), max_cut_err
    elif cut == "exclude_grid":
        # Grids are single dirs; per-experiment scoping (this grid's own defects).
        exclude = grid_cut_err_instances(source, settings)
    elif cut == "exclude_full":
        exclude = load_excluded_instances(CUT_ERR_EXCLUSIONS_FILE)
    elif cut != "none":
        raise ValueError(f"unknown cut mode {cut!r}")

    print(f"[{name or output.name}] {len(settings)} settings, baseline "
          f"{baseline or '(none)'}, cut={cut}"
          + (f", {len(exclude)} excluded" if exclude else ""))
    run_full_comparison(
        source, baseline=baseline, output_dir=output, target=target,
        df_transform=_compose(pre_transform, penalty), pd_convergence=False,
        min_buses=min_buses, min_time=min_time, style_map=style_map,
        label_map={s: labels[s] for s in settings}, exclude_instances=exclude,
        include_settings=settings, param_columns=param_columns,
        param_keep_label=(bool(param_columns) if param_keep_label is None
                          else param_keep_label),
        cut_err_drop=drop,
    )
    if rank:
        best_from_short(output / "comparison_short_all.md", max_cut_err)
