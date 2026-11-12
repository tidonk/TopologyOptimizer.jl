#!/usr/bin/env python3
"""(Re)build the CUT_ERROR exclusion list for the **full** test set.

Writes ``1-baseline/check/results/CUT_ERR_INSTANCES.txt`` — instances with a cut
error (proven ``CUT_ERROR`` status, or an objective disagreement between settings
that both claim OPTIMAL) in any full run of the pglib-opf set.  All full-run
tables share this one list, so every stage is compared on the same instances.

There is deliberately **no grid list**.  The grid set holds only ~31 instances but
is probed by ~11 grid searches of up to 44 settings, so a global union would
exclude nearly all of it — a defect in the cut-filtering grid would shrink the
unrelated partial-decomposition grid.  Grid exclusions are instead scoped **per
experiment** and computed on the fly by :func:`reader.grid_cut_err_instances`
(proven CUT_ERROR status only).  The grid and full sets reuse instance *names*
for *different* problems, so their exclusions must never be shared anyway.

SCOPE.  Only runs that are finally included in ``across_all/`` are scanned, and
within them only the settings **actually reported** — an instance must not be
discarded because two settings we never publish disagreed on it.  The scan set
is therefore derived from:

* ``across_all.PICKS`` — the canonical per-stage full/grid runs.  All of their
  settings count as reported (``_write_per_stage_final_dumps`` emits a
  ``compiled_final_<setting>.csv`` for every setting it finds there).
* the ``SOURCES`` constant exported by each ``visu/<stage>/compare_full*.py`` /
  ``compare_grid*.py`` — a list of ``(dir, reported_settings)`` pairs.  Those
  modules are import-safe (their work sits behind ``if __name__ == "__main__"``),
  so this stays automatically in sync with whatever the leaves actually report.

``EXCLUDED_STAGES`` opts a stage out of the global union.  ``9-nested`` is opted
out: its BPFI settings have real, surviving cut errors, and including them grew
the list 31 -> 42, shrinking the shared full test set to 111 for *every* other
stage — the same coupling that made a global *grid* union unusable.  Stage 9 is
not part of the paper's across_all, so its defects are scoped to its own leaf
(``9-nested/compare_full.py`` computes its exclusions from its own SOURCES).

Run:  python3 build_cut_err_exclusions.py
"""

import importlib.util
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from reader import merge_specs, update_cut_err_file, CUT_ERR_EXCLUSIONS_FILE  # noqa: E402
import across_all  # noqa: E402


def _load_module(path: Path):
    """Import a compare script by path without executing its main()."""
    spec = importlib.util.spec_from_file_location(f"_cmp_{path.stem}_{path.parent.name}", path)
    if spec is None or spec.loader is None:
        return None
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except Exception as e:  # a stage script may reference missing data
        print(f"  [warn] could not import {path}: {e}")
        return None
    return mod


# Stages whose cut errors must not shrink the shared full test set (see above).
EXCLUDED_STAGES = {"9-nested"}


def _leaf_sources(pattern: str) -> list[tuple[Path, list[str]]]:
    """Collect ``SOURCES`` from every visu/<stage>/<pattern> compare script."""
    out: list[tuple[Path, list[str]]] = []
    for script in sorted(SCRIPT_DIR.glob(f"*/{pattern}")):
        if script.parent.name in EXCLUDED_STAGES:
            print(f"  [scope] skipping {script.parent.name} (opted out of the global list)")
            continue
        mod = _load_module(script)
        src = getattr(mod, "SOURCES", None) if mod else None
        if not src:
            continue
        for d, settings in src:
            out.append((Path(d), list(settings) if settings else None))
    return out


def full_run_specs() -> list[tuple[Path, list[str] | None]]:
    specs: list[tuple[Path, list[str] | None]] = []
    for _d, _sg, _srcg, _sf, src_full in across_all.PICKS:   # canonical full runs
        for d in across_all._full_dirs(src_full):
            specs.append((d, None))                          # all settings reported
    specs += _leaf_sources("compare_full*.py")               # per-leaf reported settings
    return _dedup(specs)


# Single implementation, shared with the leaf scripts that scope their own
# exclusions (see reader.merge_specs).
_dedup = merge_specs


def main() -> None:
    # FULL test set only.  Grid exclusions are scoped per grid experiment and
    # computed on the fly by reader.grid_cut_err_instances(), so there is no
    # grid exclusion file: a global union across the ~11 grid searches would
    # wipe out most of the 31-instance grid set.
    for label, specs, path, use_obj in (
        ("FULL", full_run_specs(), CUT_ERR_EXCLUSIONS_FILE, True),
    ):
        # Clean rebuild: this scans every in-scope source, so the result is
        # complete on its own, and re-scoping actually takes effect.
        result = update_cut_err_file(specs, path, reset=True, use_obj_deviation=use_obj)
        print(f"[{label}] scanned {len(specs)} dir(s) -> {len(result)} excluded "
              f"instance(s) [{'status+objective' if use_obj else 'status only'}] "
              f"-> {path.name}")


if __name__ == "__main__":
    main()
