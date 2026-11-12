#!/usr/bin/env python3
"""Cross-stage comparison runner: collects the best setting from each algorithmic
stage and compares them side-by-side on a common instance set.

HOW IT WORKS
============
1. ``PICKS`` defines one representative setting per stage (display name, setting
   key, result directories for the 31-instance grid run and the full 174-instance
   run).
2. ``_collect(grid=True/False)`` loads each pick's data.  For the full run,
   multiple result directories are merged and the arithmetic mean of ``time_solve``
   is computed per instance across all available runs.
3. ``_write(tag, df)`` writes every table and plot to
   ``1-baseline/check/results/across_all/<tag>/``:
     - ``comparison{,_all}.{md,tex}``          long per-instance tables
     - ``comparison_short{,_all}.{md,tex}``    aggregated summary tables
     - ``performance_profile_time_solve.{pdf,png}``

OUTPUT FOLDERS
==============
* ``FINAL_GRID`` — 31-instance ``instances_grid/`` set (paper headline table).
  Stages without a dedicated grid run are skipped.
* ``FINAL``      — full 174-instance pglib-opf set; multiple partial runs are
  averaged per instance (arithmetic mean, non-OPTIMAL penalised to TIME_LIMIT).

ADDING / CHANGING PICKS
=======================
Each entry in ``PICKS`` is a 5-tuple:

    (display_name, setting_grid, src_grid, setting_full, src_full)

* ``display_name``  – label used in all output tables and plots.
* ``setting_grid``  – folder/setting name inside the grid result directory.
* ``src_grid``      – Path to grid result directory, or ``None`` to skip FINAL_GRID.
* ``setting_full``  – setting name in the full-run directories (``None`` = same as
                      ``setting_grid``).
* ``src_full``      – base-name Path, list of Paths, or ``None`` to skip FINAL.
                      A **base-name Path** (e.g. ``.../results/2-SP_solving``)
                      auto-discovers its numbered repetition directories
                      (``2-SP_solving1``/``2``/``3``, any count) via
                      :func:`reader.discover_reps`; this is the preferred form
                      and removes the old hard-coded 1/2/3 assumption.  A list
                      of explicit Paths is still accepted (used verbatim).  In
                      both cases the directories are merged and averaged.
"""

import sys
from pathlib import Path
import os

import pandas as pd

# ---------------------------------------------------------------------------
# Bootstrap: make sure the visu helpers in this folder are importable.
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from reader import (  # noqa: E402
    read_experiment, discover_reps, mean_per_instance_full, load_excluded_instances,
    drop_excluded, grid_cut_err_instances, CUT_ERR_EXCLUSIONS_FILE,
)
from tables import (  # noqa: E402
    TIME_LIMIT,
    comparison_short_all_table,
    comparison_short_table,
    comparison_short_to_latex,
    comparison_short_to_markdown,
    comparison_table,
    to_latex,
    to_markdown,
)
from perf_profile import performance_profile  # noqa: E402

# ---------------------------------------------------------------------------
# Global configuration
# ---------------------------------------------------------------------------
SHARED = Path(os.environ.get("TOPOLOGYOPTIMIZER_EXPERIMENT_ROOT", Path(__file__).resolve().parents[4]))       # root of the experiment directories
R = SHARED                              # shorthand used inside PICKS below

# All outputs are written under this directory.
OUTPUT_ROOT = SCRIPT_DIR.parents[1] / "results" / "across_all"

# Column in the experiment DataFrames to use as the performance metric.
TARGET = "time_solve"

# The setting that serves as Wilcoxon-test baseline (must match a display_name).
BASELINE = "Benders LEQRF"

# CUT_ERROR exclusions.  The FULL runs share one global list.  GRID exclusions
# are scoped per grid experiment and computed on the fly (a global union across
# the ~11 grid searches would wipe out most of the 31-instance grid set), and are
# never mixed with the full list — same instance names, different problems.
_EXCLUDED = load_excluded_instances(CUT_ERR_EXCLUSIONS_FILE)            # full runs

# ---------------------------------------------------------------------------
# PICKS: one representative entry per algorithmic stage.
# Format: (display_name, setting_grid, src_grid, setting_full, src_full)
# See module docstring for details.
# ---------------------------------------------------------------------------
PICKS = [
    # Stage 1 — compact Benders (no SP selection or cut filtering)
    ("Compact",
     "compact",
     R / "1-baseline/check/results/1-baseline-compact-grid",
     None,                                 # setting_full = same as setting_grid
     R / "1-baseline/check/results/1-baseline-compact"),   # base name → -compact1/2/3

    # Stage 2 — full-LP Benders (warm-start with full LP relaxation)
    ("Benders Full LP",
     "benders_FLLP",
     None,                                 # no dedicated grid run → skipped in FINAL_GRID
     None,
     R / "2-SP_solving/check/results/2-SP_solving"),       # base name → 2-SP_solving1/2/3

    # Stage 2 — LEQ+RF Benders (our main stage-2 baseline)
    ("Benders LEQRF",
     "benders_LEQRF",
     R / "2-SP_solving/check/results/2-SP_solving_grid",
     None,
     R / "2-SP_solving/check/results/2-SP_solving"),

    # Stage 4 — cut filtering only (best cut-filtering setting)
    ("Cut Filtering",
     "cf10_sp100",
     R / "4-cut_filtering/check/results/6-cutfiltering-abs-spsel-rel",
     "cutfiltering",                       # setting key differs in full-run dirs
     R / "4-cut_filtering/check/results/4-cut_filtering"),  # base → 4-cut_filtering1/2/3

    # Stage 3 — SP selection only (best scenario-selection setting)
    ("SP Selection",
     "rnd_abs_2",
     R / "3-SP_selection/check/results/1-n_hits_nostab",
     "spsel",
     R / "3-SP_selection/check/results/3-SP_selection"),    # base → 3-SP_selection1/2/3

    # Stage 4 — combined SP selection + cut filtering (best combined setting)
    ("SP Sel. + Cut Filt.",
     "cf05_sp005",
     R / "4-cut_filtering/check/results/6-cutfiltering-abs-spsel-rel",
     "cutfiltering_spsel",
     R / "4-cut_filtering/check/results/4-cut_filtering"),

    # Stage 5 — partial decomposition (FINAL_GRID only — full-run pending)
    ("Partial Decomp",
     "base_only",
     R / "5-partial_decomp/check/results/5-partial_decomp_num_grid",
     "partial_decomp",
     R / "5-partial_decomp/check/results/5-partial_decomp"),  # base → 5-partial_decomp1/2/3

    # Stage 82 — fractional separation (FINAL_GRID only)
    ("Frac. Sepa.",
     "all_nodes",
     R / "82-fractional_sepa/check/results/frac_sepa",
     None,
     None),

    # Stage 83 — branching priority (FINAL_GRID only)
    ("Branching",
     "cprox_num_20",
     R / "83-branching/check/results/1-bprio_grid",
     None,
     None),

    # Stage 84 — in/out strengthening (FINAL_GRID only)
    ("In/Out Streng.",
     "s005",
     R / "84-in_out_strengthening/check/results/in_out",
     None,
     None),

    # Stage 85 — Benders primal heuristic (FINAL_GRID only)
    ("Benders Heur.",
     "benders-primal-heur",
     R / "85-benders-primal-heur/check/results/benders-primal-heur",
     None,
     None),

     # Stage 9 — nested Benders (FINAL and FINAL_GRID)
    ("Nested Benders",
     "base_pf_intermediate",
     R / "9-nested/check/results/1-nested_grid",
     None,
     R / "9-nested/check/results/9-nested-bpfi"),           # base → 9-nested-bpfi1
]


# ---------------------------------------------------------------------------
# Helper: average multiple runs per instance
# ---------------------------------------------------------------------------

def _full_dirs(src_full) -> list[Path]:
    """Resolve a pick's ``src_full`` field to a concrete list of run dirs.

    * ``None``               → ``[]`` (stage skipped in FINAL).
    * ``list`` of Paths      → used verbatim (explicit repetition dirs).
    * a single base-name Path → expanded via :func:`reader.discover_reps` into
      its numbered repetition dirs (``<base>1``/``2``/``3``, any count).
    """
    if src_full is None:
        return []
    if isinstance(src_full, (list, tuple)):
        return list(src_full)
    return discover_reps(src_full)


# ---------------------------------------------------------------------------
# Helper: build time-limit–penalised copy of a DataFrame
# ---------------------------------------------------------------------------

def _penalize_tl(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy where non-OPTIMAL rows have ``TARGET`` set to TIME_LIMIT.

    Used for the ``comparison_all`` variant where unsolved instances are shown
    with their penalty value instead of being marked as errors.
    """
    out = df.copy()
    out.loc[out["status"] != "OPTIMAL", TARGET] = TIME_LIMIT
    return out


# ---------------------------------------------------------------------------
# Data collection
# ---------------------------------------------------------------------------

def _collect(grid: bool) -> pd.DataFrame:
    """Load and combine rows for all picks.

    Parameters
    ----------
    grid : if True, load the 31-instance grid run (``src_grid``);
           if False, load the full 174-instance run (``src_full``) with
           arithmetic mean aggregation over all provided directories.

    Returns
    -------
    DataFrame with a ``settings`` column containing the display name of each
    pick.  Empty picks are skipped with a printed warning.
    """
    rows: list[pd.DataFrame] = []

    for display, setting_grid, src_grid, setting_full, src_full in PICKS:

        if grid:
            # ---- Grid mode: single directory, filter by setting name ----
            if src_grid is None or not src_grid.exists():
                print(f"  [skip] {display}: no source")
                continue
            try:
                sub = read_experiment(src_grid)
            except Exception as e:
                print(f"  [skip] {display}: {e}")
                continue
            sub = sub[sub["settings"] == setting_grid]
            if sub.empty:
                print(f"  [skip] {display}: setting '{setting_grid}' missing in {src_grid.name}")
                continue
            # Per-experiment grid scoping: this grid dir's own cut errors only.
            sub = drop_excluded(sub, grid_cut_err_instances(src_grid, [setting_grid]))
            src_label = src_grid.name

        else:
            # ---- Full-run mode: merge multiple dirs, average per instance ----
            setting = setting_full if setting_full is not None else setting_grid
            srcs: list[Path] = _full_dirs(src_full)
            if not srcs:
                print(f"  [skip] {display}: no source")
                continue

            parts: list[pd.DataFrame] = []
            for s in srcs:
                if not s.exists():
                    continue
                try:
                    part = read_experiment(s)
                    filtered = part[part["settings"] == setting]
                    if not filtered.empty:
                        parts.append(filtered)
                except Exception as e:
                    print(f"  [warn] {display} {s.name}: {e}")

            if not parts:
                print(f"  [skip] {display}: setting '{setting}' not found in sources")
                continue

            # Merge all partial runs and compute per-instance arithmetic mean.
            # Each part is filtered to a single setting, so the full-schema mean
            # yields one row per instance (extra averaged columns are harmless;
            # only instance/tal/settings/time_solve/status are used downstream).
            merged = pd.concat(parts, ignore_index=True)
            if _EXCLUDED and "instance" in merged.columns:
                merged = merged[~merged["instance"].astype(str).isin(_EXCLUDED)]
            sub = mean_per_instance_full(merged)
            src_label = f"{len(parts)} dirs"

        # Rename the setting column to the human-readable display name
        sub = sub.copy()
        sub["settings"] = display
        rows.append(sub)
        print(f"  [pick] {display}: {len(sub)} rows from {src_label}")

    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


# ---------------------------------------------------------------------------
# Output writer
# ---------------------------------------------------------------------------

def _write(tag: str, df: pd.DataFrame) -> None:
    """Write all tables and the performance profile for one experiment set.

    Outputs written to ``OUTPUT_ROOT / tag /``:
    * ``comparison.md / .tex``          per-instance table; † = solved by all settings
    * ``comparison_all.md / .tex``      same but with TL-penalised unsolved instances
    * ``comparison_short.md / .tex``    summary (solved count, mean, geo. mean, Wilcoxon)
    * ``comparison_short_all.md / .tex``  as above but over all instances
    * ``performance_profile_time_solve.pdf / .png``

    Parameters
    ----------
    tag : folder name under ``OUTPUT_ROOT`` (e.g. "FINAL" or "FINAL_GRID").
    df  : combined DataFrame from ``_collect``.
    """
    if df.empty:
        print(f"  [skip] {tag}: empty")
        return

    out_dir = OUTPUT_ROOT / tag
    out_dir.mkdir(parents=True, exist_ok=True)

    # Preserve the PICKS order for columns / legend entries
    present = set(df["settings"].unique())
    order = [d for d, *_ in PICKS if d in present]

    # Choose baseline; fall back to first pick if BASELINE is not present
    baseline = BASELINE if BASELINE in present else order[0]
    if baseline != BASELINE:
        print(f"  [warn] baseline '{BASELINE}' missing, using '{baseline}'")

    # ---- Long per-instance comparison tables ----
    # comparison.md    : dagger (†) marks instances solved by *all* settings (strict)
    # comparison_all.md: dagger marks instances run by all; unsolved → TIME_LIMIT
    for suffix, df_in, mode in (
        ("",     df,               "optimal"),
        ("_all", _penalize_tl(df), "ran"),
    ):
        tbl = comparison_table(
            df_in, target=TARGET, baseline=baseline,
            include_settings=order, wilcoxon_test=True,
            mode=mode,
        )
        (out_dir / f"comparison{suffix}.md").write_text(to_markdown(tbl))
        (out_dir / f"comparison{suffix}.tex").write_text(to_latex(tbl))

    # ---- Short summary tables ----
    short = comparison_short_table(
        df, target=TARGET, baseline=baseline,
        include_settings=order, df_with_errors=df, wilcoxon_test=True,
    )
    (out_dir / "comparison_short.md").write_text(comparison_short_to_markdown(short))
    (out_dir / "comparison_short.tex").write_text(comparison_short_to_latex(short))

    short_all = comparison_short_all_table(
        df, target=TARGET, baseline=baseline,
        include_settings=order, wilcoxon_test=True,
    )
    (out_dir / "comparison_short_all.md").write_text(comparison_short_to_markdown(short_all))
    (out_dir / "comparison_short_all.tex").write_text(comparison_short_to_latex(short_all))

    # ---- Dolan-Moré performance profile ----
    # Drop instances that hit TIME_LIMIT for *every* setting before plotting,
    # so the profile reflects only instances where at least one configuration
    # produced useful data. Report the number of dropped instances.
    group_cols = ["instance", "tal"] if "tal" in df.columns else ["instance"]
    is_tl = df["status"] != "OPTIMAL"
    all_tl = df.assign(_tl=is_tl).groupby(group_cols)["_tl"].all()
    drop_keys = all_tl[all_tl].reset_index()[group_cols]
    n_dropped = len(drop_keys)
    if n_dropped:
        df_pp = df.merge(drop_keys, on=group_cols, how="left", indicator=True)
        df_pp = df_pp[df_pp["_merge"] == "left_only"].drop(columns="_merge")
        print(f"  [info] {tag}: dropped {n_dropped} all-TIME_LIMIT instance(s) "
              f"from performance profile")
    else:
        df_pp = df
    performance_profile(
        df_pp, target=TARGET,
        output=out_dir / f"performance_profile_{TARGET}.pdf",
        baseline=baseline, include_settings=order,
    )

    print(f"  [ok]   {tag}: {len(df)} rows, {df['settings'].nunique()} settings -> {out_dir}")


# ---------------------------------------------------------------------------
# Per-stage final-run dumps
# ---------------------------------------------------------------------------

def _stage_name_from_path(path: Path) -> str | None:
    """Extract the algorithmic-stage name from a result-directory path.

    Expects paths under ``<experiment-root>/<stage>/check/results/<exp>``.
    Returns ``<stage>`` (e.g. ``"2-SP_solving"``) or ``None`` if the layout
    does not match.
    """
    try:
        rel = path.resolve().relative_to(SHARED)
    except ValueError:
        return None
    parts = rel.parts
    return parts[0] if parts else None


def _write_per_stage_final_dumps() -> None:
    """For every stage that has 3 FINAL-repetition directories, write one
    pair of intermediate CSVs **per setting present in the source dirs** into
    ``OUTPUT_ROOT / <stage>/``:

    * ``compiled_final_<setting>.csv``       — full-schema concatenation of
      all repetitions for that setting; an extra ``repetition`` column tags
      the source directory (trailing digit, e.g. ``1``/``2``/``3``).
    * ``compiled_final_<setting>_mean.csv``  — per-instance mean across
      repetitions (TL-penalised, full schema preserved via
      :func:`reader.mean_per_instance_full`).

    Settings are auto-detected from the source dirs (no PICKS filter) so e.g.
    stage 2 emits CSV pairs for ``benders_FLLP``, ``benders_LEQ`` *and*
    ``benders_LEQRF`` even though PICKS only lists two of them.

    Stages with no FINAL reps are skipped silently.
    """
    # Collect: stage -> set(unique source dirs)
    stage_sources: dict[str, list[Path]] = {}
    for _display, _setting_grid, _src_grid, _setting_full, src_full in PICKS:
        srcs: list[Path] = _full_dirs(src_full)
        if not srcs:
            continue
        stage = _stage_name_from_path(srcs[0])
        if stage is None:
            continue
        bucket = stage_sources.setdefault(stage, [])
        for s in srcs:
            if s not in bucket:
                bucket.append(s)

    for stage, srcs in stage_sources.items():
        # Read every repetition once
        all_rows: list[pd.DataFrame] = []
        for s in srcs:
            if not s.exists():
                continue
            try:
                part = read_experiment(s)
            except Exception as e:
                print(f"  [warn] {stage}/{s.name}: {e}")
                continue
            if part.empty:
                continue
            tail = s.name.rstrip("0123456789")
            tag = s.name[len(tail):] or s.name
            part = part.copy()
            part["repetition"] = tag
            part["source_dir"] = s.name
            all_rows.append(part)

        if not all_rows:
            print(f"  [skip] {stage}: no FINAL repetitions on disk")
            continue

        combined = pd.concat(all_rows, ignore_index=True)
        if "settings" not in combined.columns:
            print(f"  [warn] {stage}: no 'settings' column in source data")
            continue
        # Drop globally-excluded CUT_ERROR instances from the merged full-run CSV.
        if _EXCLUDED and "instance" in combined.columns:
            combined = combined[~combined["instance"].astype(str).isin(_EXCLUDED)]

        out_dir = OUTPUT_ROOT / stage
        out_dir.mkdir(parents=True, exist_ok=True)

        for setting, sub in combined.groupby("settings", sort=True):
            safe = str(setting).replace("/", "_")
            compiled_path = out_dir / f"compiled_final_{safe}.csv"
            sub.to_csv(compiled_path, index=False)

            for_mean = sub.drop(columns=["repetition", "source_dir"],
                                errors="ignore")
            mean_df = mean_per_instance_full(for_mean)
            mean_path = out_dir / f"compiled_final_{safe}_mean.csv"
            mean_df.to_csv(mean_path, index=False)

            print(f"  [ok]   {stage}/{safe}: {len(sub)} raw / "
                  f"{len(mean_df)} mean rows")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    # Intermediate per-stage dumps of the FINAL (full-pglib) runs.
    print("\n=== per-stage FINAL dumps ===")
    _write_per_stage_final_dumps()
    # FINAL_GRID uses grid result dirs (31 instances, for paper headline table)
    # FINAL uses full-run dirs (174 instances, aggregated over 3 runs)
    for tag, grid in (("FINAL_GRID", True), ("FINAL", False)):
        print(f"\n=== {tag} ===")
        _write(tag, _collect(grid))


if __name__ == "__main__":
    main()
