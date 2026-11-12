"""Read TopologyOptimizer benchmark CSVs and compile experiment data.

Expected directory layout:
    check/results/<experiment_name>/<setting_name>/*.csv

Each CSV has columns: instance_file, settings, tal, time_*, instance, status, objective, n_*, ...
"""

import re
from pathlib import Path
import os
import pandas as pd

from selection import TIME_LIMIT

# ---------------------------------------------------------------------------
# Global CUT_ERROR exclusion list (full-run test set)
# ---------------------------------------------------------------------------
# Instances that hit a spurious ``CUT_ERROR`` in *any* run are listed, one
# instance name per line, in these files.  Every script that compiles the runs
# into a merged CSV / builds the comparison tables / performance profiles reads
# the appropriate list (via :func:`load_excluded_instances`) and drops exactly
# those instances, so the test set is consistent across all stages.
#
# NOTE: this global list covers the **full** test set only.  The grid test set is
# handled differently — see :func:`grid_cut_err_instances`.  (The two sets reuse
# instance *names* for *different* problems, so they must never share a list.)
CUT_ERR_EXCLUSIONS_FILE = Path(
    os.environ.get("TOPOLOGYOPTIMIZER_CUT_ERROR_LIST", Path(__file__).resolve().parents[2] / "results/CUT_ERR_INSTANCES.txt")        # full
)

# Relative tolerance for deciding that two OPTIMAL objective values *disagree*.
# Settings that solved the same instance must reach the same optimum, so a larger
# relative spread means some setting cut off the true optimum (a cut error).
# Below this it is treated as benign solver / numerical variation.
# Single knob shared by the exclusion scan (reader) and the majority-vote CutErr
# reporting (tables._majority_objective_errors).
CUT_ERR_OBJ_TOLERANCE = 1e-2


def drop_excluded(df: pd.DataFrame, excluded: set[str]) -> pd.DataFrame:
    """Return *df* without rows whose ``instance`` is in *excluded* (no-op if empty)."""
    if excluded and "instance" in df.columns:
        return df[~df["instance"].astype(str).isin(excluded)]
    return df


def load_excluded_instances(path: str | Path = CUT_ERR_EXCLUSIONS_FILE) -> set[str]:
    """Read the CUT_ERROR exclusion list; return the set of instance names.

    Blank lines and ``#`` comments are ignored.  Missing file → empty set.
    """
    path = Path(path)
    if not path.is_file():
        return set()
    out: set[str] = set()
    for line in path.read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            out.add(line)
    return out


def merge_specs(specs):
    """Merge ``(dir, settings)`` specs so each directory is scanned **once**.

    :func:`collect_cut_err_instances` works one directory at a time, so a spec
    list that names the same directory once per setting (the natural shape of a
    leaf script's ``SOURCES``) would only ever see one setting per pass — and
    objective *disagreement* between settings could never be detected.  Merging
    first restores it.  ``None`` (= all settings) wins over any restricted list.
    """
    merged: dict[str, list[str] | None] = {}
    for spec in specs:
        d, keep = spec if isinstance(spec, (tuple, list)) else (spec, None)
        k = str(d)
        if k not in merged:
            merged[k] = list(keep) if keep else None
        elif merged[k] is None or not keep:
            merged[k] = None
        else:
            merged[k] = sorted(set(merged[k]) | set(keep))
    return [(Path(k), v) for k, v in merged.items()]


def collect_cut_err_instances(specs, obj_thresh: float = CUT_ERR_OBJ_TOLERANCE,
                              use_obj_deviation: bool = True) -> set[str]:
    """Scan result directories for instances affected by a cut error.

    *specs* is an iterable of either a directory, or a ``(directory, settings)``
    pair where ``settings`` is the list of settings **actually reported** from
    that directory in the final across_all output (``None`` = all settings in the
    dir).  Restricting to the reported settings matters: an instance must not be
    thrown away because two settings we never publish disagreed on it.

    Two independent symptoms are detected, in a single pass per directory:

    1. **Explicit** — a reported setting has ``status == "CUT_ERROR"``.
    2. **Objective disagreement** (only when *use_obj_deviation*) — among the
       reported settings that returned ``OPTIMAL`` for an instance, the objective
       values disagree by more than *obj_thresh* (relative).  They solved the
       same problem, so they must agree; a divergent value means some setting cut
       off the true optimum.  Rather than trusting any one setting (the baseline
       is *not* privileged), the whole instance is deemed untrustworthy.

    ``use_obj_deviation=False`` restricts detection to *proven* cut errors
    (symptom 1 only).  This is used for the **grid** test set: it holds only ~31
    instances and is probed by ~11 grid searches with up to 44 settings each, so
    unioning objective-deviation flags across them would wipe out nearly the whole
    set even though a typical disagreement is one deviating setting out of 44.

    Returns the set of ``instance`` names to exclude from the test set.
    """
    found: set[str] = set()
    for spec in specs:
        if isinstance(spec, (tuple, list)):
            d, keep = spec
        else:
            d, keep = spec, None
        d = Path(d)
        if not d.is_dir():
            continue
        try:
            df = read_experiment(d)
        except FileNotFoundError:
            continue
        if keep is not None and "settings" in df.columns:
            df = df[df["settings"].isin(list(keep))]
        if df.empty:
            continue

        # 1. explicit CUT_ERROR status
        if {"status", "instance"}.issubset(df.columns):
            found |= set(df.loc[df["status"] == "CUT_ERROR", "instance"].astype(str))

        # 2. objective disagreement among OPTIMAL runs (opt-in)
        if not use_obj_deviation:
            continue
        if not {"status", "instance", "objective", "settings"}.issubset(df.columns):
            continue
        opt = df[df["status"] == "OPTIMAL"].copy()
        opt["objective"] = pd.to_numeric(opt["objective"], errors="coerce")
        opt = opt.dropna(subset=["objective"])
        if opt.empty:
            continue
        piv = opt.pivot_table(index="instance", columns="settings",
                              values="objective", aggfunc="first")
        if piv.shape[1] < 2:
            continue
        spread = (piv.max(axis=1) - piv.min(axis=1)) / piv.abs().max(axis=1).clip(lower=1e-10)
        found |= set(spread.index[spread > obj_thresh].astype(str))
    return found


def grid_cut_err_instances(src, settings=None) -> set[str]:
    """Proven cut-error instances of a **single grid experiment**, computed live.

    The grid test set has only ~31 instances but is probed by ~11 grid searches,
    so a *global* exclusion list (unioned across them) removes almost everything:
    a defect in the cut-filtering grid would shrink the unrelated
    partial-decomposition grid.  Grid exclusions are therefore scoped **per
    experiment** — each grid table drops only its own defects — and computed on
    the fly, so no grid exclusion file exists.

    Only *proven* ``CUT_ERROR`` status counts here (no objective-deviation
    flagging): with up to 44 settings per grid search, a single deviating setting
    would otherwise condemn the whole instance.
    """
    return collect_cut_err_instances([(src, settings)], use_obj_deviation=False)


def update_cut_err_file(dirs, path: str | Path = CUT_ERR_EXCLUSIONS_FILE,
                        reset: bool = False,
                        use_obj_deviation: bool = True) -> set[str]:
    """Write the CUT_ERROR instances found in *dirs* to the exclusion file.

    By default the file only grows (an instance that *ever* cut-errored stays
    excluded).  Pass ``reset=True`` for a clean rebuild — correct when the caller
    rescans *all* in-scope sources anyway, and required for a change in scanning
    scope (e.g. restricting to reported settings) to actually take effect.

    Written sorted with a short header.  Returns the full resulting set.
    """
    path = Path(path)
    result = collect_cut_err_instances(dirs, use_obj_deviation=use_obj_deviation)
    if not reset:
        result |= load_excluded_instances(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    header = (
        "# Instances excluded from the FULL test set because they hit a CUT_ERROR\n"
        "# in at least one full-run test.  One instance name per line.  Consumed by\n"
        "# reader.load_excluded_instances() in every full-run compile/table/perf script.\n"
    )
    path.write_text(header + "\n".join(sorted(result)) + ("\n" if result else ""))
    return result


def read_setting(setting_dir: str | Path) -> pd.DataFrame:
    """Read all CSVs from a single setting directory and combine them."""
    setting_dir = Path(setting_dir)
    csvs = sorted(setting_dir.glob("*.csv"))
    if not csvs:
        raise FileNotFoundError(f"No CSV files found in {setting_dir}")
    frames = [pd.read_csv(f) for f in csvs]
    return pd.concat(frames, ignore_index=True)


def read_experiment(experiment_dir: str | Path) -> pd.DataFrame:
    """Read all settings from an experiment directory.

    Each subdirectory of experiment_dir is treated as a setting.
    Returns a combined DataFrame with a 'setting' column derived from 'settings' CSV column.
    """
    experiment_dir = Path(experiment_dir)
    setting_dirs = sorted(
        d for d in experiment_dir.iterdir()
        if d.is_dir() and not d.name.startswith(".")
        and not d.name.endswith("_infeasible")
        and d.name not in {"log_metrics", "old_stats"}
    )
    if not setting_dirs:
        raise FileNotFoundError(f"No setting subdirectories in {experiment_dir}")

    frames = []
    for d in setting_dirs:
        try:
            df = read_setting(d)
            # Folder name is authoritative for the setting label so that
            # symlinked dirs (e.g. cf100_sp02_010 -> cf100_sp02_005) are
            # ingested under the linked name, not the original CSV value.
            df["settings"] = d.name
            frames.append(df)
        except FileNotFoundError:
            continue

    if not frames:
        raise FileNotFoundError(f"No CSV data found in {experiment_dir}")
    out = pd.concat(frames, ignore_index=True)

    # Sanity-check: duplicate (instance, [tal,] settings) rows would be
    # silently collapsed by ``pivot_table(aggfunc='first')`` downstream.
    key_cols = ["instance", "settings"]
    if "tal" in out.columns:
        key_cols.insert(1, "tal")
    dup_mask = out.duplicated(subset=key_cols, keep=False)
    if dup_mask.any():
        n_dup = int(dup_mask.sum())
        sample = out.loc[dup_mask, key_cols].drop_duplicates().head(5).to_dict("records")
        print(f"WARNING: {n_dup} duplicate rows in {experiment_dir} for keys {key_cols}; sample: {sample}")
    return out


# ---------------------------------------------------------------------------
# Multiple-repetition ("full run") support
# ---------------------------------------------------------------------------
# A "full run" of an experiment is dispatched several times; the repetitions
# historically lived in sibling directories sharing a common base name and a
# trailing repetition index, e.g. ``2-SP_solving1`` / ``2-SP_solving2`` /
# ``2-SP_solving3``.  The helpers below generalise that fixed 1/2/3 convention
# to an arbitrary number of repetitions and let downstream code merge + average
# them per instance.

def discover_reps(base_dir: str | Path) -> list[Path]:
    """Discover the repetition directories for a full-run experiment base name.

    Given a base path such as ``.../check/results/9-cutfiltering-full`` this
    returns every existing repetition directory ``<base>1``, ``<base>2``, …
    sorted by repetition index.  If no numbered siblings exist but the
    unsuffixed ``<base>`` directory does, ``[<base>]`` is returned (a single,
    un-replicated run).  Returns ``[]`` if nothing matches.

    Numbered repetitions deliberately take precedence over the unsuffixed base
    directory: this reproduces the historical behaviour where ``2-SP_solving1``
    /``2/``3`` are the three runs and the bare ``2-SP_solving`` directory (a
    separate single run) is *not* pulled into the average.
    """
    base_dir = Path(base_dir)
    parent, name = base_dir.parent, base_dir.name
    numbered: list[tuple[int, Path]] = []
    if parent.is_dir():
        for d in parent.iterdir():
            if not d.is_dir():
                continue
            m = re.fullmatch(re.escape(name) + r"(\d+)", d.name)
            if m:
                numbered.append((int(m.group(1)), d))
    if numbered:
        numbered.sort(key=lambda t: t[0])
        return [p for _, p in numbered]
    if base_dir.is_dir():
        return [base_dir]
    return []


def consolidate_status(statuses) -> str:
    """Combine the per-repetition statuses of one instance into a single status.

    Rules:

    * all repetitions ``OPTIMAL``              → ``"OPTIMAL"``;
    * all repetitions share one status         → **that** status (so an instance
      that is ``INFEASIBLE`` / ``CUT_ERROR`` in *every* run keeps that label and
      is subsequently excluded by :func:`selection.build_selection`, exactly as
      a single-run comparison would);
    * mixed (not all OPTIMAL, statuses differ) → ``"TIME_LIMIT"`` (treated as an
      unsolved instance).

    This deliberately does **not** collapse every non-OPTIMAL status to
    ``TIME_LIMIT``: doing so previously re-labelled unanimous ``INFEASIBLE`` runs
    as ``TIME_LIMIT``, which then leaked into the ``TL`` column of the comparison
    tables instead of being dropped as broken.
    """
    vals = set(statuses)
    if vals == {"OPTIMAL"}:
        return "OPTIMAL"
    if len(vals) == 1:
        return next(iter(vals))
    return "TIME_LIMIT"


def mean_per_instance_full(df: pd.DataFrame, target: str = "time_solve") -> pd.DataFrame:
    """Aggregate multiple runs per ``(instance[, tal], settings)`` by mean.

    Numeric columns are arithmetic-mean averaged across repetitions.  Non-numeric
    columns take the first value (assumed identical across runs for the same
    instance/setting).  ``target`` (default ``time_solve``) is penalised to
    :data:`TIME_LIMIT` for any non-OPTIMAL row *before* averaging, so a setting
    that only solved an instance in some runs is penalised appropriately.  The
    aggregated ``status`` is derived by :func:`consolidate_status` (unanimous
    status preserved, e.g. all-``INFEASIBLE`` stays ``INFEASIBLE``; mixed →
    ``"TIME_LIMIT"``).

    The full CSV schema is preserved so the result can be fed straight into
    :func:`read_experiment`-style consumers.
    """
    if df.empty:
        return df.copy()

    df = df.copy()
    if target in df.columns and "status" in df.columns:
        df.loc[df["status"] != "OPTIMAL", target] = TIME_LIMIT

    group_cols = ["instance"]
    if "tal" in df.columns:
        group_cols.append("tal")
    if "settings" in df.columns:
        group_cols.append("settings")

    numeric_cols = [
        c for c in df.select_dtypes(include="number").columns
        if c not in group_cols
    ]
    other_cols = [
        c for c in df.columns
        if c not in group_cols and c not in numeric_cols and c != "status"
    ]

    agg_spec: dict = {c: "mean" for c in numeric_cols}
    agg_spec.update({c: "first" for c in other_cols})
    if "status" in df.columns:
        agg_spec["status"] = consolidate_status

    return df.groupby(group_cols, as_index=False).agg(agg_spec)


def read_runs(dirs, average: bool = True, target: str = "time_solve") -> pd.DataFrame:
    """Read one or more repetition directories and optionally average them.

    Each directory in *dirs* is read via :func:`read_experiment` and tagged with
    a ``repetition`` column (the trailing digits of the directory name, or the
    full name if it has none) and a ``source_dir`` column.  Directories that do
    not exist or contain no CSV data are skipped with a warning.

    If *average* is True (the default) the concatenated rows are collapsed to one
    row per ``(instance[, tal], settings)`` via :func:`mean_per_instance_full`
    (the ``repetition`` / ``source_dir`` tags are dropped first).  If False, the
    raw tagged concatenation is returned.

    CUT_ERROR exclusion is applied by the caller (``run_full_comparison`` /
    ``across_all``) with the test-set-appropriate list, *not* here — the grid and
    full test sets reuse instance names for different problems, so their
    exclusion lists must never be mixed.
    """
    dirs = [Path(d) for d in dirs]
    parts: list[pd.DataFrame] = []
    for d in dirs:
        if not d.is_dir():
            continue
        try:
            part = read_experiment(d)
        except FileNotFoundError as e:
            print(f"  [warn] read_runs: {e}")
            continue
        if part.empty:
            continue
        tail = d.name.rstrip("0123456789")
        rep = d.name[len(tail):] or d.name
        part = part.copy()
        part["repetition"] = rep
        part["source_dir"] = d.name
        parts.append(part)

    if not parts:
        raise FileNotFoundError(f"No CSV data found in any of: {[str(d) for d in dirs]}")

    combined = pd.concat(parts, ignore_index=True)
    if not average:
        return combined

    for_mean = combined.drop(columns=["repetition", "source_dir"], errors="ignore")
    return mean_per_instance_full(for_mean, target=target)


def read_setting_sources(spec, target: str = "time_solve",
                         average: bool = True) -> pd.DataFrame:
    """Assemble one comparison frame from **per-setting** sources.

    *spec* maps ``setting_name -> sources``, where ``sources`` is either

    * an iterable of run directories (the setting's folder is named like the
      setting), or
    * a ``(dirs, source_folder)`` tuple, when the real runs live under a
      **different folder name** than the label we report it under.  This happens
      because settings were sometimes filed under a renamed folder (e.g. the
      runs for ``div_decay_5_005`` actually live in ``div_decay_5`` folders).

    Each setting is read from *its own* sources and aggregated with the current
    :func:`mean_per_instance_full`; the per-setting frames are then concatenated
    into a single comparison frame.  A mismatch between the folder name and the
    ``settings`` column recorded inside the CSVs is reported, so silent renames
    surface instead of being trusted blindly.

    Why this exists: some "full run" directories were hand-assembled by copying
    **pre-computed mean CSVs** in from other experiments (and other stage repos).
    Reading those freezes the aggregation at whatever code produced them, and
    silently drops settings whose reps live elsewhere.  Pointing at the original
    per-run CSVs instead means the mean is always recomputed with the current
    logic (status consolidation, TL penalisation, …) and nothing is copied.

    Settings with no readable source are skipped with a warning.
    """
    frames: list[pd.DataFrame] = []
    for setting, sources in dict(spec).items():
        # Optional (dirs, source_folder) form for renamed folders.
        if isinstance(sources, tuple) and len(sources) == 2 and isinstance(sources[1], str):
            dirs, folder = sources
        else:
            dirs, folder = sources, setting

        parts: list[pd.DataFrame] = []
        inner: set[str] = set()
        for d in dirs:
            sub = Path(d) / folder
            if not sub.is_dir():
                continue
            try:
                part = read_setting(sub)
            except FileNotFoundError:
                continue
            part = part.copy()
            if "settings" in part.columns:
                inner |= set(part["settings"].astype(str).unique())
            # We label by the *reported* setting name, not the folder.
            part["settings"] = setting
            part["source_dir"] = Path(d).name
            parts.append(part)

        # Surface renames: what the CSVs call themselves vs what we report.
        if inner and inner != {setting}:
            print(f"  [note] '{setting}' read from folder '{folder}' "
                  f"(CSV settings column says {sorted(inner)})")
        if not parts:
            print(f"  [warn] read_setting_sources: no data for setting '{setting}'")
            continue
        combined = pd.concat(parts, ignore_index=True)
        n_reps = combined["source_dir"].nunique()
        combined = combined.drop(columns=["source_dir"], errors="ignore")
        if average:
            combined = mean_per_instance_full(combined, target=target)
        print(f"  [setting] {setting}: {n_reps} run(s) -> {len(combined)} rows")
        frames.append(combined)

    if not frames:
        raise FileNotFoundError("read_setting_sources: no data for any setting")
    return pd.concat(frames, ignore_index=True)


def clean_instance_name(name: str) -> str:
    """Strip common pglib prefixes and .m suffix for display."""
    name = name.replace("pglib_opf_case", "").replace("pglib_scopf_case", "").replace(".m", "")
    return name


def instance_sort_key(name: str) -> tuple:
    """Extract numeric size from instance name for sorting (e.g., '118_ieee_01' → 118)."""
    import re
    m = re.search(r"(\d+)", clean_instance_name(name))
    return (int(m.group(1)) if m else 0, name)


# Statuses that indicate a truly broken run (exclude entire instance group)
_EXCLUDE_STATUSES = {"CUT_ERROR", "INFEASIBLE", "INFEASIBLE_OR_UNBOUNDED", "DID_NOT_RUN"}


def cleanup_data(df: pd.DataFrame, exclude_all_timelimit: bool = True,
                 baseline: str | None = None) -> tuple[pd.DataFrame, set]:
    """Remove broken rows and optionally resource-limited groups.

    If *baseline* is given, drops entire instance groups where the baseline
    setting has a bad status.  For other settings, drops individual bad rows
    but records which (instance[, tal], setting) tuples had errors so the
    table can mark them with a symbol (e.g. '--‡').

    Returns
    -------
    (cleaned_df, error_keys) where error_keys is a set of tuples
    (instance, [tal,] setting) that had an error status (non-baseline).
    """
    has_tal = "tal" in df.columns
    group_cols = ["instance", "tal"] if has_tal else ["instance"]

    bad_mask = df["status"].isin(_EXCLUDE_STATUSES)
    error_keys: set = set()

    if bad_mask.any():
        bad_rows = df[bad_mask]

        # Drop entire groups where baseline has a bad status
        if baseline is not None:
            baseline_bad = bad_rows[bad_rows["settings"] == baseline]
            if len(baseline_bad) > 0:
                bad_keys = baseline_bad[group_cols].drop_duplicates()
                before = len(df)
                df = df.merge(bad_keys, on=group_cols, how="left", indicator=True)
                df = df[df["_merge"] == "left_only"].drop(columns="_merge").reset_index(drop=True)
                print(f"cleanup: removed {before - len(df)} rows (baseline '{baseline}' had bad status)")
                # Recompute bad_mask after group removal
                bad_mask = df["status"].isin(_EXCLUDE_STATUSES)
                bad_rows = df[bad_mask]

        # Record error keys for non-baseline bad rows (for table annotation)
        for _, row in bad_rows.iterrows():
            key = tuple(row[c] for c in group_cols) + (row["settings"],)
            error_keys.add(key)

        # Drop individual bad rows
        if bad_mask.any():
            before = len(df)
            df = df[~bad_mask].reset_index(drop=True)
            print(f"cleanup: removed {before - len(df)} individual rows with bad status")

    # Drop groups where ALL settings hit TIME_LIMIT / MEMORY_LIMIT
    if exclude_all_timelimit and len(df) > 0:
        tl_statuses = {"TIME_LIMIT", "MEMORY_LIMIT"}
        grouped = df.groupby(group_cols)["status"].apply(
            lambda s: s.isin(tl_statuses).all()
        )
        all_tl_keys = grouped[grouped].reset_index()[group_cols]
        if len(all_tl_keys) > 0:
            before = len(df)
            df = df.merge(all_tl_keys, on=group_cols, how="left", indicator=True)
            df = df[df["_merge"] == "left_only"].drop(columns="_merge")
            print(f"cleanup: removed {before - len(df)} rows where all settings hit time/mem limit")

    return df.reset_index(drop=True), error_keys


def compile_experiment_csv(experiment_dir: str | Path, output: str | Path | None = None) -> pd.DataFrame:
    """Read experiment and write compiled CSV. Returns the DataFrame."""
    experiment_dir = Path(experiment_dir)
    df = read_experiment(experiment_dir)
    if output is None:
        output = experiment_dir / "compiled.csv"
    df.to_csv(output, index=False)
    print(f"Compiled {len(df)} rows to {output}")
    return df


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python reader.py <experiment_dir> [output.csv]")
        sys.exit(1)
    out = sys.argv[2] if len(sys.argv) > 2 else None
    compile_experiment_csv(sys.argv[1], out)
