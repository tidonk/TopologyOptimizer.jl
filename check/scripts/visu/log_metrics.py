"""Per-instance log-derived metrics and wall-time totals.

For every ``(instance, setting [, tal])`` row in the experiment CSVs computes
five metrics. The first four come from the Gurobi master log and include
**the full warmstart + main-run timeline** (warmstart sub-MIPs accumulate into
the wall-time offset, so an RSCOTS stage that improves over the initial SCOPF
solution counts as a non-trivial improvement). The fifth comes from the CSV:

* ``time_to_first``             — wall time of the first feasible incumbent
  across all stages.
* ``time_to_first_improving``  — wall time of the first incumbent strictly
  better than the very first feasible incumbent. For runs where the initial
  incumbent is already optimal this is scored as ``0``.
* ``time_to_optimal``           — wall time of the first incumbent that matches
  the final best objective within a small relative tolerance.
* ``pdi``                       — primal-dual integral (left-continuous step,
  trajectory extended to ``total_time``).  The **primal reference** is the
  per-instance SCOPF objective (``scopf_threshold``): any recorded primal
  worse than that threshold — including the leading no-incumbent interval —
  is replaced by ``scopf_threshold``, so all settings start the PDI from the
  same quality reference.  For warmstarted runs the **dual during warmstart
  sub-MIP stages** (SCOPF, RSCOTS, …) is zeroed out because those duals are
  not valid lower bounds for the full main-solve problem.
* ``time_total``                — total wall time including the warmstart
  (``time_total`` column if present, else ``time_solve + time_warmstart``).

Emits, under ``<experiment_dir>/log_metrics/<metric>/``:

* ``comparison.{md,tex}``           — long table, cleaned (drops all-TL groups).
* ``comparison_all.{md,tex}``       — long table, all rows.
* ``comparison_short.{md,tex}``     — solved-by-all short table.
* ``comparison_short_all.{md,tex}`` — all-instance short table.

Plus ``log_metrics/log_metrics.csv`` with raw per-row values and diagnostics.

Usage
-----
    python check/scripts/visu/log_metrics.py <experiment_dir> [baseline]
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from gurobi_log import (
    compute_pd_integral, compute_primal_integral, compute_time_to_first,
    compute_time_to_first_improving, compute_time_to_optimal,
    parse_gurobi_log, read_log,
)
from reader import cleanup_data, read_experiment
from tables import (
    comparison_short_all_table, comparison_short_table,
    comparison_short_to_latex, comparison_short_to_markdown,
    comparison_table, to_latex, to_markdown,
)


_LOG_RE = "GRB_MP_*.out"
_TS_RE = re.compile(r"_(\d{4}-\d{2}-\d{2}_\d{2}_\d{2}_\d{2})$")


def drop_incomplete_instances(
    df: pd.DataFrame,
    include_settings: list[str] | None = None,
) -> tuple[pd.DataFrame, list]:
    """Remove instances where any relevant setting is missing or has CUT_ERROR.

    An instance is *incomplete* if, for the settings considered
    (``include_settings``, or all present settings when not specified), at
    least one setting has no data row (NaN in the status pivot) or a
    ``CUT_ERROR`` status.

    Returns ``(filtered_df, dropped)`` where ``dropped`` is the list of
    instance identifiers that were removed (strings or (instance, tal) tuples).
    """
    has_tal = "tal" in df.columns
    df2 = df.copy()
    if has_tal:
        df2["tal"] = df2["tal"].fillna(0).astype(int)
    group_cols = ["instance", "tal"] if has_tal else ["instance"]

    all_settings = sorted(df2["settings"].unique())
    check_settings = [s for s in (include_settings or all_settings)
                      if s in set(all_settings)]

    status_piv = df2.pivot_table(
        index=group_cols, columns="settings", values="status", aggfunc="first",
    ).reindex(columns=check_settings)

    incomplete = (status_piv.isna().any(axis=1) |
                  (status_piv == "CUT_ERROR").any(axis=1))
    bad_indices = incomplete[incomplete].index

    if len(bad_indices) == 0:
        return df2, []

    if has_tal:
        bad_set = set(bad_indices.tolist())
        mask = ~pd.MultiIndex.from_arrays([df2["instance"], df2["tal"]]).isin(bad_set)
    else:
        bad_set = set(bad_indices.tolist())
        mask = ~df2["instance"].isin(bad_set)

    dropped = sorted(bad_set)
    return df2[mask].reset_index(drop=True), dropped

_LOG_METRICS = (
    "time_to_first",
    "time_to_first_improving",
    "time_to_optimal",
    "pdi",
    "primal_integral",
)
_CSV_METRICS = ("time_total",)
METRICS = _LOG_METRICS + _CSV_METRICS


def _instance_from_log(log_path: Path, setting: str) -> str:
    """Recover the instance name from a ``GRB_MP_<inst>_<setting>_<ts>.out`` filename."""
    stem = log_path.stem.removeprefix("GRB_MP_")
    ts_m = _TS_RE.search(stem)
    if ts_m:
        stem = stem[: ts_m.start()]
    suffix = f"_{setting}"
    if stem.endswith(suffix):
        stem = stem[: -len(suffix)]
    return stem


def parse_log_metrics(log_path: Path, scopf_obj: float = float("nan"),
                      use_csvs_if_possible: bool = True) -> dict:
    """Parse one Gurobi log and return full-timeline metrics + diagnostics.

    Loads via :func:`gurobi_log.read_log` (so a previously-exported
    ``logs_parsed/<stem>.csv`` is used when available; otherwise the raw log
    is parsed and the CSV is written for next time).  ``scopf_obj`` is the
    per-instance non-trivial threshold (see :func:`_compute_metrics`); if not
    provided, the first feasible incumbent of *this* log is used.

    When Gurobi solves the instance in presolve / at the root LP the
    node-table is empty and the trajectory-derived metrics would all be NaN.
    In that case (``optimal=True``) we fall back to ``total_time``.

    When the SCOPF threshold already matches the final best-known objective
    both ``time_to_first_improving`` and ``time_to_optimal`` are scored as
    ``0`` (no strictly-better solution can appear).
    """
    data = read_log(log_path, use_csvs_if_possible=use_csvs_if_possible)
    return _compute_metrics(data, scopf_obj)
def _first_feasible(primals: np.ndarray) -> float:
    """First non-NaN primal value (= first feasible incumbent), or NaN."""
    finite = primals[~np.isnan(primals)]
    return float(finite[0]) if len(finite) else float("nan")


def _fair_primal_dual(
    data: dict, threshold: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build fair primal/dual trajectories for PDI / primal-integral.

    Returns ``(times, primals_fair, duals_fair)`` where:

    * The leading no-incumbent interval (NaN primals) is filled with
      ``max(threshold, first_real_main_primal)`` — the same anchor logic used
      by ``_anchor_value`` in the aggregated PD plot.  Warmstarted runs
      typically load the SCOPF result as their first main-solve incumbent, so
      their anchor equals ``threshold``.  ``no_warmstart`` may start with a
      worse incumbent, so its anchor is its actual first primal — it is *not*
      artificially improved to the SCOPF quality it did not earn.
    * Primals during warmstart sub-MIP stages are replaced with the same
      anchor, because sub-MIP node-table primals are objective values of a
      *different* optimization problem (SCOPF / RSCOTS) and are not comparable
      to the main MP's objective.
    * Duals during warmstart sub-MIP stages are zeroed because those duals
      are not valid lower bounds for the main problem.
    """
    times   = np.asarray(data["times"],   dtype=float)
    primals = np.asarray(data["primals"], dtype=float).copy()
    duals   = np.asarray(data["duals"],   dtype=float).copy()

    main_start = float(data["restarts"][-1]) if data["restarts"] else 0.0
    in_warmstart = times < main_start

    # First real incumbent from the main solve (ignores warmstart sub-MIP rows).
    main_primals_nz = primals[~in_warmstart & ~np.isnan(primals)]
    first_real = float(main_primals_nz[0]) if len(main_primals_nz) > 0 else np.nan

    if not np.isnan(threshold):
        anchor = threshold if np.isnan(first_real) else max(threshold, first_real)
    else:
        anchor = first_real

    if not np.isnan(anchor):
        # Fill NaN primals (leading no-incumbent interval) with anchor.
        primals = np.where(np.isnan(primals), anchor, primals)
        # Replace sub-MIP window with anchor (those primals are on a different
        # objective scale and must not influence the integral).
        primals[in_warmstart] = anchor

    duals[in_warmstart] = 0.0
    return times, primals, duals


def _compute_metrics(data: dict, scopf_obj: float) -> dict:
    """Derive the five per-row metrics from a parsed log + per-instance SCOPF.

    ``scopf_obj`` is the **per-instance** non-trivial threshold (typically the
    first feasible incumbent produced by the SCOPF stage of a warmstart run on
    the same instance). It is used uniformly across all settings — including
    the non-warmstart baseline — so the fairness of ``time_to_first_improving``
    is preserved.

    ``time_to_first_improving`` searches the **full** timeline (including
    warmstart sub-MIP stages): if a warmstart stage produces an incumbent
    strictly better than the threshold, the timer is stopped there.

    Both ``pdi`` and ``primal_integral`` use the **same** fair primal/dual
    trajectory built by :func:`_fair_primal_dual`, so they are directly
    comparable across warmstart and no-warmstart settings.
    """
    t_f   = compute_time_to_first(data["times"], data["primals"])
    initial_obj = _first_feasible(data["primals"])
    threshold = scopf_obj if not np.isnan(scopf_obj) else initial_obj
    t_nt  = compute_time_to_first_improving(
        data["times"], data["primals"], threshold)
    t_opt = compute_time_to_optimal(
        data["times"], data["primals"], data["objective"])

    times_np, primals_fair, duals_fair = _fair_primal_dual(data, threshold)
    pdi = compute_pd_integral(
        times_np, primals_fair, duals_fair,
        total_time=data["total_time"])
    primal_int = compute_primal_integral(
        times_np, primals_fair, data["objective"],
        total_time=data["total_time"])

    obj = data["objective"]
    threshold_was_optimal = (not np.isnan(threshold) and not np.isnan(obj)
                             and abs(threshold - obj)
                                 <= max(1e-6, 1e-4 * abs(obj)))

    if data["optimal"] and len(data["times"]) == 0:
        # Solved at root/presolve — no node-table emitted.
        t_opt = float(data["total_time"])
        if np.isnan(t_f):
            t_f = float(data["total_time"])
        if np.isnan(threshold) or not threshold_was_optimal:
            t_nt = float(data["total_time"])

    if threshold_was_optimal:
        # SCOPF threshold already matches the final best-known objective — no
        # strictly-better solution can ever appear; score 0 (best possible)
        # instead of TL.
        if np.isnan(t_nt):
            t_nt = 0.0
        t_opt = float(t_f) if not np.isnan(t_f) else 0.0

    return {
        "time_to_first":            t_f,
        "time_to_first_improving": t_nt,
        "time_to_optimal":          t_opt,
        "pdi":                      pdi,
        "primal_integral":          primal_int,
        "log_total_time":           data["total_time"],
        "log_objective":            data["objective"],
        "log_warmstart_objective":  data["warmstart_objective"],
        "log_optimal":              bool(data["optimal"]),
        "n_traj_points":            int(len(data["times"])),
        "scopf_threshold":          threshold,
    }


def _collect_logs(experiment_dir: Path) -> dict[tuple[str, str], Path]:
    """Map ``(instance, setting) -> newest log file`` across the experiment."""
    _SKIP = {"log_metrics", "old_stats"}
    latest: dict[tuple[str, str], tuple[str, Path]] = {}
    for setting_dir in sorted(p for p in experiment_dir.iterdir()
                              if p.is_dir() and p.name not in _SKIP
                              and not p.name.startswith(".")):
        log_dir = setting_dir / "logs"
        if not log_dir.is_dir():
            continue
        for log_file in log_dir.glob(_LOG_RE):
            setting = setting_dir.name
            inst = _instance_from_log(log_file, setting)
            ts_m = _TS_RE.search(log_file.stem)
            ts = ts_m.group(1) if ts_m else ""
            key = (inst, setting)
            if key not in latest or ts > latest[key][0]:
                latest[key] = (ts, log_file)
    return {k: v[1] for k, v in latest.items()}


def build_log_cache(experiment_dir: Path,
                    use_csvs_if_possible: bool = True) -> tuple[dict, dict]:
    """Build the log-file map and raw parse cache for *experiment_dir*.

    Returns ``(logs, raw_cache)`` where:

    * ``logs``      — ``{(instance, setting): Path}`` newest log per row.
    * ``raw_cache`` — ``{Path: read_log(..., use_csvs_if_possible=...)}``
      parsed log data, one entry per unique log file.  Each entry is loaded
      from the sibling ``logs_parsed/<stem>.csv`` when available; otherwise
      the raw log is parsed and the CSV is written as a side effect.
      Pass ``use_csvs_if_possible=False`` to force a full re-parse that
      also refreshes the on-disk CSVs.

    Call this once and pass both to :func:`collect_log_metrics` and to
    :func:`~tables.run_full_comparison` (via its ``log_raw_cache`` /
    ``log_files_map`` params) to avoid parsing every log twice.
    """
    logs = _collect_logs(experiment_dir)
    raw_cache: dict[Path, dict] = {
        log: read_log(log, use_csvs_if_possible=use_csvs_if_possible)
        for log in set(logs.values())
    }
    return logs, raw_cache


def _scopf_threshold_per_instance(
    logs: dict[tuple[str, str], Path], raw_cache: dict[Path, dict],
) -> dict[str, float]:
    """Per-instance SCOPF threshold (lowest sub-MIP-loaded first incumbent).

    Across all logs that **loaded a user MIP start** (= ran a warmstart
    sub-MIP), pick the lowest first feasible primal value per instance.  This
    is the same reference used for both ``time_to_first_improving`` and the
    fair primal/dual trajectory in :func:`_fair_primal_dual`.
    """
    out: dict[str, float] = {}
    for (inst, _), log in logs.items():
        data = raw_cache[log]
        # Use the SCOPF-warmstart result explicitly: this is the value that the
        # RSCOTS warmstart loaded as its MIP start (first "Loaded user MIP start"
        # in the log = section 2).  It represents the starting quality of the
        # warmstart pipeline and is the fair anchor for no_warmstart comparisons.
        ws = data.get("scopf_objective", np.nan)
        if np.isnan(ws):
            continue
        if ws < out.get(inst, np.inf):
            out[inst] = ws
    return out


_EMPTY_ROW = {
    "time_to_first": np.nan, "time_to_first_improving": np.nan,
    "time_to_optimal": np.nan, "pdi": np.nan, "primal_integral": np.nan,
    "log_total_time": np.nan, "log_objective": np.nan,
    "log_warmstart_objective": np.nan,
    "log_optimal": False, "n_traj_points": 0,
    "scopf_threshold": np.nan,
}


def _per_row_metrics(
    df: pd.DataFrame, logs: dict, raw_cache: dict,
    scopf_per_inst: dict[str, float],
) -> pd.DataFrame:
    """Compute per-row metrics, caching by (log_file, threshold)."""
    cache: dict[tuple[Path, float], dict] = {}
    rows = []
    for _, row in df.iterrows():
        log = logs.get((row["instance"], row["settings"]))
        if log is None:
            rows.append(_EMPTY_ROW)
            continue
        scopf = scopf_per_inst.get(row["instance"], np.nan)
        ck = (log, scopf)
        if ck not in cache:
            cache[ck] = _compute_metrics(raw_cache[log], scopf)
        rows.append(cache[ck])
    return pd.DataFrame(rows, index=df.index)


def _add_time_total(out: pd.DataFrame) -> pd.DataFrame:
    """Ensure ``time_total`` column exists, falling back to solve + warmstart."""
    if "time_total" in out.columns:
        total = pd.to_numeric(out["time_total"], errors="coerce")
    else:
        total = pd.Series(np.nan, index=out.index)
    if total.isna().any() and {"time_solve", "time_warmstart"}.issubset(out.columns):
        fallback = (pd.to_numeric(out["time_solve"], errors="coerce").fillna(0)
                    + pd.to_numeric(out["time_warmstart"], errors="coerce").fillna(0))
        total = total.where(total.notna(), fallback)
    out["time_total"] = total
    return out


def collect_log_metrics(
    experiment_dir: Path,
    raw_cache: dict | None = None,
    logs: dict | None = None,
    use_csvs_if_possible: bool = True,
) -> pd.DataFrame:
    """Join experiment CSV rows with parsed (full-timeline) log metrics.

    Adds the CSV-derived ``time_total`` column (preferring an existing
    ``time_total`` value, else ``time_solve + time_warmstart``).

    The non-trivial threshold for every row is the per-instance SCOPF
    objective — derived as the first feasible incumbent of the first warmstart
    sub-MIP from *any* warmstart-enabled log on the same instance. This makes
    ``time_to_first_improving`` fair across the non-warmstart baseline and
    every warmstart variant (they all have to beat the same SCOPF reference).

    Pass pre-built ``logs`` and ``raw_cache`` (from :func:`build_log_cache`)
    to avoid re-parsing logs that were already read for another purpose
    (e.g. the ``gap_map`` in :func:`~tables.run_full_comparison`).
    """
    df = read_experiment(experiment_dir)
    if logs is None:
        logs = _collect_logs(experiment_dir)
    if raw_cache is None:
        raw_cache = {log: read_log(log, use_csvs_if_possible=use_csvs_if_possible)
                     for log in set(logs.values())}

    scopf_per_inst = _scopf_threshold_per_instance(logs, raw_cache)
    metrics_df = _per_row_metrics(df, logs, raw_cache, scopf_per_inst)
    out = pd.concat([df.reset_index(drop=True),
                     metrics_df.reset_index(drop=True)], axis=1)
    return _add_time_total(out)


def write_metrics_csv(df: pd.DataFrame, out_path: Path) -> Path:
    """Write per-(instance, setting [, tal]) log-metrics CSV."""
    cols = ["instance", "settings"]
    if "tal" in df.columns:
        cols.append("tal")
    cols += ["status", "objective", "time_solve", "time_total",
             "time_to_first", "time_to_first_improving", "time_to_optimal", "pdi",
             "scopf_threshold",
             "log_total_time", "log_objective", "log_warmstart_objective",
             "log_optimal", "n_traj_points"]
    keep = [c for c in cols if c in df.columns]
    df[keep].to_csv(out_path, index=False, quoting=csv.QUOTE_MINIMAL)
    return out_path


def _synthetic_status(values: pd.Series, base_status: pd.Series) -> pd.Series:
    """OPTIMAL when the metric is defined, propagate CUT_ERROR, else TIME_LIMIT.

    Lets the comparison-table helpers reuse their built-in TIME_LIMIT penalty
    for rows where the metric was never observed.
    """
    synth = pd.Series("TIME_LIMIT", index=values.index, dtype=object)
    synth[values.notna()] = "OPTIMAL"
    synth[base_status == "CUT_ERROR"] = "CUT_ERROR"
    return synth


def _prepare_metric_df(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    """Copy of ``df`` with the right ``status`` column for ``metric``.

    Log-derived metrics use a synthetic status (OPTIMAL where the metric is
    defined, TL otherwise); CSV metrics keep the original status.
    """
    sub = df.copy()
    if metric in _LOG_METRICS:
        sub["status"] = _synthetic_status(sub[metric], sub["status"])
    return sub


def _write_metric_subfolder(
    metric_dir: Path, metric: str, df_all: pd.DataFrame, baseline: str,
    label_map: dict[str, str] | None,
    include_settings: list[str] | None,
    df_all_instances: pd.DataFrame | None = None,
) -> None:
    """Write the four comparison tables for one metric.

    Both long tables (comparison.md and comparison_all.md) show ALL instances.
    - comparison.md:     † = OPTIMAL for all settings; stats over †-only.
    - comparison_all.md: † = ran for all settings without error/missing
                          (TL counts as ran); stats over †-only, TL penalized.
    For both: failed/missing instances (no †) are excluded from statistics
    for ALL settings, not just the one that failed. Instance names are noted.
    """
    metric_dir.mkdir(parents=True, exist_ok=True)
    has_baseline = baseline in df_all["settings"].unique()

    # Full dataset: all instances including errored; nullify metric for bad-status rows
    _df_for_all = df_all_instances if df_all_instances is not None else df_all
    _BAD_STATUS = {"CUT_ERROR", "INFEASIBLE", "INFEASIBLE_OR_UNBOUNDED", "DID_NOT_RUN"}
    _df_full = _df_for_all.copy()
    _df_full.loc[_df_full["status"].isin(_BAD_STATUS), metric] = np.nan

    # Error keys (--‡ annotation) and clean subset for short tables
    _, error_keys = cleanup_data(_df_for_all, baseline=baseline)
    df_clean, _ = cleanup_data(df_all, baseline=baseline)

    # comparison.md: all instances; † = OPTIMAL-for-all; stats over †-only
    long_clean = comparison_table(
        _df_full, target=metric, baseline=baseline,
        wilcoxon_test=has_baseline, label_map=label_map,
        include_settings=include_settings, error_keys=error_keys,
        mode="optimal",
    )
    # comparison_all.md: all instances; † = ran-for-all (no error/missing); TL penalized
    long_all = comparison_table(
        _df_full, target=metric, baseline=baseline,
        wilcoxon_test=has_baseline, label_map=label_map,
        include_settings=include_settings, error_keys=error_keys,
        mode="ran",
    )
    short = comparison_short_table(
        df_clean, target=metric, baseline=baseline,
        label_map=label_map, include_settings=include_settings,
        df_with_errors=df_all, wilcoxon_test=has_baseline,
    )
    # short_all: instances without CUT_ERROR/missing (same as †-semantics in comparison_all.md)
    # df_all is already filtered by the caller to exclude errored instances.
    short_all = comparison_short_all_table(
        df_all, target=metric, baseline=baseline,
        label_map=label_map, include_settings=include_settings,
        wilcoxon_test=has_baseline,
    )

    # Build footnotes with excluded-from-stats instance names
    _excl_clean = long_clean.attrs.get("excluded_from_stats", [])
    _excl_all = long_all.attrs.get("excluded_from_stats", [])

    _clean_md = to_markdown(long_clean)
    _clean_md += (
        "\n† solved by all settings within the time limit; included in summary statistics."
    )
    if _excl_clean:
        _clean_md += (
            f" {len(_excl_clean)} instance(s) excluded from statistics"
            f" (timed out or missing): {', '.join(_excl_clean)}."
        )
    _clean_md += "\n"
    (metric_dir / "comparison.md").write_text(_clean_md)
    (metric_dir / "comparison.tex").write_text(to_latex(long_clean))

    _all_md = to_markdown(long_all)
    _all_md += (
        "\n† ran for all settings without error; included in summary statistics"
        " (timed-out results penalized at the time limit)."
    )
    if _excl_all:
        _all_md += (
            f" {len(_excl_all)} instance(s) excluded from statistics"
            f" (CUT_ERROR or missing data): {', '.join(_excl_all)}."
        )
    _all_md += "\n"
    (metric_dir / "comparison_all.md").write_text(_all_md)
    (metric_dir / "comparison_all.tex").write_text(to_latex(long_all))
    (metric_dir / "comparison_short.md").write_text(
        comparison_short_to_markdown(short))
    (metric_dir / "comparison_short.tex").write_text(
        comparison_short_to_latex(short))
    (metric_dir / "comparison_short_all.md").write_text(
        comparison_short_to_markdown(short_all))
    (metric_dir / "comparison_short_all.tex").write_text(
        comparison_short_to_latex(short_all))
    print(f"Written {metric}/: comparison[_all|_short[_all]].{{md,tex}}")


def write_metric_tables(
    out_dir: Path,
    df: pd.DataFrame,
    baseline: str,
    label_map: dict[str, str] | None = None,
    include_settings: list[str] | None = None,
    df_all_instances: pd.DataFrame | None = None,
) -> None:
    """For each metric, write four tables under ``out_dir/<metric>/``."""
    for metric in METRICS:
        sub = _prepare_metric_df(df, metric)
        sub_all = (_prepare_metric_df(df_all_instances, metric)
                   if df_all_instances is not None else None)
        _write_metric_subfolder(out_dir / metric, metric, sub,
                                baseline=baseline, label_map=label_map,
                                include_settings=include_settings,
                                df_all_instances=sub_all)


_WS_IMPROVING_RTOL = 1e-4
_WS_IMPROVING_ATOL = 1e-9


def _warmstart_efficiency_rows(
    df: pd.DataFrame,
    settings_order: list[str],
    label_map: dict[str, str] | None,
) -> list[tuple[str, int, int, int]]:
    """Return ``(label, n_solved, n_ws_found, n_ws_improving)`` per setting."""
    ws_obj = pd.to_numeric(df["log_warmstart_objective"], errors="coerce")
    scopf  = pd.to_numeric(df["scopf_threshold"],         errors="coerce")
    tol    = np.maximum(_WS_IMPROVING_ATOL, _WS_IMPROVING_RTOL * scopf.abs())
    df2 = df.copy()
    df2["_ws_found"]     = ws_obj.notna()
    df2["_ws_improving"] = ws_obj < (scopf - tol)

    has_switch_col = "n_switching_actions" in df2.columns
    rows: list[tuple] = []
    for s in settings_order:
        sub = df2[df2["settings"] == s]
        if len(sub) == 0:
            continue
        label = label_map.get(s, s) if label_map else s
        n_inst      = int(len(sub))
        n_solved    = int((sub["status"] == "OPTIMAL").sum())
        n_found     = int(sub["_ws_found"].sum())
        n_improving = int(sub["_ws_improving"].sum())
        if has_switch_col:
            opt_sub = sub[sub["status"] == "OPTIMAL"]
            n_actions = pd.to_numeric(opt_sub["n_switching_actions"],
                                      errors="coerce").fillna(0)
            n_switch = int((n_actions > 0).sum())
            rows.append((label, n_inst, n_solved, n_found, n_improving, n_switch))
        else:
            rows.append((label, n_inst, n_solved, n_found, n_improving))
    return rows


def _rows_to_markdown(
    header: list[str], rows: list[tuple], alignments: list[str],
) -> str:
    """Render ``header`` + ``rows`` as a GFM Markdown table."""
    widths = [max(len(h), max((len(str(r[i])) for r in rows), default=0))
              for i, h in enumerate(header)]
    sep_parts = []
    for w, a in zip(widths, alignments):
        if a == "r":
            sep_parts.append("-" * (w - 1) + ":")
        elif a == "c":
            sep_parts.append(":" + "-" * (w - 2) + ":")
        else:
            sep_parts.append("-" * w)
    h_line  = "| " + " | ".join(h.ljust(w) for h, w in zip(header, widths)) + " |"
    s_line  = "| " + " | ".join(sep_parts) + " |"
    r_lines = ["| " + " | ".join(str(v).rjust(w) if a == "r" else str(v).ljust(w)
                                  for v, w, a in zip(row, widths, alignments))
               + " |" for row in rows]
    return "\n".join([h_line, s_line] + r_lines) + "\n"


def _rows_to_latex(header: list[str], rows: list[tuple]) -> str:
    """Render ``header`` + ``rows`` as a LaTeX booktabs tabular."""
    ncols = len(header)
    col_fmt = "l" + "r" * (ncols - 1)
    lines = [
        r"\begin{tabular}{" + col_fmt + r"}",
        r"\toprule",
        " & ".join(r"\textbf{" + h + "}" for h in header) + r" \\",
        r"\midrule",
    ]
    for row in rows:
        lines.append(" & ".join(str(v) for v in row) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines) + "\n"


def write_warmstart_efficiency_table(
    out_dir: Path,
    df: pd.DataFrame,
    label_map: dict[str, str] | None = None,
    include_settings: list[str] | None = None,
) -> None:
    """Write a per-setting warmstart efficiency summary.

    For each setting counts the number of instances where:

    * **Instances** — total instances in the experiment.
    * **Solved** — status is ``OPTIMAL``.
    * **Warmstart loaded** — main MP received a MIP start
      (``log_warmstart_objective`` is not NaN).
    * **RSCOTS improved** — the loaded warmstart is strictly better than the
      SCOPF threshold (``log_warmstart_objective < scopf_threshold - tol``).
    * **Switch in opt.** — among OPTIMAL solves, ``n_switching_actions > 0``
      (present only when the column is available in ``df``).

    Output is written to ``out_dir/warmstart_efficiency/``.
    """
    tbl_dir = out_dir / "warmstart_efficiency"
    tbl_dir.mkdir(parents=True, exist_ok=True)

    settings_order = include_settings or sorted(df["settings"].unique())
    rows = _warmstart_efficiency_rows(df, settings_order, label_map)
    has_switch_col = len(rows) > 0 and len(rows[0]) == 6

    if has_switch_col:
        header = ["Setting", "Instances", "Solved",
                  "Warmstart loaded", "RSCOTS improved", "Switch in opt."]
        aligns = ["l", "r", "r", "r", "r", "r"]
    else:
        header = ["Setting", "Instances", "Solved",
                  "Warmstart loaded", "RSCOTS improved"]
        aligns = ["l", "r", "r", "r", "r"]

    (tbl_dir / "warmstart_efficiency.md").write_text(
        _rows_to_markdown(header, rows, aligns))
    (tbl_dir / "warmstart_efficiency.tex").write_text(
        _rows_to_latex(header, rows))
    print("Written warmstart_efficiency/: warmstart_efficiency.{md,tex}")


def run(
    experiment_dir: Path,
    baseline: str = "no_warmstart",
    label_map: dict[str, str] | None = None,
    include_settings: list[str] | None = None,
    use_csvs_if_possible: bool = True,
    output_dir: Path | None = None,
) -> tuple[pd.DataFrame, dict, dict]:
    """Run the full log-metrics pipeline and return ``(df, raw_cache, logs)``.

    The returned ``raw_cache`` and ``logs`` can be passed to
    :func:`~tables.run_full_comparison` (``log_raw_cache`` /
    ``log_files_map`` params) so that Gurobi logs are not parsed a second
    time for the ``gap_map``.  Pass ``use_csvs_if_possible=False`` to force
    every ``.out`` log to be re-parsed (also overwriting the on-disk
    ``logs_parsed/<stem>.csv`` cache).
    """
    # ``output_dir`` (zero-touch) redirects derived artifacts out of the
    # read-only source repo; log parsing still reads the source logs.
    out_dir = Path(output_dir) if output_dir is not None else experiment_dir / "log_metrics"
    out_dir.mkdir(parents=True, exist_ok=True)
    logs, raw_cache = build_log_cache(experiment_dir, use_csvs_if_possible=use_csvs_if_possible)
    df = collect_log_metrics(experiment_dir, raw_cache=raw_cache, logs=logs)
    write_metrics_csv(df, out_dir / "log_metrics.csv")
    df_tables, dropped = drop_incomplete_instances(df, include_settings=include_settings)
    if dropped:
        print(f"log_metrics: excluded {len(dropped)} incomplete instance(s): "
              f"{[d[0].replace('pglib_scopf_case', '') if isinstance(d, tuple) else str(d) for d in dropped]}")
    write_metric_tables(out_dir, df_tables, baseline=baseline,
                        label_map=label_map, include_settings=include_settings,
                        df_all_instances=df)
    write_warmstart_efficiency_table(out_dir, df_tables, label_map=label_map,
                                     include_settings=include_settings)
    return df_tables, raw_cache, logs


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python log_metrics.py <experiment_dir> [baseline]",
              file=sys.stderr)
        sys.exit(1)
    exp = Path(sys.argv[1]).resolve()
    base = sys.argv[2] if len(sys.argv) > 2 else "no_warmstart"
    if not exp.is_dir():
        print(f"Not a directory: {exp}", file=sys.stderr)
        sys.exit(1)
    out_df, _, _ = run(exp, baseline=base)
    n = len(out_df)
    counts = {m: int(out_df[m].notna().sum()) for m in METRICS}
    print(f"Processed {n} rows: " +
          ", ".join(f"{c} {m}" for m, c in counts.items()) + ".")
