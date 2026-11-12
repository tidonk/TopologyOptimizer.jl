"""Gurobi log parser and primal-dual integral plotter.

Parses Gurobi MIP log files, extracts primal/dual convergence data
(including warmstart restarts), computes the primal-dual integral,
and plots publication-ready convergence figures.

Supports overlaying multiple log files on a single plot.
"""

from pathlib import Path
import json
import re
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from style import COLORS, setup_latex_style, fmt_setting


# ---------------------------------------------------------------------------
# Parsed-log on-disk format (CSV + JSON header) and read_log() dispatcher.
#
# A parsed log is dumped as
#
#     # META {"objective": ..., "restarts": [...], ...}
#     t,primal,dual,phase
#     0.0,NaN,123.4,ws0
#     ...
#
# stored at ``<setting_dir>/logs_parsed/<log_stem>.csv`` (sibling of
# ``<setting_dir>/logs/<log_stem>.out``).  :func:`read_log` is the unified
# entry point: it loads from CSV when one exists and *use_csvs_if_possible*
# is True, otherwise it parses the raw log and (always) overwrites the CSV.
# ---------------------------------------------------------------------------

_META_PREFIX = "# META "


def _parsed_csv_path(log_path: Path) -> Path:
    """Return the parsed-CSV path corresponding to *log_path*."""
    log_path = Path(log_path)
    return log_path.parent.parent / "logs_parsed" / f"{log_path.stem}.csv"


def _phase_labels(times: np.ndarray, restarts: list[float]) -> list[str]:
    """Label each trajectory point by phase: ``ws0, ws1, ..., main``."""
    if len(times) == 0:
        return []
    n_ws = len(restarts)
    labels = []
    for t in times:
        # Count how many restarts have occurred up to and including t.
        k = int(np.searchsorted(restarts, t, side="right")) if n_ws else 0
        labels.append("main" if k >= n_ws else f"ws{k}")
    return labels


def dump_parsed_log(data: dict, csv_path: Path) -> Path:
    """Write *data* (a :func:`parse_gurobi_log` result) to *csv_path*.

    Trajectory columns ``t,primal,dual,phase`` go into the CSV body; all
    scalar fields (``optimal``, ``objective``, ``warmstart_objective``,
    ``gap``, ``total_time``, ``total_nodes``, ``total_iterations``,
    ``restarts``) are JSON-serialised on a leading ``# META …`` comment line.
    """
    csv_path = Path(csv_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    meta = {
        "optimal":             bool(data["optimal"]),
        "objective":           float(data["objective"]),
        "warmstart_objective": float(data["warmstart_objective"]),
        "scopf_objective":     float(data.get("scopf_objective", np.nan)),
        "gap":                 float(data["gap"]),
        "total_time":          float(data["total_time"]),
        "total_nodes":         int(data["total_nodes"]),
        "total_iterations":    int(data["total_iterations"]),
        "restarts":            [float(r) for r in data["restarts"]],
    }
    phases = _phase_labels(data["times"], meta["restarts"])
    df = pd.DataFrame({
        "t":      data["times"],
        "primal": data["primals"],
        "dual":   data["duals"],
        "phase":  phases,
    })
    with csv_path.open("w") as fh:
        fh.write(_META_PREFIX + json.dumps(meta, allow_nan=True) + "\n")
        df.to_csv(fh, index=False)
    return csv_path


def load_parsed_log(csv_path: Path) -> dict:
    """Reverse of :func:`dump_parsed_log` — reconstruct a parse_gurobi_log dict."""
    csv_path = Path(csv_path)
    with csv_path.open() as fh:
        first = fh.readline()
        if not first.startswith(_META_PREFIX):
            raise ValueError(f"{csv_path}: missing META header")
        meta = json.loads(first[len(_META_PREFIX):])
        if "scopf_objective" not in meta:
            raise ValueError(f"{csv_path}: missing scopf_objective (stale cache — will re-parse)")
        df = pd.read_csv(fh)
    return {
        "times":               df["t"].to_numpy(dtype=float),
        "primals":             df["primal"].to_numpy(dtype=float),
        "duals":               df["dual"].to_numpy(dtype=float),
        "optimal":             bool(meta["optimal"]),
        "objective":           float(meta["objective"]),
        "warmstart_objective": float(meta["warmstart_objective"]),
        "scopf_objective":     float(meta["scopf_objective"]),
        "gap":                 float(meta["gap"]),
        "total_time":          float(meta["total_time"]),
        "total_nodes":         int(meta["total_nodes"]),
        "total_iterations":    int(meta["total_iterations"]),
        "restarts":            list(meta["restarts"]),
    }


def read_log(log_path: str | Path, use_csvs_if_possible: bool = True) -> dict:
    """Unified entry: load parsed CSV when available, else parse raw log.

    When *use_csvs_if_possible* is ``True`` and a sibling
    ``logs_parsed/<stem>.csv`` exists, that CSV is loaded (no log parsing).
    Otherwise the raw ``.out`` log is parsed via :func:`parse_gurobi_log`
    and the resulting trajectory is exported to the CSV path (overwriting
    any stale CSV).  This makes ``use_csvs_if_possible=False`` a full
    force-reparse that refreshes the on-disk cache as a side effect.
    """
    log_path = Path(log_path)
    csv_path = _parsed_csv_path(log_path)
    if use_csvs_if_possible and csv_path.is_file():
        try:
            return load_parsed_log(csv_path)
        except Exception as exc:  # noqa: BLE001 — fall back to a fresh parse
            warnings.warn(f"failed to load {csv_path} ({exc}); re-parsing")
    data = parse_gurobi_log(log_path, main_run_only=False)
    try:
        dump_parsed_log(data, csv_path)
    except OSError as exc:
        warnings.warn(f"failed to dump {csv_path} ({exc})")
    return data


def parse_gurobi_log(filepath: str | Path, main_run_only: bool = False) -> dict:
    """Parse a Gurobi log file and extract convergence data.

    Parameters
    ----------
    filepath : Gurobi .out log path.
    main_run_only : if True, drop any data from preceding solver invocations
        (e.g. warmstart sub-MIPs) and report only the final ``Gurobi Optimizer
        version`` block — times are then relative to the start of that block.

    Returns
    -------
    dict with keys:
        times, primals, duals : list[float] — convergence trajectory
        optimal : bool — whether optimal solution was found
        objective : float — best objective (NaN if not found)
        total_time, total_nodes, total_iterations : float/int
        restarts : list[float] — time offsets where solver restarted
    """
    filepath = Path(filepath)
    text = filepath.read_text()
    lines = text.split("\n")

    times, primals, duals = [], [], []
    optimal = False
    objective = np.nan
    warmstart_obj = np.nan
    gap = np.nan
    total_time = 0.0
    total_nodes = 0
    total_iterations = 0
    restarts = []

    time_offset = 0.0
    run_active = False
    scopf_obj = np.nan

    def _parse_val(s: str) -> float:
        if s == "-":
            return np.nan
        try:
            return float(s)
        except (ValueError, TypeError):
            return np.nan

    for line in lines:
        if "Gurobi Optimizer version" in line:
            if main_run_only:
                times, primals, duals = [], [], []
                restarts = []
                time_offset = 0.0
                total_time = 0.0
                total_nodes = 0
                total_iterations = 0
                objective = np.nan
                warmstart_obj = np.nan
                scopf_obj = np.nan
                gap = np.nan
                optimal = False
            elif times:
                time_offset = max(times)
                restarts.append(time_offset)
            run_active = False

        elif "Nodes    |    Current Node    |     Objective Bounds      |     Work" in line:
            run_active = True

        elif run_active and re.search(r"\d+s", line):
            tokens = line.split()
            time_match = re.search(r"(\d+)s", line)
            if not time_match:
                continue
            t = float(time_match.group(1)) + time_offset

            # Find the token ending with 's' to locate bounds
            idx = next(
                (i for i, tok in enumerate(tokens) if tok.endswith("s") and tok[:-1].isdigit()),
                None,
            )
            if idx is not None and idx >= 4:
                primal = _parse_val(tokens[idx - 4])
                dual = _parse_val(tokens[idx - 3])
                times.append(t)
                primals.append(primal)
                duals.append(max(dual, 0.0) if not np.isnan(dual) else dual)

        elif "Optimal solution found" in line:
            optimal = True

        elif "Loaded user MIP start with objective" in line:
            m = re.search(r"Loaded user MIP start with objective ([+-]?\d[\d.eE+-]*)", line)
            if m:
                v = float(m.group(1))
                warmstart_obj = v
                if np.isnan(scopf_obj):  # first occurrence = SCOPF result loaded into RSCOTS
                    scopf_obj = v

        elif "Best objective" in line:
            m = re.search(r"Best objective ([+-]?\d[\d.eE+-]*)", line)
            if m:
                objective = float(m.group(1))
            m_gap = re.search(r"gap ([\d.]+)%", line)
            if m_gap:
                gap = float(m_gap.group(1)) / 100.0

        elif "Explored" in line:
            m = re.search(
                r"Explored (\d+) nodes \((\d+) simplex iterations\) in ([\d.]+) seconds", line
            )
            if m:
                total_nodes += int(m.group(1))
                total_iterations += int(m.group(2))
                total_time += float(m.group(3))

    return {
        "times": np.array(times, dtype=float),
        "primals": np.array(primals, dtype=float),
        "duals": np.array(duals, dtype=float),
        "optimal": optimal,
        "objective": objective,
        "warmstart_objective": warmstart_obj,
        "scopf_objective": scopf_obj,
        "gap": gap,
        "total_time": total_time,
        "total_nodes": total_nodes,
        "total_iterations": total_iterations,
        "restarts": restarts,
    }


def compute_pd_integral(times: np.ndarray, primals: np.ndarray, duals: np.ndarray,
                        total_time: float | None = None) -> float:
    """Compute the primal-dual integral (area between primal and dual curves).

    Gurobi's primal/dual bounds are **piecewise-constant step functions** — a
    bound holds from one node-table emission until the next — so the integral
    is computed as ``Σ gap_i · (t_{i+1} - t_i)`` (left-continuous rectangles),
    not via trapezoidal averaging.

    For the leading interval before any incumbent has been found (NaN primals),
    the first eventually-discovered incumbent is back-filled, so the gap is
    well defined from t=0. When ``total_time`` is given the trailing interval
    from the last emission to ``total_time`` is also accumulated with the final
    gap held constant.
    """
    if len(times) == 0:
        return 0.0
    p = primals.astype(float, copy=True)
    finite_p = np.where(~np.isnan(p))[0]
    if len(finite_p) > 0:
        p[:finite_p[0]] = p[finite_p[0]]
    mask = ~(np.isnan(p) | np.isnan(duals))
    t, p, d = times[mask], p[mask], duals[mask]
    if len(t) == 0:
        return 0.0
    gap = np.abs(p - d)
    # left-continuous step: gap[i] holds over [t[i], t[i+1])
    dt   = np.diff(t)
    area = float(np.sum(gap[:-1] * dt))
    if total_time is not None and total_time > t[-1]:
        area += float(gap[-1] * (total_time - t[-1]))
    return area


def compute_primal_integral(
    times: np.ndarray,
    primals: np.ndarray,
    optimal_obj: float,
    total_time: float | None = None,
) -> float:
    """Compute the primal integral: area between the primal bound and the known optimum.

    Uses the same left-continuous step-function convention as
    :func:`compute_pd_integral`.  Leading NaN primal intervals are back-filled
    from the first feasible incumbent; if no incumbent is ever found the
    function returns ``float("nan")``.  The trailing interval from the last log
    entry to ``total_time`` (if given) is accumulated with the final primal
    value held constant.

    For a minimization problem this equals
    ``Σ max(P_i − f*, 0) · (t_{i+1} − t_i)``,
    measuring total solver-time spent above the optimal objective.
    Returns ``0.0`` for an empty trajectory (root/presolve solve with no
    node-table entries) and ``float("nan")`` when *optimal_obj* is unknown.
    """
    if np.isnan(optimal_obj):
        return float("nan")
    if len(times) == 0:
        return 0.0
    p = primals.astype(float, copy=True)
    finite_p = np.where(~np.isnan(p))[0]
    if len(finite_p) == 0:
        return float("nan")
    p[:finite_p[0]] = p[finite_p[0]]
    mask = ~np.isnan(p)
    t, p = times[mask], p[mask]
    if len(t) == 0:
        return float("nan")
    gap = np.maximum(p - optimal_obj, 0.0)
    dt = np.diff(t)
    area = float(np.sum(gap[:-1] * dt))
    if total_time is not None and total_time > t[-1]:
        area += float(gap[-1] * (total_time - t[-1]))
    return area


def compute_time_to_first(times: np.ndarray, primals: np.ndarray) -> float:
    """Time of the first feasible incumbent in the trajectory.

    Returns the earliest entry where the primal column is not NaN (Gurobi prints
    ``-`` until an incumbent is found). Returns NaN if no incumbent was logged.
    """
    if len(times) == 0:
        return float("nan")
    mask = ~np.isnan(primals)
    if not mask.any():
        return float("nan")
    return float(times[mask][0])


def compute_time_to_first_improving(
    times: np.ndarray,
    primals: np.ndarray,
    threshold_obj: float | None = None,
    rtol: float = 1e-4,
    atol: float = 1e-9,
) -> float:
    """Time of the first incumbent strictly better than ``threshold_obj``.

    For minimization the new incumbent must satisfy
    ``primal < threshold_obj - max(atol, rtol * |threshold_obj|)``. The relative
    tolerance absorbs Gurobi's truncated printout (e.g. ``Loaded user MIP start
    with objective 68674.2`` or the first sub-MIP incumbent).

    When ``threshold_obj`` is None or NaN the first incumbent in the trajectory
    is used as the reference (so the metric measures the time until the
    **initial** feasible solution is strictly improved upon — across the whole
    warmstart + main-run timeline).

    The full timeline (including warmstart sub-MIP stages) is searched: if a
    warmstart stage (SCOPF, RSCOTS) produces an incumbent strictly better than
    ``threshold_obj``, the timer is stopped there.
    """
    if len(times) == 0:
        return float("nan")
    mask = ~np.isnan(primals)
    if not mask.any():
        return float("nan")
    t, p = np.asarray(times, dtype=float)[mask], np.asarray(primals, dtype=float)[mask]
    if threshold_obj is None or np.isnan(threshold_obj):
        threshold_obj = float(p[0])
    threshold = threshold_obj - max(atol, rtol * abs(threshold_obj))
    hits = np.where(p < threshold)[0]
    if len(hits) == 0:
        return float("nan")
    return float(t[hits[0]])


def compute_time_to_optimal(
    times: np.ndarray,
    primals: np.ndarray,
    objective: float,
    rtol: float = 1e-6,
    atol: float = 1e-9,
) -> float:
    """Time at which the incumbent first reached the final objective value.

    The Gurobi log truncates printed digits (e.g. ``8026.7547538`` vs ``8026.755``
    or ``Best objective 8.026754e+03``), so we compare with a relative tolerance.
    Returns NaN if ``objective`` is NaN or no incumbent matched it.
    """
    if np.isnan(objective) or len(times) == 0:
        return float("nan")
    mask = ~np.isnan(primals)
    if not mask.any():
        return float("nan")
    t = times[mask]
    p = primals[mask]
    tol = max(atol, rtol * abs(objective))
    hits = np.where(np.abs(p - objective) <= tol)[0]
    if len(hits) == 0:
        return float("nan")
    return float(t[hits[0]])


def plot_primal_dual(
    log_files: list[str | Path],
    labels: list[str] | None = None,
    output: str | Path | None = None,
) -> plt.Figure:
    """Plot primal-dual convergence for one or more Gurobi log files.

    Parameters
    ----------
    log_files : paths to Gurobi .out log files.
    labels : display labels for each log file (defaults to filenames).
    output : save path (extension determines format).
    """
    setup_latex_style()

    if labels is None:
        labels = [Path(f).stem for f in log_files]

    fig, ax = plt.subplots(figsize=(8, 5))

    all_restarts = []
    for i, (logfile, label) in enumerate(zip(log_files, labels)):
        data = read_log(logfile)
        t, p, d = data["times"], data["primals"], data["duals"]

        if len(t) == 0:
            continue

        integral = compute_pd_integral(t, p, d)
        color = COLORS[i % len(COLORS)]

        primal_label = f"{label} primal (PDI$={integral:.1f}$)"
        dual_label = f"{label} dual"

        ax.plot(t, p, color=color, linewidth=2, marker="x", markersize=3, label=primal_label)
        ax.plot(t, d, color=color, linewidth=2, linestyle="--", label=dual_label)

        # Optimal line
        if data["optimal"] and not np.isnan(data["objective"]):
            ax.axhline(data["objective"], color=color, linestyle=":", linewidth=1, alpha=0.6)

        # Collect restarts with time offset for multi-log
        for r in data["restarts"]:
            all_restarts.append(r)

    # Draw restart markers
    for j, r in enumerate(all_restarts):
        ax.axvline(r, color="gray", linestyle=":", linewidth=0.8, alpha=0.5)
        ax.annotate(
            f"Restart {j + 1}", xy=(r, 1), xycoords=("data", "axes fraction"),
            fontsize=8, color="gray", ha="right", va="top", rotation=90,
        )

    ax.set_xscale("log")
    ax.set_xlabel(r"Time (s)")
    ax.set_ylabel(r"Objective value")
    ax.legend(frameon=False, fontsize=9)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:,.0f}"))

    fig.tight_layout()
    if output:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output)
        print(f"Saved primal-dual plot to {output}")
    return fig


def _interpolate_to_grid(times, values, grid):
    """Interpolate values onto a common time grid, forward-filling NaNs."""
    result = np.full_like(grid, np.nan)
    mask = ~np.isnan(values)
    t, v = times[mask], values[mask]
    if len(t) == 0:
        return result
    # Forward-fill: for each grid point, use the last known value
    idxs = np.searchsorted(t, grid, side="right") - 1
    valid = idxs >= 0
    result[valid] = v[idxs[valid]]
    return result


def plot_aggregated_pd(
    log_groups: dict[str, list[str | Path]],
    output: str | Path | None = None,
    label_map: dict[str, str] | None = None,
    style_map: dict[str, dict] | None = None,
    show_pdi_table: bool = True,
    relative: bool = True,
    include_settings: list[str] | None = None,
    top_n: int | None = None,
    baseline: str | None = None,
    include_warmstart: bool = True,
    scopf_thresholds: dict[str, float] | None = None,
    true_objectives: dict[str, float] | None = None,
    true_obj_rel_threshold: float = 0.05,
) -> plt.Figure:
    """Plot aggregated (mean) primal-dual convergence per setting.

    Parameters
    ----------
    log_groups : mapping from setting name → list of Gurobi log file paths.
    output : save path.
    label_map : optional mapping from setting name → display label.
    style_map : optional mapping from setting name → dict with 'color', 'linestyle' keys.
    show_pdi_table : whether to display the Mean PDI table next to the plot.
    relative : if True, normalize per instance before aggregating:
        - time axis → t / T_total (fraction of total solve time)
        - obj axis  → value / best_objective (relative to optimal)
    top_n : if set, only plot the *top_n* settings with lowest PDI plus the
        baseline (if provided and present).  Selection is based on a fast
        pre-computation pass over all log files.
    baseline : name of the baseline setting — always included when *top_n* is set.
    include_warmstart : if True (default), trajectories include the warmstart
        sub-MIP timeline.  If False, each trajectory is trimmed to start at
        ``restarts[-1]`` (start of the main solve) and its time axis is
        shifted to start at ``t=0``.
    scopf_thresholds : optional ``{instance_name: scopf_obj}`` map.  When
        given, a trajectory whose leading interval has no incumbent yet is
        anchored at the per-instance SCOPF objective (rather than at ``f*``)
        — both ``no_warmstart`` and warmstarted runs then enter the plot at
        the same SCOPF reference quality, so the aggregated curve descends
        monotonically (on average) instead of jumping up once the first
        no-warmstart incumbent arrives.  Without this map, leading
        no-incumbent intervals are anchored at ``f*`` (relative value
        ``1.0``).
    true_objectives : optional ``{instance_name: scopf_obj}`` map of the
        true compiled SCOPF optimal objectives (e.g. from compiled.csv).
        When provided it is used as f* for normalization (overriding the
        Gurobi log's ``best_obj``) and any log whose best objective exceeds
        the compiled optimum by more than *true_obj_rel_threshold* is
        excluded.  This filters out runs where the Benders auxiliary variable
        η did not converge to 0, which would otherwise inflate f* and
        produce misleadingly flat (near-optimal) trajectories.
    true_obj_rel_threshold : relative tolerance for the log-vs-compiled
        objective consistency check.  Logs with
        ``log_obj / true_obj - 1 > true_obj_rel_threshold`` are excluded.
        Default 0.05 (5 %).
    """
    setup_latex_style()

    def _anchor_value(d: dict, log_path: Path, setting: str) -> float:
        """Initial primal to use before any real incumbent (NaN-safe).

        Returns ``max(SCOPF_threshold, first_real_primal)`` so the leading
        step never *jumps up* once the first true incumbent arrives.  For
        warmstarted settings the first main-solve primal is usually no
        worse than the SCOPF result, so the anchor reduces to the SCOPF
        threshold (preserving the "warmstarted runs enter at SCOPF
        quality" semantics).  For ``no_warmstart`` the first incumbent is
        typically worse than SCOPF, so the anchor back-fills with that
        first real primal — keeping the per-instance step function
        monotonically non-increasing.
        """
        scopf_anchor = np.nan
        if scopf_thresholds is not None:
            inst = _extract_instance_name(log_path, setting)
            v = scopf_thresholds.get(inst)
            if v is not None and not np.isnan(v):
                scopf_anchor = float(v)
        if np.isnan(scopf_anchor):
            scopf_anchor = d["objective"]
        primals_nz = d["primals"][~np.isnan(d["primals"])]
        first_real = float(primals_nz[0]) if len(primals_nz) else np.nan
        if np.isnan(first_real):
            return scopf_anchor
        if np.isnan(scopf_anchor):
            return first_real
        return max(scopf_anchor, first_real)

    def _prep_data(setting: str, files):
        """Parse log files, optionally trim warmstart, anchor leading NaNs."""
        all_data = []
        for f in files:
            f = Path(f)
            d = read_log(f)
            # Skip logs where the Gurobi MP objective is inflated above the
            # compiled SCOPF optimal (indicates Benders η did not converge to 0).
            if true_objectives is not None:
                _inst = _extract_instance_name(f, setting)
                _true_opt = true_objectives.get(_inst)
                if (
                    _true_opt is not None
                    and not np.isnan(_true_opt)
                    and _true_opt > 0
                    and not np.isnan(d["objective"])
                    and d["objective"] / _true_opt - 1 > true_obj_rel_threshold
                ):
                    continue
            if not include_warmstart and d["restarts"]:
                t0 = d["restarts"][-1]
                m = d["times"] >= t0
                d = {**d,
                     "times":      d["times"][m] - t0,
                     "primals":    d["primals"][m],
                     "duals":      d["duals"][m],
                     "total_time": max(d["total_time"] - t0, 0.0),
                     "restarts":   []}
            if len(d["times"]) == 0:
                continue
            # Anchor the leading no-incumbent interval at the SCOPF threshold
            # (or, lacking one, at the optimum) so the aggregated curve has a
            # well-defined value from t=0 instead of leaving a gap.
            if np.isnan(d["primals"][0]):
                anchor = _anchor_value(d, f, setting)
                if not np.isnan(anchor):
                    d["times"]   = np.concatenate([[0.0],           d["times"]])
                    d["primals"] = np.concatenate([[anchor],        d["primals"]])
                    d["duals"]   = np.concatenate([[d["duals"][0]], d["duals"]])
            # Extend to terminal time with the final bound held constant.
            t_end = d["total_time"]
            if t_end > 0 and t_end > d["times"][-1]:
                obj = d["objective"]
                final_p = obj if (d["optimal"] and not np.isnan(obj)) else d["primals"][-1]
                final_d = obj if (d["optimal"] and not np.isnan(obj)) else d["duals"][-1]
                d["times"]   = np.append(d["times"],   t_end)
                d["primals"] = np.append(d["primals"], final_p)
                d["duals"]   = np.append(d["duals"],   final_d)
            d["_instance_name"] = _extract_instance_name(f, setting)
            all_data.append(d)
        return all_data

    # ref_times is populated after parsed_cache is built (see below)
    ref_times: dict[str, float] = {}

    def _mean_pdi(all_data):
        integrals = []
        for d in all_data:
            if relative:
                inst = d.get("_instance_name", "")
                T_ref = ref_times.get(inst, 0.0)
                t_total = T_ref if T_ref > 0 else (d["total_time"] if d["total_time"] > 0 else d["times"][-1])
                opt = (true_objectives.get(inst) if true_objectives and inst else None)
                if opt is None or np.isnan(opt):
                    opt = d["objective"]
                if t_total <= 0 or np.isnan(opt) or opt == 0:
                    continue
                rel_t = d["times"] / t_total
                integrals.append(compute_pd_integral(
                    rel_t, d["primals"] / opt, d["duals"] / opt))
            else:
                integrals.append(compute_pd_integral(d["times"], d["primals"], d["duals"]))
        return np.mean(integrals) if integrals else float("inf")

    def _line_style(i: int, s: str) -> dict:
        if style_map and s in style_map:
            return {k: v for k, v in style_map[s].items() if k in ("color", "linestyle")}
        return {"color": COLORS[i % len(COLORS)]}

    # Sort log_groups by include_settings order
    ordered_settings = sorted(
        log_groups.keys(),
        key=lambda s: include_settings.index(s)
        if include_settings and s in include_settings else 999,
    )

    # Pre-compute PDI for all settings and cache parsed data
    parsed_cache: dict[str, list] = {}
    pdi_cache: dict[str, float] = {}
    for s in ordered_settings:
        data = _prep_data(s, log_groups[s])
        if data:
            parsed_cache[s] = data
    # Build per-instance reference time from baseline for fair normalised comparison.
    # All settings are then normalised by T_ref[instance] = baseline total_time.
    if baseline and relative and baseline in parsed_cache:
        for d in parsed_cache[baseline]:
            t = d["total_time"] if d["total_time"] > 0 else (d["times"][-1] if len(d["times"]) > 0 else 0.0)
            if t > 0:
                ref_times[d["_instance_name"]] = t
    for s, data in parsed_cache.items():
        pdi_cache[s] = _mean_pdi(data)

    # Select settings to plot: top_n by lowest PDI + baseline
    if top_n is not None:
        non_base = [(s, v) for s, v in pdi_cache.items() if s != baseline]
        non_base_sorted = sorted(non_base, key=lambda x: x[1])
        selected = set(s for s, _ in non_base_sorted[:top_n])
        if baseline and baseline in pdi_cache:
            selected.add(baseline)
        plot_order = [s for s in ordered_settings if s in selected]
    else:
        plot_order = [s for s in ordered_settings if s in parsed_cache]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    n_total = 0
    pdi_lines = []

    for i, setting in enumerate(plot_order):
        all_data = parsed_cache.get(setting)
        if not all_data:
            continue

        n_total += len(all_data)

        if relative:
            grid = np.linspace(0, 1, 500)
            primal_stack, dual_stack = [], []
            for d in all_data:
                inst = d.get("_instance_name", "")
                T_ref = ref_times.get(inst, 0.0)
                t_total = T_ref if T_ref > 0 else (d["total_time"] if d["total_time"] > 0 else d["times"][-1])
                opt = (true_objectives.get(inst) if true_objectives and inst else None)
                if opt is None or np.isnan(opt):
                    opt = d["objective"]
                if t_total <= 0 or np.isnan(opt) or opt == 0:
                    continue
                rel_t = d["times"] / t_total
                primal_stack.append(_interpolate_to_grid(rel_t, d["primals"] / opt, grid))
                dual_stack.append(_interpolate_to_grid(rel_t, d["duals"] / opt, grid))
            if not primal_stack:
                continue
            primal_stack = np.array(primal_stack)
            dual_stack = np.array(dual_stack)
        else:
            t_max = max(d["times"].max() for d in all_data)
            grid = np.linspace(1, t_max, 500)
            primal_stack = np.array([
                _interpolate_to_grid(d["times"], d["primals"], grid) for d in all_data])
            dual_stack = np.array([
                _interpolate_to_grid(d["times"], d["duals"], grid) for d in all_data])

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            mean_primal = np.nanmean(primal_stack, axis=0)
            mean_dual = np.nanmean(dual_stack, axis=0)
        if np.all(np.isnan(mean_primal)):
            continue

        mean_integral = pdi_cache.get(setting, 0.0)
        ls = _line_style(i, setting)
        color = ls.get("color", "black")
        lbl = fmt_setting(setting, label_map)

        # Solid line for mean primal; shaded band down to mean dual (the PD gap)
        ax.plot(grid, mean_primal, linewidth=1.8, label=lbl, **ls)
        ax.fill_between(grid, mean_dual, mean_primal,
                        color=color, alpha=0.12, linewidth=0)

        pdi_lines.append((lbl, mean_integral))

    if relative:
        ax.set_xlabel(r"Relative time $t / T_{\mathrm{ref}}$" if ref_times else r"Relative time $t / T$")
        ax.set_ylabel(r"Relative objective $f / f^*$")
        ax.set_ylim(bottom=0.5, top=min(10.0, ax.get_ylim()[1]))
    else:
        ax.set_xscale("log")
        ax.set_xlabel(r"Time (s)")
        ax.set_ylabel(r"Objective value")
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:,.0f}"))

    ax.legend(frameon=False, fontsize="small", title="Setting",
              loc="upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0)

    # Mean PDI table (sorted by PDI ascending)
    if show_pdi_table and pdi_lines:
        pdi_lines_sorted = sorted(pdi_lines, key=lambda x: x[1])
        col_labels = ["Setting", r"$\overline{\mathrm{PDI}}$"]
        cell_text = [[name, f"${val:.2f}$"] for name, val in pdi_lines_sorted]
        n_rows = len(pdi_lines_sorted) + 1
        row_h = 0.06
        table_h = row_h * n_rows
        table = ax.table(cellText=cell_text, colLabels=col_labels,
                         loc="lower left", bbox=[1.02, 0.0, 0.50, table_h],
                         colWidths=[0.72, 0.28])
        table.auto_set_font_size(False)
        table.set_fontsize(7)
        for key, cell in table.get_celld().items():
            cell.set_edgecolor("gray")
            cell.set_linewidth(0.5)
            cell.PAD = 0.03
            if key[0] == 0:
                cell.set_facecolor("#f0f0f0")
                cell.set_text_props(weight="bold")
            else:
                cell.set_facecolor("white")

    fig.tight_layout()
    fig.subplots_adjust(right=0.65)

    if output:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, bbox_inches="tight")
        if output.suffix == ".pdf":
            fig.savefig(output.with_suffix(".png"), bbox_inches="tight")
        print(f"Saved aggregated PD plot to {output}")
    return fig


def _extract_instance_name(filepath: Path, setting: str) -> str:
    """Extract instance name from a Gurobi log filename.

    Filename pattern: GRB_MP_<instance>_<setting>_<date>.out
    """
    stem = filepath.stem
    name = stem.removeprefix("GRB_MP_")
    # Remove _<setting>_YYYY-MM-DD_HH_MM_SS suffix
    suffix = re.sub(
        rf"_{re.escape(setting)}_\d{{4}}-\d{{2}}-\d{{2}}_\d{{2}}_\d{{2}}_\d{{2}}$", "", name
    )
    return suffix


def plot_per_instance_pd(
    log_groups: dict[str, list[str | Path]],
    output_dir: str | Path,
    label_map: dict[str, str] | None = None,
    style_map: dict[str, dict] | None = None,
    include_settings: list[str] | None = None,
    include_warmstart: bool = True,
) -> None:
    """Create one PD convergence plot per instance with all settings overlaid.

    Parameters
    ----------
    log_groups : mapping from setting name → list of Gurobi log file paths.
    output_dir : directory to save per-instance plots.
    label_map : optional mapping from setting name → display label.
    style_map : optional mapping from setting name → dict with 'color' key.
    include_warmstart : if False, trim each trajectory to start at the main
        solve (``restarts[-1]``) and re-zero its time axis at ``t=0``.
    """
    from collections import defaultdict

    setup_latex_style()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    def _line_style(i: int, s: str) -> dict:
        if style_map and s in style_map:
            return {k: v for k, v in style_map[s].items() if k in ("color", "linestyle")}
        return {"color": COLORS[i % len(COLORS)]}

    def _trim_and_extend(data: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return ``(t, p, d)`` ready to plot (warmstart trimmed if requested)."""
        t, p, d = data["times"], data["primals"], data["duals"]
        if not include_warmstart and data["restarts"]:
            t0 = data["restarts"][-1]
            m = t >= t0
            t, p, d = t[m] - t0, p[m], d[m]
        if len(t) == 0:
            return t, p, d
        t_end = (data["total_time"] - (data["restarts"][-1] if
                 (not include_warmstart and data["restarts"]) else 0.0))
        if t_end > 0 and t_end > t[-1]:
            obj = data["objective"]
            final_p = obj if (data["optimal"] and not np.isnan(obj)) else p[-1]
            final_d = obj if (data["optimal"] and not np.isnan(obj)) else d[-1]
            t = np.append(t, t_end)
            p = np.append(p, final_p)
            d = np.append(d, final_d)
        return t, p, d

    # Group files by instance across settings
    # instance_data[instance_name][setting] = Path
    instance_data: dict[str, dict[str, Path]] = defaultdict(dict)
    if include_settings is not None:
        present = set(log_groups.keys())
        sorted_settings = [s for s in include_settings if s in present]
        extra = sorted(s for s in present if s not in set(include_settings))
        sorted_settings.extend(extra)
    else:
        sorted_settings = sorted(log_groups.keys())
    for setting, files in log_groups.items():
        for f in files:
            f = Path(f)
            inst = _extract_instance_name(f, setting)
            instance_data[inst][setting] = f

    for instance, setting_files in sorted(instance_data.items()):
        fig, ax = plt.subplots(figsize=(7, 4.5))
        has_data = False

        for i, setting in enumerate(sorted_settings):
            if setting not in setting_files:
                continue
            data = read_log(setting_files[setting])
            t, p, d = _trim_and_extend(data)
            if len(t) == 0:
                continue

            ls = _line_style(i, setting)
            lbl = fmt_setting(setting, label_map)
            ax.plot(t, p, linewidth=1.5, label=f"{lbl} primal", **ls)
            ax.plot(t, d, linewidth=1.5, linestyle="--",
                    color=ls.get("color", "black"), alpha=0.5, label=f"{lbl} dual")
            has_data = True

        if not has_data:
            plt.close(fig)
            continue

        ax.set_title(instance.replace("_", r"\_"), fontsize=10)
        ax.set_xscale("log")
        ax.set_xlabel(r"Time (s)")
        ax.set_ylabel(r"Objective value")
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:,.0f}"))
        ax.legend(frameon=False, fontsize="small", loc="upper left",
                  bbox_to_anchor=(1.02, 1), borderaxespad=0)
        fig.tight_layout()
        fig.subplots_adjust(right=0.65)

        out_path = output_dir / f"{instance}.pdf"
        fig.savefig(out_path)
        fig.savefig(out_path.with_suffix(".png"))
        plt.close(fig)

    print(f"Saved {len(instance_data)} per-instance PD plots to {output_dir}")


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python gurobi_log.py <logfile1> [logfile2 ...] [--output path]")
        sys.exit(1)

    args = sys.argv[1:]
    out = None
    if "--output" in args:
        idx = args.index("--output")
        out = args[idx + 1]
        args = args[:idx] + args[idx + 2:]

    plot_primal_dual(args, output=out)
