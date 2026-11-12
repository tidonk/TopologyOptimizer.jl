"""Comparison tables: instances × settings with summary statistics and optional Wilcoxon test.

Produces tables with:
  - Rows: one per instance (or instance+TAL pair)
  - Columns: one per setting
  - Cells: target metric value (e.g., time_solve)
  - Summary rows: Sum, Mean, Geometric Mean
  - Optional: Wilcoxon signed-rank test p-value and significance vs a baseline setting

All instance partitioning is delegated to :mod:`selection`, so the default
table, the ``_all`` table, the short tables and the heatmaps share **exactly**
the same set of instances in their statistics.

Two "modes" are exposed:

* ``mode='optimal'`` (default file names ``comparison_<target>.*``,
  ``comparison_short.*``).  Daggers and statistics are restricted to
  instances solved to ``OPTIMAL`` by every setting (``selection.optimal``).
* ``mode='ran'`` (file names ``comparison_<target>_all.*``,
  ``comparison_short_all.*``).  Daggers and statistics include every
  instance for which every setting produced a result CSV
  (``selection.ran``); non-OPTIMAL cells are clipped to ``TIME_LIMIT`` so
  every setting contributes equally.

Both modes list the same superset (``selection.listing``); the only thing
that changes between them is which rows are daggered and aggregated.
"""

from pathlib import Path
import re
import warnings

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from reader import read_experiment, read_runs, read_setting_sources, drop_excluded, compile_experiment_csv, instance_sort_key, clean_instance_name, cleanup_data, CUT_ERR_OBJ_TOLERANCE
from selection import (
    BROKEN_STATUSES, InstanceSelection, RESOURCE_LIMIT_STATUSES, build_selection,
)
from selection import TIME_LIMIT as _SELECTION_TIME_LIMIT

TIME_LIMIT = _SELECTION_TIME_LIMIT  # re-export for downstream `from tables import TIME_LIMIT`
GEO_SHIFT = 10.0  # shift applied to geometric means

SUMMARY_LABELS = {"Solved", "Solved by all", "Sum", "Mean", "Geo. Mean",
                  "Wilcoxon p", "Significance"}


def _shifted_geomean(col) -> float:
    """Shifted geometric mean ``exp(mean(log(v + s))) - s`` with ``s = GEO_SHIFT``.

    Drops NaN, clips negative values to 0. Returns NaN if no finite values.
    """
    v = pd.to_numeric(pd.Series(col), errors="coerce").dropna().clip(lower=0)
    if v.empty:
        return np.nan
    return float(np.exp(np.log(v + GEO_SHIFT).mean()) - GEO_SHIFT)


def _instance_key(row: pd.Series) -> str:
    """Build a unique row key from instance (+ tal if present)."""
    key = row["instance"]
    if "tal" in row.index and pd.notna(row.get("tal")):
        key = f"{key}_tal{int(row['tal'])}"
    return key


# --------------------------------------------------------------------------
# Selection-aware helpers
# --------------------------------------------------------------------------
def _group_to_display(group_key, group_cols: list[str]) -> str:
    """Map a (instance,) or (instance, tal) tuple to a display label."""
    if not isinstance(group_key, tuple):
        group_key = (group_key,)
    inst = clean_instance_name(str(group_key[0]))
    if len(group_cols) > 1 and len(group_key) > 1 and int(group_key[1]) > 0:
        return f"{inst}_tal{int(group_key[1])}"
    return inst


def _sorted_groups(index: pd.Index, group_cols: list[str]) -> list[tuple]:
    """Return the groups in ``index`` as tuples, sorted by instance size."""
    out: list[tuple] = []
    for key in index:
        if isinstance(key, tuple):
            out.append(tuple(key))
        else:
            out.append((key,))

    def _sort_key(k):
        return instance_sort_key(str(k[0])) + (k[1] if len(k) > 1 else 0,)

    return sorted(out, key=_sort_key)


def _error_marker_from(sel: InstanceSelection,
                       display_index: list[str],
                       group_order: list[tuple],
                       label_map: dict[str, str] | None) -> pd.DataFrame:
    """Build the boolean ``error_marker`` frame used by the formatters."""
    cols = [label_map.get(s, s) if label_map else s for s in sel.settings]
    em = pd.DataFrame(False, index=display_index, columns=cols)
    for disp, grp in zip(display_index, group_order):
        for raw_s, lbl_s in zip(sel.settings, cols):
            if (*grp, raw_s) in sel.error_keys:
                em.loc[disp, lbl_s] = True
    return em


def comparison_table(
    df: pd.DataFrame | None = None,
    target: str = "time_solve",
    baseline: str = "default",
    wilcoxon_test: bool = False,
    significance: float = 0.05,
    label_map: dict[str, str] | None = None,
    gap_map: dict[tuple, float] | None = None,
    error_keys: set | None = None,
    include_settings: list[str] | None = None,
    *,
    selection: InstanceSelection | None = None,
    mode: str = "optimal",
) -> pd.DataFrame:
    """Build a long instance-by-setting comparison table.

    Parameters
    ----------
    df
        Long-format experiment dataframe.  Ignored when *selection* is given.
    selection
        Pre-built :class:`~selection.InstanceSelection`.  If omitted, one is
        built from *df* with :func:`selection.build_selection`.
    mode
        ``'optimal'`` (default) — dagger and statistics restricted to
        ``selection.optimal``.  ``'ran'`` — dagger and statistics over
        ``selection.ran`` with non-OPTIMAL cells clipped to ``TIME_LIMIT``
        (the ``_all`` variant).
    target, baseline, include_settings
        Used only when *selection* is omitted, to build one on the fly.
    wilcoxon_test, significance
        One-sided Wilcoxon ``baseline > other`` test (lower-is-better) over
        the stats set.
    label_map
        ``{raw_setting -> display_label}`` applied to the column index.
    gap_map, error_keys
        Optional overrides; defaults come from *selection* when not given.

    Returns
    -------
    DataFrame with one row per instance (dagger when in the stats set), plus
    summary rows (``Solved``, ``Solved by all``, ``Sum``, ``Mean``,
    ``Geo. Mean`` and optionally ``Wilcoxon p`` / ``Significance``).
    """
    if mode not in {"optimal", "ran"}:
        raise ValueError(f"mode must be 'optimal' or 'ran', got {mode!r}")

    if selection is None:
        if df is None:
            raise ValueError("either df or selection must be provided")
        selection = build_selection(
            df, settings=include_settings, baseline=baseline, target=target,
        )

    settings = selection.settings
    group_cols = selection.group_cols
    target = selection.target

    # ---- Listing rows -------------------------------------------------
    group_order = _sorted_groups(selection.listing, group_cols)
    display_index = [_group_to_display(g, group_cols) for g in group_order]

    # Build the per-instance value table (sorted, with display labels)
    vals = selection.vals.reindex(group_order)
    vals.index = display_index

    # Stats set & dagger criterion
    stats_groups = selection.index_set(mode)
    stats_group_tuples = {tuple(g) if isinstance(g, tuple) else (g,)
                          for g in stats_groups}
    dagger_mask = pd.Series(
        [g in stats_group_tuples for g in group_order], index=display_index,
    )

    # For 'ran' mode, non-OPTIMAL cells in the stats rows are TL-clipped so
    # every setting contributes equally to the summary.
    if mode == "ran":
        stats_vals = selection.vals_for("ran", penalize_nonoptimal=True)
        stats_vals = stats_vals.reindex([g for g in group_order if g in stats_group_tuples])
        stats_vals.index = [_group_to_display(g, group_cols)
                            for g in group_order if g in stats_group_tuples]
    else:
        stats_vals = vals.loc[dagger_mask[dagger_mask].index]

    # ---- Dagger-decorate the display index ---------------------------
    daggered_index = [
        f"{idx}\u2020" if dagger_mask.loc[idx] else idx for idx in display_index
    ]
    vals_display = vals.copy()
    vals_display.index = daggered_index
    stats_vals.index = [
        f"{idx}\u2020" if not str(idx).endswith("\u2020") else idx
        for idx in stats_vals.index
    ]

    # Error marker frame (aligned to daggered index, label_map applied)
    error_marker = _error_marker_from(selection, display_index, group_order, label_map)
    error_marker.index = daggered_index
    if error_keys is not None:
        # Caller-supplied override: rebuild marker from those keys.
        em = pd.DataFrame(False, index=daggered_index,
                          columns=[label_map.get(s, s) if label_map else s
                                   for s in settings])
        for disp, grp in zip(daggered_index, group_order):
            for raw_s, lbl_s in zip(settings, em.columns):
                if (*grp, raw_s) in error_keys:
                    em.loc[disp, lbl_s] = True
        error_marker = em

    # ---- Per-setting "Solved" counts (informational) -----------------
    solved_count = pd.Series(index=settings, dtype=object)
    for s in settings:
        v = selection.vals[s]
        st = selection.status[s]
        n_solved = int(((v < TIME_LIMIT) & v.notna() & st.eq("OPTIMAL")).sum())
        solved_count[s] = str(n_solved)
    n_total = len(selection.listing)
    n_in_stats = int(dagger_mask.sum())

    counts = pd.DataFrame(index=["Solved", "Solved by all"], columns=settings)
    counts.loc["Solved"] = solved_count
    counts.loc["Solved by all"] = f"{n_in_stats}/{n_total}"

    summary = pd.DataFrame(index=["Sum", "Mean", "Geo. Mean"], columns=settings)
    summary.loc["Sum"] = stats_vals.sum()
    summary.loc["Mean"] = stats_vals.mean()
    summary.loc["Geo. Mean"] = stats_vals.apply(_shifted_geomean)

    # Wilcoxon (over the stats set, paired by instance)
    if wilcoxon_test and baseline in settings:
        baseline_vals = stats_vals[baseline].dropna()
        p_row = pd.Series(index=settings, dtype=float, name="Wilcoxon p")
        sig_row = pd.Series(index=settings, dtype=object, name="Significance")
        for s in settings:
            if s == baseline:
                p_row[s] = np.nan
                sig_row[s] = "(baseline)"
                continue
            other_vals = stats_vals[s].dropna()
            common = baseline_vals.index.intersection(other_vals.index)
            if len(common) < 3:
                p_row[s] = np.nan
                sig_row[s] = "n/a"
                continue
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", UserWarning)
                    _, p = wilcoxon(baseline_vals.loc[common], other_vals.loc[common],
                                    alternative="greater", zero_method="zsplit")
            except ValueError:
                p = np.nan
            p_row[s] = p
            if pd.isna(p):
                sig_row[s] = "n/a"
            elif p < 0.001:
                sig_row[s] = "***"
            elif p < 0.01:
                sig_row[s] = "**"
            elif p < significance:
                sig_row[s] = "*"
            else:
                sig_row[s] = ""
        summary = pd.concat([summary, pd.DataFrame([p_row, sig_row])])

    result = pd.concat([vals_display, counts, summary])

    # Apply label_map to column names
    if label_map:
        result = result.rename(columns={s: label_map.get(s, s) for s in result.columns})

    # Attrs consumed by ``to_markdown`` / ``to_latex``
    result.attrs["settings"] = settings
    result.attrs["target"] = target
    result.attrs["baseline"] = baseline
    result.attrs["gap_map"] = gap_map or {}
    result.attrs["error_keys"] = error_keys if error_keys is not None else selection.error_keys
    result.attrs["error_marker"] = error_marker
    result.attrs["mode"] = mode
    result.attrs["stats_n"] = n_in_stats
    result.attrs["listing_n"] = n_total
    # Group tuples (without dagger) that contributed to statistics
    result.attrs["stats_index"] = sorted(stats_group_tuples)
    result.attrs["excluded_from_stats"] = [
        str(idx).rstrip("\u2020") for idx in result.index
        if not str(idx).endswith("\u2020")
        and idx not in {"Solved", "Solved by all", "Sum", "Mean", "Geo. Mean",
                        "Wilcoxon p", "Significance"}
    ]
    return result


def _md_plain(md: str) -> str:
    """Passthrough: the markdown keeps the same labels as the LaTeX table.

    We used to strip ``\\texttt{X}`` -> ``X`` here, but ``to_markdown`` sizes the
    column padding/separators from the *original* label strings, so rewriting the
    cell text afterwards skewed the table alignment.  Keep the labels as-is.
    """
    return md


def to_markdown(table: pd.DataFrame) -> str:
    """Export comparison table as Markdown."""
    error_marker = table.attrs.get("error_marker", pd.DataFrame())
    _is_pdi = table.attrs.get("target") in ("pdi", "primal_integral")
    fmt = table.copy()
    for c in fmt.columns:
        fmt[c] = fmt[c].apply(
            lambda v: (f"{v:.2e}" if _is_pdi else f"{v:.2f}")
            if isinstance(v, float) and not pd.isna(v) else v
        )
    # Wilcoxon p row: 4 decimals
    if "Wilcoxon p" in fmt.index:
        for c in fmt.columns:
            v = table.loc["Wilcoxon p", c]
            if isinstance(v, float) and not pd.isna(v):
                fmt.loc["Wilcoxon p", c] = f"{v:.4f}"
    # Fill NaN: --‡ for error, -- for missing
    for idx in fmt.index:
        for c in fmt.columns:
            if pd.isna(fmt.loc[idx, c]):
                is_err = (not error_marker.empty
                          and idx in error_marker.index
                          and c in error_marker.columns
                          and error_marker.loc[idx, c])
                fmt.loc[idx, c] = "--‡" if is_err else "--"
    # Bold cells faster than baseline (first column)
    _summary_labels = {"Solved", "Solved by all", "Sum", "Mean", "Geo. Mean",
                       "Wilcoxon p", "Significance"}
    _cols = list(fmt.columns)
    if len(_cols) > 1:
        _bcol = _cols[0]
        for _idx in fmt.index:
            _is_stat = _idx in {"Sum", "Mean", "Geo. Mean"}
            _is_data = _idx not in _summary_labels
            if not (_is_stat or _is_data):
                continue
            _bv = _to_float(table.loc[_idx, _bcol])
            for _c in _cols[1:]:
                _fv = _to_float(table.loc[_idx, _c])
                if (_bv is not None and _fv is not None and _fv < _bv):
                    fmt.loc[_idx, _c] = f"**{fmt.loc[_idx, _c]}**"
    md = _md_plain(fmt.to_markdown())
    has_err = not error_marker.empty and error_marker.any().any()
    has_miss = fmt.isin(["--"]).any().any()
    notes = []
    if has_miss:
        notes.append("-- = not run")
    if has_err:
        notes.append("--‡ = solver error")
    if notes:
        md += "\n\n" + " | ".join(notes)
    return md


def to_latex(table: pd.DataFrame) -> str:
    """Export comparison table as LaTeX (booktabs), without pandas jinja2 dependency."""
    cols = table.columns.tolist()
    n_cols = len(cols)
    baseline = table.attrs.get("baseline")
    gap_map = table.attrs.get("gap_map", {})
    error_marker = table.attrs.get("error_marker", pd.DataFrame())
    _is_pdi = table.attrs.get("target") in ("pdi", "primal_integral")

    def _fmt_v(v) -> str:
        """Cell formatter; scientific notation for PDI / primal integral."""
        if pd.isna(v):
            return "--"
        if isinstance(v, float):
            return f"{v:.2e}" if _is_pdi else f"{v:.2f}"
        return _latex_escape(str(v))

    # Find baseline column (may be renamed via label_map)
    baseline_col = cols[0] if n_cols > 0 else None

    summary_labels = {"Solved", "Solved by all", "Sum", "Mean", "Geo. Mean", "Wilcoxon p", "Significance"}

    # Vertical line after baseline (first) column
    col_spec = "lr|" + "r" * (n_cols - 1) if n_cols > 1 else "lr"
    lines = [
        r"\begin{tabular}{" + col_spec + "}",
        r"\toprule",
        " & ".join([""] + [_latex_escape(c) for c in cols]) + r" \\",
        r"\midrule",
    ]

    data_rows = [idx for idx in table.index if idx not in summary_labels]
    summary_rows = [idx for idx in table.index if idx in summary_labels]

    for idx in data_rows:
        cells = [_latex_escape(str(idx))]
        baseline_val = _to_float(table.loc[idx, baseline_col]) if baseline_col else None
        for c in cols:
            v = table.loc[idx, c]
            fv = _to_float(v)
            # Gap annotation for TIME_LIMIT
            gap_key = (str(idx), c)
            gap = gap_map.get(gap_key)
            if fv is not None and fv >= TIME_LIMIT and gap is not None:
                if np.isinf(gap):
                    cell = r"($\infty$\%)"
                else:
                    cell = rf"({gap:.1f}\%)"
            elif fv is not None and fv >= TIME_LIMIT:
                cell = _fmt_v(v)
            else:
                cell = _fmt_v(v)
            # NaN: distinguish error from missing
            if pd.isna(v):
                is_err = (not error_marker.empty
                          and idx in error_marker.index
                          and c in error_marker.columns
                          and error_marker.loc[idx, c])
                cell = r"--\textsuperscript{\ddag}" if is_err else "--"
            # Bold if faster than baseline
            if (baseline_val is not None and fv is not None
                    and c != baseline_col
                    and fv < baseline_val and fv < TIME_LIMIT):
                cell = rf"\textbf{{{cell}}}"
            cells.append(cell)
        lines.append(" & ".join(cells) + r" \\")

    if summary_rows:
        count_rows = [r for r in summary_rows if r in {"Solved", "Solved by all"}]
        stat_rows = [r for r in summary_rows if r in {"Sum", "Mean", "Geo. Mean"}]
        test_rows = [r for r in summary_rows if r not in {"Solved", "Solved by all", "Sum", "Mean", "Geo. Mean"}]
        if count_rows:
            lines.append(r"\midrule")
            for idx in count_rows:
                cells = [_latex_escape(str(idx))]
                for c in cols:
                    v = table.loc[idx, c]
                    cells.append(_fmt_v(v))
                lines.append(" & ".join(cells) + r" \\")
        lines.append(r"\midrule")
        for idx in stat_rows:
            cells = [_latex_escape(str(idx))]
            baseline_val = _to_float(table.loc[idx, baseline_col]) if baseline_col else None
            for c in cols:
                v = table.loc[idx, c]
                fv = _to_float(v)
                cell = _fmt_v(v)
                if (baseline_val is not None and fv is not None
                        and c != baseline_col and fv < baseline_val):
                    cell = rf"\textbf{{{cell}}}"
                cells.append(cell)
            lines.append(" & ".join(cells) + r" \\")
        if test_rows:
            lines.append(r"\midrule")
            for idx in test_rows:
                cells = [_latex_escape(str(idx))]
                for c in cols:
                    v = table.loc[idx, c]
                    if idx == "Wilcoxon p" and not pd.isna(v):
                        cells.append(f"{float(v):.4f}")
                    else:
                        cells.append(_fmt_v(v))
                lines.append(" & ".join(cells) + r" \\")

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    return "\n".join(lines)


def _to_float(v) -> float | None:
    """Try converting *v* to float; return None on failure."""
    if pd.isna(v):
        return None
    try:
        return float(v)
    except (ValueError, TypeError):
        return None


def _wilcoxon_per_setting(
    vals: pd.DataFrame,
    settings: list[str],
    baseline: str,
    significance: float = 0.05,
) -> dict[str, tuple[float, str]]:
    """One-sided Wilcoxon (baseline > other, i.e. lower-is-better) per non-baseline setting.

    Returns ``{setting -> (p_value, significance_marker)}``. The baseline maps to
    ``(NaN, "(baseline)")``. Settings with fewer than 3 comparable pairs map to
    ``(NaN, "n/a")``.
    """
    out: dict[str, tuple[float, str]] = {}
    if baseline not in vals.columns:
        return {s: (np.nan, "n/a") for s in settings}
    baseline_vals = vals[baseline].dropna()
    for s in settings:
        if s == baseline:
            out[s] = (np.nan, "(baseline)")
            continue
        if s not in vals.columns:
            out[s] = (np.nan, "n/a")
            continue
        other = vals[s].dropna()
        common = baseline_vals.index.intersection(other.index)
        if len(common) < 3:
            out[s] = (np.nan, "n/a")
            continue
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                _, p = wilcoxon(baseline_vals.loc[common], other.loc[common],
                                alternative="greater", zero_method="zsplit")
        except ValueError:
            p = np.nan
        if pd.isna(p):
            sig = "n/a"
        elif p < 0.001:
            sig = "***"
        elif p < 0.01:
            sig = "**"
        elif p < significance:
            sig = "*"
        else:
            sig = ""
        out[s] = (float(p), sig)
    return out


def comparison_short_table(
    df: pd.DataFrame | None = None,
    target: str = "time_solve",
    baseline: str = "default",
    label_map: dict[str, str] | None = None,
    include_settings: list[str] | None = None,
    df_with_errors: pd.DataFrame | None = None,
    wilcoxon_test: bool = False,
    significance: float = 0.05,
    extra_err_by: dict | None = None,
    cut_err_drop: int | None = None,
    *,
    selection: InstanceSelection | None = None,
    mode: str = "optimal",
) -> pd.DataFrame:
    """Compact summary table (one row per setting).

    Columns: ``Solved | TL | CutErr | Sum | Mean | Geo. Mean`` (+ optional
    ``Wilcoxon p`` / ``Significance``).

    Parameters
    ----------
    selection
        Pre-built :class:`~selection.InstanceSelection`.  When omitted,
        one is built from *df*.
    mode
        ``'optimal'`` — stats over ``selection.optimal`` (strict).
        ``'ran'`` — stats over ``selection.ran`` with non-OPTIMAL clipped
        to ``TIME_LIMIT``.
    df_with_errors
        Optional pre-cleanup dataframe used only for *informational* TL /
        CutErr counts when *selection* doesn't already cover those rows.
        Normally omitted: counts come from the selection's status pivot.
    """
    if mode not in {"optimal", "ran"}:
        raise ValueError(f"mode must be 'optimal' or 'ran', got {mode!r}")

    if selection is None:
        if df is None:
            raise ValueError("either df or selection must be provided")
        selection = build_selection(
            df, settings=include_settings, baseline=baseline, target=target,
        )

    settings = selection.settings
    status = selection.status

    # Per-setting counts come from the *listing* status pivot — they describe
    # the data we actually have, not the stats subset.
    # CutErr counts *cut* errors only (a spurious/invalid Benders cut).  INFEASIBLE
    # instances are genuinely infeasible problems (e.g. the instances_grid_infeasible
    # set; the full-LP baseline also reports them infeasible) — they are excluded
    # from the comparison entirely (dropped from the `ran`/`optimal` sets by
    # BROKEN_STATUSES) and must NOT be tallied here, so they do not appear in any
    # of Solved / TL / CutErr.
    _ERR_STATUSES = {"CUT_ERROR"}
    n_solved_by = status.eq("OPTIMAL").sum()
    n_tl_by = status.isin(RESOURCE_LIMIT_STATUSES).sum()
    n_err_by = status.isin(_ERR_STATUSES).sum()

    # If the caller has a richer pre-cleanup view (e.g. instances dropped
    # earlier in run_full_comparison), prefer those counts for TL / CutErr
    # so the user can see what was filtered out.
    if df_with_errors is not None and "status" in df_with_errors.columns:
        ext = df_with_errors[df_with_errors["settings"].isin(settings)]
        ext_err = ext[ext["status"].isin(_ERR_STATUSES)].groupby("settings").size()
        ext_tl = ext[ext["status"].isin(RESOURCE_LIMIT_STATUSES)].groupby("settings").size()
        for s in settings:
            n_err_by[s] = max(int(n_err_by.get(s, 0)), int(ext_err.get(s, 0)))
            n_tl_by[s] = max(int(n_tl_by.get(s, 0)), int(ext_tl.get(s, 0)))

    stats_vals = selection.vals_for(mode, penalize_nonoptimal=(mode == "ran"))
    n_stats = len(stats_vals)

    wilcoxon_map = (_wilcoxon_per_setting(stats_vals, settings, baseline, significance)
                    if wilcoxon_test else {})

    rows = []
    baseline_label = None          # display label of the baseline row, if present
    for s in settings:
        col = stats_vals[s] if s in stats_vals.columns else pd.Series(dtype=float)
        row = {
            "Solved": str(int(n_solved_by.get(s, 0))),
            "TL": str(int(n_tl_by.get(s, 0))),
            "CutErr": str(int(n_err_by.get(s, 0)) + int((extra_err_by or {}).get(s, 0))),
            "Sum": col.sum(),
            "Mean": col.mean() if not col.empty else np.nan,
            "Geo. Mean": _shifted_geomean(col),
        }
        if wilcoxon_test:
            p, sig = wilcoxon_map.get(s, (np.nan, "n/a"))
            row["Wilcoxon p"] = p
            row["Significance"] = sig
        lbl = label_map.get(s, s) if label_map else s
        if s == baseline:
            baseline_label = lbl
        rows.append((lbl, row))

    result = pd.DataFrame(dict(rows)).T

    # Drop settings with more than ``cut_err_drop`` cut errors — too unreliable to
    # show — but never the baseline row.  The count feeds the footnote.
    n_cut_dropped = 0
    if cut_err_drop is not None and "CutErr" in result.columns:
        ce = pd.to_numeric(result["CutErr"], errors="coerce").fillna(0)
        keep = (ce <= cut_err_drop) | (result.index == baseline_label)
        n_cut_dropped = int((~keep).sum())
        result = result[keep]

    # Drop an all-zero CutErr column.  Cut-error instances are now *excluded*
    # from the test set rather than tallied (see reader.collect_cut_err_instances
    # / grid_cut_err_instances), so the column is normally empty and just adds
    # noise.  It reappears automatically if a real cut error ever survives.
    if "CutErr" in result.columns and (result["CutErr"].astype(str) == "0").all():
        result = result.drop(columns="CutErr")

    result.attrs["baseline"] = baseline
    # Display label of the baseline row, so formatters can highlight relative to
    # the *actual* baseline rather than assuming it is the first row (it often
    # is not — e.g. stage 4 deliberately places SP2 fifth).
    result.attrs["baseline_label"] = baseline_label
    result.attrs["n_instances"] = n_stats
    result.attrs["mode"] = mode
    if mode == "optimal":
        result.attrs["n_note"] = "instances (solved by all settings)"
    else:
        result.attrs["n_note"] = (
            f"instances (all; non-OPTIMAL/missing set to {TIME_LIMIT}s)"
        )
    if n_cut_dropped:
        result.attrs["n_note"] += (
            f"; {n_cut_dropped} setting(s) dropped for > {cut_err_drop} cut errors"
        )
    # Group tuples used in statistics — for cross-artifact consistency checks.
    stats_idx = selection.index_set(mode)
    result.attrs["stats_index"] = sorted(
        tuple(g) if isinstance(g, tuple) else (g,) for g in stats_idx
    )
    return result


def comparison_short_all_table(
    df: pd.DataFrame | None = None,
    target: str = "time_solve",
    baseline: str = "default",
    label_map: dict[str, str] | None = None,
    include_settings: list[str] | None = None,
    wilcoxon_test: bool = False,
    significance: float = 0.05,
    *,
    selection: InstanceSelection | None = None,
) -> pd.DataFrame:
    """Backward-compatible alias for ``comparison_short_table(mode='ran')``."""
    return comparison_short_table(
        df=df, target=target, baseline=baseline, label_map=label_map,
        include_settings=include_settings, wilcoxon_test=wilcoxon_test,
        significance=significance, selection=selection, mode="ran",
    )


def _clean_param_value(v: str) -> str:
    """Normalise a label value: ``1\\%`` -> ``0.01``, ``01`` -> ``1``, ``100\\%`` -> ``1``."""
    s = v.strip()
    pct = s.endswith(r"\%") or s.endswith("%")
    s = s.rstrip("%").rstrip("\\") if pct else s
    try:
        f = float(s)
    except ValueError:
        return v.strip()
    if pct:
        f /= 100.0
    return str(int(f)) if f == int(f) else f"{f:g}"


def _md_param_header(h: str) -> str:
    """Plain-text version of a LaTeX macro header, e.g. ``$\\spmaxn$`` -> ``spmaxn``."""
    return re.sub(r"[${}\\]", "", h).strip() or h


def label_to_param_cells(label: str, keys, default: str = r"$\infty$",
                         inf_when: dict | None = None) -> list[str]:
    """Expand a ``build_label()``-style label into ordered parameter cells.

    Grid labels encode the differing settings inline.  Two dialects occur::

        \\texttt{n: STAB=Inf, ON=F, FC=Inf, H=Inf, H\\%=0.5, PB=F}   (comma)
        \\texttt{CF=01 / SP=1\\%}                                     (slash)

    Given ``keys = ("H", r"H\\%", "PBB", "FC", "STAB")`` the first yields::

        ["$\\infty$", "0.5", "$\\infty$", "$\\infty$", "$\\infty$"]

    so the single label column can be replaced by proper parameter columns.
    Values are normalised (``1\\%`` -> ``0.01``, ``01`` -> ``1``).  ``Inf``, an
    absent key (the setting leaves that limit unset), or a value listed in
    ``inf_when[key]`` (e.g. ``CF=100`` meaning "no cut limit") renders as
    *default*.
    """
    body = label
    m = re.match(r"\s*\\texttt\{(.*)\}\s*$", body, flags=re.S)
    if m:
        body = m.group(1)
    if ":" in body:
        body = body.split(":", 1)[1]

    found: dict[str, str] = {}
    for part in re.split(r"[,/]", body):
        if "=" not in part:
            continue
        k, v = part.split("=", 1)
        found[k.strip()] = v.strip()

    inf_when = inf_when or {}
    cells = []
    for k in keys:
        v = found.get(k)
        sentinels = inf_when.get(k, ())
        if isinstance(sentinels, str):
            sentinels = (sentinels,)
        if (v is None or v in {"Inf", "inf", "Infinity", "None", ""}
                or v in tuple(sentinels)):
            cells.append(default)
        else:
            cells.append(_clean_param_value(v))
    return cells


def label_name_part(label: str) -> str:
    """Return the non-parameter part of a ``build_label()``-style label.

    ``\\texttt{CProx}, n=5`` -> ``\\texttt{CProx}``.  Used together with
    ``param_keep_label``: some paper tables keep a method-name column *and* add
    parameter columns (stage 83: ``\\texttt{CProx} & 5 & …``).  Labels without a
    parameter part are returned unchanged.

    The name itself may contain commas (``\\texttt{SP hybrid, no CF}, n=6``), so
    the split keeps every segment up to the first that contains ``=`` rather than
    just the first segment.
    """
    parts = re.split(r"[,/]", label)
    keep = []
    for part in parts:
        if "=" in part:                   # first parameter: the name ends here
            break
        keep.append(part)
    return ", ".join(k.strip() for k in keep).strip()


def comparison_short_to_latex(table: pd.DataFrame, param_columns=None,
                              inf_when: dict | None = None,
                              param_keep_label: bool = False) -> str:
    """Render the short comparison table as LaTeX (booktabs).

    ``param_columns`` optionally replaces the single label column with one column
    per setting parameter (the paper format).  It is a sequence of
    ``(header, key)`` or ``(header, key, table_format)`` tuples, where *key* is a
    ``build_label()`` abbreviation (``H``, ``H\\%``, ``FC``, ``STAB``, ``PBB``…)
    and *table_format* — if given — renders that column as a siunitx
    ``S[table-format=...]`` column.  Example::

        param_columns=[(r"$\\spmaxn$", "H"), (r"$\\spmaxrel$", r"H\\%", "3.2"),
                       (r"$\\sppb$", "PBB"), (r"$\\spminfeas$", "FC"),
                       (r"$\\spstab$", "STAB")]
    """
    cols = table.columns.tolist()
    # Split into count, stat, and test columns
    _COUNT = {"Solved", "Solved by all", "TL", "CutErr"}
    _STATS = {"Sum", "Mean", "Geo. Mean"}
    count_cols = [c for c in cols if c in _COUNT]
    stat_cols = [c for c in cols if c in _STATS]
    test_cols = [c for c in cols if c not in _COUNT and c not in _STATS]

    all_cols = count_cols + stat_cols + test_cols

    # Leading column(s): either the plain label, or one column per setting
    # parameter when *param_columns* is given (the paper format).
    if param_columns:
        _pc = [(t[0], t[1], (t[2] if len(t) > 2 else None)) for t in param_columns]
        lead_spec = "".join(f"S[table-format={fmt}]" if fmt else "r"
                            for _h, _k, fmt in _pc)
        # siunitx requires text headers of S columns to be braced.
        # Brace every parameter header (siunitx requires it for S columns;
        # the paper tables brace the others too).
        lead_head = [f"{{{h}}}" for h, _k, _fmt in _pc]
        if param_keep_label:
            # Keep the method-name column in front of the parameter columns.
            lead_spec = "l" + lead_spec
            lead_head = [""] + lead_head
    else:
        _pc = None
        lead_spec = "l"
        lead_head = [""]

    # Separate the parameter block from the count block, as in the paper tables.
    col_spec = lead_spec + ("|" if _pc else "") + "r" * len(count_cols)
    if stat_cols:
        col_spec += "|" + "r" * len(stat_cols)
    if test_cols:
        col_spec += "|" + "r" * len(test_cols)

    lines = [
        r"\begin{tabular}{" + col_spec + "}",
        r"\toprule",
        " & ".join(lead_head + [_latex_escape(c) for c in all_cols]) + r" \\",
        r"\midrule",
    ]

    # Highlight relative to the *baseline* row, not merely the first row: the
    # baseline is often not listed first (e.g. stage 4 places SP2 fifth), and
    # comparing against row 0 would bold the baseline itself against an
    # arbitrary setting.  Fall back to row 0 only if the label is unavailable.
    _baseline_idx = table.attrs.get("baseline_label")
    if _baseline_idx is None or _baseline_idx not in table.index:
        _baseline_idx = table.index[0] if len(table.index) > 0 else None
    for idx in table.index:
        if _pc:
            cells = label_to_param_cells(str(idx), [k for _h, k, _f in _pc],
                                         inf_when=inf_when)
            if param_keep_label:
                cells = [_latex_escape(label_name_part(str(idx)))] + cells
        else:
            cells = [_latex_escape(str(idx))]
        for c in count_cols:
            cells.append(str(table.loc[idx, c]))
        for c in stat_cols:
            v = table.loc[idx, c]
            cell = _fmt_cell(v)
            if _baseline_idx is not None and idx != _baseline_idx:
                bv = _to_float(table.loc[_baseline_idx, c])
                fv = _to_float(v)
                if bv is not None and fv is not None and fv < bv:
                    cell = rf"\textbf{{{cell}}}"
            cells.append(cell)
        for c in test_cols:
            v = table.loc[idx, c]
            if c == "Wilcoxon p" and not pd.isna(v):
                cells.append(f"{float(v):.4f}")
            else:
                cells.append(_fmt_cell(v))
        lines.append(" & ".join(cells) + r" \\")

    n = table.attrs.get("n_instances", "?")
    n_note = table.attrs.get("n_note", "instances (solved by all settings)")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(rf"% n = {n} {n_note}")
    return "\n".join(lines)


def comparison_short_to_markdown(table: pd.DataFrame, param_columns=None,
                                 inf_when: dict | None = None,
                                 param_keep_label: bool = False) -> str:
    """Render the short comparison table as Markdown.

    ``param_columns`` / ``inf_when`` behave exactly as in
    :func:`comparison_short_to_latex`: the single label column is replaced by one
    column per setting parameter (``$\\infty$`` rendered as ``inf`` for Markdown).
    """
    n = table.attrs.get("n_instances", "?")
    n_note = table.attrs.get("n_note", "instances (solved by all settings)")
    fmt = table.copy()
    if "Wilcoxon p" in fmt.columns:
        fmt["Wilcoxon p"] = fmt["Wilcoxon p"].apply(
            lambda v: f"{v:.4f}" if pd.notna(v) else "--")
    if "Significance" in fmt.columns:
        fmt["Significance"] = fmt["Significance"].fillna("--")
    # Bold stat cells better than the *baseline* row (which is often not the
    # first row — see comparison_short_to_latex for the rationale).
    _stat_cols = [c for c in fmt.columns if c in {"Sum", "Mean", "Geo. Mean"}]
    if len(fmt.index) > 0 and _stat_cols:
        # These columns receive formatted strings (possibly bolded); cast up
        # front so pandas does not warn about incompatible dtype assignment.
        fmt[_stat_cols] = fmt[_stat_cols].astype(object)
        _bidx = table.attrs.get("baseline_label")
        if _bidx is None or _bidx not in fmt.index:
            _bidx = fmt.index[0]
        for _c in _stat_cols:
            _bv = _to_float(table.loc[_bidx, _c])
            for _idx in fmt.index:
                _fv = _to_float(table.loc[_idx, _c])
                s = f"{_fv:.2f}" if _fv is not None else "--"
                if _idx != _bidx and _bv is not None and _fv is not None and _fv < _bv:
                    s = f"**{s}**"
                fmt.loc[_idx, _c] = s

    if param_columns:
        # Replace the label index with one column per setting parameter, so the
        # Markdown carries the same information as the LaTeX paper table.
        _pc = [(t[0], t[1]) for t in param_columns]
        _hdrs = [_md_param_header(h) for h, _k in _pc]
        _keys = [k for _h, k in _pc]
        expanded = [label_to_param_cells(str(i), _keys, default="inf",
                                         inf_when=inf_when) for i in fmt.index]
        params = pd.DataFrame(expanded, columns=_hdrs, index=fmt.index)
        if param_keep_label:
            params.insert(0, "", [label_name_part(str(i)) for i in fmt.index])
        fmt = pd.concat([params, fmt], axis=1).reset_index(drop=True)

    md = fmt.to_markdown(index=not bool(param_columns))
    md = _md_plain(md)
    md += f"\n\n*n = {n} {n_note}*"
    return md


    """Try to convert a table cell value to float, return None on failure."""
    if pd.isna(v):
        return None
    try:
        return float(v)
    except (ValueError, TypeError):
        return None


def _latex_escape(s: str) -> str:
    """Escape LaTeX special chars (safe for simple text)."""
    # Don't escape user-constructed LaTeX macros to avoid corruption.
    # For strings containing a backslash (e.g. ``\\texttt{rnd_rel_001}`` or
    # ``\\texttt{50\\%}``) we only escape *unescaped* ``%`` and ``_`` so
    # macro arguments render cleanly without producing ``\\\\%`` artefacts.
    if "\\" in s:
        import re as _re
        s = _re.sub(r"(?<!\\)%", r"\\%", s)
        s = _re.sub(r"(?<!\\)_", r"\\_", s)
        return s
    return (s.replace("%", r"\%")
            .replace("_", r"\_")
            .replace("&", r"\&")
            .replace("#", r"\#")
            .replace("$", r"\$")
            .replace("{", r"\{")
            .replace("}", r"\}")
            .replace("\u2020", r"$^{\dagger}$"))


def _fmt_cell(v) -> str:
    """Format a single table cell value."""
    if pd.isna(v):
        return "--"
    if isinstance(v, float):
        return f"{v:.2f}"
    return _latex_escape(str(v))


def write_tables(
    experiment_dir: str | Path,
    target: str = "time_solve",
    baseline: str = "default",
    wilcoxon_test: bool = False,
    output_dir: str | Path | None = None,
):
    """Read experiment, build table, write .md and .tex files."""
    experiment_dir = Path(experiment_dir)
    df = read_experiment(experiment_dir)

    has_baseline = baseline in df["settings"].unique()
    table = comparison_table(
        df, target=target, baseline=baseline,
        wilcoxon_test=wilcoxon_test and has_baseline,
    )

    if output_dir is None:
        output_dir = experiment_dir
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    md_path = output_dir / f"comparison_{target}.md"
    tex_path = output_dir / f"comparison_{target}.tex"
    md_path.write_text(to_markdown(table))
    tex_path.write_text(to_latex(table))
    print(f"Written: {md_path}, {tex_path}")
    return table


def _consistency_assert(*tables: pd.DataFrame, label: str = "") -> None:
    """Sanity-check that every supplied table aggregates over the same
    instance set.  Prints a warning (and raises in strict mode) when the
    ``stats_index`` attrs disagree.
    """
    indices = [t.attrs.get("stats_index") for t in tables if t is not None]
    indices = [i for i in indices if i is not None]
    if len(indices) < 2:
        return
    first = indices[0]
    for other in indices[1:]:
        if other != first:
            sym_diff = set(first) ^ set(other)
            print(f"WARNING [{label}]: stats_index disagreement across artifacts "
                  f"({len(sym_diff)} differing instance(s))")
            return
    print(f"consistency check [{label}]: {len(first)} instance(s), all artifacts agree.")


def _majority_objective_errors(obj_pivot: pd.DataFrame, settings: list[str],
                               baseline: str | None = None,
                               thresh: float = CUT_ERR_OBJ_TOLERANCE) -> dict[str, int]:
    """Flag per-setting objective deviations from the **majority** objective.

    For each instance the OPTIMAL objective values of all settings are clustered
    with relative tolerance *thresh*; the largest cluster is taken as the
    consensus ("majority opinion") ground truth.  Every setting outside that
    cluster is counted as a cut error for that instance.

    Unlike the previous baseline-relative rule the baseline enjoys **no special
    status**: if the baseline is the outlier, *it* gets flagged.  Ties between
    equally large clusters are broken in favour of the cluster containing the
    baseline (if any), then by the lower objective value.

    Returns ``{setting: n_instances_flagged}`` (non-zero entries only).
    """
    obj = obj_pivot.apply(pd.to_numeric, errors="coerce")
    cols = [s for s in settings if s in obj.columns]
    errs: dict[str, int] = {}
    if not cols:
        return errs

    for _idx, row in obj[cols].iterrows():
        vals = row.dropna()
        if len(vals) < 2:
            continue
        clusters: list[list] = []          # [representative_value, [settings]]
        for s, v in vals.items():
            for c in clusters:
                if abs(v - c[0]) / max(abs(c[0]), 1e-10) <= thresh:
                    c[1].append(s)
                    break
            else:
                clusters.append([float(v), [s]])
        if len(clusters) == 1:
            continue                       # unanimous — nothing to flag
        best = max(len(c[1]) for c in clusters)
        cands = [c for c in clusters if len(c[1]) == best]
        if len(cands) > 1 and baseline is not None:
            with_bl = [c for c in cands if baseline in c[1]]
            if with_bl:
                cands = with_bl
        consensus = min(cands, key=lambda c: c[0])
        for c in clusters:
            if c is consensus:
                continue
            for s in c[1]:
                errs[s] = errs.get(s, 0) + 1
    return errs


def run_full_comparison(
    experiment_dir,
    baseline: str = "default",
    target: str = "time_solve",
    output_dir: str | Path | None = None,
    exclude_instances: set | None = None,
    param_columns=None,
    inf_when: dict | None = None,
    param_keep_label: bool = False,
    pd_convergence: bool = True,
    style_map: dict[str, dict] | None = None,
    label_map: dict[str, str] | None = None,
    min_buses: int = 0,
    include_settings: list[str] | None = None,
    min_time: float = 0.0,
    log_raw_cache: dict | None = None,
    log_files_map: dict | None = None,
    plot_per_instance_time_bar: bool = False,
    plot_cuts_boxplot: bool = False,
    perf_profile_top_n: int = 5,
    pd_top_n: int | None = None,
    pd_include_warmstart: bool = True,
    pd_per_instance: bool = True,
    pd_scopf_thresholds: dict[str, float] | None = None,
    df_transform=None,
    cut_err_drop: int | None = None,
):
    """Run the complete comparison pipeline.

    ``experiment_dir`` may be:

    * a single directory (``str``/``Path``) — legacy behaviour, read via
      :func:`read_experiment`;
    * a **list of repetition directories** — read and averaged per instance via
      :func:`read_runs` (arithmetic mean, non-OPTIMAL ``target`` penalised to
      ``TIME_LIMIT``, status OPTIMAL only if every repetition was).  Use this for
      "full run" comparisons that span several dispatched repetitions; or
    * a **dict** ``{setting: [run dirs]}`` — a per-setting source map, read via
      :func:`read_setting_sources`.  Use this when the settings of one comparison
      were produced by different dispatch runs (possibly in different stage
      repos): every setting is re-aggregated from its own original per-run CSVs,
      so nothing depends on pre-computed mean CSVs copied into a hand-assembled
      directory.

    ``output_dir`` selects where every artifact is written.  When omitted it
    defaults to the (first) source directory, preserving legacy behaviour.  Pass
    an explicit ``output_dir`` to write into e.g.
    ``results/across_all/<stage>/full/<group>/`` without touching the source
    result directories.

    ``exclude_instances`` is an optional set of instance names to drop before any
    table/plot is built — pass the **test-set-appropriate** CUT_ERROR list
    (``reader.load_excluded_instances(CUT_ERR_EXCLUSIONS_FILE)`` for full runs,
    ``…_FILE_GRID`` for grid runs); never mix the two (same names, different
    problems).

    One :class:`~selection.InstanceSelection` is built and shared by every
    artifact (long tables, short tables, heatmaps when invoked downstream,
    performance profiles, time decompositions, PD-convergence plots).  This
    guarantees that the statistics in the LaTeX table, the values in a
    heatmap cell and the geometric means in the short summary are computed
    over the same set of instances.

    Two table variants are produced:

    * ``comparison_<target>.{md,tex}`` + ``comparison_short.{md,tex}`` —
      strict (``mode='optimal'``).
    * ``comparison_<target>_all.{md,tex}`` + ``comparison_short_all.{md,tex}``
      — inclusive (``mode='ran'``, non-OPTIMAL clipped to ``TIME_LIMIT``).

    Plot artifacts (status bar, performance profile, time decomposition,
    boxplots) all consume the same restricted dataframe ``selection.df``.
    """
    import re as _re
    from perf_profile import performance_profile
    from time_bar import per_instance_time_bar, status_bar_chart, cuts_boxplot, time_boxplot
    from gurobi_log import plot_aggregated_pd, plot_per_instance_pd

    # Resolve source directory(ies) and the (possibly separate) output dir.
    if isinstance(experiment_dir, dict):
        # Per-setting sources: each setting is re-aggregated from its own real
        # repetition dirs (never from copied-in mean CSVs).  A value may be a
        # plain list of dirs, or a ``(dirs, source_folder)`` pair for settings
        # whose runs live under a differently-named folder.
        def _dirs_of(v):
            if isinstance(v, tuple) and len(v) == 2 and isinstance(v[1], str):
                return v[0]
            return v
        src_dirs = [Path(d) for v in experiment_dir.values() for d in _dirs_of(v)]
        raw_df = read_setting_sources(experiment_dir, target=target)
    elif isinstance(experiment_dir, (list, tuple)):
        src_dirs = [Path(d) for d in experiment_dir]
        raw_df = read_runs(src_dirs, average=True, target=target)
    else:
        src_dirs = [Path(experiment_dir)]
        raw_df = read_experiment(src_dirs[0])
    out_dir = Path(output_dir) if output_dir is not None else src_dirs[0]
    out_dir.mkdir(parents=True, exist_ok=True)

    # Drop the caller-supplied CUT_ERROR exclusion list (test-set specific: the
    # caller passes the full-set or grid-set list, never both).
    if exclude_instances and "instance" in raw_df.columns:
        n0 = raw_df["instance"].nunique()
        raw_df = drop_excluded(raw_df, exclude_instances)
        n1 = raw_df["instance"].nunique()
        if n0 != n1:
            print(f"cut-err exclusion: dropped {n0 - n1} instance(s)")

    # Central derived column: time_with_warmstart = time_solve + time_warmstart
    # (time_warmstart is 0 for non-warmstarted settings, so it's a no-op there).
    # Defined here once so no compare script needs its own transform — just pass
    # ``target="time_with_warmstart"``.
    if "time_with_warmstart" not in raw_df.columns:
        # A missing column must default to a scalar 0 — ``DataFrame.get(col, 0)``
        # returns the int ``0`` (no ``.fillna``) when the column is absent, e.g.
        # cold FINAL data has no ``time_warmstart`` column at all.
        _ts = (pd.to_numeric(raw_df["time_solve"], errors="coerce").fillna(0)
               if "time_solve" in raw_df.columns else 0)
        _tw = (pd.to_numeric(raw_df["time_warmstart"], errors="coerce").fillna(0)
               if "time_warmstart" in raw_df.columns else 0)
        raw_df = raw_df.assign(time_with_warmstart=_ts + _tw)

    # Optional caller-supplied transform for any *other* derived columns.
    if df_transform is not None:
        raw_df = df_transform(raw_df)

    # Build the canonical instance selection.  This collapses what used to
    # be several ad-hoc filters (min_buses, missing-CSV pruning, min_time,
    # baseline-bad pruning) into one well-defined object.
    selection = build_selection(
        raw_df,
        settings=include_settings,
        baseline=baseline,
        target=target,
        min_buses=min_buses,
        min_time=min_time,
    )
    print(f"build_selection: {selection.summary()}")

    df = selection.df
    settings = selection.settings
    group_cols = selection.group_cols
    has_baseline = baseline in settings

    # ---- Side artifacts that need a status-rich pre-cleanup view -------
    # The status distribution chart shows every status (including missing
    # ones via DID_NOT_RUN), so it operates on the raw, pre-selection df
    # restricted to the chosen settings.
    raw_in_settings = raw_df[raw_df["settings"].isin(settings)].reset_index(drop=True)

    # Compiled CSV (raw, before any filtering).  For a single source dir keep
    # the legacy re-read so the file is byte-identical to before; for averaged
    # multi-run input persist the already-merged frame.
    if len(src_dirs) == 1:
        compile_experiment_csv(src_dirs[0], out_dir / "compiled.csv")
    else:
        raw_df.to_csv(out_dir / "compiled.csv", index=False)

    status_bar_chart(raw_in_settings, output=out_dir / "status_distribution.pdf",
                     label_map=label_map, baseline=baseline,
                     include_settings=include_settings)

    # ---- Objective consistency check (informational) ------------------
    obj_err_by: dict[str, int] = {}
    if "objective" in df.columns and len(settings) > 1:
        opt_df = df[df["status"] == "OPTIMAL"].copy()
        obj_pivot = opt_df.pivot_table(index=group_cols, columns="settings",
                                       values="objective", aggfunc="first")
        obj_complete = obj_pivot.dropna()
        if len(obj_complete) > 0:
            obj_vals = obj_complete.apply(pd.to_numeric, errors="coerce")
            obj_range = obj_vals.max(axis=1) - obj_vals.min(axis=1)
            rel_range = obj_range / obj_vals.abs().max(axis=1)
            thresh = CUT_ERR_OBJ_TOLERANCE
            bad = rel_range[rel_range > thresh]
            if len(bad) > 0:
                print(f"\nWARNING: {len(bad)} instance(s) have objective discrepancies (>{thresh*100:.2f}% rel):")
                short_settings = [s.replace("cutsel_", "").replace("kmeans_", "km_")
                                  for s in obj_vals.columns]
                col_w = max(12, max(len(s) for s in short_settings))
                row_w = max(len(idx if isinstance(idx, str) else "_".join(str(x) for x in idx))
                            for idx in bad.index)
                header = " " * (row_w + 2) + "".join(f"{s:>{col_w}}" for s in short_settings) + f"{'rel_range':>{col_w}}"
                print(header)
                print("-" * len(header))
                for idx in bad.index:
                    name = idx if isinstance(idx, str) else "_".join(str(x) for x in idx)
                    vals_s = "".join(f"{obj_vals.loc[idx, s]:>{col_w}.2f}" for s in obj_vals.columns)
                    vals_s += f"{rel_range.loc[idx]:>{col_w}.6f}"
                    print(f"{name:<{row_w}}  {vals_s}")
                print()
            else:
                print(f"objective check: {len(obj_complete)} all-OPTIMAL instances match (max rel range {rel_range.max():.2e})")

        # Per-setting count of instances whose OPTIMAL objective deviates from the
        # **majority** objective across settings (consensus ground truth, rather
        # than trusting the baseline).  Reported as CutErr in the summary tables
        # because a wrong objective signals a spurious cut.  The baseline itself
        # is flagged when it is the outlier.
        obj_err_by = _majority_objective_errors(obj_pivot, settings, baseline)

    # ---- Gap map (for TIME_LIMIT annotations in the long table) -------
    gap_map = {}
    tl_rows = df[df["status"] == "TIME_LIMIT"]
    if len(tl_rows) > 0:
        _parse_gurobi_log = None
        for _, row in tl_rows.iterrows():
            setting = row["settings"]
            inst = row["instance"]
            inst_key = clean_instance_name(inst)
            if "tal" in row.index and pd.notna(row.get("tal")) and int(row["tal"]) > 0:
                inst_key = f"{inst_key}_tal{int(row['tal'])}"
            col_name = label_map.get(setting, setting) if label_map else setting

            if log_raw_cache is not None and log_files_map is not None:
                log_file = log_files_map.get((inst, setting))
                if log_file is not None and log_file in log_raw_cache:
                    data = log_raw_cache[log_file]
                    gap_map[(inst_key, col_name)] = (
                        data["gap"] if not np.isnan(data["gap"]) else float("inf")
                    )
                    continue
            if _parse_gurobi_log is None:
                from gurobi_log import read_log as _parse_gurobi_log
            # Logs may live in any of the source repetition dirs; use the first
            # match found across them.
            _found = False
            for _sd in src_dirs:
                log_dir = _sd / setting / "logs"
                if not log_dir.is_dir():
                    continue
                for log_file in log_dir.glob("GRB_MP_*.out"):
                    if inst.replace("pglib_opf_", "") in log_file.name or inst in log_file.name:
                        data = _parse_gurobi_log(log_file)
                        gap_map[(inst_key, col_name)] = (
                            data["gap"] if not np.isnan(data["gap"]) else float("inf")
                        )
                        _found = True
                        break
                if _found:
                    break

    # ---- Tables: strict (mode=optimal) and inclusive (mode=ran) -------
    long_opt = comparison_table(
        selection=selection, baseline=baseline,
        wilcoxon_test=has_baseline, label_map=label_map,
        gap_map=gap_map, mode="optimal",
    )
    (out_dir / f"comparison_{target}.md").write_text(to_markdown(long_opt))
    (out_dir / f"comparison_{target}.tex").write_text(to_latex(long_opt))

    long_all = comparison_table(
        selection=selection, baseline=baseline,
        wilcoxon_test=has_baseline, label_map=label_map,
        gap_map=gap_map, mode="ran",
    )
    (out_dir / f"comparison_{target}_all.md").write_text(to_markdown(long_all))
    (out_dir / f"comparison_{target}_all.tex").write_text(to_latex(long_all))

    short_opt = comparison_short_table(
        selection=selection, baseline=baseline, label_map=label_map,
        df_with_errors=raw_in_settings, mode="optimal",
        wilcoxon_test=has_baseline, extra_err_by=obj_err_by,
        cut_err_drop=cut_err_drop,
    )
    (out_dir / "comparison_short.md").write_text(
        comparison_short_to_markdown(short_opt, param_columns=param_columns, inf_when=inf_when,
                                     param_keep_label=param_keep_label))
    (out_dir / "comparison_short.tex").write_text(
        comparison_short_to_latex(short_opt, param_columns=param_columns, inf_when=inf_when,
                                     param_keep_label=param_keep_label))

    short_all = comparison_short_table(
        selection=selection, baseline=baseline, label_map=label_map,
        df_with_errors=raw_in_settings, mode="ran",
        wilcoxon_test=has_baseline, extra_err_by=obj_err_by,
        cut_err_drop=cut_err_drop,
    )
    (out_dir / "comparison_short_all.md").write_text(
        comparison_short_to_markdown(short_all, param_columns=param_columns, inf_when=inf_when,
                                     param_keep_label=param_keep_label))
    (out_dir / "comparison_short_all.tex").write_text(
        comparison_short_to_latex(short_all, param_columns=param_columns, inf_when=inf_when,
                                     param_keep_label=param_keep_label))

    # Cross-artifact consistency: the long & short tables for the same mode
    # must aggregate over the exact same instance set.
    _consistency_assert(long_opt, short_opt, label="mode=optimal")
    _consistency_assert(long_all, short_all, label="mode=ran")

    # Performance profile: limit to baseline + top-N non-baseline by
    # shifted-geomean over the *ran* ("_all") partition (TIME_LIMIT clipped
    # at TIME_LIMIT) so the ranking reflects the same instance set users see
    # in the long/short _all tables.
    pp_settings: list[str] | None = include_settings
    if perf_profile_top_n is not None and perf_profile_top_n > 0:
        ran_vals = selection.vals_for("ran", penalize_nonoptimal=True)
        ranked: list[tuple[float, str]] = []
        for s in settings:
            col = ran_vals[s].dropna() if s in ran_vals.columns else pd.Series(dtype=float)
            gm = _shifted_geomean(col.clip(lower=0).to_numpy()) if not col.empty else float("inf")
            ranked.append((gm, s))
        non_base = [s for _, s in sorted(ranked) if s != baseline][:perf_profile_top_n]
        pp_settings = ([baseline] if has_baseline else []) + non_base

    performance_profile(df, target=target,
                        output=out_dir / f"performance_profile_{target}.pdf",
                        style_map=style_map, label_map=label_map,
                        baseline=baseline, include_settings=pp_settings)
    performance_profile(df, target=target, mode="time",
                        output=out_dir / f"performance_profile_{target}_cumulative.pdf",
                        style_map=style_map, label_map=label_map,
                        baseline=baseline, include_settings=pp_settings)

    # Per-setting helpers (only instances where every setting reached OPTIMAL).
    all_ok = df.groupby(group_cols)["status"].apply(
        lambda s: (s == "OPTIMAL").all()
    )
    ok_keys = all_ok[all_ok].reset_index()[group_cols]
    df_time = df.merge(ok_keys, on=group_cols)

    if plot_per_instance_time_bar:
        per_instance_time_bar(df, output=out_dir / "time_decomposition_per_instance.pdf",
                              label_map=label_map, axes_limits=[300, 3600],
                              baseline=baseline, include_settings=include_settings)

    if plot_cuts_boxplot:
        cuts_boxplot(df_time, output=out_dir / "cuts_boxplot.pdf",
                     label_map=label_map, baseline=baseline,
                     include_settings=include_settings)

    # Time boxplot: MP vs SP per setting (always on; cheap and informative).
    # Skip when MP/SP component columns are missing (e.g., compact-only experiments).
    if {"time_MP", "time_SP_N1"}.issubset(df_time.columns):
        time_boxplot(df_time, output=out_dir / "time_boxplot.pdf",
                     label_map=label_map, baseline=baseline,
                     include_settings=include_settings)
    else:
        print("[info] skipping time_boxplot: missing MP/SP component columns")

    # Aggregated PD convergence (master problem logs only)
    if pd_convergence:
        def _log_bus_count(p: Path) -> int:
            m = _re.search(r"case(\d+)", p.name)
            return int(m.group(1)) if m else 0

        # Only include logs for instances that survived cleanup (same as tables:
        # exclude instances where any setting has a BROKEN status).
        clean_instances = set(selection.ran.get_level_values("instance"))

        def _log_matches_clean(p: Path) -> bool:
            """Check if a log file belongs to an instance in the cleaned dataset."""
            for inst in clean_instances:
                if inst.replace("pglib_opf_", "") in p.name or inst in p.name:
                    return True
            return False

        # Per-instance SCOPF optimal objectives from compiled CSV, used to detect
        # and exclude logs where the Gurobi MP objective is inflated by the
        # Benders auxiliary variable η (i.e. η did not converge to 0 in that run).
        true_objectives: dict[str, float] = {}
        if "objective" in df.columns:
            _opt_rows = df[(df["status"] == "OPTIMAL") & df["objective"].notna()].copy()
            _opt_rows = _opt_rows.assign(_obj_num=pd.to_numeric(_opt_rows["objective"], errors="coerce"))
            for _inst, _grp in _opt_rows.groupby("instance"):
                _v = _grp["_obj_num"].min()
                if np.isfinite(_v):
                    true_objectives[_inst] = float(_v)

        # Gather MP logs for each setting across all source repetition dirs.
        log_groups = {}
        for s in settings:
            files = []
            for _sd in src_dirs:
                ld = _sd / s / "logs"
                if ld.is_dir():
                    files += [f for f in ld.glob("GRB_MP_*.out")
                              if _log_bus_count(f) >= min_buses and _log_matches_clean(f)]
            if files:
                log_groups[s] = sorted(files)
        if log_groups:
            plot_aggregated_pd(log_groups, output=out_dir / "pd_convergence_aggregated.pdf",
                               label_map=label_map, style_map=style_map,
                               include_settings=include_settings,
                               top_n=pd_top_n, baseline=baseline,
                               include_warmstart=pd_include_warmstart,
                               scopf_thresholds=pd_scopf_thresholds,
                               true_objectives=true_objectives or None)
            if pd_per_instance:
                plot_per_instance_pd(log_groups,
                                     output_dir=out_dir / "pd_convergence_plots",
                                     label_map=label_map, style_map=style_map,
                                     include_settings=include_settings,
                                     include_warmstart=pd_include_warmstart)
        else:
            print("No per-setting MP logs found, skipping PD convergence.")


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python tables.py <experiment_dir> [target] [baseline] [--wilcoxon]")
        sys.exit(1)
    exp = sys.argv[1]
    target = sys.argv[2] if len(sys.argv) > 2 else "time_solve"
    baseline = sys.argv[3] if len(sys.argv) > 3 else "default"
    do_wilcoxon = "--wilcoxon" in sys.argv
    write_tables(exp, target=target, baseline=baseline, wilcoxon_test=do_wilcoxon)
