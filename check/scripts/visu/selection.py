"""Canonical instance-selection used by every comparison artifact.

All tables (long / short, default / ``_all``), heatmaps (common-solved /
``_all``), performance profiles and time-decomposition plots derive their
"which instances do I aggregate?" decision from a single
:class:`InstanceSelection` built by :func:`build_selection`.

This guarantees that the statistics reported in one artifact (e.g. the
shifted geometric mean printed in the LaTeX table) are computed over the
exact same set of instances as those in a sibling artifact (e.g. the heatmap
cell or the short-table row), eliminating silent disagreement caused by each
script applying its own pre-cleanup.

Three nested sets are exposed:

* ``listing``  — every (instance[, tal]) group that survived the experiment-
  level filters (``min_buses``, ``min_time``).  These are the rows shown in
  long tables.  No requirement on per-setting status.
* ``ran``      — subset of ``listing`` where **every** chosen setting produced
  a valid result row (status not in ``BROKEN_STATUSES`` and not missing).
  This is the dagger set for the ``_all`` table variants and the inclusion
  set for ``_all`` heatmap statistics.
* ``optimal``  — subset of ``ran`` where **every** chosen setting reached
  ``OPTIMAL``.  This is the dagger set for the default table variants and
  the inclusion set for default heatmap statistics.

The selection also carries the pre-built status / value pivots so downstream
code never has to re-pivot.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import re

import numpy as np
import pandas as pd


TIME_LIMIT = 10800
BROKEN_STATUSES = {"INFEASIBLE", "INFEASIBLE_OR_UNBOUNDED", "DID_NOT_RUN"}
RESOURCE_LIMIT_STATUSES = {"TIME_LIMIT", "MEMORY_LIMIT"}


@dataclass(frozen=True)
class InstanceSelection:
    """Frozen view over an experiment, with canonical instance partitions.

    Attributes
    ----------
    settings
        Ordered list of setting names (baseline first when present).
    baseline
        Baseline setting name (may be absent from ``settings``).
    group_cols
        ``['instance']`` or ``['instance', 'tal']`` depending on the data.
    target
        Metric column the selection was built around (e.g. ``time_solve``).
    df
        Long-format dataframe, restricted to ``settings`` and to the
        ``listing`` instances.  All downstream code should consume this
        rather than the raw experiment dataframe.
    vals
        Pivot ``listing × settings`` of *target*.  ``TIME_LIMIT`` and
        ``MEMORY_LIMIT`` cells are clipped to :data:`TIME_LIMIT`; missing /
        broken cells are NaN.
    status
        Pivot ``listing × settings`` of the raw status column.
    listing, ran, optimal
        ``pd.Index`` over the group columns; ``optimal ⊆ ran ⊆ listing``.
    error_keys
        Set of ``(*group_key, setting)`` tuples whose status is in
        :data:`BROKEN_STATUSES`.  Used to render the ``--‡`` cell marker.
    """

    settings: list[str]
    baseline: str | None
    group_cols: list[str]
    target: str
    df: pd.DataFrame
    vals: pd.DataFrame
    status: pd.DataFrame
    listing: pd.Index
    ran: pd.Index
    optimal: pd.Index
    error_keys: set = field(default_factory=set)

    # ------------------------------------------------------------------ helpers
    def vals_for(self, mode: str, *, penalize_nonoptimal: bool = False) -> pd.DataFrame:
        """Return the *vals* slice to use for statistics.

        Parameters
        ----------
        mode
            One of ``'optimal'`` (default tables/heatmaps — strict) or
            ``'ran'`` (``_all`` tables/heatmaps — inclusive).
        penalize_nonoptimal
            For ``mode='ran'``: non-OPTIMAL cells are clipped to
            :data:`TIME_LIMIT` so every setting contributes equally.  Without
            this flag, non-OPTIMAL but resource-limited cells already carry
            ``TIME_LIMIT`` from :meth:`build_selection`.
        """
        if mode == "optimal":
            return self.vals.loc[self.optimal]
        if mode == "ran":
            sub = self.vals.loc[self.ran]
            if penalize_nonoptimal:
                # Any cell where this setting is not OPTIMAL on this row gets TL.
                st = self.status.loc[self.ran]
                sub = sub.where(st.eq("OPTIMAL"), TIME_LIMIT).fillna(TIME_LIMIT)
            return sub
        raise ValueError(f"unknown mode {mode!r}; expected 'optimal' or 'ran'")

    def index_set(self, mode: str) -> pd.Index:
        """Return the instance Index used by the given stats mode."""
        if mode == "optimal":
            return self.optimal
        if mode == "ran":
            return self.ran
        if mode == "listing":
            return self.listing
        raise ValueError(f"unknown mode {mode!r}")

    # ------------------------------------------------------------------ summary
    def summary(self) -> str:
        return (
            f"InstanceSelection(target={self.target!r}, settings={len(self.settings)}, "
            f"|listing|={len(self.listing)}, |ran|={len(self.ran)}, "
            f"|optimal|={len(self.optimal)})"
        )


# --------------------------------------------------------------------- builders
def _bus_count(name: str) -> int:
    m = re.search(r"case(\d+)", str(name))
    return int(m.group(1)) if m else 0


def _resolve_settings(df: pd.DataFrame, settings: list[str] | None,
                      baseline: str | None) -> list[str]:
    present = list(df["settings"].unique())
    if settings is not None:
        ordered = [s for s in settings if s in present]
    else:
        ordered = sorted(present)
        if baseline and baseline in ordered:
            ordered.remove(baseline)
            ordered.insert(0, baseline)
    return ordered


def _apply_min_buses(df: pd.DataFrame, min_buses: int) -> pd.DataFrame:
    if min_buses <= 0:
        return df
    mask = df["instance"].apply(_bus_count) >= min_buses
    return df[mask].reset_index(drop=True)


def _apply_min_time(df: pd.DataFrame, min_time: float,
                    target: str, group_cols: list[str]) -> pd.DataFrame:
    """Drop instance groups where no setting reaches at least *min_time*.

    Looks at ``target`` and at the sum of the standard time-decomposition
    columns when present, so trivially-fast instances are dropped no matter
    which time metric is being reported.
    """
    if min_time <= 0:
        return df
    candidates = [target]
    for c in ("time_MP", "time_SP_N1"):
        if c in df.columns and c not in candidates:
            candidates.append(c)
    grp_max = df.groupby(group_cols)[candidates].max().max(axis=1)
    keep_keys = grp_max[grp_max >= min_time].index
    if len(keep_keys) == len(grp_max):
        return df
    keep_df = pd.DataFrame(list(keep_keys), columns=group_cols)
    return df.merge(keep_df, on=group_cols).reset_index(drop=True)


def build_selection(
    df: pd.DataFrame,
    *,
    settings: list[str] | None = None,
    baseline: str | None = None,
    target: str = "time_solve",
    min_buses: int = 0,
    min_time: float = 0.0,
) -> InstanceSelection:
    """Build the canonical :class:`InstanceSelection` for an experiment.

    Steps
    -----
    1. Restrict ``df`` to the requested ``settings`` (sorted with baseline
       first when not specified).
    2. Drop instances below ``min_buses`` (case-size filter).
    3. Drop instance groups whose maximum ``target`` (or known
       time-decomposition column) is below ``min_time``.
    4. Pivot to (instance[, tal]) × settings for both ``target`` and
       ``status``.
    5. Clip TL/MEMORY_LIMIT cells to :data:`TIME_LIMIT` in the value pivot.
    6. Compute the three instance partitions and the ``error_keys`` set.

    The returned dataframe ``df`` (subset of the input) is what the rest of
    the pipeline should iterate over for plotting; the pivots in ``vals`` /
    ``status`` are what the statistics code should aggregate over.
    """
    df = df.copy()
    has_tal = "tal" in df.columns
    if has_tal:
        df["tal"] = df["tal"].fillna(0).astype(int)
    group_cols = ["instance", "tal"] if has_tal else ["instance"]

    settings_ordered = _resolve_settings(df, settings, baseline)
    df = df[df["settings"].isin(settings_ordered)].reset_index(drop=True)

    df = _apply_min_buses(df, min_buses)
    df = _apply_min_time(df, min_time, target, group_cols)

    if df.empty:
        empty_idx = pd.MultiIndex.from_tuples([], names=group_cols) if has_tal \
            else pd.Index([], name="instance")
        return InstanceSelection(
            settings=settings_ordered, baseline=baseline, group_cols=group_cols,
            target=target, df=df,
            vals=pd.DataFrame(columns=settings_ordered),
            status=pd.DataFrame(columns=settings_ordered),
            listing=empty_idx, ran=empty_idx, optimal=empty_idx,
            error_keys=set(),
        )

    vals = df.pivot_table(
        index=group_cols, columns="settings", values=target, aggfunc="first",
    ).reindex(columns=settings_ordered).apply(pd.to_numeric, errors="coerce")
    status = df.pivot_table(
        index=group_cols, columns="settings", values="status", aggfunc="first",
    ).reindex(columns=settings_ordered)

    # Align indices: pivot_table can drop index entries whose values are all
    # NaN, and ``vals`` vs ``status`` may diverge when an instance has only
    # error statuses (status non-NaN, target NaN) or vice versa. Reindex
    # both to the union so downstream boolean masks line up.
    union_idx = vals.index.union(status.index)
    vals = vals.reindex(union_idx)
    status = status.reindex(union_idx)

    # Clip resource-limited cells to TIME_LIMIT (numeric value may be missing
    # for these rows; we want a finite TL value in stats / heatmap arrays).
    tl_mask = status.isin(RESOURCE_LIMIT_STATUSES)
    vals = vals.mask(tl_mask, TIME_LIMIT)

    # Drop genuinely-infeasible instances entirely.  If the baseline (the
    # reference solver) reports a broken status (INFEASIBLE / ...), the instance
    # is out of scope for the whole comparison and is removed from *every*
    # artifact — the long per-instance table, the summary counts, plots — not
    # merely from the ran/optimal stats sets.  (Mirrors the historical
    # ``reader.cleanup_data`` behaviour, which the selection pipeline otherwise
    # no longer invokes.)
    if baseline in status.columns:
        infeasible_base = status[baseline].isin(BROKEN_STATUSES)
        if infeasible_base.any():
            keep_idx = status.index[~infeasible_base]
            n_dropped = int(infeasible_base.sum())
            vals = vals.loc[keep_idx]
            status = status.loc[keep_idx]
            df = df[df.set_index(group_cols).index.isin(keep_idx)].reset_index(drop=True)
            print(f"selection: dropped {n_dropped} genuinely-infeasible instance(s) "
                  f"(baseline '{baseline}' INFEASIBLE)")

    listing = vals.index
    # I_ran: status is OPTIMAL, TIME_LIMIT or MEMORY_LIMIT for every setting
    #        (i.e. not missing and not in BROKEN_STATUSES)
    ran_mask = status.notna().all(axis=1) & (~status.isin(BROKEN_STATUSES)).all(axis=1)
    ran = listing[ran_mask]
    optimal_mask = ran_mask & status.eq("OPTIMAL").all(axis=1)
    optimal = listing[optimal_mask]

    # Error keys: per-cell entries that are CUT_ERROR / INFEASIBLE / ...
    error_keys: set = set()
    err_pos = status.isin(BROKEN_STATUSES)
    if err_pos.any().any():
        stacked = err_pos.stack()
        # When the row index is a MultiIndex, ``stacked.index`` has
        # ``len(group_cols) + 1`` levels (last = setting).  Otherwise it has
        # exactly 2 (instance, setting).  Normalise both shapes to
        # ``(*group_key, setting)`` tuples.
        n_group = len(group_cols)
        for key in stacked[stacked].index:
            if not isinstance(key, tuple):
                key = (key,)
            group_key, s = key[:n_group], key[n_group]
            error_keys.add(tuple(group_key) + (s,))

    return InstanceSelection(
        settings=settings_ordered, baseline=baseline, group_cols=group_cols,
        target=target, df=df, vals=vals, status=status,
        listing=listing, ran=ran, optimal=optimal, error_keys=error_keys,
    )


__all__ = [
    "TIME_LIMIT",
    "BROKEN_STATUSES",
    "RESOURCE_LIMIT_STATUSES",
    "InstanceSelection",
    "build_selection",
]
