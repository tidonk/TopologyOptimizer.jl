"""Smoke test for selection.build_selection invariants.

Run with::

    python3 check/scripts/visu/test_selection.py

No pytest needed.  Exits non-zero on any failed assertion.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from selection import (  # noqa: E402
    BROKEN_STATUSES, RESOURCE_LIMIT_STATUSES, TIME_LIMIT,
    build_selection,
)


def _row(instance, settings, status, time_solve, tal=0, time_MP=None, time_SP_N1=None):
    return {
        "instance": instance, "settings": settings, "tal": tal,
        "status": status, "time_solve": time_solve,
        "time_MP": time_MP if time_MP is not None else max(time_solve / 2, 0.1),
        "time_SP_N1": time_SP_N1 if time_SP_N1 is not None else max(time_solve / 2, 0.1),
    }


def _build_test_df() -> pd.DataFrame:
    """5 instances x 3 settings, mix of statuses to exercise every branch.

    Layout (status table):

                     base   altA   altB
        case100_1   OPTIMAL OPTIMAL OPTIMAL    -> in I_opt (and I_ran, I_all)
        case100_2   OPTIMAL OPTIMAL TIME_LIMIT -> in I_ran   (not I_opt)
        case200_1   OPTIMAL TIME_L  TIME_LIMIT -> in I_ran   (not I_opt)
        case200_2   OPTIMAL CUT_ERR OPTIMAL    -> in I_all   (not I_ran)
        case300_1   OPTIMAL   --    OPTIMAL    -> in I_all   (altA missing)
    """
    rows = []
    rows.append(_row("pglib_scopf_case100_1", "base", "OPTIMAL",  120.0))
    rows.append(_row("pglib_scopf_case100_1", "altA", "OPTIMAL",  100.0))
    rows.append(_row("pglib_scopf_case100_1", "altB", "OPTIMAL",  150.0))

    rows.append(_row("pglib_scopf_case100_2", "base", "OPTIMAL",  200.0))
    rows.append(_row("pglib_scopf_case100_2", "altA", "OPTIMAL",  220.0))
    rows.append(_row("pglib_scopf_case100_2", "altB", "TIME_LIMIT", np.nan))

    rows.append(_row("pglib_scopf_case200_1", "base", "OPTIMAL",  500.0))
    rows.append(_row("pglib_scopf_case200_1", "altA", "TIME_LIMIT", np.nan))
    rows.append(_row("pglib_scopf_case200_1", "altB", "TIME_LIMIT", np.nan))

    rows.append(_row("pglib_scopf_case200_2", "base", "OPTIMAL",  800.0))
    rows.append(_row("pglib_scopf_case200_2", "altA", "CUT_ERROR", np.nan))
    rows.append(_row("pglib_scopf_case200_2", "altB", "OPTIMAL", 900.0))

    rows.append(_row("pglib_scopf_case300_1", "base", "OPTIMAL", 1100.0))
    # altA missing for case300_1
    rows.append(_row("pglib_scopf_case300_1", "altB", "OPTIMAL", 1050.0))
    return pd.DataFrame(rows)


def _check(condition: bool, msg: str) -> None:
    if not condition:
        print(f"FAIL: {msg}")
        sys.exit(1)
    print(f"  ok: {msg}")


def main() -> None:
    df = _build_test_df()
    sel = build_selection(
        df, settings=["base", "altA", "altB"], baseline="base",
        target="time_solve", min_buses=0, min_time=0.0,
    )

    print(sel.summary())

    _check(sel.settings == ["base", "altA", "altB"], "settings order preserved")
    _check(len(sel.listing) == 5, "listing has 5 instance groups (all cases)")
    _check(len(sel.ran) == 3, "ran has 3 groups (case100_1, case100_2, case200_1)")
    _check(len(sel.optimal) == 1, "optimal has 1 group (case100_1 only)")

    # Subset invariant
    listing_set = set(sel.listing.tolist())
    ran_set = set(sel.ran.tolist())
    opt_set = set(sel.optimal.tolist())
    _check(opt_set.issubset(ran_set), "optimal ⊆ ran")
    _check(ran_set.issubset(listing_set), "ran ⊆ listing")

    # error_keys: case200_2 + altA is the only CUT_ERROR.  Stored as the flat
    # (*group_key, setting) tuple — caller can `key = inst + (setting,)`.
    expected_err = {("pglib_scopf_case200_2", 0, "altA")}
    _check(sel.error_keys == expected_err, f"error_keys == {expected_err}")

    # TL clipping: case100_2/altB and case200_1/(altA,altB) should be TIME_LIMIT
    v = sel.vals
    _check(v.loc[("pglib_scopf_case100_2", 0), "altB"] == TIME_LIMIT,
           "TIME_LIMIT clipped to TIME_LIMIT scalar")
    _check(v.loc[("pglib_scopf_case200_1", 0), "altA"] == TIME_LIMIT,
           "second TIME_LIMIT cell clipped")

    # vals_for: 'optimal' must contain only case100_1
    vopt = sel.vals_for("optimal")
    _check(len(vopt) == 1 and vopt.notna().all().all(),
           "vals_for('optimal') has one row, all finite")

    # vals_for: 'ran' without penalize keeps NaN for non-resource-limit holes (none here)
    vran = sel.vals_for("ran")
    _check(len(vran) == 3, "vals_for('ran') has 3 rows")

    # vals_for: 'ran' with penalize → no NaN, max ≤ TIME_LIMIT
    vran_pen = sel.vals_for("ran", penalize_nonoptimal=True)
    _check(not vran_pen.isna().any().any(),
           "vals_for('ran', penalize_nonoptimal) has no NaN")
    _check((vran_pen <= TIME_LIMIT).all().all(),
           "all cells ≤ TIME_LIMIT after penalize")

    # status pivot present
    _check(sel.status.loc[("pglib_scopf_case200_2", 0), "altA"] == "CUT_ERROR",
           "status pivot preserves CUT_ERROR")

    # min_buses filter
    sel_big = build_selection(
        df, settings=["base", "altA", "altB"], baseline="base",
        min_buses=150,
    )
    _check(len(sel_big.listing) == 3,
           "min_buses=150 keeps only case200_1, case200_2, case300_1")

    # min_time filter
    sel_slow = build_selection(
        df, settings=["base", "altA", "altB"], baseline="base",
        min_time=300.0,
    )
    # max time_solve per instance: 100_1→150, 100_2→220, 200_1→500, 200_2→900, 300_1→1100
    _check(len(sel_slow.listing) == 3, "min_time=300 drops case100_*")

    print("\nALL SELECTION TESTS PASSED")


if __name__ == "__main__":
    main()
