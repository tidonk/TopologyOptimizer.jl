#!/usr/bin/env python3
"""Heatmap for the combined abs × rel SP-selection grid (`1-n_hits_abs_rel`).

To guarantee the cells **exactly match the table**, this does NOT recompute a
geomean from raw runs (which would diverge on the instance set and the
cut-error dropping/penalising).  Instead it reads the per-setting ``Geo. Mean``
straight out of the table's ``comparison_short_all.md`` and places each setting
at its grid coordinate:

    combined ``SP{a}SP{r}%``  -> (abs a, rel r)
    pure-abs ``SP{a}``        -> (abs a, rel 100)   (rel cap off)
    pure-rel ``SP{r}%``       -> (abs ∞, rel r)     (abs cap off; SP100% = solve all)

Axes therefore carry the full "off"/∞ edge on each side, and the no-cap corner
``(∞, 100 %)`` (``SP100%``) is the baseline — settings dropped from the table for
too many cut errors are simply absent, so a blank cell means "dropped", never a
silently different number.

Run this *after* ``compare_grid_abs_rel.py`` (the table must exist).
"""
import re
import sys
from pathlib import Path

_VISU = Path(__file__).resolve().parents[2]                      # visu/
sys.path.insert(0, str(_VISU))
from heatmap import grid_heatmap_from_values, geomeans_from_short_table   # noqa: E402

LEAF = _VISU.parents[2] / "check/results/paper/3-SP_selection/grid/1-n_hits_abs_rel"
MD = LEAF / "comparison_short_all.md"

_INF = 20            # sentinel abs level for "∞" (one tick past 10)
X_LEVELS = [1, 2, 5, 10, _INF]
Y_LEVELS = [1, 5, 10, 25, 50, 100]


def _coord(s: str):
    """Map a (de-LaTeX'd) table label to a grid (abs, rel) coordinate."""
    m = re.fullmatch(r"SP(\d+)SP(\d+)%", s)
    if m:
        return int(m.group(1)), int(m.group(2))          # combined
    m = re.fullmatch(r"SP(\d+)%", s)
    if m:
        return _INF, int(m.group(1))                     # pure rel -> abs ∞
    m = re.fullmatch(r"SP(\d+)", s)
    if m:
        return int(m.group(1)), 100                      # pure abs -> rel off (100)
    return None                                          # cf05_sp005 etc. -> not on grid


if not MD.is_file():
    print(f"warning: {MD} missing (run compare_grid_abs_rel.py first); skipping heatmap")
else:
    gm, nn, note = geomeans_from_short_table(MD, _coord)
    grid_heatmap_from_values(
        gm, nn, output=LEAF / "heatmap_time_solve_all.pdf",
        x_levels=X_LEVELS, y_levels=Y_LEVELS,
        baseline_xy=(_INF, 100),                         # the no-cap / solve-all corner
        xlabel=r"SP abs.\ hit-cap", ylabel=r"SP rel.\ hit-cap (\%)",
        xtick_labels=["1", "2", "5", "10", r"$\infty$"],
        ytick_labels=["1", "5", "10", "25", "50", "100"],
        n_note=note,
    )
    print(f"Wrote {LEAF / 'heatmap_time_solve_all.pdf'}  ({len(gm)} cells from table)")
