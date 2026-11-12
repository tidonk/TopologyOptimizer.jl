#!/usr/bin/env python3
"""Heatmap for Exp. 6 Stage 1 - 1-tal1 (pool size n x weight scheme w)."""

import sys
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _heatmap_common import render_heatmap  # noqa: E402


def _add_time_with_warmstart(df: pd.DataFrame) -> pd.DataFrame:
    """Derive time_with_warmstart = time_solve + time_warmstart.

    All settings except ``no_warmstart`` include warmstart overhead.
    """
    df = df.copy()
    t_solve = pd.to_numeric(df.get("time_solve", 0), errors="coerce").fillna(0)
    t_ws = pd.to_numeric(df.get("time_warmstart", 0), errors="coerce").fillna(0)
    # include warmstart time for EVERY run (time_warmstart is 0 for
    # non-warmstarted settings, so this is a no-op there)
    df["time_with_warmstart"] = t_solve + t_ws
    return df


render_heatmap(
    "6-seq_warmstart--1-tal1",
    baseline="no_warmstart",
    output_name="heatmap_time_with_warmstart.pdf",
    target="time_with_warmstart",
    df_transform=_add_time_with_warmstart,
)

render_heatmap(
    "6-seq_warmstart--1-tal1-pd",
    baseline="no_warmstart",
    output_name="heatmap_time_with_warmstart.pdf",
    target="time_with_warmstart",
    df_transform=_add_time_with_warmstart,
)
