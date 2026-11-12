#!/usr/bin/env python3
"""Heatmap for Exp. 4 Stage 1 - 2-cutfiltering-pb-spsel-abs (pb CF, abs SP)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _heatmap_common import render_heatmap  # noqa: E402

render_heatmap(
    "2-cutfiltering-pb-spsel-abs",
    baseline="sp10",
    is_pb=True,
)
