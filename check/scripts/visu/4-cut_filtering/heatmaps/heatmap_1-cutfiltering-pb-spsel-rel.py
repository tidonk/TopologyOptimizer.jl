#!/usr/bin/env python3
"""Heatmap for Exp. 4 Stage 1 - 1-cutfiltering-pb-spsel-rel (pb CF, rel SP)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _heatmap_common import render_heatmap  # noqa: E402

render_heatmap(
    "1-cutfiltering-pb-spsel-rel",
    baseline="sp100",
    is_pb=True,
)
