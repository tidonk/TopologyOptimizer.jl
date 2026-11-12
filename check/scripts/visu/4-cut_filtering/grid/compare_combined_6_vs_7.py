#!/usr/bin/env python3
"""Combined Exp 6 vs Exp 7 cut-filtering comparison.

Concatenates the compiled.csv of two experiments, prefixes each setting with
its experiment tag, then runs ``comparison_short_all_table`` (TL-penalized,
shifted geomean) so settings from both experiments rank on a fair common
instance set.

Outputs to ``check/results/<combined_name>/comparison_short_all.{md,tex}``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO_ROOT / "check" / "scripts" / "visu"))

from tables import (  # noqa: E402
    comparison_short_all_table,
    comparison_short_to_latex,
    comparison_short_to_markdown,
)

RESULTS = REPO_ROOT / "check" / "results"

EXPERIMENTS = [
    ("E6", "6-cutfiltering-abs-spsel-rel"),
    ("E7", "7-cutfiltering-rel-spsel-opt"),
]
OUT_NAME = "combined_6_vs_7"


def main() -> None:
    frames = []
    for tag, exp in EXPERIMENTS:
        csv = RESULTS / exp / "compiled.csv"
        if not csv.exists():
            raise SystemExit(f"missing {csv} — run the per-experiment driver first")
        df = pd.read_csv(csv, low_memory=False)
        df["settings"] = f"{tag}/" + df["settings"].astype(str)
        frames.append(df)
    both = pd.concat(frames, ignore_index=True)

    out_dir = RESULTS / OUT_NAME
    out_dir.mkdir(parents=True, exist_ok=True)
    both.to_csv(out_dir / "compiled.csv", index=False)

    short = comparison_short_all_table(both, target="time_solve", baseline="__none__")
    # sort by geo mean for readability
    short = short.sort_values("Geo. Mean")

    (out_dir / "comparison_short_all.md").write_text(comparison_short_to_markdown(short))
    (out_dir / "comparison_short_all.tex").write_text(comparison_short_to_latex(short))
    print(f"wrote {out_dir}/comparison_short_all.{{md,tex}}  ({len(short)} settings)")


if __name__ == "__main__":
    main()
