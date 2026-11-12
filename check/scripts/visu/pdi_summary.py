"""Write a per-instance PDI (primal-dual integral) summary CSV for an experiment dir.

Walks ``<experiment_dir>/<setting>/logs/`` for each setting subdirectory, parses every
Gurobi master log (``GRB_MP_*.out``), computes the PDI from the primal/dual trajectory
(including all warmstart restarts since :py:func:`parse_gurobi_log` concatenates them),
and emits ``<experiment_dir>/pdi_summary.csv`` with one row per ``(instance, setting)``.

Usage
-----
    python check/scripts/visu/pdi_summary.py <experiment_dir>

Columns: ``instance,setting,pdi,gap,objective,total_time,restarts,n_traj_points``.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

from gurobi_log import compute_pd_integral, parse_gurobi_log


_LOG_RE = "GRB_MP_*.out"


def _instance_from_log_name(log_name: str) -> str:
    """Extract the instance name from ``GRB_MP_<instance>-<timestamp>.out``."""
    stem = Path(log_name).stem
    if stem.startswith("GRB_MP_"):
        stem = stem[len("GRB_MP_"):]
    return stem.rsplit("-", 1)[0]


def collect_pdi(experiment_dir: Path) -> list[dict]:
    rows: list[dict] = []
    for setting_dir in sorted(p for p in experiment_dir.iterdir() if p.is_dir()):
        log_dir = setting_dir / "logs"
        if not log_dir.is_dir():
            continue
        for log_file in sorted(log_dir.glob(_LOG_RE)):
            data = parse_gurobi_log(log_file)
            pdi = compute_pd_integral(data["times"], data["primals"], data["duals"])
            rows.append({
                "instance":      _instance_from_log_name(log_file.name),
                "setting":       setting_dir.name,
                "pdi":           f"{pdi:.6e}",
                "gap":           "" if data["gap"] != data["gap"] else f"{data['gap']:.6e}",
                "objective":     "" if data["objective"] != data["objective"] else f"{data['objective']:.6e}",
                "total_time":    f"{data['total_time']:.4f}",
                "restarts":      len(data["restarts"]),
                "n_traj_points": len(data["times"]),
            })
    return rows


def write_pdi_summary(experiment_dir: Path, out_path: Path | None = None) -> Path:
    rows = collect_pdi(experiment_dir)
    out_path = out_path or (experiment_dir / "pdi_summary.csv")
    fieldnames = ["instance", "setting", "pdi", "gap", "objective",
                  "total_time", "restarts", "n_traj_points"]
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return out_path


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python pdi_summary.py <experiment_dir>", file=sys.stderr)
        sys.exit(1)
    exp_dir = Path(sys.argv[1]).resolve()
    if not exp_dir.is_dir():
        print(f"Not a directory: {exp_dir}", file=sys.stderr)
        sys.exit(1)
    out = write_pdi_summary(exp_dir)
    print(f"Wrote {out} ({sum(1 for _ in out.open()) - 1} rows)")
