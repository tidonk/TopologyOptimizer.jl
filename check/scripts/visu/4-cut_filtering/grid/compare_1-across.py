#!/usr/bin/env python3
"""Cross-experiment comparison of all cut-filtering settings (Exps. 1-6).

Combines the six ``{1..6}-cutfiltering-*`` experiments into a single comparison.
Labels are derived automatically by parsing each setting's TOML file (merged
onto ``settings/default.toml``) and reporting only the parameters that differ
between any two runs. A per-setting TOML dump with those differing keys is
written to ``check/results/1-across/settings_diff/``.

Drops settings that solve fewer than ``MIN_SOLVED`` instances (default 0).
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[4]
sys.path.insert(0, str(SCRIPT_DIR.parents[1]))

from reader import read_experiment, cleanup_data  # noqa: E402
from tables import (  # noqa: E402
    TIME_LIMIT,
    comparison_short_table,
    comparison_short_to_latex,
    comparison_short_to_markdown,
    comparison_table,
    to_latex,
    to_markdown,
)
from perf_profile import performance_profile  # noqa: E402
from time_bar import (  # noqa: E402
    cuts_boxplot,
    per_instance_time_bar,
    status_bar_chart,
    time_boxplot,
)

# --- Configuration ---------------------------------------------------------
RESULTS_DIR = PROJECT_ROOT / "check/results"
SETTINGS_ROOT = PROJECT_ROOT / "settings"
SETTINGS_STAGE_DIR = SETTINGS_ROOT / "4-cut_filtering"
DEFAULT_TOML = SETTINGS_ROOT / "default.toml"
OUTPUT_DIR = RESULTS_DIR / "1-across"
DEPLOY_DIR = Path.home() / "html" / "test"

# (experiment_subdir, short prefix). Settings live under SETTINGS_STAGE_DIR/<subdir>/
# TEMP (debugging): only rel-rel published to the website. Restore the full list
# below once the remaining sweeps finish on the cluster.
EXPERIMENTS = [
    # ("1-cutfiltering-pb-spsel-rel",  "pb-rel"),
    # ("2-cutfiltering-pb-spsel-abs",  "pb-abs"),
    # ("3-cutfiltering-abs-spsel-abs", "abs-abs"),
    ("4-cutfiltering-rel-spsel-rel", "rel-rel"),
    # ("5-cutfiltering-rel-spsel-abs", "rel-abs"),
    # ("6-cutfiltering-abs-spsel-rel", "abs-rel"),
]

# Baseline within the merged comparison (Wilcoxon / bolding reference).
BASELINE = "4-cutfiltering-rel-spsel-rel/cf100_sp100"
TARGET = "time_solve"
MIN_BUSES = 73
MIN_TIME = 0
MIN_SOLVED = 0

# Per-experiment linestyles; per-setting colours assigned automatically.
EXPERIMENT_LINESTYLE = {
    "pb-rel":  "solid",
    "pb-abs":  "dashed",
    "abs-abs": "dotted",
    "rel-rel": "dashdot",
    "rel-abs": (0, (3, 1, 1, 1)),
    "abs-rel": (0, (1, 1)),
}

# Optional short names for known differing keys (last path component).
KEY_ABBREV = {
    "SORTING_MODE":           "MODE",
    "SCORE_DECAY":            "DECAY",
    "STAB_ROUND_FREQ":        "STAB",
    "ETA_ROBUSTNESS":         "ETA",
    "ETA_ROBUSTNESS_FEAS":    "ETAF",
    "STOP_AFTER_TIME":        "T",
    "STOP_AFTER_N_HITS":      "H",
    "STOP_AFTER_N_HITS_REL":  r"H\%",
    "STOP_AFTER_N_MISSES":    "M",
    "STOP_AFTER_N_FEAS_CUTS": "FC",
    "STOP_ON_PB":             "PB",
    "STOP_ON_PB_REL_BUFFER":  "PBB",
    "MIN_PRIO":               "MINP",
    "ENABLED":                "ON",
}
# ---------------------------------------------------------------------------


# --- Tiny TOML reader/writer (handles the simple syntax used here) ---------
_TOML_NUM = re.compile(r"^-?\d+(\.\d+)?$")


def _strip_inline_comment(s: str) -> str:
    """Remove inline ``#`` comments outside of quoted strings."""
    in_str = False
    quote = ""
    for i, c in enumerate(s):
        if in_str:
            if c == quote and (i == 0 or s[i - 1] != "\\"):
                in_str = False
        elif c in ('"', "'"):
            in_str, quote = True, c
        elif c == "#":
            return s[:i].rstrip()
    return s.rstrip()


def _parse_value(raw: str):
    """Parse a TOML scalar; return as Python value (raw str for inline tables)."""
    raw = raw.strip()
    if not raw:
        return None
    if raw[0] in ('"', "'") and raw[-1] == raw[0]:
        return raw[1:-1]
    low = raw.lower()
    if low == "true":
        return True
    if low == "false":
        return False
    if _TOML_NUM.match(raw):
        return float(raw) if "." in raw else int(raw)
    return raw


def parse_toml_flat(path: Path) -> dict[str, object]:
    """Read a TOML file into a flat {dotted_key: value} dict (simple grammar only)."""
    out: dict[str, object] = {}
    section: list[str] = []
    if not path.is_file():
        return out
    for line in path.read_text().splitlines():
        line = _strip_inline_comment(line).strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            section = [p.strip() for p in line[1:-1].split(".")]
            continue
        if "=" not in line:
            continue
        key, _, raw = line.partition("=")
        key = key.strip()
        if not key:
            continue
        dotted = ".".join(section + [key]) if section else key
        out[dotted] = _parse_value(raw)
    return out


def _format_value(v) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return str(v)
    if isinstance(v, str):
        if v.startswith(("[", "{")) or (v.startswith('"') and v.endswith('"')):
            return v
        return f'"{v}"'
    return str(v)


def write_toml_flat(path: Path, data: dict[str, object]) -> None:
    """Write a flat dotted-key dict back to a sectioned TOML file."""
    sections: dict[tuple, list[tuple[str, object]]] = {}
    for k, v in data.items():
        parts = k.split(".")
        sec = tuple(parts[:-1])
        sections.setdefault(sec, []).append((parts[-1], v))

    lines: list[str] = []
    for sec in sorted(sections):
        if sec:
            lines.append("[" + ".".join(sec) + "]")
        for key, val in sections[sec]:
            lines.append(f"    {key} = {_format_value(val)}")
        lines.append("")
    path.write_text("\n".join(lines).rstrip() + "\n")


# --- Settings diff ---------------------------------------------------------
def collect_full_settings(experiments, results_dir, stage_dir, default_toml):
    """Return (defaults, overrides, full) merged configs for every available setting."""
    defaults = parse_toml_flat(default_toml)
    full: dict[tuple, dict[str, object]] = {}
    overrides: dict[tuple, dict[str, object]] = {}

    for subdir, prefix in experiments:
        exp_dir = results_dir / subdir
        setting_dir_root = stage_dir / subdir
        if not exp_dir.is_dir():
            continue
        setting_names = sorted(d.name for d in exp_dir.iterdir()
                               if d.is_dir() and not d.name.startswith("."))
        for s in setting_names:
            ov = parse_toml_flat(setting_dir_root / f"{s}.toml")
            overrides[(prefix, s)] = ov
            merged = dict(defaults)
            merged.update(ov)
            full[(prefix, s)] = merged
    return defaults, overrides, full


def analyse_experiments(experiments, full, overrides):
    """Per-experiment analysis.

    Returns three dicts keyed by experiment prefix:
      - ``fixed[prefix]``  : keys whose value is constant across all subsettings
                             of that experiment, with that value.
      - ``varied[prefix]`` : keys whose value differs across subsettings of
                             that experiment (the per-subsetting test knobs).
      - ``exp_diff_keys``  : list of fixed keys whose values differ between
                             at least two experiments (= inter-experiment
                             distinguishing parameters).
    """
    candidate_keys = set()
    for ov in overrides.values():
        candidate_keys.update(ov.keys())

    fixed: dict[str, dict[str, object]] = {}
    varied: dict[str, set[str]] = {}
    for _, prefix in experiments:
        cfgs = [cfg for (p, _s), cfg in full.items() if p == prefix]
        if not cfgs:
            continue
        fixed[prefix], varied[prefix] = {}, set()
        for k in candidate_keys:
            vals = {repr(c.get(k)) for c in cfgs}
            if len(vals) == 1:
                fixed[prefix][k] = cfgs[0].get(k)
            else:
                varied[prefix].add(k)

    # Inter-experiment differences (over the fixed-per-experiment portion)
    all_keys = set()
    for d in fixed.values():
        all_keys.update(d.keys())
    exp_diff_keys = []
    for k in sorted(all_keys):
        vals = {repr(fixed[p].get(k)) for p in fixed}
        if len(vals) > 1:
            exp_diff_keys.append(k)
    return fixed, varied, exp_diff_keys


def write_experiment_diffs(experiments, fixed, varied, exp_diff_keys, out_dir):
    """One TOML per experiment with inter-experiment differing keys + varied placeholders."""
    out_dir.mkdir(parents=True, exist_ok=True)
    for _, prefix in experiments:
        if prefix not in fixed:
            continue
        # Group by section
        sections: dict[tuple, list[tuple[str, object, str | None]]] = {}
        # Differing fixed values
        for k in exp_diff_keys:
            if k not in fixed[prefix]:
                continue
            parts = k.split(".")
            sections.setdefault(tuple(parts[:-1]), []).append((parts[-1], fixed[prefix][k], None))
        # Varied-per-subsetting keys: placeholder + comment
        for k in sorted(varied[prefix]):
            parts = k.split(".")
            sections.setdefault(tuple(parts[:-1]), []).append(
                (parts[-1], "<varied>", "varied across subsettings")
            )

        lines: list[str] = [f"# Experiment '{prefix}': parameters differing from"
                            f" the other 1-* experiments."]
        for sec in sorted(sections):
            if sec:
                lines.append("[" + ".".join(sec) + "]")
            for key, val, comment in sections[sec]:
                if comment:
                    lines.append(f"    {key} = {_format_value(val)}  # {comment}")
                else:
                    lines.append(f"    {key} = {_format_value(val)}")
            lines.append("")
        (out_dir / f"{prefix}.toml").write_text("\n".join(lines).rstrip() + "\n")


def _short_value(v) -> str:
    if isinstance(v, bool):
        return "T" if v else "F"
    if isinstance(v, str):
        return v.lstrip(":")
    if isinstance(v, float):
        if v.is_integer():
            return str(int(v))
        return f"{v:g}"
    return str(v)


def build_label(prefix: str, cfg: dict, exp_diff_keys: list[str],
                varied_keys: set[str]) -> str:
    """Label = inter-experiment distinguishing values + this setting's varied values."""
    parts = []
    # Skip keys that are varied within the experiment — those are reported
    # below with the per-subsetting value, avoiding duplicates in the label.
    for k in exp_diff_keys:
        if k not in cfg or k in varied_keys:
            continue
        last = k.split(".")[-1]
        short = KEY_ABBREV.get(last, last)
        parts.append(f"{short}={_short_value(cfg[k])}")
    for k in sorted(varied_keys):
        if k not in cfg:
            continue
        last = k.split(".")[-1]
        short = KEY_ABBREV.get(last, last)
        parts.append(f"{short}={_short_value(cfg[k])}")
    body = ", ".join(parts) if parts else "--"
    return rf"\texttt{{{prefix}: {body}}}"


# --- HTML data-driven page ------------------------------------------------
# The website loads instance-level data from an embedded JSON blob and
# computes per-setting statistics live in the browser. The "Solved by all"
# mode operates over the currently visible (filter-passing) settings; the
# stats are recomputed accordingly. UI state (visible columns, filters,
# sort, mode) is mirrored into the URL hash so views can be shared.

import json

_HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>__TITLE__</title>
<style>
  :root {
    --fg: #1a1a1a; --muted: #6b7280; --border: #e5e7eb;
    --bg: #ffffff; --bg-alt: #f9fafb;
    --accent: #2563eb; --accent-bg: #eff6ff;
    --baseline: #fef3c7; --baseline-border: #f59e0b;
    --mono: ui-monospace, "SF Mono", Menlo, Consolas, "Liberation Mono", monospace;
  }
  * { box-sizing: border-box; }
  body {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
                 Helvetica, Arial, sans-serif;
    margin: 0; padding: 1.5em 2em; color: var(--fg); background: var(--bg);
  }
  h1 { font-size: 1.4em; margin: 0 0 0.25em 0; font-weight: 600; }
  h2 { font-size: 1.05em; margin: 1.5em 0 0.5em 0; font-weight: 600; }
  .meta { color: var(--muted); margin-bottom: 1em; font-size: 0.9em; }
  .meta code { background: var(--bg-alt); padding: 1px 6px; border-radius: 4px;
               font-size: 0.92em; }
  .toolbar {
    display: flex; gap: 1.5em; align-items: center; flex-wrap: wrap;
    padding: 0.7em 0.9em; margin-bottom: 0.9em;
    background: var(--bg-alt); border: 1px solid var(--border);
    border-radius: 8px; font-size: 0.88em;
  }
  .toolbar label { cursor: pointer; user-select: none; }
  .toolbar .group { display: flex; align-items: center; gap: 0.5em; }
  .toolbar button {
    background: var(--bg); border: 1px solid var(--border);
    padding: 4px 10px; border-radius: 5px; cursor: pointer; font-size: inherit;
  }
  .toolbar button:hover { background: var(--accent-bg); border-color: var(--accent); }
  .cols-dropdown { position: relative; }
  .cols-menu {
    display: none; position: absolute; top: 110%; left: 0; z-index: 10;
    background: var(--bg); border: 1px solid var(--border); border-radius: 6px;
    padding: 6px 10px; box-shadow: 0 4px 12px rgba(0,0,0,0.08);
    max-height: 320px; overflow-y: auto; min-width: 180px;
  }
  .cols-menu.open { display: block; }
  .cols-menu label { display: block; padding: 2px 0; white-space: nowrap; }
  .table-wrap { overflow: auto; max-height: 70vh; border: 1px solid var(--border);
                border-radius: 8px; }
  table { border-collapse: separate; border-spacing: 0; font-size: 0.88em;
          width: 100%; background: var(--bg); }
  thead th { background: #f3f4f6; position: sticky; top: 0; z-index: 5;
             border-bottom: 2px solid #9ca3af; font-weight: 600;
             user-select: none; padding: 6px 10px; text-align: right;
             white-space: nowrap; }
  thead th.sortable { cursor: pointer; }
  thead th.sortable:hover { background: #e5e7eb; }
  thead tr.filters th { background: var(--bg-alt); top: 33px;
                        border-bottom: 1px solid var(--border); padding: 3px 5px; }
  thead tr.filters input {
    width: 100%; min-width: 50px; padding: 3px 6px; font-size: 0.9em;
    border: 1px solid var(--border); border-radius: 4px; background: var(--bg);
  }
  thead tr.filters input:focus { outline: none; border-color: var(--accent); }
  th.diff, td.diff { text-align: left; }
  th.exp, td.exp, th.name, td.name { text-align: left; font-family: var(--mono);
                                     font-size: 0.92em; }
  th.exp, th.name { font-family: inherit; }
  th.sort-asc::after  { content: " \25B2"; color: var(--accent); font-size: 0.8em; }
  th.sort-desc::after { content: " \25BC"; color: var(--accent); font-size: 0.8em; }
  tbody td { padding: 4px 10px; border-bottom: 1px solid #f3f4f6;
             white-space: nowrap; }
  tbody tr:nth-child(even) td { background: #fafafa; }
  tbody tr:hover td { background: var(--accent-bg); }
  tbody tr.baseline td { background: var(--baseline);
                         border-top: 1px solid var(--baseline-border);
                         border-bottom: 1px solid var(--baseline-border); }
  td.num { font-variant-numeric: tabular-nums; font-family: var(--mono);
           font-size: 0.92em; text-align: right; }
  .col-hidden { display: none; }
  .legend { margin-top: 0.8em; font-size: 0.82em; color: var(--muted);
            line-height: 1.6; }
  .legend code { background: var(--bg-alt); padding: 1px 5px; border-radius: 3px; }
  .footer { margin-top: 0.6em; font-size: 0.8em; color: var(--muted); }
  .profile { margin-top: 1em; padding: 0.9em 1.2em; background: var(--bg-alt);
             border: 1px solid var(--border); border-radius: 8px; }
  .profile svg { display: block; max-width: 100%; height: auto; }
  .profile .legend-box { display: flex; gap: 1.2em; flex-wrap: wrap;
                         margin-top: 0.5em; font-size: 0.85em;
                         font-family: var(--mono); }
  .profile .legend-box span { display: inline-flex; align-items: center; gap: 0.3em; }
  .profile .swatch { width: 14px; height: 3px; display: inline-block; }
</style>
</head>
<body>
<h1>__TITLE__</h1>
<div class="meta">__META__</div>
<div class="toolbar">
  <div class="group">
    <strong>Stats over:</strong>
    <label><input type="radio" name="mode" value="all" checked> All instances</label>
    <label><input type="radio" name="mode" value="solved"> Solved by all (visible settings)</label>
  </div>
  <div class="group cols-dropdown">
    <button id="cols-btn" type="button">Columns &#9662;</button>
    <div class="cols-menu" id="cols-menu"></div>
  </div>
  <div class="group">
    <button id="clear-filters" type="button">Clear filters</button>
    <button id="share-link" type="button">Copy link</button>
  </div>
  <div class="group footer" id="row-count">&nbsp;</div>
</div>
<div class="table-wrap">
<table id="t">
<thead>
<tr id="head-main"></tr>
<tr class="filters" id="head-filter"></tr>
</thead>
<tbody></tbody>
</table>
</div>
<div class="legend">
Click a column header to sort. Numeric filters accept comparators
(<code>&gt;10</code>, <code>&lt;=100</code>, <code>!=Inf</code>, <code>5..20</code>);
otherwise the filter is a case-insensitive substring match.
Statistics are recomputed for the currently visible rows.
The baseline row is highlighted.
<br>__LEGEND__
</div>

<h2 id="profile-title">Performance profile &mdash; top 5 visible settings</h2>
<div class="profile" id="profile"></div>

<script id="data" type="application/json">__DATA_JSON__</script>
<script>
(function() {
  const DATA = JSON.parse(document.getElementById("data").textContent);
  // DATA shape:
  //   {time_limit, target, baseline, settings:[{id,prefix,name,cfg:{}}, ...],
  //    diff_keys:[{key,short,last}, ...],
  //    instances: [str, ...],      // instance ids (for indexing)
  //    times: [[t_or_null, ...]]   // matrix[setting_idx][instance_idx];
  //                                 // null = missing; ≥time_limit = unsolved
  //   }
  const TL = DATA.time_limit;
  const N_SET = DATA.settings.length;
  const N_INST = DATA.instances.length;
  const TIMES = DATA.times;

  // ---- Column definitions ----------------------------------------------
  const COLS = [
    {label: "Exp", key: "exp", css: "exp", type: "text",
     value: s => s.prefix},
    {label: "Setting", key: "name", css: "name", type: "text",
     value: s => s.name},
  ];
  for (const dk of DATA.diff_keys) {
    COLS.push({label: dk.short, key: "diff:" + dk.key, css: "diff",
               type: "any", title: dk.last,
               value: s => s.cfg[dk.key]});
  }
  COLS.push({label: "Solved", key: "solved", css: "num", type: "num",
             computed: true, fmt: v => fmtInt(v)});
  COLS.push({label: "Sum", key: "sum", css: "num", type: "num",
             computed: true, fmt: v => fmtFixed(v, 1)});
  COLS.push({label: "Mean", key: "mean", css: "num", type: "num",
             computed: true, fmt: v => fmtFixed(v, 2)});
  COLS.push({label: "Geo. Mean", key: "geo", css: "num", type: "num",
             computed: true, fmt: v => fmtFixed(v, 2)});

  // ---- Formatting ------------------------------------------------------
  function fmtInt(v) {
    if (v === null || v === undefined || isNaN(v)) return "—";
    return Number(v).toLocaleString("en-US");
  }
  function fmtFixed(v, digits) {
    if (v === null || v === undefined || isNaN(v)) return "—";
    if (v === Infinity) return "∞";
    return Number(v).toLocaleString("en-US", {
      minimumFractionDigits: digits, maximumFractionDigits: digits
    });
  }
  function fmtCfg(v) {
    if (v === null || v === undefined) return "—";
    if (typeof v === "boolean") return v ? "yes" : "no";
    if (v === "Inf" || v === Infinity) return "∞";
    if (typeof v === "string") return v.replace(/^:/, "");
    return String(v);
  }
  function escapeHtml(s) {
    return String(s).replace(/[&<>"]/g, c => (
      {"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
  }

  // ---- State (with URL hash persistence) -------------------------------
  const state = {
    mode: "all",
    visibleCols: COLS.map(() => true),
    filters: COLS.map(() => ""),
    sort: {col: -1, dir: 1},
  };
  function encodeState() {
    const p = new URLSearchParams();
    p.set("mode", state.mode);
    const hidden = state.visibleCols.map((v,i)=>v?null:i).filter(i=>i!==null);
    if (hidden.length) p.set("hide", hidden.join(","));
    const f = [];
    state.filters.forEach((q,i) => {
      if (q) f.push(i + ":" + encodeURIComponent(q));
    });
    if (f.length) p.set("f", f.join(";"));
    if (state.sort.col >= 0) {
      p.set("sort", state.sort.col + ":" + (state.sort.dir === 1 ? "a" : "d"));
    }
    location.hash = p.toString();
  }
  function decodeState() {
    if (!location.hash) return;
    const p = new URLSearchParams(location.hash.slice(1));
    if (p.has("mode")) state.mode = p.get("mode");
    if (p.has("hide")) {
      const hidden = new Set(p.get("hide").split(",").map(Number));
      state.visibleCols = COLS.map((_,i) => !hidden.has(i));
    }
    if (p.has("f")) {
      for (const part of p.get("f").split(";")) {
        const [i, q] = part.split(":");
        state.filters[+i] = decodeURIComponent(q || "");
      }
    }
    if (p.has("sort")) {
      const [c, d] = p.get("sort").split(":");
      state.sort.col = +c; state.sort.dir = (d === "d" ? -1 : 1);
    }
  }
  decodeState();

  // ---- Filter parsing --------------------------------------------------
  function parseNumExpr(q) {
    q = q.trim(); if (!q) return null;
    const range = q.match(/^(-?\d+(?:\.\d+)?)\s*\.\.\s*(-?\d+(?:\.\d+)?)$/);
    if (range) {
      const lo = parseFloat(range[1]), hi = parseFloat(range[2]);
      return v => !isNaN(v) && v >= lo && v <= hi;
    }
    const m = q.match(/^(>=|<=|>|<|!=|=)?\s*(-?\d+(?:\.\d+)?|inf|-inf)$/i);
    if (!m) return null;
    const op = m[1] || "=";
    const tok = m[2].toLowerCase();
    const x = tok === "inf" ? Infinity : tok === "-inf" ? -Infinity : parseFloat(m[2]);
    return v => {
      if (isNaN(v)) return false;
      switch (op) {
        case ">":  return v >  x;
        case ">=": return v >= x;
        case "<":  return v <  x;
        case "<=": return v <= x;
        case "!=": return v != x;
        default:   return v == x;
      }
    };
  }
  function cellTextValue(setIdx, colIdx, stats) {
    const c = COLS[colIdx];
    if (c.computed) return stats[setIdx][c.key];
    const v = c.value(DATA.settings[setIdx]);
    return v;
  }
  function cellSortKey(setIdx, colIdx, stats) {
    const v = cellTextValue(setIdx, colIdx, stats);
    if (v === null || v === undefined) return -Infinity;
    if (typeof v === "number") return v;
    if (typeof v === "boolean") return v ? 1 : 0;
    if (v === "Inf") return Infinity;
    const f = parseFloat(v);
    return isNaN(f) ? String(v).toLowerCase() : f;
  }
  function matchesFilter(setIdx, colIdx, stats) {
    const q = state.filters[colIdx];
    if (!q) return true;
    const v = cellTextValue(setIdx, colIdx, stats);
    const numFn = parseNumExpr(q);
    if (numFn) {
      const num = (typeof v === "number") ? v
                 : (v === "Inf" || v === Infinity) ? Infinity
                 : parseFloat(v);
      return numFn(num);
    }
    const display = v === null || v === undefined ? "" :
                    typeof v === "boolean" ? (v?"yes":"no") :
                    String(v).replace(/^:/, "");
    return display.toLowerCase().includes(q.toLowerCase());
  }

  // ---- Stat computation ------------------------------------------------
  function computeStats(visibleSetIdx) {
    // Build per-setting stats over the chosen instance set.
    // "all"    : every instance (missing/unsolved = TL)
    // "solved" : only instances solved by ALL currently visible settings
    // For "solved by all", the eligible instance set = those for which every
    // currently visible setting has data AND that data is < TL.
    let instMask = null;
    if (state.mode === "solved" && visibleSetIdx.length) {
      instMask = new Array(N_INST).fill(true);
      for (const si of visibleSetIdx) {
        const row = TIMES[si];
        for (let j = 0; j < N_INST; j++) {
          if (row[j] === null || row[j] >= TL) instMask[j] = false;
        }
      }
    }
    const stats = new Array(N_SET);
    for (let si = 0; si < N_SET; si++) {
      const row = TIMES[si];
      let sum = 0, n = 0, sumLog = 0, nPos = 0, nSolved = 0;
      for (let j = 0; j < N_INST; j++) {
        if (instMask && !instMask[j]) continue;
        let v = row[j];
        // "All" mode: instances without data are treated as TIME_LIMIT
        // (i.e., counted as unsolved with v = TL) so all settings are
        // aggregated over the full instance set.
        // "Solved by all" mode: instMask already excludes nulls.
        if (v === null) {
          if (state.mode === "solved") continue;
          v = TL;
        }
        if (v < TL) nSolved++;
        sum += v; n++;
        if (v > 0) { sumLog += Math.log(v); nPos++; }
      }
      stats[si] = {
        solved: nSolved,
        sum: n ? sum : NaN,
        mean: n ? sum / n : NaN,
        geo: nPos ? Math.exp(sumLog / nPos) : NaN,
      };
    }
    return stats;
  }

  // ---- Render ----------------------------------------------------------
  const headMain = document.getElementById("head-main");
  const headFilter = document.getElementById("head-filter");
  const tbody = document.querySelector("#t tbody");
  const colsMenu = document.getElementById("cols-menu");
  const rowCount = document.getElementById("row-count");

  function buildHeader() {
    headMain.innerHTML = "";
    headFilter.innerHTML = "";
    colsMenu.innerHTML = "";
    COLS.forEach((c, i) => {
      const th = document.createElement("th");
      th.className = c.css + " sortable";
      if (c.title) th.title = c.title;
      th.textContent = c.label;
      if (state.sort.col === i)
        th.classList.add(state.sort.dir === 1 ? "sort-asc" : "sort-desc");
      th.addEventListener("click", () => {
        if (state.sort.col === i) state.sort.dir = -state.sort.dir;
        else { state.sort.col = i; state.sort.dir = 1; }
        encodeState(); render();
      });
      if (!state.visibleCols[i]) th.classList.add("col-hidden");
      headMain.appendChild(th);

      const tf = document.createElement("th");
      tf.className = c.css;
      if (!state.visibleCols[i]) tf.classList.add("col-hidden");
      const inp = document.createElement("input");
      inp.type = "text"; inp.placeholder = "filter…";
      inp.value = state.filters[i]; inp.dataset.col = i;
      inp.addEventListener("input", e => {
        state.filters[i] = e.target.value;
        encodeState(); render();
      });
      tf.appendChild(inp);
      headFilter.appendChild(tf);

      const lbl = document.createElement("label");
      const cb = document.createElement("input");
      cb.type = "checkbox"; cb.checked = state.visibleCols[i];
      cb.dataset.col = i;
      cb.addEventListener("change", () => {
        state.visibleCols[i] = cb.checked;
        encodeState(); render();
      });
      lbl.appendChild(cb);
      lbl.appendChild(document.createTextNode(" " + c.label));
      colsMenu.appendChild(lbl);
    });
  }

  function render() {
    // 1. Compute stats given current mode and current visible-settings set.
    //    Visible settings = those passing all (non-computed) filters.
    //    To handle the chicken-and-egg of computed-col filters, we filter in
    //    two passes: first by non-computed filters, then compute stats, then
    //    apply computed-col filters.
    const passConfig = new Array(N_SET).fill(true);
    for (let si = 0; si < N_SET; si++) {
      for (let ci = 0; ci < COLS.length; ci++) {
        if (COLS[ci].computed) continue;
        if (!matchesFilter(si, ci, null)) { passConfig[si] = false; break; }
      }
    }
    const visibleSetIdx = [];
    for (let si = 0; si < N_SET; si++) if (passConfig[si]) visibleSetIdx.push(si);
    const stats = computeStats(visibleSetIdx);

    // 2. Apply computed-column filters too.
    const finalIdx = visibleSetIdx.filter(si => {
      for (let ci = 0; ci < COLS.length; ci++) {
        if (!COLS[ci].computed) continue;
        if (!matchesFilter(si, ci, stats)) return false;
      }
      return true;
    });

    // 3. Sort.
    if (state.sort.col >= 0) {
      const ci = state.sort.col, dir = state.sort.dir;
      finalIdx.sort((a, b) => {
        const va = cellSortKey(a, ci, stats);
        const vb = cellSortKey(b, ci, stats);
        if (va < vb) return -dir;
        if (va > vb) return  dir;
        return 0;
      });
    }

    // 4. Render rows.
    tbody.innerHTML = "";
    for (const si of finalIdx) {
      const s = DATA.settings[si];
      const tr = document.createElement("tr");
      if (s.id === DATA.baseline) tr.className = "baseline";
      COLS.forEach((c, ci) => {
        const td = document.createElement("td");
        td.className = c.css;
        if (!state.visibleCols[ci]) td.classList.add("col-hidden");
        if (c.computed) {
          td.textContent = c.fmt(stats[si][c.key]);
        } else {
          const v = c.value(s);
          td.textContent = fmtCfg(v);
        }
        tr.appendChild(td);
      });
      tbody.appendChild(tr);
    }
    rowCount.textContent = `${finalIdx.length} / ${N_SET} settings`;

    // 5. Update column-header visibility (in case toggled).
    Array.from(headMain.children).forEach((th, i) =>
      th.classList.toggle("col-hidden", !state.visibleCols[i]));
    Array.from(headFilter.children).forEach((th, i) =>
      th.classList.toggle("col-hidden", !state.visibleCols[i]));

    // 6. Refresh performance profile with the top-5 currently visible settings
    //    in *current sort order*; if no explicit sort is set, fall back to
    //    ascending geometric mean.
    let ranked;
    let sortLabel;
    if (state.sort.col >= 0) {
      ranked = finalIdx.slice(0, 5);
      const dir = state.sort.dir === 1 ? "asc" : "desc";
      sortLabel = `${COLS[state.sort.col].label} ${dir}`;
    } else {
      ranked = finalIdx.slice().sort((a, b) =>
        (stats[a].geo || Infinity) - (stats[b].geo || Infinity)).slice(0, 5);
      sortLabel = "geo asc";
    }
    const title = document.getElementById("profile-title");
    if (title) title.innerHTML =
      `Performance profile &mdash; top 5 visible settings (sorted by ${sortLabel})`;
    drawProfile(ranked, stats);
  }

  // ---- Performance profile (Dolan-Moré) --------------------------------
  const PROFILE_COLORS = ["#000000","#0072B2","#D55E00","#009E73","#CC79A7",
                          "#E69F00","#56B4E9","#999999"];
  function drawProfile(setIdx, stats) {
    const host = document.getElementById("profile");
    if (!setIdx.length) { host.innerHTML = "<em>No visible settings.</em>"; return; }

    // Best time per instance across the chosen settings (treat unsolved as TL).
    const best = new Array(N_INST).fill(Infinity);
    for (const si of setIdx) {
      const row = TIMES[si];
      for (let j = 0; j < N_INST; j++) {
        let v = row[j];
        if (v === null) v = TL;
        if (v < best[j]) best[j] = v;
      }
    }
    // For each setting compute ratios; sort and build CDF on log τ.
    const TAU_MAX = 10.0;  // cap log10 τ
    const curves = setIdx.map(si => {
      const row = TIMES[si];
      const ratios = [];
      for (let j = 0; j < N_INST; j++) {
        if (best[j] === Infinity || best[j] <= 0) continue;
        let v = row[j];
        if (v === null) v = TL;
        const r = v / best[j];
        ratios.push(v >= TL ? Infinity : r);
      }
      ratios.sort((a, b) => a - b);
      const n = ratios.length;
      const pts = [];
      for (let k = 0; k < n; k++) {
        if (ratios[k] === Infinity) break;
        pts.push([ratios[k], (k + 1) / n]);
      }
      return {si, pts, n};
    });

    // SVG drawing
    const W = 720, H = 280, PAD_L = 50, PAD_R = 20, PAD_T = 10, PAD_B = 40;
    const plotW = W - PAD_L - PAD_R, plotH = H - PAD_T - PAD_B;
    const xMax = Math.log10(Math.max(2, Math.min(TAU_MAX, Math.max(
      ...curves.flatMap(c => c.pts.map(p => p[0])), 2))));
    function X(tau) { return PAD_L + (Math.log10(tau) / xMax) * plotW; }
    function Y(p)   { return PAD_T + (1 - p) * plotH; }

    let svg = `<svg viewBox="0 0 ${W} ${H}" xmlns="http://www.w3.org/2000/svg">`;
    // grid + axes
    svg += `<rect x="${PAD_L}" y="${PAD_T}" width="${plotW}" height="${plotH}" fill="white" stroke="#ccc"/>`;
    for (let g = 0; g <= 5; g++) {
      const y = PAD_T + (g / 5) * plotH;
      const p = (1 - g / 5).toFixed(1);
      svg += `<line x1="${PAD_L}" y1="${y}" x2="${PAD_L + plotW}" y2="${y}" stroke="#eee"/>`;
      svg += `<text x="${PAD_L - 6}" y="${y + 4}" text-anchor="end" font-size="10" fill="#666">${p}</text>`;
    }
    const ticks = [1, 2, 5, 10, 20, 50, 100, 1000, 10000];
    for (const t of ticks) {
      const lt = Math.log10(t);
      if (lt > xMax) break;
      const x = X(t);
      svg += `<line x1="${x}" y1="${PAD_T}" x2="${x}" y2="${PAD_T + plotH}" stroke="#eee"/>`;
      svg += `<text x="${x}" y="${PAD_T + plotH + 14}" text-anchor="middle" font-size="10" fill="#666">${t}</text>`;
    }
    svg += `<text x="${PAD_L + plotW/2}" y="${H - 5}" text-anchor="middle" font-size="11" fill="#444">performance ratio τ</text>`;
    svg += `<text x="12" y="${PAD_T + plotH/2}" font-size="11" fill="#444" transform="rotate(-90 12,${PAD_T + plotH/2})" text-anchor="middle">P(rᵢ ≤ τ)</text>`;

    // curves (step function)
    curves.forEach((c, ci) => {
      const col = PROFILE_COLORS[ci % PROFILE_COLORS.length];
      if (!c.pts.length) return;
      let d = `M ${X(1)} ${Y(0)} `;
      let prevY = 0;
      for (const [tau, p] of c.pts) {
        const xv = Math.min(tau, Math.pow(10, xMax));
        d += `L ${X(xv)} ${Y(prevY)} L ${X(xv)} ${Y(p)} `;
        prevY = p;
      }
      d += `L ${PAD_L + plotW} ${Y(prevY)}`;
      svg += `<path d="${d}" fill="none" stroke="${col}" stroke-width="1.8"/>`;
    });
    svg += `</svg>`;

    const legend = curves.map((c, ci) => {
      const col = PROFILE_COLORS[ci % PROFILE_COLORS.length];
      const s = DATA.settings[c.si];
      return `<span><span class="swatch" style="background:${col}"></span>${escapeHtml(s.id)} (mean=${fmtFixed(stats[c.si].mean, 1)})</span>`;
    }).join("");
    host.innerHTML = svg + `<div class="legend-box">${legend}</div>`;
  }

  // ---- Wire up ---------------------------------------------------------
  buildHeader();
  document.querySelectorAll("input[name=mode]").forEach(r => {
    r.checked = (r.value === state.mode);
    r.addEventListener("change", e => {
      state.mode = e.target.value; encodeState(); render();
    });
  });
  document.getElementById("cols-btn").addEventListener("click", e => {
    e.stopPropagation();
    colsMenu.classList.toggle("open");
  });
  document.addEventListener("click", () => colsMenu.classList.remove("open"));
  colsMenu.addEventListener("click", e => e.stopPropagation());
  document.getElementById("clear-filters").addEventListener("click", () => {
    state.filters = COLS.map(() => "");
    headFilter.querySelectorAll("input").forEach(i => i.value = "");
    encodeState(); render();
  });
  document.getElementById("share-link").addEventListener("click", async () => {
    try { await navigator.clipboard.writeText(location.href); }
    catch (e) { window.prompt("Copy link:", location.href); }
  });
  window.addEventListener("hashchange", () => {
    decodeState();
    // Sync header & filter widgets:
    headFilter.querySelectorAll("input").forEach(inp => {
      inp.value = state.filters[+inp.dataset.col] || "";
    });
    colsMenu.querySelectorAll("input[type=checkbox]").forEach(cb => {
      cb.checked = state.visibleCols[+cb.dataset.col];
    });
    document.querySelectorAll("input[name=mode]").forEach(r => {
      r.checked = (r.value === state.mode);
    });
    render();
  });
  render();
})();
</script>
</body>
</html>
"""


def _json_safe(v):
    """Convert TOML-parsed value to JSON-friendly form."""
    if v is None:
        return None
    if isinstance(v, bool):
        return v
    if isinstance(v, (int, float)):
        if isinstance(v, float) and (v != v):
            return None
        return v
    return str(v)


def write_comparison_html(
    df, target, baseline, include_settings, full, exp_diff_keys, varied_per_exp,
    out_path,
):
    """Data-driven HTML page with live stats, URL state, top-5 perf profile."""
    has_tal = "tal" in df.columns
    group_cols = ["instance", "tal"] if has_tal else ["instance"]

    pivoted = df.pivot_table(
        index=group_cols, columns="settings", values=target, aggfunc="first"
    ).reindex(columns=include_settings)
    for s in include_settings:
        mask = df[df["settings"] == s].set_index(group_cols)["status"] == "TIME_LIMIT"
        for idx in mask[mask].index:
            if idx in pivoted.index:
                pivoted.loc[idx, s] = TIME_LIMIT
    vals = pivoted.apply(pd.to_numeric, errors="coerce")

    # Instance identifiers (string keys, stable ordering)
    instance_ids = [
        "|".join(str(p) for p in (idx if isinstance(idx, tuple) else (idx,)))
        for idx in vals.index
    ]

    # times[setting_idx][instance_idx] = float or None
    times = []
    for s in include_settings:
        col = vals[s] if s in vals.columns else pd.Series(index=vals.index, dtype=float)
        col = col.reindex(vals.index)
        row = []
        for v in col:
            if pd.isna(v):
                row.append(None)
            else:
                row.append(float(v))
        times.append(row)

    # Diff key list
    diff_keys = list(exp_diff_keys)
    seen = set(diff_keys)
    for vk in varied_per_exp.values():
        for k in sorted(vk):
            if k not in seen:
                diff_keys.append(k)
                seen.add(k)
    # Reorder: LIMIT_MODE keys first, then their companion LIMIT_* knobs,
    # then everything else (stable within each group).
    _LIMIT_KNOB_ORDER = ("LIMIT_ABS_NUM", "LIMIT_REL_ALPHA",
                         "LIMIT_REL_ITER_ALPHA", "LIMIT_PB_RHO")
    def _diff_key_rank(k):
        last = k.split(".")[-1]
        if last == "LIMIT_MODE":
            return (0, k)
        if last in _LIMIT_KNOB_ORDER:
            return (1, _LIMIT_KNOB_ORDER.index(last), k)
        return (2, k)
    diff_keys = sorted(diff_keys, key=_diff_key_rank)
    diff_keys_info = []
    for k in diff_keys:
        last = k.split(".")[-1]
        short = KEY_ABBREV.get(last, last).replace(r"\%", "%")
        diff_keys_info.append({"key": k, "short": short, "last": last})

    # Map: LIMIT_MODE value -> the companion knob key that is *active*.
    _ACTIVE_KNOB_FOR_MODE = {
        ":abs_num":          "LIMIT_ABS_NUM",
        ":rel_num":          "LIMIT_REL_ALPHA",
        ":rel_num_per_iter": "LIMIT_REL_ITER_ALPHA",
        ":pb":               "LIMIT_PB_RHO",
    }
    _ALL_KNOBS = set(_LIMIT_KNOB_ORDER)

    # Settings metadata
    settings_meta = []
    for s in include_settings:
        prefix, name = s.split("/", 1)
        cfg = full.get((prefix, name), {})
        cfg_out = {k: _json_safe(cfg.get(k)) for k in diff_keys}
        # Blank out inapplicable LIMIT_* knobs based on this setting's LIMIT_MODE.
        # For each LIMIT_MODE key, only the matching companion knob keeps its
        # value; the other knobs sharing the same parent section are blanked.
        for k in diff_keys:
            if k.split(".")[-1] != "LIMIT_MODE":
                continue
            parent = k.rsplit(".", 1)[0]
            mode_val = cfg.get(k)
            active = _ACTIVE_KNOB_FOR_MODE.get(str(mode_val).strip() if mode_val is not None else "")
            for knob in _ALL_KNOBS:
                knob_key = f"{parent}.{knob}" if parent else knob
                if knob_key in cfg_out and knob != active:
                    cfg_out[knob_key] = ""
        settings_meta.append({
            "id": s, "prefix": prefix, "name": name, "cfg": cfg_out,
        })

    data = {
        "time_limit": float(TIME_LIMIT),
        "target": target,
        "baseline": baseline,
        "settings": settings_meta,
        "diff_keys": diff_keys_info,
        "instances": instance_ids,
        "times": times,
    }
    data_json = json.dumps(data, separators=(",", ":"), allow_nan=False)
    # Defang script close-tag inside JSON.
    data_json = data_json.replace("</", "<\\/")

    legend = "; ".join(
        f"<code>{dk['short']}</code> = <code>{dk['key']}</code>"
        for dk in diff_keys_info
    )
    title = f"Across-experiment comparison ({target})"
    meta = (f"baseline: <code>{baseline}</code> &middot; "
            f"{len(instance_ids)} instance groups &middot; "
            f"{len(include_settings)} settings &middot; "
            f"time limit = {int(TIME_LIMIT)} s")

    html = (_HTML_TEMPLATE
            .replace("__TITLE__", title)
            .replace("__META__", meta)
            .replace("__LEGEND__", legend)
            .replace("__DATA_JSON__", data_json))
    out_path.write_text(html, encoding="utf-8")


def deploy_html(src: Path, deploy_dir: Path, dest_name: str | None = None) -> Path | None:
    """Copy HTML output to ``deploy_dir`` (optionally renamed). Returns destination path or None."""
    try:
        deploy_dir.mkdir(parents=True, exist_ok=True)
        dest = deploy_dir / (dest_name or src.name)
        dest.write_text(src.read_text(encoding="utf-8"), encoding="utf-8")
        return dest
    except OSError as e:
        print(f"deploy: failed to copy {src} -> {deploy_dir}: {e}")
        return None
# ---------------------------------------------------------------------------

# --- Data assembly ---------------------------------------------------------
def _combined_name(prefix: str, setting: str) -> str:
    return f"{prefix}/{setting}"


def load_combined(experiments, results_dir):
    frames = []
    for subdir, prefix in experiments:
        exp_dir = results_dir / subdir
        if not exp_dir.is_dir():
            print(f"skip: {exp_dir} not found")
            continue
        try:
            df = read_experiment(exp_dir).copy()
        except FileNotFoundError as e:
            print(f"skip: {exp_dir}: {e}")
            continue
        df["settings"] = df["settings"].apply(lambda s: _combined_name(prefix, s))
        df["experiment"] = prefix
        frames.append(df)
    if not frames:
        raise FileNotFoundError("No experiment data found.")
    return pd.concat(frames, ignore_index=True)


def filter_by_min_solved(df, settings_order, min_solved, target):
    has_tal = "tal" in df.columns
    group_cols = ["instance", "tal"] if has_tal else ["instance"]
    pivoted = df.pivot_table(
        index=group_cols, columns="settings", values=target, aggfunc="first"
    )
    for s in pivoted.columns:
        mask = df[df["settings"] == s].set_index(group_cols)["status"] == "TIME_LIMIT"
        for idx in mask[mask].index:
            if idx in pivoted.index:
                pivoted.loc[idx, s] = TIME_LIMIT

    keep, dropped = [], []
    for s in settings_order:
        if s not in pivoted.columns:
            dropped.append((s, 0))
            continue
        col = pd.to_numeric(pivoted[s], errors="coerce")
        n_solved = int(((col < TIME_LIMIT) & col.notna()).sum())
        if n_solved >= min_solved:
            keep.append(s)
        else:
            dropped.append((s, n_solved))
    if dropped:
        print(f"min_solved={min_solved}: dropping {len(dropped)} settings")
        for s, n in dropped:
            print(f"  - {s}: solved {n}")
    return keep


def assign_colors(settings_order: list[str]) -> dict[str, str]:
    """One distinct colour per base (non-prefixed) setting name."""
    palette = [
        "#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#CC79A7",
        "#D55E00", "#6A3D9A", "#33A02C", "#1F78B4", "#B47800", "#444444",
        "#999999", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6",
    ]
    base_names = []
    for s in settings_order:
        _, base = s.split("/", 1)
        if base not in base_names:
            base_names.append(base)
    return {b: palette[i % len(palette)] for i, b in enumerate(base_names)}


# --- Main pipeline ---------------------------------------------------------
def run_across_comparison(
    experiments=EXPERIMENTS,
    results_dir=RESULTS_DIR,
    stage_dir=SETTINGS_STAGE_DIR,
    default_toml=DEFAULT_TOML,
    output_dir=OUTPUT_DIR,
    deploy_dir=DEPLOY_DIR,
    baseline=BASELINE,
    target=TARGET,
    min_buses=MIN_BUSES,
    min_time=MIN_TIME,
    min_solved=MIN_SOLVED,
):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Settings analysis -------------------------------------------------
    _, overrides, full = collect_full_settings(
        experiments, results_dir, stage_dir, default_toml
    )
    fixed_per_exp, varied_per_exp, exp_diff_keys = analyse_experiments(
        experiments, full, overrides
    )
    print(f"settings: {len(full)} settings across {len(fixed_per_exp)} experiments")
    print(f"inter-experiment differing fixed keys ({len(exp_diff_keys)}):")
    for k in exp_diff_keys:
        print(f"  - {k}")
    for prefix, vk in varied_per_exp.items():
        print(f"varied within '{prefix}' ({len(vk)}): {sorted(vk)}")
    write_experiment_diffs(experiments, fixed_per_exp, varied_per_exp,
                           exp_diff_keys, output_dir / "settings_diff")

    # 2. Combined data + labels/styles -------------------------------------
    df = load_combined(experiments, results_dir)

    if min_buses > 0:
        def _bus_count(inst: str) -> int:
            m = re.search(r"case(\d+)", inst)
            return int(m.group(1)) if m else 0
        n_before = len(df)
        df = df[df["instance"].apply(_bus_count) >= min_buses].reset_index(drop=True)
        print(f"min_buses={min_buses}: kept {len(df)}/{n_before} rows")

    df.to_csv(output_dir / "compiled.csv", index=False)
    print(f"Compiled {len(df)} rows to {output_dir / 'compiled.csv'}")

    # Setting order: per-experiment, sorted by override-file name.
    order = []
    for subdir, prefix in experiments:
        exp_dir = results_dir / subdir
        if not exp_dir.is_dir():
            continue
        for d in sorted(p for p in exp_dir.iterdir() if p.is_dir()):
            name = _combined_name(prefix, d.name)
            if name not in order:
                order.append(name)

    # Use full setting order until after cleanup; then drop low-solve settings.
    if baseline not in order:
        print(f"warning: baseline '{baseline}' not in setting set")
    base_settings_all = list(order)

    # 3. Cleanup + instance filtering (before solve-count filter) ----------
    df, error_keys = cleanup_data(df, exclude_all_timelimit=False, baseline=baseline)
    has_tal = "tal" in df.columns
    group_cols = ["instance", "tal"] if has_tal else ["instance"]

    n_settings_all = len(base_settings_all)
    min_required = min(2, n_settings_all)
    group_counts = df.groupby(group_cols)["settings"].nunique()
    complete_keys = group_counts[group_counts >= min_required].reset_index()[group_cols]
    n_before = len(df)
    df = df.merge(complete_keys, on=group_cols)
    if len(df) < n_before:
        n_dropped = len(group_counts) - len(complete_keys)
        print(f"instance filter: dropped {n_dropped} groups with <{min_required} settings")

    if min_time > 0:
        col_names = ["time_MP", "time_SP_N1"]
        present_cols = [c for c in col_names if c in df.columns]

        def _group_max_time(group):
            base = group[present_cols].sum(axis=1).max() if present_cols else 0
            tgt = group[target].max() if target in group.columns else 0
            return max(base, tgt)

        grp = df.groupby(group_cols)
        keep_mask = grp.apply(_group_max_time) >= min_time
        keep_keys = keep_mask[keep_mask].reset_index()[group_cols]
        n_before_mt = df.groupby(group_cols).ngroups
        df = df.merge(keep_keys, on=group_cols)
        n_after_mt = df.groupby(group_cols).ngroups
        if n_before_mt - n_after_mt > 0:
            print(f"min_time={min_time}: dropped {n_before_mt - n_after_mt} fast groups")

    # 4. Drop low-solve settings AFTER cleanup -----------------------------
    include_settings = filter_by_min_solved(df, base_settings_all, min_solved, target)
    df = df[df["settings"].isin(include_settings)].reset_index(drop=True)
    if baseline not in include_settings:
        print(f"warning: baseline '{baseline}' not present; "
              f"falling back to '{include_settings[0]}'")
        baseline = include_settings[0] if include_settings else baseline

    base_colors = assign_colors(include_settings)
    label_map, style_map = {}, {}
    for s in include_settings:
        prefix, name = s.split("/", 1)
        cfg = full.get((prefix, name), {})
        label_map[s] = build_label(prefix, cfg, exp_diff_keys,
                                   varied_per_exp.get(prefix, set()))
        style_map[s] = {
            "color": base_colors[name],
            "linestyle": EXPERIMENT_LINESTYLE.get(prefix, "solid"),
        }

    # 5. Tables and plots --------------------------------------------------
    status_bar_chart(df, output=output_dir / "status_distribution.pdf",
                     label_map=label_map, baseline=baseline,
                     include_settings=include_settings)

    has_baseline = baseline in include_settings

    table = comparison_table(df, target=target, baseline=baseline,
                             wilcoxon_test=has_baseline, label_map=label_map,
                             error_keys=error_keys,
                             include_settings=include_settings)
    (output_dir / f"comparison_{target}.md").write_text(to_markdown(table))
    (output_dir / f"comparison_{target}.tex").write_text(to_latex(table))

    short = comparison_short_table(df, target=target, baseline=baseline,
                                   label_map=label_map,
                                   include_settings=include_settings)
    (output_dir / "comparison_short.md").write_text(comparison_short_to_markdown(short))
    (output_dir / "comparison_short.tex").write_text(comparison_short_to_latex(short))

    write_comparison_html(
        df, target=target, baseline=baseline,
        include_settings=include_settings, full=full,
        exp_diff_keys=exp_diff_keys, varied_per_exp=varied_per_exp,
        out_path=output_dir / "comparison_short.html",
    )
    print(f"Saved sortable HTML table to {output_dir / 'comparison_short.html'}")
    deployed = deploy_html(output_dir / "comparison_short.html", deploy_dir,
                           dest_name="4-cut_filtering.html")
    if deployed:
        print(f"Deployed HTML to {deployed}")

    performance_profile(df, target=target,
                        output=output_dir / f"performance_profile_{target}.pdf",
                        style_map=style_map, label_map=label_map,
                        baseline=baseline, include_settings=include_settings)
    performance_profile(df, target=target, mode="time",
                        output=output_dir / f"performance_profile_{target}_cumulative.pdf",
                        style_map=style_map, label_map=label_map,
                        baseline=baseline, include_settings=include_settings)

    all_ok = df.groupby(group_cols)["status"].apply(lambda s: (s == "OPTIMAL").all())
    ok_keys = all_ok[all_ok].reset_index()[group_cols]
    df_time = df.merge(ok_keys, on=group_cols)
    per_instance_time_bar(df, output=output_dir / "time_decomposition_per_instance.pdf",
                          label_map=label_map, axes_limits=[300, 3600],
                          baseline=baseline, include_settings=include_settings)
    cuts_boxplot(df_time, output=output_dir / "cuts_boxplot.pdf",
                 label_map=label_map, baseline=baseline,
                 include_settings=include_settings)
    time_boxplot(df_time, output=output_dir / "time_boxplot.pdf",
                 label_map=label_map, baseline=baseline,
                 include_settings=include_settings)

    print(f"\nDone. Wrote outputs to {output_dir}")


if __name__ == "__main__":
    run_across_comparison()
