#!/usr/bin/env python3
"""Cross-experiment comparison of all ``1-*`` settings.

Combines the four ``1-n_hits*`` experiments into a single comparison. Labels
are derived automatically by parsing each setting's TOML file (merged onto
``settings/default.toml``) and reporting only the parameters that differ
between any two runs. A per-setting TOML dump with those differing keys is
written to ``check/results/1-across/settings_diff/``.

Drops settings that solve fewer than ``MIN_SOLVED`` instances (default 17).
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]              # 1-baseline (this repo)
SHARED = PROJECT_ROOT.parent                      # experiment root
STAGE_REPO = SHARED / "3-SP_selection"            # read-only data source
sys.path.insert(0, str(SCRIPT_DIR.parent))

from reader import read_experiment, cleanup_data  # noqa: E402
from tables import (  # noqa: E402
    TIME_LIMIT,
    comparison_short_all_table,
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
# Zero-touch: read the read-only 3-SP_selection repo, write into the 1-baseline
# staging tree (paper/), never back into the source repo.
RESULTS_DIR = STAGE_REPO / "check/results"
SETTINGS_ROOT = STAGE_REPO / "settings"
SETTINGS_STAGE_DIR = SETTINGS_ROOT / "3-SP_selection"
DEFAULT_TOML = SETTINGS_ROOT / "default.toml"
OUTPUT_DIR = PROJECT_ROOT / "check/results/paper/3-SP_selection/grid/1-n_hits"

# (experiment_subdir, short prefix). Settings live under SETTINGS_STAGE_DIR/<subdir>/
EXPERIMENTS = [
    ("1-n_hits",                   "n+stab"),
    ("1-n_hits_nostab",            "n"),
    ("1-n_hits_pb",                "n+pb"),
    ("1-n_hits_pb_stopafter1feas", "pb+stop"),
    ("1-n_hits_abs_rel",           "n-hyb"),
    ("6-n_since_call",             "nsc"),
]

# Baseline within the merged comparison (Wilcoxon / bolding reference).
BASELINE = "1-n_hits/rnd_rel_100"
TARGET = "time_solve"
MIN_BUSES = 73
MIN_TIME = 0
MIN_SOLVED = 17

# Per-experiment linestyles; per-setting colours assigned automatically.
EXPERIMENT_LINESTYLE = {
    "nh":  "solid",
    "ns":  "dashed",
    "pb":  "dotted",
    "pbs": "dashdot",
    "ar":  (0, (3, 1, 1, 1)),
    "nsc": (0, (1, 1)),
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

# Keys (last path component) to hide from comparison labels: either constant
# across all surviving settings (MODE=RND, T=10, DECAY) or noisy score-weight
# defaults that the reader does not need to see.
LABEL_DROP_LAST_KEYS = {
    "ON",
    "SORTING_MODE",        # always :RND in the kept experiments
    "SCORE_DECAY",         # default 0.9/1.0, not varied within RND grid
    "STOP_AFTER_TIME",     # T=10 / Inf default
    # score component weights from 6-n_since_call
    "load",
    "n_since_call",
    "reliability",
    "total_share",
    "violation_per_cut",
}

# Drop settings whose run produced at least this many CUT_ERROR rows.
CUT_ERROR_DROP_THRESHOLD = 2

# Performance profile readability cap: show only the K best settings
# (by shifted-geo-mean of TARGET, after the cut-error filter, baseline kept).
TOP_K_PROFILE = 5
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
        if last in LABEL_DROP_LAST_KEYS:
            continue
        short = KEY_ABBREV.get(last, last)
        parts.append(f"{short}={_short_value(cfg[k])}")
    for k in sorted(varied_keys):
        if k not in cfg:
            continue
        last = k.split(".")[-1]
        if last in LABEL_DROP_LAST_KEYS:
            continue
        short = KEY_ABBREV.get(last, last)
        parts.append(f"{short}={_short_value(cfg[k])}")
    body = ", ".join(parts) if parts else "--"
    return rf"\texttt{{{prefix}: {body}}}"


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
        df = read_experiment(exp_dir).copy()
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
    df_all = df.copy()  # snapshot pre-cleanup for all-instances summary
    ### TEMPORARY: remove all runs except for three top settings
    # n_before = len(df)
    # # print settings names 
    # print("All settings:")
    # for s in base_settings_all:
    #     print(f"  - {s}")

    # top_settings = [
    #     "n/rnd_abs_2",
    #     "n/rnd_abs_5",
    #     "n/rnd_abs_10",
    # ]
    # print("Keeping:")
    # for s in top_settings:
    #     print(f"  - {s}")

    # df = df[df["settings"].isin(top_settings)].reset_index(drop=True)
    # n_after = len(df)
    # print(f"Filtering to top settings: kept {n_after}/{n_before} rows")
    ###
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

    # All-instances variant: penalize missing / non-OPTIMAL with TIME_LIMIT so
    # every instance contributes to Sum/Mean/Geo. Mean. Uses the pre-cleanup
    # snapshot (df_all) restricted to the included settings.
    df_all_in = df_all[df_all["settings"].isin(include_settings)].reset_index(drop=True)
    short_all = comparison_short_all_table(df_all_in, target=target, baseline=baseline,
                                           label_map=label_map,
                                           include_settings=include_settings)
    (output_dir / "comparison_short_all.md").write_text(comparison_short_to_markdown(short_all))
    (output_dir / "comparison_short_all.tex").write_text(comparison_short_to_latex(short_all))

    # Additional filtered variant: drop any setting with >= CUT_ERROR_DROP_THRESHOLD
    # cut errors so the headline comparison_short.{md,tex} reflects only
    # numerically trustworthy runs. The full unfiltered table is preserved
    # under comparison_short_unfiltered.{md,tex}.
    from collections import Counter
    err_per_setting: Counter = Counter()
    for tup in error_keys:
        err_per_setting[tup[-1]] += 1
    dropped_for_cuterr = [
        s for s in include_settings
        if err_per_setting.get(s, 0) >= CUT_ERROR_DROP_THRESHOLD
    ]
    if dropped_for_cuterr:
        kept = [s for s in include_settings if s not in dropped_for_cuterr]
        kept_baseline = baseline if baseline in kept else (kept[0] if kept else baseline)
        # Move the unfiltered files aside
        for suffix in ("md", "tex"):
            src = output_dir / f"comparison_short.{suffix}"
            if src.exists():
                src.rename(output_dir / f"comparison_short_unfiltered.{suffix}")
            src_all = output_dir / f"comparison_short_all.{suffix}"
            if src_all.exists():
                src_all.rename(output_dir / f"comparison_short_all_unfiltered.{suffix}")
        # Re-emit the headline tables on the kept subset
        short_f = comparison_short_table(
            df[df["settings"].isin(kept)], target=target, baseline=kept_baseline,
            label_map=label_map, include_settings=kept,
        )
        short_all_f = comparison_short_all_table(
            df_all_in[df_all_in["settings"].isin(kept)],
            target=target, baseline=kept_baseline,
            label_map=label_map, include_settings=kept,
        )
        md_remark = (
            "\n\n*Note: "
            f"{len(dropped_for_cuterr)} setting(s) with \u2265 "
            f"{CUT_ERROR_DROP_THRESHOLD} CUT_ERROR rows were dropped from this "
            "table; see `comparison_short_unfiltered.md` for the full grid. "
            "Dropped: "
            + ", ".join(f"`{s}` ({err_per_setting[s]} errs)" for s in dropped_for_cuterr)
            + ".*\n"
        )
        # Emit the drop note as LaTeX comment lines so the .tex file can be
        # safely included via \resizebox{\linewidth}{!}{\input{...}} without
        # trailing prose interfering with the tabular environment.
        _dropped_list = ", ".join(
            f"{s} ({err_per_setting[s]} errs)" for s in dropped_for_cuterr
        )
        tex_remark = (
            f"\n% Note: {len(dropped_for_cuterr)} setting(s) with >= "
            f"{CUT_ERROR_DROP_THRESHOLD} CUT_ERROR rows were dropped; "
            "see comparison_short_unfiltered.tex for the full grid.\n"
            f"% Dropped: {_dropped_list}.\n"
        )
        (output_dir / "comparison_short.md").write_text(
            comparison_short_to_markdown(short_f) + md_remark)
        (output_dir / "comparison_short.tex").write_text(
            comparison_short_to_latex(short_f) + tex_remark)
        (output_dir / "comparison_short_all.md").write_text(
            comparison_short_to_markdown(short_all_f) + md_remark)
        (output_dir / "comparison_short_all.tex").write_text(
            comparison_short_to_latex(short_all_f) + tex_remark)
        print(f"comparison_short: dropped {len(dropped_for_cuterr)} setting(s) "
              f"for >= {CUT_ERROR_DROP_THRESHOLD} CUT_ERROR rows: {dropped_for_cuterr}")


    # Performance profile: too many lines becomes unreadable, so restrict to
    # the TOP_K_PROFILE best settings (lowest shifted geo-mean of TARGET on
    # the cut-error-filtered set, with baseline always included).
    pp_pool = [s for s in include_settings if s not in dropped_for_cuterr] \
        if dropped_for_cuterr else list(include_settings)
    pp_baseline = baseline if baseline in pp_pool else (pp_pool[0] if pp_pool else baseline)
    shift = 10.0
    geo = (
        df[df["settings"].isin(pp_pool)]
        .assign(_t=lambda d: d[target].astype(float) + shift)
        .groupby("settings")["_t"]
        .apply(lambda s: float(np.exp(np.log(s).mean())) - shift)
        .sort_values()
    )
    top_settings = list(geo.head(TOP_K_PROFILE).index)
    if pp_baseline not in top_settings:
        top_settings.append(pp_baseline)
    print(f"performance_profile: top {TOP_K_PROFILE} by shifted geo mean: {top_settings}")
    performance_profile(df, target=target,
                        output=output_dir / f"performance_profile_{target}.pdf",
                        style_map=style_map, label_map=label_map,
                        baseline=pp_baseline, include_settings=top_settings)
    performance_profile(df, target=target, mode="time",
                        output=output_dir / f"performance_profile_{target}_cumulative.pdf",
                        style_map=style_map, label_map=label_map,
                        baseline=pp_baseline, include_settings=top_settings)

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
