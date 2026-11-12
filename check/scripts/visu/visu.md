# Visualization Suite — Developer Reference

This document describes the architecture, call hierarchy, main functions, and key design decisions of the Python visualization suite under `check/scripts/visu/`.

---

## Directory layout

```
check/scripts/visu/
├── gurobi_log.py            # Log parsing, on-disk cache, PD-convergence plots
├── log_metrics.py           # Per-instance metrics derived from Gurobi logs
├── tables.py                # Comparison tables, summary statistics, Wilcoxon test
├── reader.py                # CSV ingestion, experiment compilation
├── selection.py             # Canonical instance-set partitioning
├── perf_profile.py          # Dolan-Moré and cumulative-over-time profiles
├── style.py                 # Colours, LaTeX style, label formatting
├── heatmap.py               # Setting × instance heatmaps
├── time_bar.py              # Status-distribution bars, time-decomposition plots
├── pdi_summary.py           # Standalone PDI summary helper
├── compare.py               # Generic experiment driver (non-seq-warmstart)
└── 6-seq_warmstarting/
    ├── _seqwarmstart_common.py   # Shared runner for seq-warmstart experiments
    ├── compare_1-tal1.py         # Entry point for the 1-tal1 sweep
    └── compare_1-tal1-pd.py      # Entry point for the 1-tal1-pd sweep
```

---

## High-level call hierarchy

```
compare_1-tal1[-pd].py
└── _seqwarmstart_common.run_experiment(experiment_dir, settings_dir, ...)
    ├── discover_settings(experiment_dir)          # lists setting folders
    ├── ordered_grid_settings(discovered)          # canonical display order
    ├── build_style_map(settings) / build_label_map(settings)
    ├── log_metrics.run(experiment_dir, ...)       # ① metrics pipeline
    │   ├── build_log_cache(experiment_dir)
    │   │   └── read_log(log_path)                # per .out file
    │   │       ├── load_parsed_log(csv_path)     # fast path: on-disk CSV
    │   │       └── parse_gurobi_log(filepath)    # slow path: raw .out
    │   │           └── dump_parsed_log(data, csv_path)  # write cache
    │   ├── collect_log_metrics(experiment_dir, raw_cache, logs)
    │   │   ├── _scopf_threshold_per_instance(logs, raw_cache)
    │   │   └── _per_row_metrics(df, logs, raw_cache, scopf_per_inst)
    │   │       └── _compute_metrics(data, scopf_obj)
    │   │           ├── _fair_primal_dual(data, threshold)
    │   │           ├── compute_pd_integral(...)
    │   │           ├── compute_primal_integral(...)
    │   │           ├── compute_time_to_first(...)
    │   │           ├── compute_time_to_first_improving(...)
    │   │           └── compute_time_to_optimal(...)
    │   ├── write_metrics_csv(df, path)
    │   ├── write_metric_tables(out_dir, df, ...)  # per-metric tables
    │   └── write_warmstart_efficiency_table(out_dir, df, ...)
    ├── _scopf_threshold_per_instance(logs, raw_cache)   # reused in plot
    └── tables.run_full_comparison(experiment_dir, ...)  # ② artifacts pipeline
        ├── build_selection(raw_df, ...)
        ├── compile_experiment_csv(experiment_dir, ...)
        ├── status_bar_chart(...)
        ├── comparison_table(selection, target, mode='optimal', ...)
        ├── comparison_table(selection, target, mode='ran', ...)
        ├── performance_profile(df, mode='factor', ...)   # Dolan-Moré
        ├── performance_profile(df, mode='cumulative', ...)
        ├── time_boxplot(df, ...)
        └── plot_aggregated_pd(log_groups, ..., scopf_thresholds)
            ├── _prep_data(setting, files)          # per setting
            │   ├── read_log(f)
            │   ├── [trim warmstart prefix if include_warmstart=False]
            │   └── _anchor_value(d, f, setting)    # leading NaN fill
            ├── _mean_pdi(all_data)                 # sidebar PDI numbers
            └── ax.fill_between(grid, mean_dual, mean_primal, ...)
```

---

## Module descriptions

### `gurobi_log.py` — log parsing and PD plots

**On-disk parsed-log cache**

Each `.out` Gurobi log is parsed once and stored as a sibling CSV at `<setting>/logs_parsed/<stem>.csv`.  The CSV body has columns `t,primal,dual,phase`; a leading `# META {...json...}` comment line holds all scalar fields.

| Function | Purpose |
|---|---|
| `read_log(log_path)` | Unified entry point. Loads from `logs_parsed/<stem>.csv` if present, otherwise calls `parse_gurobi_log` and writes the CSV. |
| `parse_gurobi_log(filepath, main_run_only)` | Parses a raw `.out` file. Returns a dict with fields: `times`, `primals`, `duals`, `restarts`, `optimal`, `objective`, `warmstart_objective`, `scopf_objective`, `gap`, `total_time`, `total_nodes`, `total_iterations`. |
| `dump_parsed_log(data, csv_path)` | Writes parsed data to the on-disk CSV cache. |
| `load_parsed_log(csv_path)` | Reads a cached CSV; raises `ValueError` for stale files missing new fields (triggers transparent re-parse). |

**`parse_gurobi_log` — key fields**

A single `.out` file may contain multiple Gurobi sessions (warmstart sub-MIPs + main solve), separated by restart boundaries. In the seq-warmstart context:

- **Section 1** — SCOPF warmstart sub-MIP (no "Loaded user MIP start" in header)
- **Section 2** — RSCOTS sub-MIP ("Loaded user MIP start with objective = SCOPF_result")
- **Section 3** — Main Benders solve ("Loaded user MIP start with objective = RSCOTS_result")

| Field | Value |
|---|---|
| `scopf_objective` | First "Loaded user MIP start" value = SCOPF result loaded into section 2 |
| `warmstart_objective` | Last "Loaded user MIP start" value = RSCOTS result loaded into section 3 |
| `restarts` | List of wall-clock times at which a new Gurobi session starts |

**Integral and timing functions**

| Function | Definition |
|---|---|
| `compute_pd_integral(times, primals, duals, total_time)` | $\int_0^T (P(t)-D(t))\,dt$ using left-continuous step functions; trajectory extended to `total_time` |
| `compute_primal_integral(times, primals, optimal_obj, total_time)` | $\int_0^T \max(P(t)-f^*, 0)\,dt$ |
| `compute_time_to_first(times, primals)` | Wall time of the first non-NaN primal |
| `compute_time_to_first_improving(times, primals, threshold)` | Wall time of the first primal strictly below `threshold` |
| `compute_time_to_optimal(times, primals, objective)` | Wall time of first primal within tolerance of `objective` |

**Aggregated PD-convergence plot**

`plot_aggregated_pd(log_groups, ..., relative=True, scopf_thresholds=..., top_n=5, baseline='no_warmstart')`

1. Calls `_prep_data(setting, files)` per setting: loads logs via `read_log`, optionally trims the warmstart prefix (`include_warmstart=False` → re-zero time to main solve start), anchors leading NaN primals via `_anchor_value`.
2. Builds `ref_times` from the baseline total times (used as $T_{\text{ref}}$ for relative-time normalization).
3. Pre-computes `pdi_cache[setting] = _mean_pdi(all_data)` for all settings; selects top-N lowest-PDI settings plus baseline for plotting.
4. Interpolates each log to a 500-point grid; computes `mean_primal` and `mean_dual` over the instance ensemble.
5. Plots mean primal as solid line; shades the PD gap with `fill_between(grid, mean_dual, mean_primal)` — **not** clipped at $f^*/f^*=1$, so the full dual lower-bound depth is visible.
6. Sets `ylim(bottom=0.5)` in relative mode.
7. Writes a sidebar table with PDI values normalized to the baseline.

**`_anchor_value` — leading-NaN fill**

For each log, the first grid point may have no incumbent yet (NaN primals). The anchor is:
```
anchor = max(scopf_thresholds[instance], first_real_main_primal)
```
This ensures the per-instance step function never jumps upward once the first real incumbent arrives. For warmstarted runs the anchor usually equals the SCOPF threshold (their first main-solve incumbent is the SCOPF result). For `no_warmstart` the anchor is the actual first B&B incumbent (which is typically worse than SCOPF — it is not artificially credited with SCOPF quality).

---

### `log_metrics.py` — per-instance metrics

**`_scopf_threshold_per_instance(logs, raw_cache)`**

For each instance, takes the minimum `scopf_objective` across all warmstarted-setting logs on that instance. This is the SCOPF warmstart result — the best quality any run received as a starting point. Used as a shared reference across all settings so metrics are directly comparable.

**`_fair_primal_dual(data, threshold)`**

Prepares the trajectory for integral computation. Key transformations:

1. **Anchor**: fills leading NaN primals with `max(threshold, first_real_main_primal)` — same logic as `_anchor_value` in the plot. Does not clip actual B&B incumbents below the SCOPF threshold.
2. **Sub-MIP window**: overrides primals during `t < main_start` with the same anchor, because sub-MIP node-table primals are objectives of a *different* optimization problem (SCOPF/RSCOTS) and must not influence the main-solve integral.
3. **Dual zeroing**: sets duals to 0 during `t < main_start`, because sub-MIP duals are not valid lower bounds for the main MP.

For `no_warmstart` there are no restarts, so steps 2 and 3 are no-ops — only step 1 applies.

**Design decision: consistent anchor everywhere**

The `_fair_primal_dual` anchor (`max(threshold, first_real)`) matches `_anchor_value` in the plot exactly. This means:
- The **visual shaded area** in the aggregated PD plot
- The **sidebar PDI numbers** in the plot legend
- The **standalone PDI table** in `log_metrics/pdi/comparison_short_all.md`
- The **standalone primal-integral table** in `log_metrics/primal_integral/`

…all measure the same quantity on the same fair trajectory.

**`_compute_metrics(data, scopf_obj)`**

Derives all five per-row metrics from a parsed log dict and the per-instance SCOPF threshold:

| Metric | Description |
|---|---|
| `time_to_first` | Wall time of first feasible (any stage) |
| `time_to_first_improving` | Wall time of first primal strictly below the SCOPF threshold |
| `time_to_optimal` | Wall time of first primal matching the final best-known objective |
| `pdi` | $\int_0^T (P_\text{fair}(t) - D_\text{fair}(t))\,dt$ via `_fair_primal_dual` |
| `primal_integral` | $\int_0^T \max(P_\text{fair}(t)-f^*,0)\,dt$ via `_fair_primal_dual` |

**`run(experiment_dir, ...)`**

Top-level entry point for the metrics pipeline. Returns `(df_tables, raw_cache, logs)`. The `raw_cache` and `logs` are passed to `run_full_comparison` so Gurobi logs are never parsed twice.

---

### `tables.py` — comparison tables and artifacts

**`run_full_comparison(experiment_dir, ...)`**

Produces all comparison artifacts in one call. Internally:
1. Reads all settings with `read_experiment(experiment_dir)`.
2. Applies an optional `df_transform` (e.g. `_add_time_with_warmstart` to add warmstart overhead to total time).
3. Builds one `InstanceSelection` shared by all downstream artifacts.
4. Writes two table variants: `mode='optimal'` (only instances where all settings solved to OPTIMAL) and `mode='ran'` (all instances with non-OPTIMAL clipped to `TIME_LIMIT`).
5. Produces performance profile, cumulative profile, time boxplot, and PD-convergence plot.

**`InstanceSelection`** (`selection.py`)

Three nested instance sets used consistently everywhere:

| Set | Meaning |
|---|---|
| `listing` | All instances surviving the experiment-level filters |
| `ran` | Subset where every setting produced a non-broken result |
| `optimal` | Subset where every setting reached OPTIMAL |

`optimal ⊆ ran ⊆ listing`. The `comparison_short.md` table aggregates over `optimal`; `comparison_short_all.md` aggregates over `ran`.

---

### `reader.py` — CSV ingestion

`read_experiment(experiment_dir)` reads every `<setting>/*.csv` under the experiment directory (skipping `log_metrics/`, `old_stats/`, dot-dirs, and `*_infeasible` dirs). The folder name is the authoritative setting label (overrides the CSV `settings` column), so symlinked directories are ingested under their link name.

---

### `_seqwarmstart_common.py` — seq-warmstart experiment runner

`run_experiment(experiment_dir, settings_dir, ...)` is the single entry point for both `compare_1-tal1.py` and `compare_1-tal1-pd.py`. It:
1. Discovers and orders settings canonically (`no_warmstart` → `scopf` → `rscots_nN_w_W` grid → `rscots_full`).
2. Builds style and label maps (colour ramps keyed by pool size `n`, line styles by weight scheme `w`).
3. Runs `log_metrics.run` to produce metric tables, then computes `scopf_thresholds` from the returned cache.
4. Runs `tables.run_full_comparison` with the shared parsed-log cache, passing `scopf_thresholds` through to the PD plot.

Key flags:
- `pd_include_warmstart=False` (default): trims sub-MIP prefix, all curves start at `t=0` of the main solve.
- `pd_top_n=5`: only the 5 lowest-PDI settings plus baseline are shown on the aggregated plot.

---

## PDI vs. primal-integral

Both metrics are computed on the **same** `_fair_primal_dual` trajectory, so they are directly comparable across settings. They differ only in what they measure:

| Metric | Formula | Typical magnitude (1-tal1-pd) |
|---|---|---|
| PDI | $\int_0^T (P(t)-D(t))\,dt$ | ~200–300M (dominated by large initial dual=0 gap) |
| Primal-integral | $\int_0^T \max(P(t)-f^*,0)\,dt$ | ~30–40M |

The shaded band in the aggregated PD plot visualizes the PDI integrand (primal − dual gap) for each setting — the sidebar numbers are the mean PDI, consistent with the standalone `log_metrics/pdi/` tables.

---

## On-disk cache invalidation

Stale parsed-log CSVs (e.g. missing the `scopf_objective` field added in v6.3.260606) are detected in `load_parsed_log` by checking for required META fields. A missing field raises `ValueError`, which `read_log` catches to fall back to a full re-parse, overwriting the CSV. To force a full re-parse of everything, pass `use_csvs_if_possible=False` to `run_experiment`.
