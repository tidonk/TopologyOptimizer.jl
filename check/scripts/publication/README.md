# InnoPTEM paper experiments

This directory runs the main evaluation and the practical evaluation at switching limits 1–5 from one implementation. The numbered main-evaluation files in `settings/` and the practical files in `settings/publication/practical/` contain only overrides of `settings/default.toml`. The default configuration is Branch Fixing (row 7). `settings/publication/origins.tsv` identifies their archived sources and SHA256 checksums. `instances.sha256` identifies the 130 paper instances; case162 is excluded.

## Environment

Use Julia 1.12.6, the pinned Julia dependency file in the repository root, Gurobi 13.0.2 with a valid license, and Python 3.11 or newer. Instantiate the Julia environment before running jobs. Install the Python packages in `requirements.txt` to export figures. Set `JULIA_DEPOT_PATH` if packages are installed outside the default depot. Use an explicit Julia executable when the system version differs.

## Run experiments

Supply the directory containing the 130 `.m` files from the archived experiment data. Inputs are checked against `instances.sha256`. Each case runs in a fresh process and records source hashes, effective parameters, solver versions, results, and logs. Existing output directories are never overwritten. `--dry-run` checks inputs and prints commands without launching a solver.

From the repository root:

```sh
python3.11 check/scripts/publication/run.py --suite main --method sequential \
  --instances check/instances/main --cases all --output check/results/publication-main \
  --julia /path/to/julia

python3.11 check/scripts/publication/run.py --suite practical --method direct --tal 2 \
  --instances check/instances/main --cases pglib_scopf_case116 \
  --output check/results/publication-practical --julia /path/to/julia
```

Main methods are `compact`, `benders`, `benders-sleu`, `spsel`, `cutfilt`, `spselcutfilt`, `partialdec`, `branchfix`, `sequential`, `hybrid`, and `benders-threads`. Practical methods are `direct`, `sequential`, `full`, `warmstart`, and `parallel`; each requires `--tal` from 1 to 5. Repeat the command for the desired methods and switching limits.

`hybrid` and `parallel` use the existing two-process launcher, with two master solver threads per process and one shared deadline. Other configurations use four master threads. The `benders-threads` configuration additionally uses four Benders threads, as in its archived run. The parallel launcher initializes a separate Julia cache for its second process; allow sufficient space and startup time.

## Controls previously selected by source copies

| Control | Main evaluation | Practical evaluation |
|---|---|---|
| `WARMSTARTING.RESTRICTED_ONLY` | `false` | `true` for `direct` and `sequential` |
| `WARMSTARTING.RSCOTS.STAGE_MODE` | `:sequential` | `:single` for `direct`, otherwise `:sequential` |
| `WARMSTARTING.RSCOTS.MIP_GAP` | `0.01` | `0.0` for restricted solves |
| `WARMSTARTING.RSCOTS.BEST_OBJ_STOP_FRACTION` | `0.9` | `0.0` disables the objective target |
| `WARMSTARTING.RSCOTS.HELPER_MIP_GAP` | negative: inherit the main solver tolerance | `0.0` for the parallel helper |
| `WARMSTARTING.RSCOTS.HELPER_DISABLE_OBJECTIVE_STOP` | `false` | `true` for the parallel helper |

Single-stage solves use only the requested switching limit. Sequential solves use stages 1 through that limit and refresh scores between stages. Standalone methods return before the unrestricted solve. The corrected dual recovery is shared by all configurations: a saved primal solution cannot be reused after additional constraints have been added.

Time limits, inclusion of initialization time, candidate selection, and solver tolerances are retained in the corresponding settings. In particular, the exact practical runs at limits 1, 2, and 5 allow six hours after SCOPF, whereas limits 3 and 4 include initialization in that budget. The figure calculations retain the common six-hour cutoff and the original phase clocks. These conventions must not be replaced by one common wall-clock calculation.

## Regenerate the current paper inputs

The reviewed experiment data remain separate from the code. Exporters read the archived publication data without changing it and write manuscript inputs to a supplied output directory. That directory may be a temporary directory with the paper's layout; it need not be the paper checkout.

Required main inputs, under `20261021-final-3h-ablation-heuristics-budget/status-audit/`:

- `snapshot.json`
- `summary-excluding-case162.json`
- `publication-overrides.json`

Also supply `20261022-final-3h-ablation-branchfix-benders-threads/status-audit/summary.json`. Practical inputs are `summary.json` and `trajectories.json` from `20261022-practical-final/`. Absolute result paths stored in those practical records are not needed by the exporters.

```sh
python3.11 check/scripts/publication/export_main.py \
  --input /path/to/main/status-audit \
  --threads-input /path/to/threads/status-audit/summary.json \
  --paper check/results/publication-export/paper \
  --diagnostics check/results/publication-export/main

python3.11 check/scripts/publication/export_figure5.py \
  --input /path/to/20261022-practical-final \
  --paper check/results/publication-export/paper

python3.11 check/scripts/publication/export_figure6.py \
  --input /path/to/20261022-practical-final \
  --paper check/results/publication-export/paper \
  --diagnostics check/results/publication-export/practical
```

The table templates retain the published formatting. Plot calculations and styles are taken from the scripts used for the current paper. The Figure 5 attainment threshold is the absolute tolerance `1e-4`; its historical output filename contains `within-one-percent` but does not define the threshold. Reviewed validity exclusions are applied by the main exporter. Diagnostic summaries are written separately from manuscript inputs.

Fresh solver runs require the same validity review before their results can replace the archived publication data. Source and settings equivalence do not guarantee identical runtimes or identical outcomes for runs stopped by a time limit. No full repeat of the experiment grid is implied by consolidation.

## Consolidation checks

All 36 files were checked against the effective settings obtained from their archived source defaults. The main table, performance profile, and three practical figure files match the current paper inputs apart from generated comments. Six selected solver runs cover the main warmstart, standalone single-stage and sequential methods, the exact baseline, parallel solving, and a two-round dual recovery. Final objectives agree with the archived results within `4.1e-7`. The full grid has not been rerun.
