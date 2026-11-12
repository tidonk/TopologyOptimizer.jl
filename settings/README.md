# Experiment settings

`default.toml` contains the final Branch Fixing configuration (row 7) from the
InnOpTEM paper's main evaluation. The numbered files contain only the changes
needed for each row, in the same order as the table.

| Row | Configuration | Settings file |
| --- | --- | --- |
| 1 | Compact | [1-compact.toml](1-compact.toml) |
| 2 | Benders | [2-benders.toml](2-benders.toml) |
| 3 | Benders SLE-U | [3-benders-sleu.toml](3-benders-sleu.toml) |
| 4a | SP Selection | [4a-sp-selection.toml](4a-sp-selection.toml) |
| 4b | Cut Filtering | [4b-cut-filtering.toml](4b-cut-filtering.toml) |
| 5 | SP+CF combined | [5-sp-cf-combined.toml](5-sp-cf-combined.toml) |
| 6 | Partial Decomp. | [6-partial-decomposition.toml](6-partial-decomposition.toml) |
| 7 | Branch Fixing | [7-branch-fixing.toml](7-branch-fixing.toml) |
| 8a | Heuristic warmstart | [8a-heuristic-warmstart.toml](8a-heuristic-warmstart.toml) |
| 8b | Parallel heuristic | [8b-parallel-heuristic.toml](8b-parallel-heuristic.toml) |
| 9 | Multithreading | [9-multithreading.toml](9-multithreading.toml) |

The main evaluation uses switching limit 1 and a three-hour time limit.
Pass the chosen file as `params_file` to `TopologyOptimizer.main`.
Row 8b requires `check/scripts/cluster/run_parallel_tal.py` to run both processes.
See [the reproduction instructions](../check/scripts/publication/README.md)
for the main and practical evaluations.
