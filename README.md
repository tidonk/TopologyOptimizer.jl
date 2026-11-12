# Fast Subproblem Solving in Benders Decomposition for Affine Potential-Based Flow Problems with Topology Switching and Scenarios

This repository contains code that can be used to replicate the results in our paper "Fast Subproblem Solving in Benders Decomposition for Affine Potential-Based Flow Problems with Topology Switching and Scenarios".
While it also contains many other features, the parts relevant for the paper are:

- Benders decomposition of the security-constrained optimal transmission switching problem (SCOTS), to which any affine potential-based flow problem can be reduced, with
  - a subproblem plugin system with different subproblem solvers

    - "classical", i.e. as LP with Gurobi (`FullLPSolver` in repo, `GRB` in paper),
    - as systems of linear equations (`LinEqSysSolver`, `SLE`),
    - as linear equation systems with factorization updates (`LinEqSysRfSolver`, `SLE-U`)

For further explanation of these, we refer to the paper. The repository is organized as follows:

 - [src](src) contains the source code for solving software
 - [check/instances/main](check/instances/main) contains all instances used for the evaluation of the paper
 - [check/scripts](check/scripts) contains scripts used to execute and visualize the experiments (see below)

Note that the contents of this repository were partially developed as part of the project "InnOpTEM" to cover a wider range of features.
While we removed most of what is not directly related to the publication, there might still remain some parts that are not relevant here.

## Installation

Use Julia 1.12.6 and a valid Gurobi license. The pinned dependencies on `paper/subproblem_solving`
include Gurobi 13.0.0. The implementation uses Gurobi's callback interface.
Make sure that `julia --version` reports the intended Julia version; use the full
path to the executable if another version is on your search path.

From the repository root:

```sh
export GRB_LICENSE_FILE=/absolute/path/to/gurobi.lic
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
```

Keep the supplied dependency versions when reproducing experiments. Install and
run with the same Julia depot. The first installation and compilation can take
several minutes. Package precompilation includes the bundled IEEE 118 example,
so the Gurobi license must already be available. No submodules are required.

For the evaluation figures and tables, use Python 3.12, LaTeX, and the
Python packages listed in `check/scripts/visu/requirements.txt`:

```sh
pip install -r check/scripts/visu/requirements.txt
```

## Manual Execution

From the repository root:

```sh
julia --startup-file=no --project=. Main.jl
```

[Main.jl](Main.jl) activates this project and solves the supplied instance
[check/instances/test/pglib_opf_case118_ieee/pglib_opf_case118_ieee_01.m](check/instances/test/pglib_opf_case118_ieee/pglib_opf_case118_ieee_01.m)
with [settings/default.toml](settings/default.toml). This is the modified IEEE 118
example with busbar data. The script prints the solution status, objective, and
number of switching actions. Solver logs are written under `check/results/logs/`.

You can specify the path to the case file, and the "application" that is executed.
The application specifies the problem that is solved and the algorithm that is used to solve it.
Relevant to the publication is:

- `DC_SCOTS_TIME_SLACK_BENDERS` - Benders decomposition

Additionally, a path to a settings file can be provided via the `params_file` argument.
The default settings can be found under [settings/default.toml](settings/default.toml).
This includes for example the subproblem solver that is used.
Custom settings can be added by copying the file and altering it as desired.
We provide a settings file for each of the experiment runs in the [settings](settings) folder.

## Reproducing the Experiments
Experiments are performed on instances synthesized from [pglib-opf](https://github.com/power-grid-lib/pglib-opf).
These can be found under [check/instances/main](check/instances/main).
The experiments in the paper were scheduled using a dispatch script for HTCondor (which likely won't work on your machine) in the command line

```bash
./check/scripts/cluster/dispatch_testset.sh --settings settings/subproblem_solving --experiment_name subproblem_solving
```

where the directory [settings/subproblem_solving](settings/subproblem_solving) contains the settings files corresponding to the different subproblem solvers evaluated in the paper.
Alternatively, you can use the local dispatch script

```bash
./check/scripts/local/dispatch_smallest.sh --n 10 --settings settings/subproblem_solving --experiment_name subproblem_solving
```

to run the n (here 10) smallest instances sequentially on your local machine.
The dispatch scripts use the file [check/scripts/cluster/solve_instance.jl](check/scripts/cluster/solve_instance.jl) to execute the solving process.

We use machines with 2 Intel Xeon L5630 quad-core processors at 2.13 GHz (8 cores, 16 GB RAM) on Debian GNU/Linux 11 across all experiments.
A three-hour time limit is imposed per instance (via the settings files).

## Reproducing Figures & Tables

All scripts used to build the figures and tables in the paper can be found under `check/scripts/visu/`.
With output files placed in the folder `check/results/<experiment>`, execute

```bash
python check/scripts/visu/compare.py <experiment>
```
All output files are placed in the same folder.
The results in the paper match to the files here as follows:

| Paper                            | Files |
| -------------------------------- | --- |
| Table 1                          | `comparison_short.md` |
| Figure 4a                        | `time_decomp_boxplot/detailed_relative.tikz` |
| Figure 4b                        | `time_decomp_boxplot/detailed.tikz` |
| Figure 5                         | `performance_profiles/total/factor.tikz` |
| Figure 6                         | `performance_profiles/SP/factor.tikz` |
| Figure 7                         | `performance_profiles/solve_time/factor.tikz` |
| Figure 8a                        | `performance_profiles/solve_time_scaled/factor.tikz` |
| Figure 8b                        | `performance_profiles/solve_time_scaled/scatter_sorted_by_sle_runtime.tikz` |
| Numbers in the text of Section 5 | `time_decomp/table.tex`, `comparison_extra.tex`, `comparison_short.md`, performance profile files |

Experiment results and historical logs are not included.

## Citation

See [CITATION.cff](CITATION.cff) for the paper citation.

## Contributors

See [CONTRIBUTORS.md](CONTRIBUTORS.md).

## License

See [LICENSE.md](LICENSE.md) and [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
