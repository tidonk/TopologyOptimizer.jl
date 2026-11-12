# TopologyOptimizer.jl

This branch provides the code, settings, instances, and reproduction scripts for
"Benders Cut Filtering for Affine Potential-Based Flow Problems with Robustness
Scenarios and Topology Switching". Use branch `paper/cut_selection` or tag
`cut_selection-2026-11-10` for this paper.

TopologyOptimizer.jl provides compact formulations and Benders decomposition for
security-constrained optimal transmission switching (SCOTS) and affine
potential-based flow problems. It builds on
[PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl).
Parts of the software were developed within the InnOpTEM project.

## Installation

Use Julia 1.12.6 and a valid Gurobi license. The pinned dependencies on `paper/cut_selection`
include Gurobi 13.0.2. The implementation uses Gurobi's callback interface.
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

## Run the bundled example

From the repository root:

```sh
julia --startup-file=no --project=. Main.jl
```

[Main.jl](Main.jl) activates this project and solves the supplied instance
[check/instances/test/pglib_opf_case118_ieee/pglib_opf_case118_ieee_01.m](check/instances/test/pglib_opf_case118_ieee/pglib_opf_case118_ieee_01.m)
with [settings/default.toml](settings/default.toml). This is the modified IEEE 118
example with busbar data. The script prints the solution status, objective, and
number of switching actions. Solver logs are written under `check/results/logs/`.

## Reproducing the cut filtering comparison

The 131 instance files are included in `check/instances/main/`. The script requires Python 3.9 or later, with no additional Python packages. From the repository root:

```sh
export GRB_LICENSE_FILE=/absolute/path/to/gurobi.lic
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
python3 check/scripts/reproduce_cutfiltering.py --instances check/instances/main --output demo-results
```

If `julia` on your path is not version 1.12.6, supply `--julia /absolute/path/to/julia-1.12.6` to the script.

The demonstration runs No filtering, SCIP cut selection, Random, Violation, Diversity, and Hybrid on `pglib_scopf_case116`, `pglib_scopf_case218_1`, and `pglib_scopf_case248`. These instances were selected from the completed comparison for short runtimes and the qualitative ordering of the methods. This is an illustrative demonstration, not a random sample or a substitute for the full benchmark. Exact rankings and runtime ratios depend on the machine.

There are 18 sequential solves, each with a 600-second limit. Their recorded solve times sum to approximately 29 minutes on the experimental machines, excluding Julia startup and compilation. Each solve starts a fresh Julia process. The settings use seed 42, four master threads, one subproblem thread, switching limit 1, no subproblem selection, no optimality-cut supplement, and feasibility priority for Violation and Hybrid. Random, Violation, Diversity, and Hybrid use a cut budget of 5% of the scenario count, whereas SCIP uses its default cut limits.

The output directory contains `runtimes.csv`, the individual results and effective settings, `performance_profile.svg`, and `performance_profile.tikz`. The profile uses solver time, excluding Julia startup, and assigns infinite ratios to unsuccessful runs. The additional `performance_profile.tex` can be compiled with pgfplots 1.18 to produce a PDF. Existing run directories are not overwritten. To regenerate only the plots:

```sh
python3 check/scripts/reproduce_cutfiltering.py --output demo-results --plot-only
```

To run the six configurations on all supplied instances with the paper's three-hour limit:

```sh
python3 check/scripts/reproduce_cutfiltering.py --instances check/instances/main --all \
    --time-limit 10800 --output full-results
```

Use only the accompanying 131 files in this directory. This command performs 786 sequential solves and can take several days. It reproduces the main performance profile, not the parameter search. For a shorter check, `--cases pglib_scopf_case116 --methods unfiltered hybrid` selects two solves. Other case and method names can be supplied in the same way.

### Docker

Build from the repository root:

```sh
docker build -f docker/Dockerfile -t topologyoptimizer-cutfiltering .
mkdir -p docker-results
docker run --rm \
    --mount type=bind,src="$(pwd)/check/instances/main",dst=/instances,readonly \
    --mount type=bind,src="$(pwd)/docker-results",dst=/results \
    --mount type=bind,src=/absolute/path/to/gurobi.lic,dst=/run/gurobi.lic,readonly \
    -e GRB_LICENSE_FILE=/run/gurobi.lic \
    topologyoptimizer-cutfiltering python3 check/scripts/reproduce_cutfiltering.py \
    --instances /instances --output /results/demo
```

The licence must support use in the container. The image uses Julia 1.12.6 and the pinned Gurobi 13.0.2 library. It does not contain a university licence configuration, development tests, or experiment logs. Compilation is deferred until the first run, when the licence is available. CI jobs likewise need their own runtime licence configuration.

### Files needed with the paper

Keep the supplied Julia dependency files, `LICENSE.md`, `LICENSES/`,
`THIRD_PARTY_NOTICES.md`, `README.md`, `Main.jl`, `src/`, `settings/default.toml`,
`check/scripts/reproduce_cutfiltering.py`, the bundled IEEE 118 example, and the
131 paper instances together. Include `docker/Dockerfile` and `.dockerignore`
for Docker use. Development tests and cluster submission scripts are not
required for this reproduction script.

## Settings and other instances

The two main application identifiers accepted by `TopologyOptimizer.main` use
direct current (DC) power flow models:

| Application | Formulation |
| --- | --- |
| `DC_SCOTS_TIME_SLACK` | Compact DC SCOTS formulation |
| `DC_SCOTS_TIME_SLACK_BENDERS` | Benders decomposition |

The `params_file` argument selects a settings file. The values in
[settings/default.toml](settings/default.toml) are loaded first, followed by any
overrides in the selected file. Use `busbars = false` for the paper instances in
`check/instances/main/` and the 31 grid instances in `check/instances/grid/`, as in
the reproduction script. The IEEE example in `Main.jl` uses the default
`busbars = true` and includes the corresponding busbar data.

The six paper configurations are defined in
`check/scripts/reproduce_cutfiltering.py`. The script applies their changes
to `settings/default.toml` and writes each run's settings with its results.

Experiment results and historical logs are not included.

## Citation

See [CITATION.cff](CITATION.cff) for the paper citation.

## Contributors

See [CONTRIBUTORS.md](CONTRIBUTORS.md).

## License

See [LICENSE.md](LICENSE.md) and [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).


