# TopologyOptimizer.jl

TopologyOptimizer.jl provides compact formulations and Benders decomposition for
security-constrained optimal transmission switching (SCOTS) and affine
potential-based flow problems. It builds on
[PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl).
Parts of the software were developed within the InnOpTEM project.

The main features are:

- SCOTS with direct current (DC) power flow models and busbar reconfiguration.
- Subproblem solvers based on linear programming, linear equation systems, and
  factorization updates (`FullLPSolver`, `LinEqSysSolver`, and `LinEqSysRfSolver`).
- Subproblem selection, Benders cut filtering, partial decomposition, and branch fixing.
- Heuristic initialization, parallel solution exchange, and subproblem multithreading.
- Loss estimation and feasibility checks using alternating current power flow models.

## Publications and branches

| Paper | Branch | Tag |
| --- | --- | --- |
| Integrated Topology and Redispatch Optimization of Real-World Power Transmission Systems using Branch-and-Benders-Cut | `paper/integrated` | `integrated-2026-11-10` |
| Benders Cut Filtering for Affine Potential-Based Flow Problems with Robustness Scenarios and Topology Switching | `paper/cut_selection` | `cut_selection-2026-11-10` |
| Fast Subproblem Solving in Benders Decomposition for Affine Potential-Based Flow Problems with Topology Switching and Scenarios | `paper/subproblem_solving` | `subproblem_solving-546be08` |

`paper/integrated` uses the same solver implementation as `main` and adds
experiment settings, reproduction scripts, and exporters for the accompanying
paper. The other paper branches retain their evaluated solver versions.
Each branch has its paper citation in `CITATION.cff`.

## Installation

Use Julia 1.12.6 and a valid Gurobi license. The pinned dependencies on `main`
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

The default configuration is Branch Fixing (row 7 of the accompanying paper's main
evaluation), with switching limit 1 and a solve limit of three hours.

## Settings and other instances

The configuration is in [settings/default.toml](settings/default.toml).
To use different settings, copy this file and pass the new path as `params_file`.
The default values are loaded first, followed by any overrides in that file.
The numbered experiment settings are on `paper/integrated`.

The two main application identifiers are:

| Application | Formulation |
| --- | --- |
| `DC_SCOTS_TIME_SLACK` | Compact DC SCOTS formulation |
| `DC_SCOTS_TIME_SLACK_BENDERS` | Benders decomposition |

For example, from Julia started with `julia --project=.` in the repository root:

```julia
using TopologyOptimizer

result = TopologyOptimizer.main(
    "check/instances/main/pglib_scopf_case85.m";
    params_file = "settings/default.toml",
    busbars = false,
    timesteps = 1,
    slack_types = [:ratea],
    plot = false,
)
display(result["info"])
```

Use `busbars = false` for the paper instances in `check/instances/main/` and
`check/instances/grid/`, as in the reproduction scripts. The IEEE example in
`Main.jl` retains the default `busbars = true`. Busbar reconfiguration itself is
disabled in `settings/default.toml`.

## Instances

`main` contains the union of 132 cases from the testsets of the accompanying paper
and the cut filtering paper in `check/instances/main/`. The paper branches retain
their exact sets: 130 cases for `paper/integrated`, 131 for `paper/cut_selection`,
and 174 for `paper/subproblem_solving`.
Every branch also includes the 31 grid instances in `check/instances/grid/` and
seven IEEE test cases under `check/instances/test/`.

Use the corresponding paper branch or tag for reproduction and follow its
instructions. On `paper/integrated`, `check/scripts/publication/README.md` describes
the main and practical evaluations. On `paper/cut_selection`, use
`check/scripts/reproduce_cutfiltering.py`. The `paper/subproblem_solving` README
covers its local and cluster execution scripts and figure generation.

## Citation

See [CITATION.cff](CITATION.cff) for the paper citation.

## Contributors

See [CONTRIBUTORS.md](CONTRIBUTORS.md).

## License

See [LICENSE.md](LICENSE.md) and [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
