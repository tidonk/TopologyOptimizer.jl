# Third-party notices

Original TopologyOptimizer.jl code is covered by [LICENSE.md](LICENSE.md).
The following upstream terms also apply to the portions identified below.

## PowerModels.jl

This project contains modified code derived from [PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl),
including the following files:

- `src/form/apo.jl`
- `src/form/acr.jl`
- `src/form/wr.jl`
- `src/form/acp.jl`
- `src/form/shared.jl`
- `src/form/dcp.jl`
- `src/io/matpower.jl`
- `src/prob/ots.jl`
- `src/core/base.jl`
- `src/core/constraint_template.jl`
- `src/core/data.jl`
- `src/core/variable.jl`
- `src/core/constraint.jl`

The PowerModels.jl code in these files has been modified for TopologyOptimizer.jl.
The original copyright, conditions, and disclaimer are retained in [LICENSES/PowerModels.txt](LICENSES/PowerModels.txt).
The pinned dependency is PowerModels.jl 0.21.5.
Keep these notices with source distributions and include them in the documentation
or other materials accompanying binary distributions. The upstream names must not
be used to imply endorsement without prior written permission.

## PGLib-OPF data

The bundled input data include modified instances derived from
[PGLib-OPF](https://github.com/power-grid-lib/pglib-opf), maintained by the IEEE PES
Task Force on Benchmarks for Validation of Emerging Power System Algorithms.
The data are provided under [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/).
The library's original license notice is retained in [LICENSES/PGLib-OPF.txt](LICENSES/PGLib-OPF.txt).
The files in `check/instances/` are modified data, rather than an unmodified copy of
the upstream library. Their existing parameter and source information is retained.
[check/instances/SOURCES.tsv](check/instances/SOURCES.tsv) identifies each input
file's source cases. Their original attribution headers, including copyright
notices and references, are retained in
[LICENSES/PGLib-OPF-sources.txt](LICENSES/PGLib-OPF-sources.txt).
See [check/instances/README.md](check/instances/README.md) for the data license
and modification information.

## Gurobi

Gurobi is a separately licensed dependency. Downloading its installer or installing
its Julia package does not provide a solver license. Use is subject to Gurobi's
[End User License Agreement](https://www.gurobi.com/EULA).
Supply a valid license at runtime. For containers, use a supported container license
and mount the license file instead of including it in the image.
