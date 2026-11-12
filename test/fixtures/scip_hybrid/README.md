The reference is unchanged SCIP 10.0.1 `src/scip/cutsel_hybrid.c`:
https://github.com/scipopt/scip/blob/v10.0.1/src/scip/cutsel_hybrid.c
SHA256: `23e87a1afb2ea568e957c7c71036b025b9dd0f39f63bbdc7d9cacc7ddeeea5ee`.
The upstream Apache license is included here.

`build_reference.py` extracts unchanged `scoring`, `selectBestCut`,
`filterWithParallelism`, and `SCIPselectCutsHybrid` bodies. The shim supplies
row efficacy, absolute parallelism, objective parallelism, integral support,
and caller-provided additive perturbations. Weights are the SCIP defaults:
efficacy 1, objective parallelism .1, integral support .1, directed cutoff 0.
Global-pool membership is false and there is no incumbent. Injected noise is
consumed in upstream order: forced rows first, then optional rows.

Build the reference library from the directory containing this README:

```sh
python3 build_reference.py /tmp/scip_hybrid_reference.so
```

Then run the test from the repository root:

```sh
julia --project=. test/cut_selection_scip_hybrid.jl
```

`SCIP_HYBRID_REFERENCE` overrides the shared-library path. Tests compare exact
ordered optional IDs, including forced filtering, caps (upstream cap zero is
unlimited), signed normals, ties, score components, random perturbations,
and adjacent floating-point score/parallelism boundaries. The Julia tests cover
real JuMP coefficients/objectives and integrality, callback physical values,
positive scaling, forced-count budgets, zero rows and nonempty fallback.
The builder also extracts unchanged `SCIPrandomSetSeed`, `randomGetRand`,
`SCIPrandomGetReal` from official v10.0.1 `misc.c` and the unchanged
`SCIPhashTwo` macro from `pub_misc.h` into a separate RNG oracle library.
Tests compare 10,000 consecutive exact draws at each of eight combinations of
seed shift and multiplier, including the production shift 42, multiplier 10,
plugin seed 0x5eed; they also check forced-first consumption and continuation.
Both upstream files are included unchanged beside the selector source.

This does not exercise SCIP presolve, row construction, global pools, directed
cutoff scoring, or a real solver callback.
