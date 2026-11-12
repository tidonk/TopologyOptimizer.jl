The reference is SCIP 10.0.1 `src/scip/cutsel_dynamic.c`, downloaded unchanged
from https://github.com/scipopt/scip/blob/v10.0.1/src/scip/cutsel_dynamic.c.
SHA256: `b02079c0fa3526c8cb9f08d7164338bd9de5a3266e3f2e57272de1f87c7c7113`.
The upstream Apache license is included here.

`build_reference.py` extracts the original scoring, projection, best-cut,
parallelism-filter, and selection functions without modifying their bodies.
C functions supply row efficacy/parallelism and SCIP's absolute
1e-9 epsilon comparisons. This tests the selection algorithm, not SCIP's
row construction, presolve, or surrounding solver.

Build the reference library from the directory containing this README:

```sh
python3 build_reference.py /tmp/scip_dynamic_reference.so
```

Then run the test from the repository root:

```sh
julia --project=. test/cut_selection_dynamic.jl
```

`SCIP_DYNAMIC_REFERENCE` can override the shared-library path. The tests cover
both modes, caps, forced rows, signed normals, ties, and epsilon boundaries.
