FROM julia:1.12.6

LABEL maintainer="Tim Donkiewicz <donkiewicz@or.rwth-aachen.de>, Oliver Gaul <gaul@or.rwth-aachen.de>"
LABEL version="0.5"

RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates python3 \
    && rm -rf /var/lib/apt/lists/*

ENV JULIA_CPU_TARGET=generic \
    JULIA_NUM_THREADS=1 \
    JULIA_PKG_PRECOMPILE_AUTO=0 \
    OPENBLAS_NUM_THREADS=1 \
    GUROBI_JL_USE_GUROBI_JLL=true

# Build from the repository root. The dependency file pins Gurobi_jll 13.0.2.
WORKDIR /opt/topologyoptimizer
COPY Project.toml Manifest.toml LICENSE.md THIRD_PARTY_NOTICES.md ./
COPY LICENSES/ LICENSES/
COPY src/ src/
COPY settings/default.toml settings/default.toml
COPY check/scripts/reproduce_cutfiltering.py check/scripts/reproduce_cutfiltering.py
COPY check/instances/README.md check/instances/SOURCES.tsv check/instances/
COPY check/instances/test/pglib_opf_case118_ieee/pglib_opf_case118_ieee_01.m check/instances/test/pglib_opf_case118_ieee/pglib_opf_case118_ieee_01.m

# The existing precompilation workload solves an example and requires a licence.
# Install packages now, but defer this workload until a licence is mounted.
RUN julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate(; allow_autoprecomp=false)'

CMD ["python3", "check/scripts/reproduce_cutfiltering.py", "--help"]
