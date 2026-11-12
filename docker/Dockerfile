FROM julia:1.12.6

LABEL maintainer="Tim Donkiewicz <donkiewicz@or.rwth-aachen.de>, Oliver Gaul <gaul@or.rwth-aachen.de>"
LABEL version="0.4"

ARG strFrom1="http:\/\/deb.debian.org\/debian-security"
ARG strTo1="https:\/\/security.debian.org\/debian-security" 
ARG strFrom2="http:\/\/deb.debian.org\/"
ARG strTo2="https:\/\/ftp.halifax.rwth-aachen.de\/"
ARG DEBIAN_FRONTEND=noninteractive

RUN apt update && apt install -y --no-install-recommends ca-certificates
RUN mkdir -p /temp/mv_sourse && cp /etc/apt/sources.list.d/debian.sources /temp/mv_sourse &&  cat /temp/mv_sourse/debian.sources | sed "s/$strFrom1/$strTo1/g" > /etc/apt/sources.list.d/debian.sources  &&  rm -rf /temp/mv_sourse
RUN mkdir -p /temp/mv_sourse && mv /etc/apt/sources.list.d/debian.sources /temp/mv_sourse &&  cat /temp/mv_sourse/debian.sources | sed "s/$strFrom2/$strTo2/g" > /etc/apt/sources.list.d/debian.sources  &&  rm -rf /temp/mv_sourse
RUN apt update && apt upgrade -y &&  apt install -y --no-install-recommends  6tunnel  wget

USER root
# Install Gurobi 13.0.2, as used in the paper experiments.
# Gurobi uses this public download in its own Dockerfiles.
# Use is subject to the Gurobi EULA; supply a valid container license at runtime.
RUN cd /opt &&  wget -q https://packages.gurobi.com/13.0/gurobi13.0.2_linux64.tar.gz &&  tar xzf gurobi13.0.2_linux64.tar.gz &&  ln -s gurobi1302 gurobi
ENV GUROBI_HOME /opt/gurobi/linux64
ENV PATH $PATH:$GUROBI_HOME/bin
ENV LD_LIBRARY_PATH $GUROBI_HOME/lib

# julia instantiate the project and precompile all packages
WORKDIR /julia
COPY Project.toml Manifest.toml LICENSE.md THIRD_PARTY_NOTICES.md ./
COPY LICENSES/ LICENSES/
# Create a minimal stub module so Julia can precompile TopologyOptimizer's dependencies.
# The real source is mounted/checked-out at CI test time.
RUN mkdir -p src && echo 'module TopologyOptimizer; end' > src/TopologyOptimizer.jl
RUN julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
RUN rm -rf src
