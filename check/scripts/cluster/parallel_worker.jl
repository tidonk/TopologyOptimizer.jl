#!/usr/bin/env julia
# Packages are provisioned before dispatch; never resolve/update in a job.
using TopologyOptimizer, JSON, LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1)

casefile, settings, output, role, limit, threads = ARGS
role in ("donor", "receiver", "control") || error("Invalid worker role")
mkpath(output)
gurobi_version_parts = (Ref{Cint}(), Ref{Cint}(), Ref{Cint}())
TopologyOptimizer.Gurobi.GRBversion(gurobi_version_parts...)
open(joinpath(output, "runtime.json"), "w") do io
    JSON.print(io, Dict("julia_version" => string(VERSION),
                        "julia_executable" => joinpath(Sys.BINDIR, Base.julia_exename()),
                        "gurobi_version" => join(getindex.(gurobi_version_parts), "."),
                        "gurobi_jl_version" => string(Base.pkgversion(TopologyOptimizer.Gurobi)),
                        "julia_threads" => Threads.nthreads(),
                        "master_threads" => parse(Int, threads),
                        "depot_paths" => DEPOT_PATH), 2)
end
helper = get(ENV, "INNOPTEM_EXCHANGE_HELPER", "tal1")
main_tal = parse(Int, get(ENV, "INNOPTEM_EXCHANGE_MAIN_TAL", "2"))
tal = role == "donor" ? (helper == "scopf" ? 0 : 1) : main_tal
params = Dict{Symbol,Any}(
    :LIMITS => Dict(:BRANCH => Dict(:TAL_TIMESTEP => tal)),
    :SOLVER => Dict(:SOLVING_TIME_LIMIT => parse(Float64, limit),
                    :MULTITHREADING => Dict(:N_THREADS_MIP_MP => parse(Int, threads),
                                            :N_THREADS_BENDERS => 1, :N_THREADS_MIP_SP => 1)),
)
if helper != "tal1"
    # Main starts immediately; only the RSCOTS helper initializes itself with SCOPF.
    params[:WARMSTARTING] = Dict(:ENABLED => role == "donor" && helper == "rscots", :WARMSTART_METHOD => :scopf)
    params[:HEURISTIC_SOLVING] = Dict(:ENABLED => false,
        :RSCOTS => Dict(:N_SW_BR_ABS => parse(Int, get(ENV, "INNOPTEM_EXCHANGE_SWITCHABLE", "16")),
            :SWITCHABLE_BRANCHES_LIMIT_MODE => :abs_num))
end
consts = Dict{Symbol,Any}(:USER_ENVIRONMENT => Dict(:DIRECTORIES => Dict(
    :RESULTS_DIR => output, :LOG_DIR => output, :SOLVER_LOG_DIR => joinpath(output, "logs"))))
result = TopologyOptimizer.main(casefile; params_file=settings, additional_params=params,
                      additional_consts=consts, busbars=false, timesteps=1, slack_types=[:ratea], plot=false)
json_safe(x) = x
json_safe(x::AbstractFloat) = isfinite(x) ? x : nothing
json_safe(x::Enum) = string(x)
json_safe(x::Symbol) = string(x)
json_safe(x::Type) = string(x)
json_safe(x::Function) = string(x)
json_safe(x::AbstractDict) = Dict(string(k) => json_safe(v) for (k, v) in x)
json_safe(x::AbstractArray) = map(json_safe, x)
open(joinpath(output, "result.json"), "w") do io
    JSON.print(io, json_safe(result), 2)
end
open(joinpath(output, "resolved-params.txt"), "w") do io
    show(io, MIME("text/plain"), TopologyOptimizer.PARAMS)
end

open(joinpath(output, "resolved-params.json"), "w") do io
    JSON.print(io, json_safe(TopologyOptimizer.PARAMS), 2)
end
