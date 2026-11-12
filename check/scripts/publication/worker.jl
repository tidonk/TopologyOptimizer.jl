#!/usr/bin/env julia
using TopologyOptimizer, JSON, LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1)
casefile, settings, output = ARGS
mkpath(output)
json_safe(x) = x
json_safe(x::AbstractFloat) = isfinite(x) ? x : nothing
json_safe(x::Enum) = string(x)
json_safe(x::Symbol) = string(x)
json_safe(x::Type) = string(x)
json_safe(x::Function) = string(x)
json_safe(x::AbstractDict) = Dict(string(k) => json_safe(v) for (k,v) in x)
json_safe(x::AbstractArray) = map(json_safe, x)
function savejson(name, value)
    path = joinpath(output, name)
    open(path * ".tmp", "w") do io
        JSON.print(io, json_safe(value), 2)
    end
    mv(path * ".tmp", path; force=true)
end
parts = (Ref{Cint}(), Ref{Cint}(), Ref{Cint}())
TopologyOptimizer.Gurobi.GRBversion(parts...)
savejson("runtime.json", Dict("julia_version" => string(VERSION),
    "gurobi_version" => join(getindex.(parts), "."), "source" => pathof(TopologyOptimizer),
    "julia_threads" => Threads.nthreads(), "depot_paths" => DEPOT_PATH))
consts = Dict{Symbol,Any}(:USER_ENVIRONMENT => Dict(:DIRECTORIES => Dict(
    :RESULTS_DIR => output, :LOG_DIR => output, :SOLVER_LOG_DIR => joinpath(output,"logs"))))
data, result = TopologyOptimizer.main(casefile; params_file=settings, additional_consts=consts,
    plot=false, busbars=false, timesteps=1, slack_types=[:ratea], return_data=true)
for key in ("admittance_warmstart", "rscots_stages", "main_time_limit")
    haskey(data, key) && (result[key] = data[key])
end
if haskey(result, "time")
    result["time"]["main_plus_warmstart"] = get(result["time"], "solve", 0.0) + get(result["time"], "warmstart", 0.0)
end
savejson("resolved-params.json", TopologyOptimizer.PARAMS)
savejson("result.json", result)
