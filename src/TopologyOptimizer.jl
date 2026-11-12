# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                           #
#  Copyright 2026 Operations Research, RWTH Aachen University               #
#                                                                           #
#  Licensed under the Apache License, Version 2.0 (the "License");          #
#  you may not use this file except in compliance with the License.         #
#  You may obtain a copy of the License at                                  #
#                                                                           #
#      http://www.apache.org/licenses/LICENSE-2.0                           #
#                                                                           #
#  Unless required by applicable law or agreed to in writing, software      #
#  distributed under the License is distributed on an "AS IS" BASIS,        #
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. #
#  See the License for the specific language governing permissions and      #
#  limitations under the License.                                           #
#                                                                           #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## Main File ##
## TopologyOptimizer Project ##
module TopologyOptimizer
using PowerModels
using PowerModels: _check_var_keys, _IM

using Dates
using Random
using JuMP
using SparseArrays
using LinearSolve, LinearAlgebra
using MetaGraphsNext
using Gurobi

using MAT

using Serialization
using Requires
using PrecompileTools

# create our module level logger
using Memento
Memento.reset!()

const _LOGGER = Memento.getlogger(@__MODULE__)

Memento.setlevel!.(Memento.getpath(_LOGGER), "emergency")
Memento.setlevel!(_LOGGER, "info")

# Register the module level logger at runtime
function __init__()
    Memento.register(_LOGGER)

    # Clear any stale handlers from precompilation, then add a fresh one with the current stdout
    empty!(_LOGGER.handlers)
    push!(_LOGGER, DefaultHandler(stdout, DefaultFormatter("{level}: {msg}"), Dict{Symbol, Any}(
      :is_colorized => false
    )))

    # Disable propagation to the root logger to prevent duplicate messages
    _LOGGER.propagate = false
    
    # Conditionally load plotting functionality if packages are available
    @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
        @require VegaLite="112f6efa-9a02-5b7d-90c0-432ed331239a" begin
            @require Colors="5ae59095-9a9b-59fe-a467-6f913c188581" begin
                include("core/plt.jl")
            end
        end
    end
end

# Create clock for printing time some steps took
global CLOCK::Float64 = time()
global _CASENAME::String = ""
push!(PowerModels._pm_global_keys, "substations")

# declare variables for TopologyOptimizer parameters and constants
PARAMS, CONSTS = Dict{Symbol, Any}(), Dict{Symbol, Any}()

function set_parameters_constants!(params::Dict{Symbol, Any}, consts::Dict{Symbol, Any})
    global PARAMS = params
    global CONSTS = consts
    global SOLVER_PRECISION_DIGITS = CONSTS[:USER_ENVIRONMENT][:PRECISION][:SOLVER_PRECISION_DIGITS]
    global SOLVER_PRECISION        = 1/10^SOLVER_PRECISION_DIGITS
end

include("core/types.jl")
include("core/data_helpers.jl")
include("core/base.jl")
include("core/graph.jl")
include("core/data.jl")
include("core/variable.jl")
include("core/constraint.jl")
include("core/constraint_template.jl")
include("core/objective.jl")
include("core/benders.jl")
include("core/cutsel.jl")
include("core/cutsel_plot.jl")

include("util/param_reader.jl")
include("util/util.jl")
include("util/warmstart.jl")
include("util/cut_carry_forward.jl")
include("util/parallel_helper_scores.jl")
include("util/admittance_scores.jl")
include("util/parallel_incumbents.jl")
include("util/common.jl")
include("util/heuristics.jl")

include("io/common.jl")
include("io/data.jl")
include("io/slack.jl")
include("io/matpower.jl")

include("ac_feas/common.jl")
include("loss_estimation/common.jl")

include("form/acp.jl")
include("form/acr.jl")
include("form/act.jl")
include("form/apo.jl")
include("form/dcp.jl")
include("form/lpac.jl")
include("form/wr.jl")
include("form/shared.jl")

include("prob/opf.jl")
include("prob/ots.jl")
include("prob/scopf.jl")
include("prob/scots.jl")

include("prob/scots_benders/data.jl")
include("prob/scots_benders/scots_benders.jl")
include("prob/scots_benders/scots_benders_nested.jl")
include("prob/scots_benders/switching_model.jl")
include("prob/scots_benders/gen_model.jl")
include("prob/scots_benders/subproblems/shared.jl")
include("prob/scots_benders/subproblems/matrix.jl")
include("prob/scots_benders/subproblems/ots_benders.jl")
include("prob/scots_benders/subproblems/ots_benders_slack.jl")
include("prob/scots_benders/subproblems/linear_solving.jl")
include("prob/scots_benders/subproblems/solvers/abstract.jl")
include("prob/scots_benders/subproblems/solvers/fulllp.jl")
include("prob/scots_benders/subproblems/solvers/lp.jl")
include("prob/scots_benders/subproblems/solvers/abstract_matrix.jl")
include("prob/scots_benders/subproblems/solvers/abstract_lineqsys.jl")
include("prob/scots_benders/subproblems/solvers/abstract_lineqsysrf.jl")
include("prob/scots_benders/subproblems/solvers/lineqsys.jl")
include("prob/scots_benders/subproblems/solvers/lineqsysrf.jl")
include("prob/scots_benders/subproblems/solvers/lineqsysbb.jl")
include("prob/scots_benders/subproblems/solvers/lineqsysbbrf.jl")
include("prob/scots_benders/subproblems/solvers/fallback.jl")


# declare variables for TopologyOptimizer parameters and constants
PARAMS, CONSTS = Dict{Symbol, Any}(), Dict{Symbol, Any}()

function load_params_constants(params_file::String, constants_file::String)
    default_params_file = joinpath(@__DIR__, "../settings/default.toml")
    global PARAMS = parse_params(default_params_file)
    if params_file != default_params_file
        extend_dict!(PARAMS, parse_params(params_file))
    end
    global CONSTS = parse_params(constants_file)
    global SOLVER_PRECISION_DIGITS = CONSTS[:USER_ENVIRONMENT][:PRECISION][:SOLVER_PRECISION_DIGITS]
    global SOLVER_PRECISION        = 1/10^SOLVER_PRECISION_DIGITS
    Random.seed!(PARAMS[:TECHNICAL][:SEED])
end

"main function to execute code from the TopologyOptimizer project"
function main(data_or_casefile::Union{String, Dict};
              application::String                       = "",                   # application identifier
              slack_types::Vector{Symbol}               = [:ratea],             # types of slack variables out of [:pg, :qg, :ratea, :dclinepf, :vad, :pd]
              busbars::Bool                             = true,                 # use formulation allowing for lossless lines (required to correctly model busbars)
              timesteps::Int64                          = 1,                    # number of time steps
              plot::Bool                                = false,                # generate plots? (requires optional plotting packages)
              params_file::String                       = joinpath(@__DIR__, "../settings/default.toml"), # file containing parameters
              additional_params::Dict{Symbol, Any}      = Dict{Symbol, Any}(),  # additional parameters (if only a few parameters are to be changed)
              additional_consts::Dict{Symbol, Any}      = Dict{Symbol, Any}(),  # additional constants (if only a few constants are to be changed)
              return_data::Bool                         = false                 # return data and result instead of just result?
            )
  
    CLOCK = time()
    t_start = time()
    start_time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    silence()

    # load parameters and constants
    load_params_constants(params_file, joinpath(@__DIR__, "constants.toml"))
    extend_dict!(PARAMS, additional_params)
    extend_dict!(CONSTS, additional_consts)
    apply_legacy_base_pf_in_master!()
    refresh_derived_params!()

    create_dirs(CONSTS[:USER_ENVIRONMENT][:DIRECTORIES])
    casename = get_casename(data_or_casefile)
    global _CASENAME = casename

    if application == "" 
        application = PARAMS[:APPLICATION]::String
    end

    # logging
    #push!(_LOGGER, DefaultHandler(joinpath(CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:LOG_DIR], "log_$casename-$(Dates.format(Dates.now(), "yymmdd-HHMMSS")).log"), DefaultFormatter("[ {level}: {msg}")))
    Memento.info(_LOGGER, "Solving $application for case $casename.\n")

    # full reading and preprocessing of data
    data, nested = read_and_preprocess_data(data_or_casefile, params_file, busbars, application, slack_types, timesteps, start_time)

    # initialize solvers
    solver_convex, solver_nlp, solver_minlp = initialize_solvers(casename, start_time, application=application)

    ### Preprocessing ###
    try
        preprocess_data!(data, application, slack_types, timesteps, nested, solver_convex, solver_nlp)
    catch error
        error isa RestrictedWarmstartFinished || rethrow()
        return return_data ? (data, error.result) : error.result
    end

    ### Solve the application ###
    result = solve_application!(data, application, slack_types, timesteps, nested, solver_convex, solver_nlp)

    ### Postprocess the results ###
    postprocess_results!(data, result, application, timesteps, nested, solver_convex, solver_nlp, t_start)

    ### Plotting (if requested and packages available) ###
    if plot
        if isdefined(TopologyOptimizer, :make_plots)
            return make_plots(data, result, casename)
        else
            @warn "Plotting requested but plotting packages not loaded. Install Plots, VegaLite, and other plotting dependencies and restart Julia."
        end
    end

    if return_data
        return data, result
    else
        return result
    end
end

case = joinpath(@__DIR__, "..", "check", "instances", "test", "pglib_opf_case118_ieee", "pglib_opf_case118_ieee_01.m")
@compile_workload begin
    redirect_stdio(stdout=devnull, stderr=devnull) do
        main(case, application = "DC_SCOTS_TIME_SLACK_BENDERS", plot = false)
    end
end

end
