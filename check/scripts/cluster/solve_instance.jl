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

include(joinpath(@__DIR__, "runtime_output.jl"))
using Printf

# define the constant for the results directory
additional_consts = Dict{Symbol, Any}(
    :USER_ENVIRONMENT => Dict{Symbol, Any}(
        :DIRECTORIES => Dict{Symbol, Any}(
            :RESULTS_DIR => joinpath(@__DIR__, "..", "..", "..", "check", "results"),
            :LOG_DIR => joinpath(@__DIR__, "..", "..", "..", "check", "results", "logs"),
            :SOLVER_LOG_DIR => joinpath(@__DIR__, "..", "..", "..", "check", "results", "solver_logs")
        )
    )
)

"function to run all optimization runs required for the TopologyOptimizer project"
function single_run(casefile, additional_consts)
    # run full optimization, including topology and RD together
    println("============================================== Running test with LEQ for $(casefile) ==============================================")
    settingsfile = joinpath(@__DIR__, "..", "..", "..", "settings", "test", "LEQ_BB.toml")
    data, result_dict = TopologyOptimizer.main(joinpath(@__DIR__, "..", casefile), additional_consts=additional_consts, params_file = settingsfile, return_data = true)
    settingsname = "$(split(basename(settingsfile), ".")[1])"
    write_runtime_data_to_csv(result_dict, TopologyOptimizer.get_casename(casefile), settingsname, TopologyOptimizer.PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP], results_dir = joinpath(TopologyOptimizer.CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:RESULTS_DIR], settingsname))
end

function multiple_runs(casefile, settingsfiles::Vector{String}, additional_consts; experiment_name="")::Vector{Dict}
    result_dicts = Vector{Dict}(undef, length(settingsfiles))
    for (i, settingsfile) in enumerate(settingsfiles)
        settingsname = "$(split(basename(settingsfile), ".")[1])"

        # set SOLVER_LOG_DIR to results/experiment_name/settingsfile/logs
        additional_consts[:USER_ENVIRONMENT][:DIRECTORIES][:SOLVER_LOG_DIR] = joinpath(@__DIR__, "..", "..", "..", "check", "results", experiment_name, split(basename(settingsfile), ".")[1], "logs")
        println("============================================== Running test for $(casefile) with settings $(settingsfile) ==============================================")
        data, result_dict = TopologyOptimizer.main(joinpath(@__DIR__, "..", casefile), additional_consts=additional_consts, params_file = settingsfile, return_data = true)

        # write results to csv
        result_dict["info"]["settings"] = settingsfile
        result_dict["info"]["casename"] = TopologyOptimizer.get_casename(casefile)
        result_dict["info"]["tal"] = TopologyOptimizer.PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP]
        write_runtime_data_to_csv(result_dict, TopologyOptimizer.get_casename(casefile), settingsname, TopologyOptimizer.PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP], results_dir = joinpath(TopologyOptimizer.CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:RESULTS_DIR], experiment_name, settingsname))
        result_dicts[i] = result_dict
    end
    return result_dicts
end

function print_comparison_table(result_dicts::Vector{Dict})
    # extract relevant data from result_dicts and print comparison table
    println("Comparison of results:")
    @printf("%-30s %-30s %-10s %-10s\n", "Case Name", "Settings File", "Runtime (s)", "Objective Value")

        for (i, result_dict) in enumerate(result_dicts)
            settingsfile = result_dict["info"]["settings"]
            casename = result_dict["info"]["casename"]
            runtime = result_dict["time"]["solve"]
            #@show(result_dict)
            objective_value = result_dict["info"]["objective"]
            #println("$(settingsfile)\t$(runtime)\t$(objective_value)")
            @printf("%-30s %-30s %-10.2f %-10.4f\n", casename, basename(settingsfile), runtime, objective_value)

    end
end

if length(ARGS) >= 2
    # casefile and settingsfile provided: run a single casefile-settingsfile combination
    casefile = ARGS[1]
    settingsfile = ARGS[2]
    experiment_name = length(ARGS) >= 3 ? ARGS[3] : "experiment"
    multiple_runs(casefile, [settingsfile], additional_consts, experiment_name=experiment_name)
elseif length(ARGS) == 1
    # only casefile provided: run single_run with default settings
    casefile = ARGS[1]
    single_run(casefile, additional_consts)
else
    base_dir = joinpath(@__DIR__, "..", "..", "..", "check", "instances")
    result_dictss = Vector{Dict}()
    # read all wrapping folders of instances
    casefolders = readdir(base_dir)
    casefolders = filter(folder -> isdir(joinpath(base_dir, folder)), casefolders)
    casefolders = ["pglib_opf_case118_ieee"]
    for casefolder in casefolders
        for casefile in readdir(joinpath(base_dir, casefolder))
            #single_run(joinpath(base_dir, casefolder, casefile), additional_consts)

            #try
                # all files in settings/
                settingsfiles = filter(file -> endswith(file, ".toml"), readdir(joinpath(@__DIR__, "..", "..", "..", "settings")))
                settingsfiles = joinpath.(joinpath(@__DIR__, "..", "..", "..", "settings"), settingsfiles)
                _result_dicts = multiple_runs(joinpath(base_dir, casefolder, casefile), settingsfiles, additional_consts, experiment_name="experiment")
                append!(result_dictss, _result_dicts)
            #catch e
            #    @error("Error occurred while running test for $(casefile): $e")
            #end
        end
    end
    print_comparison_table(result_dictss)
end