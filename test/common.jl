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

## Common Methods to Test Algorithmics ##
# is used so that test/*.jl files can be executed individually
if ! isdefined(Main, :TopologyOptimizer)
    include("../src/TopologyOptimizer.jl")
end
using .TopologyOptimizer
using PowerModels
using Test, JuMP
import Memento
import InfrastructureModels
import StatsBase
using DataFrames
using Statistics: mean
using CSV

# load parameters and constants in TopologyOptimizer
parameter_file, constant_file = joinpath(@__DIR__, "..", "settings", "default.toml"), joinpath(@__DIR__, "..", "src", "constants.toml")
TopologyOptimizer.load_params_constants(parameter_file, constant_file)
PARAMS = TopologyOptimizer.parse_params(parameter_file)
CONSTS = TopologyOptimizer.parse_params(constant_file)
SOLVER_PRECISION_DIGITS = CONSTS[:USER_ENVIRONMENT][:PRECISION][:SOLVER_PRECISION_DIGITS]
SOLVER_PRECISION = 1/10^SOLVER_PRECISION_DIGITS
TopologyOptimizer.set_parameters_constants!(PARAMS, CONSTS)

# Suppress warnings during testing.
Memento.setlevel!(Memento.getlogger(InfrastructureModels), "error")
Memento.setlevel!(Memento.getlogger(TopologyOptimizer), "error")
logger_config!("error")

"write testset to file"
function write_testset_to_file(testset::Vector{String}, instance_dir::String)
    open(joinpath(instance_dir, "testset.test"), "w") do f
        for case in testset
            println(f, case)
        end
    end
end

"write objective from test result array"
function write_expected_obj_values(runtime_data, filename)
    # write dataframe (reduced to columns :instance and :objective) with header to file in CSV format
    mkpath(dirname(filename))
    CSV.write(filename, runtime_data[:, [:instance, :objective]], writeheader=true)
end

"convert all columns to Int64, Float64 or String"
function convert_df_types!(df)
    for col in names(df)
        col_data = df[!, col]
        if all(x -> x isa Int || x isa Missing, col_data)
            df[!, col] = Int64.(col_data)
        elseif all(x -> x isa Float64 || x isa Int || x isa Missing, col_data)
            df[!, col] = Float64.(col_data)
        elseif all(x -> x isa JuMP.MOI.TerminationStatusCode, col_data)
            df[!, col] = [string(x) for x in col_data]
        else
            df[!, col] = String.(col_data)
        end
    end
end

"compile data from multiple runs of a single application into a single DataFrame"
function collect_runtime_data(results)::DataFrame
    
    # get column names
    cols_info = keys(results[1]["info"])
    cols_time = keys(results[1]["time"])
    cols_callback = keys(results[1]["callback"])

    # get data per sub-dictionary
    data_info       = [[r["info"][col] for r in results]     for col in cols_info]
    data_time       = [[r["time"][col] for r in results]     for col in cols_time] 
    data_callback   = [[r["callback"][col] for r in results] for col in cols_callback]

    # join data over "instance" column
    data = hcat(data_info..., data_time..., data_callback...)

    # append "time_" to all time columns
    cols_time = ["time_$(col)" for col in cols_time]
    colnames = [cols_info...; cols_time...; cols_callback...]

    # compile data into DataFrame
    df = DataFrame(data, Symbol.(colnames))

    # convert all columns to Int64, Float64 or String
    convert_df_types!(df)
    
    return df
end

"write runtime data to file"
function write_runtime_data(df::DataFrame; name="results", out_dir=CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:RESULTS_DIR])
    # write file
    @info("Writing results to $(joinpath(out_dir, "$(name).res"))")
    mkpath(out_dir)

    # get column names
    cols = names(df)

    # convert all columns that have "time_" as prefix to Float64
    for col in cols
        if startswith("time_", string(col)) && eltype(df[!, col]) != Float64
            df[!, col] = Float64.(df[!, col])
        end
    end

    # add statistics (sum, mean, geo mean) to table for each column that is numeric
    df_full = deepcopy(df)
    push!(df_full, ["Sum"; [eltype(df_full[!, col]) <: Number ? sum(df_full[!, col]) : NaN for col in names(df_full)[2:end]]...], promote=true)
    push!(df_full, ["Mean"; [eltype(df_full[!, col]) <: Number ? mean(df_full[!, col]) : NaN for col in names(df_full)[2:end]]...], promote=true)
    push!(df_full, ["Geo Mean"; [eltype(df_full[!, col]) <: Number ? (isempty([r for r in df_full[!, col] if r >= 0]) ? 0.0 : StatsBase.geomean([r for r in df_full[!, col] if r >= 0])) : NaN for col in names(df_full)[2:end]]...], promote=true)

    # round all Float64 columns to 6 digits
    for col in setdiff(cols, ["instance"])
        if eltype(df_full[!, col]) == Float64
            df_full[!, col] = round.(df_full[!, col], digits=6)
        end
    end

    # write table to file
    open(joinpath(out_dir, "$(name).res"), "w") do f
        println(f, df_full)
    end
    ## write to csv as well
    # remove rows sum, mean and geometric mean
    #df = df[1:end-3, :]
    # add column with index
    df_full.index = 1:nrow(df_full)
    @info("Writing results to $(joinpath(out_dir, "$(name).csv"))")
    CSV.write(joinpath(out_dir, "$(name).csv"), df_full, writeheader=true)
end


"for given cols print comparison tables for all successful applications"
function print_experiment_comparison_tables(runtime_data::Dict{String, DataFrame}, cols::Vector{String}; result_dir::String="results")
    successful_applications = [application for application in applications if !isempty(runtime_data[application])]
    for col in cols
        _compare_df = DataFrame(instance = runtime_data[applications[1]][:, :instance])
        for application in successful_applications
            # leftjoin all successful instances of all applications on the instance column, rename the column to the application name
            _compare_df = leftjoin(_compare_df, rename!(runtime_data[application][:, [:instance, Symbol(col)]], col=>application), on = :instance)
        end
        @info("Comparison of $col on all successful applications:")
        println(_compare_df)
        # write to file
        open(joinpath(result_dir, "comparison_$(col).res"), "w") do f
            println(f, _compare_df)
        end
    end
end

"get directories to set for each test run"
function get_additional_consts(path_prefix, test_name)::Dict{Symbol, Any}
    return Dict(:USER_ENVIRONMENT => Dict(:DIRECTORIES => Dict(:SOLVER_LOG_DIR    => joinpath(path_prefix, CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:RESULTS_DIR], test_name, "logs"),
                                                                :LOG_DIR           => joinpath(path_prefix, CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:RESULTS_DIR], test_name, "logs"),
                                                                :RESULTS_DIR       => joinpath(path_prefix, CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:RESULTS_DIR], test_name))))
end

"get expected results from solution file"
function get_obj_comp(application, testset, path_prefix, validate; sol_dir=CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:SOL_DIR], gen_sw=false)::Tuple{Dict{String, Float64}, String}
    obj_comp = Dict{String, Float64}()
    solu_file = ""
    if validate
        # validate against DC test results
        if startswith(application, "DC")
            solu_file = joinpath(path_prefix, sol_dir, "DC_PIPELINE_TEST$(gen_sw ? "_GEN_SW" : "").solu")
        # validate against AC test results
        elseif startswith(application, "AC")
            if startswith(application, "AC_SCOPF")
                solu_file = joinpath(path_prefix, sol_dir, "AC_SC_PIPELINE_TEST.solu")
            else
                solu_file = joinpath(path_prefix, sol_dir, "AC_PIPELINE_TEST$(gen_sw ? "_GEN_SW" : "").solu")
            end
        end
        
        # read solution file
        if isfile(solu_file)
            expected_obj = CSV.read(solu_file, DataFrame)
            obj_comp = Dict(Pair.(expected_obj.instance, expected_obj.objective))
            @info("Found expected results in $(solu_file). Will check against them.")
        else
            @warn("No expected results found in $(solu_file). Only checking termination status.")
        end
    end
    return obj_comp, solu_file
end

"main function"
function run_multi_test(input_testset::Vector{String};                       # bundled test cases 
              path_prefix::String           = joinpath(@__DIR__, ".."),         # prefix for paths
              application::String           = "DC_SCOTS_TIME_SLACK_BENDERS",    # application to run
              test_name::String             = "DC_SCOTS_BENDERS",                  # name of the test
              additional_params::Dict{Symbol, Any} = Dict{Symbol, Any}(),       # additional parameters to pass to the application
              timesteps::Int64              = 1,
              validate::Bool                = true,
              test_default::Bool            = true,                          # test without busbars and reconfiguration
              test_busbars::Bool            = true,                          # test with busbars     
              test_reconf::Bool             = true,                          # test with busbars and reconfiguration
              test_gen_sw::Bool             = true,                          # test with busbars and reconfiguration and generator switching
            )::DataFrame

    testset = [joinpath(@__DIR__, "..", "check", "instances", "test", get_casename(case), "$(get_casename(case))_01.m") for case in input_testset]

    # get expected results
    obj_comp, solu_file = get_obj_comp(application, testset, path_prefix, validate)
    obj_comp_gen_sw, solu_file_gen_sw = get_obj_comp(application, testset, path_prefix, validate, gen_sw=true)

    # run application for each instance
    allresults = []
    default_results = []
    gen_sw_results = []
    for case in testset
        additional_params[:PRINTING] = Dict(:SOLVER => Dict(:LOG_TO_CONSOLE_MP => false))
        additional_params[:SOLVER] = Dict(:MULTITHREADING => Dict(:N_THREADS_BENDERS => 1)) # ensure multithreading works!

        if test_default
            @testset verbose=true "$(TopologyOptimizer.get_casename(case))" begin
                println()
                println(" ================================ Solving instance $(TopologyOptimizer.get_casename(case)) $(test_name) ================================ ")

                result = TopologyOptimizer.main(case, application = application, plot=false, busbars=false, additional_params=additional_params, additional_consts=get_additional_consts(path_prefix, test_name))
                @test result["termination_status"] == JuMP.MOI.OPTIMAL || result["termination_status"] == JuMP.MOI.LOCALLY_SOLVED
                if haskey(obj_comp, TopologyOptimizer.get_casename(case))
                    @test isapprox(result["objective"], obj_comp[TopologyOptimizer.get_casename(case)], atol=1)
                end
                push!(allresults, result)
                push!(default_results, result)

                println(" ================================================================ ")
            end
        end

        if test_busbars
            @testset verbose=true "$(TopologyOptimizer.get_casename(case)) (with busbars)" begin
                println()
                println(" ================================ Solving instance $(TopologyOptimizer.get_casename(case)) $(test_name)_BB ================================ ")

                result = TopologyOptimizer.main(case, application = application, plot=false, busbars=true, additional_params=additional_params, additional_consts=get_additional_consts(path_prefix, "$(test_name)_BB"))
                @test result["termination_status"] == JuMP.MOI.OPTIMAL
                if haskey(obj_comp, TopologyOptimizer.get_casename(case))
                    @test isapprox(result["objective"], obj_comp[TopologyOptimizer.get_casename(case)], atol=1)
                end
                push!(allresults, result)

                println(" ================================================================ ")
            end
        end

        if test_reconf
            @testset verbose=true "$(TopologyOptimizer.get_casename(case)) (with busbars & reconf)" begin
                println()
                println(" ================================ Solving instance $(TopologyOptimizer.get_casename(case)) $(test_name)_BB_REC ================================ ")

                additional_params[:SUBSTATION_RECONF] = Dict(:BRANCH_RECONF => true, :COUPLER_RECONF => true)
                result = TopologyOptimizer.main(case, application = application, plot=false, busbars=true, additional_params=additional_params, additional_consts=get_additional_consts(path_prefix, "$(test_name)_BB_REC"))
                @test result["termination_status"] == JuMP.MOI.OPTIMAL
                if haskey(obj_comp, TopologyOptimizer.get_casename(case))
                    @test isapprox(result["objective"], obj_comp[TopologyOptimizer.get_casename(case)], atol=1)
                end
                push!(allresults, result)

                println(" ================================================================ ")
            end
        end

        if test_gen_sw
            @testset verbose=true "$(TopologyOptimizer.get_casename(case)) (with generator switching)" begin
                println()
                println(" ================================ Solving instance $(TopologyOptimizer.get_casename(case)) $(test_name)_GEN_SW ================================ ")

                delete!(additional_params, :SUBSTATION_RECONF)
                result = TopologyOptimizer.main(case, application = application, plot=false, busbars=false, additional_params=additional_params, additional_consts=get_additional_consts(path_prefix, "$(test_name)_GEN_SW"))
                @test result["termination_status"] == JuMP.MOI.OPTIMAL || result["termination_status"] == JuMP.MOI.LOCALLY_SOLVED
                if haskey(obj_comp_gen_sw, TopologyOptimizer.get_casename(case))
                    @test isapprox(result["objective"], obj_comp_gen_sw[TopologyOptimizer.get_casename(case)], atol=1)
                end
                push!(allresults, result)
                push!(gen_sw_results, result)

                println(" ================================================================ ")
            end
        end
    end

    # collect runtime data and calculate some statistics on it
    runtime_data = collect_runtime_data(allresults)
    write_runtime_data(runtime_data, name=test_name, out_dir=joinpath(path_prefix, CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:RESULTS_DIR], test_name))

    # write expected results, if file does not exist yet (should not happen, since it compromises correctness of the test)
    if ! isfile(solu_file) && validate && !isempty(default_results)
        write_expected_obj_values(collect_runtime_data(default_results), solu_file)
        @warn("Solution file does not exist yet. Writing objective function values to file $(solu_file).")
    end
    if ! isfile(solu_file_gen_sw) && validate && !isempty(gen_sw_results)
        write_expected_obj_values(collect_runtime_data(gen_sw_results), solu_file_gen_sw)
        @warn("Solution file does not exist yet. Writing objective function values to file $(solu_file_gen_sw).")
    end

    return runtime_data
end

"given a folder and a current date, search for the last index"
function get_last_run_id(dir, date)::Int64
    mkpath(dir)
    # get all files in the directory
    files = readdir(dir)
    # filter files that contain the date
    files = filter(x -> occursin(string(date), x), files)
    # get the last run id
    if isempty(files)
        return 0
    else
        return parse(Int, split(split(files[end], "_")[end], ".")[1])
    end
end

function get_function_name(funct::Function, chopped_prefix::String)
    return chopprefix(string(typeof(funct).parameters[1].name.name), "#" * chopped_prefix)
end