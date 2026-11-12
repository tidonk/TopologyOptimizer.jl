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

using DataFrames, DataStructures, CSV, Pkg
using Base.Threads
import Distributions, Random

# import solvers and JuMP
import Gurobi, Ipopt, Juniper
if haskey(Pkg.project().dependencies, "HSL_jll") 
    import HSL_jll
end
import JuMP
using Printf

export get_casename, get_formulation_from_app, print_solving_statistics, write_decomposition, _add_time_steps, adjust_ratea!
export initialize_solvers

const _mp_non_export_columns = ["BUSBARS", "slack_types", "SUBPROBLEM"]

### General ###
"Convert Any[Any[Dict{String, Any}]] to Vector{Vector{Dict{Symbol, Vector{Int64}}}}"
function apply_types_to_contingenciess(obj::Any)::Vector{Vector{Contingency}}
    # Create an empty vector to hold the final result
    result = Vector{Vector{Dict{Symbol, Vector{Int64}}}}()

    # Loop through each element in the input object
    for element in obj
        # Create an empty vector to hold the current element's result
        inner_result = Vector{Dict{Symbol, Vector{Int64}}}()

        # Loop through each dictionary in the current element
        for dict in element
            # Create a new dictionary with Symbol keys and Vector{Int64} values
            inner_dict = Dict{Symbol, Vector{Int64}}()

            # Loop through each key-value pair in the current dictionary
            for (key, value) in dict
                # Convert the key to a Symbol
                symbol_key = Symbol(key)

                # If the value is a Vector, convert it to Vector{Int64}
                if isa(value, Vector)
                    int_vector = Vector{Int64}(value)
                else
                    int_vector = Vector{Int64}()
                end

                # Add the key-value pair to the inner dictionary
                inner_dict[symbol_key] = int_vector
            end

            # Add the inner dictionary to the inner result vector
            push!(inner_result, inner_dict)
        end

        # Add the inner result vector to the final result vector
        push!(result, inner_result)
    end

    # Return the final result
    return result
end

### Solvers ###
"Initializes all solvers"
function initialize_solvers(casename, start_time; application="OTS", solver_log=PARAMS[:PRINTING][:SOLVER][:LOG_TO_CONSOLE], solver_log_dir=CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:SOLVER_LOG_DIR])
    # make sure log path exists
    mkpath(solver_log_dir)

    # define solvers
    solver_convex = optimizer_with_attributes(() -> Gurobi.Optimizer(CONSTS[:OTHER][:GUROBI_ENV]), 
                    #"OutputFlag" => solver_log, 
                    "LogToConsole" => solver_log, 
                    #"MIPgap" => 0,
                    "Method" => 3,
                    #"Crossover" => 0,
                    #"NodeMethod" => 2,
                    "Threads" => PARAMS[:SOLVER][:MULTITHREADING][:N_THREADS_MIP_MP],
                    "LogFile" => get_log_file_name(solver_log_dir, casename, "unknown", start_time),#joinpath(solver_log_dir, "GRB_$(casename)-$(replace(start_time, " " => "_")).out"),
                    # 15GB 
                    "SoftMemLimit" => 15,
                    "TimeLimit" => PARAMS[:SOLVER][:SOLVING_TIME_LIMIT])

    if haskey(Pkg.project().dependencies, "HSL_jll")
        solver_nlp = optimizer_with_attributes(Ipopt.Optimizer, 
            "hsllib" => HSL_jll.libhsl_path, 
            "print_level"=> 3,
            "linear_solver" => "ma57",
            )
    else
        solver_nlp = optimizer_with_attributes(Ipopt.Optimizer, 
            "tol"=>1e-4, 
            "print_level"=> !solver_log ? 1 : 3,
            "linear_solver" => Sys.isunix==true ? "ma27" : "mumps"
            )
    end

    if haskey(Pkg.project().dependencies, "HSL_jll")
        solver_minlp = JuMP.optimizer_with_attributes(Juniper.Optimizer, 
        "nl_solver" => JuMP.optimizer_with_attributes(Ipopt.Optimizer, 
            "hsllib" => HSL_jll.libhsl_path, 
            "print_level"=> 3,
            "linear_solver" => "ma57"
            ), 
        "log_levels"=>[])
    else
        solver_minlp = JuMP.optimizer_with_attributes(Juniper.Optimizer, 
        "nl_solver" => JuMP.optimizer_with_attributes(Ipopt.Optimizer, 
            "tol"=>1e-4, 
            "print_level"=> !solver_log ? 1 : 3,
            "linear_solver" => Sys.isunix==true ? "ma27" : "mumps"
            ), 
        "log_levels"=>[])
    end    

    return (solver_convex, solver_nlp, solver_minlp)
end

"Sets variables[i]'s attribute[i] to values[i]."
function set_variables_attributes(pm, attributes::Vector, variables::Vector{Symbol}, values::Vector)
    for (n, nw) in nws(pm)
        refs = [var(pm, n)[variable][1] for variable in variables]
        JuMP.MOI.set.([pm.model for attribute in attributes], attributes, refs, values)
    end
    return pm
end

### String Extraction ###
"Converts a casefile (including dir) to the casename"
function get_casename(casefile::String)::String
    return chopsuffix(chopsuffix(basename(casefile), ".mat"),".m")
end
function get_casename(data::Dict)::String
    return data["name"]
end

"Gets formulation from AC_SCOTS application"
function get_formulation_from_app(app)
    return CONSTS[:INTERNAL_MAPPINGS][:MODEL_SHORTNAMES_MAP][Symbol(split(app, "AC_SCOTS_")[2])]
end

### Printing ###
"Prints brief solving statistics"
function print_solving_statistics(casename, result; application="unknown")
    tal_timestep, tal_day, tal_component_and_day = PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP], PARAMS[:LIMITS][:BRANCH][:TAL_DAY], PARAMS[:LIMITS][:BRANCH][:TAL_COMPONENT_AND_DAY]
    limits = "TAL_TIMESTEP = $tal_timestep, TAL_DAY = $tal_day, TAL_COMPONENT_AND_DAY = $tal_component_and_day"
    Memento.info(_LOGGER, "Solving statistics:
 Network:                 $casename
 Application:             $application $(occursin("SCOTS", application) ? limits : "")
 Termination status:      $(result["termination_status"])
 Solving time:            $(round(result["solve_time"], digits=2)) seconds
 $("measured_runtime" in keys(result) ? "Measured time in solver: $(round(result["measured_runtime"], digits=2)) seconds" : "")
 Objective function:      $(result["objective"])
 Gurobi Version:          $(Gurobi.GRB_VERSION_MAJOR).$(Gurobi.GRB_VERSION_MINOR).$(Gurobi.GRB_VERSION_TECHNICAL)\n")
end

"for each timestep, insert a number and the respective data into a vector that is returned"
function get_timesteps_and_data(data)
    if haskey(data, "nw")
        return [(nw, data["nw"][string(nw)]) for nw in 1:length(data["nw"])]
    else
        return [(0, data)]
    end
end

"Prints redispatch. If `print_all`, then all gens (even those w/o re-dispatch) are listed with their generation."
function print_pg_delta(casename, pm_data; print_all=false)
    gen_changes = Dict{Int64, Dict{String, Vector}}()
    redispatch_cost = 0
    coeff = pm_data["cost_coeffs"]["redispatch"]
    # collect data
    data = get_timesteps_and_data(pm_data)
    log_str = "Planned re-dispatch for case $casename:\n"
    for (timestep, timestep_data) in data
        push!(gen_changes, timestep => Dict{String, Vector}())
        log_str *= "\033[3m Time step $timestep\033[0m\n"
        max_name_len = maximum([length(haskey(gen_data, "name") && PARAMS[:PRINTING][:CONFIDENTIAL] ? "$(gen_data["name"]) ($gen_id)" : "$gen_id") for (gen_id, gen_data) in timestep_data["gen"]])
        for (gen_id, gen_data) in timestep_data["gen"]
            if haskey(gen_data, "pg_delta_plus_new") # only print redispatch for non-virtual gens that were changed at all
                if gen_data["pg_delta_plus_new"] > PARAMS[:PRECISION][:POWER_PREC] || gen_data["pg_delta_minus_new"] > PARAMS[:PRECISION][:POWER_PREC] # only print if any re-dispatch took place
                    gen_id = haskey(gen_data, "name") && PARAMS[:PRINTING][:CONFIDENTIAL] ? "$(gen_data["name"]) ($gen_id)" : gen_id
                    pg_0    = round(gen_data["pg_0"], digits=2)
                    pg_new  = round(gen_data["pg_new"], digits=2)
                    pg_diff = round(pg_new - pg_0, digits=2)
                    if haskey(gen_data, "gsk_name")
                        if ! haskey(gen_changes[timestep], gen_data["gsk_name"])
                            push!(gen_changes[timestep], gen_data["gsk_name"] => [gen_data["index"], pg_0, pg_new])
                        else
                            gen_changes[timestep][gen_data["gsk_name"]][2] += pg_0
                            gen_changes[timestep][gen_data["gsk_name"]][3] += pg_new
                        end
                    else
                        push!(gen_changes[timestep], gen_id => [gen_data["index"], pg_0, pg_new])
                    end
                    cost    = round(gen_data["pg_delta_plus_new"] * (gen_data["cost_change"] + gen_data["cost_produce"]) + gen_data["pg_delta_minus_new"] * (gen_data["cost_change"] - gen_data["cost_produce"]), digits=2)
                    rdpp    = round(gen_data["RDP+"], digits=2)
                    rdpm    = round(gen_data["RDP-"], digits=2)
                    pmax    = round(gen_data["pmax"], digits=2)
                    pmin    = round(gen_data["pmin"], digits=2)
                    log_str *= "  Gen $(rpad("$gen_id:", max_name_len+1, " ")) $(lpad(pg_0, 5, " ")) => $(lpad(pg_new, 5, " ")) (re-dispatch $(lpad(pg_diff, 5, " ")) (up to: RDP- = $(lpad(rdpm, 5, " ")), RDP+ = $(lpad(rdpp, 5, " ")), pmin = $(lpad(pmin, 5, " "))), pmax = $(lpad(pmax, 5, " ")), cost $(lpad(cost, 8, " ")))\n"
                    if Memento.getlevels(_LOGGER)[Memento.getlevel(_LOGGER)] <= Memento.getlevels(_LOGGER)["debug"]
                        log_str *= "   $gen_data\n"
                    end
                    redispatch_cost += cost # do not use coeff, since real cost should be shown here (coeff printed later)
                elseif print_all # or if `print_all` is active
                    gen_id = haskey(gen_data, "name") && PARAMS[:PRINTING][:CONFIDENTIAL] ? "$(gen_data["name"]) ($gen_id)" : gen_id
                    pg_0    = round(gen_data["pg_0"], digits=2)
                    pg_new  = round(gen_data["pg_new"], digits=2)
                    pg_diff = round(pg_new - pg_0, digits=2)
                    cost    = round(gen_data["pg_delta_plus_new"] * (gen_data["cost_change"] + gen_data["cost_produce"]) + gen_data["pg_delta_minus_new"] * (gen_data["cost_change"] - gen_data["cost_produce"]), digits=2)
                    rdpp    = round(gen_data["RDP+"], digits=2)
                    rdpm    = round(gen_data["RDP-"], digits=2)
                    pmax    = round(gen_data["pmax"], digits=2)
                    pmin    = round(gen_data["pmin"], digits=2)
                    log_str *= "  Gen $(rpad("$gen_id:", max_name_len+1, " ")) $(lpad(pg_0, 5, " ")) => $(lpad(pg_new, 5, " ")) (re-dispatch $(lpad(pg_diff, 5, " ")) (up to: RDP- = $(lpad(rdpm, 5, " ")), RDP+ = $(lpad(rdpp, 5, " ")), pmin = $(lpad(pmin, 5, " "))), pmax = $(lpad(pmax, 5, " ")), cost $(lpad(cost, 8, " ")))\n"
                    if Memento.getlevels(_LOGGER)[Memento.getlevel(_LOGGER)] <= Memento.getlevels(_LOGGER)["debug"]
                        log_str *= "   $gen_data\n"
                    end
                end
            end
        end
    end
    internal_cost = round(coeff * redispatch_cost, digits=4)
    log_str *= " Cost of re-dispatch: $(round(redispatch_cost, digits=4)) (internal cost of $internal_cost)\n"
    if ! haskey(pm_data, "objective_costs")
        pm_data["objective_costs"] = Dict{String, Float64}()
    end
    pm_data["objective_costs"]["redispatch"] = internal_cost

    Memento.info(_LOGGER, log_str)
    return gen_changes
end

function print_topology_changes(casename, pm_data)::Dict{Int64, Dict{String, Pair{Int, Int}}}
    line_switches = Dict{Int64, Dict{String, Pair{Int, Int}}}()
    switching_cost = 0
    coeff = pm_data["cost_coeffs"]["br_sw"]
    # collect data
    data = get_timesteps_and_data(pm_data)
    log_str = "Planned topology actions for case $casename:\n"
    for (timestep, timestep_data) in data
        push!(line_switches, timestep => Dict{String, Pair{Int, Int}}())
        log_str *= "\033[3m Time step $timestep\033[0m\n"
        for (branch_id, branch_data) in sort!(OrderedDict(timestep_data["branch"]), by = (x -> parse(Int64, x)))
            # only print if not a virtual (_sec) line 
            if ! occursin("(", branch_id) && haskey(branch_data, "pf_new")
                branch_id = haskey(branch_data, "name") && PARAMS[:PRINTING][:CONFIDENTIAL] ? "$(branch_data["name"]) ($branch_id)" : branch_id
                if haskey(branch_data, "branch_name_intern")
                    branch_id = "$(branch_data["branch_name_intern"]) ($branch_id)"
                end
                # only print if change occurred
                if haskey(branch_data, "br_change_new") && isapprox(branch_data["br_change_new"], 1) && ! isapprox(branch_data["br_status"], branch_data["br_status_new"])
                    push!(line_switches[timestep], branch_id => (trunc(branch_data["br_status"]) => trunc(branch_data["br_status_new"])))
                    cost = round(coeff * 1, digits=2)
                    line_loading = abs(round(get(branch_data, "pf_start", branch_data["pf_new"]) / branch_data["rate_a"] * 100, digits=2))
                    log_str *= "  Branch $(rpad(branch_id, 36, " ")): $(round(branch_data["br_status"], digits=0)) => $(round(branch_data["br_status_new"], digits=0)) (cost = $cost, line loading = $line_loading %)\n"
                    switching_cost += cost
                    if Memento.getlevels(_LOGGER)[Memento.getlevel(_LOGGER)] <= Memento.getlevels(_LOGGER)["debug"]
                        log_str *= "   $branch_data\n"
                    end
                end
            end
        end
    end
    log_str *= " Internal cost of topology actions: $(round(switching_cost, digits=4))\n"
    if ! haskey(pm_data, "objective_costs")
        pm_data["objective_costs"] = Dict{String, Float64}()
    end
    pm_data["objective_costs"]["br_sw"] = switching_cost

    Memento.info(_LOGGER, log_str)
    return line_switches
end

function print_transformer_changes(casename, pm_data)::Dict{Int64, Dict{String, Pair{Float64, Float64}}}
    pst_shifts = Dict{Int64, Dict{String, Pair{Float64, Float64}}}()
    shifting_cost = 0
    coeff = pm_data["cost_coeffs"]["tr_sw"]
    # collect data
    data = get_timesteps_and_data(pm_data)
    log_str = "Planned transformer actions for case $casename:\n"
    for (timestep, timestep_data) in data
        push!(pst_shifts, timestep => Dict{String, Pair{Float64, Float64}}())
        log_str *= "\033[3m Time step $timestep\033[0m\n"
        for (branch_id, branch_data) in timestep_data["branch"]
            # only print if not a virtual (_sec) line 
            if ! occursin("(", branch_id)
                branch_id = haskey(branch_data, "name") && PARAMS[:PRINTING][:CONFIDENTIAL] ? "$(branch_data["name"]) ($branch_id)" : branch_id
                if haskey(branch_data, "branch_name_intern")
                    branch_id = "$(branch_data["branch_name_intern"]) ($branch_id)"
                end
                # only print if change occurred
                if haskey(branch_data, "tr_shift_new") && ! isapprox(branch_data["tr_shift_new"], branch_data["previous_shift_step"], atol=PARAMS[:PRECISION][:POWER_PREC])
                    push!(pst_shifts[timestep], branch_id => (branch_data["shift"] => branch_data["tr_shift_new"]))
                    cost = round(coeff * 1, digits=4)
                    old_lvl     = round((branch_data["previous_shift_step"]), digits=4)
                    new_lvl     = round((branch_data["tr_shift_new"]), digits=4)
                    old_shift   = round((branch_data["previous_shift_step"] * branch_data["shift_per_step"] * 180/π), digits=4)
                    new_shift   = round((branch_data["tr_shift_new"] * branch_data["shift_per_step"]) * 180/π, digits=4)
                    log_str *= "  Transformer $(rpad(branch_id, 36, " ")): Step $old_lvl => Step $new_lvl ($(old_shift)° => $(new_shift)°, cost = $cost)\n"
                    shifting_cost += cost
                    if Memento.getlevels(_LOGGER)[Memento.getlevel(_LOGGER)] <= Memento.getlevels(_LOGGER)["debug"]
                        log_str *= "   $branch_data\n"
                    end
                end
            end
        end
    end
    log_str *= " Internal cost of transformer shifting actions: $(round(shifting_cost, digits=4))\n"
    if ! haskey(pm_data, "objective_costs")
        pm_data["objective_costs"] = Dict{String, Float64}()
    end
    pm_data["objective_costs"]["tr_sw"] = shifting_cost

    Memento.info(_LOGGER, log_str)
    return pst_shifts
end

"print contingencies (i.e., η variables) that cause cost (i.e., η > 0)"
function print_line_overloads_benders(casename, pm_data, result; print_all=false)
    slack_cost = 0
    slack_cost_base_case = 0
    coeff = pm_data["cost_coeffs"]["overload"]
    log_str = "Planned overloaded lines for case $casename:\n"
    nws = length(result["solution"]["nw"])
    for timestep in 1:nws
        timestep = string(timestep)
        timestep_data = result["solution"]["nw"][timestep]
        log_str *= "\033[3m Time step $timestep\033[0m\n"
        if ! haskey(timestep_data, "eta")
            continue
        end
        max_name_len = maximum([length("$(get_contingency_names(pm_data, pm_data["nw"]["$timestep"]["contingencies"][parse(Int64, c_id)]))") for (c_id, _) ∈ timestep_data["eta"]])
        max_id_len   = maximum([length("$c_id") for (c_id, _) ∈ timestep_data["eta"]])
        for (c_id, eta_data) in timestep_data["eta"]
            eta_val = eta_data["η"]
            if round(eta_val, digits=SOLVER_PRECISION_DIGITS) > 0
                cost = eta_val # in Benders, we have already multiplied the coeff with the eta_val in the Benders SP obj func
                contingency = pm_data["nw"]["$timestep"]["contingencies"][parse(Int64, c_id)]
                log_str *= ("  Contingency $(rpad("$(get_contingency_names(pm_data, contingency))", max_name_len, " ")) $(rpad("($c_id)", max_id_len+2, " ")) causes overloads with cost $cost.\n")
                slack_cost += cost
                if isempty(contingency)
                    slack_cost_base_case += cost
                end
            end
        end
    end
    log_str *= " Internal cost of slack: $(round(slack_cost, digits=4))\n"
    if ! haskey(pm_data, "objective_costs")
        pm_data["objective_costs"] = Dict{String, Float64}()
    end
    pm_data["objective_costs"]["overload_all"] = slack_cost
    pm_data["objective_costs"]["overload_base_case"] = slack_cost_base_case

    Memento.info(_LOGGER, log_str)
end


"Prints lines that are overloaded, and the contingency (if applicable) in which they are"
function print_line_overloads(casename, data, result; print_all=false, for_n1=false)
    slack_cost = 0
    slack_cost_base_case = 0
    coeff = data["cost_coeffs"]["overload"]
    log_str = "Planned overloaded lines for case $casename:\n"
    nws = length(result["solution"]["nw"])
    for timestep in 1:nws
        timestep = string(timestep)
        timestep_data = result["solution"]["nw"][timestep]
        log_str *= "\033[3m Time step $timestep\033[0m\n"
        for (branch_id, branch_result_data) in timestep_data["branch"]
            # get index of branch and, if existing, contingency
            ids     = split(replace(branch_id, "(" => "", ")" => "", " " => ""), ",")
            b_id    = parse(Int, ids[1])
            c_id    = parse(Int, get(ids, 2, "-1"))
            branch_data = data["nw"][timestep]["branch"][string(b_id)]
            f_bus = branch_data["f_bus"]
            t_bus = branch_data["t_bus"]
            contingency = c_id == -1 ? Contingency() : data["nw"][timestep]["contingencies"][c_id]
            # if data has a "monitoring" list, then the name hast to be in it. If no such list exists, always print
            if branch_data["monitored"] && haskey(branch_result_data, "pf")#  haskey(data["nw"][timestep], "monitoring") || Dict(:branch => [b_id]) ∈ data["nw"][timestep]["monitoring"]
                rate_a = round(get_adjusted_rate_a(data, branch_data, for_n1=(c_id ≠ -1 || for_n1)), digits=4) + SOLVER_PRECISION
                pf     = round(abs(branch_result_data["pf"]), digits=4)
                pt     = round(abs(branch_result_data["pt"]), digits=4)
                if print_all || pf >= rate_a || pt >= rate_a
                    b_id    = haskey(branch_data, "name") ? "$(branch_data["name"]) ($b_id)" : b_id
                    if haskey(branch_data, "branch_name_intern")
                        b_id = "$(branch_data["branch_name_intern"]) $(b_id)"
                    end
                    load    = round((pf - rate_a) / rate_a * 100 + 100, digits = 2)
                    cost    = abs(round(coeff * (abs(pf) - rate_a) * (1/rate_a + 0.0001), digits=4))
                    log_str *= ("  Branch $(b_id) rating exceeded in contingency $(get_contingency_names(data, contingency)) ($(c_id))\n   $load % (p = $pf, rate_a = $rate_a, cost = $cost)\n")
                    log_str *= "$branch_data\n\n"
                    slack_cost += cost
                    if isempty(contingency)
                        slack_cost_base_case += cost
                    end
                end
            end
        end
    end
    log_str *= " Internal cost of slack: $(round(slack_cost, digits=4))\n"
    if ! haskey(data, "objective_costs")
        data["objective_costs"] = Dict{String, Float64}()
    end
    data["objective_costs"]["overload_all"] = slack_cost
    data["objective_costs"]["overload_base_case"] = slack_cost_base_case

    Memento.info(_LOGGER, log_str)
end


"Prints shutdown/startup of gens. If `print_all`, then all gens (even those w/o re-dispatch) are listed with their set status."
function print_generator_changes(casename, pm_data; print_all=false)
    gen_switching_cost = 0
    coeff = pm_data["cost_coeffs"]["gen_sw"]
    # collect data
    data = get_timesteps_and_data(pm_data)
    log_str = "Planned generator actions for case $casename:\n"
    for (timestep, timestep_data) in data
        log_str *= "\033[3m Time step $timestep\033[0m\n"
        max_name_len = maximum([length(haskey(gen_data, "name") && PARAMS[:PRINTING][:CONFIDENTIAL] ? "$(gen_data["name"]) ($gen_id)" : "$gen_id") for (gen_id, gen_data) in timestep_data["gen"]])
        for (gen_id, gen_data) in timestep_data["gen"]
            # only print redispatch for non-virtual gens
            if "gen_shutdown_new" ∉ collect(keys(gen_data))
                continue
            end
            gen_id = haskey(gen_data, "name") && PARAMS[:PRINTING][:CONFIDENTIAL] ? "$(gen_data["name"]) ($gen_id)" : gen_id
            # only print if change occurred
            if isapprox(gen_data["gen_shutdown_new"], 1.0, atol=0.1) || isapprox(gen_data["gen_startup_new"], 1.0, atol=0.1) || print_all
                gen_status_old = round(gen_data["gen_shutdown_new"], digits=0)
                gen_status_new = round(gen_data["gen_startup_new"], digits=0)
                cost = round(coeff * 1, digits=4) #* (gen_data["shutdown"] * gen_data["gen_shutdown_new"] + gen_data["startup"] * gen_data["gen_startup_new"])
                log_str *= "  Gen $(rpad("$gen_id:", max_name_len+2, " ")) $(round(gen_status_old, digits=0)) => $(round(gen_status_new, digits=0)) (cost = $cost)\n"
                gen_switching_cost += cost
                if Memento.getlevels(_LOGGER)[Memento.getlevel(_LOGGER)] <= Memento.getlevels(_LOGGER)["debug"]
                    log_str *=  "   $gen_data\n"
                end
            end
        end
    end
    log_str *= " Internal cost of generator startup and shutdown actions: $(round(gen_switching_cost, digits=4))\n"
    if ! haskey(pm_data, "objective_costs")
        pm_data["objective_costs"] = Dict{String, Float64}()
    end
    pm_data["objective_costs"]["gen_sw"] = gen_switching_cost

    Memento.info(_LOGGER, log_str)
end

function print_flow_change(casename, data)
    if PARAMS[:PRINTING][:FLOW_CHANGE]
        if :ratea ∈ data["slack_types"]
            if branch_data["rat_s_new"] > 0
                slack_info = ", with slack of $(branch_data["rat_s"]) from the limit of $(branch_data["rate_a"])"
            elseif branch_data["raf_s_new"] > 0
                slack_info = ", with slack of $(branch_data["raf_s"]) from the limit of $(branch_data["rate_a"])"
            else
                slack_info = ""
            end
        else
            slack_info = ""
        end
        log_str *= "  Branch $(lpad(branch_id, 4)) flow: $(round(branch_data["pf_start"], digits=2))"
        if round(branch_data["pf_start"], digits=1) != round(branch_data["pf_new"], digits=1)
            log_str *= " => $(round(branch_data["pf_new"], digits=2))"
        end
        log_str *= "$slack_info\n"
    end
end

function print_slack_values(casename, data, result, slack_types)
    slacks_used = false
    log_str = ""
    if "nw" ∈ collect(keys(data)) # multinetwork
        slack_data = Dict()
        for (network_id, network_data) in result["solution"]["nw"]
            slacks_nw = get_slacks_from_result(Dict("solution" => network_data), slack_types)
            # check if slacks were used at all to avoid unnecessary printing
            if size(slacks_nw, 1) > 0
                slacks_used = true
            end
            slack_data[network_id] = slacks_nw
        end
        if slacks_used
            log_str = "Slacks used for case $casename:\n"
            for (timestep, slacks) in slack_data
                log_str *= "\033[3m Time step $timestep\033[0m\n"
                for row in eachrow(slacks)
                    log_str *= "  $(row["component_type"]) $(row["component_id"]) $(row["slack_variable_name"]) = $(round(row["slack_variable_value"], digits=6))\n"
                end
            end
        end
    else # single network
        slacks = get_slacks_from_result(result, slack_types)
        # check if slacks were used at all to avoid unnecessary printing
        n_slacks_used = size(slacks, 1) 
        if n_slacks_used > 0
            log_str = "Slacks used for case $casename:\n"
            for row in eachrow(slacks)
                log_str *= "  $(row["component_type"]) $(row["component_id"]) $(row["slack_variable_name"]) = $(round(row["slack_variable_value"], digits=6))\n"
            end
        end
    end

    if log_str != ""
        Memento.warn(_LOGGER, log_str)
    end
end

"performs LODF-based k-medoids clustering for the contingencies, returns the result"
function contingency_kmedoids(nw_data, contingency_candidates, n_medoids::Int64, dist_metric)
    # compute LODF matrix
    LODF = compute_LODF(nw_data)

    return contingency_kmedoids(nw_data, LODF, contingency_candidates, n_medoids, dist_metric)
end

"performs LODF-based k-medoids clustering for the contingencies, returns the result"
function contingency_kmedoids(nw_data, LODF, contingency_candidates, n_medoids::Int64, dist_metric)
    LODF_vectors = zeros(Float64, length(contingency_candidates), size(LODF, 2))

    # sum up lodf rows for each branch of a contingency
    for (idx, contingency) in enumerate(contingency_candidates)
        branch_ids = contingency[:branch]
        valid_branch_ids = [b_id for b_id in branch_ids if 1 <= b_id <= size(LODF, 1)]
        isempty(valid_branch_ids) && continue

        # Compare contingency LODF row vectors to base case (zero vector)
        for b_id in valid_branch_ids
            LODF_vectors[idx, :] .+= LODF[b_id, :]
        end
        LODF_vectors[idx, :] ./= length(valid_branch_ids)
    end

    # calculate pairwise distances
    dist_matrix = compute_cut_distances(LODF_vectors, dist_metric)

    return Clustering.kmedoids(dist_matrix, n_medoids)
end

"diversity based contingency filtering for partial decomposition"
function get_relevant_contingencies_diversity(nw_data)::Tuple{Vector{Contingency}, Vector{Float64}}
    contingencies = Vector{Contingency}()
    violations = Vector{Float64}()

    contingency_candidates = values(nw_data["contingencies"])
    n_relevant = length(contingency_candidates)

    n_keep = n_master_contingencies(n_relevant)

    if n_relevant > 0 && n_keep > 0

        if n_keep >= n_relevant
            append!(contingencies, contingency_candidates)
            append!(violations, fill(1.0, n_relevant))
            return contingencies, violations
        end

        dist_metric = PARAMS[:BENDERS][:PARTIAL_DECOMP][:DISTANCE_METRIC]

        clustering = contingency_kmedoids(nw_data, contingency_candidates, n_keep, dist_metric)

        # k-medoids
        selected_indices = sort(clustering.medoids)

        for idx in selected_indices
            push!(contingencies, contingency_candidates[idx])
            push!(violations, 1.0)
        end
    end

    return contingencies, violations
end


"true if the partial decomposition mode derives contingency relevance from a solution instead of from the network structure"
function partial_decomp_needs_solution()::Bool
    mode = PARAMS[:BENDERS][:PARTIAL_DECOMP][:MODE]
    mode == :OVERLOAD && return true
    mode == :SW_HEUR  || return false # :DIVERSITY and :RND are structural

    weights = PARAMS[:BENDERS][:PARTIAL_DECOMP][:WEIGHTS]
    return any(weights[name] > 0 for name in (:SENS, :CONGPROX, :ETA)) # :RND and :FULL are structural
end

"true if subproblem selection or cut selection uses the LODF based contingency clustering"
function sp_diversity_clustering_active()::Bool
    sp_selection = PARAMS[:BENDERS][:PARTIAL][:SORTING_MODE] == :SCORE && PARAMS[:BENDERS][:PARTIAL][:SCORE][:diversity] > 0
    return sp_selection || PARAMS[:BENDERS][:CUTS][:SELECTION][:METHOD] == :LODF_DIVERSITY
end

"cap on the number of contingencies whose power flow is kept in the master problem"
function n_master_contingencies(n_contingencies::Int64)::Int64
    pct = PARAMS[:BENDERS][:PARTIAL_DECOMP][:N_MAX_CONTINGENCIES_PCT]
    n_keep_req = pct > 0 ? ceil(Int, pct * n_contingencies) : PARAMS[:BENDERS][:PARTIAL_DECOMP][:N_MAX_CONTINGENCIES]
    return min(n_keep_req, n_contingencies)
end

"determine the relevant contingencies, sorted by descending relevance. `result` is only required for modes for which `partial_decomp_needs_solution()` holds"
function get_relevant_contingencies(data; result=nothing, min_violation = PARAMS[:BENDERS][:PARTIAL_DECOMP][:MIN_VIOLATION])::Vector{Contingency}
    contingencies = Vector{Contingency}()
    violations = Vector{Float64}()

    mode = PARAMS[:BENDERS][:PARTIAL_DECOMP][:MODE]
    needs_solution = partial_decomp_needs_solution()

    if needs_solution && (isnothing(result) || !haskey(result["solution"], "nw"))
        @warn("Partial decomposition mode $mode derives relevance from a solution, but none is available. Skipping partial decomposition.")
        return contingencies
    end

    for timestep in 1:length(data["nw"])
        nw_data = data["nw"]["$timestep"]

        if mode == :DIVERSITY
            try
                timestep_contingencies, timestep_violations = get_relevant_contingencies_diversity(nw_data)
                append!(contingencies, timestep_contingencies)
                append!(violations, timestep_violations)
            catch e
                @warn("Error while calculating partial decomposition diversity scores: $e.")
            end
            continue
        end

        etas = Dict{String, Any}()
        if needs_solution
            timestep_data = result["solution"]["nw"]["$timestep"]
            if ! haskey(timestep_data, "eta")
                @warn("Missing solution values, skipping partial decomposition!")
                continue
            end
            etas = timestep_data["eta"]
        end

        # if SW_HEUR is the indicator for contingency relevance, we need to calculate the scores first
        if mode == :SW_HEUR
            calc_and_set_br_sw_heur_scores!(nw_data,
                                        only_switchable=false,
                                        set_switchability=false,
                                        nw_sol=needs_solution ? get(data, "result", result)["solution"]["nw"]["$timestep"] : Dict{String, Any}(),
                                        params=PARAMS[:BENDERS][:PARTIAL_DECOMP],
                                        heur_params=PARAMS[:BENDERS][:PARTIAL_DECOMP][:HEUR])
        end

        branch_data = nw_data["branch"]

        # for each contingency, determine random score, overload (eta) score or (precalculated) switching heuristic score to determine relevance
        for (c_id, contingency) in enumerate(nw_data["contingencies"])
            all(isempty, values(contingency)) && continue # base case, never a candidate

            if     mode == :RND
                push!(contingencies, contingency)
                push!(violations, rand())
            elseif mode == :OVERLOAD
                haskey(etas, "$c_id") || continue
                eta_val = etas["$c_id"]["η"]
                if round(eta_val, digits=SOLVER_PRECISION_DIGITS) > 0
                    push!(contingencies, contingency)
                    push!(violations, eta_val)
                end
            elseif mode == :SW_HEUR
                violation = 0.0
                for branch_id in contingency[:branch]
                    branch = branch_data["$branch_id"]
                    violation += branch["sw_heur_score"]
                end

                push!(contingencies, contingency)
                push!(violations, violation/length(contingency[:branch]))
            else
                @assert(false, "Unknown partial decomposition mode $mode (expected one of :RND, :OVERLOAD, :SW_HEUR, :DIVERSITY).")
            end
        end
    end

    # sort the contingencies by their violation (eta_val value) in descending order
    sorted_indices = sortperm(violations, rev=true)
    contingencies = contingencies[sorted_indices]

    # if min violation given, filter out contingencies with less violation
    if min_violation > 0 && mode != :DIVERSITY
        contingencies = [contingency for (i, contingency) in enumerate(contingencies) if violations[sorted_indices[i]] >= min_violation]
    end

    return contingencies
end

"cluster the contingencies left to the subproblems by LODF similarity, used by diversity based subproblem and cut selection"
function compute_sp_diversity_clusters!(nw_data, contingencies_relevant::Vector{Contingency})
    LODF = zeros(1,1)
    try
        LODF = compute_LODF(nw_data)
    catch e
        @warn("Error while calculating LODF for subproblem selection diversity cluster: $e.")
        return
    end

    # get the remaining contingencies to cluster on; must match the subproblem list built in build_secure_ots_time_benders
    contingency_candidates = nw_data["contingencies"]
    n_keep = min(n_master_contingencies(length(contingency_candidates)), length(contingencies_relevant))

    contingencies_sp = n_keep > 0 ? setdiff(contingency_candidates, contingencies_relevant[1:n_keep]) : copy(contingency_candidates)

    n_medoids = min(PARAMS[:BENDERS][:PARTIAL][:SCORE_DIV][:N_MEDOIDS], length(contingencies_sp))
    n_medoids > 0 || return

    # clustering once, write the result into "sp_clusters"
    if !PARAMS[:BENDERS][:PARTIAL][:SCORE_DIV][:FULL_MEDOIDS] || PARAMS[:BENDERS][:CUTS][:SELECTION][:METHOD] == :LODF_DIVERSITY
        clustering = contingency_kmedoids(nw_data, LODF, contingencies_sp, n_medoids, PARAMS[:BENDERS][:PARTIAL][:SCORE_DIV][:DISTANCE_METRIC])

        clusters = [(clustering.medoids[i], Vector{Int64}()) for i in 1:n_medoids]::Vector{Tuple{Int64, Vector{Int64}}}

        for (sp_id, cluster_id) in enumerate(clustering.assignments)
            push!(clusters[cluster_id][2], sp_id)
        end

        nw_data["sp_clusters"] = clusters

    # perform clustering iteratively, each time removing the newly found clusters
    else
        medoidss = Vector{Vector{Int64}}()
        initial_indices = [i for i in 1:length(contingencies_sp)]

        while length(contingencies_sp) > n_medoids

            # compute clustering
            clustering = contingency_kmedoids(nw_data, LODF, contingencies_sp, n_medoids, PARAMS[:BENDERS][:PARTIAL][:SCORE_DIV][:DISTANCE_METRIC])

            medoids = sort!(clustering.medoids, rev=true)

            # add intial indices of medoids to final vector
            push!(medoidss, initial_indices[medoids])

            # update initial indices: the index of the deleted medoid now points to the next element
            for medoid_id in medoids
                for i in medoid_id:(length(contingencies_sp)-1)
                    initial_indices[i] = initial_indices[i+1]
                end
            end

            # delete medoids from contingencies to be considered
            deleteat!(contingencies_sp, sort!(medoids))

        end

        !isempty(contingencies_sp) && push!(medoidss, initial_indices[1:length(contingencies_sp)])

        nw_data["sp_medoidss"] = medoidss
    end
end

"save the relevant contingencies under the timestep-dependent contingencies_relevant field and cluster the remaining ones. Does nothing if already prepared"
function prepare_contingency_selection!(data; result=nothing)
    nw_data = data["nw"]["1"]
    haskey(nw_data, "contingencies_relevant") && return

    # nested Benders keeps no contingency in the master, so it solves all of them as subproblems
    application = get(data, "application", "")
    nested = application == "DC_SCOTS_TIME_SLACK_BENDERS_NESTED" || application == "DC_SCOTS_TIME_SLACK_BENDERS_BENDERS"

    # only determine relevance if a contingency (and not just the base case) is actually kept in the master
    select = PARAMS[:BENDERS][:PARTIAL_DECOMP][:ENABLED] && !nested && n_master_contingencies(length(nw_data["contingencies"])) > 0
    relevant_contingencies = select ? get_relevant_contingencies(data, result=result) : Vector{Contingency}()

    # TODO: fix for multiple timesteps
    for (_, timestep_data) in data["nw"]
        timestep_data["contingencies_relevant"] = relevant_contingencies
        timestep_data["contingencies_relevant2"] = relevant_contingencies
    end

    sp_diversity_clustering_active() && compute_sp_diversity_clusters!(nw_data, relevant_contingencies)
end


"get the total cost of all contingencies in Benders result"
function get_total_cost_of_contingencies(data, result)::Float64
    total_cost = 0.0
    nws = length(result["solution"]["nw"])
    for timestep in 1:nws
        timestep = string(timestep)
        timestep_data = result["solution"]["nw"][timestep]
        if ! haskey(timestep_data, "eta")
            continue
        end
        for (c_id, eta_data) in timestep_data["eta"]
            eta_val = eta_data["η"]
            if round(eta_val, digits=SOLVER_PRECISION_DIGITS) > 0
                total_cost += eta_val 
            end
        end
    end
    return total_cost
end

"check overloads for all given contingencies. If `for_validation`, calculate using original obj. func., if not, calculate using line minimization to not have unconnected lines with 100% in the result"
function check_relevant_overloads(data, result, solver_convex, contingencies::Vector{Contingency}; for_validation=false)
    contingencies_copy = deepcopy(data["nw"]["1"]["contingencies"])
    data["nw"]["1"]["contingencies"] = contingencies
    
    if isempty(contingencies)
        calculate_power_flow!(data, result, solver_convex)
        data["nw"]["1"]["contingencies"] = contingencies_copy
        data["nw"]["1"]["contingencies_relevant"] = contingencies

        return data, result, "DC_SCOTS_TIME_SLACK_BENDERS"
    end

    set_dc_pf_start_values!(data, result["solution"], fix_uc=true, fix_ts=true, fix_tr=true, fix_benders=false, fix_pf=false, fix_slack=false, benders=false)
    
    @info("Solving compact model with $(length(contingencies)) contingencies to check power flows in contingencies.")
    if for_validation
        @assert ! SLACK_EQUALS_SEC_SLACK
        @warn("Use to validate Benders results only.")
        t = @elapsed begin result_compact = solve_secure_ots_time(data, DCPPowerModel, solver_convex,  weighted_obj_func_factory([objective_min_num_gen_switch, objective_min_num_tr_switch, objective_min_num_branch_switch, objective_min_redispatch_cost, objective_min_num_redispatch, objective_min_slack_relative_linear, objective_min_num_slacks], [data["cost_coeffs"]["gen_sw"], data["cost_coeffs"]["tr_sw"], data["cost_coeffs"]["br_sw"], data["cost_coeffs"]["redispatch"], data["cost_coeffs"]["num_redispatch"], data["cost_coeffs"]["overload"], data["cost_coeffs"]["num_overload"]]), slack_types=data["slack_types"]) end
        @assert(isapprox(result_compact["objective"], result["objective"], atol=SOLVER_PRECISION), "Objective function values do not match: $(result_compact["objective"]) != $(result["objective"])")
    else
        t = @elapsed begin result_compact = solve_secure_ots_time(data, DCPPowerModel, solver_convex, objective_min_line_load, slack_types=data["slack_types"]) end
    end
    result_compact["time"] = t

    print_line_overloads("case", data, result_compact)
    application = "DC_SCOTS_TIME_SLACK"

    if for_validation
        # compare result and result compact w.r.t. branch load (p)
        for (b_id, branch) ∈ data["nw"]["1"]["branch"]
            if is_branch_monitored(branch)
                pf_compact = result_compact["solution"]["nw"]["1"]["branch"][string(b_id)]["pf"]
                pf_benders = result["solution"]["nw"]["1"]["branch"][string(b_id)]["pf"]
                if abs(pf_compact - pf_benders) > POWER_PRECISION
                    @warn("Branch $b_id: Compact: $pf_compact, Benders: $pf_benders")
                end
            end
        end
    end

    result["solution"] = result_compact["solution"]

    data["nw"]["1"]["contingencies"] = contingencies_copy # restore original contingencies
    data["nw"]["1"]["contingencies_relevant"] = contingencies
    
    return data, result, application
end

"check overloads for all given contingencies by running the compact formulation for each contingency separately (fixing the br_status/status/switchable everytime)"
function check_relevant_overloads_individual(data, result, solver_convex, contingencies)
    contingencies_copy = deepcopy(data["nw"]["1"]["contingencies"])
    data["nw"]["1"]["contingencies_relevant"] = contingencies
    data["nw"]["1"]["contingencies"] = Vector{Contingency}()
    total_costs_compact = Vector{Float64}([0.0 for _ in 1:(length(contingencies) + 1)]) # +1 for the base case
    results_compact = Vector{Dict{String, Any}}([Dict{String, Any}() for _ in 1:(length(contingencies) + 1)])
    
    # get number of threads
    #n_threads = Threads.nthreads()

    #datas = [deepcopy(data) for _ in 1:n_threads]
    #results = [deepcopy(result) for _ in 1:n_threads]

    for (contingency_id, contingency) ∈ collect(enumerate([Dict(:branch => []); contingencies]))
        #result = results[threadid()]
        #data = datas[threadid()]
        # check status of contingency branches
        tmp_status = Dict{Int, Int}()
        for branch_id in contingency[:branch]
            branch_id_str = string(branch_id)
            if haskey(result["solution"]["nw"]["1"]["branch"][branch_id_str], "br_status")
                tmp_status[branch_id] = deepcopy(result["solution"]["nw"]["1"]["branch"][branch_id_str]["br_status"])
            else
                tmp_status[branch_id] = deepcopy(data["nw"]["1"]["branch"][branch_id_str]["br_status"])
            end
            result["solution"]["nw"]["1"]["branch"][branch_id_str]["br_status"] = 0
            data["nw"]["1"]["branch"][branch_id_str]["br_status"] = 0
        end

        ## solve the problem
        set_dc_pf_start_values!(data, result["solution"], fix_uc=true, fix_ts=true, fix_tr=true, fix_benders=false, fix_pf=false, fix_slack=false, benders=false)
        # compact
        #@info("Solving compact model:")
        data["SUBPROBLEM"] = true
        result_compact = solve_secure_ots_time(data, DCPPowerModel, solver_convex, weighted_obj_func_factory([objective_min_num_gen_switch, objective_min_num_tr_switch, objective_min_num_branch_switch, 
                                        objective_min_redispatch_cost, objective_min_num_redispatch, objective_min_slack_relative_linear, objective_min_num_slacks], 
                                        [data["cost_coeffs"]["gen_sw"], data["cost_coeffs"]["tr_sw"], data["cost_coeffs"]["br_sw"], data["cost_coeffs"]["redispatch"], 
                                        data["cost_coeffs"]["num_redispatch"], data["cost_coeffs"]["overload"], data["cost_coeffs"]["num_overload"]]), slack_types=data["slack_types"])
        results_compact[contingency_id] = result_compact
        total_costs_compact[contingency_id] = result_compact["objective"]
        data["SUBPROBLEM"] = false
        @info("Solved compact model for contingency $(get_contingency_names(data, contingency))")
        
        # reset br_status for the contingency branches
        for branch in contingency[:branch]
            result["solution"]["nw"]["1"]["branch"][string(branch)]["br_status"] = tmp_status[branch]
            data["nw"]["1"]["branch"][string(branch)]["br_status"] = tmp_status[branch]
        end
    end

    for (contingency_id, contingency) in enumerate([Dict(:branch => []); contingencies])
        result_compact = results_compact[contingency_id]
        for (branch_id, branch_data) in result_compact["solution"]["nw"]["1"]["branch"]
            b_c_id = "($(branch_id), $(contingency_id-1))"
            #if haskey(result["solution"]["nw"]["1"]["branch"], branch_id)
                result["solution"]["nw"]["1"]["branch"][b_c_id] = branch_data
            #end
        end
    end
    
    total_cost_compact = sum(total_costs_compact)

    # move branch solution data for base power flow into contingency id 1
    # for (branch_id, branch_data) in result["solution"]["nw"]["1"]["branch"]
        # b_c_id = "$(branch_id)"
        # result["solution"]["nw"]["1"]["branch"][b_c_id] = result["solution"]["nw"]["1"]["branch"][branch_id]
        # delete!(result["solution"]["nw"]["1"]["branch"], branch_id) # remove original branch id
    # end

    println("========================================")
    @info("Comparison of objective function values:")
    println("Total cost of compact model: $(total_cost_compact)")
    println("Total cost of Benders model: $(result["objective"])")
    if PARAMS[:LIMITS][:REDISPATCH][:PG_DELTA_NUM_DAY] == 0
        if ! isapprox(total_cost_compact, result["objective"], atol=SOLVER_PRECISION)
            @error("Objective function values do not match: $(total_cost_compact) != $(result["objective"])")
        else
            @info("Objective function value check successful!")
        end
    end
    println("========================================")
    data["nw"]["1"]["contingencies"] = contingencies_copy # restore original contingencies
end

"check overloads for all contingencies that have cost in Benders result"
function check_relevant_overloads(data, result, solver_convex; for_validation=false)
    return check_relevant_overloads(data, result, solver_convex, get_relevant_contingencies(data, result=result), for_validation=for_validation)
end

"check solution using compact formulation. If base_only, just calculate the power flow."
function calculate_power_flow!(data, result, solver_convex; base_only=true)
    @info("Calculating power flow for case $(data["name"]) with base_only = $base_only")
    data["SUBPROBLEM"] = false
    #data["nw"]["1"]["SUBPROBLEM"] = false

    # set the base power flow values
    is_nested = occursin("NESTED", data["application"]) || occursin("BENDERS_BENDERS", data["application"])
    set_dc_pf_start_values!(data, result["solution"], fix_uc=false, fix_ts=true, fix_benders=false, fix_pf=false, fix_slack=false, fix_tr=true, benders=false)

    # remove empty contingency
    if ! PARAMS[:BENDERS][:BASE_PF_IN_MASTER]
        contingencies_copy = filter!(c -> ! isempty(c), deepcopy(data["nw"]["1"]["contingencies"]))
    else
        contingencies_copy = filter(c -> ! isempty(c), deepcopy(data["nw"]["1"]["contingencies"]))
    end

    # if we want to check the base case only, set contingencies to empty 
    if base_only
        data["nw"]["1"]["contingencies"] = Vector{Contingency}()
    end
    
    # set PG_EQUALS_PGSEC = true and SLACK_EQUALS_SEC_SLACK = false
    PARAMS[:CONTINGENCIES][:PG_EQUALS_PGSEC] = true
    PARAMS[:CONTINGENCIES][:SLACK_EQUALS_SEC_SLACK] = false

    # solve the power flow to obtain values for branches and buses (non-Benders, without contingencies)
    result_pf = solve_secure_ots_time(data, DCPPowerModel, solver_convex, weighted_obj_func_factory([objective_min_gen_change, objective_min_num_gen_switch, objective_min_num_tr_switch, objective_min_num_branch_switch, objective_min_redispatch_cost, objective_min_num_redispatch, objective_min_slack_relative_linear, objective_min_num_slacks], [is_nested ? 0 : 100000, data["cost_coeffs"]["gen_sw"], data["cost_coeffs"]["tr_sw"], data["cost_coeffs"]["br_sw"], data["cost_coeffs"]["redispatch"], data["cost_coeffs"]["num_redispatch"], data["cost_coeffs"]["overload"], data["cost_coeffs"]["num_overload"]]), slack_types=data["slack_types"])
    if result_pf["termination_status"] == JuMP.MOI.INFEASIBLE || result_pf["termination_status"] == JuMP.MOI.INFEASIBLE_OR_UNBOUNDED
        @error("Power flow calculation failed: $(result_pf["termination_status"]). Result from main problem run should also have been infeasible!")
        result["termination_status"] = "CUT_ERROR"
        result["primal_status"] = "CUT_ERROR"

        return
        # JuMP.compute_conflict!(pm.model)
        # iis_model, _ = JuMP.copy_conflict(pm.model)
        # print(iis_model)
        # if PARAMS[:WRITING][:IIS]
            # JuMP.write_to_file(iis_model, "conflict_model.lp")
        # end
    end
    @assert(result_pf["termination_status"] == JuMP.MOI.OPTIMAL, "Power flow calculation failed: $(result_pf["termination_status"]). Result from main problem run should also have been infeasible!")

    if (PARAMS[:POSTPROCESSING][:CHECKS][:CHECK_WITH_COMPACT] && ! base_only)
        if ! isapprox(result["objective"], result_pf["objective"], atol=PARAMS[:POSTPROCESSING][:CHECKS][:CHECK_WITH_COMPACT_TOL])
            @error("Objective function values do not match: $(result["objective"]) != $(result_pf["objective"])")
        else
            @info("Objective function value check successful: $(result["objective"]) == $(result_pf["objective"])")
        end
    else
        if result["objective"] < result_pf["objective"] - SOLVER_PRECISION
            @error("Base power flow already more expensive than the original problem: $(result["objective"]) < $(result_pf["objective"])")
            result["termination_status"] = "CUT_ERROR"
            result["primal_status"] = "CUT_ERROR"
        else
            @info("Objective function value check successful: $(result["objective"]) >= $(result_pf["objective"])")
        end
    end

    # insert all pf values from result_pf into result
    if base_only
        for (nw, nw_result_pf) in result_pf["solution"]["nw"]
            # add branch data
            if ! haskey(result["solution"]["nw"][nw], "branch")
                result["solution"]["nw"][nw]["branch"] = nw_result_pf["branch"]
            else
                for (b_id, nw_result_pf_branch) in nw_result_pf["branch"]
                    if ! haskey(result["solution"]["nw"][nw]["branch"], b_id)
                        result["solution"]["nw"][nw]["branch"][b_id] = Dict{String, Any}()
                    end
                    nw_result_branch = result["solution"]["nw"][nw]["branch"][b_id]

                    # add electrical data (flow and rate_a) if branch is electrical
                    if haskey(nw_result_pf_branch, "pf")
                        nw_result_branch["pf"] = nw_result_pf_branch["pf"]
                        nw_result_branch["pt"] = nw_result_pf_branch["pt"]

                        if haskey(nw_result_pf_branch, "rat_s_b")
                        nw_result_branch["rat_s_b"] = nw_result_pf_branch["rat_s_b"]
                        end
                        if haskey(nw_result_pf_branch, "rat_s")
                        nw_result_branch["rat_s"] = nw_result_pf_branch["rat_s"]
                        nw_result_branch["raf_s"] = nw_result_pf_branch["raf_s"]
                        end
                    end
                end
            end

            # add bus data
            result["solution"]["nw"][nw]["bus"] = Dict{String, Dict}()
            for (bus_id, nw_result_pf_bus) in nw_result_pf["bus"]
                result["solution"]["nw"][nw]["bus"][bus_id] = Dict{String, Any}()
                # add voltage angles
                if haskey(nw_result_pf_bus, "va")
                    result["solution"]["nw"][nw]["bus"][bus_id]["va"] = nw_result_pf_bus["va"]
                end
            end

            # add gen data TODO this needs to be replaced by a proper re-computation
            result["solution"]["nw"][nw]["gen"] = nw_result_pf["gen"]
            
            # add Benders variables - set them to zero, if we set them to the objective, the objective function value of the MIPStart will be double, because the pf (and thus rate a) is used, AND the eta
            #if ! haskey(result["solution"]["nw"][nw], "eta")
            #    result["solution"]["nw"][nw]["eta"] = Dict{String, Float64}()
            #end
            #result["solution"]["nw"][nw]["eta"]["1"] = Dict("η" => 0)#result_pf["objective"]) 
        end

        # reset contingencies
        data["nw"]["1"]["contingencies"] = contingencies_copy
    else
        result["solution"] = result_pf["solution"]
    end
end

function check_pg_equals_pgsec(result)
    gen_pgs = Dict{Int64, Set}()
    for (gen_id, gen_data) in result["solution"]["gen"]
      gen_id = parse.(Int, split(replace(gen_id, r"[()]" => ""), ','))[1]
      if length(gen_id) == 2
        pg = gen_data["pg_sec"]
      else
        pg = gen_data["pg"]
      end
      if gen_id ∉ keys(gen_pgs)
        gen_pgs[gen_id] = Set()
      end
      push!(gen_pgs[gen_id], pg)
    end
    println(gen_pgs)
end

function show_system_balance(data)
    if !PARAMS[:PRINTING][:SYS_BALANCE] && !PARAMS[:PRINTING][:CONTINGENCIES] && !PARAMS[:PRINTING][:MONITORING]; return; end
    println()
    for (nw_id, nw_data) in data["nw"]
        if PARAMS[:PRINTING][:SYS_BALANCE]
            total_pd = 0
            #println("Loads:")
            for (l_id, load) in nw_data["load"]
                if load["status"] == 0 && load["pd"] > PARAMS[:PRECISION][:POWER_PREC]
                    @warn("Nonnegative deactivated load")
                end
                #if load["pd"] > PARAMS[:PRECISION][:POWER_PREC]
                #    println("   $(get_bus_by_index(nw_data, load["load_bus"])["name"]): $(round(load["pd"], digits=2))")
                #end
                total_pd += load["pd"]
            end
            total_pg = 0
            total_pg_max = 0
            #println("Gens:")
            for (g_id, gen) in nw_data["gen"]
                if gen["gen_status"] == 0 && gen["pg"] > PARAMS[:PRECISION][:POWER_PREC]
                    @warn("Nonnegative deactivated generator")
                end
                total_pg += gen["pg"]
                total_pg_max += gen["pmax"]
                #if gen["pg"] > PARAMS[:PRECISION][:POWER_PREC]
                #    println("   $(gen["name"]): $(round(gen["pg"], digits=2))")
                #end
            end
            println("  Total PD: $(round(total_pd, digits=SOLVER_PRECISION_DIGITS))\n",
                    "  Total PG: $(round(total_pg, digits=SOLVER_PRECISION_DIGITS))\n",
                    "  Diff.   : $(round(total_pg-total_pd, digits=SOLVER_PRECISION_DIGITS))\n",
                    "  Max   PG: $(round(total_pg_max, digits=SOLVER_PRECISION_DIGITS))")
            println()
        end
        
        if PARAMS[:PRINTING][:CONTINGENCIES]
            # show contingencies
            @info("Considering these contingencies:")
            for (c, contingency) in enumerate(nw_data["contingencies"])
                println("  Contingency $(lpad(c, 3, " ")) $(get_contingency_names(data, contingency)) ($(contingency[:branch]))")
            end
        end

        if PARAMS[:PRINTING][:MONITORING]
            # show monitoring branches
            if haskey(nw_data, "monitoring")
                @info("Considering these branches for monitoring:")
                for (c, contingency) in enumerate(nw_data["monitoring"])
                    println("  Montoring $(lpad(c, 3, " ")) $(get_contingency_names(data, contingency)) ($(contingency[:branch]))")
                    #println("    Type: $(nw_data["branch"][string(contingency[:branch][1])]["type"])")
                end
            end
        end
    end
    println()
end

function set_default_branch_attributes!(data)
    for (b_id, branch) in data["branch"]
        if ! haskey(branch, "status")
            branch["status"] = branch["br_status"]
            branch["br_status"] = 1
        end
        if ! haskey(branch, "switchable")
            branch["switchable"] = true
        end
        if ! haskey(branch, "shiftable")
            branch["shiftable"] = false
        end
        if ! haskey(branch, "type")
            branch["type"] = "L"
        end
        if ! haskey(branch, "monitored")
            branch["monitored"] = true
        end
        if ! haskey(branch, "electrical")
            branch["electrical"] = true
        end
        if PARAMS[:BIG_M_MULT][:SET_MAX_ANGLE_CLOSE]
            branch["angmax_close"] = PARAMS[:BIG_M_MULT][:MAX_ANG_OFF] 
            branch["angmin_close"] = -PARAMS[:BIG_M_MULT][:MAX_ANG_OFF]
        end
    end
end

function set_default_gen_attributes!(data)
    for (g_id, gen) in data["gen"]

        redispatch_model  = PARAMS[:LIMITS][:GENERATOR][:REDISPATCH_MODEL]
        gen_switching     = PARAMS[:LIMITS][:GENERATOR][:GEN_SWITCHING]
        override          = PARAMS[:LIMITS][:GENERATOR][:GEN_REDISPATCH_OVERRIDE]

        # pg_start: always records the original pg from data (before any modification)
        gen["pg_start"] = haskey(gen, "pg_0") ? gen["pg_0"] : gen["pg"]

        if !redispatch_model
            # --- No redispatch model: treat every generator as starting from zero ---
            # pg_0 = 0, RDP- = 0, RDP+ = pmax (ignores any pg_0 that may exist in data)
            gen["pg_0"] = 0.0
            gen["pg"]   = 0.0
            gen["RDP-"] = 0.0
            gen["RDP+"] = gen["pmax"]
        else
            # --- Redispatch model: use pg_0 from data if available, else pg ---
            if !haskey(gen, "pg_0")
                gen["pg_0"] = gen["pg"]
            end

            if gen_switching
                if override
                    # Generator can be switched and override is on:
                    # pmin is enforced elsewhere, so reduce it to zero here.
                    # RDP- spans down to 0, RDP+ spans up to pmax.
                    gen["pmin"] = 0.0
                    gen["RDP-"] = haskey(gen, "RDP-") ? gen["RDP-"] : gen["pg_0"]
                    gen["RDP+"] = haskey(gen, "RDP+") ? gen["RDP+"] : gen["pmax"] - gen["pg_0"]
                else
                    # Generator can be switched but override is off:
                    # pmin is enforced, so pg and pg_0 are clamped to pmin.
                    # RDP- = 0 (can't go below committed pmin), RDP+ = pmax - pmin.
                    # Warn that granting pmin for free distorts costs.
                    if redispatch_model
                        Memento.warn(_LOGGER, "GEN_SWITCHING=true, GEN_REDISPATCH_OVERRIDE=false, REDISPATCH_MODEL=true: " *
                            "pg_0 is set to pmin for gen $(g_id), which grants pmin dispatch for free and distorts marginal costs. " *
                            "Consider enabling GEN_REDISPATCH_OVERRIDE to avoid this.")
                    end
                    gen["pg_0"] = gen["pmin"]
                    gen["pg"]   = gen["pmin"]
                    gen["RDP-"] = 0.0
                    gen["RDP+"] = haskey(gen, "RDP+") ? gen["RDP+"] : gen["pmax"] - gen["pg_0"]
                end
            else
                # Generator switching is off:
                # pmin is enforced elsewhere (or override handles it); reduce pmin to zero.
                gen["pmin"] = 0.0
                gen["pg"]   = 0.0
                gen["RDP-"] = haskey(gen, "RDP-") ? gen["RDP-"] : gen["pg_0"]
                gen["RDP+"] = haskey(gen, "RDP+") ? gen["RDP+"] : gen["pmax"] - gen["pg_0"]
            end
        end

        # ramps should not be binding, if they did not exist before
        if ! haskey(gen, "ramp_60_pos")
            gen["ramp_60_pos"] = gen["RDP+"]
        end

        if ! haskey(gen, "ramp_60_neg")
            gen["ramp_60_neg"] = gen["RDP-"]
        end

        # the cost model is either piecewise linear (1) or polynomial (2) 
        # both are converted to linear, giving gen["cost_produce"] = c::Float64
        if ! haskey(gen, "cost_produce") && gen["pmax"] > 0
            cost = gen["cost"]
            model = gen["model"]
            ncost = gen["ncost"]
            if model == 1
                # piecewise linear cost function: cost = c2 * pg + c1 * z + c0 * y, where z is the auxiliary variable for the piecewise linear cost function and y is the binary variable for the piecewise linear cost function
                # we simply take the slope of the last segment as the marginal cost, which is a common approach in the literature
                gen["cost_produce"] = cost[end][2]
            elseif model == 2
                # polynomial cost function: cost = c_n * pg^n + c_(n-1) * pg^(n-1) + ... + c_1 * pg + c_0
                # we take the derivative of the cost function at pg_0 as the marginal cost, which is a common approach in the literature
                # derivative: d/dpg = n * c_n * pg^(n-1) + (n-1) * c_(n-1) * pg^(n-2) + ... + c_1
                pg_0 = gen["pg_0"]
                marginal_cost = 0.0
                for i in 1:(ncost-1)
                    # cost[i] is the coefficient for pg^(ncost-i)
                    power = ncost - i
                    marginal_cost += power * cost[i] * (pg_0^(power - 1))
                end
                gen["cost_produce"] = marginal_cost
            end
            @assert(gen["cost_produce"] >= 0, "Negative production cost for generator $(gen["name"]): $(gen["cost_produce"])")
        end

        # set changing cost: set to 0 to make redispatch basically free, set to > cost_produce if pg_0 is already not 0 (problem: we get more back than we paid if we produce less)
        if ! haskey(gen, "cost_change") && gen["pmax"] > 0
            gen["cost_change"] = gen["cost_produce"] * 1
        end


        # set status and, for now since we do not have "previous timestep" data, status_previous
        if get(gen, "pg_0", gen["pg"]) > 0
            @assert(gen["gen_status"] == 1)
            gen["status"] = 1
        else
            gen["status"] = 0
        end
        gen["status_previous"]  = gen["status"]
        # gen_status,  >  0 - machine in service
        #              <= 0 - machine out of service
        gen["gen_status"]       = 1
        gen["switchable"]       = true

        # if a gen is said to be redispatchable to an invalid area ABOVE 0 (for shutdown), fix pmin/RDP-/switchable
        if gen["pg_0"] - gen["RDP-"] > 0
            gen["switchable"] = false
            gen["RDP-"] = gen["pg_0"] - gen["pmin"]
        end
    end
end

"does all the steps to expand the graph with busbars"
function add_busbars(data)
    @info("Using busbars.")
    busbar_warn = true

    data_bus = data["bus"]::Dict{String, Any}

    if PARAMS[:SUBSTATION_RECONF][:BRANCH_RECONF]::Bool
        # Define the possible values and their probabilities
        busbars = [1, 2, 3, 4]
        probs = [0.9, 0.05, 0.03, 0.02]

        # Create a categorical distribution
        Random.seed!(1234)
        d = Distributions.Categorical(probs)

        # add a name if necessary
        for (bus_id, bus::Dict{String, Any}) in data_bus
            get!(bus, "name", bus_id)
            get!(bus, "n_busbars", busbars[rand(d)])
        end

        add_busbars_buses(data_bus)
    else
        # add a name if necessary
        for (bus_id, bus::Dict{String, Any}) in data_bus
            get!(bus, "name", bus_id)
        end

        @warn("busbar mode enabled, but BRANCH_RECONF==false. n_busbars will be ignored!")
    end

    data["substations"] = _detect_substations(data_bus, get_substation_name_filter(data))
    _rewire_graph_with_substations!(data)
    check_substation_integrity!(data["branch"]::Dict{String, Any}, data["substations"]::Dict{String, Dict{String, Vector{Any}}})
end

"add time steps, based on the vector `casefiles`. The first casefile is already in `data`"
function add_time_steps(data, casefiles, linkingfile, timesteps, slack_types, busbars)
    start_time = data["start_time"]
    data = replicate(data, timesteps)
    data["slack_types"] = slack_types
    for n in 1:timesteps
        data["nw"]["$n"] = parse_file(casefiles[n])
        t_data = data["nw"]["$n"]::Dict{String, Any}
        t_data["file"] = casefiles[n]

        # add timestep-independent data
        for key in _mp_non_export_columns
            if haskey(data, key)
                t_data[key] = data[key]
            end
        end
        
        t_data["per_unit"]  = data["per_unit"]
        t_data["start_time"] = start_time
        t_data["time"]      = n

        #add_missing_attributes(t_data, PARAMS[:TECHNICAL][:MISSING_LIST])
        set_default_branch_attributes!(t_data)
        set_default_gen_attributes!(t_data)
        
        if PARAMS[:PREPROCESSING][:REMOVE_LOOSE_ENDS]
            # multiple iteration, since a single removal might yield a new loose end
            for iter in 1:PARAMS[:PREPROCESSING][:LOOSE_ENDS_ITERS]::Int64
                set_loose_ends_to_non_electrical!(t_data)
                merge_buses!(t_data["branch"], t_data["bus"], t_data["gen"], t_data["load"], t_data["substations"])
            end
        end

        # fix non-switchable branches and calculate contingencies
        calculate_contingencies!(t_data, merge_substations=busbars)
        
        # add busbars
        if busbars
            add_busbars(t_data)
            t_data["BUSBARS"] = true
        else
            t_data["BUSBARS"] = false
        end
        
        make_bus_indices_consecutive!(t_data, t_data["branch"], t_data["bus"], t_data["gen"], t_data["load"])

        if false && PARAMS[:SUBSTATION_RECONF][:BRANCH_RECONF]
            # get dir of data["file"]
            file = t_data["file"]::String
            dir = dirname(file)
            # get file name without extension
            filename = basename(file)
            @info("Exporting busbars to file $(joinpath(dir, "$(filename)_bb.m"))")
            export_matpower(joinpath(dir, "$(filename)_bb.m"), copy(t_data))
        end
    end

    check_gen_consistency!(data)

    if timesteps > 1
        apply_linking!(data, linkingfile)
    end

    set_dc_pf_start_values!(data)
    data["multinetwork"] = true
    return data
end

"read in linkingfile and add ramping constraints to the data"
function apply_linking!(data, linkingfile::String)
    linking = CSV.read(linkingfile, DataFrame)
    for (gen_id, ramp_up, ramp_down) in eachrow(linking)
        for i in eachindex(data["nw"])
            data["nw"][string(i)]["gen"][string(gen_id)]["ramp_60_pos"] = ramp_up
            data["nw"][string(i)]["gen"][string(gen_id)]["ramp_60_neg"] = ramp_down
        end
    end
end


## All other getters ##
"gets branch voltage from voltage of connected buses. asserts bus voltages to not be different."
function get_branch_voltage(data, branch)
    kv1 = -1
    kv2 = -1
    for (bus_id, bus1) in data["bus"]
        if bus1["index"] == branch["f_bus"]
            kv1 = bus1["base_kv"]
        end
    end
    for (bus_id, bus2) in data["bus"]
        if bus2["index"] == branch["t_bus"]
            kv2 = bus2["base_kv"]
        end
    end
    @assert(kv1 ≠ -1 && kv2 ≠ -1)
    @assert(kv1 == kv2, "Branch $(branch["name"]) ($(branch["index"])) of type $(branch["type"]) has different voltages at its ends: $kv1 and $kv2")
    @assert(kv1 ∈ [27, 110, 220, 380])
    return kv1
end

"get bus by its name. returns the bus data, or if not found, an empty dict."
function get_bus_by_name(data, bus_name; include_busbar=true)::Dict
    chars = include_busbar ? 8 : 7
    for (bus_i, bus) in data["bus"]
        if bus["name"][1:chars] == bus_name
            return bus
        end
    end
    return Dict()
end

"get branch by its name. returns the branch data, or if not found, an empty dict"
function get_branch_by_name(data, branch_name; include_busbars=true)::Dict
    if include_busbars
        if haskey(data["branchname_to_id"], branch_name)
            return data["branch"][data["branchname_to_id"]]
        end
    else
        #println("  D7BUER1B_D7YWDL1O_WONNEG_O:line")
        for (b_name, b_id) in data["branchname_to_id"]
            b_name_no_busbars = "$(b_name[1:7])_$(b_name[10:16])_$(b_name[19:end])"
            if branch_name == b_name_no_busbars
                return data["branch"][b_id]
            end
        end
    end
    return Dict()
end

"get load by its host bus id. returns the load data, or if not found, an empty dict"
function get_load_by_host_bus(data, source_id)::Dict
    for (load_i, load) in data["load"]
        if load["source_id"][2] == source_id
            return load
        end
    end
    return Dict()
end

"get branch by its fbus and tbus. returns the branch data, or if not found, an empty dict"
function get_branch_by_fbus_tbus(data, fbus::Int64, tbus::Int64)::Dict
    for (branch_id, branch) in data["branch"]
        if (branch["f_bus"] == fbus && branch["t_bus"] == tbus) || (branch["f_bus"] == tbus && branch["t_bus"] == fbus)
          return branch
        end
      end
    return Dict()
end

function get_contingency_names(data, contingency)
    names = Vector{String}()
    for (c_type, c_ids) in contingency
        if c_type == :branch
            for c_id in c_ids
                branch = data["nw"]["1"]["branch"]["$c_id"]
                if haskey(branch, "branch_name_intern")
                    push!(names, branch["branch_name_intern"])
                #(f_bus, t_bus) = get_buses_of_branch(data, branch)
                else
                    push!(names, "$(get(branch, "name", ""))")
                end
            end
        end
        if c_type == :gen && length(c_ids) > 0
            @warn("Not implemented")
        end
    end
    return names
end

"get total active load at a bus"
function get_total_bus_load_p(data, bus)::Float64
    total_load = 0.0
    for (_, load) in data["load"]
        if load["source_id"][2] == bus["index"]
            total_load += load["pd"]
        end
    end
    return total_load
end

"get total reactive load at a data"
function get_total_bus_load_q(data, bus)::Float64
    total_load = 0.0
    for (_, load) in data["load"]
        if load["source_id"][2] == bus["index"]
            total_load += load["qd"]
        end
    end
    return total_load
end

"get total active generation at a bus"
function get_total_bus_gen_p(data, bus)::Float64
    total_gen = 0.0
    for (_, gen) in data["gen"]
        if gen["source_id"][2] == bus["index"]
            total_gen += gen["pg_new"]
        end
    end
    return total_gen
end

"get total reactive generation at a bus"
function get_total_bus_gen_q(data, bus)::Float64
    total_gen = 0.0
    for (_, gen) in data["gen"]
        if gen["source_id"][2] == bus["index"]
            total_gen += get(gen, "qg_new", 0.0)
        end
    end
    return total_gen
end

function write_solving_data(casename, pm_data, line_switches::Dict{Int64, Dict{String, Pair{Int, Int}}}, gen_changes::Dict{Int64, Dict{String, Vector}}, results_dir::String=CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:RESULTS_DIR])
    run_id      = "$(casename)_$(Dates.now())"
    scenario_id = casename
    @info("Writing solving data for case $casename to $results_dir")
    mkpath(results_dir)
  
    # RA_TOPO
    column_order = ["Run ID", "Scenario ID", "Timestamp", "RA ID", "RA Name", "TSO",
      "Line description mode", "UCT Node From", "UCT Node To", "Element Name",
      "Status before", "Status after", "Busbar before", "Busbar after", "Penalty Cost"]
    n_topo = 1
    for (time_step, nw_line_switches) in line_switches
      for (line_name, (status_before, status_after)) in nw_line_switches
        data = OrderedDict{String, Any}()
  
        #line_data = pm_data["nw"][string(time_step)]["branch"][string(line_name)]
        region = length(line_name) ≥ 2 ? line_name[1:2] : ""
  
        data["Run ID"]                = run_id
        data["Scenario ID"]           = scenario_id
        data["Timestamp"]             = "$(lpad(time_step-1, 2, "0")):30"
        data["RA ID"]                 = "$(region)_$(lpad(n_topo, 5, "0"))"
        data["RA Name"]               = ""
        data["TSO"]                   = region
        data["Line description mode"] = "ElementName"
        if length(line_name) < 19
            line_name = "UNKNOWN1_UNKNOWN2_L$line_name"
        end
        data["UCT Node From"]         = line_name[1:8]
        data["UCT Node To"]           = line_name[10:17]
        data["Element Name"]          = line_name[19:end]
        data["Status before"]         = status_before
        data["Status after"]          = status_after
        data["Busbar before"]         = ""
        data["Busbar after"]          = ""
        data["Penalty Cost"]          = 0
  
        transposed_data = Dict(key => [value] for (key, value) in data)
        if n_topo == 1
          CSV.write(joinpath(results_dir, "RA_TOPO_$(scenario_id).csv"), transposed_data, header=column_order)
        else
          CSV.write(joinpath(results_dir, "RA_TOPO_$(scenario_id).csv"), transposed_data, append=true, header=column_order)
        end
        n_topo += 1
      end
    end
    
    # RA_PST
    # TBD
    
    # RA_RD
    column_order = ["Run ID", "Scenario ID", "Timestamp", "RA RD ID", "TSO",
      "Generator Name", "RD description mode", 
      "P before Optimization", "P after Optimization", "P Redispatch",
      "Startup Cost", "Redispatch Cost", "Penalty Cost"]
    n_rd = 1
    for (time_step, nw_gen_changes) in gen_changes
      for (gen_name, (gen_id, pg_before, pg_after)) in nw_gen_changes
        data = OrderedDict{String, Any}()
        
        time        = pm_data["nw"][string(time_step)]["time"]
        gen_data    = pm_data["nw"][string(time_step)]["gen"][string(trunc(Int, gen_id))]
        region      = get(gen_data, "name", "R1")[1:2]
  
        data["Run ID"]                = run_id
        data["Scenario ID"]           = scenario_id
        data["Timestamp"]             = time
        data["RA RD ID"]              = get(gen_data, "ra_rd_id", -1)
        data["TSO"]                   = region
        data["Generator Name"]        = gen_name
        data["RD description mode"]   = occursin("GSK", gen_name) ? "GSK" : "NODE"
        # ...
        data["P before Optimization"] = pg_before
        data["P after Optimization"]  = pg_after
        data["P Redispatch"]          = pg_after - pg_before
        data["Startup Cost"]          = 0
        data["Redispatch Cost"]       = 0
        data["Penalty Cost"]          = 0
        
        transposed_data = OrderedDict(key => [value] for (key, value) in data)
        if n_rd == 1
          CSV.write(joinpath(results_dir, "RA_RD_$(scenario_id).csv"), transposed_data, header=column_order)
        else
          CSV.write(joinpath(results_dir, "RA_RD_$(scenario_id).csv"), transposed_data, append=true, header=column_order)
        end
        n_rd += 1
      end
    end
end

"given a multiplier, adjust rate a to new value"
function adjust_ratea!(t_data; ub_mult=PARAMS[:LIMIT_MULT][:RATEA_MULT])
    for (b_id, branch) in t_data["branch"]
        branch["rate_a"] *= ub_mult
    end
end

function reduce_max_loading!(t_data; ub_mult=PARAMS[:LIMIT_MULT][:RATEA_MULT])
    for (b_id, branch) in t_data["branch"]
        branch["max_loading_n1"] = get(branch, "max_loading_n1", 1.0) * ub_mult
        branch["max_loading"] = get(branch, "max_loading", 1.0) * ub_mult
    end
end

function reset_max_loading!(t_data; ub_mult=PARAMS[:LIMIT_MULT][:RATEA_MULT])
    for (b_id, branch) in t_data["branch"]
        branch["max_loading_n1"] = branch["max_loading_n1"] / ub_mult
        branch["max_loading"] = branch["max_loading"] / ub_mult
    end
end

"add dictionaries 'info', 'time' and 'callback' to result dictionary. 'time' contains callback times and 'callback' exists only if optional args 'callback_stats' and 'callback_times' are provided"
function add_solving_data_to_results!(data, result, casename::String; callback_stats::Dict=Dict(), callback_times::Dict=Dict())
    # add general info
    result["info"] = Dict()
    result["info"]["instance"] = casename
    result["info"]["status"] = result["termination_status"]
    result["info"]["objective"] = result["objective"]
    result["info"]["n_switching_actions"] = get(result, "n_switching_actions", 0)

    # add time info
    t_tmp = result["time"] # this is the time for building the model and solving it
    result["time"] = Dict()
    result["time"]["total"] = result["total_time"]
    result["time"]["build"] = t_tmp - result["solve_time"]
    result["time"]["solve"] = result["solve_time"]
    
    if ! isempty(callback_times)
        # total time for subproblems (preprocessing, N-1, connectedness)
        result["time"]["SP_conn"] = round(callback_times[:connect], digits=4)
        result["time"]["SP_N1"] = round(callback_times[:contingency], digits=4)
        result["time"]["SP_pre"] = round(callback_times[:pre], digits=4)
        result["time"]["SP"] = result["time"]["SP_pre"] + result["time"]["SP_N1"] + result["time"]["SP_conn"]

        # total time for master problem
        result["time"]["MP"] = result["time"]["solve"] - result["time"]["SP"]

        # more detailed N-1 times, detailed by action type (if parallelized, these are CPU time, not wall clock time!)
        result["time"]["SP_N1_primal"] = round(callback_times[:contingency_primal], digits=4)
        result["time"]["SP_N1_dual"] = round(callback_times[:contingency_dual], digits=4)
        result["time"]["SP_N1_primal_check"] = round(callback_times[:contingency_primal_check], digits=4)
        result["time"]["SP_N1_dual_check"] = round(callback_times[:contingency_dual_check], digits=4)
        result["time"]["SP_N1_cut_build"] = round(callback_times[:contingency_cut_build], digits=4)
        result["time"]["SP_N1_cut_check"] = round(callback_times[:contingency_cut_check], digits=4)
        result["time"]["SP_N1_fallback"] = round(callback_times[:contingency_fallback], digits=4)

        if PARAMS[:SOLVER][:MULTITHREADING][:N_THREADS_BENDERS] == 1
            result["time"]["SP_N1_misc"] = result["time"]["SP_N1"] - result["time"]["SP_N1_primal"] - result["time"]["SP_N1_dual"] - result["time"]["SP_N1_primal_check"] - result["time"]["SP_N1_dual_check"] - result["time"]["SP_N1_cut_build"] - result["time"]["SP_N1_cut_check"] - result["time"]["SP_N1_fallback"]
        else
            result["time"]["SP_N1_misc"] = NaN
        end

        # more detailed N-1 times, detailed by cut type (if parallelized, these are CPU time, not wall clock time!) 
        result["time"]["SP_N1_cut_opt"] = round(callback_times[:contingency_cut_opt], digits=4)
        result["time"]["SP_N1_cut_feas"] = round(callback_times[:contingency_cut_feas], digits=4)
        result["time"]["SP_N1_cut_none"] = round(callback_times[:contingency_cut_none], digits=4)

        result["time"]["SP_N1_refactorization"] = round(callback_times[:refactorization], digits=4)

    end

    if PARAMS[:POSTPROCESSING][:AC_FEAS][:CHECK] && result["info"]["status"] == JuMP.MOI.OPTIMAL
        result["time"]["ac_feas_check"] = 0
        ac_feas_result = Dict(
            "solved" => true,
            "objective" => 0.0,
            "angle_violations_sum" => 0.0,
            "rate_violations_sum" => 0.0
        )
        for nw_id in keys(result["solution"]["nw"])
            result["time"]["ac_feas_check"] += result["solution"]["nw"][nw_id]["ac_feasibility"]["run_time"]
            ac_feas_result["solved"] = ac_feas_result["solved"] && result["solution"]["nw"][nw_id]["ac_feasibility"]["solved"]
            ac_feas_result["objective"] += result["solution"]["nw"][nw_id]["ac_feasibility"]["objective"]
            ac_feas_result["angle_violations_sum"] += result["solution"]["nw"][nw_id]["ac_feasibility"]["angle_violations_sum"]
            ac_feas_result["rate_violations_sum"] += result["solution"]["nw"][nw_id]["ac_feasibility"]["rate_violations_sum"]
        end
        result["ac_feas"] = ac_feas_result
    end

    if PARAMS[:WARMSTARTING][:ENABLED] && haskey(data, "warmstart_time")
        result["time"]["warmstart"] = data["warmstart_time"]
        delete!(data, "warmstart_time")
    end

    if PARAMS[:WARMSTARTING][:ENABLED]
        result["info"]["n_carried_benders_cuts"] = get(data, "n_carried_benders_cuts", 0)
    end

    if PARAMS[:LOSS_ESTIMATION][:ESTIMATE_LOSSES]
        result["time"]["loss_estimation_time"] = data["loss_estimation_time"]
        delete!(result, "loss_estimation_time")
    end

    # add callback info
    result["callback"] = Dict()
    if ! isempty(callback_stats)
        result["callback"]["SP_optimizer_type"] = string(PARAMS[:BENDERS][:SUBPROBLEM_SOLVER_TYPE]) 
        result["callback"]["n_callbacks"] = callback_stats[:n_callbacks]
        result["callback"]["n_solved_subs_total"] = callback_stats[:n_solved_subs_total]

        # all found Benders cuts
        result["callback"]["n_benders_feas_cuts_found"] = callback_stats[:n_benders_feas_cuts_found]
        result["callback"]["n_benders_opt_cuts_found"] = callback_stats[:n_benders_opt_cuts_found]
        result["callback"]["n_benders_cuts_found"] = callback_stats[:n_benders_cuts_found]

        # all added Benders cuts (disagrees if cut selection takes place)
        #result["callback"]["n_benders_feas_cuts_added"] = callback_stats[:n_benders_feas_cuts_added] # not yet tracked separately
        #result["callback"]["n_benders_opt_cuts_added"] = callback_stats[:n_benders_opt_cuts_added] # not yet tracked separately
        result["callback"]["n_benders_cuts_added"] = callback_stats[:n_benders_cuts_added]

        # all cuts
        result["callback"]["n_benders_cut_value_disagree"] = callback_stats[:n_cut_value_disagree]
        result["callback"]["n_benders_cut_value_agree"] = callback_stats[:n_cut_value_agree]
        result["callback"]["n_benders_cut_value_fail"] = callback_stats[:n_cut_value_fail]
        result["callback"]["n_connect_cuts_added"] = callback_stats[:n_connect_cuts_added]

        # all solver runs and failures
        result["callback"]["n_leq_fails"] = callback_stats[:n_leq_fails]
        result["callback"]["n_leq_residual_failures"] = callback_stats[:n_leq_residual_failures]
        result["callback"]["n_primal_lp"] = callback_stats[:n_primal_lp]
        result["callback"]["n_fallback_primal_lp"] = callback_stats[:n_fallback_primal_lp]
        result["callback"]["n_dual_solver"] = callback_stats[:n_dual_solver]
        result["callback"]["n_dual_lp"] = callback_stats[:n_dual_lp]
        result["callback"]["n_questionable_sol"]            = callback_stats[:n_questionable_sol]

        result["callback"]["n_rf_new"] = callback_stats[:n_rf_new]
        result["callback"]["n_rf_update"] = callback_stats[:n_rf_update]
        if result["callback"]["n_questionable_sol"] > 0.3 * result["callback"]["n_solved_subs_total"]
            Memento.warn(_LOGGER, "Easy SP failed often: $(result["callback"]["n_questionable_sol"]) fails along $(result["callback"]["n_solved_subs_total"]) subproblems ($(result["callback"]["n_questionable_sol"] / result["callback"]["n_solved_subs_total"] * 100) %). Consider using a different solver.")
        end
    end
end

"set electrical=false/monitored=false/switchable=false (and remove from contingencies) for all branches that, at at least one end, have a bus that has no branches, gens or loads connected to it"
function set_loose_ends_to_non_electrical!(t_data)
    #throw(InvalidStateException("Untested feature: set_loose_ends_to_non_electrical!()", :untested))
    for (b_id, branch) ∈ t_data["branch"]
        for end_bus ∈ [branch["t_bus"], branch["f_bus"]]
            # check if one of the end buses has no branches, gens or loads
            end_bus_loads = get_connected_loads(t_data, end_bus)
            end_bus_gens  = get_connected_gens(t_data, end_bus)
            end_bus_branches::Vector{Int64} = get_connected_branches(t_data, end_bus)
            
            # remove this branch itself from the list of branches
            end_bus_branches = filter(b -> b != parse(Int64, b_id), end_bus_branches)

            # if nothing connected at this end, set branch to non-electrical, non-monitored, non-switchable and remove it from contingencies if applicable. break since only one bus has to be a loose end
            if isempty(end_bus_loads) && isempty(end_bus_gens) && isempty(end_bus_branches)
                # set branch to non-electrical
                t_data["branch"][b_id]["electrical"] = false
                # set branch to non-monitored
                t_data["branch"][b_id]["monitored"] = false
                # set branch to non-switchable
                t_data["branch"][b_id]["switchable"] = false
                # set branch to not in contingencies
                t_data["branch"][b_id]["in_contingencies"] = false

                # remove from contingencies if applicable
                #for contingency in t_data["contingencies"]
                #    if b_id ∈ contingency[:branch]
                #        delete!(t_data["contingencies"], contingency)
                #    end
                #end
                break
            end
        end
    end
end

"check sum over all data[objective_costs] (except overload_all) against result[objective], comparing base case cost"
function check_obj(data_sol_cost, result_obj; tol = PARAMS[:POSTPROCESSING][:CHECKS][:BASE_COST_CHECK_TOL], base_only::Bool=true)
    # check sum over all data[objective_costs] (except overload_all) against result[objective]
    remove_which_overload = base_only ? "overload_all" : "overload_base_case"
    objective_sum = sum(data_sol_cost[key] for key in keys(data_sol_cost) if key != remove_which_overload)
    if ! isapprox(objective_sum, result_obj, atol=tol)
        str = @sprintf("[Base cost check] Objective components do not match: obj_check = %12.2f ≠ %12.2f = obj (absolute errror: %12.2f; relative error: %6.2f%%)", objective_sum, result_obj, abs(objective_sum - result_obj), abs(objective_sum - result_obj) / max(abs(objective_sum), abs(result_obj)) * 100)
        @error("Objective components do not match: $(str)")
        display(data_sol_cost)
    end
end

"
Removes all breaker branches and busses from given network.
Reconnects the remaining components soundly by creating new branches that connect the buses that were formerly connected through the breaker.
Returns a mapping of the removed (non-breaker) branches to the new branches and corresponding breaker branches.
"
function contract_breakers!(data)

    # maps the connecting branch to the replacement branches and the corresponding breaker branches
    connecting_branches_mapping = Dict{String, Vector{Tuple{String, String}}}()

    if !haskey(data, "substations")
        return connecting_branches_mapping
    end

    for (_, substation) in data["substations"]
        for breaker_bus in substation["aux_bus"]
            breaker_branches = breaker_bus["breaker_branches"]

            # get the buses connected to the breaker with breaker branches
            connecting_buses::Vector{Int} = []
            # mapping the connecting buses to the breaker branches that connect them with the breaker
            connecting_bus_breaker_mapping = Dict{Int, String}()
            for branch_id in breaker_branches
                branch = data["branch"][branch_id]
                if branch["f_bus"] == breaker_bus["bus_i"]
                    push!(connecting_buses, branch["t_bus"])
                    connecting_bus_breaker_mapping[branch["t_bus"]] = branch_id
                else
                    push!(connecting_buses, branch["f_bus"])
                    connecting_bus_breaker_mapping[branch["f_bus"]] = branch_id
                end
            end

            # if connecting buses is empty, the breaker bus is not connected (continue with the next breaker)
            if isempty(connecting_buses)
                continue
            end

            # get the non-breaker branches connecting the breaker with other buses
            connecting_branches::Vector{String} = []
            for branch_id in keys(data["branch"])
                if !(branch_id in breaker_branches) && (data["branch"][branch_id]["f_bus"] == breaker_bus["bus_i"] || data["branch"][branch_id]["t_bus"] == breaker_bus["bus_i"])
                    push!(connecting_branches, branch_id)
                end
            end

            # create new branches for each connecting bus and each connecting branch, connecting the busses that were formerly connected through the breaker
            index = maximum([parse(Int, key) for key in keys(data["branch"])]) + 1
            for connecting_branch in connecting_branches
                branch = data["branch"][connecting_branch]
                breaker_bus_dir = branch["f_bus"] == breaker_bus["bus_i"] ? "f_bus" : "t_bus"
                connecting_branches_mapping[connecting_branch] = []
                for connecting_bus in connecting_buses
                    new_branch = deepcopy(branch)
                    new_branch["index"] = index
                    new_branch[breaker_bus_dir] = connecting_bus
                    push!(data["branch"], string(index) => new_branch)
                    push!(connecting_branches_mapping[connecting_branch], (string(index), connecting_bus_breaker_mapping[connecting_bus]))
                    index += 1
                end
                delete!(data["branch"], connecting_branch)
            end

            # remove the breaker branches and the breaker bus
            for breaker_branch in breaker_branches
                delete!(data["branch"], breaker_branch)
            end
            delete!(substation, "aux_bus")
            delete!(data["bus"], breaker_bus["bus_i"])
        end
    end
    correct_bus_types!(data)
    return connecting_branches_mapping
end

"create directories for logs, results, and solver logs"
function create_dirs(dirs::Dict{Symbol, Any})
    mkpath(dirs[:LOG_DIR])
    mkpath(dirs[:RESULTS_DIR])
    mkpath(dirs[:SOLVER_LOG_DIR])
end

"main function to solve an application with the given data and solver"
function solve_application!(data::Dict{String, Any}, 
                            application::String, 
                            slack_types::Vector{Symbol},
                            timesteps::Int64, 
                            nested::Bool, 
                            solver_convex::JuMP.MOI.OptimizerWithAttributes,
                            solver_nlp::JuMP.MOI.OptimizerWithAttributes
                           )::Dict{String, Any}

    # start optimization
    result = []
    t = @elapsed begin
    if application == "DC_LODF"
        result = solve_pf(data, DCPPowerModel, solver_convex)
        PTDF = compute_PTDF(data)
        LODF = compute_LODF(PTDF, data)
        @show(LODF)
    elseif application == "DC_OPF"
        result = solve_opf(data, DCPPowerModel, solver_convex)
    elseif application == "AC_OPF"
        result = solve_opf(data["nw"]["1"], ACRPowerModel, solver_convex)
    elseif application == "AC_OPF_TIME_SLACK" || application == "AC_SCOPF_TIME_SLACK"
        _data = data["nw"]["1"]
        _data["per_unit"] = data["per_unit"]
        _data["cost_coeffs"] = data["cost_coeffs"]
        _data["application"] = application
        _data["multinetwork"] = false
            
        
        if application == "AC_OPF_TIME_SLACK"
            limits_generator, tal_timestep = override_ac_feas_params()
            slack_types = _data["slack_types"]
            _data["contingencies"] = Vector{Contingency}()
            _result = solve_opf_slack(_data, ACRPowerModel, solver_nlp)#, weighted_obj_func_factory([objective_min_fuel_and_flow_cost, objective_min_slack_relative_linear], [data["cost_coeffs"]["redispatch"], data["cost_coeffs"]["overload"]]), slack_types=slack_types)
            reset_ac_feas_params!(_data, limits_generator, tal_timestep, slack_types)
        elseif application == "AC_SCOPF_TIME_SLACK"
            limits_generator, tal_timestep, slack_types = set_ac_feas_params!(_data)
            _data["chosen_contingencies"] = _data["contingencies"]
            _result = solve_secure_opf_slack(_data, ACRPowerModel, solver_nlp)
            reset_ac_feas_params!(_data, limits_generator, tal_timestep, slack_types)
        end
        result = Dict{String, Any}("solution"=>Dict{String, Any}("nw" => Dict{String, Any}("1" => _result["solution"])), "termination_status" => _result["termination_status"], "objective" => _result["objective"], "solve_time" => _result["solve_time"])
        
    elseif application == "DC_OTS"
        result = solve_ots(data, DCMPPowerModel, solver_convex)
    elseif application == "AC_OTS"
        result = solve_ots(data, ACPPowerModel, solver_minlp)
    elseif application == "DC_SCOTS"
        result = solve_secure_ots(data, DCPPowerModel, solver_convex, objective_min_fuel_and_flow_cost)
    elseif application == "DC_SCOTS_LAZY"
        result = solve_secure_ots_lazy(data, DCPPowerModel, solver_convex, objective_min_fuel_and_flow_cost)
    elseif application == "DC_SCOTS_DECOMP"
        result = solve_secure_ots_decomp(data, DCPPowerModel, solver_convex, objective_min_fuel_and_flow_cost)
    elseif application == "AC_SCOTS"
        result = solve_secure_ots(data, ACPPowerModel, solver_minlp, objective_min_fuel_and_flow_cost)
    elseif application == "DC_SCOTS_TIME"
        result = solve_secure_ots_time(data, DCPPowerModel, solver_convex, objective_min_redispatch, slack_types=Vector{Symbol}())
    elseif application == "AC_SCOTS_TIME"
        Memento.error(_LOGGER, "AC_SCOTS_TIME not implemented."); return
    elseif application == "DC_SCOTS_TIME_SLACK"
        result = solve_secure_ots_time(data, DCPPowerModel, solver_convex, weighted_obj_func_factory([objective_min_num_gen_switch, objective_min_num_tr_switch, objective_min_num_branch_switch, objective_min_redispatch_cost, objective_min_num_redispatch, objective_min_slack_relative_linear, objective_min_num_slacks], [data["cost_coeffs"]["gen_sw"], data["cost_coeffs"]["tr_sw"], data["cost_coeffs"]["br_sw"], data["cost_coeffs"]["redispatch"], data["cost_coeffs"]["num_redispatch"], data["cost_coeffs"]["overload"], data["cost_coeffs"]["num_overload"]]), slack_types=slack_types)
    elseif application == "DC_SCOTS_TIME_SLACK_BENDERS"
        result = solve_secure_ots_time_benders(data, DCPPowerModel, solver_convex, weighted_obj_func_factory([objective_min_num_gen_switch, objective_min_num_tr_switch, objective_min_num_branch_switch, objective_min_redispatch_cost, objective_min_num_redispatch, objective_min_slack_relative_linear, objective_min_num_slacks], [data["cost_coeffs"]["gen_sw"], data["cost_coeffs"]["tr_sw"], data["cost_coeffs"]["br_sw"], data["cost_coeffs"]["redispatch"], data["cost_coeffs"]["num_redispatch"], data["cost_coeffs"]["overload"], data["cost_coeffs"]["num_overload"]]), slack_types=slack_types)    elseif application == "AC_SCOTS_TIME_SLACK"
    elseif application == "DC_SCOTS_TIME_SLACK_BENDERS_NESTED" || application == "DC_SCOTS_TIME_SLACK_BENDERS_BENDERS"
        result = solve_secure_ots_time_benders_nested(data, DCPPowerModel, solver_convex, slack_types=slack_types)
    elseif application == "DC_SCOTS_SLACK"
        result = solve_secure_ots_slack(data, DCPPowerModel, solver_convex, objective_min_slack_relative_linear, slack_types=slack_types)
    elseif application == "DC_SCOTS_SLACK_MULTI"
        result = solve_secure_ots_slack(data, DCPPowerModel, solver_convex, ([objective_min_slack_relative_linear, objective_min_slack_relative_linear], [0.5, 1]), slack_types=slack_types)
    elseif application == "AC_SCOTS_SLACK"
        result = solve_secure_ots_slack(data, ACPPowerModel, solver_minlp, objective_min_slack_relative_linear, slack_types=slack_types)
    elseif application == "DATA_ONLY"
        return data
    elseif occursin("AC_SCOTS_", application)
        network_formulation_model = get_formulation_from_app(application)
        result = solve_secure_ots(data, network_formulation_model, solver_convex, objective_min_fuel_and_flow_cost)
    elseif application == "PLOT_ONLY"
        data["nw"]["1"]["per_unit"] = data["per_unit"]
        result = solve_plot_problem(data["nw"]["1"], DCPPowerModel, solver_convex, objective_zero, slack_types=slack_types)
    elseif application ∈ ["DC_SCOTS_TIME_SLACK_BENDERS_FULL_LP", "DC_SCOTS_TIME_SLACK_BENDERS_LP", "DC_SCOTS_TIME_SLACK_BENDERS_DUAL_JL", "DC_SCOTS_TIME_SLACK_BENDERS_PRIMAL_JL"]
        tmp_solver_type = deepcopy(PARAMS[:BENDERS][:SUBPROBLEM_SOLVER_TYPE])
        if application     == "DC_SCOTS_TIME_SLACK_BENDERS_FULL_LP"
            PARAMS[:BENDERS][:SUBPROBLEM_SOLVER_TYPE] = FullLPSolver
        elseif application == "DC_SCOTS_TIME_SLACK_BENDERS_LP"
            PARAMS[:BENDERS][:SUBPROBLEM_SOLVER_TYPE] = LPSolver
        elseif occursin("JL", application)
            PARAMS[:BENDERS][:SUBPROBLEM_SOLVER_TYPE] = LinEqSysSolver
        end

        tmp_dual_solve = deepcopy(PARAMS[:BENDERS][:SUBPROBLEM_SOLVER_DUAL_SOLVE])
        if occursin("PRIMAL", application)
            PARAMS[:BENDERS][:SUBPROBLEM_SOLVER_DUAL_SOLVE] = false
        else
            PARAMS[:BENDERS][:SUBPROBLEM_SOLVER_DUAL_SOLVE] = true
        end

        result = solve_secure_ots_time_benders(data, DCPPowerModel, solver_convex, weighted_obj_func_factory([objective_min_num_gen_switch, objective_min_num_tr_switch, objective_min_num_branch_switch, objective_min_redispatch_cost, objective_min_num_redispatch, objective_min_slack_relative_linear, objective_min_num_slacks], [data["cost_coeffs"]["gen_sw"], data["cost_coeffs"]["tr_sw"], data["cost_coeffs"]["br_sw"], data["cost_coeffs"]["redispatch"], data["cost_coeffs"]["num_redispatch"], data["cost_coeffs"]["overload"], data["cost_coeffs"]["num_overload"]]), slack_types=slack_types)

        PARAMS[:BENDERS][:SUBPROBLEM_SOLVER_TYPE] = tmp_solver_type
        PARAMS[:BENDERS][:SUBPROBLEM_SOLVER_DUAL_SOLVE] = tmp_dual_solve
    else
        Memento.error(_LOGGER, "Functionality $application with $timesteps time step(s) is not implemented.")
    end
    end
    if result != []
        result["time"] = t
    else
        Memento.error(_LOGGER, "Problem could not be solved.")
    end
    
    return result
end 