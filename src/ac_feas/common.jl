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

if Sys.which("matlab") !== nothing && haskey(Pkg.project().dependencies, "MATLAB") # binary and package MATLAB is installed
    include("matpower.jl")
end

"prints information about the ac-feasibility of the given result dict and network id"
function print_info_acfeas(result::Dict{String, Any}, nw_id::String)
    @assert(haskey(result["solution"]["nw"][nw_id], "ac_feasibility"), "No AC feasibility information available for network $nw_id.")
    ac_feasibility = result["solution"]["nw"][nw_id]["ac_feasibility"]
    println("AC-feasibility check for network $nw_id:")
    println("Runtime: ", ac_feasibility["run_time"])
    if ac_feasibility["solved"]
        println("AC feasibility solved")
        println("Objective: ", ac_feasibility["objective"])
        println("Angle Violation sum: ", ac_feasibility["angle_violations_sum"])
        println("Rate Violation sum: ", ac_feasibility["rate_violations_sum"])
        if length(ac_feasibility["angle_violations"]) > 0
            println("Branches with angle violations: ", join(ac_feasibility["angle_violations"], " "))
        end
        if length(ac_feasibility["rate_violations_from"]) > 0
            println("Branches with rate violations from: ", join(ac_feasibility["rate_violations_from"], " "))
        end
        if length(ac_feasibility["rate_violations_to"]) > 0
            println("Branches with rate violations to: ", join(ac_feasibility["rate_violations_to"], " "))
        end
    else
        println("AC feasibility not solved")
    end
    println()
end

"returns the index of the replacement branch (determined recursively) of the given connecting branch, whose breaker branch is switched on and who appears in the contracted network"
function get_mapped_branch(actual_branches::Dict{String, Any}, calculated_indices::Vector{String}, connecting_branches_mapping::Dict{String, Vector{Tuple{String, String}}}, connecting_branch::String)
    for (replacement_branch, breaker_branch) in connecting_branches_mapping[connecting_branch]
        i = findfirst(x -> x == replacement_branch, calculated_indices)
        if isnothing(i)
            return get_mapped_branch(actual_branches, calculated_indices, connecting_branches_mapping, replacement_branch)
        end
        return i
    end
    return nothing
end

"modifies the result dict with the given ac-feasibility information: result and runtime of the ac-feasibility check, ac power flows and violations of the branches and bus voltages"
function modify_result!(
        result::Dict{String, Any}, nw_id::String, objectives::Tuple{Bool, Float64, Float64, Float64}, run_time::Float64,
        rate_violations_from::Vector{Float64}, rate_violations_to::Vector{Float64}, angle_violations::Vector{Float64},
        pfs::Vector{Float64}, qfs::Vector{Float64}, pts::Vector{Float64}, qts::Vector{Float64}, vms::Vector{Float64}, vas::Vector{Float64},
        indices_branches::Vector{String} = map(string, 1:length(pfs)), indices_buses::Vector{String} = map(string, 1:length(vms)),
        connecting_branches_mapping::Dict{String, Vector{Tuple{String, String}}} = Dict{String, Vector{Tuple{String, String}}}()
    ) 

    # get the indices of all branches that have violations above 0
    rate_violations_from_indices::Vector{String}, rate_violations_to_indices::Vector{String}, angle_violations_indices::Vector{String} = [], [], []

    function set_branch_data!(branch::Dict{String, Any}, pf::Float64, qf::Float64, pt::Float64, qt::Float64, rate_violation_from::Float64, rate_violation_to::Float64, angle_violation::Float64)
        branch["ac_p_fr"] = pf
        branch["ac_q_fr"] = qf
        branch["ac_p_to"] = pt
        branch["ac_q_to"] = qt
        branch["ac_rate_violation_from"] = rate_violation_from
        branch["ac_rate_violation_to"] = rate_violation_to
        branch["ac_angle_violation"] = angle_violation
    end

    if objectives[1]

        # write pfs, qfs, pts, qts and violations to the result dictionary
        for (i, index) in enumerate(indices_branches)
            
            if !(haskey(result["solution"]["nw"][nw_id]["branch"], string(index)))
                continue
            end

            set_branch_data!(result["solution"]["nw"][nw_id]["branch"][index], pfs[i], qfs[i], pts[i], qts[i], rate_violations_from[i], rate_violations_to[i], angle_violations[i])   
            
            rate_violations_from[i] > 0 && push!(rate_violations_from_indices, index)
            rate_violations_to[i] > 0 && push!(rate_violations_to_indices, index)
            angle_violations[i] > 0 && push!(angle_violations_indices, index)
        end

        # write ac power flows and violations to the connecting branches (removed in the contracted network)
        for connecting_branch in keys(connecting_branches_mapping)
            
            if !haskey(result["solution"]["nw"][nw_id]["branch"], connecting_branch)
                continue
            end

            i = get_mapped_branch(result["solution"]["nw"][nw_id]["branch"], indices_branches, connecting_branches_mapping, connecting_branch)
            set_branch_data!(result["solution"]["nw"][nw_id]["branch"][index], pfs[i], qfs[i], pts[i], qts[i], rate_violations_from[i], rate_violations_to[i], angle_violations[i])
            
            rate_violations_from[i] > 0 && push!(rate_violations_from_indices, connecting_branch)
            rate_violations_to[i] > 0 && push!(rate_violations_to_indices, connecting_branch)
            angle_violations[i] > 0 && push!(angle_violations_indices, connecting_branch)
        end

        # write vms and vas to the result dictionary
        for (i, index) in enumerate(indices_buses)
            result["solution"]["nw"][nw_id]["bus"][index]["vm_ac"] = vms[i]
            result["solution"]["nw"][nw_id]["bus"][index]["va_ac"] = vas[i]
        end
    end
    
    # provide general ac-feasibility information and give violation indices
    result["solution"]["nw"][nw_id]["ac_feasibility"] = Dict(
        "solved" => objectives[1], 
        "objective" => objectives[2],
        "angle_violations_sum" => objectives[3],
        "rate_violations_sum" => objectives[4],
        "run_time" => run_time,
        "angle_violations" => angle_violations_indices,
        "rate_violations_from" => rate_violations_from_indices,
        "rate_violations_to" => rate_violations_to_indices
    )

end

"given power flows, rate limits and active status, calculates the violations of the branches"
function calculate_rate_violations(pfs::Vector{Float64}, qfs::Vector{Float64}, pts::Vector{Float64}, qts::Vector{Float64}, rates::Vector{Float64}, active::Vector{Bool})

    # calculate the apparent power flows
    Sfs = sqrt.(pfs.^2 + qfs.^2)
    Sts = sqrt.(pts.^2 + qts.^2)

    # obtain violations by calculating the diff of the apparent power flows and the rate limits
    rate_violations_from, rate_violations_to = (Sfs .* active .- rates), (Sts .* active .- rates)

    # set the diff to 0 if it is smaller than the precision
    rate_violations_from .= @. ifelse(rate_violations_from < SOLVER_PRECISION, 0, rate_violations_from) 
    rate_violations_to .= @. ifelse(rate_violations_to < SOLVER_PRECISION, 0, rate_violations_to) 
    
    return rate_violations_from, rate_violations_to
end

"Given the voltage angles and the min/max voltage angle differences, calculates the angle difference violations of the branches"
function calculate_angle_diff_violations(va_frs::Vector{Float64}, va_tos::Vector{Float64}, angmins::Vector{Float64}, angmaxs::Vector{Float64})

    # calculate the angle difference
    angle_diffs = va_tos .- va_frs

    # obtain violations by calculating the diff of the apparent power flows and the rate limits
    angle_violations = @. ifelse(angle_diffs < angmins, 
                                 angmins - angle_diffs, 
                                 ifelse(angle_diffs > angmaxs, angle_diffs - angmaxs, 0.0))
                                 
    return angle_violations
end

"calculates ac power flows and voltages for the given solution and network data"
function calculate_branch_data(
        solution_pf::Dict{String, Any},
        data::Dict{String, Any},
        power_model::Type,
        vms::Vector{Float64}, 
        vas::Vector{Float64},
        indices_buses::Vector{String} = map(string, 1:length(vms))
    )

    indices_branches, pfs, qfs, pts, qts, rates, va_frs, va_tos, angmins, angmaxs, active = String[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Bool[]
        
    if power_model == ACRPowerModel
        for (branch_id, branch) in data["branch"]

            push!(indices_branches, branch_id)
            push!(rates, get_adjusted_rate_a(data, branch))
            push!(active, get_status(branch) > 0.5)

            if get_status(branch) < 0.5
                push!(pfs, 0.0)
                push!(qfs, 0.0)
                push!(pts, 0.0)
                push!(qts, 0.0)
                push!(va_frs, 0.0)
                push!(va_tos, 0.0)
                push!(angmins, branch["angmin"])
                push!(angmaxs, branch["angmax"])
                continue
            end

            f_bus, t_bus = string(branch["f_bus"]), string(branch["t_bus"])
            vm_fr, vm_to = vms[findfirst(x -> x == f_bus, indices_buses)], vms[findfirst(x -> x == t_bus, indices_buses)]
            va_fr, va_to = vas[findfirst(x -> x == f_bus, indices_buses)], vas[findfirst(x -> x == t_bus, indices_buses)]
            push!(va_frs, va_fr)
            push!(va_tos, va_to)
            push!(angmins, branch["angmin"])
            push!(angmaxs, branch["angmax"])

            push!(pfs, solution_pf["branch"][branch_id]["pf"])
            push!(qfs, solution_pf["branch"][branch_id]["qf"])

            push!(pts, solution_pf["branch"][branch_id]["pt"])
            push!(qts, solution_pf["branch"][branch_id]["qt"])
        end
    else
        for (branch_id, branch) in data["branch"]

            push!(indices_branches, branch_id)
            push!(rates, get_adjusted_rate_a(data, branch))
            push!(active, branch["br_status"] > 0.5)

            f_bus, t_bus = string(branch["f_bus"]), string(branch["t_bus"])
            vm_fr, vm_to = solution_pf["bus"][f_bus]["vm"], solution_pf["bus"][t_bus]["vm"]
            va_fr, va_to = solution_pf["bus"][f_bus]["va"], solution_pf["bus"][t_bus]["va"]
            push!(va_frs, va_fr)
            push!(va_tos, va_to)
            push!(angmins, branch["angmin"])
            push!(angmaxs, branch["angmax"])

            g, b = calc_branch_y(branch)
            tr, ti = calc_branch_t(branch)
            g_fr, b_fr = branch["g_fr"], branch["b_fr"]
            tm = branch["tap"]
            g_to, b_to = branch["g_to"], branch["b_to"]

            p_fr = (g+g_fr)/tm^2*vm_fr^2 + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to))
            q_fr = -(b+b_fr)/tm^2*vm_fr^2 - (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to))
            push!(pfs, p_fr)
            push!(qfs, q_fr)

            p_to = (g+g_to)*vm_to^2 + (-g*tr-b*ti)/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/tm^2*(vm_to*vm_fr*sin(va_to-va_fr))
            q_to = -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-g*tr-b*ti)/tm^2*(vm_to*vm_fr*sin(va_to-va_fr))
            push!(pts, p_to)
            push!(qts, q_to)
        end
    end

    return pfs, qfs, pts, qts, rates, va_frs, va_tos, angmins, angmaxs, active, indices_branches
end

function calculate_bus_data(solution_pf::Dict{String, Any}, data::Dict{String, Any}, power_model::Type)
    indices_buses, vms, vas = String[], Float64[], Float64[]
    if power_model == ACRPowerModel 
        for bus_id in keys(data["bus"])
            push!(indices_buses, bus_id)
            vr = solution_pf["bus"][bus_id]["vr"]
            vi = solution_pf["bus"][bus_id]["vi"]
            push!(vms, sqrt(vr^2 + vi^2))
            push!(vas, atan(vi / vr))
        end
        return vms, vas, indices_buses
    else
        indices_buses, vms, vas = String[], Float64[], Float64[]
        for bus_id in keys(data["bus"])
            push!(indices_buses, bus_id)
            push!(vms, solution_pf["bus"][bus_id]["vm"])
            push!(vas, solution_pf["bus"][bus_id]["va"])
        end
    end
    return vms, vas, indices_buses
end

"sets params neccessary for the ac_feas check"
function set_ac_feas_params!(data::Dict{String, Any})::Tuple{Dict{Symbol, Any}, Int64, Vector{Symbol}}
    
    limit_generator, tal_timestep = override_ac_feas_params()

    slack_types = deepcopy(data["slack_types"])
    !(:pd in data["slack_types"]) && push!(data["slack_types"], :pd)
    !(:qd in data["slack_types"]) && push!(data["slack_types"], :qd)

    return limit_generator, tal_timestep, slack_types
end

function override_ac_feas_params()::Tuple{Dict{Symbol, Any}, Int64}
    limit_generator = deepcopy(PARAMS[:LIMITS][:GENERATOR])
    tal_timestep = PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP]
    tal_timestep = PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP]
    extend_dict!(
        PARAMS, 
        Dict{Symbol, Any}(
            :LIMITS => Dict(
                :BRANCH => Dict(
                    :TAL_TIMESTEP => 0
                ),
                :GENERATOR => Dict(
                    :GEN_SWITCHING => false,
                    :GEN_REDISPATCH_OVERRIDE => true,
                    :GEN_REDISPATCH_OVERRIDE_DATA => true
                )
            )
        )
    )
    @info(_LOGGER, "Overwrote LIMITS.BRANCH.TAL_TIMESTEP, LIMITS.GENERATOR.GEN_SWITCHING and LIMITS.GENERATOR.GEN_REDISPATCH_OVERRIDE.")
    return limit_generator, tal_timestep
    @info(_LOGGER, "Overwrote LIMITS.BRANCH.TAL_TIMESTEP, LIMITS.GENERATOR.GEN_SWITCHING and LIMITS.GENERATOR.GEN_REDISPATCH_OVERRIDE.")
    return limit_generator, tal_timestep
end

"resets params neccessary for the ac_feas check"
function reset_ac_feas_params!(data::Dict{String, Any}, limit_generator::Dict{Symbol, Any}, tal_timestep::Int64, slack_types::Vector{Symbol})
    PARAMS[:LIMITS][:GENERATOR] = limit_generator
    PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP] = tal_timestep
    data["slack_types"] = slack_types
end

"setting all values in the data dict neccessary for the ac power flow calculation. contingency selection can be :all for all contingencies, :relevant for only contingencies relevant to the dc solution or :none for no contingencies"
function set_ac_data_values!(data::Dict{String, Any}, result::Dict{String, Any}, contingency_selection::Symbol = :all)

    convert_to_generic_cost_data!(data)
    # set start values for AC power flow if known
    set_dc_pf_start_values!(data, result["solution"], fix_uc=false, fix_ts=false, fix_benders=false, fix_pf=false, fix_slack=false, fix_pg=false)
    
    if contingency_selection == :all
        @warn("All contingencies are included in the AC feasibility check. This might lead to long computation times.")
    elseif contingency_selection == :relevant
        @warn("Only contingencies relevant to the DC solution are included in the AC feasibility check. This might lead to discrepencies between true AC feasibility and reported AC feasibility.")
        relevant_contingencies = get_relevant_contingencies(data, result=result)
    else
        @warn("No contingencies are included in the AC feasibility check. Quality of AC feasibility check will be reduced significantly.")
    end

    connecting_branches_mappings = Dict{String, Any}()
    for nw_id in keys(result["solution"]["nw"])

        # contract network graph to remove breakers not handled by solve_opf_slack
        connecting_branches_mappings[nw_id] = contract_breakers!(data["nw"][nw_id])

        # select contingencies to be included
        if contingency_selection == :all
            data["nw"][nw_id]["chosen_contingencies"] = data["nw"][nw_id]["contingencies"]
        elseif contingency_selection == :relevant
            data["nw"][nw_id]["chosen_contingencies"] = relevant_contingencies
        else
            data["nw"][nw_id]["chosen_contingencies"] = [] 
        end
        
        # solve a slack optimal power flow model
        data["nw"][nw_id]["cost_coeffs"] = data["cost_coeffs"]
    end
    return connecting_branches_mappings
end

"""Get result vector of ac scopf slack"""
function set_result_vector!(result::Dict{String, Any}, result_opf::Dict{String, Any}, data::Dict{String, Any}, nw_id::String, start_time::Float64, connecting_branches_mappings::Dict{String, Any}, set_result::Bool=true)
    result_vector::Tuple{Bool, Float64, Float64, Float64} = (false, 0.0, 0.0, 0.0)
    if result_opf["termination_status"] == JuMP.MOI.OPTIMAL || result_opf["termination_status"] == JuMP.MOI.LOCALLY_SOLVED
            
        # calculate the branch data, the violations and obtain the objective value
            vms, vas, indices_buses = calculate_bus_data(result_opf["solution"], data["nw"][nw_id], PARAMS[:POSTPROCESSING][:AC_FEAS][:POWER_MODEL])
            pfs, qfs, pts, qts, rates, va_frs, va_tos, angmins, angmaxs, active, indices_branches = calculate_branch_data(result_opf["solution"], data["nw"][nw_id], PARAMS[:POSTPROCESSING][:AC_FEAS][:POWER_MODEL], vms, vas, indices_buses)
            rate_violations_from, rate_violations_to = calculate_rate_violations(pfs, qfs, pts, qts, rates, active)
            angle_violations::Vector{Float64} = calculate_angle_diff_violations(va_frs, va_tos, angmins, angmaxs)
            
            branch_changes = sum([branch["br_change"] for (branch_id, branch) in result["solution"]["nw"][nw_id]["branch"] if is_branch_switchable(data["nw"][nw_id]["branch"][branch_id])])
            tr_changes = PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_SHIFTING] ? sum([branch["tr_shift_change"] for (branch_id, branch) in result["solution"]["nw"][nw_id]["branch"] if is_branch_shiftable(data["nw"][nw_id]["branch"][branch_id])]) : 0

            result_opf["objective"] += branch_changes * data["cost_coeffs"]["br_sw"] + tr_changes * data["cost_coeffs"]["tr_sw"]
        
            run_time = time() - start_time
            result_vector = (true, result_opf["objective"], sum(angle_violations), sum(rate_violations_from) + sum(rate_violations_to))
            set_result && modify_result!(result, nw_id, result_vector, run_time, rate_violations_from, rate_violations_to, angle_violations, pfs, qfs, pts, qts, vms, vas, indices_branches, indices_buses, connecting_branches_mappings[nw_id])
        else
            run_time = time() - start_time
            set_result && modify_result!(result, nw_id, result_vector, run_time, Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{String}(), Vector{String}(), connecting_branches_mappings[nw_id])
        end      
    return result_vector
end

"checks the ac feasibility of a given network and result dict using a linear programming solver"
function check_ac_feasibility_lp!(result::Dict{String, Any}, data::Dict{String, Any}, solver_nlp::JuMP.MOI.OptimizerWithAttributes)

    # check if a solver is provided
    @assert(solver_nlp !== nothing, "No solver provided for AC feasibility check.")

    result_opfs = Dict{String, Any}()
    connecting_branches_mappings = set_ac_data_values!(data, result, :all) 
    limit_generator, tal_timestep, slack_types = set_ac_feas_params!(data)
    
    for nw_id in keys(result["solution"]["nw"])

        start_time = time()
        
        result_opf = solve_secure_opf_slack(data["nw"][nw_id], PARAMS[:POSTPROCESSING][:AC_FEAS][:POWER_MODEL], solver_nlp)
        set_result_vector!(result, result_opf, data, nw_id, start_time, connecting_branches_mappings)
        result_opfs[nw_id] = result_opf

        # print the according information
        print_info_acfeas(result, nw_id)

    end
    reset_ac_feas_params!(data, limit_generator, tal_timestep, slack_types)
    return result_opfs
end

"updates the network data with the given result solution"
function update_network!(data::Dict{String, Any}, result_solution::Dict{String, Any}, nw_id::String)

    # delete branches which are switched off from data_new
    for (branch_id, branch) in data["nw"][nw_id]["branch"]
        if is_branch_switchable(branch) && get_status(result_solution["branch"][branch_id]) <= 0.5
            data["nw"][nw_id]["branch"][branch_id]["br_status"] = 0.0
            data["nw"][nw_id]["branch"][branch_id]["status"] = 0.0
        end
    end

    # set generator power in data_new to values of the DC solution
    for gen_id in keys(result_solution["gen"])
        if haskey(result_solution["gen"][gen_id], "gen_status") && result_solution["gen"][gen_id]["gen_status"] <= 0.5
            data["nw"][nw_id]["gen"][gen_id]["gen_status"] = 0.0
        end
        if haskey(result_solution["gen"][gen_id], "status") && result_solution["gen"][gen_id]["status"] <= 0.5
            data["nw"][nw_id]["gen"][gen_id]["status"] = 0.0
        end
        if (haskey(result_solution["gen"][gen_id], "status") && result_solution["gen"][gen_id]["status"] >= 0.5) || (haskey(result_solution["gen"][gen_id], "gen_status") && result_solution["gen"][gen_id]["gen_status"] >= 0.5)
            data["nw"][nw_id]["gen"][gen_id]["pg"] = result_solution["gen"][gen_id]["pg"]
        end
    end

    # set initial voltage angle
    for bus_id in keys(result_solution["bus"])
        if haskey(result_solution["bus"][bus_id], "va")
            data["nw"][nw_id]["bus"][bus_id]["va"] = result_solution["bus"][bus_id]["va"]
        end
    end

    # remove the load of the loss estimation
    if haskey(data["nw"][nw_id], "load")
        for load_id in keys(data["nw"][nw_id]["load"])
            if get(data["nw"][nw_id]["load"][load_id], "loss_load", false)
                delete!(data["nw"][nw_id]["load"], load_id)
            end
        end
    end
end

"checks the ac feasibility of a given network with respect to the given result dict and writes the results to the result dict"
function check_ac_feasibility!(result::Dict{String, Any}, data::Dict{String, Any}, solver_nlp::JuMP.MOI.OptimizerWithAttributes = nothing)
    
    if result["termination_status"] != JuMP.MOI.OPTIMAL
        println("No optimal DC solution found, skipping AC feasibility check.")
        return
    end

    # update the network according to the result
    data_copy = deepcopy(data)
    for nw_id in keys(result["solution"]["nw"])
        update_network!(data_copy, result["solution"]["nw"][nw_id], nw_id)
    end

    # check for AC feasibility (method depending on PARAMS[:POSTPROCESSING][:AC_FEAS][:METHOD])
    if PARAMS[:POSTPROCESSING][:AC_FEAS][:METHOD] == :MATPOWER
        check_ac_feasibility_matpower!(result, data_copy)
    elseif PARAMS[:POSTPROCESSING][:AC_FEAS][:METHOD] == :LP
        return check_ac_feasibility_lp!(result, data_copy, solver_nlp)
    else
        error("AC feasibility method not supported.")
    end

end