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

export estimate_losses!

import JuMP

if Sys.which("matlab") !== nothing && haskey(Pkg.project().dependencies, "MATLAB") # binary and package MATLAB is installed
    include("matpower.jl")
    include("../ac_feas/matpower.jl")
end

"add the given loss to the given network and add a bus load accordingly"
function add_loss!(network::Dict{String, Any}, loss::Float64, branch_index::String)
    # add the losses to the branch
    network["branch"][branch_index]["loss"] = loss

    # add a new load to the bus of the branch
    loss == 0 && return
    bus_index = loss >= 0 ? network["branch"][branch_index]["t_bus"] : network["branch"][branch_index]["f_bus"]
    load_index = length(network["load"]) + 1
    @assert(!haskey(network["load"], load_index), "Load index already exists.")
    load = create_bus_load(abs(loss), 0.0, bus_index, 1, load_index)
    load["loss_load"] = true
    network["load"][string(load_index)] = load
end

"add the given losses to the given network and add bus loads accordingly"
function add_losses!(
    network::Dict{String, Any}, 
    losses::Vector{Float64}, 
    indices_branches::Vector{String} = map(string, 1:length(losses)),
    connecting_branches_mapping::Dict{String, Vector{Tuple{String, String}}} = Dict{String, Vector{Tuple{String, String}}}()
)
    
    # add load entry to network if not present
    if !haskey(network, "load")
        network["load"] = Dict{String, Dict{String, Any}}()
    end

    # add loads to existing branches
    for (i, index) in enumerate(indices_branches)
        if haskey(network["branch"], index)
            add_loss!(network, losses[i], index)
        end        
    end

    # if the network was contraced, maps loads to the original branches
    for connecting_branch in keys(connecting_branches_mapping)
        if !haskey(network["branch"], connecting_branch)
            continue
        end
        i = get_mapped_branch(network["branch"], indices_branches, connecting_branches_mapping, connecting_branch)
        add_loss!(network, losses[i], connecting_branch)
    end

end

"estimates the losses of the branches of a given network using a linear programming solver"
function estimate_losses_lp!(data::Dict{String, Any}, solver_nlp::JuMP.MOI.OptimizerWithAttributes)
    
    # check if a solver is provided
    @assert(solver_nlp !== nothing, "No solver provided for Constant Loss Estimation.")

    overall_run_time = 0.0

    function estimate_losses!(network, cost_coeffs = nothing)
        start_time = time()

        # contract network graph to remove breakers not handled by solve_pf
        contracted_network = deepcopy(network)
        connecting_branches_mapping = contract_breakers!(contracted_network)

        # solve a optimal power flow model
        contracted_network["cost_coeffs"] = cost_coeffs
        if PARAMS[:LOSS_ESTIMATION][:MODEL] == :OPF
            result_opf = solve_opf_slack(contracted_network, PARAMS[:LOSS_ESTIMATION][:POWER_MODEL], solver_nlp)
        elseif PARAMS[:LOSS_ESTIMATION][:MODEL] == :SCOPF
            contracted_network["chosen_contingencies"] = [] # TODO replace by actual chosen /relevant contingencies
            @warn("Contingency selection for loss estimation not implemented yet. Running SCOPF without contingencies.")
            result_opf = solve_scopf_slack(contracted_network, PARAMS[:LOSS_ESTIMATION][:POWER_MODEL], solver_nlp)
        else
            error("Model not supported for loss estimation.")
        end

        if result_opf["termination_status"] == JuMP.MOI.OPTIMAL || result_opf["termination_status"] == JuMP.MOI.LOCALLY_SOLVED
            # calculate and add the losses
            vms, vas, indices_buses = calculate_bus_data(result_opf["solution"], contracted_network, PARAMS[:LOSS_ESTIMATION][:POWER_MODEL])
            pfs, _, pts, _, _, _, _, _, _, _, indices_branches = calculate_branch_data(result_opf["solution"], contracted_network, PARAMS[:LOSS_ESTIMATION][:POWER_MODEL], vms, vas, indices_buses)      
            losses = pfs + pts
            add_losses!(network, losses, indices_branches, connecting_branches_mapping)
        end

        return time() - start_time
    end

    if haskey(data, "nw")
        # iterate over all networks
        for nw_id in keys(data["nw"])
            convert_to_generic_cost_data!(data)
            run_time = estimate_losses!(data["nw"][nw_id], data["cost_coeffs"])
            overall_run_time += run_time
        end
    else
        add_missing_attributes(data, PARAMS[:TECHNICAL][:MISSING_LIST])
        overall_run_time = estimate_losses!(data, data["cost_coeffs"])
    end

    data["loss_estimation_time"] = overall_run_time
end

"estimates the losses of a given network and adds the losses to the branches and adds a bus load accordingly"
function estimate_losses!(data::Dict{String, Any}, solver_nlp::JuMP.MOI.OptimizerWithAttributes = nothing)
    # estimate losses (method depending on LOSS_ESTIMATION_METHOD)
    if PARAMS[:LOSS_ESTIMATION][:METHOD] == :MATPOWER
        return estimate_losses_matpower!(data)
    elseif PARAMS[:LOSS_ESTIMATION][:METHOD] == :LP
        return estimate_losses_lp!(data, solver_nlp)
    else
        error("Loss estimation method not supported.")
    end
end