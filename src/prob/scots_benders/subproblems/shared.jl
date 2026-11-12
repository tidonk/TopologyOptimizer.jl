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

# functions used by different subproblem solvers

"extracts and preprocesses data needed in the solving and evaluation process to reduce look-up effort & time spent on repeating computations"
function collect_branch_data(pm_sub::AbstractPowerModel; nw::Int64=nw_id_default)::Dict{Int64, Dict{Symbol, Any}}
    branch_data = Dict{Int64, Dict{Symbol, Any}}() # [branch_id] => ...
    
    obj_func_terms = JuMP.objective_function(pm_sub.model).terms

    # collect data for each branch, mostly we just copy the values from pm_sub.data. This increases storage needs, but
    # we later save the look-ups to get `ref(pm_sub, nw, :branch)`, Symbol look-ups are faster, and each `branch` Dict
    # has fewer entries
    for (branch_id, branch) in ref(pm_sub, nw, :branch)

        branch_data[branch_id] = Dict(:f_bus => branch["f_bus"],
                                      :t_bus => branch["t_bus"],
                                      :switchable => branch["switchable"],
                                      :shiftable => branch["shiftable"],
                                      :status => branch["status"],
                                      :type => branch["type"],
                                      :monitored => branch["monitored"],
                                      :angmin_close => get(branch, "angmin_close", ref(pm_sub, nw, :off_angmin)),
                                      :angmax_close => get(branch, "angmax_close", ref(pm_sub, nw, :off_angmax)),
                                      :electrical => branch["electrical"])

        # some data is only relevant if the branch is electrical -> continue otherwise
        if pm_sub.data["BUSBARS"] && ! branch["electrical"]
            continue
        end
        
        for pair in [:angmin => branch["angmin"],
                     :angmax => branch["angmax"],
                     :b => branch["b"],
                     :index => branch["index"],
                     :shift => branch["shift"],
                     :rate_a => branch["rate_a"],
                     :adjusted_rate_a => get_adjusted_rate_a(pm_sub, branch),
                     :max_loading => get(branch, "max_loading", 100),       # TODO: remove and only use :adjusted_rate_a
                     :max_loading_n1 => get(branch, "max_loading_n1", 100)  # TODO: remove and only use :adjusted_rate_a
                     ]
            push!(branch_data[branch_id], pair)
        end

        # shift_per_step is only needed if the branch is actually shiftable
        if get(branch, "shiftable", false)
            branch_data[branch_id][:shift_per_step] = branch["shift_per_step"]
        end

        branch_data[branch_id][:obj_coeff] = obj_func_terms[var(pm_sub, nw, :rate_a_slack, (branch_id, branch["f_bus"], branch["t_bus"]))]
    end

    return branch_data
end

"extends power-balance duals to detached breaker buses on active electrical lines"
function _update_detached_bus_duals!(solver::AbstractMatrixSolver, dualss, master_var_vals, outage_branches)
    for (branch_id, branch) in solver.branch_data
        is_branch_electrical(branch) || continue
        branch_id in outage_branches && continue
        (is_branch_switchable(branch) ? master_var_vals[:z_branch][branch_id] > 0.5 : branch[:status] == 1) || continue

        f_bus, t_bus = branch[:f_bus], branch[:t_bus]
        f_col, t_col = solver.bus_to_col[f_bus], solver.bus_to_col[t_bus]
        if f_col <= 0 < t_col
            bus_id, value = f_bus, solver.duals[t_col]
        elseif t_col <= 0 < f_col
            bus_id, value = t_bus, solver.duals[f_col]
        else
            continue
        end

        # This zero-injection bus was eliminated from A. Matching its neighbour's
        # dual keeps the dangling line's flow stationarity equation satisfied.
        dualss[value >= 0 ? :power_balance_ub : :power_balance_lb][bus_id] = value
    end
end

function _update_bb_component_duals!(solver::AbstractMatrixSolver, dualss::Dict{Symbol, Dict{Int, Float64}})
    "sets dualss `value` for the branch, and then applies the `value` to the other buses' busbar_sum"
    function _update_vad_single!(bus_id, branch_id, value, busbar_sums)
        branch = solver.branch_data[branch_id]
        direction = branch[:f_bus] == bus_id ? -1 : 1
        
        dual_val = direction * value
        
        # switchable -> on_off
        if is_branch_switchable(branch)
            
            # positive dual value -> ub, else lb
            if dual_val < 0
                dualss[:vad_on_off_ub][branch_id] = dual_val
            else
                dualss[:vad_on_off_lb][branch_id] = dual_val
            end
    
        else
            # positive dual value -> ub, else lb
            if dual_val < 0
                dualss[:vad_ub][branch_id] = dual_val
            else
                dualss[:vad_lb][branch_id] = dual_val
            end

        end
        
        other_bus_id = direction == -1 ? branch[:t_bus] : branch[:f_bus]

        # update the surplus of the corresponding busbar
        busbar_sums[other_bus_id] += value
    end

    # pre-allocate
    bus_to_connected_branches =     [Vector{Tuple{Int64, Float64}}() for _ in 1:solver.n_buses] # store branch_id, dual value
    bus_to_connected_aux_branches = [Vector{Int64}() for _ in 1:solver.n_buses]

    dualss_ohm_fr_on_off_lb = dualss[:ohm_fr_on_off_lb]
    dualss_ohm_fr_on_off_ub = dualss[:ohm_fr_on_off_ub]
    dualss_ohm_fr_lb = dualss[:ohm_fr_lb]
    dualss_ohm_fr_ub = dualss[:ohm_fr_ub]

    for (branch_id, branch) in solver.branch_data

        f_bus_id = branch[:f_bus]::Int64
        t_bus_id = branch[:t_bus]::Int64

        type = branch[:type]::String
        if type == "L"
            b = branch[:b]::Float64

            if is_branch_switchable(branch)
                dual_value = b * (dualss_ohm_fr_on_off_lb[branch_id] + dualss_ohm_fr_on_off_ub[branch_id])
                dual_value += dualss[:vad_on_off_lb][branch_id] + dualss[:vad_on_off_ub][branch_id]
            else
                dual_value = b * (dualss_ohm_fr_lb[branch_id] + dualss_ohm_fr_ub[branch_id])
                dual_value += dualss[:vad_lb][branch_id] + dualss[:vad_ub][branch_id]
            end

            push!(bus_to_connected_branches[f_bus_id], (branch_id,  dual_value))
            push!(bus_to_connected_branches[t_bus_id], (branch_id, -dual_value))
        elseif type == "breaker" # the breaker bus is always the t bus
            push!(bus_to_connected_aux_branches[t_bus_id], branch_id)
        else
            push!(bus_to_connected_aux_branches[f_bus_id], branch_id)
            push!(bus_to_connected_aux_branches[t_bus_id], branch_id)
        end
    end

    # we need to solve an uncapacited flow problem in each substation component. Smaller overall values are better (maybe), but any solution will do
    for station_component ∈ solver.components

        busbar_sums = Dict(bus_id => 0.0 for bus_id in station_component.busbar_buses)

        # breakers first, since they are always leaves
        for bus_id ∈ station_component.breaker_buses

            # calculate bus surplus
            busbar_sum = 0.0
            for (branch_id, dual_value) in bus_to_connected_branches[bus_id]
                busbar_sum += dual_value
            end
            
            # there should only be one active breaker branch, send all of the surplus over it
            for branch_id in bus_to_connected_aux_branches[bus_id]
                if branch_id in station_component.breaker_branches || branch_id in station_component.coupler_branches
                    _update_vad_single!(bus_id, branch_id, busbar_sum, busbar_sums)
                end
            end

        end

        # calculate surplus at each busbar
        for bus_id ∈ station_component.busbar_buses

            busbar_sum = busbar_sums[bus_id]
            for (branch_id, dual_value) in bus_to_connected_branches[bus_id]
                busbar_sum += dual_value
            end

            busbar_sums[bus_id] =  busbar_sum
        end

        # handle all busbars - we construct a spanning tree and then iteratively eliminate its leaves
        # start with the element with maximum abs sum, this ensures a min flow for ≤ 3 busbars
        visited_busbars = [argmax((bus_id->abs(busbar_sums[bus_id])), station_component.busbar_buses)]
        predecessor_branches = Vector{Int}()
        
        # we do a breadth-first search
        pos=1
        while length(visited_busbars) < length(station_component.busbar_buses)
            @assert(pos <= length(visited_busbars), "Some busbars are not reachable from the starting busbar. Substation: $station_component")
            # get the next busbar
            bus_id = visited_busbars[pos]
            
            # check for each incident branch whether the other side is unvisited, and add it if so
            for branch_id in bus_to_connected_aux_branches[bus_id]
                branch = solver.branch_data[branch_id]
                other_bus_id = branch[:f_bus] == bus_id ? branch[:t_bus] : branch[:f_bus]
                
                if other_bus_id ∉ visited_busbars && other_bus_id ∈ station_component.busbar_buses
                    push!(visited_busbars, other_bus_id)
                    push!(predecessor_branches, branch_id)
                end
            end
            
            pos += 1
        end

        # assign values - we iterate backwards and thus always have a leaf (the last (pos=1) busbar will have surplus 0)
        for pos = length(visited_busbars):-1:2
            bus_id = visited_busbars[pos]
            busbar_sum = busbar_sums[bus_id]

            # send all of the surplus to the "predecessor"
            _update_vad_single!(bus_id, predecessor_branches[pos-1], busbar_sum, busbar_sums)
        end

    end
end

"builds a column used to extend the solvers matrix A into its feasibility cut equivalent from branch"
function build_feas_col(solver::AbstractMatrixSolver, branch::Dict{Symbol, Any}, m::Int64)
    # 1.0 at f_bus, -1.0 at t_bus, else 0.0
    return sparse([solver.bus_to_col[branch[:f_bus]], solver.bus_to_col[branch[:t_bus]]], [1, 1], [1.0, -1.0], m, 1)
end

"gets a diagonal value used to extend the solvers matrix A into its feasibility cut equivalent from branch"
function get_feas_diag(master_var_vals::Dict{Symbol, AbstractVector}, outage_branches::Vector{Int64}, branch_id::Int64, branch::Dict{Symbol, Any}, vad::Float64)::Float64
    # we need the the violated angle limit
    is_off = is_branch_switchable(branch) ? (master_var_vals[:z_branch][branch_id]::Float64 < 0.5) : (branch[:status]::Int == 0)
        
    if branch_id ∈ outage_branches || is_off
        angmin = branch[:angmin_close]::Float64
        angmax = branch[:angmax_close]::Float64
    else
        angmin = branch[:angmin]::Float64
        angmax = branch[:angmax]::Float64
    end

    if vad < angmin - SOLVER_PRECISION::Float64
        diag = angmin
    else
        diag = angmax
    end

    return diag
end