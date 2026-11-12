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

# busbar version of the linear equation system solver
"intializes the solver"
function init_lineqsysbb_solvers(pm, pm_subs, n_threads::Int64, sub_cons_for_dualss; nw::Int64=nw_id_default)::Vector{LinEqSysBBSolver}
    components, detached_components, branch_to_comp = calc_full_busbar_components(pm_subs[1], [branch_id for (branch_id, branch) in ref(pm, nw, :branch) if get(branch, "status", branch["br_status"]) == 0], nw=nw)
    #@assert(components.busbar_buses, "Substation $component has no busbars")
    A, power_surplus, bus_to_col = build_base_pf_matrix(pm_subs[1], components, detached_components, nw=nw)

    # add FullLP solver as fallback solver (one per thread)
    fallback_solver_solvers = init_fulllp_solvers(pm_subs, n_threads, sub_cons_for_dualss)
    fallback_solvers = [FallbackSolver(solver=fallback_solver_solvers[t], pm_sub=pm_subs[t], nw=nw) for t in 1:n_threads]

    branch_data = collect_branch_data(pm_subs[1], nw=nw)

    return [[LinEqSysBBSolver(A=A, power_surplus=power_surplus, bus_to_col=bus_to_col, components=components, detached_components=detached_components, branch_to_comp=branch_to_comp, fallback_solver=fallback_solvers[1], branch_data=branch_data)]; [LinEqSysBBSolver(A=deepcopy(A), power_surplus=deepcopy(power_surplus), bus_to_col=deepcopy(bus_to_col), components=deepcopy(components), detached_components=deepcopy(detached_components), branch_to_comp=deepcopy(branch_to_comp), fallback_solver=fallback_solvers[t], branch_data=branch_data) for t in 2:n_threads]]
end

"applies the current master solution to solve"
function apply_master!(solver::LinEqSysBBSolver, pm::AbstractPowerModel, master_var_vals::Dict{Symbol, AbstractVector}; nw::Int64=nw_id_default)
    add_master_trafo_vals!(solver, master_var_vals)

    solver.A, solver.power_surplus = apply_master_vals!(solver.A, pm, solver.power_surplus, solver.bus_to_col, solver.components, solver.detached_components, solver.branch_to_comp, master_var_vals, nw=nw)
end

"resets the solver back to its base state"
function reset_master!(solver::LinEqSysBBSolver, pm::AbstractPowerModel, master_var_vals::Dict{Symbol, AbstractVector}; nw::Int64=nw_id_default)
    solver.A, solver.power_surplus = reset_master_vals!(solver.A, pm, solver.power_surplus, solver.bus_to_col, solver.components, solver.detached_components, solver.branch_to_comp, master_var_vals, nw=nw)
end

"applies one contingency to the solver"
function apply_contingency!(solver::LinEqSysBBSolver, pm::AbstractPowerModel, branches::Vector{Int64}, master_var_vals::Dict{Symbol, AbstractVector}, int_sol::Bool; nw::Int64=nw_id_default)::Dict{Int64, Dict{Symbol, Float64}}
    @assert(int_sol, "Only integer solutions are supported for the matrix solver")
    for branch_id in branches
        branch = ref(pm, nw, :branch, branch_id)

        if solver.bus_to_col[branch["f_bus"]::Int64] ≤ 0 || solver.bus_to_col[branch["t_bus"]::Int64] ≤ 0
            continue
        end
        
        # remove outaged transformer from solver
        if branch["transformer"]
            reset_trafo_val!(solver.power_surplus, branch, solver.bus_to_col, get_transformer_shift(solver, branch))   
        end

        remove_branch!(solver.A, pm, branch_id, solver.bus_to_col, nw=nw)
        solver.A, solver.power_surplus, _ = update_outage_line!(solver.A, pm, solver.power_surplus, branch_id, solver.bus_to_col, solver.components, solver.detached_components, solver.branch_to_comp, false, nw=nw)
    end

    # update fallback solver as well and returns its info
    return apply_contingency!(solver.fallback_solver, pm, branches, master_var_vals, int_sol, nw=nw)
end

"reverts the removal of given (contingency) branches"
function reset_contingency!(solver::LinEqSysBBSolver, pm::AbstractPowerModel, branches::Vector{Int64}, old_rhss::Dict{Int64, Dict{Symbol, Float64}}; nw::Int64=nw_id_default)
    for branch_id in branches
        branch = ref(pm, nw, :branch, branch_id)
        solver.A, solver.power_surplus, _ = update_outage_line!(solver.A, pm, solver.power_surplus, branch_id, solver.bus_to_col, solver.components, solver.detached_components, solver.branch_to_comp, true, nw=nw)

        if solver.bus_to_col[branch["f_bus"]] ≤ 0 || solver.bus_to_col[branch["t_bus"]] ≤ 0
            continue
        end
        
        add_branch!(solver.A, pm, branch_id, solver.bus_to_col, nw=nw)

        # apply values of outaged transformer again
        if branch["transformer"]
            apply_trafo_val!(solver.power_surplus, branch, solver.bus_to_col, get_transformer_shift(solver, branch))   
        end
    end

    # update fallback solver as well 
    reset_contingency!(solver.fallback_solver, pm, branches, old_rhss, nw=nw)
end

"updates the dual problem to fit the current master values, contingency & primal solution for a feasibility cut !after the primal was solved!"
function update_dual_feas!(solver::LinEqSysBBSolver, pm::AbstractPowerModel, master_var_vals::Dict{Symbol, AbstractVector}, outage_branches::Vector{Int64}, betas::Dict{Int64, Float64}; nw::Int64=nw_id_default)
    ## we are constructing the leq (  A    1    | 0 )
    ##                             (  A   -1    | 0 )  
    ##                             (  dᵀ vadᵐᵃˣ | 1 )
    ## where the 1 and -1 are at the positions of fbus and tbus off a branch with vad violation, and dᵀ is the objective function coefficient of the dual subproblem

    # we are solving something like Ax = 0
    solver.dual_rhs = zeros(Float64, length(solver.power_surplus))
    
    # add obj func > 0 constraint (we need a positive ray)
    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{Float64}()

    for col_id in unique(solver.bus_to_col)
        col_id ≤ 0 && continue

        obj_coeff = solver.power_surplus[col_id]
        if obj_coeff != 0.0
            push!(I, 1)
            push!(J, col_id)
            push!(V, obj_coeff)
        end
    end
    
    # find a branch with vad violation
    solver.A_feas_cut = deepcopy(solver.A)
    m, n = size(solver.A_feas_cut)
    for (branch_id, vad) in betas
        branch = solver.branch_data[branch_id]
        
        # add col to A
        solver.A_feas_cut = [solver.A_feas_cut build_feas_col(solver, branch, m)]

        # get diag for A
        n += 1
        push!(I, 1)
        push!(J, n)
        push!(V, get_feas_diag(master_var_vals, outage_branches, branch_id, branch, vad))
    end

    solver.A_feas_cut = [solver.A_feas_cut ; sparse(I,J,V, 1, size(solver.A_feas_cut)[2])]
    push!(solver.dual_rhs, 6.66)
end

"updates dualss with the solver's current dual values"
function update_dualss_opt!(dualss::Dict{Symbol, Dict{Int, Float64}}, pm_sub::AbstractPowerModel, master_var_vals::Dict{Symbol, AbstractVector}, outage_branches::Vector{Int64}, solver::LinEqSysBBSolver, betas::Dict{Int64, Float64}; nw::Int64=nw_id_default)
    master_var_vals_branch = master_var_vals[:z_branch]::JuMP.Containers.DenseAxisArray{Float64, 1, Tuple{Vector{Int64}}, Tuple{JuMP.Containers._AxisLookup{Dict{Int64, Int64}}}}

    solver_duals = solver.duals::Vector{Float64}
    normalize_dual_solution!(solver_duals, false)
    
    # set δ
    for bus_id in eachindex(dualss[:power_balance_lb])
        # ub if the value is positive, lb else
        b_col = solver.bus_to_col[bus_id]
        if b_col > 0
            if solver_duals[b_col] ≥ 0
                dualss[:power_balance_ub][bus_id] = solver_duals[b_col]
            else
                dualss[:power_balance_lb][bus_id] = solver_duals[b_col]
            end
        end
    end

    _update_detached_bus_duals!(solver, dualss, master_var_vals, outage_branches)

    # compute & set μ (ohm_fr_on_off_lb/ub), π (rating_from_slack_on_off_lb/ub), β (rating_from_slack_lb/ub)
    for (branch_id, branch) ∈ solver.branch_data

        β = get(betas, branch_id, 0.0)

        f_bus::Int64 = branch[:f_bus]
        t_bus::Int64 = branch[:t_bus]

        # for μ₍ᵤᵥ₎, π₍ᵤᵥ₎, the value is -β₍ᵤᵥ₎ + (-δᵤ + δᵥ)
        dual_val = -β
        f_col, t_col = solver.bus_to_col[f_bus], solver.bus_to_col[t_bus]
        if f_col > 0 && t_col > 0 # otherwise one of the buses is detached and should be ignored here
            dual_val += - solver_duals[f_col] + solver_duals[t_col]
        end

        # detect the correct variable to use: 
        #   μ if the branch is off, π else
        #     lb if the value is positive, ub else
        switchable = is_branch_switchable(branch)
        electrical = is_branch_electrical(branch)
        monitored  = is_branch_monitored(branch)
        if branch_id ∉ outage_branches
            if electrical && switchable && master_var_vals_branch[branch_id]::Float64 > SOLVER_PRECISION::Float64 # switchable and on
                if dual_val ≥ 0
                    dualss[:ohm_fr_on_off_lb][branch_id] = dual_val
                else
                    dualss[:ohm_fr_on_off_ub][branch_id] = dual_val
                end
            elseif electrical && !switchable && branch[:status]==1 # not switchable
                if dual_val ≥ 0
                    dualss[:ohm_fr_lb][branch_id] = dual_val
                else
                    dualss[:ohm_fr_ub][branch_id] = dual_val
                end
            elseif switchable # and off
                switch_val = master_var_vals_branch[branch_id]::Float64

                # turned off non electrical branches have the difference in power balances
                if dual_val == 0 && switch_val ≤ SOLVER_PRECISION::Float64 && !electrical
                    dual_val = get(dualss[:power_balance_lb], t_bus, 0.0) + get(dualss[:power_balance_ub], t_bus, 0.0) - get(dualss[:power_balance_lb], f_bus, 0.0) - get(dualss[:power_balance_ub], f_bus, 0.0)
                end

                if dual_val ≥ 0
                    dualss[:rating_from_slack_on_off_lb][branch_id] = dual_val
                else
                    dualss[:rating_from_slack_on_off_ub][branch_id] = dual_val
                end
            end
        end

        # lb if the value is positive, ub else
        if monitored
            if β ≥ 0
                dualss[:rating_from_slack_lb][branch_id] = β
            else
                dualss[:rating_from_slack_ub][branch_id] = β
            end
        end
    end

    _update_bb_component_duals!(solver, dualss)
end

"updates dualss with the solver's current dual values"
function update_dualss_feas!(dualss::Dict{Symbol, Dict{Int, Float64}}, pm_sub::AbstractPowerModel, master_var_vals::Dict{Symbol, AbstractVector}, outage_branches::Vector{Int64}, solver::LinEqSysBBSolver, betas::Dict{Int64, Float64}; nw::Int64=nw_id_default)

    # avoid gc
    solver_duals = solver.duals::Vector{Float64}
    bus_to_col = solver.bus_to_col::Vector{Int64}
    master_var_vals_branch = master_var_vals[:z_branch]::JuMP.Containers.DenseAxisArray{Float64, 1, Tuple{Vector{Int64}}, Tuple{JuMP.Containers._AxisLookup{Dict{Int64, Int64}}}}

    normalize_dual_solution!(solver_duals, true)

    # set δ
    for bus_id in eachindex(dualss[:power_balance_lb])
        col = bus_to_col[bus_id]
        if col > 0
            # ub if the value is positive, lb else
            if solver_duals[col] ≥ SOLVER_PRECISION/1e4
                dualss[:power_balance_ub][bus_id] = solver_duals[col]
            elseif solver_duals[col] ≤ -SOLVER_PRECISION/1e4
                dualss[:power_balance_lb][bus_id] = solver_duals[col]
            end
        end
    end

    _update_detached_bus_duals!(solver, dualss, master_var_vals, outage_branches)

    # compute & set μ (ohm_fr_on_off_lb/ub), π (rating_from_slack_on_off_lb/ub), β (rating_from_slack_lb/ub)
    for (branch_id, branch::Dict{String, Any}) ∈ ref(pm_sub, nw, :branch)::Dict{Int64, Any}

        f_bus = branch["f_bus"]::Int64
        t_bus = branch["t_bus"]::Int64

        # for μ₍ᵤᵥ₎, π₍ᵤᵥ₎, the value is -β₍ᵤᵥ₎ + (-δᵤ + δᵥ)
        if bus_to_col[f_bus] ≤ 0 || bus_to_col[t_bus] ≤ 0
            dual_val = 0
        else
            dual_val = - solver_duals[bus_to_col[f_bus]] + solver_duals[bus_to_col[t_bus]]
        end

        # detect the correct variable to use: 
        #   μ if the branch is off, π else
        #     lb if the value is positive, ub else
        switchable = is_branch_switchable(branch)
        electrical = is_branch_electrical(branch)
        monitored  = is_branch_monitored(branch)
        if branch_id ∉ outage_branches
            switch_val = switchable ? master_var_vals_branch[branch_id]::Float64 : 0.0

            if switchable
                if dual_val == 0 && switch_val ≤ SOLVER_PRECISION::Float64 && !electrical
                    dual_val = get(dualss[:power_balance_lb], t_bus, 0.0) + get(dualss[:power_balance_ub], t_bus, 0.0) - get(dualss[:power_balance_lb], f_bus, 0.0) - get(dualss[:power_balance_ub], f_bus, 0.0)
                end
            end
            
            # skip very small values
            if -SOLVER_PRECISION/1e4 < dual_val && dual_val < SOLVER_PRECISION/1e4
                continue
            end

            if electrical && switchable && master_var_vals_branch[branch_id]::Float64 > SOLVER_PRECISION::Float64 # switchable and on
                if dual_val ≥ 0
                    dualss[:ohm_fr_on_off_lb][branch_id] = dual_val
                else
                    dualss[:ohm_fr_on_off_ub][branch_id] = dual_val
                end
            elseif electrical && !switchable # not switchable
                if dual_val ≥ 0
                    dualss[:ohm_fr_lb][branch_id] = dual_val
                else
                    dualss[:ohm_fr_ub][branch_id] = dual_val
                end
            elseif switchable # and off
                
                # turned off non electrical branches have the difference in power balances
                if dual_val == 0 && switch_val ≤ SOLVER_PRECISION::Float64 && !electrical
                    dual_val = get(dualss[:power_balance_lb], t_bus, 0.0) + get(dualss[:power_balance_ub], t_bus, 0.0) - get(dualss[:power_balance_lb], f_bus, 0.0) - get(dualss[:power_balance_ub], f_bus, 0.0)
                end

                if dual_val ≥ 0
                    dualss[:rating_from_slack_on_off_lb][branch_id] = dual_val
                else
                    dualss[:rating_from_slack_on_off_ub][branch_id] = dual_val
                end
            end
        end
    end

    # set χ
    for (branch_id, vad) in betas
        branch = solver.branch_data[branch_id]
        if branch_id ∉ outage_branches 
            if is_branch_switchable(branch)
                if master_var_vals_branch[branch_id]::Float64 > 0.5
                    if vad < branch[:angmin] - SOLVER_PRECISION
                        dualss[:vad_on_off_lb][branch_id] = solver.duals[end]
                    else
                        dualss[:vad_on_off_ub][branch_id] = solver.duals[end]
                    end
                else
                    if vad < branch[:angmin_close] - SOLVER_PRECISION
                        dualss[:vad_on_off_lb][branch_id] = -solver.duals[end]
                    else
                        dualss[:vad_on_off_ub][branch_id] = -solver.duals[end]
                    end
                end
            else
                if vad < branch[:angmin] - SOLVER_PRECISION
                    dualss[:vad_lb][branch_id] = solver.duals[end]
                else
                    dualss[:vad_ub][branch_id] = solver.duals[end]
                end
            end
        end
    end

    _update_bb_component_duals!(solver, dualss)
end