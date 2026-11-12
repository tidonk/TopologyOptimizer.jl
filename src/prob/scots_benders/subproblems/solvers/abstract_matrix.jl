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

## general AbstractMatrixSolver
"returns calculated value of vad and pf at the branch"
function get_vad(solver::AbstractMatrixSolver, branch::Dict)::Float64
    return solver.vas[solver.bus_to_col[branch["f_bus"]]] - solver.vas[solver.bus_to_col[branch["t_bus"]]]
end
function get_pf(solver::AbstractMatrixSolver, branch::Dict, vad::Float64)::Float64
    return -branch["b"] * (vad - get_transformer_shift(solver, branch))
end
function get_vad(solver::AbstractMatrixSolver, branch::Dict{Symbol, Any})::Float64
    return solver.vas[solver.bus_to_col[branch[:f_bus]]] - solver.vas[solver.bus_to_col[branch[:t_bus]]]
end
function get_vad(solver::AbstractMatrixSolver, f_bus::Int64, t_bus::Int64)::Float64
    return solver.vas[solver.bus_to_col[f_bus]] - solver.vas[solver.bus_to_col[t_bus]]
end
function get_pf(solver::AbstractMatrixSolver, branch::Dict{Symbol, Any}, vad::Float64)::Float64
    b::Float64 = branch[:b]
    return -b * (vad - get_transformer_shift(solver, branch))
end
function get_vad_pf(solver::AbstractMatrixSolver, branch::Dict; nw::Int64=nw_id_default)::Tuple{Float64, Float64}
    vad = get_vad(solver, branch)
    return (vad, get_pf(solver, branch, vad))
end

"applies the current master solution to solver"
function apply_master!(solver::AbstractMatrixSolver, pm::AbstractPowerModel, master_var_vals::Dict{Symbol, AbstractVector}; nw::Int64=nw_id_default)
    add_master_trafo_vals!(solver, master_var_vals)
    
    solver.A, solver.power_surplus = apply_master_vals!(solver.A, pm, solver.power_surplus, solver.bus_to_col, master_var_vals, nw=nw)
end

"resets the solver back to its base state"
function reset_master!(solver::AbstractMatrixSolver, pm::AbstractPowerModel, master_var_vals::Dict{Symbol, AbstractVector}; nw::Int64=nw_id_default)
    solver.A, solver.power_surplus = reset_master_vals!(solver.A, pm, solver.power_surplus, solver.bus_to_col, master_var_vals, nw=nw)
end

"applies one contingency to the solver"
function apply_contingency!(solver::AbstractMatrixSolver, pm::AbstractPowerModel, branches::Vector{Int64}, master_var_vals::Dict{Symbol, AbstractVector}, int_sol::Bool; nw::Int64=nw_id_default)::Dict{Int64, Dict{Symbol, Float64}}
    @assert(int_sol, "Only integer solutions are supported for the matrix solver")
    for branch_id in branches
        remove_branch!(solver.A, pm, branch_id, solver.bus_to_col, nw=nw)
        branch = ref(pm, nw, :branch, branch_id)
        reset_trafo_val!(solver.power_surplus, branch, solver.bus_to_col, get_transformer_shift(solver, branch))
    end
    
    # update fallback solver as well and returns its info
    return apply_contingency!(solver.fallback_solver, pm, branches, master_var_vals, int_sol, nw=nw)
end

"reverts the removal of given (contingency) branches"
function reset_contingency!(solver::AbstractMatrixSolver, pm::AbstractPowerModel, branches::Vector{Int64}, old_rhss::Dict{Int64, Dict{Symbol, Float64}}; nw::Int64=nw_id_default)
    for branch_id in branches
        add_branch!(solver.A, pm, branch_id, solver.bus_to_col, nw=nw)
        branch = ref(pm, nw, :branch, branch_id)
        apply_trafo_val!(solver.power_surplus, branch, solver.bus_to_col, get_transformer_shift(solver, branch))
    end

    # update fallback solver as well
    reset_contingency!(solver.fallback_solver, pm, branches, old_rhss, nw=nw)
end

"updates the dual problem to fit the current master values, contingency & primal solution for an optimality cut !after the primal was solved!"
function update_dual_opt!(solver::AbstractMatrixSolver, pm::AbstractPowerModel, master_var_vals::Dict{Symbol, AbstractVector}, betas::Dict{Int64, Float64}; nw::Int64=nw_id_default)
    if typeof(solver) == LinEqSysBBSolver || typeof(solver) == LinEqSysBBRfSolver || length(solver.dual_rhs) ≠ length(solver.power_surplus)
        solver.dual_rhs = zeros(Float64, length(solver.power_surplus))
    else
        fill!(solver.dual_rhs, 0.0)
    end

    for (branch_id, _) in betas
        branch = solver.branch_data[branch_id]
        _, pf = get_vad_pf(solver, branch, nw=nw)


        ## fill rhs for branches with positive slack
        # we need to add to the rhs of node u
        # for each outgoing branch (u,v) from u
        #  pf>0 -> - cˢˡᵃᶜᵏ₍ᵤᵥ₎ * B₍ᵤᵥ₎
        #  pf<0 ->   cˢˡᵃᶜᵏ₍ᵤᵥ₎ * B₍ᵤᵥ₎
        # for each incoming branch (v,u) to   u
        #  pf>0 ->   cˢˡᵃᶜᵏ₍ᵥᵤ₎ * B₍ᵥᵤ₎
        #  pf<0 -> - cˢˡᵃᶜᵏ₍ᵥᵤ₎ * B₍ᵥᵤ₎
        # we are iterating over branches, so instead we have
        #  pf>0 -> rhsᵤ += - cˢˡᵃᶜᵏ₍ᵤᵥ₎ * B₍ᵤᵥ₎
        #          rhsᵥ +=   cˢˡᵃᶜᵏ₍ᵤᵥ₎ * B₍ᵤᵥ₎
        #  pf<0 -> rhsᵤ +=   cˢˡᵃᶜᵏ₍ᵤᵥ₎ * B₍ᵤᵥ₎
        #          rhsᵥ += - cˢˡᵃᶜᵏ₍ᵤᵥ₎ * B₍ᵤᵥ₎

        #                         B₍ᵤᵥ₎      *         cˢˡᵃᶜᵏ₍ᵤᵥ₎
        value = (pf > 0 ? -1 : 1)*branch[:b]*branch[:obj_coeff]

        f_col = solver.bus_to_col[branch[:f_bus]]
        t_col = solver.bus_to_col[branch[:t_bus]]

        (f_col <= 0 || t_col <= 0) && continue

        solver.dual_rhs[f_col] += value
        solver.dual_rhs[t_col] -= value
    end
end

"updates the dual problem to fit the current master values, contingency & primal solution for a feasibility cut !after the primal was solved!"
function update_dual_feas!(solver::AbstractMatrixSolver, pm::AbstractPowerModel, master_var_vals::Dict{Symbol, AbstractVector}, outage_branches::Vector{Int64}, betas::Dict{Int64, Float64}; nw::Int64=nw_id_default)
    ## we are constructing the leq (  A    1    )   ( 0 )
    ##                             (  |   -1    ) = ( 0 )  
    ##                             (  dᵀ vadᵐᵃˣ )   ( 1 )
    ## where the 1 and -1 are at the positions of fbus and tbus of a branch with vad violation (the rest of the column is 0), and dᵀ is the objective function coefficient of the dual subproblem

    # we are solving something like Ax = 0
    if typeof(solver) == LinEqSysBBSolver || typeof(solver) == LinEqSysBBRfSolver || length(solver.dual_rhs) ≠ length(solver.power_surplus)
        solver.dual_rhs = zeros(Float64, length(solver.power_surplus))
    else
        fill!(solver.dual_rhs, 0.0)
    end
    
    # add obj func constraint (we need a positive ray)
    I = Vector{Int64}()   # row index
    J = Vector{Int64}()   # col index
    V = Vector{Float64}() # value
    for bus_id in 1:solver.n_buses
        bus_col = solver.bus_to_col[bus_id]
        bus_col <= 0 && continue
        
        obj_coeff = solver.power_surplus[bus_col]
        if obj_coeff != 0.0
            push!(I, 1)
            push!(J, bus_col)
            push!(V, obj_coeff)
        end
    end
    
    # find a branch with vad violation
    solver.A_feas_cut = deepcopy(solver.A)
    for (branch_id, vad) in betas
        branch = solver.branch_data[branch_id]

        f_col = solver.bus_to_col[branch[:f_bus]]
        t_col = solver.bus_to_col[branch[:t_bus]]

        (f_col <= 0 || t_col <= 0) && continue
        
        # add col to A
        a = sparse([f_col, t_col], [1, 1], [1.0, -1.0], solver.n_buses, 1)
        solver.A_feas_cut = [solver.A_feas_cut a]

        is_off = is_branch_switchable(branch) ? (master_var_vals[:z_branch][branch_id]::Float64 < 0.5) : (branch[:status]::Int == 0)

        if branch_id ∈ outage_branches || is_off
            angmin = branch[:angmin_close]::Float64
            angmax = branch[:angmax_close]::Float64
        else
            angmin = branch[:angmin]::Float64
            angmax = branch[:angmax]::Float64
        end

        push!(I, 1)
        push!(J, size(solver.A_feas_cut)[2])
        if vad < angmin - SOLVER_PRECISION::Float64
            push!(V, angmin)
        else
            push!(V, angmax)
        end
    end

    solver.A_feas_cut = [solver.A_feas_cut ; sparse(I,J,V, 1, size(solver.A_feas_cut)[2])]
    push!(solver.dual_rhs, 6.66)
end

"updates dualss with the solver's current dual values"
function update_dualss_opt!(dualss::Dict{Symbol, Dict{Int, Float64}}, pm_sub::AbstractPowerModel, master_var_vals::Dict{Symbol, AbstractVector}, outage_branches::Vector{Int64}, solver::AbstractMatrixSolver, betas::Dict{Int64, Float64}; nw::Int64=nw_id_default)

    # avoid gc
    solver_duals = solver.duals::Vector{Float64}
    bus_to_col = solver.bus_to_col::Vector{Int64}
    master_var_vals_branch = master_var_vals[:z_branch]::JuMP.Containers.DenseAxisArray{Float64, 1, Tuple{Vector{Int64}}, Tuple{JuMP.Containers._AxisLookup{Dict{Int64, Int64}}}}

    # compute & set μ (ohm_fr_on_off_lb/ub), π (rating_from_slack_on_off_lb/ub), β (rating_from_slack_lb/ub)
    for (branch_id, branch::Dict{Symbol, Any}) ∈ solver.branch_data

        β = get(betas, branch_id, 0.0)

        f_bus = branch[:f_bus]::Int64
        t_bus = branch[:t_bus]::Int64

        # for μ₍ᵤᵥ₎, π₍ᵤᵥ₎, the value is -β₍ᵤᵥ₎ + (-δᵤ + δᵥ)
        dual_val = -β - solver_duals[bus_to_col[f_bus]] + solver_duals[bus_to_col[t_bus]]

        # detect the correct variable to use: 
        #   μ if the branch is off, π else
        #     lb if the value is positive, ub else
        if is_branch_switchable(branch)
            switch_val = master_var_vals_branch[branch_id]::Float64
            if branch_id ∉ outage_branches && switch_val > SOLVER_PRECISION::Float64
                if dual_val >= 0.0
                    dualss[:ohm_fr_on_off_lb][branch_id] = dual_val
                else
                    dualss[:ohm_fr_on_off_ub][branch_id] = dual_val
                end
            else
                if dual_val >= 0.0
                    dualss[:rating_from_slack_on_off_lb][branch_id] = dual_val
                else
                    dualss[:rating_from_slack_on_off_ub][branch_id] = dual_val
                end
            end
        elseif branch_id ∉ outage_branches && branch[:status] == 1
            # Fixed active lines still have Ohm constraints. Their RHS can
            # contribute to the cut when a transformer has a phase shift.
            dualss[dual_val >= 0.0 ? :ohm_fr_lb : :ohm_fr_ub][branch_id] = dual_val
        end

        # lb if the value is positive, ub else
        if β >= 0
            dualss[:rating_from_slack_lb][branch_id] = β
        else
            dualss[:rating_from_slack_ub][branch_id] = β
        end
    end

    # set δ
    for bus_id in eachindex(dualss[:power_balance_lb])
        # ub if the value is positive, lb else
        if solver_duals[bus_to_col[bus_id]] >= 0.0
            dualss[:power_balance_ub][bus_id] = solver_duals[bus_to_col[bus_id]]
        else
            dualss[:power_balance_lb][bus_id] = solver_duals[bus_to_col[bus_id]]
        end
    end
end

"updates dualss with the solver's current dual values"
function update_dualss_feas!(dualss::Dict{Symbol, Dict{Int, Float64}}, pm_sub::AbstractPowerModel, master_var_vals::Dict{Symbol, AbstractVector}, outage_branches::Vector{Int64}, solver::AbstractMatrixSolver, betas::Dict{Int64, Float64}; nw::Int64=nw_id_default)

    # avoid gc
    solver_duals = solver.duals::Vector{Float64}
    bus_to_col = solver.bus_to_col::Vector{Int64}
    master_var_vals_branch = master_var_vals[:z_branch]::JuMP.Containers.DenseAxisArray{Float64, 1, Tuple{Vector{Int64}}, Tuple{JuMP.Containers._AxisLookup{Dict{Int64, Int64}}}}

    # normalize va values
    dual_val = solver_duals[end-1]
    for other_bus_id in 1:solver.n_buses
        solver_duals[bus_to_col[other_bus_id]] -= dual_val
    end

    # compute & set μ (ohm_fr_on_off_lb/ub), π (rating_from_slack_on_off_lb/ub), β (rating_from_slack_lb/ub)
    for (branch_id, branch::Dict{Symbol, Any}) ∈ solver.branch_data
        f_bus = branch[:f_bus]::Int64
        t_bus = branch[:t_bus]::Int64

        # for μ₍ᵤᵥ₎, π₍ᵤᵥ₎, the value is -β₍ᵤᵥ₎ + (-δᵤ + δᵥ)
        dual_val = - solver_duals[bus_to_col[f_bus]] + solver_duals[bus_to_col[t_bus]]

        # skip very small values
        if -SOLVER_PRECISION < dual_val && dual_val < SOLVER_PRECISION
            continue
        end

        # detect the correct variable to use: 
        #   μ if the branch is off, π else
        #     lb if the value is positive, ub else
        if is_branch_switchable(branch)
            switch_val = master_var_vals_branch[branch_id]::Float64
            if branch_id ∉ outage_branches && switch_val > SOLVER_PRECISION::Float64
                if dual_val >= 0.0
                    dualss[:ohm_fr_on_off_lb][branch_id] = dual_val
                else
                    dualss[:ohm_fr_on_off_ub][branch_id] = dual_val
                end
            else
                if dual_val >= 0.0
                    dualss[:rating_from_slack_on_off_lb][branch_id] = dual_val
                else
                    dualss[:rating_from_slack_on_off_ub][branch_id] = dual_val
                end
            end
        elseif branch_id ∉ outage_branches && branch[:status] == 1
            # Fixed active lines still have Ohm constraints. Their RHS can
            # contribute to the cut when a transformer has a phase shift.
            dualss[dual_val >= 0.0 ? :ohm_fr_lb : :ohm_fr_ub][branch_id] = dual_val
        end
    end

    # set δ
    for bus_id in eachindex(dualss[:power_balance_lb])
        # ub if the value is positive, lb else
        if solver_duals[bus_to_col[bus_id]] >= SOLVER_PRECISION/1e4
            dualss[:power_balance_ub][bus_id] = solver_duals[bus_to_col[bus_id]]
        elseif solver_duals[bus_to_col[bus_id]] <= -SOLVER_PRECISION/1e4
            dualss[:power_balance_lb][bus_id] = solver_duals[bus_to_col[bus_id]]
        end
    end

    # set χ
    for (branch_id, vad) in betas
        branch = solver.branch_data[branch_id]
        if is_branch_switchable(branch)
            if branch_id ∉ outage_branches 
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

            end
        elseif branch_id ∉ outage_branches && branch[:status] == 1
            # Connectedness preprocessing can fix this line on. Omitting its
            # angle-bound multiplier would turn the Farkas ray into an invalid cut.
            key = vad < branch[:angmin] - SOLVER_PRECISION ? :vad_lb : :vad_ub
            dualss[key][branch_id] = solver.duals[end]
        end
    end
end

"updates dualss with the solver's current dual values, for when the dual solution has been obtained via early solve (no separate dual solve)"
function update_dualss_feas_early_solve!(dualss::Dict{Symbol, Dict{Int, Float64}}, solver::AbstractMatrixSolver, master_var_vals::Dict{Symbol, AbstractVector}, outage_branches::Vector{Int64}, betas::Dict{Int64, Float64})
    for (branch_id, _) = betas
        branch = solver.branch_data[branch_id]

        f_bus_id = branch[:f_bus]
        t_bus_id = branch[:t_bus]

        disc_bus_id = abs(solver.vas[solver.bus_to_col[branch[:f_bus]]]) > abs(solver.vas[solver.bus_to_col[branch[:t_bus]]]) ? f_bus_id : t_bus_id

        dualss[:power_balance_lb][disc_bus_id]          = -1.0
        dualss[:rating_from_slack_on_off_ub][branch_id] = -1.0

        # we need to find all other branches connected to this one and set `dualss[:power_balance_lb][other_branch[:*_bus]] = -1.0; and `dualss[:rating_from_slack_on_off_lb][other_branch_id] = 1.0` for all turned off but not outaged branches
        visited_branches = [branch_id]
        buses_to_visit = [disc_bus_id]

        pos = 1
        while pos ≤ length(buses_to_visit)
            cur_bus_id = buses_to_visit[pos]

            # search through all connected on branches
            for (other_branch_id, other_branch) in solver.branch_data
                if other_branch_id ∉ visited_branches

                    is_on = is_branch_switchable(other_branch) ? (master_var_vals[:z_branch][other_branch_id]::Float64 > 0.5) : (other_branch[:status]::Int == 1)
                    if is_on

                        if other_branch[:f_bus] == cur_bus_id
                            push!(buses_to_visit, other_branch[:t_bus])
                            dualss[:power_balance_lb][other_branch[:t_bus]] = -1.0
                        elseif other_branch[:t_bus] == cur_bus_id
                            push!(buses_to_visit, other_branch[:f_bus])
                            dualss[:power_balance_lb][other_branch[:f_bus]] = -1.0
                        else
                            continue
                        end

                    elseif other_branch_id ∉ outage_branches
                        if other_branch[:f_bus] == cur_bus_id || other_branch[:t_bus] == cur_bus_id

                            if haskey(dualss[:rating_from_slack_on_off_lb], other_branch_id)
                                dualss[:rating_from_slack_on_off_lb][other_branch_id] = 1.0
                            end

                        else
                            continue
                        end

                    end

                    push!(visited_branches, other_branch_id) # this includes off && outaged branches
                end
            end

            pos += 1
        end

        break
    end
end

"check dual solution of a solver for inconsistencies"
function check_dual_solution(solver::AbstractMatrixSolver)
    if PARAMS[:BENDERS][:CHECKS][:CHECK_QUESTIONABLE_SOLS]
        # if a majority of the dual variables have the same values, something probably went wrong => solve the LP
        for (_, val) in countmap(round.(solver.duals, digits=SOLVER_PRECISION_DIGITS-1))
            if val >= 0.75*solver.n_buses
                @warn("Subproblem solver returned questionable solution (questionability score: $(round(val/solver.n_buses * 100, digits=2)) %). Checking with fallback solver.")
                callback_stats[:n_questionable_sol] += 1
                return false
            end
        end
    end
    return true
end

"normalizes the dual solution using an arbitrary element of the solution vector"
function normalize_dual_solution!(solver_duals::Vector{Float64}, is_feas::Bool)
    # omit the very last entry if this is a feasibility cut, its value is fixed
    last_pos = length(solver_duals)
    if is_feas
        last_pos -= 1
    else
        return
    end

    ref_row_id = argmin(solver_duals[1:last_pos])

    dual_val = solver_duals[ref_row_id]
    for dual_row_id in 1:last_pos
        solver_duals[dual_row_id] -= dual_val + 1.0
    end
end
