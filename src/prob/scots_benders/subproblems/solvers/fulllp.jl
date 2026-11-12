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

## Full LP solver
"returns calculated value of vad and pf at the branch"
function get_vad_pf(solver::FullLPSolver, branch::Dict; nw::Int64=nw_id_default)::Tuple{Float64, Float64}
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (branch["index"], f_bus, t_bus)

    pf  = JuMP.value(var(solver.pm, nw, :p,  f_idx))
    vad = JuMP.value(var(solver.pm, nw, :va, f_bus)) - JuMP.value(var(solver.pm, nw, :va, t_bus))

    return (vad, pf)
end
function get_vad_pf(solver::FullLPSolver, branch::Dict{Symbol, Any}; nw::Int64=nw_id_default)::Tuple{Float64, Float64}
    f_bus = branch[:f_bus]
    t_bus = branch[:t_bus]
    f_idx = (branch[:index], f_bus, t_bus)

    pf  = JuMP.value(var(solver.pm, nw, :p,  f_idx))
    vad = JuMP.value(var(solver.pm, nw, :va, f_bus)) - JuMP.value(var(solver.pm, nw, :va, t_bus))

    return (vad, pf)
end


"returns calculated value of va at the bus"
function get_va(solver::FullLPSolver, bus_id::Int64; nw::Int64=nw_id_default)::Float64
    return JuMP.value(var(solver.pm, nw, :va, bus_id))
end

"intializes the solvers"
function init_fulllp_solvers(pms, n_threads::Int64, sub_cons_for_dualss::Vector{Dict{Symbol, Dict{Int, JuMP.ConstraintRef}}})::Vector{FullLPSolver}
    return [FullLPSolver(pm=pms[i], sub_cons_for_duals=sub_cons_for_dualss[i]) for i in 1:n_threads]
end

"solves problem in its current state, returns success flag"
function solve_primal!(solver::FullLPSolver)::Bool
    JuMP.optimize!(solver.pm.model)

    # for debugging - infeasible subproblems might hint at data errors
    if PARAMS[:BENDERS][:CHECKS][:CHECK_WITH_IIS] && JuMP.termination_status(solver.pm.model) == JuMP.MOI.INFEASIBLE
        @info("Optimizing IIS of Subproblem:")
        JuMP.compute_conflict!(solver.pm.model)
        iis_model, _ = JuMP.copy_conflict(solver.pm.model)
        print(iis_model)     
    end
    
    # return true even if the solver failed (i.e., if JuMP.termination_status(solver.pm.model) ∈ [JuMP.MOI.INFEASIBLE, JuMP.MOI.OPTIMAL] does not hold), because we handle deviating termination statuses later in the code
    return true
end

"solves dual problem in its current state !after the primal was solved! to obtain a optimality cut, returns success flag"
function solve_dual_opt!(solver::FullLPSolver)::Bool
    return true # if the primal was solved, we already have the dual values
end

"applies one contingency to the solver"
function apply_contingency!(solver::FullLPSolver, pm_sub::AbstractPowerModel, branches::Vector{Int64}, master_var_vals::Dict{Symbol, AbstractVector}, int_sol::Bool; nw::Int64=nw_id_default)::Dict{Int64, Dict{Symbol, Float64}}
    fix_master_variables!(solver.pm, nw, master_var_vals, int_sol)

    # prepare to store right hand sides of constraints
    tmp_rhsss = Dict{Int64, Dict{Symbol, Float64}}()

    for branch_id in branches
        branch = ref(pm_sub, nw, :branch, branch_id)

        # if branch is switchable, it has on/off constraints linked to z_branch, which was fixed in `fix_master_variables!` above
        if ! is_branch_switchable(branch)
            # prepare to save right hand sides of constraints
            tmp_rhsss[branch["index"]] = Dict{Symbol, Float64}()
            tmp_rhss = tmp_rhsss[branch["index"]]

            if is_branch_electrical(branch)
                # get thermal limit constraints (only exist for electrical branches) and store rhss, then relax constraints by setting rhs to -/+ "large"
                ohm_fr_ub = solver.sub_cons_for_duals[:ohm_fr_ub][branch["index"]]
                ohm_fr_lb = solver.sub_cons_for_duals[:ohm_fr_lb][branch["index"]]

                tmp_rhss[:ohm_fr_lb] = JuMP.normalized_rhs(ohm_fr_ub)
                tmp_rhss[:ohm_fr_ub] = JuMP.normalized_rhs(ohm_fr_lb)

                JuMP.set_normalized_rhs(ohm_fr_lb, -100000.0)
                JuMP.set_normalized_rhs(ohm_fr_ub,  100000.0)
            end

            # get voltage angle difference constraints and store rhss, then relax constraints by setting rhs to -/+ "large"
            vad_lb    = solver.sub_cons_for_duals[:vad_lb][branch["index"]]
            vad_ub    = solver.sub_cons_for_duals[:vad_ub][branch["index"]]

            tmp_rhss[:vad_lb]    = JuMP.normalized_rhs(vad_lb)
            tmp_rhss[:vad_ub]    = JuMP.normalized_rhs(vad_ub)

            JuMP.set_normalized_rhs(vad_lb,    -100000.0)
            JuMP.set_normalized_rhs(vad_ub,     100000.0)

            # fix p=0
            JuMP.fix(var(pm_sub, nw, :p, (branch["index"], branch["f_bus"], branch["t_bus"])), 0.0, force=true)
        end
    end

    return tmp_rhsss
end

"reverts the removal of constraints of given (contingency) branches"
function reset_contingency!(solver::FullLPSolver, pm_sub::AbstractPowerModel, branches::Vector{Int64}, tmp_rhsss::Dict{Int64, Dict{Symbol, Float64}}; nw::Int64=nw_id_default)
    for branch_id in branches
        branch = ref(pm_sub, nw, :branch, branch_id)

        if ! is_branch_switchable(branch)
            # get stored rhss
            tmp_rhss = tmp_rhsss[branch["index"]]

            if is_branch_electrical(branch)
                # get thermal limit constraints (only exist for electrical branches) and apply stored rhss
                ohm_fr_ub = solver.sub_cons_for_duals[:ohm_fr_ub][branch["index"]]
                ohm_fr_lb = solver.sub_cons_for_duals[:ohm_fr_lb][branch["index"]]

                JuMP.set_normalized_rhs(ohm_fr_lb, tmp_rhss[:ohm_fr_lb])
                JuMP.set_normalized_rhs(ohm_fr_ub, tmp_rhss[:ohm_fr_ub]) 
            end

            # get voltage angle difference constraints and apply stored rhss
            vad_lb    = solver.sub_cons_for_duals[:vad_lb][branch["index"]]
            vad_ub    = solver.sub_cons_for_duals[:vad_ub][branch["index"]]

            JuMP.set_normalized_rhs(vad_lb,    tmp_rhss[:vad_lb])
            JuMP.set_normalized_rhs(vad_ub,    tmp_rhss[:vad_ub])

            # unfix p
            JuMP.unfix(var(pm_sub, nw, :p, (branch["index"], branch["f_bus"], branch["t_bus"])))
        end
    end
end

"updates dualss with the solver's current dual values"
function update_dualss_opt!(dualss::Dict{Symbol, Dict{Int, Float64}}, pm_sub::AbstractPowerModel, master_var_vals::Dict{Symbol, AbstractVector}, outage_branches::Vector{Int64}, solver::FullLPSolver, betas::Dict{Int64, Float64}; nw::Int64=nw_id_default)
    for (con_symbol, cons_dict) in solver.sub_cons_for_duals
        for (comp_id, con) in cons_dict
            dualss[con_symbol][comp_id] = JuMP.dual(con)
        end
    end
end
function update_dualss_feas!(dualss::Dict{Symbol, Dict{Int, Float64}}, pm_sub::AbstractPowerModel, master_var_vals::Dict{Symbol, AbstractVector}, outage_branches::Vector{Int64}, solver::FullLPSolver, betas::Dict{Int64, Float64}; nw::Int64=nw_id_default)
    update_dualss_opt!(dualss, pm_sub, master_var_vals, outage_branches, solver, betas, nw=nw)
end
function update_dualss_opt!(dualss::Dict{Symbol, Dict{Int, Float64}}, sub_cons_for_duals_gen_model)
    for (con_symbol, cons_dict) in sub_cons_for_duals_gen_model
        for (comp_id, con) in cons_dict
            dualss[con_symbol][comp_id] = JuMP.dual(con)
        end
    end
end
function update_dualss_feas!(dualss::Dict{Symbol, Dict{Int, Float64}}, sub_cons_for_duals_gen_model)
    return update_dualss_opt!(dualss, sub_cons_for_duals_gen_model)
end