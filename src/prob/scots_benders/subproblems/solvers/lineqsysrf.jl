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

## linear equation system with re-factorization
"intializes the solvers"
function init_lineqsysrf_solvers(pm, pm_subs, n_threads::Int64, sub_cons_for_dualss; nw::Int64=nw_id_default)::Vector{LinEqSysRfSolver}
    A, power_surplus, bus_to_col = build_base_pf_matrix(pm, nw=nw)

    # add FullLP solver as fallback solver (one per thread)
    fallback_solver_solvers = init_fulllp_solvers(pm_subs, n_threads, sub_cons_for_dualss)
    fallback_solvers = [FallbackSolver(solver=fallback_solver_solvers[t], pm_sub=pm_subs[t], nw=nw) for t in 1:n_threads]

    branch_data = collect_branch_data(pm_subs[1], nw=nw)

    solvers = Vector{LinEqSysRfSolver}()

    v = LinearSolve.LinearVerbosity(
        solver_failure = LinearSolve.SciMLLogging.Silent()
    )

    for t in 1:n_threads
        A_cp             = t==1 ? A             : deepcopy(A)
        power_surplus_cp = t==1 ? power_surplus : deepcopy(power_surplus)
        bus_to_col_cp    = t==1 ? bus_to_col    : deepcopy(bus_to_col)
        
        linear_problem = LinearSolve.LinearProblem(A_cp, power_surplus_cp)
        linsolve = LinearSolve.init(linear_problem, verbose = v)

        push!(solvers, LinEqSysRfSolver(A=A_cp, power_surplus=power_surplus_cp, linsolve=linsolve, bus_to_col=bus_to_col_cp, fallback_solver=fallback_solvers[t], branch_data = branch_data))

    end

    return solvers
end

"applies the current master solution to solver"
function apply_master!(solver::LinEqSysRfSolver, pm::AbstractPowerModel, master_var_vals::Dict{Symbol, AbstractVector}; nw::Int64=nw_id_default)
    add_master_trafo_vals!(solver, master_var_vals)

    if solver.is_factorized
        solver.L, solver.U, solver.power_surplus = apply_master_vals!(solver.L, solver.U, pm, solver.power_surplus, solver.bus_to_col, master_var_vals, solver.p, solver.q, solver.Rs, nw=nw)
        solver.A, solver.power_surplus = apply_master_vals!(solver.A, pm, solver.power_surplus, solver.bus_to_col, master_var_vals, nw=nw, gens=false)
    else
        solver.A, solver.power_surplus = apply_master_vals!(solver.A, pm, solver.power_surplus, solver.bus_to_col, master_var_vals, nw=nw)
    end
end

"resets the solver back to its base state"
function reset_master!(solver::LinEqSysRfSolver, pm::AbstractPowerModel, master_var_vals::Dict{Symbol, AbstractVector}; nw::Int64=nw_id_default)
    if solver.is_factorized
        solver.L, solver.U, solver.power_surplus = reset_master_vals!(solver.L, solver.U, pm, solver.power_surplus, solver.bus_to_col, master_var_vals, solver.p, solver.q, solver.Rs, nw=nw)
        solver.A, solver.power_surplus = reset_master_vals!(solver.A, pm, solver.power_surplus, solver.bus_to_col, master_var_vals, nw=nw, gens=false)
    else
        solver.A, solver.power_surplus = reset_master_vals!(solver.A, pm, solver.power_surplus, solver.bus_to_col, master_var_vals, nw=nw)
    end
end

"applies one contingency to the solver"
function apply_contingency!(solver::LinEqSysRfSolver, pm::AbstractPowerModel, branches::Vector{Int64}, master_var_vals::Dict{Symbol, AbstractVector}, int_sol::Bool; nw::Int64=nw_id_default)::Dict{Int64, Dict{Symbol, Float64}}
    @assert(int_sol, "Only integer solutions are supported for the matrix solver")

    if solver.is_factorized
        for branch_id in branches
            remove_branch!(solver.L, solver.U, pm, branch_id, solver.bus_to_col, solver.p, solver.q, solver.Rs, nw=nw)
            remove_branch!(solver.A, pm, branch_id, solver.bus_to_col, nw=nw)
        end
    else
        for branch_id in branches
            remove_branch!(solver.A, pm, branch_id, solver.bus_to_col, nw=nw)
        end
    end
    
    # update fallback solver as well and returns its info
    return apply_contingency!(solver.fallback_solver, pm, branches, master_var_vals, int_sol, nw=nw)
end

"reverts the removal of given (contingency) branches"
function reset_contingency!(solver::LinEqSysRfSolver, pm::AbstractPowerModel, branches::Vector{Int64}, old_rhss::Dict{Int64, Dict{Symbol, Float64}}; nw::Int64=nw_id_default)
    if solver.is_factorized
        for branch_id in branches
            add_branch!(solver.L, solver.U, pm, branch_id, solver.bus_to_col, solver.p, solver.q, solver.Rs, nw=nw)
            add_branch!(solver.A, pm, branch_id, solver.bus_to_col, nw=nw)
        end
    else
        for branch_id in branches
            add_branch!(solver.A, pm, branch_id, solver.bus_to_col, nw=nw)
        end
    end

    # update fallback solver as well
    reset_contingency!(solver.fallback_solver, pm, branches, old_rhss, nw=nw)
end

"updates the dual problem to fit the current master values, contingency & primal solution for a feasibility cut !after the primal was solved!"
function update_dual_feas!(solver::LinEqSysRfSolver, pm::AbstractPowerModel, master_var_vals::Dict{Symbol, AbstractVector}, outage_branches::Vector{Int64}, betas::Dict{Int64, Float64}; nw::Int64=nw_id_default)
    ## we are constructing the leq (  A    1    | 0 )
    ##                             (  A   -1    | 0 )  
    ##                             (  dᵀ vadᵐᵃˣ | 1 )
    ## where the 1 and -1 are at the positions of fbus and tbus off a branch with vad violation, and dᵀ is the objective function coefficient of the dual subproblem

    # we are solving something like Ax = 0
    fill!(solver.dual_rhs, 0.0)
    push!(solver.dual_rhs, 6.66)
    
    # add obj func constraint (we need a positive direction)
    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{Float64}()
    for bus_id in 1:solver.n_buses
        obj_coeff = solver.power_surplus[solver.bus_to_col[bus_id]]
        if obj_coeff != 0.0
            push!(I, 1)
            push!(J, solver.bus_to_col[bus_id])
            push!(V, obj_coeff)
        end
    end

    solver.feas_row = sparse(I,J,V, 1, size(solver.A)[2])
    
    # get a branch with vad violation
    for (branch_id, vad) in betas
        branch = solver.branch_data[branch_id]
        
        # get feas col and diag
        solver.feas_col = build_feas_col(solver, branch, solver.n_buses)

        solver.feas_diag = get_feas_diag(master_var_vals, outage_branches, branch_id, branch, vad)

        break
    end
end