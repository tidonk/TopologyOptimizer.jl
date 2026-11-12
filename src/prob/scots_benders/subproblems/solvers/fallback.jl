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

# buffers update info so that apply_* is called only when needed
## abstract
"returns calculated value of vad and pf at the branch"
function get_vad_pf(solver::FallbackSolver, branch::Dict; nw::Int64=nw_id_default)::Tuple{Float64, Float64}
    return get_vad_pf(solver.solver, branch, nw=nw)
end

"solves problem in its current state, returns success flag"
function solve_primal!(solver::FallbackSolver)::Bool
    return solve_primal!(solver.solver)
end

"solves dual problem in its current state !after the primal was solved! to obtain a optimality cut, returns success flag"
function solve_dual_opt!(solver::FallbackSolver)::Bool
    return solver_dual!(solver.solver)
end

"applies the current master solution to solver"
function apply_master!(solver::FallbackSolver, pm::AbstractPowerModel, master_var_vals::Dict{Symbol, AbstractVector}; nw::Int64=nw_id_default)
    return apply_master!(solver.solver, pm, master_var_vals, nw=nw)
end

"resets the solver back to its base state"
function reset_master!(solver::FallbackSolver, pm::AbstractPowerModel, master_var_vals::Dict{Symbol, AbstractVector}; nw::Int64=nw_id_default)
    return reset_master!(solver.solver, pm, master_var_vals, nw=nw)
end

"applies one contingency to the solver"
function apply_contingency!(solver::FallbackSolver, pm::AbstractPowerModel, branches::Vector{Int64}, master_var_vals::Dict{Symbol, AbstractVector}, is_int::Bool; nw::Int64=nw_id_default)::Dict{Int64, Dict{Symbol, Float64}}
    solver.contingency_branches, solver.master_var_vals, solver.is_int = branches, master_var_vals, is_int
    return Dict{Int64, Dict{Symbol, Float64}}()
end

"resets the solver to the current master solution"
function reset_contingency!(solver::FallbackSolver, pm::AbstractPowerModel, branches::Vector{Int64}, old_rhss::Dict{Int64, Dict{Symbol, Float64}}; nw::Int64=nw_id_default)
end

"updates the dual problem to fit the current master values, contingency & primal solution for an optimality cut !after the primal was solved!"
function update_dual_opt!(solver::FallbackSolver, pm::AbstractPowerModel, master_var_vals::Dict{Symbol, AbstractVector}, betas::Dict{Int64, Float64}; nw::Int64=nw_id_default)
    return update_dual_opt!(solver.solver, pm, master_var_vals, betas, nw=nw)
end

"updates the dual problem to fit the current master values, contingency & primal solution for a feasibility cut !after the primal was solved!"
function update_dual_feas!(solver::FallbackSolver, pm::AbstractPowerModel, master_var_vals::Dict{Symbol, AbstractVector}, outage_branches::Vector{Int64}, betas::Dict{Int64, Float64}; nw::Int64=nw_id_default)
    return update_dual_feas!(solver.solver, pm, master_var_vals, outage_branches, betas, nw=nw)
end

"applies the buffered update of one contingency to the fallback solver (solver.solver) by calling `apply_contingency!` with the buffered data"
function apply_fallback_update!(solver::FallbackSolver)::Dict{Int64, Dict{Symbol, Float64}}
    return apply_contingency!(solver.solver, solver.pm_sub, solver.contingency_branches, solver.master_var_vals, solver.is_int, nw=solver.nw)
end

"resets the buffered update of one contingency in the fallback solver (solver.solver) by calling `reset_contingency!` with the buffered data"
function reset_fallback_update!(solver::FallbackSolver, old_rhss::Dict{Int64, Dict{Symbol, Float64}})
    return reset_contingency!(solver.solver, solver.pm_sub, solver.contingency_branches, old_rhss, nw=solver.nw)
end