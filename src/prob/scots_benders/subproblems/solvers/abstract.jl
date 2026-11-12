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

## abstract
"returns calculated value of vad and pf at the branch"
function get_vad_pf(solver::AbstractSubproblemSolver, branch::Dict; nw::Int64=nw_id_default)::Tuple{Float64, Float64}
end

"solves problem in its current state, returns success flag"
function solve_primal!(solver::AbstractSubproblemSolver)::Bool
    return false
end

"solves dual problem in its current state !after the primal was solved! to obtain a optimality cut, returns success flag"
function solve_dual_opt!(solver::AbstractSubproblemSolver)::Bool
    return false
end

"solves dual problem in its current state !after the primal was solved! to obtain a feasibility cut, returns success flag"
function solve_dual_feas!(solver::AbstractSubproblemSolver)::Bool
    return false
end

"applies the current master solution to solver"
function apply_master!(solver::AbstractSubproblemSolver, pm::AbstractPowerModel, master_var_vals::Dict{Symbol, AbstractVector}; nw::Int64=nw_id_default)
end

"resets the solver back to its base state"
function reset_master!(solver::AbstractSubproblemSolver, pm::AbstractPowerModel, master_var_vals::Dict{Symbol, AbstractVector}; nw::Int64=nw_id_default)
end

"applies one contingency to the solver"
function apply_contingency!(solver::AbstractSubproblemSolver, pm::AbstractPowerModel, branches::Vector{Int64}, master_var_vals::Dict{Symbol, AbstractVector}, is_int::Bool; nw::Int64=nw_id_default)::Dict{Int64, Dict{Symbol, Float64}}
end

"resets the solver to the current master solution"
function reset_contingency!(solver::AbstractSubproblemSolver, pm::AbstractPowerModel, branches::Vector{Int64}, old_rhss::Dict{Int64, Dict{Symbol, Float64}}; nw::Int64=nw_id_default)
end

"updates the dual problem to fit the current master values, contingency & primal solution for an optimality cut !after the primal was solved!"
function update_dual_opt!(solver::AbstractSubproblemSolver, pm::AbstractPowerModel, master_var_vals::Dict{Symbol, AbstractVector}, betas::Dict{Int64, Float64}; nw::Int64=nw_id_default)
end

"updates the dual problem to fit the current master values, contingency & primal solution for a feasibility cut !after the primal was solved!"
function update_dual_feas!(solver::AbstractSubproblemSolver, pm::AbstractPowerModel, master_var_vals::Dict{Symbol, AbstractVector}, outage_branches::Vector{Int64}, betas::Dict{Int64, Float64}; nw::Int64=nw_id_default)
end

"check dual solution of a solver for inconsistencies"
function check_dual_solution(solver::AbstractSubproblemSolver)
    return true
end