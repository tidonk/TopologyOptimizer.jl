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

## lp
"returns calculated value of vad and pf at the branch"
function get_vad_pf(solver::LPSolver, branch::Dict; nw::Int64=nw_id_default)::Tuple{Float64, Float64}
    vad = JuMP.value(var(solver.pm, nw, :va, branch["f_bus"])) - JuMP.value(var(solver.pm, nw, :va, branch["t_bus"]))
    pf  = JuMP.value(var(solver.pm, nw, :p, (branch["index"], branch["f_bus"], branch["t_bus"])))
    return (vad, pf)
end

"intializes the solvers"
function init_lp_solvers(pm, n_threads::Int64; nw::Int64=nw_id_default)::Vector{LPSolver}
    # TODO: is it really required to copy the full `pm` just for the lp?
    pm_lps = create_lp!.([deepcopy(pm) for _ in 1:n_threads], nw=nw)

    branch_data = collect_branch_data(pm_lps[1], nw=nw)

    return [LPSolver(pm=pm_lp, nw=nw, branch_data = branch_data) for pm_lp in pm_lps]
end

"solves problem in its current state, returns success flag"
function solve_primal!(solver::LPSolver)::Bool
    JuMP.optimize!(solver.pm.model)
    return JuMP.termination_status(solver.pm.model) == JuMP.MOI.OPTIMAL
end

"applies one contingency to the solver"
function apply_contingency!(solver::LPSolver, pm::AbstractPowerModel, branches::Vector{Int64}, master_var_vals::Dict{Symbol, AbstractVector}, int_sol::Bool; nw::Int64=nw_id_default)::Dict{Int64, Dict{Symbol, Float64}}
    fix_master_variables!(solver.pm, nw, master_var_vals, int_sol, ignore_var_names = [:z_gen]) # this also applies the master
    return Dict{Symbol, Float64}()
end
