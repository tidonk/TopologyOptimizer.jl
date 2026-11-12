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

# shared functions for LinEqSysSolver and LinEqSysBBSolver
"solves problem in its current state, returns success flag"
function solve_primal!(solver::AbstractLinEqSysSolver)::Bool
    if !_is_matrix_connected(solver.A)
        return false
    end
    try
        solver.vas = solver.A\solver.power_surplus
    catch e
        if isa(e, LinearAlgebra.SingularException)
            return false
        else
            rethrow(e)
        end
    end
    return true
end

"solves dual problem in its current state !after the primal was solved! to obtain a optimality cut, returns success flag"
function solve_dual_opt!(solver::AbstractLinEqSysSolver)::Bool
    try
        solver.duals = solver.A\solver.dual_rhs
    catch e
        if isa(e, LinearAlgebra.SingularException)
            callback_stats[:n_leq_fails] += 1
            return false
        else
            rethrow(e)
        end
    end
    return true
end

"solves dual problem in its current state !after the primal was solved! to obtain a feasibility cut, returns success flag"
function solve_dual_feas!(solver::AbstractLinEqSysSolver)::Bool
    if !_is_matrix_connected(solver.A_feas_cut)
        callback_stats[:n_leq_fails] += 1
        solver.dual_rhs=zeros(solver.n_buses)
        return false
    end
    try
        solver.duals = solver.A_feas_cut\solver.dual_rhs
    catch e
        if isa(e, LinearAlgebra.SingularException)
            callback_stats[:n_leq_fails] += 1
            solver.dual_rhs=zeros(solver.n_buses)
            return false
        else
            rethrow(e)
        end
    end
    solver.dual_rhs=zeros(solver.n_buses)
    return true
end