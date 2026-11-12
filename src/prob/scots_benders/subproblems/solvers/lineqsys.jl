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

## linear equation system
"intializes the solvers"
function init_lineqsys_solvers(pm, pm_subs, n_threads::Int64, sub_cons_for_dualss; nw::Int64=nw_id_default)::Vector{LinEqSysSolver}
    A, power_surplus, bus_to_col = build_base_pf_matrix(pm, nw=nw)

    # add FullLP solver as fallback solver (one per thread)
    fallback_solver_solvers = init_fulllp_solvers(pm_subs, n_threads, sub_cons_for_dualss)
    fallback_solvers = [FallbackSolver(solver=fallback_solver_solvers[t], pm_sub=pm_subs[t], nw=nw) for t in 1:n_threads]

    branch_data = collect_branch_data(pm_subs[1], nw=nw)

    return [[LinEqSysSolver(A=A, power_surplus=power_surplus, bus_to_col=bus_to_col, fallback_solver=fallback_solvers[1], branch_data=branch_data)]; [LinEqSysSolver(A=deepcopy(A), power_surplus=deepcopy(power_surplus), bus_to_col=bus_to_col, fallback_solver=fallback_solvers[t], branch_data = branch_data) for t in 2:n_threads]]
end