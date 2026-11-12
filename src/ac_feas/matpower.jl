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

using MATLAB

"extracts ac power flows of branches from the matpower results"
function get_branch_data_matpower(baseMVA::Float64)

    # obtain the (re)active power flows
    mat"$pfs = results.branch(:, PF);"
    mat"$qfs = results.branch(:, QF);"
    mat"$pts = results.branch(:, PT);"
    mat"$qts = results.branch(:, QT);"

    return pfs / baseMVA, qfs / baseMVA, pts / baseMVA, qts / baseMVA
end

"extracts bus voltages from the matpower results"
function get_bus_data_matpower(baseMVA::Float64)
    # obtain the voltage angles and magnitudes
    mat"$vms = results.bus(:, VM);"
    mat"$vas = results.bus(:, VA);"
    return vms, vas / baseMVA
end

"checks the ac feasibility of a given network and result dict using matpower"
function check_ac_feasibility_matpower!(result::Dict{String, Any}, data::Dict{String, Any})

    for nw_id in keys(result["solution"]["nw"])
    
        start_time = time()

        # export data to matpower temporary case file
        case_file_path = "temp_case.m"
        export_file(case_file_path, data["nw"][nw_id])

        # read network in matpower
        mat"define_constants;"
        mat"$mpc = loadcase($case_file_path);" # TODO potentially move this to outer function
        
        # remove the temporary case file
        rm(case_file_path)

        # check ac feasibility
        mat"results = runpf($mpc);"
        mat"$feasible = results.success;"

        # parse violations
        baseMVA = data["nw"][nw_id]["baseMVA"]
        pfs, qfs, pts, qts = get_branch_data_matpower(baseMVA)
        vms, vas = get_bus_data_matpower(baseMVA)
        active = collect(mpc["branch"][:, 11] .> 0.5)
        rates = mpc["branch"][:, 6] / baseMVA
        violations_from, violations_to = calculate_violations(pfs, qfs, pts, qts, rates, active)

        run_time = time() - start_time

        # write ac-feasibility results to result dict and print the according information
        modify_result!(result, nw_id, feasible, run_time, violations_from, violations_to, pfs, qfs, pts, qts, vms, vas)
        print_info_acfeas(result, nw_id)
    end

end