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

"estimates the losses of the branches of a given network using matpower"
function estimate_losses_matpower!(data::Dict{String, Any})

    overall_run_time = 0.0

    for network in values(data["nw"])

        start_time = time()

        # export data to matpower temporary case file
        case_file_path = "temp_case.m"
        export_file(case_file_path, network)

        # read network in matpower
        mat"define_constants;"
        mat"$mpc = loadcase($case_file_path);"

        # remove the temporary case file
        rm(case_file_path)

        # run optimal power flow
        mat"results = runopf($mpc);"

        # check if opf converged
        mat"$feasible = results.success;"
        if !feasible
            println("OPF did not converge for network $(network["name"]). Constant loss estimation not possible.")
            return
        end

        # get branch data
        baseMVA = network["baseMVA"]
        mat"$pfs = results.branch(:, PF);"
        mat"$pts = results.branch(:, PT);"
        pfs, pts = pfs / baseMVA, pts / baseMVA

        # calculate losses
        losses = pts - pfs

        # calculate the run time
        overall_run_time += time() - start_time

        # write losses to branches
        add_losses!(network, losses)
    end
    
    data["loss_estimation_time"] = overall_run_time
end