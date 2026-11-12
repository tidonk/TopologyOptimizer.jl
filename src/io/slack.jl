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

export write_matpower_with_slacks

"apply found slacks and export as matpower file"
function write_matpower_with_slacks(data, result, casename; slacked_cases_dir=CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:SLACKED_CASES_DIR], slack_cases_filetype=CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:SLACKED_CASE_FILETYPE], slack_map=CONSTS[:INTERNAL_MAPPINGS][:SLACK_VARIABLES][:SLACK_VARIABLE_CALC_MAP])
    
    data        = deepcopy(data)
    result      = deepcopy(result)
    slack_types = data["slack_types"]

    # if multinetwork, use first network and first network's solution
    if "nw" in collect(keys(data))
        # only use the first network
        for key in keys(data["nw"]["1"])
            data[key] = data["nw"]["1"][key]
        end
        delete!(data, "nw")
        delete!(data, "multinetwork")

        # only use result of first network solving
        for key in keys(result["solution"]["nw"]["1"])
            result["solution"][key] = result["solution"]["nw"]["1"][key]
        end
        delete!(result["solution"], "nw")
    end

    # get slack variable values from provided solution
    slacks = get_slacks_from_result(result, slack_types)

    # iterate over component types, ids and slack variables
    for c_type_data in groupby(slacks, :component_type)
        c_type = string(c_type_data.component_type[1])
        for c_id_data in groupby(c_type_data, :component_id)
            c_id = string(c_id_data.component_id[1])
            for c_var_data in groupby(c_id_data, :slack_variable_name)
                sl_var = Symbol(c_var_data.slack_variable_name[1])

                # get old param and value to add or subtract (whatever slack_map[sl_var][1] says)
                old_param = data[c_type][c_id][slack_map[sl_var][2]]
                slack_val = c_var_data.slack_variable_value[1]
                
                # calculate and set new param value
                new_param = slack_map[sl_var][1](old_param, slack_val)
                data[c_type][c_id][slack_map[sl_var][2]] = new_param
            end
        end
    end

    # export modified case data to case file
    Memento.info(_LOGGER, "Writing case file including modifications according to slack optimal values to '$(joinpath(slacked_cases_dir, "$casename.$slack_cases_filetype"))'.\n")
    mkpath(slacked_cases_dir)
    export_matpower(joinpath(slacked_cases_dir, "$casename.$slack_cases_filetype"), data)
end