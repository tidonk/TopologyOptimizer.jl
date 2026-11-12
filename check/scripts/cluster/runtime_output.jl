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

using TopologyOptimizer

function extract_keys(result_dict::Dict)
    keys_vec::Vector{String} = []
    for field in ["time", "info", "callback", "ac_feas"]
        if haskey(result_dict, field) && isa(result_dict[field], Dict)
            if field == "time"
                prefix = "time_"
            elseif field == "ac_feas"
                prefix = "ac_feas_"
            else
                prefix = ""
            end
            for subfield in keys(result_dict[field])
                push!(keys_vec, prefix * string(subfield))
            end
        end
    end
    return keys_vec
end

# write results
function write_runtime_data_to_csv(result_dict, casename::String, settingsname::String, tal::Int; results_dir = CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:LOG_DIR])
    # write results
    header = "instance_file,settings,tal," * join(extract_keys(result_dict), ",") * "\n"

    # create a row with the results
    results = []
    for field in ["time", "info", "callback"]
        if haskey(result_dict, field) && isa(result_dict[field], Dict)
            for subfield in keys(result_dict[field])
                if typeof(result_dict[field][subfield]) == Bool
                    push!(results, Int(result_dict[field][subfield]))
                else
                    push!(results, result_dict[field][subfield])
                end
            end
        end
    end

    row = vcat([casename, settingsname, tal], results)

    @info("Writing runtime data for $(casename) to results file $(joinpath(results_dir, "$(casename).$(settingsname).csv"))")
    mkpath(results_dir)
    open(joinpath(results_dir, "$(casename).$(settingsname).csv"), "w") do f
        write(f, header * join(row, ",") * "\n")
    end
end

