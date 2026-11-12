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

using TOML

"given two dictionaries, extend the first dictionary with the second dictionaries entries (adhere to nested structure)"
function extend_dict!(dict1::Dict, dict2::Dict)
    for (key, value) in dict2
        if haskey(dict1, key) && value isa AbstractDict && dict1[key] isa AbstractDict
            extend_dict!(dict1[key], value)
        else
            dict1[key] = value
        end
    end
end

"convert parameter dict by replacing string keys with symbols and recursively applying the convert method to all values"
function convert_params(value::AbstractDict)
    return Dict{Symbol, Any}(Symbol(k) => convert_params(v) for (k, v) in value)
end

"convert lists of values by applying the convert method to all elements"
function convert_params(value::AbstractVector)
    return map(convert_params, value)
end

"convert values of type string by checking if they are a reference to a symbol, function, variable or simply a string"
function convert_params(value::String)

    # empty strings remain as-is
    if isempty(value)
        return value
    end

    # check if value is a reference to a symbol
    if startswith(value, ":")
        return Symbol(value[2:end])
    end

    # check if value is a reference to a filetype (handled extra, as e.g. json is a function as well as a filetype)
    if startswith(value, ".")
        return value[2:end]
    end

    # check if value is reference to numerical value infinity
    if value == "Inf" || value == "-Inf"
        return parse(Float64, value)
    end
    
    try
        # try to parse the value as a julia expression
        parsed_value = Meta.parse(value)
        evaluated_value = eval(parsed_value)
        # if the value is a function, replace the value by base invoking the function to avoid errors due to the functions age
        return isa(evaluated_value, Function) ? 
            (args...) -> Base.invokelatest(evaluated_value, args...) : 
            evaluated_value
    catch e
        # if the value is not a valid julia expression, return the value as string
        return value
    end
end

"convert all other values by returning them as is"
function convert_params(value::Any)
    return value
end

"parse the parameter file at the given file path and return the converted parameters (replace string keys with symbols and convert values)"
function parse_params(file_path::String)
    toml_data::Dict{String, Any} = TOML.parsefile(file_path)
    converted_toml_data::Dict{Symbol, Any} = convert_params(toml_data)
    return converted_toml_data
end

"backwards compatibility: if a legacy `[BENDERS] BASE_PF_IN_MASTER` value is present (e.g. from older TOMLs or `additional_params`), move it into `[BENDERS.PARTIAL_DECOMP]` (set `ENABLED=true`, `INCLUDE_BASE_PF=<legacy>`) and remove the legacy key. Call once after all TOML/extend operations and before `refresh_derived_params!()`."
function apply_legacy_base_pf_in_master!()
    if haskey(PARAMS[:BENDERS], :BASE_PF_IN_MASTER)
        legacy = pop!(PARAMS[:BENDERS], :BASE_PF_IN_MASTER)
        PARAMS[:BENDERS][:PARTIAL_DECOMP][:ENABLED]         = true
        PARAMS[:BENDERS][:PARTIAL_DECOMP][:INCLUDE_BASE_PF] = legacy
    end
    return nothing
end

"refresh derived flags in PARAMS that combine multiple TOML fields (call after any load/mutation of `[BENDERS.PARTIAL_DECOMP]`)."
function refresh_derived_params!()
    pd = PARAMS[:BENDERS][:PARTIAL_DECOMP]
    PARAMS[:BENDERS][:BASE_PF_IN_MASTER] = pd[:ENABLED] && pd[:INCLUDE_BASE_PF]
    return nothing
end