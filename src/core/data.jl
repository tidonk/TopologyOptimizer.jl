# This file contains modified code derived from PowerModels.jl.
# Copyright (c) 2016, Los Alamos National Security, LLC. All rights reserved.
# These portions retain the license and disclaimer in LICENSES/PowerModels.txt.
# PowerModels.jl code in this file has been modified for TopologyOptimizer.jl.

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

"applies certain parameters directly to `data`: sets shift/switchability "
function apply_params!(data::Dict{String, <:Any}, busbars::Bool)
    if haskey(data, "nw")
        for (_, data_nw) in data["nw"]
            apply_params!(data_nw, busbars)
        end
        return
        
    end
    if !PARAMS[:LIMITS][:TRANSFORMER][:USE_TRAFOS]
        for (_, branch) in data["branch"]
            if get(branch, "transformer", false)
                branch["shift"] = 0.0
                branch["tap"] = 1.0
            end
        end
    end
    if PARAMS[:DEBUG][:REMOVE_SHUNTS]
        # remove shunts from data
        data["shunt"] = Dict{String, Any}()
    end
    if PARAMS[:DEBUG][:REMOVE_DC_LINES]
        # remove dc lines from data
        data["dcline"] = Dict{String, Any}()
    end
    if !PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_SHIFTING]
        for (_, branch) in data["branch"]
            branch["shiftable"] = false
        end
    else
        for (_, branch) in data["branch"]
            if !haskey(branch, "shiftable")
                branch["shiftable"] = false
            end
        end
    end
    if !PARAMS[:LIMITS][:GENERATOR][:GEN_SWITCHING]
        for (_, gen) in data["gen"]
            gen["switchable"] = false
        end
    end
    if !PARAMS[:LIMITS][:BRANCH][:BRANCH_SWITCHING]
        for (_, branch) in data["branch"]
            branch["switchable"] = false
        end
    end
    if !busbars
        for (_, branch) in data["branch"]
            branch["electrical"] = true
        end
    end
end

component_table(data::Dict{String,<:Any}, component::String, args...) = _IM.component_table(data, component, args...)

"recursively applies new_data to data, not overwriting information, for a single timestep"
function update_data_safe!(data, result, timestep)
    # add all entries from `result["solution"]` to data, and append "_new" to each
    for component in ["branch", "gen", "bus", "load"]
        for (c_id, c_data) in get(result["solution"], component, [])
            if ! occursin("(", c_id)
                #existing_keys = keys(data[component][first(c)])
                for key in keys(result["solution"][component][c_id])
                    new_key = "$(key)_new" # * (key in existing_keys ? "_new" : "")
                    data[component][c_id][new_key] = c_data[key]
                end
            end
        end
    end
end

"wrapper function to recursively applies new_data to data, not overwriting information, for one or multiple time steps"
function update_data_safe!(data, result)
    # if multinetwork, call for each network
    if haskey(data, "nw")
        for timestep ∈ 1:length(data["nw"])
            update_data_safe!(data["nw"][string(timestep)], Dict("solution" => result["solution"]["nw"][string(timestep)]), timestep)
        end
    else
        update_data_safe!(data, result, 0)
    end
end

"Deletes all components specified in inactive_components from data"
function remove_inactive_components!(data::Dict{String,<:Any}, inactive_components::Dict{Symbol,Vector{Int64}})
    for (component_type, component_ids) in inactive_components
        delete!.([data[String(component_type)] for id in component_ids], ["$c" for c in component_ids])
    end
end

"Deletes all components specified in inactive_components from data"
function remove_inactive_components!(data::Dict{String,<:Any}, inactive_component_type::String, inactive_component_id::String)
    remove_inactive_components!(data, Dict(inactive_component_type => [inactive_component_id]))
end

"Transforms network data into per-unit"
function make_per_unit!(data::Dict{String,<:Any})
    pm_data = get_pm_data(data)

    if !haskey(pm_data, "per_unit") || pm_data["per_unit"] == false
        pm_data["per_unit"] = true
        apply_pm!(PowerModels._make_per_unit!, pm_data)
    end

    if haskey(pm_data, "gen")
        mva_base = pm_data["baseMVA"]

        # to be consistent with matpower's opf.flow_lim= 'I' with current magnitude
        # limit defined in MVA at 1 p.u. voltage
        ka_base = mva_base

        rescale        = x -> x/mva_base
        
        _apply_func! = PowerModels._apply_func!
    end
end