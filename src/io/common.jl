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

using FileIO, JLD2, JSON, DataFrames, JuMP

export check_casefile, check_cachefile, write_slack_values, read_contingenciess, write_contingenciess

"set start values for solver"
function set_dc_pf_start_values!(data::Dict{String,<:Any})
    for (n_id, network) in data["nw"]
        for (g_id, gen) in network["gen"]
            gen["pg_start"] = gen["pg"]
            if gen["pg"] > gen["pmax"]
                @error("Initial solution inconsistency: gen $(gen["name"]) pmax < pg <=> $(branch["pmax"]) < $(branch["pg"])")
            end
            gen["z_gen_start"] = gen["status"]
            #gen["qg_start"] = gen["qg"]
            #gen["z_gen_start"] = gen["status"]
        end
        for (bus_id, bus) in network["bus"]
            bus["va_start"] = bus["va"]
            #bus["va_start"] -= add[parse(Int64, bus_id)]
            #@show(bus_id, bus["va_start"])
            #bus["vm_start"] = bus["vm"]
            #bus["vr_start"] = bus["vr"]
            #bus["vi_start"] = bus["vi"]
        end
        for (b_id, branch) in network["branch"]
            if !haskey(branch, "type") || branch["type"] ∉ ["K", "breaker"]
                if haskey(branch, "pf")
                    branch["p_start"]  = branch["pf"]
                    branch["pf_start"] = branch["pf"]
                    branch["pt_start"] = branch["pt"]
                    branch["q_start"] = branch["qf"]
                    branch["qf_start"] = branch["qf"]
                    branch["qt_start"] = branch["qt"]
                end
            end
            branch["z_branch_start"] = haskey(branch, "status") ? branch["status"] : 1
        end
    end
end

"set start values. round for integers, ceil for limits"
function set_dc_pf_start_values!(data::Dict{String, <:Any}, solution; fix_uc=false, fix_ts=false, fix_pf=false, fix_tr=false, fix_benders=false, fix_slack=false, fix_pg=false, benders=true)
    @assert(!fix_pg || fix_uc, "power generation must be fixed if uc is fixed")
    data["WARMSTART"]   = true
    data["FIX_TS"]      = fix_ts
    data["FIX_UC"]      = fix_uc
    data["FIX_PF"]      = fix_pf
    data["FIX_TR"]      = fix_tr
    data["FIX_BENDERS"] = fix_benders
    data["FIX_SLACK"]   = fix_slack
    data["FIX_PG"]      = fix_pg
    for (n_id, network) in data["nw"]
        network["WARMSTART"]   = true
        network["FIX_TS"]      = fix_ts
        network["FIX_UC"]      = fix_uc
        network["FIX_PF"]      = fix_pf
        network["FIX_TR"]      = fix_tr
        network["FIX_BENDERS"] = fix_benders
        network["FIX_SLACK"]   = fix_slack
        network["FIX_PG"]      = fix_pg
        _solution = solution["nw"][n_id]
        if haskey(_solution, "gen") 
            for (g_id, gen) in network["gen"]
            # there doesn't have to be a solution for the gen, if we solve BENDERS_NESTED
                if haskey(_solution["gen"], g_id)
                    if haskey(_solution["gen"][g_id], "pg")
                        if fix_pg 
                            gen["pg_start"]                 = _solution["gen"][g_id]["pg"]
                            pg_delta                        = _solution["gen"][g_id]["pg"] - gen["pg_0"]
                            gen["RDP+"]                     += pg_delta
                            gen["RDP-"]                     -= pg_delta
                            if pg_delta < 0
                                gen["pg_delta_plus_start"]      = pg_delta
                                gen["pg_delta_minus_start"]     = 0
                            else
                                gen["pg_delta_plus_start"]      = 0
                                gen["pg_delta_minus_start"]     = pg_delta
                            end

                        else
                            gen["pg_start"]                 = _solution["gen"][g_id]["pg"]
                            gen["qg_start"]                 = _solution["gen"][g_id]["qg"]
                        end
                    end

                    if haskey(_solution["gen"][g_id], "pg_delta_plus")
                        gen["pg_delta_plus_start"]      = _solution["gen"][g_id]["pg_delta_plus"]
                        gen["pg_delta_minus_start"]     = _solution["gen"][g_id]["pg_delta_minus"]
                        gen["pg_delta_on_off_start"]    = round(get(_solution["gen"][g_id], "pg_delta_on_off", 1))
                    end

                    if haskey(_solution["gen"][g_id], "gen_startup")
                        gen["z_gen_startup_start"]      = round(get(_solution["gen"][g_id], "gen_startup", 1))
                        gen["z_gen_shutdown_start"]     = round(get(_solution["gen"][g_id], "gen_shutdown", 1))
                        gen["z_gen_start"]              = round(get(_solution["gen"][g_id], "gen_status", 1))
                    end                    
                end
            end
        end

        # there doesn't have to be a solution for the bus, if PARAMS[:BENDERS][:BASE_PF_IN_MASTER] = false
        if haskey(_solution, "bus")
            for (bus_id, bus) in network["bus"]
                if ! haskey(_solution["bus"], bus_id)
                    continue
                end
                if haskey(_solution["bus"][bus_id], "va")
                    bus["va_start"] = _solution["bus"][bus_id]["va"]
                end
            end
        end

        n_added_ratea_slack = 0
        n_added_z_branch = 0
        n_added_tr_shift = 0
        n_added_eta = 0

        if haskey(_solution, "branch")
            for (b_id, branch) in network["branch"]
                if ! haskey(_solution["branch"], b_id)
                    continue
                end
                # there doesn't have to be a solution for the branch, if PARAMS[:BENDERS][:BASE_PF_IN_MASTER] = false
                if haskey(_solution["branch"][b_id], "pf")
                    branch["p_start"]               = _solution["branch"][b_id]["pf"]
                    branch["pf_start"]              = _solution["branch"][b_id]["pf"]
                    branch["pt_start"]              = _solution["branch"][b_id]["pt"]
                end

                if haskey(_solution["branch"][b_id], "rat_s")
                    n_added_ratea_slack += 1
                    branch["rate_a_slack_on_off_start"] = round(get(_solution["branch"][b_id], "rat_s_b", 0))
                    branch["rat_s_start"] = _solution["branch"][b_id]["rat_s"]
                    branch["raf_s_start"] = _solution["branch"][b_id]["raf_s"]
                end

                if haskey(_solution["branch"][b_id], "br_status")
                    n_added_z_branch += 1
                    branch["z_branch_start"]        = round(_solution["branch"][b_id]["br_status"])
                end
                   
                if haskey(_solution["branch"][b_id], "br_change")
                    branch["z_branch_change_start"] = round(_solution["branch"][b_id]["br_change"])
                    if branch["z_branch_change_start"] == 1
                        @info("Switched branch: $b_id $(get(branch, "branch_name_intern", ""))")
                    end
                end

                if haskey(_solution["branch"][b_id], "tr_shift")
                    n_added_tr_shift += 1
                    branch["tr_shift_start"]        = _solution["branch"][b_id]["tr_shift"]
                    branch["tr_shift_change_start"] = get(_solution["branch"][b_id], "tr_shift_change", branch["cur_shift_step"])
                end
            end
        end
        # assign η variable values as well
        if haskey(_solution, "eta") && benders
            n_added_eta += 1
            eta_dict = Dict{String, Dict{String, Any}}(Dict("$i" => Dict("contingency" => contingency, "η_start" => ceil(get(_solution["eta"], "$i", Dict("η" => 0.0))["η"], digits=PARAMS[:PRECISION][:SOLVER_PREC_DIGITS]-2)) for (i, contingency) in enumerate(network["contingencies"])))
            network["η"] = eta_dict
            @assert(length(_solution["eta"]) == length(network["contingencies"]))
        end
        @debug("Set start values for network $n_id: $(length(network["gen"])) gens, $(length(network["bus"])) buses, $(length(network["branch"])) branches, $n_added_ratea_slack rate_a_slack_on_off, $n_added_z_branch z_branch, $n_added_tr_shift tr_shift, $n_added_eta η.")
    end
end

"check if casefile (including dir) exists"
function check_casefile(casefile)
    # check if case file
    if ! isfile(casefile)
        Memento.error(_LOGGER, "Case file $casefile not found.")
    end
end

"helper function to extract slack variable values from result DataFrame"
function get_slacks_from_result(result, slack_types::Vector{Symbol})
    # initialize dataframe to save slack variables to 
    data = DataFrame(component_type=Symbol[], component_id=String[], slack_variable_name=String[], slack_variable_value=Float64[])

    # obtain data from result vector
    for slack_type in slack_types
        slack_component_type = CONSTS[:INTERNAL_MAPPINGS][:SLACK_VARIABLES][:SLACK_COMPONENT_TYPE][slack_type]
        for component in get(result["solution"], string(slack_component_type), [])
            for slack_variable_name in CONSTS[:INTERNAL_MAPPINGS][:SLACK_VARIABLES][:RESULT_ARRAY_SLACK_VARIABLE_NAME_MAP][slack_type]
                slack_variable_name = string(slack_variable_name)
                if slack_variable_name in keys(component[2])
                    slack_val = component[2][slack_variable_name]
                    slack_val = ceil(slack_val, digits=SOLVER_PRECISION_DIGITS)
                    @assert(slack_val >= 0, "Slack value $slack_variable_name = $slack_val for component $slack_component_type $component is negative.")
                    if slack_val > SOLVER_PRECISION
                        push!(data, (slack_component_type, component[1], slack_variable_name, slack_val))
                    end
                end
            end
        end
    end
    return data
end

"write slack variable values"
function write_slack_values(slack_types::Vector{Symbol}, result, casename::String; slack_cache=CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:SLACK_CACHE_DIR], filetype=slack_cache=CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:SLACK_CACHE_FILETYPE])
    slacks = get_slacks_from_result(result, slack_types)
    # write file
    mkpath(slack_cache)
    CSV.write(joinpath(slack_cache,"slack_$(casename).$(filetype)"), slacks)
end

"adds buses to all existing buses according to their `n_busbars` to turn them into substations"
function add_busbars_buses(bus_data::Dict)::Nothing
    max_bus_id = maximum([parse(Int64, bus_id) for (bus_id, bus) in bus_data])
    current_bus_id = max_bus_id + 1
    
    # add new buses to each bus according to their `n_busbars` 
    for bus_id in deepcopy(eachindex(bus_data)) # deepcopy is used to ensure that we don't visit the same element twice due to adding stuff to the Dict
        bus = bus_data[bus_id]

        # add a name if necessary
        if !haskey(bus, "name")
            bus["name"] = bus_id
        end
        
        # copy the bus, then replace its `name`, `index` and `bus_i`
        for i in 2:bus["n_busbars"]
            new_bus = deepcopy(bus)
            new_bus["name"] = new_bus["name"] * "_$i"
            new_bus["bus_i"] = current_bus_id
            new_bus["index"] = current_bus_id

            bus_data["$current_bus_id"] = new_bus

            current_bus_id += 1
            
        end

        # the existing bus is the "first" busbar of a substation
        bus["name"] = bus["name"] * "_1"
    end
end

"add missing attributes with values following the provided functions to components in data"
function add_missing_attributes(data, missing_list::Dict{Symbol})
    for (component_type, attribute_list) in missing_list
        for (_, component) in data[string(component_type)]
            unresolved_attributes::Vector{Tuple{Symbol, Function}} = []
            for (attribute, value_func) in attribute_list
                if !haskey(component, "$attribute")
                    try
                        component[string(attribute)] = value_func(component)
                    catch
                        push!(unresolved_attributes, (attribute, value_func))
                    end
                end
            end
            for (attribute, value_func) in unresolved_attributes
                component[string(attribute)] = value_func(component)
            end
        end
    end
end

"set bus index to be consecutive"
function make_bus_indices_consecutive!(data::Dict, data_branch::Dict{String, Any}, data_bus::Dict{String, Any}, data_gen::Dict{String, Any}, data_load::Dict{String, Any})
    # create a mapping from the old bus indices to the new, consecutive ones
    bus_ids = sort([bus["index"] for (bus_id, bus) in data_bus])
    if bus_ids[end] == length(bus_ids) # if the indices are already consecutive
        return
    else
        bus_id_map = Dict{Int64, Int64}(bus_id => i for (i, bus_id) in enumerate(bus_ids))
    end

    # replace data_bus
    for (_, bus) in data_bus
        new_index = bus_id_map[bus["index"]]
        bus["index"] = new_index
        bus["bus_i"] = new_index
        bus["source_id"] = ["bus", new_index]
    end
    data["bus"] = Dict{String, Any}(string(bus["index"]) => bus for (b_id, bus) in data_bus)

    for (_, branch) in data_branch
        branch["f_bus"] = bus_id_map[branch["f_bus"]]
        branch["t_bus"] = bus_id_map[branch["t_bus"]]
    end
    for (_, load) in data_load
        load["load_bus"] = bus_id_map[load["load_bus"]]
        load["source_id"] = ["bus", load["load_bus"]]
    end
    for (_, gen) in data_gen
        gen["gen_bus"] = bus_id_map[gen["gen_bus"]]
        gen["source_id"] = ["bus", gen["gen_bus"]]
    end
    for (_, shunt) in data["shunt"]
        shunt["shunt_bus"] = bus_id_map[shunt["shunt_bus"]]
        shunt["source_id"] = ["bus", shunt["shunt_bus"]]
    end

    # TODO: other components (???,...)

    # update substation data structure
    if haskey(data, "substations")
        for substation in values(data["substations"])
            map!(x -> string(bus_id_map[parse(Int64, x)]), substation["busbar"], substation["busbar"])

            for aux_bus in substation["aux_bus"]
                aux_bus["bus_i"] = bus_id_map[aux_bus["bus_i"]]
            end
        end
    end

    @assert(maximum([parse(Int64, b_id) for (b_id, bus) in data["bus"]]) == length(data["bus"]), "Bus indices are not consecutive.")
end


"set branch index to be consecutive"
function make_branch_indices_consecutive!(data::Dict{String, Any}, data_branch::Dict{String, Any}, data_bus::Dict{String, Any}, data_substations::Dict{String, Dict{String, Vector{Any}}}, data_contingencies::Vector{Dict{Symbol, Vector{Int64}}})
    # create a mapping from the old branch indices to the new, consecutive ones
    branch_id_map = Dict{Int64, Int64}(branch["index"] => i for (i, (branch_id, branch)) in enumerate(data_branch))

    # replace data_branch
    for (_, branch) in data_branch
        new_index = branch_id_map[branch["index"]]
        if haskey(branch, "breaker_branches")
           map!(branch_id -> branch_id_map[branch_id], branch["breaker_branches"], branch["breaker_branches"])
        end
        branch["index"] = new_index
        branch["source_id"] = ["branch", new_index]
    end
    data["branch"] = Dict{String, Any}(string(branch["index"]) => branch for (_, branch) in data_branch)

    # TODO: other components (shunts, ...)?

    # update substation data structure ("K" and "aux_bus" -> "breaker_branches")
    for substation in values(data_substations)
        # map coupler indices
        map!(branch_id -> string(branch_id_map[parse(Int64, branch_id)]), substation["K"], substation["K"])

        for aux_bus in substation["aux_bus"]
            # update in substation
            map!(branch_id -> branch_id_map[branch_id], aux_bus["breaker_branches"], aux_bus["breaker_branches"])
        end
    end

    # update contingency list
    for contingency in data_contingencies
        map!(branch_id -> branch_id_map[branch_id], contingency[:branch], contingency[:branch])
    end

    @assert(maximum([parse(Int64, b_id) for (b_id, branch) in data["branch"]]) == length(data["branch"]), "Branch indices are not consecutive. Length of data_branch: $(length(data["branch"])). Maximum branch index: $(maximum([parse(Int64, b_id) for (b_id, branch) in data["branch"]])). Length of branch_id_map: $(length(branch_id_map)).")
end

"round component attributes to a given number of digits to avoid numerical problems"
function round_component_data!(data; digits=PARAMS[:PRECISION][:SOLVER_PREC_DIGITS])
    if PARAMS[:PRECISION][:SOLVER_PREC_DIGITS] < 0 || ! PARAMS[:ROUNDING][:INPUT]; return; end

    for (i, branch) in data["branch"]
        branch["angmin"] = floor(branch["angmin"], digits=digits)
        branch["angmax"] = ceil(branch["angmax"], digits=digits)
        branch["rate_a"] = ceil(branch["rate_a"], digits=digits)
    end
    
    for (i, load) in data["load"]
        load["pd"] = round(load["pd"], digits=digits)
    end
    #for (i, bus) in data["bus"]
        #bus["pd"] = round(bus["pd"], digits=digits)
    #end
end
