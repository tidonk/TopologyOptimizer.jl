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

export _rewire_graph_with_substations!, _detect_substations
"get bus by its name. returns the bus id, or if not found,-1."
function get_bus_id_by_name(data, bus_name; include_busbar=true)::String
    chars = include_busbar ? 8 : 7
    for (bus_id, bus) in data["bus"]
        if bus["name"][1:chars] == bus_name
            return bus_id
        end
    end
    return "-1"
end

function get_substation_name_filter(data)
    return CONSTS[:USER_ENVIRONMENT][:BUSBAR_RECONFIGURATION][:SUBSTATION_NAME_FILTERS][:pglib]
end

"groups individual buses into substations based on their names"
function _detect_substations(data_bus::Dict{String, Any}, substation_name_filter)::Dict{String, Dict{String,Vector{Any}}}
    # id => {field_name => [value, ...], ...}
    substations = Dict{String, Dict{String,Vector{Any}}}()
    for (b_id, bus) in data_bus
        substation_name = substation_name_filter(bus["name"]::String)::SubString{String}
        if substation_name ∉ keys(substations)
            substations[substation_name] = Dict("busbar" => [b_id], "aux_bus" => Vector{Dict{String, Any}}(), "K" => Vector{String}())
        else
            push!(substations[substation_name]["busbar"], b_id)
        end

        bus["substation"] = substation_name
    end

    # count of number of busbars in substation: 1,2,3,4,>4
    counts = [0,0,0,0,0]
    (substation_name -> counts[min(length(substations[substation_name]["busbar"]),5)] += 1).(keys(substations)) # for each substation, increase counter corresponding to its length

    @info("Detected $(length(keys(substations))) substations.
    #lengths 1:$(counts[1]), 2:$(counts[2]), 3:$(counts[3]), 4:$(counts[4]), >4:$(counts[5])")
    return substations
end


"adds buses and branches representing couplers and breakers of substations, and connects them with the remaining network components"
function _rewire_graph_with_substations!(data)
    substation_name_filter = CONSTS[:USER_ENVIRONMENT][:BUSBAR_RECONFIGURATION][:SUBSTATION_NAME_FILTERS][:pglib]
    current_bus_index = maximum(parse.(Int64, [bus_id for (bus_id, bus) in data["bus"]])) + 1

    "builds a new auxiliary branch object"
    function _build_aux_branch(f_bus_i             ::String, 
                               t_bus_i             ::String,
                               type                ::String,
                               t_substation_name   ::String, 
                               aux_index           ::Int64, 
                               current_branch_index::Int64; 
                               br_status           ::Int64  =1, 
                               status              ::Int64  =1, 
                               switchable          ::Bool   =true, 
                               rate_a              ::Float64=10000.0, 
                               in_contingencies             =false, 
                               monitored                    =false, 
                               orig_branch_id      ::Int    =-1, 
                               branch_name_intern           =""      )::Dict{String, Any}

        type_name = type == "breaker" ? "brk" : "line"
        second_name_entry = type == "breaker" ? t_substation_name : data["bus"][t_bus_i]["name"]
        @assert(type ∈ ["breaker", "K"])
        return Dict{String, Any}(
            "name" => "$(data["bus"][f_bus_i]["name"])_$(second_name_entry)_$(aux_index):$(type_name)",
            "f_bus" => parse(Int64, f_bus_i),
            "t_bus" => parse(Int64, t_bus_i),
            "switchable" => switchable,
            "br_status" => br_status,
            "status" => status,
            "status_previous" => status,
            "z_branch_start" => status,
            "index" => current_branch_index,
            "in_contingencies" => in_contingencies,
            "monitored" => monitored,
            "orig_index" => orig_branch_id,
            "type" => type,
            "angmin" => 0,
            "angmax" => 0,
            "electrical" => false,
            "tap" => 0,
            "rate_a" => rate_a,
            "transformer" => false # TODO: sometimes true?
        )
    end

    "returns true iff the bus is switchable (ie should be connected to switchable branches)"
    function is_bus_switchable(bus)::Bool
        return true
    end

    # non-switchable branches that are off can be deleted
    for (b_id, branch) in copy(data["branch"])
        t_bus_id = string(branch["t_bus"])
        f_bus_id = string(branch["f_bus"])

        if (!is_bus_switchable(data["bus"][t_bus_id]) && !is_bus_switchable(data["bus"][f_bus_id])) && branch["status"] == 0
            delete!(data["branch"], b_id)
            continue
        end

        if "type" ∉ keys(branch)
            branch["type"] = "L"
        end
        if "name" ∉ keys(branch)
            branch["name"] = "$(t_bus_id)_$(f_bus_id)"
        end
    end

    branch_data_copy = copy(data["branch"]) # this is useful as we will be manipulating the branches

    # connect all busbars within a substation with couplers, if the busbars are supposed to be switchable
    n_couplers = 0
    current_branch_index = maximum(parse.(Int64, [branch_id for (branch_id, branch) in data["branch"]])) + 1
    if PARAMS[:SUBSTATION_RECONF][:COUPLER_RECONF]
        for (substation_name, substation) in data["substations"]
            for (busbar_i_idx1, busbar_i1) in enumerate(substation["busbar"][begin:end-1])
                if !is_bus_switchable(data["bus"][busbar_i1])
                    continue
                end

                for busbar_i2 in substation["busbar"][(busbar_i_idx1+1):end]
                    if !is_bus_switchable(data["bus"][busbar_i2])
                        continue
                    end

                    # if coupler is not switchable, we have continued before already
                    coupler_branch = _build_aux_branch(busbar_i1, busbar_i2, "K", substation_name, n_couplers, switchable=true, current_branch_index, status=0, br_status=1)
                    n_couplers += 1
                    current_branch_index += 1

                    push!(data["branch"], "$(coupler_branch["index"])" => coupler_branch)
                    push!(substation["K"], "$(coupler_branch["index"])")
                end
            end
        end
    end

    # replace all branch ends with auxiliary breaker buses, and connect those with the corresponding busbars
    n_breaker_branches = 0
    n_breaker_buses = 0
    n_removed_couplers = 0
    n_mislabeled_couplers = 0
    if PARAMS[:SUBSTATION_RECONF][:BRANCH_RECONF]
        for (b_id, branch) in branch_data_copy

            t_bus_id = string(branch["t_bus"])
            f_bus_id = string(branch["f_bus"])

            t_bus_name = data["bus"][t_bus_id]["name"]
            f_bus_name = data["bus"][f_bus_id]["name"]

            branch["breaker_branches"] = Vector{Int64}()

            # check if branch is in same substation (ie coupler branch), and remove it if so, but use its status
            if chop(t_bus_name) == chop(f_bus_name)
                if !is_bus_switchable(data["bus"][t_bus_id]) 
                    continue
                end

                for coupler_id in data["substations"][chop(t_bus_name)]["K"]
                    coupler = data["branch"][coupler_id]

                    # if ends are the same, independent of direction, this is the correct coupler
                    if string(coupler["f_bus"]) ∈ [t_bus_id, f_bus_id] && string(coupler["t_bus"]) ∈ [t_bus_id, f_bus_id]
                        if branch["transformer"]
                            data["branch"][coupler_id] = branch
                            data["branch"][coupler_id]["index"] = parse(Int64, coupler_id)
                            # coupler + transformer
                            data["branch"][coupler_id]["type"] = "K"
                            data["branch"][coupler_id]["transformer"] = true
                        else
                            if branch["transformer"]
                                @show(branch)
                            end
                            coupler[   "status"] = branch[   "status"]
                            coupler["br_status"] = branch["br_status"]
                            coupler["rate_a"]    = branch["rate_a"]
                        end

                        break
                    end

                end

                if branch["type"] ∉ ["K", "coupler", "T"]
                    n_mislabeled_couplers += 1
                    PARAMS[:PRINTING][:CONFIDENTIAL] && @warn("Removing intra-substation branch $(branch["name"]) of type $(branch["type"]) that is not labeled as coupler.")
                elseif ! branch["transformer"]
                    @debug("Removing intra-substation branch $(branch["name"]) of type $(branch["type"]).")
                end
                delete!(data["branch"], b_id)
                
                n_removed_couplers += 1
                continue
            end

            if is_branch_switchable(branch)
                was_replaced = false
                # copy f_bus and t_bus into a new auxiliary bus each, which represent breakers
                for (bus_id, other_bus_id) ∈ zip([t_bus_id, f_bus_id], [f_bus_id, t_bus_id])
                    bus       = data["bus"][      bus_id]
                    other_bus = data["bus"][other_bus_id]

                    substation_name       = substation_name_filter(      bus["name"])
                    other_substation_name = substation_name_filter(other_bus["name"])

                    if length(data["substations"][substation_name]["busbar"]) > 1
                        n_breaker_buses += 1
                        
                        ## build new bus
                        breaker_bus = deepcopy(bus)

                        # replace fields of original bus to fit the new bus
                        breaker_bus["bus_i"] = current_bus_index
                        breaker_bus["index"] = current_bus_index
                        breaker_bus["name"] = "$(substation_name)_$(length(get(data["substations"][substation_name], "aux_bus", []))+1):brk"
                        breaker_bus["source_id"] = ["aux", bus["index"]]
                        breaker_bus["type"] = "breaker"
                        breaker_bus["breaker_branches"] = Vector{Int64}()

                        breaker_bus_i = current_bus_index
                        
                        current_bus_index += 1

                        ## set ends of branch to new aux node
                        branch[(bus_id == t_bus_id) ? "t_bus" : "f_bus"] = breaker_bus_i

                        ## add to data and to substation
                        push!(data["bus"], "$breaker_bus_i"=>breaker_bus)
                        push!(data["substations"][substation_name]["aux_bus"], Dict("bus_i" => breaker_bus_i, "breaker_branches" => breaker_bus["breaker_branches"]))

                        ## add one auxiliary branch from the new breaker to each busbar of the corresponding substation
                        # these branches are only used to indicate connections, so we can omit most of the info contained in normal branches
                        for busbar_i in data["substations"][substation_name]["busbar"]
                            n_breaker_branches += 1

                            # if this is the connection that was connected for the original branch, use its status, otherwise 0
                            status = busbar_i ∈ [f_bus_id, t_bus_id] ? branch["status"] : 0
                            breaker_branch = _build_aux_branch("$busbar_i", "$breaker_bus_i", "breaker", String(other_substation_name), n_breaker_branches, current_branch_index, br_status = branch["br_status"], status=status, switchable=branch["switchable"], rate_a = branch["rate_a"], orig_branch_id=parse(Int64, b_id))
                            breaker_branch["original_branch"] = branch["index"]
                            current_branch_index += 1

                            branch_i = breaker_branch["index"]
                            
                            push!(branch["breaker_branches"], branch_i)
                            push!(breaker_bus["breaker_branches"], branch_i)
                            push!(data["branch"], "$branch_i" => breaker_branch)
                        end

                        was_replaced = true
                    end 
                end

                if was_replaced
                    branch["switchable"] = false
                    branch["status"] = 1
                end
            else
                branch["status"] = branch["br_status"]
            end
            if length(branch["breaker_branches"]) == 0
                branch["breaker_branches"] = [branch["index"]]
            end
        end
    end
    if n_mislabeled_couplers > 0
        @warn("Removed $n_mislabeled_couplers intra-substation branches that are not labeled as couplers.")
    end
    @info("Removed $n_removed_couplers existing couplers.
      Added $n_breaker_buses breaker buses, $n_breaker_branches breaker branches and $n_couplers couplers.\n\n")
end


# Modified from PowerModels.jl (https://github.com/lanl-ansi/PowerModels.jl), see LICENSE
"Export power network data to a file in the matpower format"
function export_matpower(file::AbstractString, data::Dict{String,Any})
    open(file, "w") do io
        export_matpower(io, data)
    end
end

# Modified from PowerModels.jl (https://github.com/lanl-ansi/PowerModels.jl), see LICENSE
"Export power network data in the matpower format"
function export_matpower(io::IO, data::Dict{String,Any})
    filter!(x -> x[1] ∉ _mp_non_export_columns, data)
    if _IM.ismultinetwork(data)
        Memento.error(_LOGGER, "export_matpower does not yet support multinetwork data")
    end

    data = deepcopy(data)

    #convert data to mixed unit
    if data["per_unit"]
       make_mixed_units!(data)
    end

    # make all costs have the name number of items (to prepare for table export)
    standardize_cost_terms!(data)

    # create some useful maps and data structures
    buses = Dict{Int, Dict}()
    for (idx,bus) in data["bus"]
        buses[bus["index"]] = bus
    end
    generators = Dict{Int, Dict}()
    for (idx,gen) in data["gen"]
        generators[gen["index"]] = gen
    end
    storage = Dict{Int, Dict}()
    for (idx,strg) in data["storage"]
        storage[strg["index"]] = strg
    end
    branches = Dict{Int, Dict}()
    for (idx,branch) in data["branch"]
       branches[branch["index"]] = branch
    end
    dclines = Dict{Int, Dict}()
    for (idx,dcline) in data["dcline"]
       dclines[dcline["index"]] = dcline
    end
    switches = Dict{Int, Dict}()
    for (idx,switch) in data["switch"]
       switches[switch["index"]] = switch
    end
    ne_branches = Dict{Int, Dict}()
    if haskey(data, "ne_branch")
        for (idx,branch) in data["ne_branch"]
            ne_branches[branch["index"]] = branch
        end
    end

    pd = Dict{Int, Float64}()
    qd = Dict{Int, Float64}()
    gs = Dict{Int, Float64}()
    bs = Dict{Int, Float64}()

    # collect all the loads
    for (idx,bus) in data["bus"]
        pd[bus["index"]] = 0
        qd[bus["index"]] = 0
    end
    for (idx,load) in data["load"]
        bus = buses[load["load_bus"]]
        pd[bus["index"]] += load["pd"]
        qd[bus["index"]] += load["qd"]
    end

    # collect all the shunts
    for (idx,bus) in data["bus"]
        gs[bus["index"]] = 0
        bs[bus["index"]] = 0
    end
    for (idx, shunt) in data["shunt"]
        bus = buses[shunt["shunt_bus"]]
        bs[bus["index"]] = bs[bus["index"]] + shunt["bs"]
        gs[bus["index"]] = gs[bus["index"]] + shunt["gs"]
    end

    mvaBase = data["baseMVA"]

    # Print the header information
    println(io, "%% MATPOWER Case Format : Version 2")
    println(io, "function mpc = ", data["name"])
    println(io, "mpc.version = '2';")
    println(io)
    println(io, "%%-----  Power Flow Data  -----%%")
    println(io, "%% system MVA base")
    println(io, "mpc.baseMVA = ", mvaBase, ";")

    # Print the bus data
    println(io, "%% bus data")
    println(io, "%    bus_i    type    Pd    Qd    Gs    Bs    area    Vm    Va    baseKV    zone    Vmax    Vmin")
    println(io, "mpc.bus = [")
    for (idx,bus) in sort(collect(buses), by=(x) -> x.first)
        println(io,
            "\t", bus["index"],
            "\t", bus["bus_type"],
            "\t", PowerModels._get_default(pd, bus["index"]),
            "\t", PowerModels._get_default(qd, bus["index"]),
            "\t", PowerModels._get_default(gs, bus["index"]),
            "\t", PowerModels._get_default(bs, bus["index"]),
            "\t", PowerModels._get_default(bus, "area"),
            "\t", PowerModels._get_default(bus, "vm"),
            "\t", PowerModels._get_default(bus, "va"),
            "\t", PowerModels._get_default(bus, "base_kv"),
            "\t", PowerModels._get_default(bus, "zone"),
            "\t", PowerModels._get_default(bus, "vmax"),
            "\t", PowerModels._get_default(bus, "vmin"),
            "\t", PowerModels._get_default(bus, "lam_p", ""),
            "\t", PowerModels._get_default(bus, "lam_q", ""),
            "\t", PowerModels._get_default(bus, "mu_vmax", ""),
            "\t", PowerModels._get_default(bus, "mu_vmin",""),
        )
    end
    println(io, "];")
    println(io)

    # Print the bus names
    if haskey(collect(values(buses))[1], "bus_name")
        println(io, "%% bus names")
        println(io, "mpc.bus_name = {")
        for (idx,bus) in sort(buses)
            println(io, "\t'", bus["bus_name"], "'")
        end
        println(io, "};")
        println(io)
    end

    # Print the generator data
    println(io, "%% generator data")
    println(io, "%    bus    Pg    Qg    Qmax    Qmin    Vg    mBase    status    Pmax    Pmin    Pc1    Pc2    Qc1min    Qc1max    Qc2min    Qc2max    ramp_agc    ramp_10    ramp_30    ramp_q    apf")
    println(io, "mpc.gen = [")
    i = 1
    for (idx,gen) in sort(collect(generators), by=(x) -> x.first)
        if idx != gen["index"]
            Memento.warn(_LOGGER, "The index of the generator does not match the matpower assigned index. Any data that uses generator indexes for reference is corrupted.");
        end
        println(io, 
            "\t", gen["gen_bus"],
            "\t", PowerModels._get_default(gen, "pg"),
            "\t", PowerModels._get_default(gen, "qg"),
            "\t", PowerModels._get_default(gen, "qmax"),
            "\t", PowerModels._get_default(gen, "qmin"),
            "\t", PowerModels._get_default(gen, "vg"),
            "\t", PowerModels._get_default(gen, "mbase"),
            "\t", PowerModels._get_default(gen, "gen_status"),
            "\t", PowerModels._get_default(gen, "pmax"),
            "\t", PowerModels._get_default(gen, "pmin"),
            "\t", PowerModels._get_default(gen, "pc1", ""),
            "\t", PowerModels._get_default(gen, "pc2", ""),
            "\t", PowerModels._get_default(gen, "qc1min", ""),
            "\t", PowerModels._get_default(gen, "qc1max", ""),
            "\t", PowerModels._get_default(gen, "qc2min", ""),
            "\t", PowerModels._get_default(gen, "qc2max", ""),
            "\t", PowerModels._get_default(gen, "ramp_agc", ""),
            "\t", (haskey(gen, "ramp_10") ? gen["ramp_10"] : haskey(gen, "ramp_30") ? 0 : ""),
            "\t", PowerModels._get_default(gen, "ramp_30", ""),
            "\t", PowerModels._get_default(gen, "ramp_q", ""),
            "\t", PowerModels._get_default(gen, "apf", ""),
            "\t", PowerModels._get_default(gen, "mu_pmax", ""),
            "\t", PowerModels._get_default(gen, "mu_pmin", ""),
            "\t", PowerModels._get_default(gen, "mu_qmax", ""),
            "\t", PowerModels._get_default(gen, "mu_qmin", ""),
        )
        i = i+1
    end
    println(io,"];")
    println(io)


    if length(storage) > 0
        # Print the storage data
        println(io, "%% storage data")
        println(io, "%    storage_bus    ps    qs    energy    energy_rating    charge_rating    discharge_rating    charge_efficiency    discharge_efficiency    thermal_rating    qmin    qmax    r    x    p_loss    q_loss    status")
        println(io, "mpc.storage = [")
        i = 1
        for (idx,strg) in sort(collect(storage), by=(x) -> x.first)
            if idx != strg["index"]
                Memento.warn(_LOGGER, "The index of the storage does not match the matpower assigned index. Any data that uses storage indexes for reference is corrupted.");
            end
            println(io, "\t", strg["storage_bus"],
                "\t", PowerModels._get_default(strg, "ps"),
                "\t", PowerModels._get_default(strg, "qs"),
                "\t", PowerModels._get_default(strg, "energy"),
                "\t", PowerModels._get_default(strg, "energy_rating"),
                "\t", PowerModels._get_default(strg, "charge_rating"),
                "\t", PowerModels._get_default(strg, "discharge_rating"),
                "\t", PowerModels._get_default(strg, "charge_efficiency"),
                "\t", PowerModels._get_default(strg, "discharge_efficiency"),
                "\t", PowerModels._get_default(strg, "thermal_rating"),
                "\t", PowerModels._get_default(strg, "qmin"),
                "\t", PowerModels._get_default(strg, "qmax"),
                "\t", PowerModels._get_default(strg, "r"),
                "\t", PowerModels._get_default(strg, "x"),
                "\t", PowerModels._get_default(strg, "p_loss"),
                "\t", PowerModels._get_default(strg, "q_loss"),
                "\t", PowerModels._get_default(strg, "status"),
            )
            i = i+1
        end
        println(io,"];")
        println(io)
    end


    # Print the branch data
    println(io, "%% branch data")
    println(io, "%    f_bus    t_bus    r    x    b    rateA    rateB    rateC    ratio    angle    status    angmin    angmax")
    println(io, "mpc.branch = [")
    i = 1
    for (idx,branch) in sort(collect(branches), by=(x) -> x.first)
        if get(branch, "type", "") ∈ ["breaker", "K"]; continue; end
        if idx != branch["index"]
            Memento.warn(_LOGGER, "The index of the branch does not match the matpower assigned index. Any data that uses branch indexes for reference is corrupted.");
        end
        println(io,
            "\t", PowerModels._get_default(branch, "f_bus"),
            "\t", PowerModels._get_default(branch, "t_bus"),
            "\t", PowerModels._get_default(branch, "br_r"),
            "\t", PowerModels._get_default(branch, "br_x"),
            "\t", (branch["b_to"] + branch["b_fr"]) != NaN ? (branch["b_to"] + branch["b_fr"]) : 0,
            "\t", PowerModels._get_default(branch, "rate_a"),
            "\t", PowerModels._get_default(branch, "rate_b"),
            "\t", PowerModels._get_default(branch, "rate_c"),

            "\t", (branch["transformer"] ? branch["tap"] : 0),
            "\t", (branch["transformer"] ? branch["shift"] : 0),
            "\t", PowerModels._get_default(branch, "br_status"),
            "\t", PowerModels._get_default(branch, "angmin"),
            "\t", PowerModels._get_default(branch, "angmax"),
            "\t", PowerModels._get_default(branch, "pf", ""),
            "\t", PowerModels._get_default(branch, "qf", ""),
            "\t", PowerModels._get_default(branch, "pt", ""),
            "\t", PowerModels._get_default(branch, "qt", ""),
            "\t", PowerModels._get_default(branch, "mu_sf", ""),
            "\t", PowerModels._get_default(branch, "mu_st", ""),
            "\t", PowerModels._get_default(branch, "mu_angmin", ""),
            "\t", PowerModels._get_default(branch, "mu_angmax", ""),
        )
        i = i+1
    end
    println(io, "];")
    println(io)

    if length(dclines) > 0 
        # print the dcline data
        println(io, "%% dcline data")
        println(io, "%    f_bus    t_bus    status    Pf    Pt    Qf    Qt    Vf    Vt    Pmin    Pmax    QminF    QmaxF    QminT    QmaxT    loss0    loss1")
        println(io, "mpc.dcline = [")
        for (idx, dcline) in sort(collect(dclines), by=(x) -> x.first)
            println(io,
                "\t", PowerModels._get_default(dcline, "f_bus"),
                "\t", PowerModels._get_default(dcline, "t_bus"),
                "\t", PowerModels._get_default(dcline, "br_status"),
                "\t", PowerModels._get_default(dcline, "pf"),
                "\t", -PowerModels._get_default(dcline, "pt"),  # opposite convention
                "\t", -PowerModels._get_default(dcline, "qf"),  # opposite convention
                "\t", -PowerModels._get_default(dcline, "qt"),  # opposite convention
                "\t", PowerModels._get_default(dcline, "vf"),
                "\t", PowerModels._get_default(dcline, "vt"),
                "\t", (haskey(dcline, "mp_pmin") ? dcline["mp_pmin"] : min(dcline["pmaxt"], dcline["pmint"], dcline["pmaxf"], dcline["pminf"])),
                "\t", (haskey(dcline, "mp_pmax") ? dcline["mp_pmax"] : max(dcline["pmaxt"], dcline["pmint"], dcline["pmaxf"], dcline["pminf"])),
                "\t", PowerModels._get_default(dcline, "qminf"),
                "\t", PowerModels._get_default(dcline, "qmaxf"),
                "\t", PowerModels._get_default(dcline, "qmint"),
                "\t", PowerModels._get_default(dcline, "qmaxt"),
                "\t", PowerModels._get_default(dcline, "loss0"),
                "\t", PowerModels._get_default(dcline, "loss1"),
                "\t", PowerModels._get_default(dcline, "mu_pmin", ""),
                "\t", PowerModels._get_default(dcline, "mu_pmax", ""),
                "\t", PowerModels._get_default(dcline, "mu_qminf", ""),
                "\t", PowerModels._get_default(dcline, "mu_qmaxf", ""),
                "\t", PowerModels._get_default(dcline, "mu_qmint", ""),
                "\t", PowerModels._get_default(dcline, "mu_qmaxt", ""),
            )
        end
        println(io, "];")
        println(io)
    end


    if length(switches) > 0
        # print the switch data
        println(io, "%% switch data")
        println(io, "%  f_bus    t_bus    psw    qsw    state    thermal_rating    status")
        println(io, "mpc.switch = [")
        for (idx, switch) in sort(collect(switches), by=(x) -> x.first)
            println(io,
                "\t", PowerModels._get_default(switch, "f_bus"),
                "\t", PowerModels._get_default(switch, "t_bus"),
                "\t", PowerModels._get_default(switch, "psw"),
                "\t", PowerModels._get_default(switch, "qsw"),
                "\t", PowerModels._get_default(switch, "state"),
                "\t", PowerModels._get_default(switch, "thermal_rating"),
                "\t", PowerModels._get_default(switch, "status")
            )
        end
        println(io, "];")
        println(io)
    end


    # Print the gen cost data
    PowerModels._export_cost_data(io, generators, "mpc.gencost")

    # Print the dcline cost data
    PowerModels._export_cost_data(io, dclines, "mpc.dclinecost")

    # ne branch is not part of the matpower specs. However, it is treated as a special case by the matpower parser
    # for example, br_b is converted into b_to and b_fr
    if haskey(data, "ne_branch")
        println(io, "%column_names%    f_bus    t_bus    br_r    br_x    br_b    rate_a    rate_b    rate_c    tap    shift    br_status    angmin    angmax    construction_cost")
        println(io, "mpc.ne_branch = [")
        i = 1
        for (idx,branch) in sort(collect(ne_branches), by=(x) -> x.first)
            if idx != branch["index"]
                Memento.warn(_LOGGER, "The index of the ne_branch does not match the matpower assigned index. Any data that uses branch indexes for reference is corrupted.");
            end
            println(io,
                "\t", branch["f_bus"],
                "\t", branch["t_bus"],
                "\t", PowerModels._get_default(branch, "br_r"),
                "\t", PowerModels._get_default(branch, "br_x"),
                "\t", (haskey(branch,"b_to") ? branch["b_to"] + branch["b_fr"]  : 0),
                "\t", PowerModels._get_default(branch, "rate_a"),
                "\t", PowerModels._get_default(branch, "rate_b"),
                "\t", PowerModels._get_default(branch, "rate_c"),
                "\t", (branch["transformer"] ? branch["tap"] : 0),
                "\t", (branch["transformer"] ? branch["shift"] : 0),
                "\t", PowerModels._get_default(branch, "br_status"),
                "\t", PowerModels._get_default(branch, "angmin"),
                "\t", PowerModels._get_default(branch, "angmax"),
                "\t", PowerModels._get_default(branch, "construction_cost"),
            )
            i = i+1
        end
        println(io, "];");
        println(io)
    end

    if all(haskey(bus, "name") for (i,bus) in buses)
        # Print the bus name data
        println(io, "%% bus names")
        println(io, "mpc.bus_name = {")
        for (idx,bus) in sort(collect(buses), by=(x) -> x.first)
            println(io,
                "\t", "'", bus["name"], "'"
            )
        end
        println(io, "};")
        println(io)
    end

    # Print the extra bus data
    PowerModels._export_extra_data(io, data, "bus", Set(["index", "source_id", "gs", "bs", "zone", "bus_i", "bus_type", "qd",  "vmax", "area", "vmin", "va", "vm", "base_kv", "pd", "name", "lam_p", "lam_q", "mu_vmax", "mu_vmin"]); postfix="_data")

    # Print the extra generator data
    PowerModels._export_extra_data(io, data, "gen", Set(["index", "source_id", "gen_bus", "pg", "qg", "qmax", "qmin", "vg", "mbase", "gen_status", "pmax", "pmin", "pc1", "pc2", "qc1min", "qc1max", "qc2min", "qc2max", "ramp_agc", "ramp_10", "ramp_30", "ramp_q", "apf", "ncost", "model", "shutdown", "startup", "cost", "mu_pmax", "mu_pmin", "mu_qmax", "mu_qmin"]); postfix="_data")

    # Print the extra storage data
    PowerModels._export_extra_data(io, data, "storage", Set(["index", "source_id", "storage_bus", "ps", "qs", "energy", "energy_rating", "charge_rating", "discharge_rating", "charge_efficiency", "discharge_efficiency", "thermal_rating", "qmin", "qmax", "r", "x", "p_loss", "q_loss", "status"]); postfix="_data")

    # Print the extra branch data
    PowerModels._export_extra_data(io, data, "branch", Set(["index", "source_id", "f_bus", "t_bus", "br_r", "br_x", "br_b", "b_to", "b_fr", "rate_a", "rate_b", "rate_c", "tap", "shift", "br_status", "angmin", "angmax", "transformer", "g_to", "g_fr", "pf", "qf", "pt", "qt", "mu_sf", "mu_st", "mu_angmin", "mu_angmax"]); postfix="_data")

    # Print the extra dcline data
    PowerModels._export_extra_data(io, data, "dcline", Set(["index", "source_id", "mu_qmaxt", "mu_qmint", "mu_qmaxf", "mu_qminf", "mu_pmax", "mu_pmin", "loss0", "loss1", "qmint", "qmaxt", "pmin", "pmax", "qminf", "qmaxf", "f_bus", "t_bus", "br_status", "pf", "pt", "qf", "qt", "vf", "vt", "ncost", "model", "shutdown", "pmaxt", "startup", "pmint", "cost", "pminf", "pmaxf", "mp_pmax", "mp_pmin"]); postfix="_data")

    # Print the extra switch data
    PowerModels._export_extra_data(io, data, "switch", Set(["index", "source_id", "f_bus", "t_bus", "psw", "qsw", "state", "thermal_rating", "status"]); postfix="_data")


    # Print the extra ne_branch data
    if haskey(data, "ne_branch")
        PowerModels._export_extra_data(io, data, "ne_branch", Set(["index", "source_id", "f_bus", "t_bus", "br_r", "br_x", "br_b", "b_to", "b_fr", "rate_a", "rate_b", "rate_c", "tap", "shift", "br_status", "angmin", "angmax", "transformer", "construction_cost", "g_to", "g_fr"]); postfix="_data")
    end

    # print the extra load data
    PowerModels._export_extra_data(io, data, "load", Set(["index", "source_id", "load_bus", "status", "qd", "pd"]); postfix="_data")

    # print the extra shunt data
    PowerModels._export_extra_data(io, data, "shunt", Set(["index", "source_id", "shunt_bus", "status", "gs", "bs"]); postfix="_data")

    # print the extra component data
    for (key, value) in data
        if key != "bus" && key != "gen" && key != "branch" && key != "load" && key != "shunt" && key != "storage" && key != "dcline" && key != "switch" && key != "ne_branch" && key != "version" && key != "baseMVA" && key != "per_unit" && key != "name" && key != "source_type" && key != "source_version"
            PowerModels._export_extra_data(io, data, key)
        end
    end
end