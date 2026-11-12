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

using Graphs, SimpleGraphs, SimpleGraphAlgorithms
import SimpleGraphAlgorithms.min_edge_cut
export get_critical_gens

### Getters ###
"get all branches that are connected to a bus"
function get_connected_branches(data, bus_id::Int64; only_active=false)::Vector{Int64}
  branch_ids = Vector{Int64}()
  # for each branch, check its f_bus and t_bus
  for (_, branch) in data["branch"]
    if ((branch["f_bus"] == bus_id) || (branch["t_bus"] == bus_id)) && (!only_active || get(branch, "status", branch["br_status"]) != 0)
      # this branch is connected to bus_id, so add it
      if data["BUSBARS"] && branch["type"] == "breaker"
        branch_ids = [branch_ids; branch["original_branch"]]
      else
        branch_ids = [branch_ids; branch["index"]]
      end
    end
  end
  
  #@assert length(branch_ids) > 0 || only_active
  return branch_ids
end

"get all (real) branches that are connected to each bus, for every bus. TODO: call get_connected_branches?"
function get_connected_branchess(pm::AbstractPowerModel; nw=nw_id_default, only_active=false, btc=Vector())::Dict{Int64, Set{Int64}}
  branch_idss = Dict(bus_id => Set{Int64}() for (bus_id, _) in ref(pm, nw, :bus))

  # for each branch, check its f_bus and t_bus
  for (_, branch) in ref(pm, nw, :branch)
    if !only_active || get_status(branch) != 0
      # this branch is active, so add it to connected branches of both ends
      considered_buses = [branch["f_bus"], branch["t_bus"]]
      for bus_id in considered_buses
        if pm.data["BUSBARS"] && branch["type"] == "breaker"
          push!(branch_idss[bus_id], branch["original_branch"])
        elseif get(branch, "type", "") ≠ "K"
          push!(branch_idss[bus_id], branch["index"])
        end
      end
    end
  end
  
  return branch_idss
end

"get all loads that are connected to a bus"
function get_connected_loads(data, bus_id::Int64; only_active=false)::Vector{Int64}
  load_ids = Vector{Int64}()
  # for each load, check its load_bus
  for load in data["load"]
    if (load[2]["load_bus"] == bus_id) && (!only_active || load[2]["pd"] > 0)
      # this load is connected to bus_id, so add it
      load_ids = [load_ids; load[2]["index"]]
    end
  end
  
  return load_ids
end

"get all gens that are connected to a bus"
function get_connected_gens(data, bus_id::Int64; only_active=false)::Vector{Int64}
  gen_ids = Vector{Int64}()
  # for each gen, check its gen_bus
  for gen in data["gen"]
    if (gen[2]["gen_bus"] == bus_id) && (!only_active || gen[2]["pg"] > 0)
      # this gen is connected to bus_id, so add it
      gen_ids = [gen_ids; gen[2]["index"]]
    end
  end
  return gen_ids
end

function get_branch_by_index(data, index::Int64)
  for (br_id, branch) in data["branch"]
    if branch["index"] == index
      return branch
    end
  end
  @warn("Could not find branch with index $index.")
end

"given a branch multicontingency, check if there exist any buses that have exactly these branches connected"
function get_induced_bus_contingencies(data, contingency; only_active=false)::Vector{Int64}
  exclude_buses = Vector{Int64}()
  if haskey(contingency, :branch) && length(contingency[:branch]) > 1 # branch multicontingency might include a bus
    connected_buses = Set()
    # get all buses that are at either end of all of the multicontingency branches
    for br_index in contingency[:branch]
      branch = get_branch_by_index(data, br_index)
      if only_active && get(branch, "status", branch["status"]) == 0
        continue
      end
      push!(connected_buses, branch["f_bus"])
      push!(connected_buses, branch["t_bus"])
    end
    for bus_index in connected_buses
      # check if the considered bus only has the branches of the multicontingency connected to it
      if get_connected_branches(data, bus_index, only_active=only_active) ⊆ contingency[:branch]
        push!(exclude_buses, bus_index)
      end
    end
  end
  return exclude_buses
end

"function to calculate radial branches and branches that would generate islands"
function get_critical_branches(data; nw=nw_id_default, check_bus_autonomy=false, contingency::Contingency=Contingency(), only_active=true)::Set{Contingency}
  if "nw" ∈ keys(data)
    data = data["nw"]["$nw"]
  end
  
  # create nodes, representing buses, and save what buses are considered
  g = Graph()
  contingency_branches = deepcopy(get(contingency, :branch, []))
  for contingency_branch in get(contingency, :branch, [])
    append!(contingency_branches, get(get_branch_by_index(data, contingency_branch), "breaker_branches", []))
  end
  induced_bus_contingencies = get_induced_bus_contingencies(data, Dict(:branch => contingency_branches))
  considered_buses = Dict{String, Any}()
  for (b_id::String, bus::Dict{String, Any}) in data["bus"]
    if bus["index"] ∉ get(contingency, :bus, []) && bus["index"] ∉ induced_bus_contingencies
      push!(considered_buses, b_id => bus)
      add_vertex!(g)
    end
  end

  # Graphs does not support naming nodes, so we have to create a map from the data bus ids (might be 1,2,100,9002,...) to the graph node ids (1,2,3,4,...)
  bus_id_data_to_graph = Dict(bus["index"] => new_bus_id for (new_bus_id, (bus_id, bus)) in enumerate(considered_buses))
  bus_id_graph_to_data = Dict(new_bus_id => bus["index"] for (new_bus_id, (bus_id, bus)) in enumerate(considered_buses))

  # create edges, representing branches
  for (br_index, branch) in data["branch"]
    if branch["index"] ∉ contingency_branches && (!only_active || get(branch, "status", branch["br_status"]) == 1)
      add_edge!(g, bus_id_data_to_graph[branch["f_bus"]], bus_id_data_to_graph[branch["t_bus"]])
    end
  end

  # generate plot of graph
  #@show(g, contingency)
  #display(graphplot(g, curves=false))

  # if graph is not connected, there already exist islands
  @assert(MetaGraphsNext.is_connected(g) || only_active, "Graph is not connected, i.e., case has islands.")

  # get bridges of the graph (edges that would create islands if removed)
  graph_bridges = [[src(n), dst(n)] for n in bridges(g)] 

  # TODO: check if radial branches can exist outside of the set of graph bridges
  radial_branches_candidates_edges = [(bus_id_graph_to_data[bridge[1]], bus_id_graph_to_data[bridge[2]]) for bridge in graph_bridges]
  radial_branches_candidates = [get_branch_id_from_fbus_tbus(data, fbus, tbus) for (fbus, tbus) in radial_branches_candidates_edges]

  # initialize list of true radial branches (and island-creating branches)
  radial_branches_ids = Set{Contingency}()

  # check if the outer end of the radial branch cannot provide for themselves
  for candidate_branch in radial_branches_candidates
    if ! check_bus_autonomy
      # define branch as radial
      push!(radial_branches_ids, Dict(:branch => [candidate_branch]))
    else
      # check the degree of the two connected nodes
      f_bus_id = data["branch"][string(candidate_branch)]["f_bus"]    
      t_bus_id = data["branch"][string(candidate_branch)]["t_bus"]
      f_bus_branches = get_connected_branches(data, f_bus_id)
      t_bus_branches = get_connected_branches(data, t_bus_id)
      f_bus_deg = length(f_bus_branches)
      t_bus_deg = length(t_bus_branches)
      # if two connected components with at least two nodes per component would be created ("islanding")
      # TODO: validate this assumption by checking the exact definition of islanding
      if (f_bus_deg > 1) && (t_bus_deg > 1)
        # case: this candidate would create two islands, and thus should not be considered as basic contingency
        Memento.trace(_LOGGER, "Branch $candidate_branch creates islands.")
        push!(radial_branches_ids, Dict(:branch => [candidate_branch]))
      elseif f_bus_deg == 1
        # case: only the currently viewed line is connected to f_bus
        # check if f_bus can provide for itself 
        if is_bus_autonomous(data, f_bus_id)
          Memento.trace(_LOGGER, "f_bus $f_bus_id of branch $candidate_branch can provide for itself.")
        else 
          Memento.trace(_LOGGER, "f_bus $f_bus_id of branch $candidate_branch cannot provide for itself, thus branch $candidate_branch is a radial line.")
          push!(radial_branches_ids, Dict(:branch => [candidate_branch]))
        end
      elseif t_bus_deg == 1
        # case: only the currently viewed line is connected to t_bus
        # check if t_bus can provide for itself 
        if is_bus_autonomous(data, t_bus_id)
          Memento.trace(_LOGGER, "t_bus $t_bus_id of branch $candidate_branch can provide for itself.")
        else 
          Memento.trace(_LOGGER, "t_bus $t_bus_id of branch $candidate_branch cannot provide for itself, thus branch $candidate_branch is a radial line.")
          push!(radial_branches_ids, Dict(:branch => [candidate_branch]))
        end
      end
    end
  end 

  return radial_branches_ids
end

"get total pd and qd over all loads and pmax and qmax over all gens and compare"
function is_bus_autonomous(data, bus_id)::Bool
  # get connected loads and gens for given bus
  load_ids = get_connected_loads(data, bus_id)
  gen_ids = get_connected_gens(data, bus_id)
  
  # get total pd and qd from loads
  pd_total, qd_total = 0, 0
  for load_id in load_ids
    pd_total += data["load"][load_id]["pd"]
    qd_total += data["load"][load_id]["qd"]
  end

  # get total max possible pg and qg from gens
  pg_max_total, qg_max_total = 0, 0
  for gen_id in gen_ids
    pg_max_total += data["gen"][gen_id]["pmax"]
    qg_max_total += data["gen"][gen_id]["qmax"]
  end
  
  # calculate balances 
  # TODO: check validity for reactive power flow (qg/qd)
  p_balance = pg_max_total - pd_total
  q_balance = qg_max_total - qd_total

  # check balance
  if (p_balance < 0) && (q_balance < 0)
    # demand cannot be covered by gens at this node
    return false
  else
    # demand can be covered by gens at this node
    return true
  end
end

"returns a selection of gens that cannot be turned off while remaining feasible"
function get_critical_gens(data; nw=nw_id_default, contingency::Contingency=Contingency(), active_only=true)::Set{Contingency}
  if "nw" ∈ keys(data)
    data = data["nw"]["$nw"]
  end

  # get total load
  total_balance_p     = 0
  total_balance_q_max = 0
  for (load_id, load) in data["load"]
    if load["source_id"][2] ∉ get(contingency, :bus, [])
      total_balance_p     -= load["pd"]
      total_balance_q_max -= load["qd"]
    end
  end

  # add maximum generation
  total_balance_q_min = total_balance_q_max
  for (gen_id, gen) in data["gen"]
    if gen["index"] ∉ get(contingency, :gen, [])
      total_balance_p     += gen["pmax"]
      total_balance_q_min += gen["qmin"]
      total_balance_q_max += gen["qmax"]
    end
  end

  # check if removing the gen will put us beyond the bound
  critical_gen_ids = Set{Contingency}()
  for (gen_id, gen) in data["gen"]
    if gen["index"] ∉ get(contingency, :gen, [])
      if  total_balance_p     - gen["pmax"] < 0                  ||
         (total_balance_q_min - gen["qmin"] > 0 && !active_only) ||
         (total_balance_q_max - gen["qmax"] < 0 && !active_only)
         
         push!(critical_gen_ids, Dict(:gen => [gen["index"]]))
      end
    end
  end

  return critical_gen_ids
end

"gets all supported critical components and checks their validity for the current implementation"
function get_critical_components(data_t, contingencies; n=nw_id_default)::Set{Contingency}
  critical_components = Set{Contingency}()

  # calculate critical components in case of any of the determined "contingencies that can be considered in the base case"
  for contingency in append!(contingencies, Contingency())
    union!(critical_components, get_critical_branches(data_t, nw=n, contingency=contingency)) # TODO: check if this works with multi-contingencies
    union!(critical_components, get_critical_gens(data_t, nw=n, contingency=contingency))
  end

  # check correctness of all critical components (compatibility with current, tested functionality)
  for component in critical_components
    @assert(length(component) == 1, "Contingencies with multiple component types not yet filtered out of decisions.")
    (component_type, component_ids) = first(component)
    @assert(length(component_ids) == 1, "Contingencies with multiple components per type not yet filtered out of decisions.")
  end

  return critical_components
end

"wrapper function to convert vector of critical branches to atomic contingencies"
function remove_branch_contingencies!(data, critical_branch_contingencies::Set{Int64})
    remove_branch_contingencies!(data, [Dict(:branch => [branch]) for branch in critical_branch_contingencies])
end

"given a list of contingencies with branches, remove critical branch contingencies from the data"
function remove_branch_contingencies!(data, critical_branch_contingencies::Vector{Contingency})
    n_before_removal = length(data["contingencies"])
    filter!(contingency -> !haskey(contingency, :branch) || contingency ∉ critical_branch_contingencies, data["contingencies"])
    @debug("[Critical Contingencies] Removed $(n_before_removal - length(data["contingencies"])) unfulfillable branch contingencies.")
    sort!(data["contingencies"], by=c->sum(values(c[:branch])))
end

"given a list of contingencies with generators, remove critical generator contingencies from the data"
function remove_gen_contingencies!(data, critical_gens_contingencies::Set{Contingency})
    n_before_removal = length(data["contingencies"])
    filter!(contingency -> !haskey(contingency, :gen) || contingency ∉ critical_gens_contingencies, data["contingencies"])
    @debug("[Critical Contingencies] Removed $(n_before_removal - length(data["contingencies"])) unfulfillable generator contingencies.")
end

"given a list of critical branches, make them not-switchable"
function fix_branches!(data, splitting_contingencies::Vector{Contingency})
    # the branches of one-element splitting contingencies can not be switched, otherwise the graph would be split
    n_line_fixed = 0
    n_breaker_fixed = 0

    # check each branch contained in a splitting contingency of size 1
    for contingency in splitting_contingencies
        length(contingency[:branch]) > 1 && continue

        for branch_id in contingency[:branch]
            branch = data["branch"]["$branch_id"]

            # fix if switchable
            if get(branch, "switchable", false) == true # explicit type conversion to Bool here in case the value is non-boolean TODO remove me later
                n_line_fixed += 1
                branch["switchable"] = false
            end

            # the associated breaker branches cannot be switched either
            if haskey(branch, "breaker_branches")
                for breaker_branch_id in branch["breaker_branches"]
                    breaker_branch = data["branch"]["$breaker_branch_id"]

                    # fix if switchable
                    if breaker_branch["switchable"] == true # explicit type conversion to Bool here in case the value is non-boolean TODO remove me later
                        n_breaker_fixed += 1
                        breaker_branch["switchable"] = false
                    end
                end
            end
        end
    end
    @debug("[Critical Contingencies] Fixed $n_line_fixed lines & $n_breaker_fixed breaker branches.")
end

"given the data, fix generation limits and status of generators"
function fix_gens!(data)
  # TODO: fix_critical_gens (fix required gens to remain on, reduce pmax based on single connected lines, ...)
  @warn("Not implemented")
end

"returns all basic, fulfillable contingencies. Also fixes splitting branches and critical gens."
function get_basic_contingencies(data;
                                 component_types::Vector{Symbol}= PARAMS[:CONTINGENCIES][:COMPONENT_TYPES], 
                                 exclude_critical_branches::Bool= true, 
                                 exclude_critical_gens::Bool    = true, 
                                 merge_substations::Bool        = false, 
                                 fix_splitting_branches::Bool   = true,
                                 fix_critical_gens::Bool        = false
                                )::Vector{Contingency}

  contingencies = Vector{Contingency}()
  if !haskey(data, "contingencies")
    for component_type in component_types
      contingencies_component = [Dict(component_type => [parse(Int64, component_id)]) for (component_id, component_data) in data[string(component_type)] if get(component_data, "type", "") ∉ ["breaker", "K"]]
      contingencies = [contingencies; contingencies_component]
      data["contingencies"] = contingencies
    end
  else
    contingencies = data["contingencies"]
    @info("Using given contingencies ($(length(contingencies)) of $(length(data["branch"])) branch contingencies).")
  end

  return preprocess_contingencies_and_switchability!(data, contingencies; 
                                 component_types=component_types, 
                                 exclude_critical_branches=exclude_critical_branches, 
                                 exclude_critical_gens=exclude_critical_gens, 
                                 merge_substations=merge_substations, 
                                 fix_splitting_branches=fix_splitting_branches,
                                 fix_critical_gens=fix_critical_gens
                                )
end

"returns all basic, fulfillable contingencies. Also fixes splitting branches and critical gens."
function preprocess_contingencies_and_switchability!(data,
                                 contingencies::Vector{Contingency}; 
                                 component_types::Vector{Symbol}= PARAMS[:CONTINGENCIES][:COMPONENT_TYPES], 
                                 exclude_critical_branches::Bool= true, 
                                 exclude_critical_gens::Bool    = true, 
                                 merge_substations::Bool        = false, 
                                 fix_splitting_branches::Bool   = true,
                                 fix_critical_gens::Bool        = false
                                )::Vector{Contingency}
    if exclude_critical_branches || fix_splitting_branches
        _, splitting_contingencies = get_splitting_contingencies(data, contingencies, merge_substations=merge_substations)
    
        if exclude_critical_branches
            remove_branch_contingencies!(data, splitting_contingencies)
            remove_branch_contingencies!(data, collect(get_critical_branches(data)))
        end


        if fix_splitting_branches
            fix_branches!(data, splitting_contingencies)
        end
    end

    if exclude_critical_gens
        critical_gens_contingencies = get_critical_gens(data)
        remove_gen_contingencies!(data, critical_gens_contingencies)
    end

    if fix_critical_gens
        fix_gens!(data)
    end

    return contingencies
end

"set flags for branches and gens such that they are forced on and cannot be deactivated"
function set_critical_component_status_switchable!(data)
  # initialize set of critical components
  critical_components = get_critical_components(data, data["contingencies"])

  # set switchable and status for all components that cannot be turned off (since a contingency based on the changed topology would then create islands)
  for component in critical_components
    (component_type, (component_id, _...)) = first(component)
    # set switchable and current status
    data[string(component_type)][string(component_id)]["switchable"]     = false
    data[string(component_type)][string(component_id)]["status"] = 1
  end

  # print info on how many gens and branches were determined to be N-2 critical
  Memento.info(_LOGGER, "Fixed $(count(x -> haskey(x, :branch), critical_components)) branches (of $(length(data["branch"])) branches) to on and non-switchable.")
  Memento.info(_LOGGER, "Fixed $(count(x -> haskey(x, :gen), critical_components)) gens (of $(length(data["gen"])) gens) to on and non-switchable.")
end

"helper function to set flags for multinetworks"
function set_critical_component_status_switchable_over_time!(data, timesteps; kwargs...)
  for time_step in 1:timesteps
    # set flags for branches and gens such that they are forced on and cannot be deactivated
    set_critical_component_status_switchable!(data["nw"]["$time_step"]; kwargs...)
  end
end

"find all critical components and remove them from the set of contingencies and set their status and switchable flags"
function calculate_contingencies!(data; kwargs...)
  # set contingencies to be considered in that timestep (excludes critical branches and gens!)
  data["contingencies"] = get_basic_contingencies(data; kwargs...)
  
  return
end

"helper function to calculate contigencies for multinetworks"
function calculate_contingencies_over_time!(data, timesteps; kwargs...)
  for time_step in 1:timesteps
    # get correct data and calculate basic contingencies for that timestep
    calculate_contingencies!(data["nw"]["$time_step"]; kwargs...)
  end
end

### connectivity
"returns a min cut for `graph`"
function _get_min_cut(graph::AbstractGraph)::Vector{Tuple{String, String}}
  # create simple graph from graph (needed for min_cut...)
  s_graph = IntGraph()
  for label in labels(graph)
    add!(s_graph, code_for(graph, label))
  end
  add_edges!(s_graph, [(code_for(graph, l1), code_for(graph, l2)) for (l1, l2) in edge_labels(graph)])
  cost_func = get_weight_function(graph)
  costs = Dict(((code_for(graph, l1), code_for(graph, l2)), cost_func(graph[(l1,l2)...])) for (l1,l2) in edge_labels(graph))

  return [(label_for(graph, c1), label_for(graph, c2)) for (c1,c2) in min_edge_cut(s_graph, costs)]
end

"return all buses that are not connected to the main network (the largest connected component)"
function get_island_buses(graph::AbstractGraph)
  island_buses        = Vector{String}()
  island_buses_labels = Vector{String}()
  components          = connected_components(graph)

  largest_component_idx = argmax((c -> length(components[c])), 1:length(components))
  for (c_id, component) in enumerate(components)
    if c_id == largest_component_idx
      continue
    end
    for v_id in component
      push!(island_buses_labels, label_for(graph, v_id))
      append!(island_buses, graph[label_for(graph, v_id)])
    end
  end

  return island_buses_labels, island_buses
end

"removes a given vector of buses, and a given set of branches from the substation data"
function remove_branches_and_buses_from_substation_data!(data, island_buses::Vector{String}, deleted_branches::Set{Int64})
  deleted_substations = 0
  for (s, substation) in data["substations"]
      #@show(s, substation["busbar"], substation["K"], substation["aux_bus"])
      filter!(busbar  -> busbar ∉ island_buses, substation["busbar"])
      filter!(coupler -> parse(Int64, coupler) ∉ deleted_branches, substation["K"])
      filter!(aux_bus -> string(aux_bus["bus_i"]) ∉ island_buses, substation["aux_bus"])
      for aux_bus in substation["aux_bus"]
          filter!(branch_id -> branch_id ∉ deleted_branches, aux_bus["breaker_branches"])
      end
      if isempty(substation["aux_bus"]) && isempty(substation["K"])
          @debug("Deleted substation $s")
          delete!(data["substations"], s)
          deleted_substations += 1
      end
  end

  for (branch_id, branch) ∈ data["branch"]
    if haskey(branch, "breaker_branches")
        filter!(branch_id -> branch_id ∉ deleted_branches, branch["breaker_branches"])
    end
  end
  @debug("[Critical Contingencies] Removed $deleted_substations substations.")
end

"return the graph consisting of only the largest island, and all buses that were not connected to it"
function _reduce_to_largest_island(data::Dict{<:Union{String, Symbol}, Any})::Tuple{AbstractGraph, Vector{String}}
  # get graph (with all connected components) and detect islands
  off_branches = [branch_id for (branch_id::String, branch::Dict{String, Any}) in data["branch"] if branch["status"] ≠ 1]
  graph = _create_underlying_graph(data, merge_substations=false, off_branches=off_branches, only_largest_component=false)
  island_buses_labels, island_buses_b_ids = get_island_buses(graph)

  # remove island_buses from graph
  for island_bus_label in island_buses_labels
    @assert(rem_vertex!(graph, code_for(graph, island_bus_label)), "Could not remove bus $island_bus_code from graph.")
  end
  @assert(MetaGraphsNext.is_connected(graph))
  @debug("[Largest Island] Identified $(length(island_buses_labels)) buses not connected to the main network.")
  return graph, island_buses_b_ids
end

"remove buses and loads, gens and branches connected to them"
function _remove_buses!(pm_data, buses_to_remove::Vector{String}; remove_outside_D=false)
    removed_buses       = 0
    removed_loads       = 0
    removed_gens        = 0
    removed_branches    = 0
    deleted_branches    = Vector{String}()

    # delete buses that should be removed
    for (bus_id, bus) in pm_data["bus"]
        if bus_id ∈ buses_to_remove
            delete!(pm_data["bus"], bus_id)
            removed_buses += 1
            # if the bus had a load connected to it, remove the load as well
            for (load_id, load) in pm_data["load"]
                if load["load_bus"] == bus["index"]
                    delete!(pm_data["load"], load_id)
                    removed_loads += 1
                end
            end
            # if the bus had a gen connected to it, remove the gen as well
            for (gen_id, gen) in pm_data["gen"]
                if gen["gen_bus"] == bus["index"]
                    #@info("Removing gen $(gen["name"]) connected to bus $(bus["name"]): $gen")
                    delete!(pm_data["gen"], gen_id)
                    removed_gens += 1
                end
            end
            # if the bus had a branch connected to it, remove the branch as well
            for (branch_id, branch) in pm_data["branch"]
                if branch["f_bus"] == bus["index"] || branch["t_bus"] == bus["index"]
                    #@info("Removing branch $(branch["branch_name_intern"]) connected to bus $(bus["name"]): $branch")
                    delete!(pm_data["branch"], branch_id)
                    push!(deleted_branches, branch_id)
                    removed_branches += 1
                end
            end
        end
    end

    @debug("[Bus Removal] Removed $removed_buses buses, $removed_loads loads, $removed_gens gens, and $removed_branches branches.")
    return Set(parse.(Int64, deleted_branches))
end

"searches for all `contingencies` that would split the network given by `data` by computing min cuts and then iteratively contracting violating edges.

`merge_substations`: merges connected parts of busbars (only useful if busbar expansion is applied)
`only_active`: only branches that are switched on are added to the graph
`only_powered`: splitting of components that are not inherently powered (ie do not contain loads) is ignored, it does not lead to contingencies to be considered splitting"
function get_splitting_contingencies(data::Dict{String, Any}, contingencies::Vector; merge_substations::Bool=false, only_active::Bool=true, only_powered::Bool=false)::Tuple{Vector{Contingency}, Vector{Contingency}}
  if only_active
    off_branches = [branch_id for (branch_id::String, branch::Dict{String, Any}) in data["branch"] if isapprox(branch["status"], 0.0, atol =SOLVER_PRECISION)]
  else
    off_branches = Vector{String}()
  end
  graph = _create_underlying_graph(data, merge_substations=merge_substations, off_branches=off_branches)
  
  return get_splitting_contingencies(graph, contingencies, data, only_powered=only_powered)
end

"searches for all `contingencies` that would split `graph` by computing min cuts and then iteratively contracting violating edges.

`only_powered`: splitting of components that are not inherently powered (ie do not contain loads) is ignored, it does not lead to contingencies to be considered splitting"
function get_splitting_contingencies(in_graph::AbstractGraph, contingencies::Vector{Contingency}, data; only_powered::Bool=false)::Tuple{Vector{Contingency}, Vector{Contingency}}
  @assert(MetaGraphsNext.is_connected(in_graph))

  open_contingencies::Vector{Contingency} = copy(contingencies)
  splitting_contingencies = Vector{Contingency}()
  graph = in_graph

  # while open contingencies remain
  # 1. contract all edges not contained in any contingency
  # 2. compute a min cut
  # 3. remove all contingencies with size smaller than the cut
  # 4. all contingencies that contain the cut are splitting -> remove them as well
  # 5. if none of the edges could be contracted in 1., we need to remove some contigencies:
  #   a. select the edge in the cut that is contained in the fewest number of contingencies
  #   b. for each of these contingencies, explicitly check whether they are splitting and then remove them
  # 6. go to 1.
  # the first iteration, we check all contingencies of size 1 by computing all bridges instead, since this is somewhat cheaper than the min cuts
  any_update=true
  bridges_checked = false
  while !isempty(open_contingencies) && any_update
    any_update = false

    if !bridges_checked
      # compute all bridges, as they are equal to all min cuts of size 1
      cuts = [[(label_for(graph, src(edge)), label_for(graph, dst(edge)))] for edge in MetaGraphsNext.bridges(graph)]::Vector{Vector{Tuple{String, String}}} # converts the output format into our label format

    else
      # 1.

      update = true
      while update
        update = false
        for edge_label in edge_labels(graph)
          new_branches = Vector() # branches that are not contained

          # check if contained
          for branch_id in graph[edge_label...]

            contained = false
            for contingency in open_contingencies
              if parse(Int64, branch_id) ∈ contingency[:branch]
                contained = true
                break
              end
            end

            # add to new_branches if contained
            if contained
              push!(new_branches, branch_id)
            end
          end

          # update, if empty contract
          if length(new_branches) < length(graph[edge_label...])
            graph[edge_label...] = new_branches
            any_update = true
          end
          if isempty(new_branches)
            update = true
            any_update = true

            # we contract l2 into l1
            l1, l2 = edge_label

            # remove the edge
            rem_edge!(graph, code_for(graph, l1), code_for(graph, l2))

            # re-attach all edges attached to l2 to l1
            for n_label in deepcopy(neighbor_labels(graph, l2)) # deepcopy is necessary
              if haskey(graph, (l1, n_label)...) # if there also is an edge between l1 & n_label, merge the edges
                append!(graph[l1, n_label], graph[l2, n_label])
              else # re-attach
                graph[l1, n_label] = graph[l2, n_label]
              end

              rem_edge!(graph, code_for(graph, l2), code_for(graph, n_label))
            end

            # add the to l2 associated buses to l1
            append!(graph[l1], graph[l2])

            rem_vertex!(graph, code_for(graph, l2))
            break
          end
        end
      end
      
      # 2.
      cuts = [_get_min_cut(graph)]
    end
      
    for cut in cuts

      # 3.
      cut_length = sum([length(graph[edge_label...]) for edge_label in cut])

      if cut_length > 1
        still_open_contingencies = copy(open_contingencies)
        filter!(contingency -> length(contingency[:branch]) ≥ cut_length, open_contingencies)
      else
        still_open_contingencies = Vector{Contingency}()  
      end

      # 4.
      for contingency in open_contingencies
        
        # each edge in the cut needs to be contained
        contained = true
        for edge_label in cut
          
          # the edge is contained if at least one of its associated branches is contained in the contingency
          for branch_id in (graph[edge_label...])::Vector{String} # all branches associated with this edge
            if parse(Int64, branch_id::String) ∉ contingency[:branch]
              contained = false
              break
            end
          end

          if !contained
            break
          end

        end

        # if the cut is contained in the contingency, it is splitting, otherwise we still don't know
        push!((contained ? splitting_contingencies : still_open_contingencies), contingency)
      end

      if length(still_open_contingencies) < length(open_contingencies)
        any_update = true
      end
      open_contingencies = still_open_contingencies
    end

    # all contingencies of size 1 that are not contained in a bridge can be removed
    if !bridges_checked
      open_contingencies = filter!(x->length(x[:branch]) ≥ 2, open_contingencies)

      if !isempty(open_contingencies)
        bridges_checked = true
        graph = deepcopy(in_graph)
      end
    end

    ## 5.
    if !any_update && !isempty(open_contingencies)
      # a
      min_count = length(open_contingencies) + 1 # the current min count of contingencies covering a label
      min_contingencies = Vector() # the contingencies containing the edge with the min count

      # find an edge_label with the lowest count of occurences
      for edge_label in cuts[1]
        cur_count = 0
        cur_contingencies = Vector()

        # check for each contingency whether it contains at least one of the associated branch_ids
        for contingency in open_contingencies

          for branch_id in graph[edge_label...]
            if parse(Int64, branch_id) ∈ contingency[:branch]
              cur_count += 1
              push!(cur_contingencies, contingency)

              break
            end
          end

        end

        # update min_...
        if count < min_count
          min_count = count
          min_contingencies = cur_contingencies

          # there is always at least one
          if min_count == 1
            break
          end
        end

      end

      @assert(min_count > 0)

      # b
      for contingency in min_contingencies
        removed_edges = Vector{Tuple{Tuple{String, String}, Vector{String}}}() # ((label1, label2), [branch_ids])

        # remove all edges of this contingency from the graph
        for edge_label in edge_labels(graph)
          # only if all branches of the edge are contained in the contingency, remove the edge
          if length([branch_id for branch_id in graph[edge_label...] if (parse(Int64, branch_id) ∈ contingency[:branch])]) == length(graph[edge_label...])
            push!(removed_edges, (edge_label, graph[edge_label...]))
            rem_edge!(graph, edge_label)
          end
        end

        # if the remaining graph is no longer connected, this contingency is splitting
        if !MetaGraphsNext.is_connected(graph)
          push!(splitting_contingencies, contingency)
        end

        # re-add the removed edges
        for (edge_label, edge_data) in removed_edges
          graph[edge_label] = edge_data
        end
      end

      # all contingencies here have been visited
      open_contingencies = [contingency for contingency in open_contingencies if contingency ∉ min_contingencies]
      any_update = true
    end
  end


  # remove contingencies where the split does not lead to at least two powered components
  powered_splitting_contingencies = Vector{Contingency}()
  if only_powered
    for contingency in splitting_contingencies
      powered_splitting_contingencies = Vector{Contingency}()
      removed_edges = Vector{Tuple{Tuple{String, String}, Vector{String}}}() # ((label1, label2), [branch_ids])

      # remove all edges in the contingency
      for (e1,e2) in edge_labels(graph)
        if graph[e1,e2] ⊆ string.(contingency[:branch])
          push!(removed_edges, ((e1,e2), graph[e1,e2]))
          rem_edge!(graph, code_for(graph, e1), code_for(graph, e2))
        end
      end

      components = connected_components(graph)

      # check for each component whether it is powered until we found at least two
      n_powered = 0
      for (c, connected_component) in enumerate(components)
        if _network_powered(data, bus_list=append!([parse.(Int64, graph[label]) for label in label_for.([graph], connected_component)]...))
            n_powered += 1

            if n_powered ≥ 2
              push!(powered_splitting_contingencies, contingency)
              break
            end
        end

        # if we need more powered components than are unchecked, we can stop
        if c + 2 - n_powered > length(components)
          break
        end
      end

      # re-add the removed edges
      for (edge_label, edge_data) in removed_edges
        graph[edge_label...] = edge_data
      end

    end
  end

  return powered_splitting_contingencies, splitting_contingencies
end

"calculates connected components within each substation"
function calc_busbar_components(data, off_branches::Vector{String}=Vector{String}())
    componentss = Dict{String, Vector{Vector{String}}}()
    for (substation_name, substation) in data["substations"]::Dict{String, Dict{String, Vector{Any}}}
        # get all coupler branches and breaker branches in that substation
        open_branches = [branch for branch::String in substation["K"]::Vector{Any} if branch ∉ off_branches]::Vector{String} # unvisited couplers
        for aux_bus in substation["aux_bus"] # unvisited breakers
            for branch_id in aux_bus["breaker_branches"]
                if "$branch_id" ∉ off_branches
                    push!(open_branches, "$branch_id")
                end
            end
        end
        
        # each component contains buses of the substation that are connected
        if length(open_branches) > 0
            components = Vector{Vector{String}}()
        else # no branches ⇒ no connection
            components::Vector{Vector{String}} = [substation["busbar"]::Vector{Any}]
        end

        # fill up each component until it is complete, then proceed to the next one
        while !isempty(open_branches)
            # initialize the component with an arbitrary branch
            branch_id = pop!(open_branches)
            branch = data["branch"][branch_id]

            push!(components, ["$(branch["f_bus"])", "$(branch["t_bus"])"])
            
            visited_branches = Vector{String}() # branches visited in the loop

            # while there are updates to the current component, continue to look for branches incident to it
            new = true
            while new
                new = false
                # check if the branch is incident, and if it introduces a new vertex into the component
                for branch_id in open_branches
                    branch = data["branch"][branch_id]

                    if "$(branch["f_bus"])" ∈ components[end]
                        if "$(branch["t_bus"])" ∉ components[end] # both are contained ⇒ connected, but no update; only f_bus incident ⇒ update
                            push!(components[end], "$(branch["t_bus"])")
                            
                            new = true
                        end

                        push!(visited_branches, branch_id)
                    elseif "$(branch["t_bus"])" ∈ components[end] # only t_bus incident ⇒ update
                        push!(components[end], "$(branch["f_bus"])")

                        push!(visited_branches, branch_id)
                        new = true
                    end
                end

                # remove visited branches from open_branches
                open_branches = [branch_id for branch_id in open_branches if branch_id ∉ visited_branches]
            end
        end

        componentss[substation_name] = collect(Set(components))
    end

    return componentss
end

"create underlying graph of the power system"
function _create_underlying_graph(data::Dict{<:Union{String, Symbol}, Any}; merge_substations::Bool=false, merge_non_switchable::Bool=false, off_branches::Vector{String}=Vector{String}(), only_largest_component::Bool=true)::AbstractGraph
    # in this function, branch & bus refer to entries in pm.data and edge & vertex to components of the created graph

    graph = MetaGraph(
        Graph();
        label_type=String,
        edge_data_type=Vector{String},     # associated `br_id`s  (since we can't include multi-edges)
        vertex_data_type=Vector{String},   # associated `bus_id`s (each vertex correponds to a connected part of a busbar, or buses that are connected by non-switchable branches) 
        weight_function=ed -> length(ed), # the weight of an edge is equal to the number of associated branches,
        default_weight=0
    )

    bus_to_label = Dict{String, String}() # bus_id -> label of vertex

    ## merge the connected busbars in each substation
    if merge_substations && haskey(data, "substations")
        componentss = calc_busbar_components(data, off_branches)
        for (substation_name, _) in data["substations"]
            components = componentss[substation_name]

            # add components as vertices
            for (c, component) in enumerate(components)
                @assert(MetaGraphsNext.add_vertex!(graph, "$(substation_name)_$c", component), "Could not add vertex $component for substation $substation_name.")
                
                for bus_id in component
                    bus_to_label["$bus_id"] = "$(substation_name)_$c"
                end
            end
        end
    end

    # add all buses that are not already covered as vertices
    for (bus_id, bus) ∈ (haskey(data, "bus") ? data["bus"] : data[:bus])
        if !merge_substations || !haskey(bus_to_label, "$bus_id") # don't add if already covered
            @assert(MetaGraphsNext.add_vertex!(graph, "$bus_id", ["$bus_id"]), "Could not add bus $bus_id: $bus")
            
            bus_to_label["$bus_id"] = "$bus_id"
        end
    end

    # merge vertices connected by non-switchable branches
    if merge_non_switchable
        for (br_id, branch) ∈ (haskey(data, "branch") ? data["branch"] : data[:branch])
            if branch["br_status"] == 1 && "$br_id" ∉ off_branches && !branch["switchable"]
                f_vertex = bus_to_label["$(branch["f_bus"])"]
                t_vertex = bus_to_label["$(branch["t_bus"])"]

                # merge vertices if not already merged
                if f_vertex != t_vertex
                    # transfer associated buses from second to first vertex
                    append!(graph[f_vertex], graph[t_vertex])

                    for bus_id in graph[t_vertex]
                        bus_to_label[bus_id] = f_vertex
                    end

                    # remove the second vertex
                    @assert(rem_vertex!(graph, code_for(graph, t_vertex)))
                end
            end
        end
    end

    # add edges
    for (br_id::Union{String, Int}, branch::Dict) in (haskey(data, "branch") ? data["branch"] : data[:branch])::Dict
        if branch["br_status"]::Int64 == 1 && "$br_id" ∉ off_branches && !(merge_substations && branch["type"]::String ∈ ["breaker", "K"]) && !(merge_non_switchable && !branch["switchable"]) # ignore branches already used for merging
            f_vertex = bus_to_label["$(branch["f_bus"])"]
            t_vertex = bus_to_label["$(branch["t_bus"])"]

            if f_vertex != t_vertex # might be a single vertex due to merging
                if MetaGraphsNext.has_edge(graph, code_for(graph, f_vertex), code_for(graph, t_vertex))
                    @debug("Graph already has edge $(branch["f_bus"]) -> $(branch["t_bus"]).")
                    push!(graph[f_vertex, t_vertex], "$br_id") # multi-edge
                else
                    @assert(MetaGraphsNext.add_edge!(graph, f_vertex, t_vertex, ["$br_id"]), "Could not add edge $(branch["f_bus"]) -> $(branch["t_bus"]) ($br_id).")
                end
            end
        end
    end

    # remove all but a biggest connected components
    if only_largest_component
        con_components = [label_for.([graph], component) for component in (connected_components(graph))]
        _, maxidx = findmax((x->length(x)).(con_components))
        for (idx, component) in enumerate(con_components)
            if idx ≠ maxidx
                for bus in component
                    @assert(rem_vertex!(graph, code_for(graph, bus)))
                end
            end
        end
    end
    
    return graph
end

"use steiner tree to check number of edges needed to connect the underlying graph as induced by the current master switching decisions given by `disabled_branches`"
function _underlying_graph_steiner_tree_approx(pm, master_var_vals, graph::AbstractGraph, disabled_branches_contingency::Vector{String}, disabled_branches::Vector{Vector{String}}; only_powered::Bool=true, nw::Int64=nw_id_default)::Int64
    # TODO: remember removed edges and re-add them after checking connectedness
    g = deepcopy(graph)

    for (e1,e2) in edge_labels(g)
      for b_id in disabled_branches_contingency
        if b_id ∈ graph[e1,e2]
          @assert(rem_edge!(g, code_for(g, e1), code_for(g, e2)), "Could not remove edge $((e1,e2)).")
        end
      end
    end

    # check all connected components except for the main network for power demand. If there is none, this is just an artifact and all demand is still covered
    components = connected_components(g)
    islands = Vector{Vector{Int64}}()
    connected_components_wo_islands = Vector{Vector{Int64}}() # components of size > 1 that are not powered
    if length(components) > 1
      for connected_component in components
            if ! only_powered || _network_powered(pm.data["nw"]["$nw"]["bus"], pm.data["nw"]["$nw"]["load"], pm.data["nw"]["$nw"]["gen"], bus_list=append!([parse.(Int64, g[label]) for label in label_for.([g], connected_component)]...))
                push!(islands, connected_component)

            elseif length(connected_component) > 1
                push!(connected_components_wo_islands, connected_component)
            end

        end
    end
    

    # build the steiner tree graph
    if length(islands) > 1
        return 1 # FIXME: see #215
        st_graph = MetaGraph(
            Graph();
            label_type=Int64,
            edge_data_type=String,
            vertex_data_type=Int64
        )

        # add a vertex for each island/other connected component, and keep track of the contained buses
        bus_id_to_v_id = Dict()
        vertex_counter = 1
        for connected_component_vector in [islands, connected_components_wo_islands]
            for connected_component in connected_component_vector
                @assert(MetaGraphsNext.add_vertex!(st_graph, vertex_counter, vertex_counter))

                # for each bus of this connected component, map from the original bus id in graph `g` to the connected component's id (which is `vertex_counter`)
                for bus in connected_component
                    bus_id_to_v_id[label_for(g, bus)] = vertex_counter

                    for ass_bus_id in g[label_for(g, bus)]
                      bus_id_to_v_id[ass_bus_id] = vertex_counter
                    end
                end

                vertex_counter += 1
            end
        end

        # for each disabled branch, add an edge between the nodes, or the respective island
        for (b_id, fbus, tbus) in disabled_branches 
            # add keys to dict
            # for each bus included as endpoint of any disabled branch, check if it was already (partly) in one of the connected components. If not, add it as vertex to the st_graph.
            for bus in filter((x->!haskey(bus_id_to_v_id, x)), [fbus, tbus])
                # the bus is not contained in a connected component, so we create a new vertex
                @assert(MetaGraphsNext.add_vertex!(st_graph, vertex_counter, vertex_counter))
                bus_id_to_v_id[bus] = vertex_counter
                vertex_counter += 1
            end

            # v1 is the smaller one
            v1, v2 = bus_id_to_v_id[fbus] < bus_id_to_v_id[tbus] ? (bus_id_to_v_id[fbus], bus_id_to_v_id[tbus]) : (bus_id_to_v_id[tbus], bus_id_to_v_id[fbus])

            # add an edge with between the mapped buses, if it does not exist yet
            if (!has_edge(st_graph, code_for(st_graph, v1), code_for(st_graph, v2)))
                @assert(MetaGraphsNext.add_edge!(st_graph, v1, v2, b_id))
            end
        end

        # vertices 1:length(islands) are per construction corresponding to the islands
        st = steiner_tree(st_graph, [code_for(st_graph, i) for i in (1 : length(islands))])
    else
        return 0
    end

    # approximation factor is 2-2/t where t is the number of terminals
    return ceil(ne(st)/(2-2/length(islands)))
end

"(deprecated?) check connectedness of the underlying graph as induced by the current master switching decisions given by `disabled_branches`"
function _underlying_graph_is_connected(pm, master_var_vals, graph::AbstractGraph, disabled_branches::Vector{Vector{Int64}}, enabled_branches::Vector{Vector{Int64}}, induced_isolated_nodes::Vector{Int64}, true_branch_endpoints::Set{Vector{Int64}}; nw::Int64=nw_id_default)::Bool
    # TODO: remember removed edges and re-add them after checking connectedness
    g = deepcopy(graph)

    # add branches that are to be enabled
    for (b_id, fbus, tbus) in enabled_branches
        if ! haskey(g, fbus)
            #println("Adding vertex $fbus.")
            @assert(MetaGraphsNext.add_vertex!(g, fbus, fbus))
            @assert(isinteger(code_for(g, fbus)))
        end
        if ! haskey(g, tbus)
            #println("Adding vertex $tbus.")
            @assert(MetaGraphsNext.add_vertex!(g, tbus, tbus))
            @assert(isinteger(code_for(g, tbus)))
        end
        if (!has_edge(g, code_for(g, fbus), code_for(g, tbus)) && !has_edge(g, code_for(g, tbus), code_for(g, fbus)) && !haskey(g, fbus, tbus))
            @assert(add_edge!(g, fbus, tbus, b_id), "Could not add edge $(fbus) -> $(tbus).")
        end
    end

    # remove branches that are to be disabled
    for (b_id, fbus, tbus) in disabled_branches
        if haskey(g, fbus, tbus)
            @assert(rem_edge!(g, code_for(g, fbus), code_for(g, tbus)), "Could not remove edge $(fbus) -> $(tbus) ($b_id).")
        end
    end

    # remove nodes that were previously only connected to branches that are now turned off
    for bus in induced_isolated_nodes
        @assert(rem_vertex!(g, code_for(g, bus)))
    end

    # check all connected components except for the main network for power demand. If there is none, this is just an artifact and all demand is still covered
    components = connected_components(g)
    islands = []
    if length(components) > 1
        for connected_component in components 
            if _network_powered(pm.data["nw"]["$nw"]["bus"], pm.data["nw"]["$nw"]["load"], pm.data["nw"]["$nw"]["gen"], bus_list=label_for.([g], connected_component))
                push!(islands, connected_component)
            end
        end
    end

    return length(islands) <= 1
end

#--------------------------------------------------------
"""
`min_edge_cut(G)` returns a minimum size set of edges whose removal
disconnects `G`. The graph must have at least two vertices.

`min_edge_cut(G,s,t)` is a minimum size set of edges whose removal
separates vertices `s` and `t`.
"""
function min_edge_cut(G::UG{T}, s::T, t::T, costs::Dict{Tuple{T,T},Int64}, flag::Bool = true)::Set{Tuple{T,T}} where {T}
    n = NV(G)
    n > 1 || error("Graph must have at least two vertices")

    if flag
        has(G, s) || error("$s is not a vertex of this graph")
        has(G, t) || error("$t is not a vertex of this graph")
        s != t || error("vertices must not be the same")
    else
        if cache_check(G, :min_edge_cut)
            return cache_recall(G, :min_edge_cut)
        end

        if !SimpleGraphAlgorithms.is_connected(G)
            X = Set{Tuple{T,T}}()
            cache_save(G, :min_edge_cut, X)
            return X
        end
    end

    VV = vlist(G)
    EE = elist(G)

    MOD = SimpleGraphAlgorithms.Model(SimpleGraphAlgorithms.get_solver())

    SimpleGraphAlgorithms.@variable(MOD, a[VV], Bin)  # in part 1
    SimpleGraphAlgorithms.@variable(MOD, b[VV], Bin)  # in part 2
    SimpleGraphAlgorithms.@variable(MOD, c[EE], Bin)  # span between the parts

    for v in VV
      SimpleGraphAlgorithms.@constraint(MOD, a[v] + b[v] == 1)
    end

    SimpleGraphAlgorithms.@constraint(MOD, sum(a[v] for v in VV) >= 1)
    SimpleGraphAlgorithms.@constraint(MOD, sum(b[v] for v in VV) >= 1)

    if flag
      SimpleGraphAlgorithms.@constraint(MOD, a[s] == 1)
      SimpleGraphAlgorithms.@constraint(MOD, b[t] == 1)
    end


    for e in EE
      u, v = e
      SimpleGraphAlgorithms.@constraint(MOD, a[u] + b[v] - 1 <= c[e])
      SimpleGraphAlgorithms.@constraint(MOD, a[v] + b[u] - 1 <= c[e])
    end

    SimpleGraphAlgorithms.@objective(MOD, Min, sum(costs[e]*c[e] for e in EE))

    SimpleGraphAlgorithms.optimize!(MOD)
    status = Int(SimpleGraphAlgorithms.termination_status(MOD))

    C = SimpleGraphAlgorithms.value.(c)

    X = Set(e for e in EE if C[e] > 0.5)

    if !flag
      SimpleGraphAlgorithms.cache_save(G, :min_edge_cut, X)
    end
    return X
end

function min_edge_cut(G::UG, costs::Dict{Tuple{T, T},Int64})::Set where {T}
    n = NV(G)
    n > 1 || error("Graph must have at least two vertices")
    s = first(G.V)
    return min_edge_cut(G, s, s, costs, false)
end