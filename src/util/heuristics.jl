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

import StatsBase: sample, countmap

"alternative function to compute PTDF matrix"
function compute_PTDF(data)#, br_ids::Vector{Int64} = collect(keys(data["branch"])))
    bus_count = length(data["bus"])
    branch_count = length(data["branch"])
    
    # Identify the reference bus
    slack_bus = parse(Int64, findfirst(bus -> bus["bus_type"] == 3, data["bus"])[1])
    
    # Initialize B' matrix
    B_prime = zeros(bus_count, bus_count)
    
    # Populate B' matrix based on branch data
    for (b_id, branch) in data["branch"]
        if branch["status"] == 0
            continue # skip branches that are not in service
        end
        from_bus = branch["f_bus"]
        to_bus = branch["t_bus"]
        x = branch["br_x"]
        
        B_prime[from_bus, from_bus] += 1/x
        B_prime[to_bus, to_bus] += 1/x
        B_prime[from_bus, to_bus] -= 1/x
        B_prime[to_bus, from_bus] -= 1/x
    end
    
    # Exclude the reference bus from B' matrix
    B_prime = B_prime[setdiff(1:bus_count, slack_bus), setdiff(1:bus_count, slack_bus)]
    
    # Invert B' matrix
    B_prime_inv = inv(B_prime)
    
    # Initialize PTDF matrix
    PTDF = zeros(branch_count, bus_count)
    
    for (b_id, branch) in data["branch"]
        from_bus = branch["f_bus"]
        to_bus = branch["t_bus"]
        x = branch["br_x"]
        
        e_f = zeros(bus_count - 1)
        e_t = zeros(bus_count - 1)
        
        if from_bus != slack_bus
            e_f[from_bus > slack_bus ? from_bus - 1 : from_bus] = 1
        end
        if to_bus != slack_bus
            e_t[to_bus > slack_bus ? to_bus - 1 : to_bus] = 1
        end
        
        # TODO: this only has to be calculated for lines for which we want the LODF later on
        PTDF[parse(Int64, b_id), setdiff(1:bus_count, slack_bus)] = (e_f - e_t)' * B_prime_inv / x
    end
    
    return PTDF
end

"computes one entry of the LODF matrix"
function compute_LODF(PTDF, branch_data_of_i::Dict, i::Int, j::Int)
    # Calculate LODF based on PTDF
    from_bus_i = branch_data_of_i["f_bus"]
    to_bus_i = branch_data_of_i["t_bus"]
    
    PTDF_ij = PTDF[i, from_bus_i] - PTDF[i, to_bus_i]
    LODF = PTDF[j, from_bus_i] - PTDF[j, to_bus_i]
    LODF /= (1 - PTDF_ij)
    return LODF
end
"compute full LODF matrix for all branch pairs in the network"
function compute_LODF(PTDF, data)
    branch_count = length(data["branch"])
    @assert(branch_count == length(data["branch"]), "Maximum branch index does not match the number of branches")

    LODF = zeros(branch_count, branch_count)
    for i in 1:branch_count
        for j in 1:branch_count
            if i != j
                LODF[i,j] = compute_LODF(PTDF, data["branch"]["$i"], i, j)
            end
        end
    end
    
    return LODF
end

"get all branches that are at most n_hops (buses) away. Exclude the branch itself"
function get_neighborhood_branches(nw_data, branch, n_hops)::Vector{Int64}
    neighbors = Set{String}()
    cur_hop_set = Set([branch["f_bus"], branch["t_bus"]])
    for _ in 1:n_hops
        next_hop_set = Set{Int64}()
        for bus_id in cur_hop_set
            branch_ids = ["$b_id" for b_id in get_connected_branches(nw_data, bus_id) if "$b_id" ∉ neighbors]
            #@info("Branches connected to bus $bus_id: $branch_ids")
            union!(next_hop_set, Set([branch["f_bus"] for branch in [nw_data["branch"][b_id] for b_id in branch_ids] if branch["f_bus"] ≠ bus_id]), Set([branch["t_bus"] for branch in [nw_data["branch"][b_id] for b_id in branch_ids] if branch["t_bus"] ≠ bus_id]))
            union!(neighbors, Set(branch_ids))
        end
        cur_hop_set = next_hop_set
    end
    return collect(parse.(Int64, neighbors))
end


### Score Calculation ###
"Assign random switching relevance scores to each switchable branch in the network."
function calc_br_sw_heur_scores_rnd(nw_data; kwargs...)::Dict{Int64, Float64}
    br_sw_heur_scores = Dict{Int64, Float64}()

    for (_, branch) in nw_data["branch"]
        if branch["switchable"]
            br_sw_heur_scores[branch["index"]] = rand()
        end
    end
    return br_sw_heur_scores
end

"Assign the same switching relevance scores to each switchable branch in the network."
function calc_br_sw_heur_scores_full(nw_data; kwargs...)::Dict{Int64, Float64}
    br_sw_heur_scores = Dict{Int64, Float64}()

    for (_, branch) in nw_data["branch"]
        if branch["switchable"]
            br_sw_heur_scores[branch["index"]] = 1.0
        end
    end
    return br_sw_heur_scores
end

"collect pfs. If there are power flows for contingencies and `aggregate_over_contingencies`, aggregate them using max abs function"
function collect_pfs(nw_data, nw_sol; contingency_id::Int64=-1)::Dict{Int64, Float64}
    if contingency_id == -1
        pfs = Dict{Int64, Float64}([branch["index"] => nw_sol["branch"][string(branch["index"])]["pf"] for (_, branch) in nw_data["branch"]])
    else
        pfs = Dict{Int64, Float64}([branch["index"] => nw_sol["branch"]["($b_id, $contingency_id)"]["pf"] for (_, branch) in nw_data["branch"]])
    end
    return pfs
end

"
    calc_br_sw_heur_scores_sens(nw_data; nw_sol::Dict{String, Any}=Dict(), LOAD_LIMIT::Float64=0.6)::Dict{Int64, Float64}

Assign switching relevance scores to each switchable branch in the network based on the sensitivity heuristic (Liu et al., 2012). 

## Input
- `nw_data` -- the network data dictionary
- `nw_sol` -- the network solution dictionary
- `LOAD_LIMIT` -- the line loading at which a branch is considered congested
"
function calc_br_sw_heur_scores_sens(nw_data; nw_sol::Dict{String, Any}=Dict(), LOAD_LIMIT=0.9, N_MIN_CONGESTED_BRANCHES::Int64=PARAMS[:WARMSTARTING][:RSCOTS][:HEUR][:SENS][:N_MIN_CONGESTED_BRANCHES])::Dict{Int64, Float64}
    @assert(! haskey(nw_data["branch"]["1"], "pf_new"), "Data was overwritten by the solution. Provide the original data.")
    @assert(haskey(nw_sol, "branch"), "Solution empty or does not contain branch power flows")
    @assert(LOAD_LIMIT > 0.0, "Load limit must be greater than 0")
    #@assert(! nw_data["BUSBARS"], "BUSBARS not yet supported.") # TODO: to support them, simply deactivate switchability, do `merge_buses()`, all low-impedance branches will have been merged, run score, map back to original branches

    # calc PTDF and LODF 
    PTDF = compute_PTDF(nw_data)
    @assert(! any(isnan.(PTDF)), "PTDF matrix contains $(sum(isnan.(PTDF))) NaN values")
        
    # collect power flow and branch indices of overloaded branches
    pfs = collect_pfs(nw_data, nw_sol)
    congested_branches_idxs = [branch["index"] for (_, branch) ∈ nw_data["branch"] if is_branch_monitored(branch) && (pfs[branch["index"]] > LOAD_LIMIT * get_adjusted_rate_a(nw_data, branch))]
    @info("Number of (nearly) overloaded branches with load limit $LOAD_LIMIT: $(length(congested_branches_idxs))")

    # if too few branches are congested, set congested_branches_idxs to the N_MIN_CONGESTED_BRANCHES branch ids
    if length(congested_branches_idxs) < N_MIN_CONGESTED_BRANCHES
        @warn("No branches are congested. Using the $N_MIN_CONGESTED_BRANCHES most loaded branches instead.")
        sorted_branches_idxs = sort(collect(keys(pfs)), by = x -> abs(pfs[x]), rev = true)
        congested_branches_idxs = [idx for idx in sorted_branches_idxs[1:min(N_MIN_CONGESTED_BRANCHES, length(sorted_branches_idxs))]]
    end
    @info("Number of (nearly) overloaded branches after check: $(length(congested_branches_idxs))")

    # calculate sensitivities
    br_sw_heur_scores = Dict{Int64, Float64}()
    for (_, branch) in nw_data["branch"]
        if branch["switchable"]
            k  = branch["index"]
            s  = minimum([(abs((pfs[k] * compute_LODF(PTDF, nw_data["branch"][string(b)], b, k)) + pfs[b]) - abs(pfs[b])) / abs(pfs[b]) for b ∈ congested_branches_idxs if pfs[b] != 0], init = 0.0)
            br_sw_heur_scores[branch["index"]]  = -s # turn around the sign, since we want higher scores for branches with large, negative sensitivity
        end
    end

    # print table with branch names and scores
    # @info("Branch sensitivities:")
    # for (branch_idx, score) in sort(collect(br_sw_heur_scores), by = x -> x[2])
    #    branch = nw_data["branch"][string(branch_idx)]
    #    @printf("  %-36s (%4d): %-8.3f\n", branch["branch_name_intern"], branch_idx, score)
    # end

    # show distribution of scores, binned into 0.1 intervals, sort and show
    score_distribution = Dict{Float64, Int64}(reverse(sort(collect(countmap(round.(values(br_sw_heur_scores), digits=1))))))
    @info("Distribution of scores: $(join(["$(k): $(v)" for (k, v) in sort(collect(score_distribution), by = x -> x[1])], ", "))")

    return br_sw_heur_scores
end

"
    calc_br_sw_heur_scores_congprox(nw_data; nw_sol::Dict{String, Any}=Dict(), LOAD_LIMIT::Float64=0.6, N_HOPS::Int64=2)::Dict{Int64, Float64}

Assign switching relevance scores to each branch in the network based on the congestion proximity heuristic.
Congested branches get 1.0, branches that are at most `N_HOPS` away get 1/hop_dist, all others get 0.

## Input
- `nw_data` -- the network data dictionary
- `nw_sol` -- the network solution dictionary
- `LOAD_LIMIT` -- the line loading at which a branch is considered congested
- `N_HOPS` -- the number of hops until a branch is assigned a non-zero score
"
function calc_br_sw_heur_scores_congprox(nw_data::Dict{String, Any}; 
                                         nw_sol::Dict{String, Any}=Dict(), 
                                         LOAD_LIMIT::Float64=0.8, 
                                         N_HOPS::Int64=2,
                                         kwargs...)::Dict{Int64, Float64}
    @assert(! haskey(nw_data["branch"]["1"], "pf_new"), "Data was overwritten by the solution. Provide the original data.")
    @assert(haskey(nw_sol, "branch"), "Solution empty or does not contain branch power flows")
    @assert(LOAD_LIMIT > 0.0, "Load limit must be greater than 0")
    
    br_sw_heur_scores = Dict{Int64, Float64}()
    
    # set congested branches to switchable
    congested_branches::Vector{Dict} = [branch for (b_id, branch) in nw_data["branch"] if is_branch_monitored(branch) && abs(nw_sol["branch"][b_id]["pf"]) > get_adjusted_rate_a(nw_data, branch) * LOAD_LIMIT]
    @info("Number of (nearly) overloaded branches: $(length(congested_branches))")

    # initialize scores for each switchable branch to 0.0
    for (_, branch) in nw_data["branch"]
        if branch["switchable"]
            br_sw_heur_scores[branch["index"]] = 0.0
        end
    end

    # for each congested branch and all branches that are at most `N_HOPS` away, set score to 1/(1+hopdist)
    for branch in congested_branches
        for hop_dist in 1:N_HOPS
            br_idxs = get_neighborhood_branches(nw_data, branch, hop_dist)
            for br_idx in br_idxs
                if is_branch_switchable(nw_data["branch"][string(br_idx)])
                    br_sw_heur_scores[br_idx] = max(br_sw_heur_scores[br_idx], 1/(1+hop_dist))
                end
            end
        end
        if is_branch_switchable(branch)
            br_sw_heur_scores[branch["index"]] = 1.0
        end
    end

    # show distribution of scores 
    @info("Distribution of scores: $(Dict{Float64, Int64}(reverse(sort(collect(countmap(values(br_sw_heur_scores)))))))")
    return br_sw_heur_scores
end


### Score Calculation Wrappers ###
"
    calc_br_sw_heur_scores(nw_data; heuristic::Symbol=:CONGPROX, kwargs...)::Dict{Int64, Float64}

Calculate branch switching relevance scores based on the selected heuristic.

## Input
- `nw_data` -- the data dictionary for a single network
- `heuristic` -- the heuristic to use. Possible values are:
    - `rand`: randomly assign a score to each branch
    - `full`: each branch gets the same score of
    - `sens`: assign a score based on the sensitivity of the branch
    - `congprox`: assign a score based on the proximity to congested branches
- `nw` -- the network index to calculate the scores for
- `kwargs` -- additional keyword arguments to pass to the heuristic calculation function. Should usually contain nw_sol, since most reasonable heuristics require it

## Output
- a dictionary with branch indices as keys and the assigned scores as values. Scores are in the range [0, 1].
"
function calc_br_sw_heur_scores(nw_data::Dict{String, Any}; heuristic::Symbol=:CONGPROX, kwargs...)::Dict{Int64, Float64}
    if heuristic == :RND
        br_sw_heur_scores = calc_br_sw_heur_scores_rnd(nw_data; kwargs...)
    elseif heuristic == :FULL
        br_sw_heur_scores = calc_br_sw_heur_scores_full(nw_data; kwargs...)
    elseif heuristic == :SENS
        br_sw_heur_scores = calc_br_sw_heur_scores_sens(nw_data; kwargs...)
    elseif heuristic == :CONGPROX
        br_sw_heur_scores = calc_br_sw_heur_scores_congprox(nw_data; kwargs...)
    else
        @error("Unknown heuristic: $heuristic")
    end

    # check scores
    n_replaced = 0
    for (k, v) in br_sw_heur_scores
        if isnan(v)
            br_sw_heur_scores[k] = 0.0
            n_replaced += 1
        end
    end
    if n_replaced > 0
        @info("Replaced $n_replaced NaNs by 0 ($(round(n_replaced / length(br_sw_heur_scores) * 100, digits=2))%)")
        if n_replaced / length(br_sw_heur_scores) > 0.2
            @warn("NaN ratio too high (> 20%), running again with congprox heuristic.")
            br_sw_heur_scores = calc_br_sw_heur_scores_congprox(nw_data; kwargs..., PARAMS[:WARMSTARTING][:RSCOTS][:HEUR][:CONGPROX]...)
        end
    end

    return br_sw_heur_scores
end

"
    calc_and_set_br_sw_heur_scores!(nw_data; set_switchability::Bool=true, percentile::Float64=0.2, heuristic::Symbol=:RND, kwargs...)

Calculate and set branch switching relevance scores based on the selected heuristic.

## Input
- `nw_data` -- the data dictionary for a single network
- `set_switchability` -- whether to set the `switchable` flag for each branch based on the heuristic
- `percentile` -- the percentile of the branch switching priorities to use as the cutoff for setting the `switchable` flag
- `force_percentile` -- whether to enforce the number of switchable branches to be at most `percentile` (if `true`) or to be exactly `percentile` (if `false`) by sampling randomly from the branches with the same score
- `force_mode` -- the mode to use for enforcing the number of switchable branches. Possible values are:
    - `sample`: randomly sample branches with the same score until the number of switchable branches is at most `percentile`
    - `first`: remove all branches with the same score until the number of switchable branches is at most `percentile`
- `heuristic` -- the heuristic to use
- `nw` -- the network index to calculate the scores for
- `kwargs` -- additional keyword arguments to pass to the heuristic calculation function. Should usually contain nw_sol, since most reasonable heuristics require it
"
function calc_and_set_br_sw_heur_scores!(nw_data; 
                                         set_switchability::Bool=true, 
                                         params::Dict{Symbol, Any},
                                         heur_params::Dict{Symbol, Any},
                                         kwargs...)

    heuristic = params[:BR_SW_HEUR]
    branch_switching_scores::Dict{Int64, Float64} = calc_br_sw_heur_scores(nw_data; heuristic=heuristic, heur_params..., kwargs...)
    
    if set_switchability
        percentile = params[:N_SW_BR_REL]
        force_percentile = params[:N_SW_BR_FORCE]
        force_mode = params[:N_SW_BR_FORCE_MODE]
        cutoff_mode = params[:SWITCHABLE_BRANCHES_LIMIT_MODE]
        br_sw_heur_score_values::Vector{Float64} = sort(collect(values(branch_switching_scores)))

        if cutoff_mode == :abs_num
            br_sw_heur_score_cutoff = br_sw_heur_score_values[max(1, length(br_sw_heur_score_values) - params[:N_SW_BR_ABS] + 1)]
        elseif cutoff_mode == :abs_score
            br_sw_heur_score_cutoff = params[:SW_BR_ABS_SCORE_LIM]
        elseif cutoff_mode == :rel
            br_sw_heur_score_cutoff = br_sw_heur_score_values[round(Int, length(br_sw_heur_score_values) * percentile)]
        else
            @error("Unknown cutoff mode: $cutoff_mode. Use :abs or :rel.")
        end

        # if mode is :FULL, set score cutoff to -10000 to make sure all branches are switchable
        if heuristic == :FULL
            br_sw_heur_score_cutoff = -10000.0
            force_percentile = false
        end
        @info("Cutoff for switchable branches: $br_sw_heur_score_cutoff")
    
        sw_branches_after_score::Vector{Int64} = [k for (k, v) in branch_switching_scores if v ≥ br_sw_heur_score_cutoff]
        if isempty(sw_branches_after_score)
            # if the set is empty, the MIP might become an LP, not calling the callback, and thus making the remaining pipeline fail
            @warn("No switchable branches found with score above cutoff $br_sw_heur_score_cutoff. Setting the first branch as switchable.")
            sw_branches_after_score = [first(keys(branch_switching_scores))]
        end

        # filter out branches with the same score randomly
        if force_percentile && length(sw_branches_after_score) > round(Int, length(br_sw_heur_score_values) * percentile)
            # get all branches with the same score (highest score still above cutoff) 
            max_score_above_cutoff = minimum([v for (k, v) in branch_switching_scores if v ≥ br_sw_heur_score_cutoff])
            sw_branches_tied::Vector{Int64} = [k for (k, v) in branch_switching_scores if v == max_score_above_cutoff]

            # select random branches to remove
            if force_mode == :sample
                n_to_remove = length(sw_branches_after_score) - round(Int, length(br_sw_heur_score_values) * percentile)
                sw_branches_after_score = setdiff(sw_branches_after_score, sample(sw_branches_tied, n_to_remove, replace=false))
            # remove all branches with the same score (highest score still below cutoff)
            elseif force_mode == :first
                sw_branches_after_score = setdiff(sw_branches_after_score, sw_branches_tied)
            end
        end
    end

    # sort the branch switching scores and give each branch a priority (integer) 
    branch_switching_priorities = Dict{Int64, Int64}()
    for (i, branch_idx) in enumerate(sort(collect(keys(branch_switching_scores)), by = x -> branch_switching_scores[x], rev = true))
        branch_switching_priorities[branch_idx] = i
    end


    for (_, branch) in nw_data["branch"]
        branch["old_switchable"] = branch["switchable"]

        # if branch is not switchable, no score had to be assigned
        if ! branch["switchable"]
            continue
        end

        # set score (float)
        branch["sw_heur_score"] = branch_switching_scores[branch["index"]]

        # set priority (int)
        branch["sw_heur_priority"] = branch_switching_priorities[branch["index"]]

        # set switchability
        if set_switchability
            branch["switchable"] = branch["index"] ∈ sw_branches_after_score || get(branch, "heur_switched", false)
        end
    end
    @info("Switchable branches: $(length(filter(x -> x[2]["switchable"], nw_data["branch"])))")
end

### Helper Functions ###
"reset the switchability flags of all branches in the network to their original values"
function reset_switchability!(nw_data)
    for (_, branch) in nw_data["branch"]
        branch["heur_switchable"] = branch["switchable"]
        branch["switchable"] = branch["old_switchable"]
        branch["heur_switched"] = Bool(get(branch, "z_branch_change_start", 0.0))
        if branch["heur_switched"]
            @info("Branch $(branch["index"]) was switched in the heuristic. Leaving it as switchable.")
        end
    end
end

### Primal Heuristics (embedded in solving) ###
"given some indices of z_branch variables that should be set to 0, build a list of solutions where each solution has exactly one of these variables set to 0"
function build_benders_solutions(master_vars, master_var_vals, cb_data, z_branch_vars_not_one_idxs)::Tuple{Vector, Vector{Vector}}
    sols = Vector{Tuple{Vector{VariableRef}, Vector{Float64}}}()
    master_var_vals_modifieds = Vector{Dict{Symbol, AbstractVector}}() # store modified master variable values for each solution
    for idx in z_branch_vars_not_one_idxs
        master_var_vals_modified = copy(master_var_vals) # make a copy of the master variable values
        master_var_vals_modified[:z_branch][idx] = 0.0 # set the variable to 0

        variables, values = build_benders_solution(master_vars, 
                                                   master_var_vals[:η], 
                                                   cb_data, 
                                                   :η, 
                                                   master_var_vals = master_var_vals)
        push!(sols, (variables, values))
        push!(master_var_vals_modifieds, master_var_vals_modified)
    end
    return sols, master_var_vals_modifieds
end

"provide a list of all solutions where in each, exactly one of the z_branch variables that are < 1 in master_var_vals is set to 0"
function scots_rounding_heuristic_all_frac(pm, master_vars, master_var_vals, cb_data)
    sols = Vector{Tuple{Vector{VariableRef}, Vector{Float64}}}()
    # get all branch variables that are not set to 1 in the master solution
    z_branch_vars_not_one_idxs = Vector{Int64}()
    for var in master_vars[:z_branch]
        idx = get_var_idx(var)
        val = JuMP.callback_value(cb_data, var)
        if isapprox(val, 1.0, atol=SOLVER_PRECISION)
            continue # skip variables that are set to 1
        else
            push!(z_branch_vars_not_one_idxs, idx)
        end
    end
    @debug("HEUR: Using branch variables with indices $(z_branch_vars_not_one_idxs) for heuristic solution.")
    return build_benders_solutions(master_vars, master_var_vals, cb_data, z_branch_vars_not_one_idxs)
end

"use each branch's `sw_heur_score` for the list of branches to switch in the heuristic. Not recommended on fractional solutions"
function scots_rounding_heuristic_sw_score(pm, master_vars, master_var_vals, cb_data)
    sols = Vector{Tuple{VariableRef, Float64}}()
    # get all branch variables that are not set to 1 in the master solution
    z_branch_vars_not_one_idxs = Vector{Int64}()
    for var in master_vars[:z_branch]
        idx = get_var_idx(var)
        sw_heur_score = pm.data["nw"]["1"]["branch"][string(idx)]["sw_heur_score"]
        if sw_heur_score > 0.2
            val = JuMP.callback_value(cb_data, var)
            # skip variables already switched off
            if isapprox(val, 1.0, atol=SOLVER_PRECISION)
                push!(z_branch_vars_not_one_idxs, idx)
            end
        end
    end
    @debug("HEUR: Using branch variables with indices $(z_branch_vars_not_one_idxs) for heuristic solution.")
    return build_benders_solutions(master_vars, master_var_vals, cb_data, z_branch_vars_not_one_idxs)
end

"wrapper function for the rounding heuristic, which selects the appropriate heuristic based on the parameter settings"
function scots_rounding_heuristic(pm, master_vars, master_var_vals, cb_data)
    if PARAMS[:PRIMAL_HEURISTICS][:HEURISTIC] == :SW_SCORE
        return scots_rounding_heuristic_sw_score(pm, master_vars, master_var_vals, cb_data)
    elseif PARAMS[:PRIMAL_HEURISTICS][:HEURISTIC] == :ALL_FRAC
        return scots_rounding_heuristic_all_frac(pm, master_vars, master_var_vals, cb_data)
    else
        @error("Unknown heuristic type $(PARAMS[:PRIMAL_HEURISTICS][:HEURISTIC]) for Benders callback.")
        return [], master_var_vals
    end
end

"function to add branch switching priorities and set switchability based on the heuristic scores for the main solving run"
function enable_heuristic_solving!(data, application::String, timesteps::Int)
    if ! occursin("SCOTS", application) || timesteps > 1
        return
    end

    @warn("Solving $application heuristically using switching heuristic '$(PARAMS[:HEURISTIC_SOLVING][:RSCOTS][:BR_SW_HEUR])'.")
    @assert(haskey(data, "result"), "No result data available for heuristic solving. Try enabling WARMSTART_MODE to calculate powerflow for branch switching relevance score calculation.")
    calc_and_set_br_sw_heur_scores!(data["nw"]["1"], 
                                    nw_sol=data["result"]["solution"]["nw"]["1"], 
                                    params=PARAMS[:HEURISTIC_SOLVING][:RSCOTS], 
                                    heur_params=PARAMS[:HEURISTIC_SOLVING][:RSCOTS][:HEUR][PARAMS[:HEURISTIC_SOLVING][:RSCOTS][:BR_SW_HEUR]])

    if PARAMS[:HEURISTIC_SOLVING][:KEEP_WARMSTART_SWITCHABILITY]
        for (_, nw_data) in data["nw"]
            for (b_id, branch) in nw_data["branch"]
                if haskey(branch, "heur_switchable") && branch["heur_switchable"]
                    branch["switchable"] = true
                end
            end
        end
    end
end