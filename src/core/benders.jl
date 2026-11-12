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

import JuMP

##### prio list handling #####

###### init functions #####
"sets the priority of any subproblem to 0.0"
function _prio_init_zero(c::Int64)::Float64
    return 0.0
end

##### pre prio functions #####
"
    _prio_pre_weighted_sum!(subproblem_prio_list::Vector{Dict{Symbol, Any}}, counters::Dict{Symbol, Real})

For a given list of score names and their weights, calculate the weighted sum and write it into entry[:priority] for each subproblem.

## Arguments
- `subproblem_prio_list`: list of subproblems
- `counters`: dictionary of counters (not yet used)
"
function _prio_pre_weighted_sum!(subproblem_prio_list::Vector{Dict{Symbol, Any}}, counters::Dict{Symbol, Real})
    if PARAMS[:BENDERS][:PARTIAL][:SORTING_MODE] == :SCORE
        score_names = [score_name for (score_name, score_weight) in PARAMS[:BENDERS][:PARTIAL][:SCORE] if score_weight > 0.0]
        if :total_share ∈ score_names
            n_cuts_total = max(sum(e[:n_cuts_used] for e in subproblem_prio_list), 1)
        end

        #@info("Calculating scores for $(length(subproblem_prio_list)) subproblems")
        for entry in subproblem_prio_list
            # calc reliability
            if :reliability ∈ score_names
                entry[:reliability] = entry[:n_calls] == 0 ? 0 : entry[:n_cuts_used] / entry[:n_calls]
            end

            # calc total_share
            if :total_share ∈ score_names
                entry[:total_share] = entry[:n_cuts_used] / n_cuts_total
            end

            if :n_since_call ∈ score_names
                entry[:n_since_call] = entry[:round_called] == -1 ? 0 : callback_stats[:n_callbacks] - entry[:round_called]
            end

            if :violation_per_cut ∈ score_names
                entry[:violation_per_cut] = entry[:n_cuts_used] == 0 ? 0 : entry[:sum_violations] / entry[:n_calls]
            end

            # calc weighted sum of scores to get priority
            entry[:priority] = sum(isnan(entry[score_name]) ? 0.0 : entry[score_name] * PARAMS[:BENDERS][:PARTIAL][:SCORE][score_name] for score_name ∈ score_names)

            #@info("Score components for subproblem $(entry[:contingency_id]): $(join(["$(score_name): $(entry[score_name]) * $(weight) = $(entry[score_name] * weight)" for (score_name, weight) in zip(score_names, weights)], ", "))")
        end

        # scale priorities to [0, 1]
        max_prio = maximum(entry[:priority] for entry in subproblem_prio_list)
        min_prio = minimum(entry[:priority] for entry in subproblem_prio_list)
        if (max_prio - min_prio) == 0
            return
        else
            for entry in subproblem_prio_list
                entry[:priority] = (entry[:priority] - min_prio) / (max_prio - min_prio)
            end
        end
    elseif PARAMS[:BENDERS][:PARTIAL][:SORTING_MODE] == :RND
        # randomize priorities from 0 to 1
        for entry in subproblem_prio_list
            entry[:priority] = rand()
        end
    else
        @error("Unknown score method: $(PARAMS[:BENDERS][:PARTIAL][:SCORE_METHOD])")
    end
end

##### post prio functions #####
function _prio_post_set_call(subproblem_prio_list_item::Dict{Symbol, Any})::Float64
    subproblem_prio_list_item[:round_called] = callback_stats[:n_callbacks]
end

"returns `1` if at least one constraint was found by the subproblem, else `0`"
function _prio_post_set_on_hit(old_prio::Float64, con::Union{JuMP.ScalarConstraint, Nothing})::Float64
    return (isnothing(con) ? 0 : 1)
end

"sets priority to `1` if at least one constraint was found by the subproblem, else depreciates it by `depreciation_rate`"
function _prio_post_depreciation_set(old_prio::Float64, con::Union{JuMP.ScalarConstraint, Nothing}, depreciation_rate=0.9)::Float64
    return (isnothing(con) ? old_prio*depreciation_rate : 1)
end

"adds `1` to the priority if at least one constraint was found by the subproblem, else depreciates it by `depreciation_rate`"
function _prio_post_depreciation_add(old_prio::Float64, con::Union{JuMP.ScalarConstraint, Nothing}, depreciation_rate=0.9)::Float64
    return (isnothing(con) ? old_prio*depreciation_rate : old_prio+1)
end

"appreciate priority if a cut is found by `appreciation_rate` (can be <1) and depreciate by `depreciation_rate` if no cut was found."
function _prio_post_reliability_mult(old_prio::Float64, con::Union{JuMP.ScalarConstraint, Nothing}, depreciation_rate=0.999, appreciation_rate=1.001)::Float64
    return (isnothing(con) ? old_prio*depreciation_rate : old_prio*appreciation_rate)
end

"adds 0.5 if a cut was found, subtracts 0.5 if no cut was found. capped at -2/+2."
function _prio_post_reliability_add(old_prio::Float64, con::Union{JuMP.ScalarConstraint, Nothing})::Float64
    return (isnothing(con) ? max(old_prio - 0.5, -2) : min(old_prio + 0.5, 2))
end


##### setters #####
"branches that fail in the current subproblem. If `filter=true`, only return branches that are not already turned off."
function get_outage_branches(pm_sub::AbstractPowerModel, contingency::Contingency, master_var_vals::Dict{Symbol, AbstractVector}; filter::Bool=false, nw::Int64=nw_id_default)::Vector{Int64}
    outage_branches = Vector{Int64}(copy(contingency[:branch]))

    # filter out branches that are already off (if they are switchable, check z_branch, else check status)
    if filter
        for branch_id in outage_branches
            branch = ref(pm_sub, 1, :branch, branch_id)
            # if branch is inactive (if possible, by choice, else by data)
            if (is_branch_switchable(branch) ? isapprox(master_var_vals[:z_branch][branch_id], 0.0, atol=SOLVER_PRECISION) : branch["status"] == 0)
                filter!(b_id -> b_id ≠ branch_id, outage_branches)
            end
        end
    end

    return outage_branches
end

"based on `branches`, turn off the corresponding branch in `master_var_vals`"
function switch_on_outage_branches!(master_var_vals::Dict{Symbol, AbstractVector}, branches::Vector{Int64})
    if ! PARAMS[:LIMITS][:BRANCH][:BRANCH_SWITCHING]
        return
    end
    for branch_id in branches
        if branch_id ∈ axes(master_var_vals[:z_branch], 1)
            master_var_vals[:z_branch][branch_id] = 1.0
        end
    end
end

"based on branches, turn off (if not already off in `master_var_vals` and return master_var_vals and all switching actions)"
function switch_off_outage_branches!(master_var_vals::Dict{Symbol, AbstractVector}, branches::Vector{Int64})
    if ! PARAMS[:LIMITS][:BRANCH][:BRANCH_SWITCHING]
        return
    end
    for branch_id in branches
        if branch_id ∈ axes(master_var_vals[:z_branch], 1)
            master_var_vals[:z_branch][branch_id] = 0.0
        end
    end
end


##### getters #####
"helper function to categorize subproblem constraints for easier access in the callback function"
function get_sub_cons_for_duals(pm_sub::AbstractPowerModel)::Dict{Symbol, Dict{Int, JuMP.ConstraintRef}}
    # create cons dict
    sub_cons_for_duals = Dict(
        :power_balance_lb => Dict(),
        :power_balance_ub => Dict(),
        :rating_from_slack_lb => Dict(),
        :rating_from_slack_ub => Dict(),
        :rating_from_slack_on_off_lb => Dict(),
        :rating_from_slack_on_off_ub => Dict(),
        :ohm_fr_lb => Dict(),
        :ohm_fr_ub => Dict(),
        :ohm_fr_on_off_lb => Dict(),
        :ohm_fr_on_off_ub => Dict(),
        :vad_lb => Dict(),
        :vad_ub => Dict(),
        :vad_on_off_lb => Dict(),
        :vad_on_off_ub => Dict(),
        :breaker_branch_limit => Dict(),
        #:theta_ref => Dict(),
        :pg_on_off_lb => Dict(),
        :pg_on_off_ub => Dict()
    )

    # categorize subproblem constraints
    conss = JuMP.all_constraints(pm_sub.model, include_variable_in_set_constraints = false)

    n_other = 0
    for con in conss
        added = false
        for cons_symbol in keys(sub_cons_for_duals)

            # check if constraint is of matching type, then add to dict with the corresponding comp(onent)_id
            if occursin(string(cons_symbol), JuMP.name(con))
                comp_id = parse(Int64, split(JuMP.name(con), "_")[end])
                push!(sub_cons_for_duals[cons_symbol], comp_id => con)
                added = true
                break
            end

        end
        if ! added
            n_other += 1
        end
    end
    @assert(n_other <= 1, "We did not use all dual variables! Remaining: $n_other") # 1 for theta_ref
    return sub_cons_for_duals
end

"helper function to categorize middle problem constraints (gen master problem) for easier access in the callback function"
function get_sub_cons_for_duals_gen_model(pm_sub::AbstractPowerModel; with_pf::Bool=false)::Dict{Symbol, Dict{Int, JuMP.ConstraintRef}}
    # create cons dict
    sub_cons_for_duals = Dict(
        #:z_branch_fix => Dict(),
        :power_balance_ub => Dict(),
        :power_balance_lb => Dict(),
        :pg_delta_minus_on_off_linking => Dict(),
        :pg_delta_plus_on_off_linking => Dict(),
        :pg_delta_minus => Dict(),
        :pg_delta_plus => Dict(),
        :pg_delta_plus_lb => Dict(),
        :pg_delta_plus_ub => Dict(),
        :pg_delta_minus_lb => Dict(),
        :pg_delta_minus_ub => Dict(),
        :pg_delta_linking => Dict(),
        :pg_gsk_linking_minus => Dict(),
        :pg_gsk_linking_plus => Dict(),
        :pg_on_off_lb => Dict(),
        :pg_on_off_ub => Dict(),
        :benders_cut => Dict() 
    )

    if with_pf
        # add constraints for power flow
        push!(sub_cons_for_duals,
            :ohm_fr_on_off_ub => Dict(),
            :ohm_fr_on_off_lb => Dict(),
            :ohm_fr_ub => Dict(),
            :ohm_fr_lb => Dict(),
            :vad_on_off_ub => Dict(),
            :vad_on_off_lb => Dict(),
            :vad_ub => Dict(),
            :vad_lb => Dict(),
            :rating_from_slack_ub => Dict(),
            :rating_from_slack_lb => Dict(),
            :rating_from_slack_on_off_ub => Dict(),
            :rating_from_slack_on_off_lb => Dict(),
        )
    end

    # categorize subproblem constraints
    conss = JuMP.all_constraints(pm_sub.model, include_variable_in_set_constraints = false)

    n_other = 0
    for con in conss
        added = false
        for cons_symbol in keys(sub_cons_for_duals)

            # check if constraint is of matching type, then add to dict with the corresponding comp(onent)_id
            if occursin(string(cons_symbol), JuMP.name(con))
                comp_id = parse(Int64, split(JuMP.name(con), "_")[end])
                push!(sub_cons_for_duals[cons_symbol], comp_id => con)
                added = true
                break
            end

        end
        if ! added && ! occursin("theta_ref", JuMP.name(con)) # theta_ref never has a nonzero dual variable, so we can ignore it
            @info("Not added: $(JuMP.name(con))")
            n_other += 1
        end
    end
    @assert(n_other == 0, "We did not use all dual variables! Remaining: $n_other")
    return sub_cons_for_duals
end

"get the basename and index of a variable as symbol, e.g. '1_pg[5]' => :pg and '1_z_branch[13]' => :z_branch"
function get_var_basename_idx(var::JuMP.VariableRef)::Tuple{Symbol, Int64}
    name = JuMP.name(var)

    name_first = findfirst('_', name)
    name_last  = findlast('[', name)
    idx_last   = findlast(']', name)

    return Symbol(name[name_first+1:name_last-1]), parse(Int64, name[name_last+1:idx_last-1])
end

"get the basename of a variable as symbol, e.g. '1_pg[5]' => :pg and '1_z_branch[13]' => :z_branch"
function get_var_basename(var::JuMP.VariableRef)::Symbol
    name = JuMP.name(var)

    first = findfirst('_', name)
    last  = findlast('[', name)

    return Symbol(name[first+1:last-1])
end

"get the index of a variable as symbol, e.g. '1_pg[5]' => 5 and '1_z_branch[13]' => 13"
function get_var_idx(var::JuMP.VariableRef)::Int64
    return parse(Int64, match(r"\[(\d+)\]", JuMP.name(var)).captures[1])
end

##### helper functions #####
"fix variables from the master in the subproblem"
function fix_master_variables!(pm::AbstractPowerModel, nw::Int64, master_var_vals::Dict{Symbol, AbstractVector}, int_sol::Bool; ignore_var_names=Vector{Symbol}())
    # fix variables and unset their binary status
    for (var_name, vals) in master_var_vals
        if var_name == :η || var_name == :η_rd || var_name == :z_branch_change || var_name == :tr_shift_change || var_name ∈ ignore_var_names
            continue
        end
        vars = var(pm, nw, var_name)

        for var_idx in axes(vars, 1)
            x = vars[var_idx]
            
            is_int = true
            if JuMP.is_binary(x)
                JuMP.unset_binary(x)
            elseif JuMP.is_integer(x)
                JuMP.unset_integer(x)
            else
                is_int = false
            end

            if is_int && int_sol
                @assert(isapprox(round(vals[var_idx], digits=0), vals[var_idx], atol =SOLVER_PRECISION*10), "Master variable $var_name is not integer: $(vals[var_idx])")
                JuMP.fix(x, round(vals[var_idx], digits=0), force=true)
            else
                vl = isapprox(vals[var_idx], 0.0, atol =SOLVER_PRECISION) ? 0 : vals[var_idx]             
                JuMP.fix(x, vl, force=true)
            end
        end
    end
end

"get the average line load in the current solution for a branch, i.e., the power flow on the branch divided by the adjusted rate_a"
function get_line_load_in_cur_sol(pm::AbstractPowerModel, nw::Int64, branch_id::Int64, cb_data::Gurobi.CallbackData)::Float64
    p = JuMP.callback_value(cb_data, var(pm, nw, :p, (branch_id, ref(pm, nw, :branch, branch_id)["f_bus"], ref(pm, nw, :branch, branch_id)["t_bus"])))
    return abs(p) / get_adjusted_rate_a(pm, ref(pm, nw, :branch, branch_id))
end
function get_line_load_in_cur_sol(pm::AbstractPowerModel, nw::Int64, branch_id::Int64, var_obj::Dict{Any, Any})::Float64
    p = JuMP.value(var_obj[(branch_id, ref(pm, nw, :branch, branch_id)["f_bus"], ref(pm, nw, :branch, branch_id)["t_bus"])])
    return abs(p) / get_adjusted_rate_a(pm, ref(pm, nw, :branch, branch_id))
end

"for later score calculation, add information of the current solution to the subproblem_prio_list"
function _add_data_to_subproblem_prio_list!(entry::Dict{Symbol, Any}, pm::AbstractPowerModel, lookup::Union{Gurobi.CallbackData, Dict{Any, Any}}; nw::Int=1)
    outage_branches = entry[:contingency][:branch] # we need the power flow, so we do not consider breaker branches
    
    if ! isempty(outage_branches) && PARAMS[:BENDERS][:BASE_PF_IN_MASTER]
        # add average line load in master solution amongst outage branches
        #entry[:load]    = sum(abs(JuMP.callback_value(cb_data, var(pm, nw, :p, (b, ref(pm, nw, :branch, b)["f_bus"], ref(pm, nw, :branch, b)["t_bus"])))) / get_adjusted_rate_a(pm, ref(pm, nw, :branch, b)) for b in outage_branches)/length(outage_branches)
        entry[:load] = sum(get_line_load_in_cur_sol(pm, nw, b, lookup) for b in outage_branches)/length(outage_branches)

        # add sw_heur_score (note: based on warm-start solution, not current master solution, as the calculation can be rather expensive)
        entry[:sw_score] = sum(get(ref(pm, nw, :branch, b), "sw_heur_score", 0.0) for b in outage_branches)/length(outage_branches)

        # TODO: add dual values to use for alpha-criterion?
        #entry[:alpha] = JuMP.dual(pm.model, ...)
    else
        entry[:load] = 0.0
        entry[:sw_score] = 0.0
        #entry[:alpha] = 0.0
    end
end

"print switching decisions made by the master problem. returns branches switched off in the problem."
function get_master_decisions(pm::AbstractPowerModel, cb_data; nw::Int64=nw_id_default)::Tuple{Vector{Int64}, Vector{Int64}}
    switched_off_by_master = Vector{Int64}()
    switched_on_by_master = Vector{Int64}()

    if haskey(var(pm, nw), :z_branch_change)
        for z_change_var in var(pm, nw, :z_branch_change)
            z_change_val = JuMP.callback_value(cb_data, z_change_var)
            if isapprox(z_change_val, 1, atol = SOLVER_PRECISION) # z_branch_change is 1, so a switch took place
                br_id       = get_var_idx(z_change_var)
                branch      = ref(pm, nw, :branch, br_id)
                z_val       = JuMP.callback_value(cb_data, var(pm, nw, :z_branch, br_id))
                
                # if there is no cost on z_branch_change, z_branch_change can be 1, but no switch took place
                if ! isapprox(z_val, branch["status"], atol =SOLVER_PRECISION) 
                    if branch["status"] == 1
                        push!(switched_off_by_master, br_id)
                        #@info("Switching OFF branch $(lpad(br_id, 4, ' ')) $(branch["name"])")
                    elseif branch["status"] == 0
                        push!(switched_on_by_master, br_id)
                        #@info("Switching  ON branch $(lpad(br_id, 4, ' ')) $(branch["name"])")
                    end
                end
            end
        end
    end
    return switched_off_by_master, switched_on_by_master
end

"check if the network induced by the `bus_list` has any loads that need to be satisfied"
function _network_powered(data::Dict; bus_list::Vector{Int64}=[])::Bool
    return _network_powered(data[:bus], data[:load], data[:gen], bus_list = bus_list)
end

"check if the network induced by the `bus_list` has any loads that need to be satisfied"
function _network_powered(bus_data::Dict, load_data::Dict, gen_data::Dict; bus_list::Vector{Int64}=[])::Bool
    if isempty(bus_list)
        bus_list = keys(bus_data)
    end
    total_pd = 0
    for (_, load) in load_data
        if load["load_bus"] ∈ bus_list && load["status"] == 1
            total_pd += load["pd"]::Float64
        end
    end
    total_pg_min, total_pg_max = 0.0, 0.0
    for (_, gen) in gen_data
        if gen["gen_bus"] ∈ bus_list && gen["status"] == 1
            total_pg_min += gen["switchable"] ? min(0, gen["pmin"]::Float64) : gen["pmin"]::Float64
            total_pg_max += gen["pmax"]::Float64
        end
    end
    
    return total_pd > 0 || total_pg_min > 0
end

"get the f_bus and t_bus of a branch"
function _get_branch_fbus_tbus(branch_id::Int64, pm::AbstractPowerModel; nw::Int64=nw_id_default)::Vector{Int64}
    branch = ref(pm, nw, :branch, branch_id)
    return [branch["f_bus"], branch["t_bus"]]
end

"get the true endpoints of branches, i.e., get the original branch, and then extend by the breaker branches at either end, if applicable"
function _get_branch_endpoints(branch_id::Int64, pm::AbstractPowerModel; nw::Int64=nw_id_default)::Vector{Int64}
    branch = ref(pm, nw, :branch, branch_id)
    endpoints = Vector{Int64}()
    # obtain an original branch
    if get(branch, "type", "") == "breaker"
        orig_b_id = branch["original_branch"]
    else
        orig_b_id = branch_id
    end
    orig_branch = ref(pm, nw, :branch, orig_b_id)

    # add orig branch endpoints, in case no or just partial breaker extension took place
    push!(endpoints, orig_branch["f_bus"])
    push!(endpoints, orig_branch["t_bus"])

    # obtain all _active_ breakers of the original branch
    for breaker_id ∈ get(orig_branch, "breaker_branches", [])
        breaker = ref(pm, nw, :branch, breaker_id)
        # use only active breaker
        if breaker["status"] == 1
            # get outward-facing endpoint and replace the endpoint of the original line with it 
            if breaker["f_bus"] ∈ endpoints
                filter!(e->e≠breaker["f_bus"], endpoints)
                push!(endpoints, breaker["t_bus"])
            else
                filter!(e->e≠breaker["t_bus"], endpoints)
                push!(endpoints, breaker["f_bus"])
            end
        end
    end
    @assert(length(endpoints) == 2, "Branch $branch_id has $(length(endpoints)) endpoints ($(endpoints)), but should have 2.")
    return endpoints
end

##### cut building #####
"returns a benders expression to be used in a Benders cut `expr <= η` following the current `dualss`and `outage_branches` values"
function build_benders_expression(pm_sub::AbstractPowerModel, constant::Float64, lhss::Dict{Symbol, Dict{Int, JuMP.GenericAffExpr}}, dualss::Dict{Symbol, Dict{Int, Float64}}, master_vars::Dict{Symbol, AbstractVector}, outage_branches::Vector{Int64})::JuMP.AffExpr
    expr = JuMP.@expression(pm_sub.model, constant + sum(sum(lhss[con_symbol][comp_id] * dual_val for (comp_id, dual_val) in duals if dual_val != 0.0) for (con_symbol, duals) in dualss))

    remove_var_from_expr!(expr, master_vars, outage_branches)

    return expr
end
"updates a benders expression to be used in a Benders cut `expr <= η` following the current `dualss`and `outage_branches` values"
function update_benders_expression!(expr::AffExpr, constant::Float64, lhss::Dict{Symbol, Dict{Int, JuMP.GenericAffExpr}}, dualss::Dict{Symbol, Dict{Int, Float64}})::JuMP.AffExpr
    terms = expr.terms::JuMP.OrderedCollections.OrderedDict{VariableRef, Float64}

    # reset
    for (var_ref, _) in terms
        terms[var_ref] = 0.0
    end

    # update cut
    for (con_symbol, duals) in dualss
        lhs = lhss[con_symbol]
        
        for (comp_id, dual_val) in duals
            if dual_val != 0.0
                lhs_val = lhs[comp_id]

                constant += lhs_val.constant::Float64 * dual_val
                
                for (var_ref, coeff) in lhs_val.terms::JuMP.OrderedCollections.OrderedDict{VariableRef, Float64}
                    terms[var_ref] += coeff * dual_val    
                end

            end
        end

    end

    expr.constant = constant

    return expr
end
function update_benders_expression!(expr::AffExpr, constant::Float64, lhss::Dict{Symbol, Dict{Int, JuMP.GenericAffExpr}}, dualss::Dict{Symbol, Dict{Int, Float64}}, master_vars::Dict{Symbol, AbstractVector}, outage_branches::Vector{Int64})::JuMP.AffExpr
    update_benders_expression!(expr, constant, lhss, dualss)

    remove_var_from_expr!(expr, master_vars, outage_branches)

    return expr
end

"reset all variable coefficients in the cut expression to 0"
function reset_cut_expr!(cut_expr::AffExpr, master_vars::Dict{Symbol, AbstractVector})
    for (_, var_vector) in master_vars
        for var in var_vector
            cut_expr.terms[var] = 0.0
        end
    end
end 
function reset_cut_exprs!(cut_exprs::Vector{AffExpr}, master_vars::Dict{Symbol, AbstractVector})
    for cut_expr in cut_exprs
        reset_cut_expr!(cut_expr, master_vars)
    end
end

"getter function for the dictionary of variable refs to their current solution value"
function get_master_var_vals_vars(master_vars::Dict{Symbol, AbstractVector}, cb_data::Gurobi.CallbackData, n_threads::Int64)::Vector{Dict{JuMP.VariableRef, Float64}}
    master_var_vals_var = Dict{VariableRef, Float64}()
    for (_, vars) in master_vars
        for var in vars
            master_var_vals_var[var] = JuMP.callback_value.(cb_data, var)
        end
    end
    return [copy(master_var_vals_var) for _ in 1:n_threads]::Vector{Dict{VariableRef, Float64}}
end
function get_master_var_vals_vars(master_vars::Dict{Symbol, AbstractVector}, master_var_vals::Dict{Symbol, AbstractVector}, n_threads::Int64)::Vector{Dict{VariableRef, Float64}}
    master_var_vals_var = Dict{VariableRef, Float64}()
    for (_, vars) in master_vars
        for var in vars
            var_basename = get_var_basename(var)
            if var_basename == :eta
                var_basename = :η
            end
            master_var_vals_var[var] = master_var_vals[var_basename][get_var_idx(var)]
        end
    end
    return [copy(master_var_vals_var) for _ in 1:n_threads]::Vector{Dict{VariableRef, Float64}}
end