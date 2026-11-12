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

"basic redispatch minimization"
function objective_min_redispatch(pm::AbstractPowerModel; kwargs...)
    #FIXME: tries to perform negative redispatch on negative gens (i.e., loads).
    return JuMP.@objective(pm.model, Min,
        sum(
            sum( var(pm, n,   :pg_delta_minus, i) for (i,gen) in nw_ref[:gen]) + 
            sum( var(pm, n,   :pg_delta_plus, i) for (i,gen) in nw_ref[:gen])
        for (n, nw_ref) in nws(pm))
    )
end

"minimize number of branch switching actions"
function objective_min_num_branch_switch(pm::AbstractPowerModel; kwargs...)
    return JuMP.@objective(pm.model, Min,
        sum(
            sum( z_branch_change for z_branch_change in var(pm, n,   :z_branch_change))
        for (n, nw_ref) in nws(pm))
    )
end

"minimize number of generator startup/shutdown actions"
function objective_min_num_gen_switch(pm::AbstractPowerModel; kwargs...)
    return JuMP.@objective(pm.model, Min,
        sum(
            sum( var(pm, n,   :z_gen_shutdown, i) for (i,br) in nw_ref[:gen]) +
            sum( var(pm, n,   :z_gen_startup, i) for (i,br) in nw_ref[:gen])
        for (n, nw_ref) in nws(pm))
    )
end

"minimize overall line load"
function objective_min_line_load(pm::AbstractPowerModel; kwargs...)
    p_abs = var(pm, 1)[:p_abs] = JuMP.@variable(pm.model,
        [(l,i,j) in ref(pm, 1, :arcs)], base_name="p_abs",
        start = comp_start_value(ref(pm, 1, :branch, l), "p_abs_start")
    )

    for (l, branch) in ref(pm, 1, :branch)
        i = branch["f_bus"]
        j = branch["t_bus"]
        JuMP.@constraint(pm.model, p_abs[(l,i,j)] >= var(pm, 1, :p, (l,i,j)))
        JuMP.@constraint(pm.model, p_abs[(l,i,j)] >= -var(pm, 1, :p, (l,i,j)))
    end

    return JuMP.@objective(pm.model, Min,
        sum(
            sum( p_abs[(i, branch["f_bus"], branch["t_bus"])] for (i, branch) in nw_ref[:branch])
        for (n, nw_ref) in nws(pm))
    )
end

"minimize number of transformer shifting actions"
function objective_min_num_tr_switch(pm::AbstractPowerModel; kwargs...)
    return JuMP.@objective(pm.model, Min,
        sum(
            sum( var(pm, n,   :tr_shift_change, i) for (i,br) in nw_ref[:branch] if is_branch_shiftable(pm, n, i))
        for (n, nw_ref) in nws(pm))
    )
end

"minimize deviation from initial solution"
function objective_min_change(pm::AbstractPowerModel; nw::Int=nw_id_default, kwargs...)
    println(" -> Obj min change")
    nw = 1
    p_dev = var(pm, nw)[:p_abs] = JuMP.@variable(pm.model,
        [(l,i,j) in ref(pm, nw, :arcs_from)], base_name="$(nw)_p_abs",
        start = comp_start_value(ref(pm, nw, :branch, l), "p_abs_start")
    )
    p = var(pm, nw, :p)

    for (l, branch) in ref(pm, nw, :branch)
        i = branch["f_bus"]
        j = branch["t_bus"]
        JuMP.@constraint(pm.model, p[(l,i,j)] - branch["pf_start"] <= p_dev[(l,i,j)])
        JuMP.@constraint(pm.model, branch["pf_start"] - p[(l,i,j)] <= p_dev[(l,i,j)])
    end

    

    return JuMP.@objective(pm.model, Min,
        sum(
            #sum( p_abs[(i, branch["f_bus"], branch["t_bus"])] - abs(branch["pf_start"]) for (i, branch) in nw_ref[:branch]) +
            sum(p_dev[(i, branch["f_bus"], branch["t_bus"])] for (i, branch) in nw_ref[:branch])
        for (n, nw_ref) in nws(pm))
    )
end

"minimize deviation of generator active power from initial solution (pg_start)"
function objective_min_gen_change(pm::AbstractPowerModel; nw::Int=nw_id_default, kwargs...)
    nw = 1
    pg_dev = var(pm, nw)[:pg_dev] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :gen)], base_name="$(nw)_pg_dev",
        lower_bound = 0.0
    )

    for (i, gen) in ref(pm, nw, :gen)
        pg_i = var(pm, nw, :pg, i)
        pg_start = gen["pg_start"]
        JuMP.@constraint(pm.model,  pg_i - pg_start <= pg_dev[i])
        JuMP.@constraint(pm.model,  pg_start - pg_i <= pg_dev[i])
    end

    return JuMP.@objective(pm.model, Min,
        sum(
            sum(var(pm, n, :pg_dev, i) for (i, gen) in nw_ref[:gen])
        for (n, nw_ref) in nws(pm))
    )
end

"for each gen, our costs are equal to `pg * cost_produce` + `(pg_delta_plus + pg_delta_minus) * cost_change"
function _expr_min_redispatch_cost(pm::AbstractPowerModel)
    return JuMP.@expression(pm.model,
        sum(
            sum(
                (var(pm, nw, :pg, i) - gen["pg_0"])                                * gen["cost_produce"] +
                (var(pm, nw, :pg_delta_minus, i) + var(pm, nw, :pg_delta_plus, i)) * gen["cost_change"]
            for (i, gen) in nw_ref[:gen] if((haskey(gen, "RDP-") && gen["RDP-"]>0)
                                          ||(haskey(gen, "RDP+") && gen["RDP+"]>0)))
        for (nw, nw_ref) in nws(pm)
        )
    )
end
function objective_min_redispatch_cost(pm::AbstractPowerModel)
    return JuMP.@objective(pm.model, Min,
        _expr_min_redispatch_cost(pm)
    )
end

function objective_min_switching(pm::AbstractPowerModel)
    return JuMP.@objective(pm.model, Min,
        sum(
            sum(
                (var(pm, nw, :z_branch_change, i)
            for (i, gen) in nw_ref[:branch]))
        for (nw, nw_ref) in nws(pm)
        )
    )
end

"objective function that minimizes the number of redispatched generators"
function objective_min_num_redispatch(pm::AbstractPowerModel; kwargs...)
    return JuMP.@objective(pm.model, Min,
        sum(
            sum( var(pm, n,   :pg_delta_on_off, i) for (i,gen) in nw_ref[:gen])
        for (n, nw_ref) in nws(pm))
    )
end

"objective function that minimizes the cost for starting and shutting down generators"
function objective_min_gen_start_stop_cost(pm::AbstractPowerModel; kwargs...)
    return JuMP.@objective(pm.model, Min,
        sum(
            sum( gen["shutdown"] * var(pm, n,   :z_gen_shutdown, i) + gen["startup"] * var(pm, n,   :z_gen_startup, i) for (i,gen) in nw_ref[:gen])
        for (n, nw_ref) in nws(pm))
    )
end

"objective function that minimizes the number of active branches"
function objective_min_branch(pm::AbstractPowerModel; kwargs...)

    return JuMP.@objective(pm.model, Min,
        sum(
            sum( var(pm, n,   :z_branch, i) for (i,branch) in nw_ref[:branch])
        for (n, nw_ref) in nws(pm))
    )

end

"returns a list of all on_off variables types associated with the `slack_types`"
function get_slack_var_symbols_on_off(pm::AbstractPowerModel, slack_types::Vector{Symbol})::Vector{Symbol}
    return (symbol -> Symbol("$(symbol)_on_off")).(get_slack_var_symbols(pm, slack_types))
end

"returns a list of all variables types associated with the `slack_types`"
function get_slack_var_symbols(pm::AbstractPowerModel, slack_types::Vector{Symbol})::Vector{Symbol}   
    var_symbols = append!(:pg                        ∈ slack_types ? [:pmin_slack,               :pmax_slack      ] : [],
                          :dclinepf                  ∈ slack_types ? [:p_dc_slack_lower,         :p_dc_slack_upper] : [], 
                          :ratea                     ∈ slack_types ? [                           :rate_a_slack    ] : [],
                          :ratea ∈ slack_types && !PARAMS[:CONTINGENCIES][:SLACK_EQUALS_SEC_SLACK] && haskey(var(pm,[n for (n, _) in nws(pm)][1]), :rate_a_slack_sec) ? [:rate_a_slack_sec] : [], # the key does not exist in Benders mode
                          :vad                       ∈ slack_types ? [:vad_slack_lower,          :vad_slack_upper ] : [],
                          :pd                        ∈ slack_types ? [:pd_slack_upper, :pd_slack_lower            ] : [])
    if typeof(pm) <: AbstractACPModel
                  append!(var_symbols,
                          :vm                        ∈ slack_types ? [:vm_slack_lower,           :vm_slack_upper  ] : [],
                          :qg                        ∈ slack_types ? [:qmin_slack,               :qmax_slack      ] : [])
    end
    
    return var_symbols
end

"returns the to the variable corresponding slack coefficient for relative coefficients"
function _get_slack_coeff(pm::AbstractPowerModel, n::Int, var_symbol::Symbol, vr, quadratic::Bool=false)
    split_name = split("$vr", "_")

    variable_name_slack_type = Dict(v => k for (k, vs) in CONSTS[:INTERNAL_MAPPINGS][:SLACK_VARIABLES][:SLACK_VARIABLE_INTERNAL_NAME_MAP] for v in vs)
    slack_type = variable_name_slack_type[var_symbol]
    slack_component_type = CONSTS[:INTERNAL_MAPPINGS][:SLACK_VARIABLES][:SLACK_COMPONENT_TYPE][slack_type]

    if slack_type == :vad && !quadratic
        return 10000
    end

    # get index
    idx = split(split_name[end][begin:end-1], "[")[2]
    if slack_component_type ∈ [:gen, :bus, :load]
        idx = [parse(Int, idx)]
    else
        if occursin("sec", "$split_name")
            idx = split(idx[begin+1:end-1], ",")
            idx = replace.(idx, "(" => "")
            #idx = replace.(idx, ")" => "")
            idx = (parse(Int, idx[1]), parse(Int, idx[3]), parse(Int, idx[4]))
        else
            idx = split(idx[begin+1:end-1], ",")
            idx = (parse(Int, idx[1]), parse(Int, idx[2]), parse(Int, idx[3]))
        end
    end

    if slack_type == :vad
        # compute coeff based on rhs of `constraint_ohms_yt_*_on_off_slack``
        branch = ref(pm, n, :branch)[idx[1]]
        g, b = calc_branch_y(branch)
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]

        va_fr = var(pm, n, :va, f_bus)
        va_to = var(pm, n, :va, t_bus)
        
        return 1/-b*(va_fr - va_to)
    else
        # get component
        component = ref(pm, n, slack_component_type)[idx[1]]
    end

    # get ub
    lb1, ub1, lb2, ub2 = _get_slack_bounds(pm, n, idx[1], component, slack_type)

    # return ub, for dcline we need to find out which ub       
    if slack_component_type == :dcline 
        if idx == (l, component["f_bus"], component["t_bus"]) # f_bus
            ub = ub1
        else # t_bus
            assert(idx == (l, component["t_bus"], component["f_bus"]))
            ub = ub2
        end
    else
        ub = ub1 # TODO - this should depend on the solver's float precision
    end

    return 1/(ub + 0.0001)
end

"objective function that minimizes the usage of slack variables, with larger individual deviations being punished more severely"
function objective_min_slack_absolute_quadratic(pm::AbstractPowerModel; kwargs...)
    return _objective_min_slack(pm::AbstractPowerModel; linear=false, relative=false, kwargs...)
end

"objective function that minimizes the usage of slack variables, with any amount of deviation punished equally"
function objective_min_slack_absolute_linear(pm::AbstractPowerModel; kwargs...)
    return _objective_min_slack(pm::AbstractPowerModel; linear=true, relative=false, kwargs...)
end

"objective function that minimizes the usage of slack variables, with larger individual deviations being punished more severely"
function objective_min_slack_relative_quadratic(pm::AbstractPowerModel; kwargs...)
    return _objective_min_slack(pm::AbstractPowerModel; linear=false, relative=true, kwargs...)
end

"objective function that minimizes the usage of slack variables, with any amount of deviation punished equally"
function objective_min_slack_relative_linear(pm::AbstractPowerModel; kwargs...)
    return _objective_min_slack(pm::AbstractPowerModel; linear=true, relative=true, kwargs...)
end

"objective function that minimizes maximum slack over all components"
function objective_min_max_slack(pm::AbstractPowerModel; kwargs...)
    var_max_slack = JuMP.@variable(pm.model, base_name="max_slack")
    
    # link var_max_slack with all slack variables
    slack_types = pm.data["slack_types"]
    for slack_type in slack_types
        var_symbols = get_slack_var_symbols(pm, [slack_type])

        for (n,nw) in nws(pm)
            varss = [get(var(pm, n), var_symbol, Dict()) for var_symbol in var_symbols]

            # vad and vm slack variables aren't also component types
            slack_component_type = CONSTS[:INTERNAL_MAPPINGS][:SLACK_VARIABLES][:SLACK_COMPONENT_TYPE][slack_type]

            # add constraints for all components corresponding to the slack type
            for i in ids(pm, n, slack_component_type)
                
                # get component ids
                if slack_component_type ∈ [:bus, :gen]
                    idx = i
                elseif slack_component_type ∈ [:branch]
                    idx = ref(pm, n, :arcs)[i]
                else #dcline
                    dcline = ref(pm, n, :dcline)[i]
                    idx = (i, dcline["f_bus"], dcline["t_bus"]) # f_idx
                    t_idx = (i, dcline["t_bus"], dcline["f_bus"])
                end

                # add constraints
                for vars in varss
                    if isempty(vars)
                        continue
                    end
                    
                    JuMP.@constraint(pm.model, vars[idx] <= var_max_slack)

                    if slack_component_type == :dcline
                        JuMP.@constraint(pm.model, vars[t_idx] <= var_max_slack)
                    end
                end
            end
        end
    end

    return JuMP.@objective(pm.model, Min, var_max_slack)
end


"objective function that minimizes the usage of slack variables"
function _objective_min_slack(pm::AbstractPowerModel; linear=true, relative=false, kwargs...)
    # choose variables that appear in the obj func depending on the types of slack variables we are considering
    slack_types = ref(pm, ismultinetwork(pm) ? 1 : 0, :slack_types)
    var_symbols = get_slack_var_symbols(pm, slack_types)

    # relative = relative to maximum, else absolute value
    if relative
        coeffs = [[[_get_slack_coeff(pm, n, symbol, v, false) for (i,v) in enumerate(var(pm,n,symbol))] for symbol in var_symbols] for (n,nw) in nws(pm)]
    else
        coeffs = [[[1.0 for (i,v) in enumerate(var(pm,n,symbol))] for symbol in var_symbols] for (n,nw) in nws(pm)]
    end

    # sum over all slack variables, with coeffients depending on `relative`, and power on `linear`
    return JuMP.@objective(pm.model, Min,
        sum( 
            sum( 
                sum(
                    (coeffs[n + (ismultinetwork(pm) ? 0 : 1)][s][v] * (typeof(vr) <: JuMP.VariableRef ? vr : vr[2]))^(linear ? 1 : 2) for (v,vr) in enumerate(var(pm,n,symbol))
                )
                for (s, symbol) in enumerate(var_symbols)
            )
            for (n,nw) in nws(pm) 
        )
    )
end

"minimizes the number of non-zero slack variables"
function objective_min_num_slacks(pm::AbstractPowerModel; kwargs...)
    slack_types = pm.data["slack_types"]
    #n_networks = length(nws(pm))
    if ! ismultinetwork(pm)
        networks = [0] 
    else
        networks = 1:length(nws(pm))
    end

    for n in networks
        variable_component_slack_on_off(pm, slack_types, nw=n)

        # link slack_on_off and slack
        for slack_type in slack_types
            slack_component_type = CONSTS[:INTERNAL_MAPPINGS][:SLACK_VARIABLES][:SLACK_COMPONENT_TYPE][slack_type]
                for i in ids(pm, n, slack_component_type)
                    constraint_component_slack_on_off(pm, i, slack_type, nw=n)
            end
        end
    end

    var_symbols = get_slack_var_symbols_on_off(pm, slack_types)

    return JuMP.@objective(pm.model, Min,
        sum( 
            sum( 
                sum(
                    vr for vr in var(pm,n,symbol)
                )
                for (s, symbol) in enumerate(var_symbols)
            )
            for (n,nw) in nws(pm) 
        )
    )
end

"empty objective"
function objective_zero(pm::AbstractPowerModel; kwargs...)
    return JuMP.@objective(pm.model, Min, 0)
end

"wrapper for single objective functions"
function _weighted_obj_func_factory(obj_func::Function)
    return obj_func
end

##### multi-objective #####

"wrapper for weighted objective functions"
function _weighted_obj_func_factory(obj_funcs::Tuple{Vector{<:Function}, Vector{<:Real}})
    return weighted_obj_func_factory(obj_funcs...)
end

"creates a weighted objective function without the need for additional functionality in the objective function functions"
function weighted_obj_func_factory(obj_funcs::Vector{<:Function}, weights::Vector{<:Real})
    weighted_objective_function = function (pm::AbstractPowerModel)
        w_obj_func = JuMP.@expression(pm.model, 0)
        for (obj_func, weight) in zip(obj_funcs, weights)
            if weight != 0
                # apply objective function to model, and overwrite the previous one
                obj_func(pm)

                # add to weighted objective function
                w_obj_func += weight * JuMP.objective_function(pm.model)
            end
        end

        if PARAMS[:CUTS][:NO_NEGATIVE_OBJ]  
            JuMP.@constraint(pm.model, w_obj_func >= 0)
        end
        
        return JuMP.@objective(pm.model, Min, w_obj_func)
    end
    return weighted_objective_function
end