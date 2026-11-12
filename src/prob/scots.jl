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

#### General Assumptions of these OTS Models ####
#
# - if the branch status is 0 in the input, it is out of service and forced to 0 in OTS
# - the network will be maintained as one connected component (i.e. at least n-1 edges)
#
using Gurobi
import JuMP
""
function solve_secure_ots(input::Any, model_type::Type, optimizer::Any, obj_funcs::Union{Function, Tuple{Vector{<:Function}, Vector{<:Real}}}; solver_specific_build_function=ssbf(pm)=pm, kwargs...)
    return solve_secure_model(input, model_type, optimizer, _weighted_obj_func_factory(obj_funcs) ∘ solver_specific_build_function ∘ build_secure_ots; ref_extensions=[ref_add_on_off_va_bounds!], kwargs...)
end

""
function solve_secure_ots_decomp(input::Any, model_type::Type, optimizer::Any, obj_funcs::Union{Function, Tuple{Vector{<:Function}, Vector{<:Real}}}; solver_specific_build_function=ssbf(pm)=pm, kwargs...)
    return solve_secure_model(input, model_type, optimizer, _weighted_obj_func_factory(obj_funcs) ∘ solver_specific_build_function ∘ build_secure_ots_decomp; ref_extensions=[ref_add_on_off_va_bounds!], kwargs...)
end

""
function solve_secure_ots_lazy(input::Any, model_type::Type, optimizer::Any, obj_funcs::Union{Function, Tuple{Vector{<:Function}, Vector{<:Real}}}; solver_specific_build_function=ssbf(pm)=pm, kwargs...)
    return solve_secure_model(input, model_type, optimizer, _weighted_obj_func_factory(obj_funcs) ∘ solver_specific_build_function ∘ build_secure_ots_lazy; ref_extensions=[ref_add_on_off_va_bounds!], kwargs...)
end

""
function solve_secure_ots_time(input::Any, model_type::Type, optimizer::Any, obj_funcs::Union{Function, Tuple{Vector{<:Function}, Vector{<:Real}}}; solver_specific_build_function=ssbf(pm)=pm, multinetwork=true, slack_types=[:pg, :ratea, :dclinepf, :vad], kwargs...)
    return solve_secure_model(input, model_type, optimizer, _weighted_obj_func_factory(obj_funcs) ∘ solver_specific_build_function ∘ build_secure_ots_time, slack_types; ref_extensions=[ref_add_on_off_va_bounds!], multinetwork=true, kwargs...)
end

""
function solve_secure_ots_slack(input::Any, model_type::Type, optimizer::Any, obj_funcs::Union{Function, Tuple{Vector{<:Function}, Vector{<:Real}}}; solver_specific_build_function=ssbf(pm)=pm, slack_types=[:pg, :qg, :ratea, :dclinepf, :vad], kwargs...)
    return solve_secure_model(input, model_type, optimizer, _weighted_obj_func_factory(obj_funcs) ∘ solver_specific_build_function ∘ build_secure_ots_slack, slack_types; ref_extensions=[ref_add_on_off_va_bounds!], kwargs...)
end

"hands input data to plotter without solving power flows etc"
function solve_plot_problem(input::Any, model_type::Type, optimizer::Any, obj_funcs::Union{Function, Tuple{Vector{<:Function}, Vector{<:Real}}}; solver_specific_build_function=ssbf(pm)=pm, slack_types=[:pg, :qg, :ratea, :dclinepf, :vad], kwargs...)
    return solve_secure_model(input, model_type, optimizer, _weighted_obj_func_factory(obj_funcs) ∘ solver_specific_build_function ∘ build_plot_problem, slack_types; ref_extensions=[ref_add_on_off_va_bounds!], kwargs...)
end
 
"Security-Constrained OTS for a single timestep"
function build_secure_ots(pm::AbstractPowerModel; n=nw_id_default)
    contingencies::Vector{Contingency} = ref(pm, n, :contingencies)

    bridges = Vector{Int}()
    for contingency in append!(contingencies, Contingency())
        for bridge in get_critical_branches(pm.data, contingency=contingency, nw=n) # todo: check if this works with multi-contingencies
            append!(bridges, get(bridge, :branch, []))
        end
    end
    pm.data["bridges"] = bridges
    
    variable_branch_indicator(pm, nw=n)
    variable_gen_indicator(pm, nw=n)
    variable_bus_voltage_on_off(pm, nw=n)
    variable_gen_power(pm, nw=n)
    variable_branch_power(pm, nw=n)
    variable_dcline_power(pm, nw=n)

    constraint_model_voltage_on_off(pm, nw=n)

    for i in ids(pm, n, :gen)
        gen = ref(pm, :gen, i, nw=n)
        # TODO: move this into a wrapper in `constraint_template.jl`
        constraint_gen_power_on_off(pm, n, i, gen["pmin"], gen["pmax"], gen["qmin"], gen["qmax"])
    end

    for i in ids(pm, n, :ref_buses)
        constraint_theta_ref(pm, i, nw=n)
    end

    for i in ids(pm, n, :bus)
        constraint_power_balance(pm, i, nw=n)
    end

    for i in ids(pm, n, :branch)
        if pm.data["BUSBARS"] && ! is_branch_electrical(pm, n, i)
                # for coupler or breaker, there is no thermal limit
        elseif is_branch_switchable(pm, n, i)
            constraint_ohms_yt_from_on_off(pm, i, nw=n)
            constraint_ohms_yt_to_on_off(pm, i, nw=n)
        else
            constraint_ohms_yt_from(pm, i, nw=n)
            constraint_ohms_yt_to(pm, i, nw=n)
        end

        if is_branch_switchable(pm, n, i)
            constraint_voltage_angle_difference_on_off(pm, i, nw=n)

            constraint_thermal_limit_from_on_off(pm, i, nw=n)
            constraint_thermal_limit_to_on_off(pm, i, nw=n)
        else
            constraint_voltage_angle_difference(pm, i, nw=n)

            if is_branch_monitored(pm, n, i)
                constraint_thermal_limit_from(pm, i, nw=n)
                constraint_thermal_limit_to(pm, i, nw=n)
            end
        end
    end

    for i in ids(pm, n, :dcline)
        constraint_dcline_power_losses(pm, i, nw=n)
    end

    # add constraints for all contingencies
    variable_bus_voltage_sec_on_off(pm, contingencies, nw=n)
    variable_gen_power_sec(pm, contingencies, nw=n)
    variable_branch_power_sec(pm, contingencies, nw=n)
    variable_dcline_power_sec(pm, contingencies, nw=n)

    constraint_model_voltage_sec_on_off(pm, contingencies, nw=n)

    # set subproblem flag in order to consider N-1 rating for branches. If SUBPROBLEM was true already, then there should be no additional contingencies
    subproblem_tmp = get(pm.data, "SUBPROBLEM", false)
    pm.data["SUBPROBLEM"] = true
    for (c, contingency) in enumerate(contingencies)
        for i in ids(pm, n, :ref_buses)
            #if i ∉ get(contingency, :bus, [])
                #constraint_theta_ref_sec(pm, i, c, contingency)
            #end
            # TODO: if i ∉ get(contingency, :bus, [])
            # what to do if a ref bus is in the contingency?
        end

        for i in ids(pm, n, :bus)
            if i ∉ get(contingency, :bus, [])
                constraint_power_balance_sec(pm, i, c, contingency, nw=n)
            end
        end

        for i in ids(pm, n, :branch)
            if i ∉ get(contingency, :branch, [])
                if pm.data["BUSBARS"] && ! is_branch_electrical(pm, n, i)
                    # no electrical properties to be enforced for auxiliary branches
                elseif is_branch_switchable(pm, n, i)
                    constraint_ohms_yt_from_sec_on_off(pm, i, c, contingency, nw=n)
                    constraint_ohms_yt_to_sec_on_off(pm, i, c, contingency, nw=n)
                else
                    constraint_ohms_yt_from_sec(pm, i, c, contingency, nw=n)
                    constraint_ohms_yt_to_sec(pm, i, c, contingency, nw=n)
                end
        
                if is_branch_switchable(pm, n, i)
                    constraint_voltage_angle_difference_sec_on_off(pm, i, c, contingency, nw=n)

                    constraint_thermal_limit_from_sec_on_off(pm, i, c, contingency, nw=n)
                    constraint_thermal_limit_to_sec_on_off(pm, i, c, contingency, nw=n)
                else
                    constraint_voltage_angle_difference_sec(pm, i, c, contingency, nw=n)

                    if is_branch_monitored(pm, n, i)
                        constraint_thermal_limit_from_sec(pm, i, c, contingency, nw=n)
                        constraint_thermal_limit_to_sec(pm, i, c, contingency, nw=n)
                    end
                end
            end
        end
        
        for i in ids(pm, n, :dcline)
            if i ∉ get(contingency, :dcline, [])
               constraint_dcline_power_losses_sec(pm, i, c, contingency, nw=n)
            end
        end
    end
    pm.data["SUBPROBLEM"] = subproblem_tmp

    if ! ismultinetwork(pm)
        # add constraints to limit the total number of possible actions
        constraint_top_actions_limit(pm)
    end

    return pm
end

"Security-Constrained OTS with slack variables for a subset of {gen, branch, dcline, vad}"
function build_secure_ots_slack(pm::AbstractPowerModel; n=nw_id_default, with_base_pf::Bool=true)
    contingencies::Vector{Contingency} = ref(pm, n, :contingencies)
    slack_types = ref(pm, n, :slack_types)

    PARAMS[:LIMITS][:BRANCH][:BRANCH_SWITCHING] && variable_branch_indicator(pm, nw=n)
    PARAMS[:LIMITS][:GENERATOR][:GEN_SWITCHING] && variable_gen_indicator(pm, nw=n)
    # only the vm can be (un)bounded here, not the va
    with_base_pf && variable_bus_voltage_on_off(pm, bounded=:vm       ∉ slack_types, nw=n)
    PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_SHIFTING] && variable_transformer_shift(  pm, nw=n)

    variable_gen_power(         pm,                                  nw=n) # bounded handeled separately for pg and qg
    if with_base_pf
        variable_branch_power(      pm, bounded=:ratea    ∉ slack_types, nw=n)
        variable_dcline_power(      pm, bounded=:dclinepf ∉ slack_types, nw=n)

        # impose limits for all parameters that were previously made unbounded
        variable_component_slack(pm, slack_types, nw=n)
        for slack_type in slack_types
            if slack_type ∉ [:ratea, :vad, :pd]
                for i in ids(pm, n, CONSTS[:INTERNAL_MAPPINGS][:SLACK_VARIABLES][:SLACK_COMPONENT_TYPE][slack_type])
                    constraint_component_slack(pm, i, slack_type, nw=n)
                end
            end
        end
    end

    # only one active breaker branch per breaker bus
    if PARAMS[:LIMITS][:BRANCH][:BRANCH_SWITCHING]
        for i in ids(pm, n, :bus)
            if pm.data["BUSBARS"] && get(ref(pm, n, :bus, i), "type", "") ∈ ["breaker"]
                constraint_breaker_branch_limit(pm, i, nw=n)
            end
        end
    end

    if with_base_pf
        constraint_model_voltage_on_off(pm, n)

        for i in ids(pm, n, :ref_buses)
            constraint_theta_ref(pm, i, nw=n)
        end

        for i in ids(pm, n, :bus)
            constraint_power_balance(pm, i, nw=n)
        end

        for i in ids(pm, n, :branch)
            if pm.data["BUSBARS"] && ! is_branch_electrical(pm, n, i)
                    # for coupler or breaker, there is no thermal limit
            elseif is_branch_switchable(pm, n, i)
                constraint_ohms_yt_from_on_off(pm, i, nw=n)
                constraint_ohms_yt_to_on_off(pm, i, nw=n)
            else
                constraint_ohms_yt_from(pm, i, nw=n)
                constraint_ohms_yt_to(pm, i, nw=n)
            end

            if is_branch_switchable(pm, n, i)
                constraint_voltage_angle_difference_slack_on_off(pm, i, nw=n)

                constraint_thermal_limit_from_slack_on_off(pm, i, nw=n)
                constraint_thermal_limit_to_slack_on_off(pm, i, nw=n)
            else
                constraint_voltage_angle_difference_slack(pm, i, nw=n)

                if is_branch_monitored(pm, n, i)
                    constraint_thermal_limit_from_slack(pm, i, nw=n)
                    constraint_thermal_limit_to_slack(pm, i, nw=n)
                end
            end
        end

        for i in ids(pm, n, :dcline)
            constraint_dcline_power_losses(pm, i, nw=n)
        end
    else
        # no base PF: aggregate power-balance over all components (replaces per-bus constraint_power_balance)
        constraint_power_balance_all(pm, nw=n)
    end
    
    if PARAMS[:LIMITS][:GENERATOR][:GEN_SWITCHING] 
        for i in ids(pm, n, :gen)
            # ensure that pg is only != 0 if the gen is on
            constraint_gen_power_on_off(pm, i, nw=n)
        end
    end

    if ! ismultinetwork(pm)
        # add constraints to limit the total number of possible actions
        constraint_top_actions_limit(pm)
    end
    
    if isempty(contingencies)
        return pm
    end
    # add constraints for all contingencies
    variable_bus_voltage_sec_on_off(pm, contingencies, nw=n)
    variable_gen_power_sec(   pm, contingencies,                                  nw=n)
    variable_branch_power_sec(pm, contingencies, bounded=:ratea    ∉ slack_types, nw=n)
    variable_dcline_power_sec(pm, contingencies, bounded=:dclinepf ∉ slack_types, nw=n)

    if :ratea ∈ slack_types
        variable_rate_a_slack_sec(pm, contingencies, nw=n)
    end

    # set subproblem flag in order to consider N-1 rating for branches. If SUBPROBLEM was true already, then there should be no additional contingencies
    subproblem_tmp = get(pm.data, "SUBPROBLEM", false)
    pm.data["SUBPROBLEM"] = true
    # impose limits for all parameters that were previously made unbounded
    for (c, contingency) in enumerate(contingencies)

        for i in ids(pm, n, :ref_buses)
            constraint_theta_ref_sec(pm, i, c, contingency, nw=n)
        end

        for slack_type in slack_types
            if slack_type ∉ [:ratea, :vad, :pd]
                slack_component_type = CONSTS[:INTERNAL_MAPPINGS][:SLACK_VARIABLES][:SLACK_COMPONENT_TYPE][slack_type]
                for i in ids(pm, n, slack_component_type)
                    if i ∉ get(contingency, slack_component_type, [])
                        constraint_component_slack_sec(pm, i, c, contingency, slack_type, nw=n)
                    end
                end
            end
        end

        for i in ids(pm, n, :bus)
            if i ∉ get(contingency, :bus, [])
                constraint_power_balance_sec(pm, i, c, contingency, nw=n)
            end
        end

        for i in ids(pm, n, :branch)
            if i ∉ get(contingency, :branch, [])
                if pm.data["BUSBARS"] && ! is_branch_electrical(pm, n, i)
                    # no electrical properties to be enforced for auxiliary branches
                elseif is_branch_switchable(pm, n, i)
                    constraint_ohms_yt_from_sec_on_off(pm, i, c, contingency, nw=n)
                    constraint_ohms_yt_to_sec_on_off(pm, i, c, contingency, nw=n)
                else
                    constraint_ohms_yt_from_sec(pm, i, c, contingency, nw=n)
                    constraint_ohms_yt_to_sec(pm, i, c, contingency, nw=n)
                end

                if is_branch_switchable(pm, n, i)
                    constraint_voltage_angle_difference_slack_sec_on_off(pm, i, c, contingency, nw=n)

                    constraint_thermal_limit_from_slack_sec_on_off(pm, i, c, contingency, nw=n)
                    constraint_thermal_limit_to_slack_sec_on_off(pm, i, c, contingency, nw=n)
                else
                    constraint_voltage_angle_difference_slack_sec(pm, i, c, contingency, nw=n)

                    if is_branch_monitored(pm, n, i)
                        constraint_thermal_limit_from_slack_sec(pm, i, c, contingency, nw=n)
                        constraint_thermal_limit_to_slack_sec(pm, i, c, contingency, nw=n)
                    end
                end
            end
        end

        for i in ids(pm, n, :dcline)
            if i ∉ get(contingency, :dcline, [])
               constraint_dcline_power_losses_sec(pm, i, c, contingency, nw=n)
            end
        end
    end
    pm.data["SUBPROBLEM"] = subproblem_tmp

    return pm
end

"Security-Constrained OTS with time expansion to account for multiple time steps in the optimization. `with_base_pf = false` omits the base-case powerflow variables (pf & va) and their constraints. Per-contingency PFs are added for whatever is currently in `ref(pm, t, :contingencies)` (empty → no contingency PFs)."
function build_secure_ots_time(pm::AbstractPowerModel; with_base_pf::Bool=true)
    n_timesteps = length(nws(pm))
    # for each timestep, extend the secure OTS 
    for timestep in 1:n_timesteps
        pm = build_secure_ots_slack(pm::AbstractPowerModel; n=timestep, with_base_pf=with_base_pf)

        # binary decision variable whether a gen is active in a timestep or not. incurs shutdown/startup costs.
        #variable_gen_indicator(pm, nw=timestep)
        # for each branch l in each time step, add z_branch_change[l] that is 1 iff the value of z_branch[l] changed in the current time step
        PARAMS[:LIMITS][:BRANCH][:BRANCH_SWITCHING] && variable_branch_change_indicator(pm, nw=timestep)
        # for each branch l that is a transformer in each time step, add tr_shift_change[l] that is 1 iff the value of tr_shift[l] changed in the current time step
        PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_SHIFTING] && variable_transformer_shift_change_indicator(pm, nw=timestep)
        # for each gen l in each time step, add z_gen_change[l] that is 1 iff the value of z_gen[l] changed in the current time step
        PARAMS[:LIMITS][:GENERATOR][:GEN_SWITCHING] && variable_gen_change_indicator(pm, nw=timestep, relax=true)
        # for each gen g in each time step, add pg_delta[g] that is between (1-\delta_g) * pg_0 and (1+\delta_g) * pg_0.
        variable_gen_power_change(pm, nw=timestep, bounded=false)
        PARAMS[:LIMITS][:REDISPATCH][:PG_DELTA_NUM_DAY]  < Inf && variable_gen_power_change_on_off(pm, nw=timestep)
        
        for i in ids(pm, timestep, :gen)
            # link z_gen to z_gen_startup and z_gen_shutdown
            PARAMS[:LIMITS][:GENERATOR][:GEN_SWITCHING]  && constraint_gen_change_linking(pm, timestep, i)
            # ensure that pg is only != 0 if the gen is on
            #constraint_gen_power_on_off(pm, i, nw=timestep)
            # that pg is exactly pg_0 + pg_delta for all gens
            constraint_gen_power_delta_linking(pm, i, nw=timestep)
            # ensure that pg_delta is only != 0 iff pg_delta_on_off == 1
            PARAMS[:LIMITS][:REDISPATCH][:PG_DELTA_NUM_DAY] < Inf && constraint_gen_power_change_on_off(pm, timestep, i)
            # ensure that the previous pg is not too different from the current pg (obey gradients/ramps)
            #! PARAMS[:LIMITS][:GENERATOR][:GEN_REDISPATCH_OVERRIDE] && constraint_gen_power_change_linking(pm, timestep, i)
        end

        for i in ids(pm, timestep, :branch)
            if PARAMS[:LIMITS][:BRANCH][:BRANCH_SWITCHING] && is_branch_switchable(pm, timestep, i)
                constraint_branch_change_linking(pm, i, nw=timestep)
            end
            if PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_SHIFTING] && is_branch_shiftable(pm, timestep, i)
                constraint_transformer_shift_change(pm, i, nw=timestep)
                constraint_transformer_shift_change_linking(pm, i, nw=timestep)
            end
        end

        # restrict the number of topological changes per time step
        PARAMS[:LIMITS][:BRANCH][:BRANCH_SWITCHING] && constraint_branch_max_changes_timestep(pm, timestep)
        # restrict the number of topological changes per time step
        PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_SHIFTING] && constraint_transformer_max_changes_timestep(pm, timestep)
        
        # ensure that auxiliary gens to model GSKs are linked
        constraint_gen_power_gsk_linking(pm, nw=timestep)
        
    end

    # restrict the total number of topological changes
    PARAMS[:LIMITS][:BRANCH][:BRANCH_SWITCHING] && constraint_branch_max_changes_total(pm, n_timesteps)
    # restrict the total number of transformer shift changes
    PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_SHIFTING]     && constraint_transformer_max_changes_total(pm, n_timesteps)
    # restrict the total number of re-dispatches
    PARAMS[:LIMITS][:REDISPATCH][:PG_DELTA_NUM_DAY]     < Inf && constraint_gen_max_changes_total(pm, n_timesteps)
    # restrict the total number of startups and shutdowns
    PARAMS[:LIMITS][:GENERATOR][:GEN_STARTUP_NUM_DAY]  ≥ 0 && constraint_gen_max_startup_total(pm, n_timesteps)
    PARAMS[:LIMITS][:GENERATOR][:GEN_SHUTDOWN_NUM_DAY] ≥ 0 && constraint_gen_max_shutdown_total(pm, n_timesteps)

    # restrict the total number of changes per component
    for i in ids(pm, 1, :branch)
        if PARAMS[:LIMITS][:BRANCH][:BRANCH_SWITCHING] && is_branch_switchable(pm, 1, i)
            constraint_branch_max_changes_component(pm, n_timesteps, i)
            if length(pm.data["nw"]) > 1
                constraint_switch_duration(pm, n_timesteps, i)
            end
        end
        if PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_SHIFTING] && is_branch_shiftable(pm, 1, i)
            constraint_transformer_max_changes_component(pm, n_timesteps, i)
        end
    end

    # minimum startup / shutdown durations
    for (i, gen) in ref(pm, 1, :gen)
        if is_gen_switchable(gen)
            PARAMS[:LIMITS][:GENERATOR][:GEN_MIN_UP]   ≥ 1 && constraint_startup_minup(pm, n_timesteps, i)
            PARAMS[:LIMITS][:GENERATOR][:GEN_MIN_DOWN] ≥ 1 && constraint_shutdown_mindown(pm, n_timesteps, i)
        end
    end

    # link all slacks. can be disabled if pmin/pmax/ratea/... are different in each step and thus could require a different slack value
    #constraint_slack_linking(pm, n_timesteps)
    return pm
end

"Master problem (gen) for nested decomposition"
function build_gen_model_time(pm::AbstractPowerModel; with_pf::Bool=true)
    n_timesteps = length(nws(pm))
    # for each timestep, extend the secure OTS 
    for timestep in 1:n_timesteps
        pm = build_gen_model(pm::AbstractPowerModel, n=timestep, with_pf=with_pf)
    end

    # restrict the total number of re-dispatches
    PARAMS[:LIMITS][:REDISPATCH][:PG_DELTA_NUM_DAY]     < Inf && constraint_gen_max_changes_total(pm, n_timesteps)

    return pm
end

"Master problem (switching) for nested decomposition"
function build_switching_model_time(pm::AbstractPowerModel; with_pf::Bool=true)
    n_timesteps = length(nws(pm))
    # for each timestep, extend the secure OTS 
    for timestep in 1:n_timesteps
        build_switching_model(pm::AbstractPowerModel; n=timestep, with_pf=with_pf)
    end

    if n_timesteps > 1
        # restrict the total number of topological changes
        constraint_branch_max_changes_total(pm, n_timesteps)
        # restrict the total number of transformer shift changes
        PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_SHIFTING]     && constraint_transformer_max_changes_total(pm, n_timesteps)

        # restrict the total number of changes per component
        for i in ids(pm, 1, :branch)
            if is_branch_switchable(pm, 1, i)
                constraint_branch_max_changes_component(pm, n_timesteps, i)
                if length(pm.data["nw"]) > 1
                    constraint_switch_duration(pm, n_timesteps, i)
                end
            end
            if is_branch_shiftable(pm, 1, i)
                constraint_transformer_max_changes_component(pm, n_timesteps, i)
            end
        end
    end

    return pm
end

"Its only purpose is to plot :c"
function build_plot_problem(pm::AbstractPowerModel; n=nw_id_default)
    contingencies::Vector{Contingency} = Vector{Contingency}()#ref(pm, n, :contingencies)

    #bridges = Vector{Int}()
    #for contingency in append!(contingencies, Contingency())
    #    for bridge in get_critical_branches(pm.data, contingency=contingency, nw=n) # todo: check if this works with multi-contingencies
    #        append!(bridges, get(bridge, :branch, []))
    #    end
    #end
    #pm.data["bridges"] = bridges
    
    variable_branch_indicator(pm, nw=n)
    variable_gen_indicator(pm, nw=n)
    variable_bus_voltage_on_off(pm, nw=n)
    variable_gen_power(pm, nw=n)
    variable_branch_power(pm, nw=n)
    variable_dcline_power(pm, nw=n)

    # add constraints for all contingencies
    variable_bus_voltage_sec_on_off(pm, contingencies, nw=n)
    variable_gen_power_sec(pm, contingencies, nw=n)
    variable_branch_power_sec(pm, contingencies, nw=n)
    variable_dcline_power_sec(pm, contingencies, nw=n)

    return pm
end