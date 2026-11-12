function solve_secure_opf_slack(file, model_type::Type, optimizer; kwargs...)
    return solve_model(file, model_type, optimizer, build_secure_opf_slack; kwargs...)
end

"AC OPf with slack and (selected) contengencies (AC-SCOPF-SLACK) [neccessary to pass selected contingencies]"
function build_secure_opf_slack(pm::AbstractPowerModel; nw::Int=nw_id_default)

    variable_bus_voltage(pm, nw=nw)
    variable_gen_power(pm, nw=nw)
    variable_branch_power(pm, nw=nw)
    PARAMS[:LIMITS][:GENERATOR][:GEN_SWITCHING] && variable_gen_indicator(pm, nw=nw)
    variable_rate_a_slack(pm, nw=nw)
    variable_dcline_power(pm, nw=nw)
    variable_gen_power_change(pm, nw=nw, bounded=false)

    contingencies = pm.data["chosen_contingencies"]
    slack_types = ref(pm, ismultinetwork(pm) ? 1 : 0, :slack_types)
    if :pd in slack_types
        variable_pd_slack(pm, nw=nw)
        if !isempty(contingencies)
            variable_pd_slack_sec(pm, contingencies, nw=nw)
        end
    end
    if :qd in slack_types
        variable_qd_slack(pm, nw=nw)
        if !isempty(contingencies)
            variable_qd_slack_sec(pm, contingencies, nw=nw)
        end
    end
    
    weighted_obj_func_factory( 
        [objective_min_redispatch_cost, objective_min_num_redispatch, objective_min_slack_relative_linear, objective_min_num_slacks],
        [pm.data["cost_coeffs"]["redispatch"], pm.data["cost_coeffs"]["num_redispatch"], pm.data["cost_coeffs"]["overload"], pm.data["cost_coeffs"]["num_overload"]]
    )(pm)
   
    if PARAMS[:LIMITS][:GENERATOR][:GEN_SWITCHING] 
        for i in ids(pm, nw, :gen)
            # ensure that pg is only != 0 if the gen is on
            constraint_gen_power_on_off(pm, i, nw=nw)
        end
    end

    constraint_model_voltage(pm, nw)

    for i in ids(pm, nw, :ref_buses)
        constraint_theta_ref(pm, i, nw=nw)
    end

    for i in ids(pm, nw, :bus)
        constraint_power_balance(pm, i, nw=nw)
    end

    for i in ids(pm, nw, :branch)
        if pm.data["BUSBARS"] && ! _is_branch_electrical(pm, nw, i)
            # for coupler or breaker, there is no thermal limit
        else
            @assert(ref(pm, nw, :branch, i)["angmin"] != 0, "angmin is not defined for branch $i: $(ref(pm, nw, :branch, i))")
            constraint_ohms_yt_from(pm, i, nw=nw)
            constraint_ohms_yt_to(pm, i, nw=nw)
        end

        constraint_voltage_angle_difference(pm, i, nw=nw)

        if is_branch_monitored(pm, nw, i)
            constraint_thermal_limit_from_slack(pm, i, nw=nw)
            constraint_thermal_limit_to_slack(pm, i, nw=nw)
        end
    end

    for i in ids(pm, nw, :dcline)
        constraint_dcline_power_losses(pm, i, nw=nw)
    end

    for i in ids(pm, nw, :gen)
        constraint_gen_power_delta_linking(pm, i, nw=nw)
    end

    # contingeny constraints
    if isempty(contingencies)
        return pm
    end
    variable_bus_voltage_sec(pm, contingencies, nw=nw)
    variable_gen_power_sec(pm, contingencies, nw=nw)
    variable_branch_power_sec(pm, contingencies, nw=nw)
    variable_dcline_power_sec(pm, contingencies, nw=nw)
    variable_rate_a_slack_sec(pm, contingencies, nw=nw)

    constraint_model_voltage_sec(pm, nw, contingencies)
    
    for (c, contingency) in enumerate(contingencies)

        for i in ids(pm, nw, :bus)
            if i ∉ get(contingency, :bus, [])
                constraint_power_balance_sec(pm, i, c, contingency, nw=nw)
            end
        end
        
        for i in ids(pm, nw, :branch)
            if i ∉ get(contingency, :branch, [])
                if pm.data["BUSBARS"] && ! is_branch_electrical(pm, nw, i)
                    # no electrical properties to be enforced for auxiliary branches
                else
                    constraint_ohms_yt_from_sec(pm, i, c, contingency, nw=nw)
                    constraint_ohms_yt_to_sec(pm, i, c, contingency, nw=nw)
                end
        
                constraint_voltage_angle_difference_sec(pm, i, c, contingency, nw=nw)

                if is_branch_monitored(pm, nw, i)
                    constraint_thermal_limit_from_slack_sec(pm, i, c, contingency, nw=nw)
                    constraint_thermal_limit_to_slack_sec(pm, i, c, contingency, nw=nw)
                end
            end
        end
        
        for i in ids(pm, nw, :dcline)
            if i ∉ get(contingency, :dcline, [])
               constraint_dcline_power_losses_sec(pm, i, c, contingency, nw=nw)
            end
        end
    end
    
    return pm 
end