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

"Gen model (for nested Benders decomposition). Includes branch, trafo, gen switching vars that are unconstrained and will be fixed by the nested Benders algorithm"
function build_gen_model(pm::AbstractPowerModel; n::Int=nw_id_default, with_pf::Bool=false, include_switching_vars::Bool=true)

    if include_switching_vars
        if PARAMS[:LIMITS][:BRANCH][:BRANCH_SWITCHING]
            variable_branch_indicator(pm, nw=n, relax=true)
        end

        if PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_SHIFTING]
            variable_transformer_shift(pm, nw=n, discrete=false)
        end

        if PARAMS[:LIMITS][:GENERATOR][:GEN_SWITCHING] 
            variable_gen_indicator(pm, nw=n, relax=true)
        end
    end

    if true#REDISPATCHING
        variable_gen_power(pm, nw=n, bounded=false)
        variable_gen_power_change(pm, nw=n)
    
        constraint_gen_power_gsk_linking(pm, nw=n) # ensure that auxiliary gens to model GSKs are linked
        for i in ids(pm, n, :gen)
            constraint_gen_power_delta_linking(pm, i, nw=n) # that pg is exactly pg_0 + pg_delta for all gens
            constraint_gen_power_change_linking(pm, n, i) # obey gradients
        end

        if PARAMS[:LIMITS][:GENERATOR][:GEN_SWITCHING] 
            for i in ids(pm, n, :gen)
                constraint_gen_power_on_off(pm, i, nw=n) # ensure that pg is only != 0 if the gen is on
            end
        end
    end
    
    if ! with_pf
        constraint_power_balance_all(pm, nw=n)
    else
        variable_bus_voltage_on_off(pm, nw=n)
        variable_branch_power(pm, bounded=:ratea ∉ ref(pm, n, :slack_types), nw=n)
        variable_dcline_power(pm, nw=n)

        # impose limits for all parameters that were previously made unbounded
        variable_component_slack(pm, ref(pm, n, :slack_types), nw=n)

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
    end



    return pm
end