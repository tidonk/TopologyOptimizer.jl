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

"Security-Constrained OTS with slack variables for a subset of {gen, branch, dcline, vad}"
function build_ots_slack_benders(pm::AbstractPowerModel; n=nw_id_default)
    slack_types = ref(pm, n, :slack_types)

    PARAMS[:LIMITS][:BRANCH][:BRANCH_SWITCHING] && variable_branch_indicator(  pm, nw=n)
    PARAMS[:LIMITS][:GENERATOR][:GEN_SWITCHING] && variable_gen_indicator(pm, nw=n)
    # only the vm can be (un)bounded here, not the va
    variable_bus_voltage_on_off(pm, bounded=:vm       ∉ slack_types, nw=n)
    PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_SHIFTING] && variable_transformer_shift(  pm, nw=n)

    variable_gen_power(         pm, bounded =false,                  nw=n) # bounded handeled separately for pg and qg
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

    constraint_model_voltage_on_off(pm, n)

    for i in ids(pm, n, :ref_buses)
        constraint_theta_ref(pm, i, nw=n)
    end

    for i in ids(pm, n, :bus)
        constraint_power_balance(pm, i, nw=n)
    end

    for i in ids(pm, n, :branch)
        if pm.data["BUSBARS"] && ! is_branch_electrical(pm, n, i)
            # no electrical properties to be enforced for auxiliary branches
        elseif is_branch_switchable(pm, n, i)
            constraint_ohms_yt_from_on_off(pm, i, nw=n)
            constraint_ohms_yt_to_on_off(pm, i, nw=n)
        else
            constraint_ohms_yt_from(pm, i, nw=n)
            constraint_ohms_yt_to(pm, i, nw=n)
        end

        if is_branch_switchable(pm, n, i)
            constraint_thermal_limit_from_slack_on_off(pm, i, nw=n)
            constraint_thermal_limit_to_slack_on_off(pm, i, nw=n)

            constraint_voltage_angle_difference_slack_on_off(pm, i, nw=n)
        else
            if is_branch_monitored(pm, n, i)
                constraint_thermal_limit_from_slack(pm, i, nw=n)
                constraint_thermal_limit_to_slack(pm, i, nw=n)
            end

            constraint_voltage_angle_difference_slack(pm, i, nw=n)
        end
    end

    for i in ids(pm, n, :dcline)
        constraint_dcline_power_losses(pm, i, nw=n)
    end
    
    if PARAMS[:LIMITS][:GENERATOR][:GEN_SWITCHING] 
        for i in ids(pm, n, :gen)
            # ensure that pg is only != 0 if the gen is on
            constraint_gen_power_on_off(pm, i, nw=n)
        end
    end

    return pm
end