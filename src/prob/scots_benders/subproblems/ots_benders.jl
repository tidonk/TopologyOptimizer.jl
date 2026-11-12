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

"for benders subproblems. does not allow for switching gens."
function build_ots_benders(pm::AbstractPowerModel; relax=false, nw=nw_id_default)
    variable_branch_indicator(pm, nw=nw, relax=relax)
    variable_bus_voltage_on_off(pm, nw=nw)
    variable_gen_power(pm, nw=nw)
    variable_branch_power(pm, bounded=false, nw=nw)
    variable_dcline_power(pm, nw=nw)
    variable_transformer_shift(  pm, nw=nw)

    constraint_model_voltage_on_off(pm, nw=nw)

    for i in ids(pm, nw, :ref_buses)
        constraint_theta_ref(pm, i, nw=nw)
    end

    for i in ids(pm, nw, :bus)
        constraint_power_balance(pm, i, nw=nw)
    end

    for i in ids(pm, nw, :branch)
        if pm.data["BUSBARS"] && ! is_branch_electrical(pm, nw, i)
            # no electrical properties to be enforced for auxiliary branches
        elseif is_branch_switchable(pm, nw, i)
            constraint_ohms_yt_from_on_off(pm, i, nw=nw)
            constraint_ohms_yt_to_on_off(pm, i, nw=nw)
        else
            constraint_ohms_yt_from(pm, i, nw=nw)
            constraint_ohms_yt_to(pm, i, nw=nw)
        end

        # thermal and angle limits are calculated after solving the MIP, so special constraints on them are not needed (except for on/off)
        constraint_branch_from_on_off(pm, i, nw=nw)
        constraint_branch_to_on_off(pm, i, nw=nw)
    end

    for i in ids(pm, nw, :dcline)
        constraint_dcline_power_losses(pm, i, nw=nw)
    end

    return pm
end