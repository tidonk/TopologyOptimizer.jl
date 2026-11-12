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

"OTS without powerflow (to be used with benders)"
function build_ots_no_pf(pm::AbstractPowerModel; n=nw_id_default)
    PARAMS[:LIMITS][:BRANCH][:BRANCH_SWITCHING] && variable_branch_indicator(pm, nw=n)
    PARAMS[:LIMITS][:GENERATOR][:GEN_SWITCHING] && variable_gen_indicator(pm, nw=n)
    # only the vm can be (un)bounded here, not the va
    PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_SHIFTING] && variable_transformer_shift(pm, nw=n)

    variable_gen_power(pm, nw=n) # bounded handeled separately for pg and qg

    # only one active breaker branch per breaker bus
    if PARAMS[:LIMITS][:BRANCH][:BRANCH_SWITCHING]
        for i in ids(pm, n, :bus)
            if pm.data["BUSBARS"] && get(ref(pm, n, :bus, i), "type", "") ∈ ["breaker"]
                constraint_breaker_branch_limit(pm, i, nw=n)
            end
        end
    end

    constraint_power_balance_all(pm, nw=n)

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

    return pm
end