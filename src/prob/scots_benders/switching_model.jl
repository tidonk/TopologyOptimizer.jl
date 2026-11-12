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

"Switching model (for nested Benders decomposition)"
function build_switching_model(pm::AbstractPowerModel; n::Int=nw_id_default, with_pf::Bool=false)
    if PARAMS[:LIMITS][:BRANCH][:BRANCH_SWITCHING]
        variable_branch_indicator(pm, nw=n)
        variable_branch_change_indicator(pm, nw=n)
        constraint_branch_max_changes_timestep(pm, n)

        for i in ids(pm, n, :branch)
            if is_branch_switchable(pm, n, i)
                constraint_branch_change_linking(pm, i, nw=n)
            end
        end

        # only one active breaker branch per breaker bus
        for i in ids(pm, n, :bus)
            if pm.data["BUSBARS"] && get(ref(pm, n, :bus, i), "type", "") ∈ ["breaker"]
                constraint_breaker_branch_limit(pm, i, nw=n)
            end
        end
    end

    if PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_SHIFTING]
        variable_transformer_shift(pm, nw=n)
        variable_transformer_shift_change_indicator(pm, nw=n)
        constraint_transformer_max_changes_timestep(pm, n)  
        
        for i in ids(pm, n, :branch)
            if is_branch_shiftable(pm, n, i)
                constraint_transformer_shift_change(pm, i, nw=n)
                constraint_transformer_shift_change_linking(pm, i, nw=n)
            end 
        end
    end

    if PARAMS[:LIMITS][:GENERATOR][:GEN_SWITCHING] 
        variable_gen_indicator(pm, nw=n)
        variable_gen_change_indicator(pm, nw=n, relax=true)

        # restrict the total number of startups and shutdowns
        PARAMS[:LIMITS][:GENERATOR][:GEN_STARTUP_NUM_DAY]  ≥ 0 && constraint_gen_max_startup_total(pm, length(nws(pm)))
        PARAMS[:LIMITS][:GENERATOR][:GEN_SHUTDOWN_NUM_DAY] ≥ 0 && constraint_gen_max_shutdown_total(pm, length(nws(pm)))

        # minimum startup / shutdown durations
        for i in ids(pm, n, :gen)
            constraint_gen_change_linking(pm, n, i) # link z_gen to z_gen_startup and z_gen_shutdown
            PARAMS[:LIMITS][:GENERATOR][:GEN_MIN_UP]   ≥ 1 && constraint_startup_minup(pm, length(nws(pm)), i)
            PARAMS[:LIMITS][:GENERATOR][:GEN_MIN_DOWN] ≥ 1 && constraint_shutdown_mindown(pm, length(nws(pm)), i)
        end

        # cuts
        constraint_power_balance_switching(pm, nw=n)
    end

    # on top of switching_model, add gen_model (without switching vars) to determine generation in switching master, too
    if PARAMS[:BENDERS][:NESTED][:GENERATION_IN_MP]
        build_gen_model(pm, n=n, include_switching_vars=false, with_pf=with_pf)
    end

    return pm
end