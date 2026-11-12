# This file contains modified code derived from PowerModels.jl.
# Copyright (c) 2016, Los Alamos National Security, LLC. All rights reserved.
# These portions retain the license and disclaimer in LICENSES/PowerModels.txt.
# PowerModels.jl code in this file has been modified for TopologyOptimizer.jl.

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

### thermal limit
## from

"`p[f_idx]^2 + q[f_idx]^2 <= rate_a^2`"
function constraint_thermal_limit_from_slack(pm::AbstractPowerModel, n::Int, f_idx, rate_a)
    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)
    rate_a_slack = var(pm, n, :rate_a_slack, f_idx)

    JuMP.@constraint(pm.model, p_fr^2 + q_fr^2 <= (rate_a + rate_a_slack)^2, base_name = "$(n)_rating_from_slack_$(f_idx)") 
end

"`p[f_idx]^2 + q[f_idx]^2 <= (rate_a * z_branch[i])^2`"
function constraint_thermal_limit_from_slack(pm::AbstractPowerModel, n::Int, i, f_idx, rate_a)
    @assert(:ratea ∈ ref(pm, n, :slack_types))

    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)

    rate_a_slack = var(pm, n, :rate_a_slack, f_idx)
    JuMP.@constraint(pm.model, p_fr^2 + q_fr^2 <= (rate_a + rate_a_slack)^2, base_name = "$(n)_rating_from_slack_$(i)")
end

"`p[f_idx]^2 + q[f_idx]^2 <= (rate_a * z_branch[i])^2`"
function constraint_thermal_limit_from_slack_on_off(pm::AbstractPowerModel, n::Int, i, f_idx, rate_a)
    @assert(:ratea ∈ ref(pm, n, :slack_types))

    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)
    z = var(pm, n, :z_branch, i)
    
    p_max = ref(pm, n, :p_max)

    rate_a_slack = var(pm, n, :rate_a_slack, f_idx)
    JuMP.@constraint(pm.model, p_fr^2 + q_fr^2 <= (rate_a)^2*z*(p_max)) # switching
    JuMP.@constraint(pm.model, p_fr^2 + q_fr^2 <= (rate_a + rate_a_slack)^2) # bound
end

"`p[f_idx]^2 + q[f_idx]^2 <= rate_a^2`"
function constraint_thermal_limit_from_sec(pm::AbstractPowerModel, n::Int, i::Int, f_idx, rate_a)
    p_fr_sec = var(pm, n, :p_sec, f_idx)
    q_fr_sec = var(pm, n, :q_sec, f_idx)
    
    JuMP.@constraint(pm.model, p_fr_sec^2 + q_fr_sec^2 <= rate_a^2, base_name = "$(n)_rating_from_sec_$(f_idx)")
end

"`p[f_idx]^2 + q[f_idx]^2 <= (rate_a * z_branch[i])^2`"
function constraint_thermal_limit_from_sec_on_off(pm::AbstractPowerModel, n::Int, i, f_idx, rate_a)
    p_sec_fr = var(pm, n, :p_sec, f_idx)
    q_sec_fr = var(pm, n, :q_sec, f_idx)
    z = var(pm, n, :z_branch, i)

    JuMP.@constraint(pm.model, p_sec_fr^2 + q_sec_fr^2 <= rate_a^2*z^2)
end

"`p[f_idx]^2 + q[f_idx]^2 <= (rate_a * z_branch[i])^2`"
function constraint_thermal_limit_from_slack_sec(pm::AbstractPowerModel, n::Int, i, f_idx, rate_a)
    @assert(:ratea ∈ ref(pm, n, :slack_types))

    p_sec_fr = var(pm, n, :p_sec, f_idx)
    q_sec_fr = var(pm, n, :q_sec, f_idx)

    (l,c),i,j = f_idx

    rate_a_slack = var(pm, n, :rate_a_slack_sec, f_idx)
    JuMP.@constraint(pm.model, p_sec_fr^2 + q_sec_fr^2 <= (rate_a + rate_a_slack)^2, base_name = "$(n)_rating_from_slack_$(f_idx)") # bound
end

"`p[f_idx]^2 + q[f_idx]^2 <= (rate_a * z_branch[i])^2`"
function constraint_thermal_limit_from_slack_sec_on_off(pm::AbstractPowerModel, n::Int, i, f_idx, rate_a)
    @assert(:ratea ∈ ref(pm, n, :slack_types))

    p_sec_fr = var(pm, n, :p_sec, f_idx)
    q_sec_fr = var(pm, n, :q_sec, f_idx)
    z = var(pm, n, :z_branch, i)
    
    p_max = ref(pm, n, :p_max)

    (l,c),i,j = f_idx

    rate_a_slack = var(pm, n, :rate_a_slack_sec, f_idx)
    JuMP.@constraint(pm.model, p_sec_fr^2 + q_sec_fr^2 <= (rate_a)^2*z*(p_max)) # switching
    JuMP.@constraint(pm.model, p_sec_fr^2 + q_sec_fr^2 <= (rate_a + rate_a_slack)^2) # bound
end

## to

"`p[t_idx]^2 + q[t_idx]^2 <= rate_a^2`"
function constraint_thermal_limit_to_slack(pm::AbstractPowerModel, n::Int, t_idx, rate_a)
    p_to = var(pm, n, :p, t_idx)
    q_to = var(pm, n, :q, t_idx)
    rate_a_slack = var(pm, n, :rate_a_slack, t_idx)

    JuMP.@constraint(pm.model, p_to^2 + q_to^2 <= (rate_a + rate_a_slack)^2, base_name = "$(n)_rating_to_slack_$(t_idx)")
end

"`p[t_idx]^2 + q[t_idx]^2 <= (rate_a * z_branch[i])^2`"
function constraint_thermal_limit_to_slack(pm::AbstractPowerModel, n::Int, i, t_idx, rate_a)
    @assert(:ratea ∈ ref(pm, n, :slack_types))

    p_to = var(pm, n, :p, t_idx)
    q_to = var(pm, n, :q, t_idx)
    
    rate_a_slack = var(pm, n, :rate_a_slack, t_idx)
    JuMP.@constraint(pm.model, p_to^2 + q_to^2 <= (rate_a + rate_a_slack)^2, base_name = "$(n)_rating_to_slack_$(i)")
end

"`p[t_idx]^2 + q[t_idx]^2 <= (rate_a * z_branch[i])^2`"
function constraint_thermal_limit_to_slack_on_off(pm::AbstractPowerModel, n::Int, i, t_idx, rate_a)
    @assert(:ratea ∈ ref(pm, n, :slack_types))
    
    p_to = var(pm, n, :p, t_idx)
    q_to = var(pm, n, :q, t_idx)
    z = var(pm, n, :z_branch, i)

    p_max = ref(pm, n, :p_max)
    
    rate_a_slack = var(pm, n, :rate_a_slack, t_idx)
    JuMP.@constraint(pm.model, p_to^2 + q_to^2 <= (rate_a)^2*z*p_max, base_name = "$(n)_rating_to_switch_$(i)") # switching
    JuMP.@constraint(pm.model, p_to^2 + q_to^2 <= (rate_a + rate_a_slack)^2, base_name = "$(n)_rating_to_bound_$(i)") # bound
end

"`p[t_idx]^2 + q[t_idx]^2 <= rate_a^2`"
function constraint_thermal_limit_to_sec(pm::AbstractPowerModel, n::Int, i::Int, t_idx, rate_a)
    p_to_sec = var(pm, n, :p_sec, t_idx)
    q_to_sec = var(pm, n, :q_sec, t_idx)

    JuMP.@constraint(pm.model, p_to_sec^2 + q_to_sec^2 <= rate_a^2, base_name = "$(n)_rating_to_sec_$(t_idx)")
end



"`p[t_idx]^2 + q[t_idx]^2 <= (rate_a * z_branch[i])^2`"
function constraint_thermal_limit_to_sec_on_off(pm::AbstractPowerModel, n::Int, i, t_idx, rate_a)
    p_sec_to = var(pm, n, :p_sec, t_idx)
    q_sec_to = var(pm, n, :q_sec, t_idx)
    z = var(pm, n, :z_branch, i)

    JuMP.@constraint(pm.model, p_sec_to^2 + q_sec_to^2 <= rate_a^2*z^2)
end

"`p[t_idx]^2 + q[t_idx]^2 <= (rate_a * z_branch[i])^2`"
function constraint_thermal_limit_to_slack_sec(pm::AbstractPowerModel, n::Int, i, t_idx, rate_a)
    @assert(:ratea ∈ ref(pm, n, :slack_types))
    
    p_sec_to = var(pm, n, :p_sec, t_idx)
    q_sec_to = var(pm, n, :q_sec, t_idx)

    (l,c),i,j = t_idx
    
    rate_a_slack = var(pm, n, :rate_a_slack_sec, t_idx)
    JuMP.@constraint(pm.model, p_sec_to^2 + q_sec_to^2 <= (rate_a + rate_a_slack)^2, base_name = "$(n)_rating_to_slack_$(t_idx)") # bound
end

"`p[t_idx]^2 + q[t_idx]^2 <= (rate_a * z_branch[i])^2`"
function constraint_thermal_limit_to_slack_sec_on_off(pm::AbstractPowerModel, n::Int, i, t_idx, rate_a)
    @assert(:ratea ∈ ref(pm, n, :slack_types))

    p_sec_to = var(pm, n, :p_sec, t_idx)
    q_sec_to = var(pm, n, :q_sec, t_idx)
    z = var(pm, n, :z_branch, i)

    p_max = ref(pm, n, :p_max)

    (l,c),i,j = t_idx

    rate_a_slack = var(pm, n, :rate_a_slack_sec, t_idx)
    JuMP.@constraint(pm.model, p_sec_to^2 + q_sec_to^2 <= (rate_a)^2*z*p_max) # switching
    JuMP.@constraint(pm.model, p_sec_to^2 + q_sec_to^2 <= (rate_a + rate_a_slack)^2) # bound
end

### dcline

function constraint_dcline_power_losses_sec(pm::AbstractPowerModel, n::Int, f_bus, t_bus, f_idx, t_idx, loss0, loss1)
    p_sec_fr = var(pm, n, :p_dc_sec, f_idx)
    p_sec_to = var(pm, n, :p_dc_sec, t_idx)

    JuMP.@constraint(pm.model, (1-loss1) * p_sec_fr + (p_sec_to - loss0) == 0)
end

"""
do nothing, most models to not require any model-specific on/off voltage constraints
"""
function constraint_model_voltage_sec(pm::AbstractPowerModel, n::Int, contingencies::Vector{Contingency})
end

#### these constraints are kept here for their name
# Generic thermal limit constraint
"`p[f_idx]^2 + q[f_idx]^2 <= rate_a^2`"
### Only adds constraint name
function constraint_thermal_limit_from(pm::AbstractPowerModel, n::Int, f_idx, rate_a)
    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)

    JuMP.@constraint(pm.model, p_fr^2 + q_fr^2 <= rate_a^2, base_name = "$(n)_rating_from_$(f_idx)")
end

"`p[t_idx]^2 + q[t_idx]^2 <= rate_a^2`"
function constraint_thermal_limit_to(pm::AbstractPowerModel, n::Int, t_idx, rate_a)
    p_to = var(pm, n, :p, t_idx)
    q_to = var(pm, n, :q, t_idx)

    JuMP.@constraint(pm.model, p_to^2 + q_to^2 <= rate_a^2, base_name = "$(n)_rating_to_$(t_idx)")
end

"on/off constraint for generators"
# TODO: the case that it does not have the key :qg should be handled by calling the constraint with the same name, but in apo.jl
function constraint_gen_power_slack_on_off(pm::AbstractPowerModel, n::Int, i::Int, pmin, pmax, qmin, qmax)
    pg = var(pm, n, :pg, i)

    JuMP.set_lower_bound(pg, min(0, ref(pm, n, :gen, i)["pmin"]))
    #JuMP.set_upper_bound(pg, pmax)

    z = var(pm, n, :z_gen, i)
    pmax_slack = get(var(pm, n), :pmax_slack, Dict(i => 0))[i]
    pmin_slack = get(var(pm, n), :pmin_slack, Dict(i => 0))[i]
    qmax_slack = get(var(pm, n), :qmax_slack, Dict(i => 0))[i]
    qmin_slack = get(var(pm, n), :qmin_slack, Dict(i => 0))[i]

    JuMP.@constraint(pm.model, pg <= pmax*PARAMS[:BIG_M_MULT][:PG]  *z, base_name = "$(n)_pg_slack_on_off_ub_$(i)")
    JuMP.@constraint(pm.model, pg >= pmin*PARAMS[:BIG_M_MULT][:PG]  *z, base_name = "$(n)_pg_slack_on_off_lb_$(i)")
    JuMP.@constraint(pm.model, pg <= pmax + pmax_slack, base_name = "$(n)_pg_slack_ub_$(i)")
    JuMP.@constraint(pm.model, pg >= pmin - pmin_slack, base_name = "$(n)_pg_slack_lb_$(i)")

    if haskey(var(pm, n), :qg)
        qg = var(pm, n, :qg, i)
        JuMP.@constraint(pm.model, qg <= qmax*PARAMS[:BIG_M_MULT][:PG]  *z, base_name = "$(n)_qg_slack_on_off_ub_$(i)")
        JuMP.@constraint(pm.model, qg >= qmin*PARAMS[:BIG_M_MULT][:PG]  *z, base_name = "$(n)_qg_slack_on_off_lb_$(i)")
        JuMP.@constraint(pm.model, qg <= qmax + qmax_slack, base_name = "$(n)_qg_slack_ub_$(i)")
        JuMP.@constraint(pm.model, qg >= qmin - qmin_slack, base_name = "$(n)_qg_slack_lb_$(i)")
    end
end

"on/off constraint for generators"
# TODO: the case that it does not have the key :qg should be handled by calling the constraint with the same name, but in apo.jl
function constraint_gen_power_on_off(pm::AbstractPowerModel, n::Int, i::Int, pmin, pmax, qmin, qmax)
    pg = var(pm, n, :pg, i)
    z = var(pm, n, :z_gen, i)

    JuMP.@constraint(pm.model, pg <= pmax*z, base_name = "$(n)_pg_on_off_ub_$(i)")
    JuMP.@constraint(pm.model, pg >= pmin*z, base_name = "$(n)_pg_on_off_lb_$(i)")
    if haskey(var(pm, n), :qg)
        qg = var(pm, n, :qg, i)
        JuMP.@constraint(pm.model, qg <= qmax*z, base_name = "$(n)_qg_on_off_ub_$(i)")
        JuMP.@constraint(pm.model, qg >= qmin*z, base_name = "$(n)_qg_on_off_lb_$(i)")
    end
end

"in case of contingencies, gens cannot be changed"
function constraint_pg_equal_pgsec(pm::AbstractPowerModel; nw::Int=nw_id_default)
    pg      = get(var(pm, nw),   :pg, Dict())
    qg      = get(var(pm, nw),   :qg, Dict())
    pg_sec  = get(var(pm, nw),   :pg_sec, Dict())
    qg_sec  = get(var(pm, nw),   :qg_sec, Dict())

    for i in ids(pm, nw, :gen)
        for (c, contingency) in enumerate(ref(pm, nw, :contingencies))
            if i ∉ get(contingency, :gen, [])
                JuMP.@constraint(pm.model, pg[i] == pg_sec[(i, c)])
                JuMP.@constraint(pm.model, qg[i] == qg_sec[(i, c)])
            end
        end
    end
end

"""
Creates Line Flow constraint for DC Lines (Matpower Formulation)

```
p_fr + p_to == loss0 + p_fr * loss1
```
"""
### Only adds constraint name
function constraint_dcline_power_losses(pm::AbstractPowerModel, n::Int, f_bus, t_bus, f_idx, t_idx, loss0, loss1)
    p_fr = var(pm, n, :p_dc, f_idx)
    p_to = var(pm, n, :p_dc, t_idx)

    JuMP.@constraint(pm.model, (1-loss1) * p_fr + (p_to - loss0) == 0, base_name = "$(n)_dcline_power_losses_$(f_idx)_$(t_idx)")
end