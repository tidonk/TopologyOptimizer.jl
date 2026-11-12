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

""
function variable_bus_voltage_sec_on_off(pm::AbstractWRModel, contingencies::Vector{Contingency}; kwargs...)
    variable_bus_voltage_magnitude_sqr_sec(pm, contingencies; kwargs...)

    variable_branch_voltage_magnitude_fr_sqr_sec_on_off(pm, contingencies; kwargs...)
    variable_branch_voltage_magnitude_to_sqr_sec_on_off(pm, contingencies; kwargs...)
    variable_branch_voltage_product_sec_on_off(pm, contingencies; kwargs...)
end

""
function constraint_model_voltage_sec_on_off(pm::AbstractWRModel, n::Int, contingencies::Vector{Contingency})
    w_sec  = var(pm, n, :w_sec)
    wr_sec = var(pm, n, :wr_sec)
    wi_sec = var(pm, n, :wi_sec)
    z  = var(pm, n, :z_branch)
    
    w_fr_sec = var(pm, n, :w_fr_sec)
    w_to_sec = var(pm, n, :w_to_sec)

    constraint_voltage_magnitude_sqr_from_sec_on_off(pm, n, contingencies)
    constraint_voltage_magnitude_sqr_to_sec_on_off(pm, n, contingencies)
    constraint_voltage_product_sec_on_off(pm, n, contingencies)

    
    for (c, contingency) in enumerate(contingencies)
        for (l,i,j) in ref(pm, n, :arcs_from)
            if l ∉ get(contingency, :branch, [])
                _IM.relaxation_complex_product_on_off(pm.model, w_sec[(i,c)], w_sec[(j,c)], wr_sec[(l,c)], wi_sec[(l,c)], z[l])
                _IM.relaxation_equality_on_off(pm.model, w_sec[(i,c)], w_fr_sec[(l,c)], z[l])
                _IM.relaxation_equality_on_off(pm.model, w_sec[(j,c)], w_to_sec[(l,c)], z[l])
            end
        end
    end
end

""
function constraint_voltage_magnitude_from_sec_on_off(pm::AbstractWRModel, n::Int, contingencies::Vector{Contingency})
    buses = ref(pm, n, :bus)
    branches = ref(pm, n, :branch)

    vm_fr_sec = var(pm, n, :vm_fr_sec)
    z = var(pm, n, :z_branch)

    for (c, contingency) in enumerate(contingencies)
        for (i, branch) in ref(pm, n, :branch)
            if i ∉ get(contingency, :branch, [])

                JuMP.@constraint(pm.model, vm_fr_sec[(i,c)] <= z[i]*buses[branch["f_bus"]]["vmax"])
                JuMP.@constraint(pm.model, vm_fr_sec[(i,c)] >= z[i]*buses[branch["f_bus"]]["vmin"])
            end
        end
    end
end

""
function constraint_voltage_magnitude_to_sec_on_off(pm::AbstractWRModel, n::Int, contingencies::Vector{Contingency})
    buses = ref(pm, n, :bus)
    branches = ref(pm, n, :branch)

    vm_to_sec = var(pm, n, :vm_to_sec)
    z = var(pm, n, :z_branch)


    for (c, contingency) in enumerate(contingencies)
        for (i, branch) in ref(pm, n, :branch)
            if i ∉ get(contingency, :branch, [])

                JuMP.@constraint(pm.model, vm_to_sec[(i,c)] <= z[i]*buses[branch["t_bus"]]["vmax"])
                JuMP.@constraint(pm.model, vm_to_sec[(i,c)] >= z[i]*buses[branch["t_bus"]]["vmin"])
            end
        end
    end
end

""
function constraint_voltage_magnitude_sqr_from_sec_on_off(pm::AbstractWRModel, n::Int, contingencies::Vector{Contingency})
    buses = ref(pm, n, :bus)
    branches = ref(pm, n, :branch)

    w_fr_sec = var(pm, n, :w_fr_sec)
    z = var(pm, n, :z_branch)

    for (c, contingency) in enumerate(contingencies)
        for (i, branch) in ref(pm, n, :branch)
            if i ∉ get(contingency, :branch, [])

                JuMP.@constraint(pm.model, w_fr_sec[(i,c)] <= z[i]*buses[branch["f_bus"]]["vmax"]^2)
                JuMP.@constraint(pm.model, w_fr_sec[(i,c)] >= z[i]*buses[branch["f_bus"]]["vmin"]^2)
            end
        end
    end
end

""
function constraint_voltage_magnitude_sqr_to_sec_on_off(pm::AbstractWRModel, n::Int, contingencies::Vector{Contingency})
    buses = ref(pm, n, :bus)
    branches = ref(pm, n, :branch)

    w_to_sec = var(pm, n, :w_to_sec)
    z = var(pm, n, :z_branch)

    for (c, contingency) in enumerate(contingencies)
        for (i, branch) in ref(pm, n, :branch)
            if i ∉ get(contingency, :branch, [])

                JuMP.@constraint(pm.model, w_to_sec[(i,c)] <= z[i]*buses[branch["t_bus"]]["vmax"]^2)
                JuMP.@constraint(pm.model, w_to_sec[(i,c)] >= z[i]*buses[branch["t_bus"]]["vmin"]^2)
            end
        end
    end
end

""
function constraint_voltage_product_sec_on_off(pm::AbstractWRModel, n::Int, contingencies)
    wr_min, wr_max, wi_min, wi_max = ref_calc_voltage_product_bounds(ref(pm, n, :buspairs))

    bi_bp = Dict((i, (b["f_bus"], b["t_bus"])) for (i,b) in ref(pm, n, :branch))

    wr_sec = var(pm, n, :wr_sec)
    wi_sec = var(pm, n, :wi_sec)
    z  = var(pm, n, :z_branch)

    for (c, contingency) in enumerate(contingencies)
        for b in ids(pm, n, :branch)
            if b ∉ get(contingency, :branch, [])
                JuMP.@constraint(pm.model, wr_sec[(b,c)] <= z[b]*wr_max[bi_bp[b]])
                JuMP.@constraint(pm.model, wr_sec[(b,c)] >= z[b]*wr_min[bi_bp[b]])
                JuMP.@constraint(pm.model, wi_sec[(b,c)] <= z[b]*wi_max[bi_bp[b]])
                JuMP.@constraint(pm.model, wi_sec[(b,c)] >= z[b]*wi_min[bi_bp[b]])
            end
        end
    end
end

"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[f_idx] ==        g/tm*w_fr[i] + (-g*tr+b*ti)/tm*(wr[i]) + (-b*tr-g*ti)/tm*(wi[i])
q[f_idx] == -(b+c/2)/tm*w_fr[i] - (-b*tr-g*ti)/tm*(wr[i]) + (-g*tr+b*ti)/tm*(wi[i])
```
"""
function constraint_ohms_yt_from_sec(pm::AbstractWRModel, n::Int, i, c::ContingencyID, contingency::Contingency, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_fr_sec = var(pm, n, :p_sec, f_idx)
    q_fr_sec = var(pm, n, :q_sec, f_idx)
    w_fr_sec = var(pm, n, :w_fr_sec, (i,c))
    wr_sec   = var(pm, n, :wr_sec, (i,c))
    wi_sec   = var(pm, n, :wi_sec, (i,c))

    JuMP.@constraint(pm.model, p_fr_sec ==  (g+g_fr)/tm^2*w_fr_sec + (-g*tr+b*ti)/tm^2*wr_sec + (-b*tr-g*ti)/tm^2*wi_sec, base_name="$(n)_C$(c)_ohm_p_fr_$(i)" )
    JuMP.@constraint(pm.model, q_fr_sec == -(b+b_fr)/tm^2*w_fr_sec - (-b*tr-g*ti)/tm^2*wr_sec + (-g*tr+b*ti)/tm^2*wi_sec, base_name="$(n)_C$(c)_ohm_q_fr_$(i)" )
end

"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[f_idx] ==        g/tm*w_fr[i] + (-g*tr+b*ti)/tm*(wr[i]) + (-b*tr-g*ti)/tm*(wi[i])
q[f_idx] == -(b+c/2)/tm*w_fr[i] - (-b*tr-g*ti)/tm*(wr[i]) + (-g*tr+b*ti)/tm*(wi[i])
```
"""
function constraint_ohms_yt_from_sec_on_off(pm::AbstractWRModel, n::Int, i, c::ContingencyID, contingency::Contingency, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max)
    p_fr_sec = var(pm, n, :p_sec, f_idx)
    q_fr_sec = var(pm, n, :q_sec, f_idx)
    w_fr_sec = var(pm, n, :w_fr_sec, (i,c))
    wr_sec   = var(pm, n, :wr_sec, (i,c))
    wi_sec   = var(pm, n, :wi_sec, (i,c))

    JuMP.@constraint(pm.model, p_fr_sec ==  (g+g_fr)/tm^2*w_fr_sec + (-g*tr+b*ti)/tm^2*wr_sec + (-b*tr-g*ti)/tm^2*wi_sec )
    JuMP.@constraint(pm.model, q_fr_sec == -(b+b_fr)/tm^2*w_fr_sec - (-b*tr-g*ti)/tm^2*wr_sec + (-g*tr+b*ti)/tm^2*wi_sec )
end

"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[t_idx] ==        g*w_to[i] + (-g*tr-b*ti)/tm*(wr[i]) + (-b*tr+g*ti)/tm*(-wi[i])
q[t_idx] == -(b+c/2)*w_to[i] - (-b*tr+g*ti)/tm*(wr[i]) + (-g*tr-b*ti)/tm*(-wi[i])
```
"""
function constraint_ohms_yt_to_sec(pm::AbstractWRModel, n::Int, i, c::ContingencyID, contingency::Contingency, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    p_to_sec = var(pm, n, :p_sec, t_idx)
    q_to_sec = var(pm, n, :q_sec, t_idx)
    w_to_sec = var(pm, n, :w_to_sec, (i,c))
    wr_sec = var(pm, n, :wr_sec, (i,c))
    wi_sec = var(pm, n, :wi_sec, (i,c))

    JuMP.@constraint(pm.model, p_to_sec ==  (g+g_to)*w_to_sec + (-g*tr-b*ti)/tm^2*wr_sec + (-b*tr+g*ti)/tm^2*-wi_sec, base_name="$(n)_C$(c)_ohm_p_to_$(i)" )
    JuMP.@constraint(pm.model, q_to_sec == -(b+b_to)*w_to_sec - (-b*tr+g*ti)/tm^2*wr_sec + (-g*tr-b*ti)/tm^2*-wi_sec, base_name="$(n)_C$(c)_ohm_q_to_$(i)" )
end

"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[t_idx] ==        g*w_to[i] + (-g*tr-b*ti)/tm*(wr[i]) + (-b*tr+g*ti)/tm*(-wi[i])
q[t_idx] == -(b+c/2)*w_to[i] - (-b*tr+g*ti)/tm*(wr[i]) + (-g*tr-b*ti)/tm*(-wi[i])
```
"""
function constraint_ohms_yt_to_sec_on_off(pm::AbstractWRModel, n::Int, i, c::ContingencyID, contingency::Contingency, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max)
    p_to_sec = var(pm, n, :p_sec, t_idx)
    q_to_sec = var(pm, n, :q_sec, t_idx)
    w_to_sec = var(pm, n, :w_to_sec, (i,c))
    wr_sec = var(pm, n, :wr_sec, (i,c))
    wi_sec = var(pm, n, :wi_sec, (i,c))

    JuMP.@constraint(pm.model, p_to_sec ==  (g+g_to)*w_to_sec + (-g*tr-b*ti)/tm^2*wr_sec + (-b*tr+g*ti)/tm^2*-wi_sec )
    JuMP.@constraint(pm.model, q_to_sec == -(b+b_to)*w_to_sec - (-b*tr+g*ti)/tm^2*wr_sec + (-g*tr-b*ti)/tm^2*-wi_sec )
end

"`angmin*wr[i] <= wi[i] <= angmax*wr[i]`"
function constraint_voltage_angle_difference_slack(pm::AbstractWRModel, n::Int, f_idx, angmin, angmax, vad_min, vad_max)
    @error("Voltage angle difference slack variables not implemented for this formulation.")

    i, f_bus, t_bus = f_idx
    wr = var(pm, n, :wr, i)
    wi = var(pm, n, :wi, i)

    JuMP.@constraint(pm.model, wi <= tan(angmax)*wr)
    JuMP.@constraint(pm.model, wi >= tan(angmin)*wr)
end

"`angmin*wr[i] <= wi[i] <= angmax*wr[i]`"
function constraint_voltage_angle_difference_slack_on_off(pm::AbstractWRModel, n::Int, f_idx, angmin, angmax, vad_min, vad_max)
    @error("Voltage angle difference slack variables not implemented for this formulation.")

    i, f_bus, t_bus = f_idx
    wr = var(pm, n, :wr, i)
    wi = var(pm, n, :wi, i)

    JuMP.@constraint(pm.model, wi <= tan(angmax)*wr)
    JuMP.@constraint(pm.model, wi >= tan(angmin)*wr)
end

"`angmin*wr[i] <= wi[i] <= angmax*wr[i]`"
function constraint_voltage_angle_difference_sec(pm::AbstractWRModel, n::Int, c::ContingencyID, contingency::Contingency, f_idx, angmin, angmax)
    (i,c), f_bus, t_bus = f_idx
    wr_sec = var(pm, n, :wr_sec, (i,c))
    wi_sec = var(pm, n, :wi_sec, (i,c))

    JuMP.@constraint(pm.model, wi_sec <= tan(angmax)*wr_sec, base_name="$(n)_C$(c)_vad_ub_$(i)")
    JuMP.@constraint(pm.model, wi_sec >= tan(angmin)*wr_sec, base_name="$(n)_C$(c)_vad_lb_$(i)")
end

"`angmin*wr[i] <= wi[i] <= angmax*wr[i]`"
function constraint_voltage_angle_difference_sec_on_off(pm::AbstractWRModel, n::Int, c::ContingencyID, contingency::Contingency, f_idx, angmin, angmax, vad_min, vad_max)
    (i,c), f_bus, t_bus = f_idx
    wr_sec = var(pm, n, :wr_sec, (i,c))
    wi_sec = var(pm, n, :wi_sec, (i,c))

    JuMP.@constraint(pm.model, wi_sec <= tan(angmax)*wr_sec)
    JuMP.@constraint(pm.model, wi_sec >= tan(angmin)*wr_sec)
end

"`angmin*wr[i] <= wi[i] <= angmax*wr[i]`"
function constraint_voltage_angle_difference_slack_sec(pm::AbstractWRModel, n::Int, c::ContingencyID, contingency::Contingency, f_idx, angmin, angmax, vad_min, vad_max)
    @error("Voltage angle difference slack variables not implemented for this formulation.")
    
    (i,c), f_bus, t_bus = f_idx
    wr_sec = var(pm, n, :wr_sec, (i,c))
    wi_sec = var(pm, n, :wi_sec, (i,c))

    JuMP.@constraint(pm.model, wi_sec <= tan(angmax)*wr_sec)
    JuMP.@constraint(pm.model, wi_sec >= tan(angmin)*wr_sec)
end

"`angmin*wr[i] <= wi[i] <= angmax*wr[i]`"
function constraint_voltage_angle_difference_slack_sec_on_off(pm::AbstractWRModel, n::Int, c::ContingencyID, contingency::Contingency, f_idx, angmin, angmax, vad_min, vad_max)
    @error("Voltage angle difference slack variables not implemented for this formulation.")
    
    (i,c), f_bus, t_bus = f_idx
    wr_sec = var(pm, n, :wr_sec, (i,c))
    wi_sec = var(pm, n, :wi_sec, (i,c))

    JuMP.@constraint(pm.model, wi_sec <= tan(angmax)*wr_sec)
    JuMP.@constraint(pm.model, wi_sec >= tan(angmin)*wr_sec)
end

"`t[ref_bus] == 0`"
function constraint_theta_ref_sec(pm::AbstractQCWRModel, n::Int, i::Int, c::ContingencyID, contingency::Dict)
    JuMP.@constraint(pm.model, var(pm, n, :va_sec)[(i,c)] == 0)
end

""
function variable_bus_voltage_sec_on_off(pm::AbstractQCWRModel, contingencies::Vector{Contingency}; kwargs...)
    variable_bus_voltage_angle_sec(pm, contingencies; kwargs...)
    variable_bus_voltage_magnitude_sec(pm, contingencies; kwargs...)
    variable_branch_voltage_magnitude_fr_sec_on_off(pm, contingencies; kwargs...)
    variable_branch_voltage_magnitude_to_sec_on_off(pm, contingencies; kwargs...)

    variable_bus_voltage_magnitude_sqr_sec(pm, contingencies; kwargs...)
    variable_branch_voltage_magnitude_fr_sqr_sec_on_off(pm, contingencies; kwargs...)
    variable_branch_voltage_magnitude_to_sqr_sec_on_off(pm, contingencies; kwargs...)

    variable_branch_voltage_product_sec_on_off(pm, contingencies; kwargs...)

    variable_branch_voltage_product_angle_sec_on_off(pm, contingencies; kwargs...)
    variable_buspair_voltage_product_magnitude_sec_on_off(pm, contingencies; kwargs...)
    variable_branch_cosine_sec_on_off(pm, contingencies; kwargs...)
    variable_branch_sine_sec_on_off(pm, contingencies; kwargs...)
    variable_branch_current_magnitude_sqr_sec_on_off(pm, contingencies; kwargs...) # includes 0, but needs new indexes
end

""
function variable_branch_voltage_product_angle_sec_on_off(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, report::Bool=true)
    sec_branches = [(i,c) for i in ids(pm, nw, :branch) for (c, contingency) in enumerate(contingencies) if i ∉ get(contingency, :branch, [])]
    
    td_sec = var(pm, nw)[:td_sec] = JuMP.@variable(pm.model,
        [(i,c) in sec_branches], base_name="$(nw)_td_sec",
        lower_bound = min(0, ref(pm, nw, :branch, i, "angmin")),
        upper_bound = max(0, ref(pm, nw, :branch, i, "angmax")),
        start = comp_start_value(ref(pm, nw, :branch, i), "td_start")
    )

    report && sol_component_value(pm, nw, :branch, :td_sec, sec_branches, td_sec)
end

""
function variable_buspair_voltage_product_magnitude_sec_on_off(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, report::Bool=true)
    branches = ref(pm, nw, :branch)
    buses = ref(pm, nw, :bus)
    sec_branches = [(i,c) for i in ids(pm, nw, :branch) for (c, contingency) in enumerate(contingencies) if i ∉ get(contingency, :branch, [])]

    vv_min = Dict((l, buses[branch["f_bus"]]["vmin"]*buses[branch["t_bus"]]["vmin"]) for (l, branch) in branches)
    vv_max = Dict((l, buses[branch["f_bus"]]["vmax"]*buses[branch["t_bus"]]["vmax"]) for (l, branch) in branches)

    vv_sec = var(pm, nw)[:vv_sec] = JuMP.@variable(pm.model,
        [(l,c) in sec_branches], base_name="$(nw)_vv_sec",
        lower_bound = min(0, vv_min[l]),
        upper_bound = max(0, vv_max[l]),
        start = comp_start_value(ref(pm, nw, :branch, l), "vv_start", 1.0)
    )

    report && sol_component_value(pm, nw, :branch, :vv_sec, sec_branches, vv_sec)
end

""
function variable_branch_cosine_sec_on_off(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, report::Bool=true)
    sec_branches = [(i,c) for i in ids(pm, nw, :branch) for (c, contingency) in enumerate(contingencies) if i ∉ get(contingency, :branch, [])]

    cos_min = Dict(((i,c), -Inf) for (i,c) in sec_branches)
    cos_max = Dict(((i,c),  Inf) for (i,c) in sec_branches)

    for (i,c) in sec_branches
        branch = ref(pm,nw,:branch)[i]
        angmin = branch["angmin"]
        angmax = branch["angmax"]
        if angmin >= 0
            cos_max[(i,c)] = cos(angmin)
            cos_min[(i,c)] = cos(angmax)
        end
        if angmax <= 0
            cos_max[(i,c)] = cos(angmax)
            cos_min[(i,c)] = cos(angmin)
        end
        if angmin < 0 && angmax > 0
            cos_max[(i,c)] = 1.0
            cos_min[(i,c)] = min(cos(angmin), cos(angmax))
        end
    end

    cs_sec = var(pm, nw)[:cs_sec] = JuMP.@variable(pm.model,
        [(i,c) in sec_branches], base_name="$(nw)_cs_sec",
        lower_bound = min(0, cos_min[(i,c)]),
        upper_bound = max(0, cos_max[(i,c)]),
        start = comp_start_value(ref(pm, nw, :branch, i), "cs_start", 1.0)
    )

    report && sol_component_value(pm, nw, :branch, :cs_sec, sec_branches, cs_sec)
end

""
function variable_branch_sine_sec_on_off(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, report::Bool=true)
    sec_branches = [(i,c) for i in ids(pm, nw, :branch) for (c, contingency) in enumerate(contingencies) if i ∉ get(contingency, :branch, [])]
    
    si_sec = var(pm, nw)[:si_sec] = JuMP.@variable(pm.model,
        [(l,c) in sec_branches], base_name="$(nw)_si_sec",
        lower_bound = min(0, sin(ref(pm, nw, :branch, l, "angmin"))),
        upper_bound = max(0, sin(ref(pm, nw, :branch, l, "angmax"))),
        start = comp_start_value(ref(pm, nw, :branch, l), "si_start")
    )

    report && sol_component_value(pm, nw, :branch, :si_sec, sec_branches, si_sec)
end

""
function variable_branch_current_magnitude_sqr_sec_on_off(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, report::Bool=true)
    ccm_min = Dict((l, 0) for l in ids(pm, nw, :branch))

    branches = ref(pm, nw, :branch)
    ccm_max = Dict()
    for (l, branch) in branches
        vm_fr_min = ref(pm, nw, :bus, branch["f_bus"], "vmin")
        ccm_max[l] = ((branch["rate_a"]*branch["tap"])/vm_fr_min)^2
    end

    sec_branches = [(i,c) for i in ids(pm, nw, :branch) for (c, contingency) in enumerate(contingencies) if i ∉ get(contingency, :branch, [])]

    ccm_sec = var(pm, nw)[:ccm_sec] = JuMP.@variable(pm.model,
        [(l,c) in sec_branches], base_name="$(nw)_ccm_sec",
        lower_bound = ccm_min[l],
        upper_bound = ccm_max[l],
        start = comp_start_value(ref(pm, nw, :branch, l), "ccm_start")
    )

    report && sol_component_value(pm, nw, :branch, :ccm_sec, sec_branches, ccm_sec)
end

""
function constraint_model_voltage_sec_on_off(pm::AbstractQCWRModel, n::Int, contingencies::Vector{Contingency})
    v_sec = var(pm, n, :vm_sec)
    t_sec = var(pm, n, :va_sec)
    vm_fr_sec = var(pm, n, :vm_fr_sec)
    vm_to_sec = var(pm, n, :vm_to_sec)

    td_sec = var(pm, n, :td_sec)
    si_sec = var(pm, n, :si_sec)
    cs_sec = var(pm, n, :cs_sec)
    vv_sec = var(pm, n, :vv_sec)

    w_sec = var(pm, n, :w_sec)
    w_fr_sec = var(pm, n, :w_fr_sec)
    w_to_sec = var(pm, n, :w_to_sec)

    wr_sec = var(pm, n, :wr_sec)
    wi_sec = var(pm, n, :wi_sec)

    z = var(pm, n, :z_branch)

    td_lb = ref(pm, n, :off_angmin)
    td_ub = ref(pm, n, :off_angmax)
    td_max = max(abs(td_lb), abs(td_ub))

    for (c, contingency) in enumerate(contingencies)
        for i in ids(pm, n, :bus)
            if i ∉ get(contingency, :bus, [])
                _IM.relaxation_sqr(pm.model, v_sec[(i,c)], w_sec[(i,c)])
            end
        end
    end

    constraint_voltage_magnitude_from_sec_on_off(pm, n, contingencies) # bounds on vm_fr
    constraint_voltage_magnitude_to_sec_on_off(pm, n, contingencies) # bounds on vm_to
    constraint_voltage_magnitude_sqr_from_sec_on_off(pm, n, contingencies) # bounds on w_fr
    constraint_voltage_magnitude_sqr_to_sec_on_off(pm, n, contingencies) # bounds on w_to
    constraint_voltage_product_sec_on_off(pm, n, contingencies) # bounds on wr, wi

    for (c, contingency) in enumerate(contingencies)
        for (l, branch) in ref(pm, n, :branch)
            if l ∉ get(contingency, :branch, [])
                i = branch["f_bus"]
                j = branch["t_bus"]

                JuMP.@constraint(pm.model, t_sec[(i,c)] - t_sec[(j,c)] >= td_sec[(l,c)] + td_lb*(1-z[l]))
                JuMP.@constraint(pm.model, t_sec[(i,c)] - t_sec[(j,c)] <= td_sec[(l,c)] + td_ub*(1-z[l]))

                relaxation_sin_on_off(pm.model, td_sec[(l,c)], si_sec[(l,c)], z[l], td_max)
                relaxation_cos_on_off(pm.model, td_sec[(l,c)], cs_sec[(l,c)], z[l], td_max)
                _IM.relaxation_product_on_off(pm.model, vm_fr_sec[(l,c)], vm_to_sec[(l,c)], vv_sec[(l,c)], z[l])
                _IM.relaxation_product_on_off(pm.model, vv_sec[(l,c)], cs_sec[(l,c)], wr_sec[(l,c)], z[l])
                _IM.relaxation_product_on_off(pm.model, vv_sec[(l,c)], si_sec[(l,c)], wi_sec[(l,c)], z[l])

                # this constraint is redudant and useful for debugging
                #_IM.relaxation_complex_product(pm.model, w[_seci], w[j], wr[l], wi[l])

                #cs4 = _IM.relaxation_complex_product_on_off(pm.model, w[i], w[j], wr[l], wi[l], z[l])
                _IM.relaxation_equality_on_off(pm.model, v_sec[(i,c)], vm_fr_sec[(l,c)], z[l])
                _IM.relaxation_equality_on_off(pm.model, v_sec[(j,c)], vm_to_sec[(l,c)], z[l])
                _IM.relaxation_equality_on_off(pm.model, w_sec[(i,c)], w_fr_sec[(l,c)], z[l])
                _IM.relaxation_equality_on_off(pm.model, w_sec[(j,c)], w_to_sec[(l,c)], z[l])

                # to prevent this constraint from being posted on multiple parallel branchs
                # TODO needs on/off variant
                constraint_power_magnitude_sqr_sec_on_off(pm, l, nw=n, c, contingency)
                constraint_power_magnitude_link_sec_on_off(pm, l, nw=n, c, contingency) # different index set
            end
        end
    end
end

"`p[arc_from]^2 + q[arc_from]^2 <= w[f_bus]/tm*cm[i]`"
function constraint_power_magnitude_sqr_sec_on_off(pm::AbstractQCWRModel, n::Int, i, c, f_bus, arc_from, tm)
    w_sec    = var(pm, n, :w_sec, (f_bus,c))
    p_fr_sec = var(pm, n, :p_sec, arc_from)
    q_fr_sec = var(pm, n, :q_sec, arc_from)
    ccm_sec   = var(pm, n, :ccm_sec, (i,c))
    z = var(pm, n, :z_branch, i)

    # TODO see if there is a way to leverage relaxation_complex_product_on_off here
    w_lb, w_ub = _IM.variable_domain(w_sec)
    ccm_lb, ccm_ub = _IM.variable_domain(ccm_sec)
    z_lb, z_ub = _IM.variable_domain(z)

    JuMP.@constraint(pm.model, p_fr_sec^2 + q_fr_sec^2 <= w_sec*ccm_sec*z_ub/tm^2)
    JuMP.@constraint(pm.model, p_fr_sec^2 + q_fr_sec^2 <= w_ub*ccm_sec*z/tm^2)
    JuMP.@constraint(pm.model, p_fr_sec^2 + q_fr_sec^2 <= w_sec*ccm_ub*z/tm^2)
end

"`ccm[f_bus,t_bus] == (g^2 + b^2)*(w[f_bus]/tm + w[t_bus] - 2*(tr*wr[f_bus,t_bus] + ti*wi[f_bus,t_bus])/tm) - c*q[f_idx] - ((c/2)/tm)^2*w[f_bus]`"
function constraint_power_magnitude_link_sec_on_off(pm::AbstractQCWRModel, n::Int, i, c, arc_from, g, b, g_fr, b_fr, g_to, b_to, tr, ti, tm)
    w_fr_sec = var(pm, n, :w_fr_sec, (i,c))
    w_to_sec = var(pm, n, :w_to_sec, (i,c))
    p_fr_sec = var(pm, n, :p_sec, arc_from)
    q_fr_sec = var(pm, n, :q_sec, arc_from)
    wr_sec   = var(pm, n, :wr_sec, (i,c))
    wi_sec   = var(pm, n, :wi_sec, (i,c))
    ccm_sec   = var(pm, n, :ccm_sec, (i,c))

    ym_sh_sqr = g_fr^2 + b_fr^2

    JuMP.@constraint(pm.model, ccm_sec == (g^2 + b^2)*(w_fr_sec/tm^2 + w_to_sec - 2*(tr*wr_sec + ti*wi_sec)/tm^2) - ym_sh_sqr*(w_fr_sec/tm^2) + 2*(g_fr*p_fr_sec - b_fr*q_fr_sec))
end