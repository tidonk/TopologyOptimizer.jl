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

"`angmin <= z_branch[i]*(t[f_bus] - t[t_bus]) <= angmax` + slack linking"
function constraint_voltage_angle_difference_slack(pm::AbstractDCPModel, n::Int, f_idx, angmin, angmax, vad_min, vad_max)
    i, f_bus, t_bus = f_idx
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)
    
    vad_slack_lower = get(var(pm, n), :vad_slack_lower, Dict(f_idx => 0))[f_idx]
    vad_slack_upper = get(var(pm, n), :vad_slack_upper, Dict(f_idx => 0))[f_idx]

    # switched on?
    JuMP.@constraint(pm.model, va_fr - va_to <= angmax + vad_slack_upper, base_name = "$(n)_vad_slack_ub_$(i)")
    JuMP.@constraint(pm.model, va_fr - va_to >= angmin - vad_slack_lower, base_name = "$(n)_vad_slack_lb_$(i)")
end

# Modified from PowerModels.jl (https://github.com/lanl-ansi/PowerModels.jl), see LICENSE
"`angmin*z_branch[i] + vad_min*(1-z_branch[i]) <= t[f_bus] - t[t_bus] <= angmax*z_branch[i] + vad_max*(1-z_branch[i])`"
function constraint_voltage_angle_difference_on_off(pm::AbstractDCPModel, n::Int, f_idx, angmin, angmax, vad_min, vad_max)
    i, f_bus, t_bus = f_idx

    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)
    z = var(pm, n, :z_branch, i)

    JuMP.@constraint(pm.model, va_fr - va_to <= angmax*z + vad_max*(1-z), base_name = "$(n)_vad_on_off_ub_$(i)")
    JuMP.@constraint(pm.model, va_fr - va_to >= angmin*z + vad_min*(1-z), base_name = "$(n)_vad_on_off_lb_$(i)")
end

"`angmin <= z_branch[i]*(t[f_bus] - t[t_bus]) <= angmax` + slack linking"
function constraint_voltage_angle_difference_slack_on_off(pm::AbstractDCPModel, n::Int, f_idx, angmin, angmax, vad_min, vad_max)
    i, f_bus, t_bus = f_idx
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)
    z = var(pm, n, :z_branch, i)
    
    vad_slack_lower = get(var(pm, n), :vad_slack_lower, Dict(f_idx => 0))[f_idx]
    vad_slack_upper = get(var(pm, n), :vad_slack_upper, Dict(f_idx => 0))[f_idx]

    # switched on?
    @assert(angmin < 0)
    JuMP.@constraint(pm.model, va_fr - va_to <= angmax + vad_slack_upper + (1-z) * (angmax * PARAMS[:BIG_M_MULT][:VAD]  ), base_name = "$(n)_vad_slack_on_off_ub_$(i)")
    JuMP.@constraint(pm.model, va_fr - va_to >= angmin - vad_slack_lower + (1-z) * (angmin * PARAMS[:BIG_M_MULT][:VAD]  ), base_name = "$(n)_vad_slack_on_off_lb_$(i)")
end

"`angmin*z_branch[i] + vad_min*(1-z_branch[i]) <= t[f_bus] - t[t_bus] <= angmax*z_branch[i] + vad_max*(1-z_branch[i])`"
function constraint_voltage_angle_difference_sec(pm::AbstractDCPModel, n::Int, c::ContingencyID, contingency::Dict, f_idx, angmin, angmax)
    (i,c), f_bus, t_bus = f_idx

    va_sec_fr = var(pm, n, :va_sec, (f_bus,c))
    va_sec_to = var(pm, n, :va_sec, (t_bus,c))

    # vad_max / vad_min are big M constants
    JuMP.@constraint(pm.model, va_sec_fr - va_sec_to <= angmax, base_name = "$(n)_C$(c)_vad_ub_$(i)")
    JuMP.@constraint(pm.model, va_sec_fr - va_sec_to >= angmin, base_name = "$(n)_C$(c)_vad_lb_$(i)")
end

"`angmin*z_branch[i] + vad_min*(1-z_branch[i]) <= t[f_bus] - t[t_bus] <= angmax*z_branch[i] + vad_max*(1-z_branch[i])`"
function constraint_voltage_angle_difference_sec_on_off(pm::AbstractDCPModel, n::Int, c::ContingencyID, contingency::Dict, f_idx, angmin, angmax, vad_min, vad_max)
    (i,c), f_bus, t_bus = f_idx

    va_sec_fr = var(pm, n, :va_sec, (f_bus,c))
    va_sec_to = var(pm, n, :va_sec, (t_bus,c))
    z = var(pm, n, :z_branch, i)

    # vad_max / vad_min are big M constants
    JuMP.@constraint(pm.model, va_sec_fr - va_sec_to <= angmax*z + vad_max*(1-z), base_name = "$(n)_C$(c)_vad_on_off_ub_$(i)")
    JuMP.@constraint(pm.model, va_sec_fr - va_sec_to >= angmin*z + vad_min*(1-z), base_name = "$(n)_C$(c)_vad_on_off_lb_$(i)")
end

"`angmin <= z_branch[i]*(t[f_bus] - t[t_bus]) <= angmax` +  slack linking"
function constraint_voltage_angle_difference_slack_sec(pm::AbstractDCPModel, n::Int, c::ContingencyID, contingency::Contingency, f_idx, angmin, angmax, vad_min, vad_max)
    (i, c), f_bus, t_bus = f_idx

    va_sec_fr = var(pm, n, :va_sec, (f_bus,c))
    va_sec_to = var(pm, n, :va_sec, (t_bus,c))

    vad_slack_lower = get(var(pm, n), :vad_slack_lower, Dict((i, f_bus, t_bus) => 0))[(i, f_bus, t_bus)]
    vad_slack_upper = get(var(pm, n), :vad_slack_upper, Dict((i, f_bus, t_bus) => 0))[(i, f_bus, t_bus)]

     # slack linking
    JuMP.@constraint(pm.model, va_sec_fr - va_sec_to <= angmax + vad_slack_upper, base_name = "$(n)_C$(c)_vad_slack_ub_$(i)")
    JuMP.@constraint(pm.model, va_sec_fr - va_sec_to >= angmin - vad_slack_lower, base_name = "$(n)_C$(c)_vad_slack_lb_$(i)")
end

"`angmin <= z_branch[i]*(t[f_bus] - t[t_bus]) <= angmax` +  slack linking"
function constraint_voltage_angle_difference_slack_sec_on_off(pm::AbstractDCPModel, n::Int, c::ContingencyID, contingency::Contingency, f_idx, angmin, angmax, vad_min, vad_max)
    (i, c), f_bus, t_bus = f_idx

    va_sec_fr = var(pm, n, :va_sec, (f_bus,c))
    va_sec_to = var(pm, n, :va_sec, (t_bus,c))
    z = var(pm, n, :z_branch, i)

    vad_slack_lower = get(var(pm, n), :vad_slack_lower, Dict((i, f_bus, t_bus) => 0))[(i, f_bus, t_bus)]
    vad_slack_upper = get(var(pm, n), :vad_slack_upper, Dict((i, f_bus, t_bus) => 0))[(i, f_bus, t_bus)]

     # slack linking
    @assert(angmin < 0)
    JuMP.@constraint(pm.model, va_sec_fr - va_sec_to <= angmax + vad_slack_upper + (1-z) * (angmax * PARAMS[:BIG_M_MULT][:VAD]  ), base_name = "$(n)_C$(c)_vad_slack_on_off_ub_$(i)")
    JuMP.@constraint(pm.model, va_sec_fr - va_sec_to >= angmin - vad_slack_lower + (1-z) * (angmin * PARAMS[:BIG_M_MULT][:VAD]  ), base_name = "$(n)_C$(c)_vad_slack_on_off_lb_$(i)")
end

"nothing to do, this model is symmetric"
function constraint_ohms_yt_to_sec(pm::AbstractDCPModel, n::Int, i, c, contingency, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
end

"nothing to do, this model is symmetric"
function constraint_ohms_yt_to_sec_on_off(pm::AbstractDCPModel, n::Int, i, c, contingency, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max)
end

"nothing to do, this model is symmetric"
function constraint_thermal_limit_to_slack_on_off(pm::AbstractDCPModel, n::Int, i, t_idx, rate_a)
end

"nothing to do, this model is symmetric"
function constraint_thermal_limit_to_slack(pm::AbstractDCPModel, n::Int, i, t_idx, rate_a)
end

"nothing to do, this model is symmetric"
function constraint_thermal_limit_to_sec(pm::AbstractDCPModel, n::Int, i::Int, t_idx, rate_a)
end

"nothing to do, this model is symmetric"
function constraint_thermal_limit_to_sec_on_off(pm::AbstractDCPModel, n::Int, i, t_idx, rate_a)
end

"nothing to do, this model is symmetric"
function constraint_thermal_limit_to_slack_sec_on_off(pm::AbstractDCPModel, n::Int, i, t_idx, rate_a)
end

""
function variable_bus_voltage_sec_on_off(pm::AbstractDCPModel, contingencies::Vector{Contingency}; kwargs...)
    variable_bus_voltage_angle_sec(pm, contingencies; kwargs...)
end

""
function variable_bus_voltage_sec_complete_on_off(pm::AbstractDCPModel, contingencies::Vector{Contingency}; kwargs...)
    variable_bus_voltage_angle_sec_complete(pm, contingencies; kwargs...)
end

"do nothing, this model does not have complex voltage variables"
function constraint_model_voltage_sec_on_off(pm::AbstractDCPModel, n::Int, contingencies::Vector{Contingency})
end

"`-b*(t[f_bus] - t[t_bus] + vad_min*(1-z_branch[i])) <= p[f_idx] <= -b*(t[f_bus] - t[t_bus] + vad_max*(1-z_branch[i]-cont))`"
function constraint_ohms_yt_from_sec(pm::AbstractDCPModel, n::Int, i, c::ContingencyID, contingency::Contingency, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm) 
    p_sec_fr  = var(pm, n,  :p_sec, f_idx) # power flow variable for this contingency
    va_sec_fr = var(pm, n, :va_sec, (f_bus,c)) # voltage angle variable for this contingency
    va_sec_to = var(pm, n, :va_sec, (t_bus,c)) # "

    shift = get_transformer_shift(pm, n, i)

    JuMP.@constraint(pm.model, p_sec_fr == -b*(va_sec_fr - va_sec_to - shift), base_name = "$(n)_C$(c)_ohm_fr_$(i)")
end

"`-b*(t[f_bus] - t[t_bus] + vad_min*(1-z_branch[i])) <= p[f_idx] <= -b*(t[f_bus] - t[t_bus] + vad_max*(1-z_branch[i]-cont))`"
function constraint_ohms_yt_from_sec_on_off(pm::AbstractDCPModel, n::Int, i, c::ContingencyID, contingency::Contingency, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max) 
    p_sec_fr  = var(pm, n,  :p_sec, f_idx) # power flow variable for this contingency
    va_sec_fr = var(pm, n, :va_sec, (f_bus,c)) # voltage angle variable for this contingency
    va_sec_to = var(pm, n, :va_sec, (t_bus,c)) # "
    z = var(pm, n, :z_branch, i)
    shift = get_transformer_shift(pm, n, i) 

    if isapprox(b, 0; atol=1e-5)
        @warn("Too small b value for branch $i, setting to 0.")
        b = 0
    end
    if b <= 0
        JuMP.@constraint(pm.model, p_sec_fr <= -b*(va_sec_fr - va_sec_to - shift + vad_max*(1-z)), base_name = "$(n)_C$(c)_ohm_fr_on_off_ub_$(i)")
        JuMP.@constraint(pm.model, p_sec_fr >= -b*(va_sec_fr - va_sec_to - shift + vad_min*(1-z)), base_name = "$(n)_C$(c)_ohm_fr_on_off_lb_$(i)")
    else # account for bound reversal when b is positive
        JuMP.@constraint(pm.model, p_sec_fr >= -b*(va_sec_fr - va_sec_to - shift + vad_max*(1-z)), base_name = "$(n)_C$(c)_ohm_fr_on_off_lb_$(i)")
        JuMP.@constraint(pm.model, p_sec_fr <= -b*(va_sec_fr - va_sec_to - shift + vad_min*(1-z)), base_name = "$(n)_C$(c)_ohm_fr_on_off_ub_$(i)")
    end
end

""
function variable_branch_power_real_sec(pm::AbstractDCPModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, pf_ub_mult_sec=PARAMS[:LIMIT_MULT][:PF_UB_MULT_SEC])
    # one variable for each arc in each contingency where it is not disabled
    sec_arcs = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs_from) for (c, contingency) in enumerate(contingencies) if l ∉ get(contingency, :branch, [])]
    
    p_sec = var(pm, nw)[:p_sec] = JuMP.@variable(pm.model,
        [((l,c),i,j) in sec_arcs], base_name="$(nw)_p_sec",
        start = comp_start_value(ref(pm, nw, :branch, l), "p_start")
    )

    if bounded
        flow_lb, flow_ub = ref_calc_branch_flow_bounds(ref(pm, nw, :branch), ref(pm, nw, :bus))

        for sec_arc in sec_arcs
            (l,c),i,j = sec_arc
            if !isinf(flow_lb[l])
                JuMP.set_lower_bound(p_sec[sec_arc], flow_lb[l])
            end
            if !isinf(flow_ub[l])
                JuMP.set_upper_bound(p_sec[sec_arc], max(pf_ub_mult_sec * flow_ub[l], flow_lb[l]))
            end
        end
    end

    for (l,branch) in ref(pm, nw, :branch)
        if haskey(branch, "pf_start")
            for (c, contingency) in enumerate(contingencies)
                if l ∉ get(contingency, :branch, [])
                    f_idx = ((l,c), branch["f_bus"], branch["t_bus"])
                    JuMP.set_start_value(p_sec[f_idx], branch["pf_start"])
                end
            end
        end
    end

    # this explicit type erasure is necessary
    p_sec_expr = Dict{Any,Any}( (((l,c),i,j), p_sec[((l,c),i,j)]) for ((l,c),i,j) in sec_arcs )
    merge!(p_sec_expr, Dict( (((l,c),j,i), -1.0*p_sec[((l,c),i,j)]) for ((l,c),i,j) in sec_arcs))
    var(pm, nw)[:p_sec] = p_sec_expr

    arcs_from_sec = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs_from) for (c, contingency) in enumerate(contingencies) if l ∉ get(contingency, :branch, [])]
    arcs_to_sec = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs_to) for (c, contingency) in enumerate(contingencies) if l ∉ get(contingency, :branch, [])]

    report && sol_component_value_edge(pm, nw, :branch, :pf, :pt, arcs_from_sec, arcs_to_sec, p_sec_expr)
end

""
function variable_branch_power_real_sec_complete(pm::AbstractDCPModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    # one variable for each arc in each contingency where it is not disabled
    sec_arcs = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs_from) for (c, contingency) in enumerate(contingencies)]
    
    p_sec = var(pm, nw)[:p_sec] = JuMP.@variable(pm.model,
        [((l,c),i,j) in sec_arcs], base_name="$(nw)_p_sec",
        start = comp_start_value(ref(pm, nw, :branch, l), "p_start")
    )

    if bounded
        flow_lb, flow_ub = ref_calc_branch_flow_bounds(ref(pm, nw, :branch), ref(pm, nw, :bus))

        for sec_arc in sec_arcs
            (l,c),i,j = sec_arc
            if !isinf(flow_lb[l])
                JuMP.set_lower_bound(p_sec[sec_arc], flow_lb[l])
            end
            if !isinf(flow_ub[l])
                JuMP.set_upper_bound(p_sec[sec_arc], flow_ub[l])
            end
        end
    end

    for (l,branch) in ref(pm, nw, :branch)
        if haskey(branch, "pf_start")
            for (c, contingency) in enumerate(contingencies)
                f_idx = ((l,c), branch["f_bus"], branch["t_bus"])
                JuMP.set_start_value(p_sec[f_idx], branch["pf_start"])
            end
        end
    end

    # this explicit type erasure is necessary
    p_sec_expr = Dict{Any,Any}( (((l,c),i,j), p_sec[((l,c),i,j)]) for ((l,c),i,j) in sec_arcs )
    merge!(p_sec_expr, Dict( (((l,c),j,i), -1.0*p_sec[((l,c),i,j)]) for ((l,c),i,j) in sec_arcs))
    var(pm, nw)[:p_sec] = p_sec_expr

    arcs_from_sec = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs_from) for (c, contingency) in enumerate(contingencies)]
    arcs_to_sec = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs_to) for (c, contingency) in enumerate(contingencies)]

    report && sol_component_value_edge(pm, nw, :branch, :pf, :pt, arcs_from_sec, arcs_to_sec, p_sec_expr)
end

"nothing to do, no voltage angle variables"
function constraint_theta_ref_sec(pm::AbstractNFAModel, contingencies::Vector{Contingency}, n::Int, i::Int)
end

""
function constraint_ohms_yt_to_sec(pm::AbstractDCPLLModel, n::Int, i, c::ContingencyID, contingency::Contingency, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    p_sec_fr = var(pm, n,  :p_sec, f_idx)
    p_sec_to = var(pm, n,  :p_sec, t_idx)
    va_fr = var(pm, n, :va_sec, (f_bus,c))
    va_to = var(pm, n, :va_sec, (t_bus,c))

    r = g/(g^2 + b^2)
    JuMP.@constraint(pm.model, p_sec_fr + p_sec_to >= r*(p_sec_fr^2), base_name="$(n)_C$(c)_ohm_to_$(i)")
end

""
function constraint_ohms_yt_to_sec_on_off(pm::AbstractDCPLLModel, n::Int, i, c::ContingencyID, contingency::Contingency, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max)
    p_fr  = var(pm, n,  :p_sec, f_idx)
    p_to  = var(pm, n,  :p_sec, t_idx)
    va_fr = var(pm, n, :va_sec, (f_bus,c))
    va_to = var(pm, n, :va_sec, (t_bus,c))
    z = var(pm, n, :z_branch, i)

    r = g/(g^2 + b^2)
    t_m = max(abs(vad_min),abs(vad_max))
    JuMP.@constraint(pm.model, p_fr + p_to >= r*( (-b*(va_fr - va_to))^2 - (-b*(t_m))^2*(1-z) ) )
    JuMP.@constraint(pm.model, p_fr + p_to >= 0)
end

"links imaginary slack variables to original variables for gens"
function constraint_gen_slack_imaginary(pm::AbstractDCPModel, n::Int, i::Int, qmin, qmax)
    # apo models ignore reactive power flow
end

"links imaginary slack variables to original variables for gens"
function constraint_gen_slack_sec_imaginary(pm::AbstractDCPModel, n::Int, i::Int, c::ContingencyID, contingency::Contingency, qmin, qmax)
    # apo models ignore reactive power flow
end

"slack variable: `v[i]` for `i` in `bus`es"
function variable_bus_voltage_magnitude_slack(pm::AbstractDCPModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    # no vm in this model
end

"slack variable: `v[i]` for `i` in `bus`es"
function variable_bus_voltage_magnitude_slack_on_off(pm::AbstractDCPModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    # no vm in this model
end

"links slack variables to original variables for bus voltage magnitude"
function constraint_bus_voltage_magnitude_slack(pm::AbstractDCPModel, n::Int, i::Int, bus, vmin, vmax)
    # no vm in this model
end

"links slack variables to original variables for bus voltage magnitude"
function constraint_bus_voltage_magnitude_slack_sec(pm::AbstractDCPModel, n::Int, i::Int, c::ContingencyID, contingency::Contingency, bus, vmin, vmax)
    # no vm in this model
end

"variable: `qq[j]` for `j` in `gen`"
function variable_gen_power_imaginary_slack_on_off(pm::AbstractDCPModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    # no reactive power in this model
end

# Modified from PowerModels.jl (https://github.com/lanl-ansi/PowerModels.jl), see LICENSE
"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[f_idx] == -b*(t[f_bus] - t[t_bus])
```
"""
function constraint_ohms_yt_from(pm::AbstractDCPModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_fr  = var(pm, n,  :p, f_idx)
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)
    (i, f_bus, t_bus) = f_idx
    shift = get_transformer_shift(pm, n, i)

    JuMP.@constraint(pm.model, p_fr <= -b*(va_fr - va_to - shift), base_name = "$(n)_ohm_fr_ub_$(i)")
    JuMP.@constraint(pm.model, p_fr >= -b*(va_fr - va_to - shift), base_name = "$(n)_ohm_fr_lb_$(i)")
    # omit reactive constraint
end

# Modified from PowerModels.jl (https://github.com/lanl-ansi/PowerModels.jl), see LICENSE
"`-b*(t[f_bus] - t[t_bus] + vad_min*(1-z_branch[i])) <= p[f_idx] <= -b*(t[f_bus] - t[t_bus] + vad_max*(1-z_branch[i]))`"
function constraint_ohms_yt_from_on_off(pm::AbstractDCPModel, n::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max)
    p_fr  = var(pm, n,  :p, f_idx)
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)

    z = var(pm, n, :z_branch, i)

    shift = get_transformer_shift(pm, n, i)
    if isapprox(b, 0; atol=1e-5)
        @warn("Too small b value for branch $i, setting to 0.")
        b = 0
    end
    if b <= 0
        JuMP.@constraint(pm.model, p_fr <= -b*(va_fr - va_to - shift + vad_max*(1-z)), base_name = "$(n)_ohm_fr_on_off_ub_$(i)")
        JuMP.@constraint(pm.model, p_fr >= -b*(va_fr - va_to - shift + vad_min*(1-z)), base_name = "$(n)_ohm_fr_on_off_lb_$(i)")
    else # account for bound reversal when b is positive
        JuMP.@constraint(pm.model, p_fr >= -b*(va_fr - va_to - shift + vad_max*(1-z)), base_name = "$(n)_ohm_fr_on_off_lb_$(i)")
        JuMP.@constraint(pm.model, p_fr <= -b*(va_fr - va_to - shift + vad_min*(1-z)), base_name = "$(n)_ohm_fr_on_off_ub_$(i)")
    end
end

# Modified from PowerModels.jl (https://github.com/lanl-ansi/PowerModels.jl), see LICENSE
"nothing to do, this model is symmetric"
function constraint_ohms_yt_to(pm::AbstractDCPModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
end


# Modified from PowerModels.jl (https://github.com/lanl-ansi/PowerModels.jl), see LICENSE
"nothing to do, this model is symmetric"
function constraint_ohms_yt_to_on_off(pm::AbstractDCPModel, n::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max)
end