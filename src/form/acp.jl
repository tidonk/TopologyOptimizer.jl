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
function constraint_power_balance_sec(pm::AbstractACPModel, n::Int, i::Int, c::ContingencyID, contingency::Contingency, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    vm_sec   = var(pm, n, :vm_sec, (i,c))
    p    = get(var(pm, n),    :p, Dict()); _check_var_keys(p, bus_arcs, "active power", "branch")

    p_sec    = get(var(pm, n),    :p_sec, Dict()); _check_var_keys(p_sec, bus_arcs, "active power", "branch")
    q_sec    = get(var(pm, n),    :q_sec, Dict()); _check_var_keys(q_sec, bus_arcs, "reactive power", "branch")
    pg_sec   = get(var(pm, n),   :pg_sec, Dict()); _check_var_keys(pg_sec, bus_gens, "active power", "generator")
    qg_sec   = get(var(pm, n),   :qg_sec, Dict()); _check_var_keys(qg_sec, bus_gens, "reactive power", "generator")
    # TODO: add ps_sec/....
    ps   = get(var(pm, n),   :ps, Dict()); _check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(var(pm, n),   :qs, Dict()); _check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(var(pm, n),  :psw, Dict()); _check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(var(pm, n),  :qsw, Dict()); _check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    p_dc = get(var(pm, n), :p_dc, Dict()); _check_var_keys(p_dc, bus_arcs_dc, "active power", "dcline")
    q_dc = get(var(pm, n), :q_dc, Dict()); _check_var_keys(q_dc, bus_arcs_dc, "reactive power", "dcline")

    # TODO these contingecies should be filtered out in advance
    if length([p_sec[((l,c),i,j)] for (l,i,j) in bus_arcs if l ∉ get(contingency, :branch, [])]) == 0
        return
    end

    # not implemented
    @assert length(ps) == length(qs) == length(psw) == length(qsw) == 0


    # the check "typeof(p[arc]) <: JuMP.NonlinearExpression" is required for the
    # case when p/q are nonlinear expressions instead of decision variables
    # once NLExpressions are first order in JuMP it should be possible to
    # remove this.
    nl_form = length(bus_arcs) > 0 && (typeof(p[iterate(bus_arcs)[1]]) <: JuMP.NonlinearExpression)

    if !nl_form
        cstr_p_sec = JuMP.@constraint(pm.model,
            sum(p_sec[((l,c),i,j)] for (l,i,j) in bus_arcs if l ∉ get(contingency, :branch, []))
            + sum(p_dc_sec[((l_dc,c),i_dc,j_dc)] for (l_dc,i_dc,j_dc) in bus_arcs_dc if l_dc ∉ get(contingency, :dcline, []))
            + sum(psw[a_sw] for a_sw in bus_arcs_sw)
            ==
            sum(pg_sec[(g, c)] for g in bus_gens if g ∉ get(contingency, :gen, []))
            - sum(ps_sec[(s,c)] for s in bus_storage if s ∉ get(contingency, :storage, []))
            - sum(pd for (i,pd) in bus_pd)
            - sum(gs for (i,gs) in bus_gs)*vm_sec^2
        )
    else
        cstr_p_sec = JuMP.@constraint(pm.model,
            sum(p_sec[((l,c),i,j)] for (l,i,j) in bus_arcs if l ∉ get(contingency, :branch, []))
            + sum(p_dc_sec[((l_dc,c),i_dc,j_dc)] for (l_dc,i_dc,j_dc) in bus_arcs_dc if l_dc ∉ get(contingency, :dcline, []))
            + sum(psw[a_sw] for a_sw in bus_arcs_sw)
            ==
            sum(pg_sec[(g, c)] for g in bus_gens if g ∉ get(contingency, :gen, []))
            - sum(ps_sec[(s,c)] for s in bus_storage if s ∉ get(contingency, :storage, []))
            - sum(pd for (i,pd) in bus_pd)
            - sum(gs for (i,gs) in bus_gs)*vm_sec^2
        )
    end

    if !nl_form
        cstr_q_sec = JuMP.@constraint(pm.model,
            sum(q_sec[((l,c),i,j)] for (l,i,j) in bus_arcs if l ∉ get(contingency, :branch, []))
            + sum(q_dc_sec[((l_dc,c),i_dc,j_dc)] for (l_dc,i_dc,j_dc) in bus_arcs_dc if l_dc ∉ get(contingency, :dcline, []))
            + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
            ==
            sum(qg_sec[(g, c)] for g in bus_gens if g ∉ get(contingency, :gen, []))
            - sum(qs_sec[(s,c)] for s in bus_storage if s ∉ get(contingency, :storage, []))
            - sum(qd for (i,qd) in bus_qd)
            - sum(bs for (i,bs) in bus_bs)*vm_sec^2
        )
    else
        cstr_q = JuMP.@NLconstraint(pm.model,
            sum(q[a] for a in bus_arcs)
            + sum(q_dc[a_dc] for a_dc in bus_arcs_dc)
            + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
            ==
            sum(qg[g] for g in bus_gens)
            - sum(qs[s] for s in bus_storage)
            - sum(qd for (i,qd) in bus_qd)
            + sum(bs for (i,bs) in bus_bs)*vm_sec^2
        )
    end

    if _IM.report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end

""
function variable_bus_voltage_sec_on_off(pm::AbstractACPModel, contingencies::Vector{Contingency}; kwargs...)
    variable_bus_voltage_angle_sec(pm, contingencies; kwargs...)
    variable_bus_voltage_magnitude_sec(pm, contingencies; kwargs...)
end

"do nothing, this model does not have complex voltage constraints"
function constraint_model_voltage_sec_on_off(pm::AbstractACPModel, n::Int, contingencies::Vector{Contingency})
    return
end


function constraint_ohms_yt_from_sec(pm::AbstractACPModel, n::Int, i, c::ContingencyID, contingency::Contingency, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_sec_fr  = var(pm, n,  :p_sec,  f_idx   )  
    q_sec_fr  = var(pm, n,  :q_sec,  f_idx   )
    vm_sec_fr = var(pm, n, :vm_sec, (f_bus,c))
    vm_sec_to = var(pm, n, :vm_sec, (t_bus,c))
    va_sec_fr = var(pm, n, :va_sec, (f_bus,c))
    va_sec_to = var(pm, n, :va_sec, (t_bus,c))

    JuMP.@NLconstraint(pm.model, p_sec_fr == (g+g_fr)/tm^2*vm_sec_fr^2 + (-g*tr+b*ti)/tm^2*(vm_sec_fr*vm_sec_to*cos(va_sec_fr-va_sec_to)) + (-b*tr-g*ti)/tm^2*(vm_sec_fr*vm_sec_to*sin(va_sec_fr-va_sec_to)) )
    JuMP.@NLconstraint(pm.model, q_sec_fr == -(b+b_fr)/tm^2*vm_sec_fr^2 - (-b*tr-g*ti)/tm^2*(vm_sec_fr*vm_sec_to*cos(va_sec_fr-va_sec_to)) + (-g*tr+b*ti)/tm^2*(vm_sec_fr*vm_sec_to*sin(va_sec_fr-va_sec_to)) )
end

function constraint_ohms_yt_from_sec_on_off(pm::AbstractACPModel, n::Int, i, c::ContingencyID, contingency::Contingency, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max)
    p_sec_fr  = var(pm, n,  :p_sec,  f_idx   )  
    q_sec_fr  = var(pm, n,  :q_sec,  f_idx   )
    vm_sec_fr = var(pm, n, :vm_sec, (f_bus,c))
    vm_sec_to = var(pm, n, :vm_sec, (t_bus,c))
    va_sec_fr = var(pm, n, :va_sec, (f_bus,c))
    va_sec_to = var(pm, n, :va_sec, (t_bus,c))
    z = var(pm, n, :z_branch, i)

    JuMP.@NLconstraint(pm.model, p_sec_fr == z*( (g+g_fr)/tm^2*vm_sec_fr^2 + (-g*tr+b*ti)/tm^2*(vm_sec_fr*vm_sec_to*cos(va_sec_fr-va_sec_to)) + (-b*tr-g*ti)/tm^2*(vm_sec_fr*vm_sec_to*sin(va_sec_fr-va_sec_to))) )
    JuMP.@NLconstraint(pm.model, q_sec_fr == z*(-(b+b_fr)/tm^2*vm_sec_fr^2 - (-b*tr-g*ti)/tm^2*(vm_sec_fr*vm_sec_to*cos(va_sec_fr-va_sec_to)) + (-g*tr+b*ti)/tm^2*(vm_sec_fr*vm_sec_to*sin(va_sec_fr-va_sec_to))) )
end

""
function constraint_ohms_yt_to_sec(pm::AbstractACPModel, n::Int, i, c::ContingencyID, contingency::Contingency, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    p_sec_to  = var(pm, n,  :p_sec,  t_idx   )
    q_sec_to  = var(pm, n,  :q_sec,  t_idx   )
    vm_sec_fr = var(pm, n, :vm_sec, (f_bus,c))
    vm_sec_to = var(pm, n, :vm_sec, (t_bus,c))
    va_sec_fr = var(pm, n, :va_sec, (f_bus,c))
    va_sec_to = var(pm, n, :va_sec, (t_bus,c))

    JuMP.@NLconstraint(pm.model, p_sec_to ==  (g+g_to)*vm_sec_to^2 + (-g*tr-b*ti)/tm^2*(vm_sec_to*vm_sec_fr*cos(va_sec_to-va_sec_fr)) + (-b*tr+g*ti)/tm^2*(vm_sec_to*vm_sec_fr*sin(va_sec_to-va_sec_fr)) )
    JuMP.@NLconstraint(pm.model, q_sec_to == -(b+b_to)*vm_sec_to^2 - (-b*tr+g*ti)/tm^2*(vm_sec_to*vm_sec_fr*cos(va_sec_to-va_sec_fr)) + (-g*tr-b*ti)/tm^2*(vm_sec_to*vm_sec_fr*sin(va_sec_to-va_sec_fr)) )
end

""
function constraint_ohms_yt_to_sec_on_off(pm::AbstractACPModel, n::Int, i, c::ContingencyID, contingency::Contingency, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max)
    p_sec_to  = var(pm, n,  :p_sec,  t_idx   )
    q_sec_to  = var(pm, n,  :q_sec,  t_idx   )
    vm_sec_fr = var(pm, n, :vm_sec, (f_bus,c))
    vm_sec_to = var(pm, n, :vm_sec, (t_bus,c))
    va_sec_fr = var(pm, n, :va_sec, (f_bus,c))
    va_sec_to = var(pm, n, :va_sec, (t_bus,c))
    z = var(pm, n, :z_branch, i)

    JuMP.@NLconstraint(pm.model, p_sec_to == z*( (g+g_to)*vm_sec_to^2 + (-g*tr-b*ti)/tm^2*(vm_sec_to*vm_sec_fr*cos(va_sec_to-va_sec_fr)) + (-b*tr+g*ti)/tm^2*(vm_sec_to*vm_sec_fr*sin(va_sec_to-va_sec_fr))) )
    JuMP.@NLconstraint(pm.model, q_sec_to == z*(-(b+b_to)*vm_sec_to^2 - (-b*tr+g*ti)/tm^2*(vm_sec_to*vm_sec_fr*cos(va_sec_to-va_sec_fr)) + (-g*tr-b*ti)/tm^2*(vm_sec_to*vm_sec_fr*sin(va_sec_to-va_sec_fr))) )
end

"links imaginary slack variables to original variables for gens"
function constraint_gen_slack_imaginary(pm::AbstractACPModel, n::Int, i::Int, qmin, qmax)
    qg = get(var(pm, n), :qg, Dict())
    qmax_slack = get(var(pm, n), :qmax_slack, Dict())
    qmin_slack = get(var(pm, n), :qmin_slack, Dict())

    JuMP.@constraint(pm.model, qg[i] <= qmax + qmax_slack[i])
    JuMP.@constraint(pm.model, qg[i] >= qmin - qmin_slack[i])
end

"links imaginary slack variables to original variables for gens"
function constraint_gen_slack_sec_imaginary(pm::AbstractACPModel, n::Int, i::Int, c::ContingencyID, contingency::Contingency, qmin, qmax)
    qg_sec = get(var(pm, n), :qg_sec, Dict())
    qmax_slack = get(var(pm, n), :qmax_slack, Dict())
    qmin_slack = get(var(pm, n), :qmin_slack, Dict())

    JuMP.@constraint(pm.model, qg_sec[(i,c)] <= qmax + qmax_slack[i])
    JuMP.@constraint(pm.model, qg_sec[(i,c)] >= qmin - qmin_slack[i])
end

"variable: `pg_delta_percentage_slack[j]` for `j` in `gen`"
function variable_pg_delta_percentage_slack_on_off(pm::AbstractPowerModel; nw::Int=nw_id_default, report::Bool=true)
    pg_delta_percentage_slack_on_off = var(pm, nw)[:pg_delta_percentage_slack_on_off] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :gen)], base_name="$(nw)_pg_delta_percentage_slack_on_off",
        binary = true,
        start = 0
    )

    report && sol_component_value(pm, nw, :gen, :pg_delta_percentage_slack_on_off, ids(pm, nw, :gen), pg_delta_percentage_slack_on_off)
end

"variable: `qq[j]` for `j` in `gen`"
function variable_gen_power_imaginary_slack_on_off(pm::AbstractACPModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    qmin_slack_on_off = var(pm, nw)[:qmin_slack_on_off] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :gen)], base_name="$(nw)_qmin_slack_on_off",
        binary = true,
        start = 0
    )

    qmax_slack_on_off = var(pm, nw)[:qmax_slack_on_off] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :gen)], base_name="$(nw)_qmax_slack_on_off",
        binary = true,
        start = 0
    )

    report && sol_component_value(pm, nw, :gen, :qmin_slack_on_off, ids(pm, nw, :gen), qmin_slack_on_off)
    report && sol_component_value(pm, nw, :gen, :qmax_slack_on_off, ids(pm, nw, :gen), qmax_slack_on_off)
end

"`angmin <= z_branch[i]*(t[f_bus] - t[t_bus]) <= angmax` + slack linking"
function constraint_voltage_angle_difference_slack(pm::AbstractACPModel, n::Int, f_idx, angmin, angmax, vad_min, vad_max)
    i, f_bus, t_bus = f_idx
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)
    
    vad_slack_lower = var(pm, n, :vad_slack_lower, f_idx)
    vad_slack_upper = var(pm, n, :vad_slack_upper, f_idx)

    # slack linking
    JuMP.@constraint(pm.model, (va_fr - va_to) <= angmax + vad_slack_upper)
    JuMP.@constraint(pm.model, (va_fr - va_to) >= angmin - vad_slack_lower)
end

"`angmin <= z_branch[i]*(t[f_bus] - t[t_bus]) <= angmax` + slack linking"
function constraint_voltage_angle_difference_slack_on_off(pm::AbstractACPModel, n::Int, f_idx, angmin, angmax, vad_min, vad_max)
    i, f_bus, t_bus = f_idx
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)
    z = var(pm, n, :z_branch, i)
    
    vad_slack_lower = var(pm, n, :vad_slack_lower, f_idx)
    vad_slack_upper = var(pm, n, :vad_slack_upper, f_idx)

    # switched on?
    JuMP.@constraint(pm.model, z*(va_fr - va_to) <= z*vad_max)
    JuMP.@constraint(pm.model, z*(va_fr - va_to) >= z*vad_min)

    # slack linking
    JuMP.@constraint(pm.model, z*(va_fr - va_to) <= angmax + vad_slack_upper)
    JuMP.@constraint(pm.model, z*(va_fr - va_to) >= angmin - vad_slack_lower)
end

"`angmin <= z_branch[i]*(t[f_bus] - t[t_bus]) <= angmax`"
function constraint_voltage_angle_difference_sec(pm::AbstractACPModel, n::Int, c::ContingencyID, contingency::Contingency, f_idx, angmin, angmax)
    (i, c), f_bus, t_bus = f_idx

    va_sec_fr = var(pm, n, :va_sec, (f_bus,c))
    va_sec_to = var(pm, n, :va_sec, (t_bus,c))

    JuMP.@constraint(pm.model, va_sec_fr - va_sec_to <= angmax)
    JuMP.@constraint(pm.model, va_sec_fr - va_sec_to >= angmin)
end

"`angmin <= z_branch[i]*(t[f_bus] - t[t_bus]) <= angmax`"
function constraint_voltage_angle_difference_sec_on_off(pm::AbstractACPModel, n::Int, c::ContingencyID, contingency::Contingency, f_idx, angmin, angmax, vad_min, vad_max)
    (i, c), f_bus, t_bus = f_idx

    va_sec_fr = var(pm, n, :va_sec, (f_bus,c))
    va_sec_to = var(pm, n, :va_sec, (t_bus,c))
    z = var(pm, n, :z_branch, i)

    JuMP.@constraint(pm.model, z*(va_sec_fr - va_sec_to) <= z*angmax)
    JuMP.@constraint(pm.model, z*(va_sec_fr - va_sec_to) >= z*angmin)
end

"`angmin <= z_branch[i]*(t[f_bus] - t[t_bus]) <= angmax` +  slack linking"
function constraint_voltage_angle_difference_slack_sec(pm::AbstractACPModel, n::Int, c::ContingencyID, contingency::Contingency, f_idx, angmin, angmax, vad_min, vad_max)
    (i, c), f_bus, t_bus = f_idx

    va_sec_fr = var(pm, n, :va_sec, (f_bus,c))
    va_sec_to = var(pm, n, :va_sec, (t_bus,c))
    
    vad_slack_lower = var(pm, n, :vad_slack_lower, (i, f_bus, t_bus))
    vad_slack_upper = var(pm, n, :vad_slack_upper, (i, f_bus, t_bus))

     # slack linking
    JuMP.@constraint(pm.model, (va_sec_fr - va_sec_to) <= angmax + vad_slack_upper)
    JuMP.@constraint(pm.model, (va_sec_fr - va_sec_to) >= angmin - vad_slack_lower)
end

"`angmin <= z_branch[i]*(t[f_bus] - t[t_bus]) <= angmax` +  slack linking"
function constraint_voltage_angle_difference_slack_sec_on_off(pm::AbstractACPModel, n::Int, c::ContingencyID, contingency::Contingency, f_idx, angmin, angmax, vad_min, vad_max)
    (i, c), f_bus, t_bus = f_idx

    va_sec_fr = var(pm, n, :va_sec, (f_bus,c))
    va_sec_to = var(pm, n, :va_sec, (t_bus,c))
    z = var(pm, n, :z_branch, i)
    
    vad_slack_lower = var(pm, n, :vad_slack_lower, (i, f_bus, t_bus))
    vad_slack_upper = var(pm, n, :vad_slack_upper, (i, f_bus, t_bus))
    
    # switched on?
    JuMP.@constraint(pm.model, z*(va_sec_fr - va_sec_to) <= z*vad_max)
    JuMP.@constraint(pm.model, z*(va_sec_fr - va_sec_to) >= z*vad_min)

     # slack linking
    JuMP.@constraint(pm.model, z*(va_sec_fr - va_sec_to) <= angmax + vad_slack_upper)
    JuMP.@constraint(pm.model, z*(va_sec_fr - va_sec_to) >= angmin - vad_slack_lower)
end

"slack variable: `v[i]` for `i` in `bus`es"
function variable_bus_voltage_magnitude_slack(pm::AbstractACPModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    vm_slack_lower = var(pm, nw)[:vm_slack_lower] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :bus)], base_name="$(nw)_vm_slack_lower",
        start = 0.0
    )

    vm_slack_upper = var(pm, nw)[:vm_slack_upper] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :bus)], base_name="$(nw)_vm_slack_upper",
        start = 0.0
    )
    
    if bounded
        for (i, bus) in ref(pm, nw, :bus)
            JuMP.set_lower_bound(vm_slack_lower[i], 0.0)
            JuMP.set_upper_bound(vm_slack_lower[i], bus["vmin"])

            JuMP.set_lower_bound(vm_slack_upper[i], 0.0)
            #JuMP.set_upper_bound(vm_slack_upper[i], bus["vmax"])
        end
    end

    report && sol_component_value(pm, nw, :bus, :vm_slack_lower, ids(pm, nw, :bus), vm_slack_lower)
    report && sol_component_value(pm, nw, :bus, :vm_slack_upper, ids(pm, nw, :bus), vm_slack_upper)
end

"slack variable: `v[i]` for `i` in `bus`es"
function variable_bus_voltage_magnitude_slack_on_off(pm::AbstractACPModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    vm_slack_lower_on_off = var(pm, nw)[:vm_slack_lower_on_off] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :bus)], base_name="$(nw)_vm_slack_lower_on_off",
        binary = true,
        start = 0.0
    )

    vm_slack_upper_on_off = var(pm, nw)[:vm_slack_upper_on_off] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :bus)], base_name="$(nw)_vm_slack_upper_on_off",
        binary = true,
        start = 0.0
    )

    report && sol_component_value(pm, nw, :bus, :vm_slack_lower_on_off, ids(pm, nw, :bus), vm_slack_lower_on_off)
    report && sol_component_value(pm, nw, :bus, :vm_slack_upper_on_off, ids(pm, nw, :bus), vm_slack_upper_on_off)
end

"links slack variables to original variables for bus voltage magnitude"
function constraint_bus_voltage_magnitude_slack(pm::AbstractACPModel, n::Int, i::Int, bus, vmin, vmax)
    vm = get(var(pm, n), :vm, Dict())
    vm_slack_upper = get(var(pm, n), :vm_slack_upper, Dict())
    vm_slack_lower = get(var(pm, n), :vm_slack_lower, Dict())

    JuMP.@constraint(pm.model, vm[i] <= vmax + vm_slack_upper[i])
    JuMP.@constraint(pm.model, vm[i] >= vmin - vm_slack_lower[i])
end

"links slack variables to original variables for bus voltage magnitude"
function constraint_bus_voltage_magnitude_slack_sec(pm::AbstractACPModel, n::Int, i::Int, c::ContingencyID, contingency::Contingency, bus, vmin, vmax)
    vm_sec = get(var(pm, n), :vm_sec, Dict())
    vm_slack_upper = get(var(pm, n), :vm_slack_upper, Dict())
    vm_slack_lower = get(var(pm, n), :vm_slack_lower, Dict())

    JuMP.@constraint(pm.model, vm_sec[(i,c)] <= vmax + vm_slack_upper[i])
    JuMP.@constraint(pm.model, vm_sec[(i,c)] >= vmin - vm_slack_lower[i])
end

# Modified from PowerModels.jl (https://github.com/lanl-ansi/PowerModels.jl), see LICENSE
""
function constraint_power_balance(pm::AbstractACPModel, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    vm   = var(pm, n, :vm, i)
    p    = get(var(pm, n),    :p, Dict()); _check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(var(pm, n),    :q, Dict()); _check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(var(pm, n),   :pg, Dict()); _check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(var(pm, n),   :qg, Dict()); _check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(var(pm, n),   :ps, Dict()); _check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(var(pm, n),   :qs, Dict()); _check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(var(pm, n),  :psw, Dict()); _check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(var(pm, n),  :qsw, Dict()); _check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    p_dc = get(var(pm, n), :p_dc, Dict()); _check_var_keys(p_dc, bus_arcs_dc, "active power", "dcline")
    q_dc = get(var(pm, n), :q_dc, Dict()); _check_var_keys(q_dc, bus_arcs_dc, "reactive power", "dcline")

    if :pd in get(pm.data, "slack_types", [])
        pd_slack_lower = get(var(pm, n), :pd_slack_lower, Dict()); _check_var_keys(pd_slack_lower, bus_arcs_dc, "pd slack lower", "load")
        pd_slack_upper = get(var(pm, n), :pd_slack_upper, Dict()); _check_var_keys(pd_slack_upper, bus_arcs_dc, "pd slack upper", "load")
        bus_loads_pd = [l_id for l_id in ids(pm, n, :load) if ref(pm, n, :load, l_id)["load_bus"] == i]
    else
        pd_slack_lower = Dict()
        pd_slack_upper = Dict()
        bus_loads_pd = []
    end

    if :qd in get(pm.data, "slack_types", [])
        qd_slack_lower = get(var(pm, n), :qd_slack_lower, Dict()); _check_var_keys(qd_slack_lower, bus_arcs_dc, "qd slack lower", "load")
        qd_slack_upper = get(var(pm, n), :qd_slack_upper, Dict()); _check_var_keys(qd_slack_upper, bus_arcs_dc, "qd slack upper", "load")
        bus_loads_qd = [l_id for l_id in ids(pm, n, :load) if ref(pm, n, :load, l_id)["load_bus"] == i]
    else
        qd_slack_lower = Dict()
        qd_slack_upper = Dict()
        bus_loads_qd = []
    end

    # the check "typeof(p[arc]) <: JuMP.NonlinearExpression" is required for the
    # case when p/q are nonlinear expressions instead of decision variables
    # once NLExpressions are first order in JuMP it should be possible to
    # remove this.
    nl_form = length(bus_arcs) > 0 && (typeof(p[iterate(bus_arcs)[1]]) <: JuMP.NonlinearExpression)

    if !nl_form
        cstr_p = JuMP.@constraint(pm.model,
            sum(p[a] for a in bus_arcs)
            + sum(p_dc[a_dc] for a_dc in bus_arcs_dc)
            + sum(psw[a_sw] for a_sw in bus_arcs_sw)
            ==
            sum(pg[g] for g in bus_gens)
            - sum(ps[s] for s in bus_storage)
            - sum(pd for (i,pd) in bus_pd)
            - sum(pd_slack_upper[load] for load in bus_loads_pd)
            + sum(pd_slack_lower[load] for load in bus_loads_pd)
            - sum(gs for (i,gs) in bus_gs)*vm^2
        )
    else
        cstr_p = JuMP.@NLconstraint(pm.model,
            sum(p[a] for a in bus_arcs)
            + sum(p_dc[a_dc] for a_dc in bus_arcs_dc)
            + sum(psw[a_sw] for a_sw in bus_arcs_sw)
            ==
            sum(pg[g] for g in bus_gens)
            - sum(ps[s] for s in bus_storage)
            - sum(pd for (i,pd) in bus_pd)
            - sum(pd_slack_upper[load] for load in bus_loads_pd)
            + sum(pd_slack_lower[load] for load in bus_loads_pd)
            - sum(gs for (i,gs) in bus_gs)*vm^2
        )
    end

    if !nl_form
        cstr_q = JuMP.@constraint(pm.model,
            sum(q[a] for a in bus_arcs)
            + sum(q_dc[a_dc] for a_dc in bus_arcs_dc)
            + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
            ==
            sum(qg[g] for g in bus_gens)
            - sum(qs[s] for s in bus_storage)
            - sum(qd for (i,qd) in bus_qd)
            - sum(qd_slack_upper[load] for load in bus_loads_qd)
            + sum(qd_slack_lower[load] for load in bus_loads_qd)
            + sum(bs for (i,bs) in bus_bs)*vm^2
        )
    else
        cstr_q = JuMP.@NLconstraint(pm.model,
            sum(q[a] for a in bus_arcs)
            + sum(q_dc[a_dc] for a_dc in bus_arcs_dc)
            + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
            ==
            sum(qg[g] for g in bus_gens)
            - sum(qs[s] for s in bus_storage)
            - sum(qd for (i,qd) in bus_qd)
            - sum(qd_slack_upper[load] for load in bus_loads_qd)
            + sum(qd_slack_lower[load] for load in bus_loads_qd)
            + sum(bs for (i,bs) in bus_bs)*vm^2
        )
    end

    if _IM.report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end

# Modified from PowerModels.jl (https://github.com/lanl-ansi/PowerModels.jl), see LICENSE
"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[f_idx] ==  (g+g_fr)/tm*v[f_bus]^2 + (-g*tr+b*ti)/tm^2*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-b*tr-g*ti)/tm^2*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus]))
q[f_idx] == -(b+b_fr)/tm*v[f_bus]^2 - (-b*tr-g*ti)/tm^2*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr+b*ti)/tm^2*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus]))
```
"""
function constraint_ohms_yt_from(pm::AbstractACPModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_fr  = var(pm, n,  :p, f_idx)
    q_fr  = var(pm, n,  :q, f_idx)
    vm_fr = var(pm, n, :vm, f_bus)
    vm_to = var(pm, n, :vm, t_bus)
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)

    JuMP.@NLconstraint(pm.model, p_fr ==  (g+g_fr)/tm^2*vm_fr^2 + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to)) )
    JuMP.@NLconstraint(pm.model, q_fr == -(b+b_fr)/tm^2*vm_fr^2 - (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to)) )
end

# Modified from PowerModels.jl (https://github.com/lanl-ansi/PowerModels.jl), see LICENSE
"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)

```
p[t_idx] ==  (g+g_to)*v[t_bus]^2 + (-g*tr-b*ti)/tm^2*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm^2*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
q[t_idx] == -(b+b_to)*v[t_bus]^2 - (-b*tr+g*ti)/tm^2*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm^2*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
```
"""
function constraint_ohms_yt_to(pm::AbstractACPModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    p_to  = var(pm, n,  :p, t_idx)
    q_to  = var(pm, n,  :q, t_idx)
    vm_fr = var(pm, n, :vm, f_bus)
    vm_to = var(pm, n, :vm, t_bus)
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)

    JuMP.@NLconstraint(pm.model, p_to ==  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/tm^2*(vm_to*vm_fr*sin(va_to-va_fr)) )
    JuMP.@NLconstraint(pm.model, q_to == -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-g*tr-b*ti)/tm^2*(vm_to*vm_fr*sin(va_to-va_fr)) )
end