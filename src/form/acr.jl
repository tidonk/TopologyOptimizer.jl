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

"reference bus angle constraint"
function constraint_theta_ref_sec(pm::AbstractACRModel, n::Int, i::Int, c::ContingencyID, contingency::Dict)
    JuMP.@constraint(pm.model, var(pm, n, :vi_sec)[(i,c)] == 0)
end

function variable_bus_voltage_sec(pm::AbstractACRModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, bounded::Bool=true, kwargs...)
    variable_bus_voltage_real_sec(pm, contingencies; nw=nw, bounded=bounded, kwargs...)
    variable_bus_voltage_imaginary_sec(pm, contingencies; nw=nw, bounded=bounded, kwargs...)

    if bounded
        for (c, contingency) in enumerate(contingencies)
            for (i, bus) in ref(pm, nw, :bus)
                constraint_voltage_magnitude_bounds_sec(pm, nw, i, c, contingency, bus["vmin"], bus["vmax"])    
            end
        end
    end
end

"`vmin <= vm[i] <= vmax`"
### Only adds constraint name
function constraint_voltage_magnitude_bounds(pm::AbstractACRModel, n::Int, i, vmin, vmax)
    @assert vmin <= vmax
    vr = var(pm, n, :vr, i)
    vi = var(pm, n, :vi, i)

    JuMP.@constraint(pm.model, vmin^2 <= (vr^2 + vi^2), base_name = "$(n)_voltage_magnitude_lb_$(i)")
    JuMP.@constraint(pm.model, vmax^2 >= (vr^2 + vi^2), base_name = "$(n)_voltage_magnitude_ub_$(i)")
end

"`vmin <= vm[i] <= vmax`"
function constraint_voltage_magnitude_bounds_sec(pm::AbstractACRModel, n::Int, i, c, contingency, vmin, vmax)
    @assert vmin <= vmax
    vr_sec = var(pm, n, :vr_sec, (i,c))
    vi_sec = var(pm, n, :vi_sec, (i,c))

    JuMP.@constraint(pm.model, vmin^2 <= (vr_sec^2 + vi_sec^2), base_name = "$(n)_voltage_magnitude_lb_($(i),$(c))")
    JuMP.@constraint(pm.model, vmax^2 >= (vr_sec^2 + vi_sec^2), base_name = "$(n)_voltage_magnitude_ub_($(i),$(c))")
end

"reference bus angle constraint"
function constraint_theta_ref(pm::AbstractACRModel, n::Int, i::Int)
    JuMP.@constraint(pm.model, var(pm, n, :vi)[i] == 0, base_name = "$(n)_theta_ref_$(i)")
end

function constraint_power_balance(pm::AbstractACRModel, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    vr = var(pm, n, :vr, i)
    vi = var(pm, n, :vi, i)
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

    cstr_p = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(p_dc[a_dc] for a_dc in bus_arcs_dc)
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        ==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pd for pd in values(bus_pd))
        - sum(pd_slack_upper[load] for load in bus_loads_pd)
        + sum(pd_slack_lower[load] for load in bus_loads_pd)
        - sum(gs for gs in values(bus_gs))*(vr^2 + vi^2),
        base_name = "$(n)_power_balance_p_$(i)"
    )
    cstr_q = JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        + sum(q_dc[a_dc] for a_dc in bus_arcs_dc)
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        ==
        sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        - sum(qd for qd in values(bus_qd))
        - sum(qd_slack_upper[load] for load in bus_loads_qd)
        + sum(qd_slack_lower[load] for load in bus_loads_qd)
        + sum(bs for bs in values(bus_bs))*(vr^2 + vi^2),
        base_name = "$(n)_power_balance_q_$(i)"
    )

    if _IM.report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end

function constraint_power_balance_sec(pm::AbstractACRModel, n::Int, i::Int, c::ContingencyID, contingency::Contingency, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    vr_sec = var(pm, n, :vr_sec, (i,c))
    vi_sec = var(pm, n, :vi_sec, (i,c))
    p_sec    = get(var(pm, n),  :p_sec, Dict()); _check_var_keys(p_sec, bus_arcs, "active power", "branch")
    q_sec    = get(var(pm, n),  :q_sec, Dict()); _check_var_keys(q_sec, bus_arcs, "reactive power", "branch")
    pg_sec   = get(var(pm, n),  :pg_sec, Dict()); _check_var_keys(pg_sec, bus_gens, "active power", "generator")
    qg_sec   = get(var(pm, n),  :qg_sec, Dict()); _check_var_keys(qg_sec, bus_gens, "reactive power", "generator")
    ps_sec   = get(var(pm, n),  :ps_sec, Dict()); _check_var_keys(ps_sec, bus_storage, "active power", "storage")
    qs_sec   = get(var(pm, n),  :qs_sec, Dict()); _check_var_keys(qs_sec, bus_storage, "reactive power", "storage")
    psw  = get(var(pm, n),  :psw, Dict()); _check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(var(pm, n),  :qsw, Dict()); _check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    p_dc_sec = get(var(pm, n),  :p_dc_sec, Dict()); _check_var_keys(p_dc_sec, bus_arcs_dc, "active power", "dcline")
    q_dc_sec = get(var(pm, n),  :q_dc_sec, Dict()); _check_var_keys(q_dc_sec, bus_arcs_dc, "reactive power", "dcline")

    if :pd in get(pm.data, "slack_types", [])
        pd_slack_lower_sec = get(var(pm, n), :pd_slack_lower_sec, Dict()); _check_var_keys(pd_slack_lower_sec, bus_arcs_dc, "pd slack lower", "load")
        pd_slack_upper_sec = get(var(pm, n), :pd_slack_upper_sec, Dict()); _check_var_keys(pd_slack_upper_sec, bus_arcs_dc, "pd slack upper", "load")
        bus_loads_pd_sec = [l_id for l_id in ids(pm, n, :load) if ref(pm, n, :load, l_id)["load_bus"] == i]
    else
        pd_slack_lower_sec = Dict()
        pd_slack_upper_sec = Dict()
        bus_loads_pd_sec = []
    end

    if :qd in get(pm.data, "slack_types", [])
        qd_slack_lower_sec = get(var(pm, n), :qd_slack_lower_sec, Dict()); _check_var_keys(qd_slack_lower_sec, bus_arcs_dc, "qd slack lower", "load")
        qd_slack_upper_sec = get(var(pm, n), :qd_slack_upper_sec, Dict()); _check_var_keys(qd_slack_upper_sec, bus_arcs_dc, "qd slack upper", "load")
        bus_loads_qd_sec = [l_id for l_id in ids(pm, n, :load) if ref(pm, n, :load, l_id)["load_bus"] == i]
    else
        qd_slack_lower_sec = Dict()
        qd_slack_upper_sec = Dict()
        bus_loads_qd_sec = []
    end

    cstr_p_sec = JuMP.@constraint(pm.model,
        sum(p_sec[((l,c),i,j)] for (l,i,j) in bus_arcs if l ∉ get(contingency, :branch, []))
        + sum(p_dc_sec[((l_dc,c),i_dc,j_dc)] for (l_dc,i_dc,j_dc) in bus_arcs_dc if l_dc ∉ get(contingency, :dcline, []))
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        ==
        sum(pg_sec[(g, c)] for g in bus_gens if g ∉ get(contingency, :gen, []))
        - sum(ps_sec[(s,c)] for s in bus_storage if s ∉ get(contingency, :storage, []))
        - sum(pd for pd in values(bus_pd))
        - sum(pd_slack_upper_sec[(load,c)] for load in bus_loads_pd_sec)
        + sum(pd_slack_lower_sec[(load,c)] for load in bus_loads_pd_sec)
        - sum(gs for gs in values(bus_gs))*(vr_sec^2 + vi_sec^2),
        base_name = "$(n)_power_balance_p_sec_($(i),$(c))"
    )

    cstr_q_sec = JuMP.@constraint(pm.model,
        sum(q_sec[((l,c),i,j)] for (l,i,j) in bus_arcs if l ∉ get(contingency, :branch, []))
        + sum(q_dc_sec[((l_dc,c),i_dc,j_dc)] for (l_dc,i_dc,j_dc) in bus_arcs_dc if l_dc ∉ get(contingency, :dcline, []))
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        ==
        sum(qg_sec[(g, c)] for g in bus_gens if g ∉ get(contingency, :gen, []))
        - sum(qs_sec[(s,c)] for s in bus_storage if s ∉ get(contingency, :storage, []))
        - sum(qd_sec for qd_sec in values(bus_qd))
        - sum(qd_slack_upper_sec[(load,c)] for load in bus_loads_qd_sec)
        + sum(qd_slack_lower_sec[(load,c)] for load in bus_loads_qd_sec)
        + sum(bs for bs in values(bus_bs))*(vr_sec^2 + vi_sec^2),
        base_name = "$(n)_power_balance_q_sec_($(i),$(c))"
    )

    if _IM.report_duals(pm)
        sol(pm, n, :bus, (i,c))[:lam_kcl_r] = cstr_p_sec
        sol(pm, n, :bus, (i,c))[:lam_kcl_i] = cstr_q_sec
    end
end

"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)
"""
### Only adds constraint name
function constraint_ohms_yt_from(pm::AbstractACRModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)
    vr_fr = var(pm, n, :vr, f_bus)
    vr_to = var(pm, n, :vr, t_bus)
    vi_fr = var(pm, n, :vi, f_bus)
    vi_to = var(pm, n, :vi, t_bus)

    JuMP.@constraint(pm.model, p_fr ==  (g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to), base_name = "$(n)_ohm_fr_p_$(f_idx)")
    JuMP.@constraint(pm.model, q_fr == -(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to), base_name = "$(n)_ohm_fr_q_$(f_idx)")
end

"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)
"""
### Only adds constraint name
function constraint_ohms_yt_to(pm::AbstractACRModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    p_to = var(pm, n, :p, t_idx)
    q_to = var(pm, n, :q, t_idx)
    vr_fr = var(pm, n, :vr, f_bus)
    vr_to = var(pm, n, :vr, t_bus)
    vi_fr = var(pm, n, :vi, f_bus)
    vi_to = var(pm, n, :vi, t_bus)

    JuMP.@constraint(pm.model, p_to ==  (g+g_to)*(vr_to^2 + vi_to^2) + (-g*tr-b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr+g*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to)), base_name = "$(n)_ohm_to_p_$(t_idx)")
    JuMP.@constraint(pm.model, q_to == -(b+b_to)*(vr_to^2 + vi_to^2) - (-b*tr+g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr-b*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to)), base_name = "$(n)_ohm_to_q_$(t_idx)")
end

"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)
"""
function constraint_ohms_yt_from_sec(pm::AbstractACRModel, n::Int, i, c, contingency, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_fr_sec = var(pm, n, :p_sec, f_idx)
    q_fr_sec = var(pm, n, :q_sec, f_idx)
    vr_fr_sec = var(pm, n, :vr_sec, (f_bus,c))
    vr_to_sec = var(pm, n, :vr_sec, (t_bus,c))
    vi_fr_sec = var(pm, n, :vi_sec, (f_bus,c))
    vi_to_sec = var(pm, n, :vi_sec, (t_bus,c))

    JuMP.@constraint(pm.model, p_fr_sec ==  (g+g_fr)/tm^2*(vr_fr_sec^2 + vi_fr_sec^2) + (-g*tr+b*ti)/tm^2*(vr_fr_sec*vr_to_sec + vi_fr_sec*vi_to_sec) + (-b*tr-g*ti)/tm^2*(vi_fr_sec*vr_to_sec - vr_fr_sec*vi_to_sec), base_name = "$(n)_ohm_fr_p_sec_$(f_idx)")
    JuMP.@constraint(pm.model, q_fr_sec == -(b+b_fr)/tm^2*(vr_fr_sec^2 + vi_fr_sec^2) - (-b*tr-g*ti)/tm^2*(vr_fr_sec*vr_to_sec + vi_fr_sec*vi_to_sec) + (-g*tr+b*ti)/tm^2*(vi_fr_sec*vr_to_sec - vr_fr_sec*vi_to_sec), base_name = "$(n)_ohm_fr_q_sec_$(f_idx)")
end

"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)
"""
function constraint_ohms_yt_to_sec(pm::AbstractACRModel, n::Int, i, c, contingency, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    p_to_sec = var(pm, n, :p_sec, t_idx)
    q_to_sec = var(pm, n, :q_sec, t_idx)
    vr_fr_sec = var(pm, n, :vr_sec, (f_bus,c))
    vr_to_sec = var(pm, n, :vr_sec, (t_bus,c))
    vi_fr_sec = var(pm, n, :vi_sec, (f_bus,c))
    vi_to_sec = var(pm, n, :vi_sec, (t_bus,c))

    JuMP.@constraint(pm.model, p_to_sec ==  (g+g_to)*(vr_to_sec^2 + vi_to_sec^2) + (-g*tr-b*ti)/tm^2*(vr_fr_sec*vr_to_sec + vi_fr_sec*vi_to_sec) + (-b*tr+g*ti)/tm^2*(-(vi_fr_sec*vr_to_sec - vr_fr_sec*vi_to_sec)), base_name = "$(n)_ohm_to_p_sec_$(t_idx)")
    JuMP.@constraint(pm.model, q_to_sec == -(b+b_to)*(vr_to_sec^2 + vi_to_sec^2) - (-b*tr+g*ti)/tm^2*(vr_fr_sec*vr_to_sec + vi_fr_sec*vi_to_sec) + (-g*tr-b*ti)/tm^2*(-(vi_fr_sec*vr_to_sec - vr_fr_sec*vi_to_sec)), base_name = "$(n)_ohm_to_q_sec_$(t_idx)")
end

"""
branch voltage angle difference bounds
"""
### Only adds constraint name
function constraint_voltage_angle_difference(pm::AbstractACRModel, n::Int, f_idx, angmin, angmax)
    i, f_bus, t_bus = f_idx

    vr_fr = var(pm, n, :vr, f_bus)
    vr_to = var(pm, n, :vr, t_bus)
    vi_fr = var(pm, n, :vi, f_bus)
    vi_to = var(pm, n, :vi, t_bus)

    JuMP.@constraint(pm.model, (vi_fr*vr_to - vr_fr*vi_to) <= tan(angmax)*(vr_fr*vr_to + vi_fr*vi_to), base_name = "$(n)_vad_ub_$(f_idx)")
    JuMP.@constraint(pm.model, (vi_fr*vr_to - vr_fr*vi_to) >= tan(angmin)*(vr_fr*vr_to + vi_fr*vi_to), base_name = "$(n)_vad_lb_$(f_idx)")
end

"""
branch voltage angle difference bounds
"""
function constraint_voltage_angle_difference_sec(pm::AbstractACRModel, n::Int, c::ContingencyID, contingency::Contingency, f_idx, angmin, angmax)
    i, f_bus, t_bus = f_idx

    vr_fr_sec = var(pm, n, :vr_sec, (f_bus,c))
    vr_to_sec = var(pm, n, :vr_sec, (t_bus,c))
    vi_fr_sec = var(pm, n, :vi_sec, (f_bus,c))
    vi_to_sec = var(pm, n, :vi_sec, (t_bus,c))

    JuMP.@constraint(pm.model, (vi_fr_sec*vr_to_sec - vr_fr_sec*vi_to_sec) <= tan(angmax)*(vr_fr_sec*vr_to_sec + vi_fr_sec*vi_to_sec), base_name = "$(n)_vad_ub_sec_$(f_idx)")
    JuMP.@constraint(pm.model, (vi_fr_sec*vr_to_sec - vr_fr_sec*vi_to_sec) >= tan(angmin)*(vr_fr_sec*vr_to_sec + vi_fr_sec*vi_to_sec), base_name = "$(n)_vad_lb_sec_$(f_idx)")
end

"variable: `qq[j]` for `j` in `gen`"
function PowerModels.variable_gen_power_imaginary(pm::ACRPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    qg = var(pm, nw)[:qg] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :gen)], base_name="$(nw)_qg",
        start = comp_start_value(ref(pm, nw, :gen, i), "qg_start")
    )

    if bounded
        for (i, gen) in ref(pm, nw, :gen)
            JuMP.set_lower_bound(qg[i], gen["qmin"])
            JuMP.set_upper_bound(qg[i], gen["qmax"])
        end
    end

    report && sol_component_value(pm, nw, :gen, :qg, ids(pm, nw, :gen), qg)
end
