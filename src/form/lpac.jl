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
function variable_bus_voltage_magnitude_sec(pm::AbstractLPACModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    sec_buses = [(i,c) for i in ids(pm, nw, :bus) for (c, contingency) in enumerate(contingencies) if i ∉ get(contingency, :bus, [])]

    phi_sec = var(pm, nw)[:phi_sec] = JuMP.@variable(pm.model,
        [(i,c) in sec_buses], base_name="$(nw)_phi_sec",
        start = comp_start_value(ref(pm, nw, :bus, i), "phi_start")
    )

    if bounded
        for ((i, bus), c) in ref(pm, nw, :bus, contingencies)
            if i ∉ get(contingencies[c], :bus, [])
                JuMP.set_lower_bound(phi_sec[(i,c)], bus["vmin"] - 1.0)
                JuMP.set_upper_bound(phi_sec[(i,c)], bus["vmax"] - 1.0)
            end
        end
    end

    report && sol_component_value(pm, nw, :bus, :phi_sec, sec_buses, phi_sec)
end

""
function constraint_power_balance_sec(pm::AbstractLPACModel, n::Int, i::Int, c::ContingencyID, contingency::Contingency, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    phi_sec  = var(pm, n, :phi_sec, (i,c))
    p_sec    = get(var(pm, n),    :p_sec, Dict()); _check_var_keys(p_sec, bus_arcs, "active power", "branch")
    q_sec    = get(var(pm, n),    :q_sec, Dict()); _check_var_keys(q_sec, bus_arcs, "reactive power", "branch")
    pg_sec   = get(var(pm, n),   :pg_sec, Dict()); _check_var_keys(pg_sec, bus_gens, "active power", "generator")
    qg_sec   = get(var(pm, n),   :qg_sec, Dict()); _check_var_keys(qg_sec, bus_gens, "reactive power", "generator")
    ps   = get(var(pm, n),   :ps, Dict()); _check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(var(pm, n),   :qs, Dict()); _check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(var(pm, n),  :psw, Dict()); _check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(var(pm, n),  :qsw, Dict()); _check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    p_dc_sec = get(var(pm, n), :p_dc_sec, Dict()); _check_var_keys(p_dc_sec, bus_arcs_dc, "active power", "dcline")
    q_dc_sec = get(var(pm, n), :q_dc_sec, Dict()); _check_var_keys(q_dc_sec, bus_arcs_dc, "reactive power", "dcline")

    # not implemented
    @assert length(ps) == length(qs) == length(psw) == length(qsw) == 0

    # TODO these contingecies should be filtered out in advance
    if length([p_sec[((l,c),i,j)] for (l,i,j) in bus_arcs if l ∉ get(contingency, :branch, [])]) == 0
        return
    end

    cstr_p = JuMP.@constraint(pm.model,
        sum(p_sec[((l,c),i,j)] for (l,i,j) in bus_arcs if l ∉ get(contingency, :branch, []))
        + sum(p_dc_sec[((l_dc,c),i_dc,j_dc)] for (l_dc,i_dc,j_dc) in bus_arcs_dc if l_dc ∉ get(contingency, :dcline, []))
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        ==
        sum(pg_sec[(g, c)] for g in bus_gens if g ∉ get(contingency, :gen, []))
        - sum(ps_sec[(s,c)] for s in bus_storage if s ∉ get(contingency, :storage, []))
        - sum(pd for pd in values(bus_pd))
        - sum(gs for gs in values(bus_gs))*(1.0 + 2*phi_sec)
    )
    cstr_q = JuMP.@constraint(pm.model,
        sum(q_sec[((l,c),i,j)] for (l,i,j) in bus_arcs if l ∉ get(contingency, :branch, []))
        + sum(q_dc_sec[((l_dc,c),i_dc,j_dc)] for (l_dc,i_dc,j_dc) in bus_arcs_dc if l_dc ∉ get(contingency, :dcline, []))
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        ==
        sum(qg_sec[(g, c)] for g in bus_gens if g ∉ get(contingency, :gen, []))
        - sum(qs[(s,c)] for s in bus_storage if s ∉ get(contingency, :storage, []))
        - sum(qd for qd in values(bus_qd))
        + sum(bs for bs in values(bus_bs))*(1.0 + 2*phi_sec)
    )

    if _IM.report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


"`angmin <= z_branch[i]*(t[f_bus] - t[t_bus]) <= angmax` + slack linking"
function constraint_voltage_angle_difference_slack(pm::AbstractLPACCModel, n::Int, f_idx, angmin, angmax, vad_min, vad_max)
    i, f_bus, t_bus = f_idx
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)
    
    vad_slack_lower = var(pm, n, :vad_slack_lower, f_idx)
    vad_slack_upper = var(pm, n, :vad_slack_upper, f_idx)

    # slack linking
    JuMP.@constraint(pm.model, z*(va_fr - va_to) <= angmax + vad_slack_upper, base_name="$(n)_vad_slack_ub_$(i)")
    JuMP.@constraint(pm.model, z*(va_fr - va_to) >= angmin - vad_slack_lower, base_name="$(n)_vad_slack_lb_$(i)")
end

"`angmin <= z_branch[i]*(t[f_bus] - t[t_bus]) <= angmax` + slack linking"
function constraint_voltage_angle_difference_slack_on_off(pm::AbstractLPACCModel, n::Int, f_idx, angmin, angmax, vad_min, vad_max)
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

"`angmin*z_branch[i] + vad_min*(1-z_branch[i]) <= t[f_bus] - t[t_bus] <= angmax*z_branch[i] + vad_max*(1-z_branch[i])`
~ the voltage angle difference at each branch must lie between the branch-specific minimum and maximum (angmin, angmax) if the branch is active,
and between pre-determined minimum/maximum values (big M) if the branch is not active
"
function constraint_voltage_angle_difference_sec(pm::AbstractLPACCModel, n::Int, c::ContingencyID, contingency::Contingency, f_idx, angmin, angmax)
    (i, c), f_bus, t_bus = f_idx

    va_fr_sec = var(pm, n, :va_sec, (f_bus,c))
    va_to_sec = var(pm, n, :va_sec, (t_bus,c))

    JuMP.@constraint(pm.model, va_fr_sec - va_to_sec <= angmax, base_name="$(n)_C$(c)_vad_ub_$(i)")
    JuMP.@constraint(pm.model, va_fr_sec - va_to_sec >= angmin, base_name="$(n)_C$(c)_vad_lb_$(i)")
end


"`angmin*z_branch[i] + vad_min*(1-z_branch[i]) <= t[f_bus] - t[t_bus] <= angmax*z_branch[i] + vad_max*(1-z_branch[i])`
~ the voltage angle difference at each branch must lie between the branch-specific minimum and maximum (angmin, angmax) if the branch is active,
and between pre-determined minimum/maximum values (big M) if the branch is not active
"
function constraint_voltage_angle_difference_sec_on_off(pm::AbstractLPACCModel, n::Int, c::ContingencyID, contingency::Contingency, f_idx, angmin, angmax, vad_min, vad_max)
    (i, c), f_bus, t_bus = f_idx

    va_fr_sec = var(pm, n, :va_sec, (f_bus,c))
    va_to_sec = var(pm, n, :va_sec, (t_bus,c))
    z = var(pm, n, :z_branch, i)

    JuMP.@constraint(pm.model, va_fr_sec - va_to_sec <= angmax*z + vad_max*(1-z))
    JuMP.@constraint(pm.model, va_fr_sec - va_to_sec >= angmin*z + vad_min*(1-z))
end

"`angmin <= z_branch[i]*(t[f_bus] - t[t_bus]) <= angmax` +  slack linking"
function constraint_voltage_angle_difference_slack_sec(pm::AbstractLPACCModel, n::Int, c::ContingencyID, contingency::Contingency, f_idx, angmin, angmax, vad_min, vad_max)
    (i, c), f_bus, t_bus = f_idx

    va_sec_fr = var(pm, n, :va_sec, (f_bus,c))
    va_sec_to = var(pm, n, :va_sec, (t_bus,c))

    vad_slack_lower = var(pm, n, :vad_slack_lower, (i, f_bus, t_bus))
    vad_slack_upper = var(pm, n, :vad_slack_upper, (i, f_bus, t_bus))

     # slack linking
    JuMP.@constraint(pm.model, z*(va_sec_fr - va_sec_to) <= angmax + vad_slack_upper, base_name="$(n)_C$(c)_vad_slack_ub_$(i)")
    JuMP.@constraint(pm.model, z*(va_sec_fr - va_sec_to) >= angmin - vad_slack_lower, base_name="$(n)_C$(c)_vad_slack_ub_$(i)")
end

"`angmin <= z_branch[i]*(t[f_bus] - t[t_bus]) <= angmax` +  slack linking"
function constraint_voltage_angle_difference_slack_sec_on_off(pm::AbstractLPACCModel, n::Int, c::ContingencyID, contingency::Contingency, f_idx, angmin, angmax, vad_min, vad_max)
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

""
function constraint_model_voltage_sec_on_off(pm::AbstractLPACCModel, n::Int, contingencies::Vector{Contingency})
    phi  = var(pm, n, :phi_sec)
    t = var(pm, n, :va_sec)
    phi_fr = var(pm, n, :phi_fr_sec)
    phi_to = var(pm, n, :phi_to_sec)

    td = var(pm, n, :td_sec)
    cs = var(pm, n, :cs_sec)

    z = var(pm, n, :z_branch)

    td_lb = ref(pm, n, :off_angmin)
    td_ub = ref(pm, n, :off_angmax)
    td_max = max(abs(td_lb), abs(td_ub))


    for (c, contingency) in enumerate(contingencies)
        for (l, branch) in ref(pm, n, :branch)
            if l ∉ get(contingency, :branch, [])
                i = branch["f_bus"]
                j = branch["t_bus"]
                
                JuMP.@constraint(pm.model, t[(i,c)] - t[(j,c)] >= td[(l,c)] + td_lb*(1-z[l]))
                JuMP.@constraint(pm.model, t[(i,c)] - t[(j,c)] <= td[(l,c)] + td_ub*(1-z[l]))
                relaxation_cos_on_off(pm.model, td[(l,c)], cs[(l,c)], z[l], td_max)

                _IM.constraint_bounds_on_off(pm.model, td[(l,c)], z[l])
                _IM.constraint_bounds_on_off(pm.model, phi_fr[(l,c)], z[l])
                _IM.constraint_bounds_on_off(pm.model, phi_to[(l,c)], z[l])
                _IM.relaxation_equality_on_off(pm.model, phi[(i,c)], phi_fr[(l,c)], z[l])
                _IM.relaxation_equality_on_off(pm.model, phi[(j,c)], phi_to[(l,c)], z[l])
            end
        end
    end
end

""
function variable_bus_voltage_sec_on_off(pm::AbstractLPACCModel, contingencies::Vector{Contingency}; kwargs...)
    variable_bus_voltage_angle_sec(pm, contingencies; kwargs...)
    variable_bus_voltage_magnitude_sec(pm, contingencies; kwargs...)

    variable_branch_voltage_magnitude_fr_sec_on_off(pm, contingencies; kwargs...)
    variable_branch_voltage_magnitude_to_sec_on_off(pm, contingencies; kwargs...)
    variable_branch_voltage_product_angle_sec_on_off(pm, contingencies; kwargs...)
    variable_branch_cosine_sec_on_off(pm, contingencies; kwargs...)
end


""
function variable_branch_voltage_magnitude_fr_sec_on_off(pm::AbstractLPACCModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, report::Bool=true)
    buses = ref(pm, nw, :bus)
    branches = ref(pm, nw, :branch)
    sec_branches = [(i,c) for i in ids(pm, nw, :branch) for (c, contingency) in enumerate(contingencies) if i ∉ get(contingency, :branch, [])]

    phi_fr_sec = var(pm, nw)[:phi_fr_sec] = JuMP.@variable(pm.model,
        [(i,c) in sec_branches], base_name="$(nw)_phi_fr_sec",
        lower_bound = min(0, buses[branches[i]["f_bus"]]["vmin"] - 1.0),
        upper_bound = max(0, buses[branches[i]["f_bus"]]["vmax"] - 1.0),
        start = comp_start_value(ref(pm, nw, :bus, branches[i]["f_bus"]), "phi_fr_start")
    )

    report && sol_component_value(pm, nw, :branch, :phi_fr_sec, sec_branches, phi_fr_sec)
end

""
function variable_branch_voltage_magnitude_to_sec_on_off(pm::AbstractLPACCModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, report::Bool=true)
    buses = ref(pm, nw, :bus)
    branches = ref(pm, nw, :branch)
    sec_branches = [(i,c) for i in ids(pm, nw, :branch) for (c, contingency) in enumerate(contingencies) if i ∉ get(contingency, :branch, [])]
    
    phi_to_sec = var(pm, nw)[:phi_to_sec] = JuMP.@variable(pm.model,
        [(i,c) in sec_branches], base_name="$(nw)_phi_to",
        lower_bound = min(0, buses[branches[i]["t_bus"]]["vmin"] - 1.0),
        upper_bound = max(0, buses[branches[i]["t_bus"]]["vmax"] - 1.0),
        start = comp_start_value(ref(pm, nw, :bus, branches[i]["t_bus"]), "phi_to_start")
    )

    report && sol_component_value(pm, nw, :branch, :phi_to_sec, sec_branches, phi_to_sec)
end

""
function constraint_ohms_yt_from_sec(pm::AbstractLPACCModel, n::Int, i, c::ContingencyID, contingency::Contingency, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_fr_sec = var(pm, n,    :p_sec, f_idx)
    q_fr_sec = var(pm, n,    :q_sec, f_idx)
    phi_fr_sec = var(pm, n, :phi_fr_sec, (i,c))
    phi_to_sec = var(pm, n, :phi_to_sec, (i,c))
    td_sec = var(pm, n, :td_sec, (i,c))
    cs_sec = var(pm, n, :cs_sec, (i,c))

    JuMP.@constraint(pm.model, p_fr_sec ==  (g+g_fr)/tm^2*(2*phi_fr_sec) + (-g*tr+b*ti)/tm^2*(cs_sec + phi_fr_sec + phi_to_sec) + (-b*tr-g*ti)/tm^2*(td_sec), base_name="$(n)_C$(c)_ohm_p_fr_$(i)")
    JuMP.@constraint(pm.model, q_fr_sec == -(b+b_fr)/tm^2*(2*phi_fr_sec) - (-b*tr-g*ti)/tm^2*(cs_sec + phi_fr_sec + phi_to_sec) + (-g*tr+b*ti)/tm^2*(td_sec), base_name="$(n)_C$(c)_ohm_q_fr_$(i)")
end

""
function constraint_ohms_yt_from_sec_on_off(pm::AbstractLPACCModel, n::Int, i, c::ContingencyID, contingency::Contingency, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max)
    p_fr_sec = var(pm, n,    :p_sec, f_idx)
    q_fr_sec = var(pm, n,    :q_sec, f_idx)
    phi_fr_sec = var(pm, n, :phi_fr_sec, (i,c))
    phi_to_sec = var(pm, n, :phi_to_sec, (i,c))
    td_sec = var(pm, n, :td_sec, (i,c))
    cs_sec = var(pm, n, :cs_sec, (i,c))
    z = var(pm, n, :z_branch, i)

    JuMP.@constraint(pm.model, p_fr_sec ==  (g+g_fr)/tm^2*(z + 2*phi_fr_sec) + (-g*tr+b*ti)/tm^2*(cs_sec + phi_fr_sec + phi_to_sec) + (-b*tr-g*ti)/tm^2*(td_sec))
    JuMP.@constraint(pm.model, q_fr_sec == -(b+b_fr)/tm^2*(z + 2*phi_fr_sec) - (-b*tr-g*ti)/tm^2*(cs_sec + phi_fr_sec + phi_to_sec) + (-g*tr+b*ti)/tm^2*(td_sec))
end

""
function constraint_ohms_yt_to_sec(pm::AbstractLPACCModel, n::Int, i, c::ContingencyID, contingency::Contingency, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    p_to_sec = var(pm, n,    :p_sec, t_idx)
    q_to_sec = var(pm, n,    :q_sec, t_idx)
    phi_fr_sec = var(pm, n, :phi_fr_sec, (i,c))
    phi_to_sec = var(pm, n, :phi_to_sec, (i,c))
    td_sec = var(pm, n, :td_sec, (i,c))
    cs_sec = var(pm, n, :cs_sec, (i,c))

    JuMP.@constraint(pm.model, p_to_sec ==  (g+g_to)*(2*phi_to_sec) + (-g*tr-b*ti)/tm^2*(cs_sec + phi_fr_sec + phi_to_sec) + (-b*tr+g*ti)/tm^2*-(td_sec), base_name="$(n)_C$(c)_ohm_p_to_$(i)" )
    JuMP.@constraint(pm.model, q_to_sec == -(b+b_to)*(2*phi_to_sec) - (-b*tr+g*ti)/tm^2*(cs_sec + phi_fr_sec + phi_to_sec) + (-g*tr-b*ti)/tm^2*-(td_sec), base_name="$(n)_C$(c)_ohm_q_to_$(i)" )
end

""
function constraint_ohms_yt_to_sec_on_off(pm::AbstractLPACCModel, n::Int, i, c::ContingencyID, contingency::Contingency, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max)
    p_to_sec = var(pm, n,    :p_sec, t_idx)
    q_to_sec = var(pm, n,    :q_sec, t_idx)
    phi_fr_sec = var(pm, n, :phi_fr_sec, (i,c))
    phi_to_sec = var(pm, n, :phi_to_sec, (i,c))
    td_sec = var(pm, n, :td_sec, (i,c))
    cs_sec = var(pm, n, :cs_sec, (i,c))
    z = var(pm, n, :z_branch, i)

    JuMP.@constraint(pm.model, p_to_sec ==  (g+g_to)*(z + 2*phi_to_sec) + (-g*tr-b*ti)/tm^2*(cs_sec + phi_fr_sec + phi_to_sec) + (-b*tr+g*ti)/tm^2*-(td_sec) )
    JuMP.@constraint(pm.model, q_to_sec == -(b+b_to)*(z + 2*phi_to_sec) - (-b*tr+g*ti)/tm^2*(cs_sec + phi_fr_sec + phi_to_sec) + (-g*tr-b*ti)/tm^2*-(td_sec) )
end