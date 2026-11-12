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

"variable: `p[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_rate_a_slack(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    # aux branches don't have slack
    flow_arcs = [(l,i,j) for (l,i,j) in ref(pm, nw, :arcs) if get(ref(pm, nw, :branch, l), "type", "") ∉ ["breaker"]]

    rate_a_slack = var(pm, nw)[:rate_a_slack] = JuMP.@variable(pm.model,
        [(l,i,j) in flow_arcs], base_name="$(nw)_rate_a_slack",
        start = comp_start_value(ref(pm, nw, :branch, l), "rate_a_slack_start")
    )

    if bounded
        for arc in flow_arcs
            JuMP.set_lower_bound(rate_a_slack[arc], 0)
        end
    end

    for (l,branch) in ref(pm, nw, :branch)
        if haskey(branch, "raf_s_start")
            f_idx = (l, branch["f_bus"], branch["t_bus"])
            JuMP.set_start_value(rate_a_slack[f_idx], branch["raf_s_start"])
        end
        if haskey(branch, "rat_s_start")
            t_idx = (l, branch["t_bus"], branch["f_bus"])
            JuMP.set_start_value(rate_a_slack[t_idx], branch["rat_s_start"])
        end
    end

    if get(pm.data, "FIX_SLACK", false) && ! get(pm.data, "SUBPROBLEM", false)
        for (l,branch) in ref(pm, nw, :branch)
            f_idx = (l, branch["f_bus"], branch["t_bus"])
            t_idx = (l, branch["t_bus"], branch["f_bus"])
            JuMP.fix(rate_a_slack[f_idx], JuMP.start_value(rate_a_slack[f_idx]), force=true)
            JuMP.fix(rate_a_slack[t_idx], JuMP.start_value(rate_a_slack[t_idx]), force=true)
        end
    end

    flow_arcs_from = [(l,i,j) for (l,i,j) in ref(pm, nw, :arcs_from) if get(ref(pm, nw, :branch, l), "type", "") ∉ ["breaker"]]
    flow_arcs_to   = [(l,i,j) for (l,i,j) in ref(pm, nw, :arcs_to)   if get(ref(pm, nw, :branch, l), "type", "") ∉ ["breaker"]]

    report && sol_component_value_edge(pm, nw, :branch, :raf_s, :rat_s, flow_arcs_from, flow_arcs_to, rate_a_slack)
end

"variable: `p[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_rate_a_slack_sec(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    # aux branches don't have slack
    sec_arcs = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs) for (c, contingency) in enumerate(contingencies) if (l ∉ get(contingency, :branch, []) && get(ref(pm, nw, :branch, l), "type", "") ∉ ["breaker"])]

    if ! PARAMS[:CONTINGENCIES][:SLACK_EQUALS_SEC_SLACK]
        rate_a_slack_sec = var(pm, nw)[:rate_a_slack_sec] = JuMP.@variable(pm.model,
        [((l,c),i,j) in sec_arcs], base_name="$(nw)_rate_a_slack_sec",
        start = 0
        )
    else
        var(pm, nw)[:rate_a_slack_sec] = Dict(((l,c),i,j) => var(pm, nw, :rate_a_slack, (l,i,j)) for ((l,c),i,j) in sec_arcs)
        rate_a_slack_sec = var(pm, nw)[:rate_a_slack_sec]
    end


    if bounded
        for arc in sec_arcs
            JuMP.set_lower_bound(rate_a_slack_sec[arc], 0)
        end
    end

    flow_arcs_from = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs_from) for (c, contingency) in enumerate(contingencies) if (l ∉ get(contingency, :branch, []) && get(ref(pm, nw, :branch, l), "type", "") ∉ ["breaker"])]
    flow_arcs_to   = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs_to)   for (c, contingency) in enumerate(contingencies) if (l ∉ get(contingency, :branch, []) && get(ref(pm, nw, :branch, l), "type", "") ∉ ["breaker"])]

    report && sol_component_value_edge(pm, nw, :branch, :raf_s_sec, :rat_s_sec, flow_arcs_from, flow_arcs_to, rate_a_slack_sec)
end

""
function variable_rate_a_slack_on_off(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    # aux branches don't have slack
    flow_arcs = [(l,i,j) for (l,i,j) in ref(pm, nw, :arcs) if get(ref(pm, nw, :branch, l), "type", "") ∉ ["breaker"]]

    rate_a_slack_on_off = var(pm, nw)[:rate_a_slack_on_off] = JuMP.@variable(pm.model,
        [(l,i,j) in flow_arcs], base_name="$(nw)_rate_a_slack_on_off",
        binary = true,
        start = comp_start_value(ref(pm, nw, :branch, l), "rate_a_slack_on_off_start")
    )

    flow_arcs_from = [(l,i,j) for (l,i,j) in ref(pm, nw, :arcs_from) if get(ref(pm, nw, :branch, l), "type", "") ∉ ["breaker"]]
    flow_arcs_to   = [(l,i,j) for (l,i,j) in ref(pm, nw, :arcs_to)   if get(ref(pm, nw, :branch, l), "type", "") ∉ ["breaker"]]

    report && sol_component_value_edge(pm, nw, :branch, :raf_s_b, :rat_s_b, flow_arcs_from, flow_arcs_to, rate_a_slack_on_off)
end

function variable_rate_a_slack_sec_on_off(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    # aux branches don't have slack
    sec_arcs = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs) for (c, contingency) in enumerate(contingencies) if (l ∉ get(contingency, :branch, []) && get(ref(pm, nw, :branch, l), "type", "") ∉ ["breaker"])]

    rate_a_slack_sec_on_off = var(pm, nw)[:rate_a_slack_on_off] = JuMP.@variable(pm.model,
        [(l,i,j) in flow_arcs], base_name="$(nw)_rate_a_slack_sec_on_off",
        binary = true,
        start = 0
    )

    flow_arcs_from = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs_from) for (c, contingency) in enumerate(contingencies) if (l ∉ get(contingency, :branch, []) && get(ref(pm, nw, :branch, l), "type", "") ∉ ["breaker"])]
    flow_arcs_to   = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs_to)   for (c, contingency) in enumerate(contingencies) if (l ∉ get(contingency, :branch, []) && get(ref(pm, nw, :branch, l), "type", "") ∉ ["breaker"])]

    report && sol_component_value_edge(pm, nw, :branch, :raf_s_b, :rat_s_b, flow_arcs_from, flow_arcs_to, rate_a_slack_sec_on_off)
end

"variable: `t[i]` for `i` in `bus`es"
function variable_bus_voltage_angle_sec(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    sec_buses = [(i,c) for i in ids(pm, nw, :bus) for (c, contingency) in enumerate(contingencies) if i ∉ get(contingency, :bus, [])]

    va_sec = var(pm, nw)[:va_sec] = JuMP.@variable(pm.model,
        [(i,c) in sec_buses], base_name="$(nw)_va_sec",
        start = comp_start_value(ref(pm, nw, :bus, i), "va_start")
    )

    report && sol_component_value(pm, nw, :bus, :va_sec, sec_buses, va_sec)
end

"real part of the voltage variable `i` in `bus`es"
function variable_bus_voltage_real_sec(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    sec_buses = [(i,c) for i in ids(pm, nw, :bus) for (c, contingency) in enumerate(contingencies) if i ∉ get(contingency, :bus, [])]
    vr_sec = var(pm, nw)[:vr_sec] = JuMP.@variable(pm.model,
        [(i, c) in sec_buses], base_name="$(nw)_vr_sec",
        start = comp_start_value(ref(pm, nw, :bus, i), "vr_sec_start", 1.0)
    )

    if bounded
        for c in eachindex(contingencies)
            for (i, bus) in ref(pm, nw, :bus)
                JuMP.set_lower_bound(vr_sec[(i,c)], -bus["vmax"])
                JuMP.set_upper_bound(vr_sec[(i,c)],  bus["vmax"])
            end
        end
    end

    report && sol_component_value(pm, nw, :bus, :vr_sec, sec_buses, vr_sec)
end

"real part of the voltage variable `i` in `bus`es"
function variable_bus_voltage_imaginary_sec(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    sec_buses = [(i,c) for i in ids(pm, nw, :bus) for (c, contingency) in enumerate(contingencies) if i ∉ get(contingency, :bus, [])]
    vi_sec = var(pm, nw)[:vi_sec] = JuMP.@variable(pm.model,
        [(i, c) in sec_buses], base_name="$(nw)_vi_sec",
        start = comp_start_value(ref(pm, nw, :bus, i), "vi_sec_start", 1.0)
    )

    if bounded
        for c in eachindex(contingencies)
            for (i, bus) in ref(pm, nw, :bus)
                JuMP.set_lower_bound(vi_sec[(i,c)], -bus["vmax"])
                JuMP.set_upper_bound(vi_sec[(i,c)],  bus["vmax"])
            end
        end
    end

    report && sol_component_value(pm, nw, :bus, :vi_sec, sec_buses, vi_sec)
end

"variable: `v[i]` for `i` in `bus`es"
function variable_bus_voltage_magnitude_sec(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    sec_buses = [(i,c) for i in ids(pm, nw, :bus) for (c, contingency) in enumerate(contingencies) if i ∉ get(contingency, :bus, [])]

    vm_sec = var(pm, nw)[:vm_sec] = JuMP.@variable(pm.model,
        [(i, c) in sec_buses], base_name="$(nw)_vm_sec",
        start = comp_start_value(ref(pm, nw, :bus, i), "vm_start", 1.0)
    )

    if bounded
        for ((i, bus), c) in ref(pm, nw, :bus, contingencies)
            if i ∉ get(contingencies[c], :bus, [])
                JuMP.set_lower_bound(vm_sec[(i,c)], bus["vmin"])
                JuMP.set_upper_bound(vm_sec[(i,c)], bus["vmax"])
            end
        end
    end

    report && sol_component_value(pm, nw, :bus, :vm_sec, sec_buses, vm_sec)
end

"variable: `0 <= vm_fr[l] <= buses[branches[l][\"f_bus\"]][\"vmax\"]` for `l` in `branch`es"
function variable_branch_voltage_magnitude_fr_sec_on_off(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, report::Bool=true)
    buses = ref(pm, nw, :bus)
    branches = ref(pm, nw, :branch)
    sec_branches = [(i,c) for i in ids(pm, nw, :branch) for (c, contingency) in enumerate(contingencies) if i ∉ get(contingency, :branch, [])]

    vm_fr_sec = var(pm, nw)[:vm_fr_sec] = JuMP.@variable(pm.model,
        [(i,c) in sec_branches], base_name="$(nw)_vm_fr_sec",
        lower_bound = 0,
        upper_bound = buses[branches[i]["f_bus"]]["vmax"],
        start = comp_start_value(ref(pm, nw, :branch, i), "vm_fr_start", 1.0)
    )

    report && sol_component_value(pm, nw, :branch, :vm_fr_sec, sec_branches, vm_fr_sec)
end

"variable: `0 <= vm_to[l] <= buses[branches[l][\"t_bus\"]][\"vmax\"]` for `l` in `branch`es"
function variable_branch_voltage_magnitude_to_sec_on_off(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, report::Bool=true)
    buses = ref(pm, nw, :bus)
    branches = ref(pm, nw, :branch)
    sec_branches = [(i,c) for i in ids(pm, nw, :branch) for (c, contingency) in enumerate(contingencies) if i ∉ get(contingency, :branch, [])]

    vm_to_sec = var(pm, nw)[:vm_to_sec] = JuMP.@variable(pm.model,
        [(i,c) in sec_branches], base_name="$(nw)_vm_to_sec",
        lower_bound = 0,
        upper_bound = buses[branches[i]["t_bus"]]["vmax"],
        start = comp_start_value(ref(pm, nw, :branch, i), "vm_to_start", 1.0)
    )

    report && sol_component_value(pm, nw, :branch, :vm_to_sec, sec_branches, vm_to_sec)
end

"variable: `w[i] >= 0` for `i` in `bus`es"
function variable_bus_voltage_magnitude_sqr_sec(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    sec_buses = [(i,c) for i in ids(pm, nw, :bus) for (c, contingency) in enumerate(contingencies) if i ∉ get(contingency, :bus, [])]

    w_sec = var(pm, nw)[:w_sec] = JuMP.@variable(pm.model,
        [(i,c) in sec_buses], base_name="$(nw)_w_sec",
        lower_bound = 0.0,
        start = comp_start_value(ref(pm, nw, :bus, i), "w_start", 1.001)
    )

    if bounded
        for ((i, bus), c) in ref(pm, nw, :bus, contingencies)
            JuMP.set_lower_bound(w_sec[(i,c)], bus["vmin"]^2)
            JuMP.set_upper_bound(w_sec[(i,c)], bus["vmax"]^2)
        end
    end

    report && sol_component_value(pm, nw, :bus, :w_sec, sec_buses, w_sec)
end

"variable: `0 <= w_fr[l] <= buses[branches[l][\"f_bus\"]][\"vmax\"]^2` for `l` in `branch`es"
function variable_branch_voltage_magnitude_fr_sqr_sec_on_off(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, report::Bool=true)
    buses = ref(pm, nw, :bus)
    branches = ref(pm, nw, :branch)
    sec_branches = [(i,c) for i in ids(pm, nw, :branch) for (c, contingency) in enumerate(contingencies) if i ∉ get(contingency, :branch, [])]

    w_fr_sec = var(pm, nw)[:w_fr_sec] = JuMP.@variable(pm.model,
        [(i,c) in sec_branches], base_name="$(nw)_w_fr_sec",
        lower_bound = 0,
        upper_bound = buses[branches[i]["f_bus"]]["vmax"]^2,
        start = comp_start_value(ref(pm, nw, :branch, i), "w_fr_start", 1.001)
    )

    report && sol_component_value(pm, nw, :branch, :w_fr_sec, sec_branches, w_fr_sec)
end

"variable: `0 <= w_to[l] <= buses[branches[l][\"t_bus\"]][\"vmax\"]^2` for `l` in `branch`es"
function variable_branch_voltage_magnitude_to_sqr_sec_on_off(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, report::Bool=true)
    buses = ref(pm, nw, :bus)
    branches = ref(pm, nw, :branch)
    sec_branches = [(i,c) for i in ids(pm, nw, :branch) for (c, contingency) in enumerate(contingencies) if i ∉ get(contingency, :branch, [])]

    w_to_sec = var(pm, nw)[:w_to_sec] = JuMP.@variable(pm.model,
        [(i,c) in sec_branches], base_name="$(nw)_w_to_sec",
        lower_bound = 0,
        upper_bound = buses[branches[i]["t_bus"]]["vmax"]^2,
        start = comp_start_value(ref(pm, nw, :branch, i), "w_to_start", 1.001)
    )

    report && sol_component_value(pm, nw, :branch, :w_to_sec, sec_branches, w_to_sec)
end

""
function variable_branch_voltage_product_sec_on_off(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, report::Bool=true)
    wr_min, wr_max, wi_min, wi_max = ref_calc_voltage_product_bounds(ref(pm, nw, :buspairs))
    bi_bp = Dict((i, (b["f_bus"], b["t_bus"])) for (i,b) in ref(pm, nw, :branch))
    sec_branches = [(i,c) for i in ids(pm, nw, :branch) for (c, contingency) in enumerate(contingencies) if i ∉ get(contingency, :branch, [])]

    wr_sec = var(pm, nw)[:wr_sec] = JuMP.@variable(pm.model,
        [(b,c) in sec_branches], base_name="$(nw)_wr_sec",
        lower_bound = min(0, wr_min[bi_bp[b]]),
        upper_bound = max(0, wr_max[bi_bp[b]]),
        start = comp_start_value(ref(pm, nw, :buspairs, bi_bp[b]), "wr_start", 1.0)
    )
    wi_sec = var(pm, nw)[:wi_sec] = JuMP.@variable(pm.model,
        [(b,c) in sec_branches], base_name="$(nw)_wi_sec",
        lower_bound = min(0, wi_min[bi_bp[b]]),
        upper_bound = max(0, wi_max[bi_bp[b]]),
        start = comp_start_value(ref(pm, nw, :buspairs, bi_bp[b]), "wi_start")
    )

    report && sol_component_value(pm, nw, :branch, :wr_sec, sec_branches, wr_sec)
    report && sol_component_value(pm, nw, :branch, :wi_sec, sec_branches, wi_sec)
end

"slack variables for voltage angle difference"
function variable_vad_slack(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    if pm.data["BUSBARS"]
        # aux branches don't have slack
        flow_arcs = [(l,i,j) for (l,i,j) in ref(pm, nw, :arcs) if get(ref(pm, nw, :branch, l), "type", "") ∉ ["breaker", "K"] || ref(pm, nw, :branch, l)["transformer"]]
    else
        flow_arcs = ref(pm, nw, :arcs)
    end
    
    vad_slack_lower = var(pm, nw)[:vad_slack_lower] = JuMP.@variable(pm.model,
        [(l,i,j) in flow_arcs], base_name="$(nw)_vad_slack_lower",
        start = 0
    )

    vad_slack_upper = var(pm, nw)[:vad_slack_upper] = JuMP.@variable(pm.model,
        [(l,i,j) in flow_arcs], base_name="$(nw)_vad_slack_upper",
        start = 0
    )

    if bounded
        for arc in ref(pm, nw, :arcs)
            
            JuMP.set_lower_bound(vad_slack_lower[arc], 0)
            JuMP.set_lower_bound(vad_slack_upper[arc], 0)
        end
    end

    if pm.data["BUSBARS"]
        flow_arcs_from = [(l,i,j) for (l,i,j) in ref(pm, nw, :arcs_from) if get(ref(pm, nw, :branch, l), "type", "") ∉ ["breaker", "K"] || ref(pm, nw, :branch, l)["transformer"]]
        flow_arcs_to   = [(l,i,j) for (l,i,j) in ref(pm, nw, :arcs_to)   if get(ref(pm, nw, :branch, l), "type", "") ∉ ["breaker", "K"] || ref(pm, nw, :branch, l)["transformer"]]
    else
        flow_arcs_from = ref(pm, nw, :arcs_from)
        flow_arcs_to   = ref(pm, nw, :arcs_to)
    end

    report && sol_component_value_edge(pm, nw, :branch, :vadf_sl, :vadt_sl, flow_arcs_from, flow_arcs_to, vad_slack_lower)
    report && sol_component_value_edge(pm, nw, :branch, :vadf_su, :vadt_su, flow_arcs_from, flow_arcs_to, vad_slack_upper)
end

"slack variables for voltage angle difference"
function variable_vad_slack_on_off(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    if pm.data["BUSBARS"]
        # aux branches don't have slack
        flow_arcs = [(l,i,j) for (l,i,j) in ref(pm, nw, :arcs) if get(ref(pm, nw, :branch, l), "type", "") ∉ ["breaker", "K"] || ref(pm, nw, :branch, l)["transformer"]]
    else
        flow_arcs = ref(pm, nw, :arcs)
    end

    vad_slack_lower_on_off = var(pm, nw)[:vad_slack_lower_on_off] = JuMP.@variable(pm.model,
        [(l,i,j) in flow_arcs], base_name="$(nw)_vad_slack_lower_on_off",
        binary=true,
        start = 0
    )

    vad_slack_upper_on_off = var(pm, nw)[:vad_slack_upper_on_off] = JuMP.@variable(pm.model,
        [(l,i,j) in flow_arcs], base_name="$(nw)_vad_slack_upper_on_off",
        binary=true,
        start = 0
    )

    if pm.data["BUSBARS"]
        flow_arcs_from = [(l,i,j) for (l,i,j) in ref(pm, nw, :arcs_from) if get(ref(pm, nw, :branch, l), "type", "") ∉ ["breaker", "K"] || ref(pm, nw, :branch, l)["transformer"]]
        flow_arcs_to   = [(l,i,j) for (l,i,j) in ref(pm, nw, :arcs_to)   if get(ref(pm, nw, :branch, l), "type", "") ∉ ["breaker", "K"] || ref(pm, nw, :branch, l)["transformer"]]
    else
        flow_arcs_from = ref(pm, nw, :arcs_from)
        flow_arcs_to   = ref(pm, nw, :arcs_to)
    end

    report && sol_component_value_edge(pm, nw, :branch, :vadf_sl_b, :vadt_sl_b, flow_arcs_from, flow_arcs_to, vad_slack_lower_on_off)
    report && sol_component_value_edge(pm, nw, :branch, :vadf_su_b, :vadt_su_b, flow_arcs_from, flow_arcs_to, vad_slack_upper_on_off)
end

"generates variables for both `active` and `reactive` generation"
function variable_gen_power_slack(pm::AbstractPowerModel; kwargs...)
    variable_gen_power_real_slack(pm; kwargs...)
    variable_gen_power_imaginary_slack(pm; kwargs...)
end

"variable: `pg[j]` for `j` in `gen`"
function variable_gen_power_real_slack(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    pmin_slack = var(pm, nw)[:pmin_slack] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :gen)], base_name="$(nw)_pmin_slack",
        start = 0
    )
    pmax_slack = var(pm, nw)[:pmax_slack] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :gen)], base_name="$(nw)_pmax_slack",
        start = 0
    )

    if bounded
        for (i, gen) in ref(pm, nw, :gen)
            JuMP.set_lower_bound(pmin_slack[i], 0)
            JuMP.set_lower_bound(pmax_slack[i], 0)
            
            if gen["pmin"] < 0
                # pmin can be decreased arbitrarily, pmax can be increased arbitrarily
            else
                # pmin can only be reduced to 0, pmax can be increased arbitrarily
                JuMP.set_upper_bound(pmin_slack[i], gen["pmin"])
            end
        end
    end

    report && sol_component_value(pm, nw, :gen, :pmin_slack, ids(pm, nw, :gen), pmin_slack)
    report && sol_component_value(pm, nw, :gen, :pmax_slack, ids(pm, nw, :gen), pmax_slack)
end

function variable_pd_slack(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    pd_slack_lower = var(pm, nw)[:pd_slack_lower] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :load)], base_name="$(nw)_pd_slack_lower",
        start = 0,
        lower_bound = 0
    )

    pd_slack_upper = var(pm, nw)[:pd_slack_upper] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :load)], base_name="$(nw)_pd_slack_upper",
        start = 0,
        lower_bound = 0
    )
    
    for (i, load) in ref(pm, nw, :load)

        pd_val = load["pd"]
        max_load = abs(pd_val)

        if pd_val > 0
            JuMP.set_upper_bound(pd_slack_lower[i], max_load)
            JuMP.fix(pd_slack_upper[i], 0.0; force=true)
        elseif pd_val < 0
            JuMP.set_upper_bound(pd_slack_upper[i], max_load)
            JuMP.fix(pd_slack_lower[i], 0.0; force=true)
        else
            JuMP.fix(pd_slack_lower[i], 0.0; force=true)
            JuMP.fix(pd_slack_upper[i], 0.0; force=true)
        end
    end

    report && sol_component_value(pm, nw, :load, :pd_slack_lower, ids(pm, nw, :load), pd_slack_lower)
    report && sol_component_value(pm, nw, :load, :pd_slack_upper, ids(pm, nw, :load), pd_slack_upper)
end

function variable_qd_slack(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    qd_slack_lower = var(pm, nw)[:qd_slack_lower] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :load)], base_name="$(nw)_qd_slack_lower",
        start = 0,
        lower_bound = 0
    )

    qd_slack_upper = var(pm, nw)[:qd_slack_upper] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :load)], base_name="$(nw)_qd_slack_upper",
        start = 0,
        lower_bound = 0
    )
    
    for (i, load) in ref(pm, nw, :load)

        pd_val = load["pd"]
        max_load = abs(pd_val)

        if pd_val > 0
            JuMP.set_upper_bound(qd_slack_lower[i], max_load)
            JuMP.fix(qd_slack_upper[i], 0.0; force=true)
        elseif pd_val < 0
            JuMP.set_upper_bound(qd_slack_upper[i], max_load)
            JuMP.fix(qd_slack_lower[i], 0.0; force=true)
        else
            JuMP.fix(qd_slack_lower[i], 0.0; force=true)
            JuMP.fix(qd_slack_upper[i], 0.0; force=true)
        end
    end

    report && sol_component_value(pm, nw, :load, :qd_slack_lower, ids(pm, nw, :load), qd_slack_lower)
    report && sol_component_value(pm, nw, :load, :qd_slack_upper, ids(pm, nw, :load), qd_slack_upper)
end

function variable_pd_slack_sec(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    sec_loads = [(i, c) for i in ids(pm, nw, :load) for (c, _) in enumerate(contingencies)]
    pd_slack_lower_sec = var(pm, nw)[:pd_slack_lower_sec] = JuMP.@variable(pm.model,
        [(i, c) in sec_loads], base_name="$(nw)_pd_slack_lower_sec",
        start = 0,
        lower_bound = 0
    )

    pd_slack_upper_sec = var(pm, nw)[:pd_slack_upper_sec] = JuMP.@variable(pm.model,
        [(i, c) in sec_loads], base_name="$(nw)_pd_slack_upper_sec",
        start = 0,
        lower_bound = 0
    )

    for (i, load) in ref(pm, nw, :load)

        pd_val = load["pd"]
        max_load = abs(pd_val)

        for (c, _) in enumerate(contingencies)
            if pd_val > 0
                JuMP.set_upper_bound(pd_slack_lower_sec[(i,c)], max_load)
                JuMP.fix(pd_slack_upper_sec[(i,c)], 0.0; force=true)
            elseif pd_val < 0
                JuMP.set_upper_bound(pd_slack_upper_sec[(i,c)], max_load)
                JuMP.fix(pd_slack_lower_sec[(i,c)], 0.0; force=true)
            else
                JuMP.fix(pd_slack_lower_sec[(i,c)], 0.0; force=true)
                JuMP.fix(pd_slack_upper_sec[(i,c)], 0.0; force=true)
            end
        end
    end

    report && sol_component_value(pm, nw, :load, :pd_slack_lower_sec, [(i,c) for i in ids(pm, nw, :load) for (c, _) in enumerate(contingencies)], pd_slack_lower_sec)
    report && sol_component_value(pm, nw, :load, :pd_slack_upper_sec, [(i,c) for i in ids(pm, nw, :load) for (c, _) in enumerate(contingencies)], pd_slack_upper_sec)
end

function variable_qd_slack_sec(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    sec_loads = [(i, c) for i in ids(pm, nw, :load) for (c, _) in enumerate(contingencies)]
    qd_slack_lower_sec = var(pm, nw)[:qd_slack_lower_sec] = JuMP.@variable(pm.model,
        [(i,c) in sec_loads], base_name="$(nw)_qd_slack_lower_sec",
        start = 0,
        lower_bound = 0
    )

    qd_slack_upper_sec = var(pm, nw)[:qd_slack_upper_sec] = JuMP.@variable(pm.model,
        [(i, c) in sec_loads], base_name="$(nw)_qd_slack_upper_sec",
        start = 0,
        lower_bound = 0
    )

    for (i, load) in ref(pm, nw, :load)
        
        pd_val = load["pd"]
        max_load = abs(pd_val)

        for (c, _) in enumerate(contingencies)
            if pd_val > 0
                JuMP.set_upper_bound(qd_slack_lower_sec[(i,c)], max_load)
                JuMP.fix(qd_slack_upper_sec[(i,c)], 0.0; force=true)
            elseif pd_val < 0
                JuMP.set_upper_bound(qd_slack_upper_sec[(i,c)], max_load)
                JuMP.fix(qd_slack_lower_sec[(i,c)], 0.0; force=true)
            else
                JuMP.fix(qd_slack_lower_sec[(i,c)], 0.0; force=true)
                JuMP.fix(qd_slack_upper_sec[(i,c)], 0.0; force=true)
            end
        end
    end

    report && sol_component_value(pm, nw, :load, :qd_slack_lower_sec, [(i,c) for i in ids(pm, nw, :load) for (c, _) in enumerate(contingencies)], qd_slack_lower_sec)
    report && sol_component_value(pm, nw, :load, :qd_slack_upper_sec, [(i,c) for i in ids(pm, nw, :load) for (c, _) in enumerate(contingencies)], qd_slack_upper_sec)
end

"variable: `qq[j]` for `j` in `gen`"
function variable_gen_power_imaginary_slack(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    qmin_slack = var(pm, nw)[:qmin_slack] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :gen)], base_name="$(nw)_qmin_slack",
        start = 0
    )

    qmax_slack = var(pm, nw)[:qmax_slack] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :gen)], base_name="$(nw)_qmax_slack",
        start = 0
    )

    if bounded
        for (i, gen) in ref(pm, nw, :gen)
            JuMP.set_lower_bound(qmin_slack[i], 0)
            #JuMP.set_upper_bound(qmin_slack[i], gen["qmin"])

            JuMP.set_lower_bound(qmax_slack[i], 0)
            #JuMP.set_upper_bound(qmin_slack[i], gen["qmax"])
        end
    end

    report && sol_component_value(pm, nw, :gen, :qmin_slack, ids(pm, nw, :gen), qmin_slack)
    report && sol_component_value(pm, nw, :gen, :qmax_slack, ids(pm, nw, :gen), qmax_slack)
end

"generates variables for both `active` and `reactive` generation"
function variable_gen_power_slack_on_off(pm::AbstractPowerModel; kwargs...)
    variable_gen_power_real_slack_on_off(pm; kwargs...)
    variable_gen_power_imaginary_slack_on_off(pm; kwargs...)
end

"variable: `pg[j]` for `j` in `gen`"
function variable_gen_power_real_slack_on_off(pm::AbstractPowerModel; nw::Int=nw_id_default, report::Bool=true)
    pmin_slack_on_off = var(pm, nw)[:pmin_slack_on_off] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :gen)], base_name="$(nw)_pmin_slack_on_off",
        binary = true,
        start = 0
    )
    pmax_slack_on_off = var(pm, nw)[:pmax_slack_on_off] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :gen)], base_name="$(nw)_pmax_slack_on_off",
        binary = true,
        start = 0
    )

    report && sol_component_value(pm, nw, :gen, :pmin_slack_on_off, ids(pm, nw, :gen), pmin_slack_on_off)
    report && sol_component_value(pm, nw, :gen, :pmax_slack_on_off, ids(pm, nw, :gen), pmax_slack_on_off)
end

"generates variables for `active` and `reactive` generation"
function variable_gen_power_sec(pm::AbstractPowerModel, contingencies::Vector{Contingency}; kwargs...)
    slack_types = get(ismultinetwork(pm) ? pm.data["nw"]["1"] : pm.data, "slack_types", [])
    variable_gen_power_real_sec(     pm, bounded=:pg ∉ slack_types, contingencies; kwargs...)
    variable_gen_power_imaginary_sec(pm, bounded=:qg ∉ slack_types, contingencies; kwargs...)
end

"variable: `pg[j, c]` for `j` in `gen` for `c` in `contingencies`"
function variable_gen_power_real_sec(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    sec_gens = [(i,c) for i in ids(pm, nw, :gen) for (c, contingency) in enumerate(contingencies) if i ∉ get(contingency, :gen, [])]

    if ! PARAMS[:CONTINGENCIES][:PG_EQUALS_PGSEC]
        pg_sec = var(pm, nw)[:pg_sec] = JuMP.@variable(pm.model,
            [(i, c) in sec_gens], base_name="$(nw)_pg_sec",
            start = comp_start_value(ref(pm, nw, :gen, i), "pg_start")
        )
        if bounded
            for ((i, gen), c) in ref(pm, nw, :gen, contingencies)
                if i ∉ get(contingencies[c], :gen, [])
                    JuMP.set_lower_bound(pg_sec[(i, c)], gen["pmin"])
                    JuMP.set_upper_bound(pg_sec[(i, c)], gen["pmax"])
                    #@show(i, JuMP.lower_bound(pg_sec[(i,c)]), JuMP.upper_bound(pg_sec[(i,c)]))
                end
            end
        end
    else
        var(pm, nw)[:pg_sec] = Dict((i, c) => var(pm, nw, :pg, i) for (i, c) in sec_gens)
        pg_sec = var(pm, nw)[:pg_sec]
    end
    
    report && sol_component_value(pm, nw, :gen, :pg_sec, sec_gens, pg_sec)
end

"variable: `qq[j]` for `j` in `gen`"
function variable_gen_power_imaginary_sec(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    sec_gens = [(i,c) for i in ids(pm, nw, :gen) for (c, contingency) in enumerate(contingencies) if i ∉ get(contingency, :gen, [])]
    
    if ! PARAMS[:CONTINGENCIES][:PG_EQUALS_PGSEC]
        qg_sec = var(pm, nw)[:qg_sec] = JuMP.@variable(pm.model,
            [(i,c) in sec_gens], base_name="$(nw)_qg_sec",
            start = comp_start_value(ref(pm, nw, :gen, i), "qg_start")
        )

        if bounded
            for ((i, gen), c) in ref(pm, nw, :gen, contingencies)
                if i ∉ get(contingencies[c], :gen, [])
                    JuMP.set_lower_bound(qg_sec[(i,c)], gen["qmin"])
                    JuMP.set_upper_bound(qg_sec[(i,c)], gen["qmax"])
                end
            end
        end
    else
        var(pm, nw)[:qg_sec] = Dict((i, c) => var(pm, nw, :qg, i) for (i, c) in sec_gens)
        qg_sec = var(pm, nw)[:qg_sec]
    end

    report && sol_component_value(pm, nw, :gen, :qg_sec, sec_gens, qg_sec)
end

""
function variable_pd_slack_on_off(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    pd_slack_lower_on_off = var(pm, nw)[:pd_slack_lower_on_off] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :load)], base_name="$(nw)_pd_slack_lower_on_off",
        binary = true,
        start = 0
    )

    pd_slack_upper_on_off = var(pm, nw)[:pd_slack_upper_on_off] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :load)], base_name="$(nw)_pd_slack_upper_on_off",
        binary = true,
        start = 0
    )

    report && sol_component_value(pm, nw, :load, :pd_slack_upper_on_off, ids(pm, nw, :load), pd_slack_upper_on_off)
    report && sol_component_value(pm, nw, :load, :pd_slack_lower_on_off, ids(pm, nw, :load), pd_slack_lower_on_off)
end

"variable: ..."
function variable_branch_power_sec(pm::AbstractPowerModel, contingencies::Vector{Contingency}; kwargs...)
    variable_branch_power_real_sec(pm, contingencies; kwargs...)
    variable_branch_power_imaginary_sec(pm, contingencies; kwargs...)
end

"variable: `p[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_branch_power_real_sec(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, pf_ub_mult_sec=PARAMS[:LIMIT_MULT][:PF_UB_MULT_SEC])
    # one variable for each arc in each contingency where it is not disabled
    sec_arcs = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs) for (c, contingency) in enumerate(contingencies) if (l ∉ get(contingency, :branch, []) && get(ref(pm, nw, :branch, l), "type", "") ∉ ["breaker"])]

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
        for (c, contingency) in enumerate(contingencies)
            if l ∉ get(contingency, :branch, [])
                if haskey(branch, "pf_start")
                    f_idx = ((l,c), branch["f_bus"], branch["t_bus"])
                    JuMP.set_start_value(p_sec[f_idx], branch["pf_start"])
                end
                if haskey(branch, "pt_start")
                    t_idx = ((l,c), branch["t_bus"], branch["f_bus"])
                    JuMP.set_start_value(p_sec[t_idx], branch["pt_start"])
                end
            end
        end
    end

    # this explicit type erasure is necessary
    p_expr = Dict{Any,Any}( (((l,c),i,j), p_sec[((l,c),i,j)]) for ((l,c),i,j) in sec_arcs )
    merge!(p_expr, Dict( (((l,c),j,i), -1.0*p_sec[((l,c),i,j)]) for ((l,c),i,j) in sec_arcs))
    var(pm, nw)[:p_sec] = p_expr

    arcs_from_sec = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs_from) for (c, contingency) in enumerate(contingencies) if l ∉ get(contingency, :branch, [])]
    arcs_to_sec = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs_to) for (c, contingency) in enumerate(contingencies) if l ∉ get(contingency, :branch, [])]

    report && sol_component_value_edge(pm, nw, :branch, :pf, :pt, arcs_from_sec, arcs_to_sec, p_sec)
end

"variable: `q_sec[l,i,j]` for `(l,i,j)` in (non-contingencied) `arcs`"
function variable_branch_power_imaginary_sec(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    # one variable for each arc in each contingency where it is not disabled
    sec_arcs = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs) for (c, contingency) in enumerate(contingencies) if l ∉ get(contingency, :branch, [])]

    q_sec = var(pm, nw)[:q_sec] = JuMP.@variable(pm.model,
        [((l,c),i,j) in sec_arcs], base_name="$(nw)_q_sec",
        start = comp_start_value(ref(pm, nw, :branch, l), "q_start")
    )

    if bounded
        flow_lb, flow_ub = ref_calc_branch_flow_bounds(ref(pm, nw, :branch), ref(pm, nw, :bus))

        for sec_arc in sec_arcs
            (l,c),i,j = sec_arc
            if !isinf(flow_lb[l])
                JuMP.set_lower_bound(q_sec[sec_arc], flow_lb[l])
            end
            if !isinf(flow_ub[l])
                JuMP.set_upper_bound(q_sec[sec_arc], flow_ub[l])
            end
        end
    end

    for (l,branch) in ref(pm, nw, :branch)
        for (c, contingency) in enumerate(contingencies)
            if l ∉ get(contingency, :branch, [])
                if haskey(branch, "qf_start")
                    f_idx = ((l, c), branch["f_bus"], branch["t_bus"])
                    JuMP.set_start_value(q_sec[f_idx], branch["qf_start"])
                end
                if haskey(branch, "qt_start")
                    t_idx = ((l, c), branch["t_bus"], branch["f_bus"])
                    JuMP.set_start_value(q_sec[t_idx], branch["qt_start"])
                end
            end
        end
    end

    # this explicit type erasure is necessary
    q_expr = Dict{Any,Any}( (((l,c),i,j), q_sec[((l,c),i,j)]) for ((l,c),i,j) in sec_arcs )
    merge!(q_expr, Dict( (((l,c),j,i), -1.0*q_sec[((l,c),i,j)]) for ((l,c),i,j) in sec_arcs))
    var(pm, nw)[:q_sec] = q_expr

    arcs_from_sec = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs_from) for (c, contingency) in enumerate(contingencies) if l ∉ get(contingency, :branch, [])]
    arcs_to_sec = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs_to) for (c, contingency) in enumerate(contingencies) if l ∉ get(contingency, :branch, [])]

    report && sol_component_value_edge(pm, nw, :branch, :qf, :qt, arcs_from_sec, arcs_to_sec, q_expr)
end

"variable: ..."
function variable_branch_power_sec_complete(pm::AbstractPowerModel, contingencies::Vector{Contingency}; kwargs...)
    variable_branch_power_real_sec_complete(pm, contingencies; kwargs...)
    variable_branch_power_imaginary_sec_complete(pm, contingencies; kwargs...)
end

function variable_dcline_power_slack(pm::AbstractPowerModel; kwargs...)
    variable_dcline_power_real_slack(pm; kwargs...)
    #variable_dcline_power_imaginary(pm; kwargs...)
end

function variable_dcline_power_slack_on_off(pm::AbstractPowerModel; kwargs...)
    variable_dcline_power_real_slack_on_off(pm; kwargs...)
    #variable_dcline_power_imaginary_on_off(pm; kwargs...)
end

"variable: `p_dc[l,i,j]` for `(l,i,j)` in `arcs_dc`"
function variable_dcline_power_real_slack(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    p_dc_slack_lower = var(pm, nw)[:p_dc_slack_lower] = JuMP.@variable(pm.model,
        [arc in ref(pm, nw, :arcs_dc)], base_name="$(nw)_p_dc_slack_lower",
    )
    p_dc_slack_upper = var(pm, nw)[:p_dc_slack_upper] = JuMP.@variable(pm.model,
        [arc in ref(pm, nw, :arcs_dc)], base_name="$(nw)_p_dc_slack_upper",
    )

    if bounded
        for (l,dcline) in ref(pm, nw, :dcline)
            f_idx = (l, dcline["f_bus"], dcline["t_bus"])
            t_idx = (l, dcline["t_bus"], dcline["f_bus"])

            JuMP.set_lower_bound(p_dc_slack_lower[f_idx], 0)
            JuMP.set_upper_bound(p_dc_slack_lower[f_idx], dcline["pminf"])

            JuMP.set_lower_bound(p_dc_slack_lower[t_idx], 0)
            JuMP.set_upper_bound(p_dc_slack_lower[t_idx], dcline["pmint"])
            ##
            JuMP.set_lower_bound(p_dc_slack_upper[f_idx], 0)
            #JuMP.set_upper_bound(p_dc_slack_upper[f_idx], dcline["pmaxf"])

            JuMP.set_lower_bound(p_dc_slack_upper[t_idx], 0)
            #JuMP.set_upper_bound(p_dc_slack_upper[t_idx], dcline["pmaxt"])
        end
    end

    for (l,dcline) in ref(pm, nw, :dcline)
        if haskey(dcline, "pf")
            f_idx = (l, dcline["f_bus"], dcline["t_bus"])
            JuMP.set_start_value(p_dc_slack_lower[f_idx], 0)
            JuMP.set_start_value(p_dc_slack_upper[f_idx], 0)
        end

        if haskey(dcline, "pt")
            t_idx = (l, dcline["t_bus"], dcline["f_bus"])
            JuMP.set_start_value(p_dc_slack_lower[t_idx], 0)
            JuMP.set_start_value(p_dc_slack_upper[t_idx], 0)
        end
    end

    report && sol_component_value_edge(pm, nw, :dcline, :pf, :pt, ref(pm, nw, :arcs_from_dc), ref(pm, nw, :arcs_to_dc), p_dc_slack_lower)
    report && sol_component_value_edge(pm, nw, :dcline, :pf, :pt, ref(pm, nw, :arcs_from_dc), ref(pm, nw, :arcs_to_dc), p_dc_slack_upper)
end

"variable: `p_dc[l,i,j]` for `(l,i,j)` in `arcs_dc`"
function variable_dcline_power_real_slack_on_off(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    p_dc_slack_lower_on_off = var(pm, nw)[:p_dc_slack_lower_on_off] = JuMP.@variable(pm.model,
        [arc in ref(pm, nw, :arcs_dc)], base_name="$(nw)_p_dc_slack_lower_on_off",
        binary=true
    )
    p_dc_slack_upper_on_off = var(pm, nw)[:p_dc_slack_upper_on_off] = JuMP.@variable(pm.model,
        [arc in ref(pm, nw, :arcs_dc)], base_name="$(nw)_p_dc_slack_upper_on_off",
        binary=true
    )

    report && sol_component_value_edge(pm, nw, :dcline, :pf_b, :pt_b, ref(pm, nw, :arcs_from_dc), ref(pm, nw, :arcs_to_dc), p_dc_slack_lower_on_off)
    report && sol_component_value_edge(pm, nw, :dcline, :pf_b, :pt_b, ref(pm, nw, :arcs_from_dc), ref(pm, nw, :arcs_to_dc), p_dc_slack_upper_on_off)
end

function variable_dcline_power_sec(pm::AbstractPowerModel, contingencies::Vector{Contingency}; kwargs...)
    variable_dcline_power_real_sec(pm, contingencies; kwargs...)
    variable_dcline_power_imaginary_sec(pm, contingencies; kwargs...)
end

"variable: `p_dc[l,i,j]` for `(l,i,j)` in `arcs_dc`"
function variable_dcline_power_real_sec(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    # one variable for each arc in each contingency where it is not disabled
    sec_arcs = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs_dc) for (c, contingency) in enumerate(contingencies) if l ∉ get(contingency, :dcline, [])]
    
    p_dc_sec = var(pm, nw)[:p_dc_sec] = JuMP.@variable(pm.model,
        [sec_arc in sec_arcs], base_name="$(nw)_p_dc_sec",
    )

    for (c, contingency) in enumerate(contingencies)
        if bounded
            for (l,dcline) in ref(pm, nw, :dcline)
                if l ∉ contingency["dcline"]
                    f_idx = ((l,c), dcline["f_bus"], dcline["t_bus"])
                    t_idx = ((l,c), dcline["t_bus"], dcline["f_bus"])

                    JuMP.set_lower_bound(p_dc_sec[f_idx], dcline["pminf"])
                    JuMP.set_upper_bound(p_dc_sec[f_idx], dcline["pmaxf"])

                    JuMP.set_lower_bound(p_dc_sec[t_idx], dcline["pmint"])
                    JuMP.set_upper_bound(p_dc_sec[t_idx], dcline["pmaxt"])
                end
            end
        end

        for (l,dcline) in ref(pm, nw, :dcline)
            if haskey(dcline, "pf")
                f_idx = ((l,c), dcline["f_bus"], dcline["t_bus"])
                JuMP.set_start_value(p_dc_sec[f_idx], dcline["pf"])
            end

            if haskey(dcline, "pt")
                t_idx = ((l,c), dcline["t_bus"], dcline["f_bus"])
                JuMP.set_start_value(p_dc_sec[t_idx], dcline["pt"])
            end
        end
    end

    arcs_from_sec = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs_from_dc) for (c, contingency) in enumerate(contingencies) if l ∉ get(contingency, :dcline, [])]
    arcs_to_sec = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs_to_dc) for (c, contingency) in enumerate(contingencies) if l ∉ get(contingency, :dcline, [])]

    report && sol_component_value_edge(pm, nw, :dcline, :pf, :pt, arcs_from_sec, arcs_to_sec, p_dc_sec)
end

"variable: `q_dc[l,i,j]` for `(l,i,j)` in `arcs_dc`"
function variable_dcline_power_imaginary_sec(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    q_dc_sec = var(pm, nw)[:q_dc_sec] = JuMP.@variable(pm.model,
        [arc in ref(pm, nw, :arcs_dc)], base_name="$(nw)_q_dc",
    )

    if bounded
        for (l,dcline) in ref(pm, nw, :dcline)
            f_idx = (l, dcline["f_bus"], dcline["t_bus"])
            t_idx = (l, dcline["t_bus"], dcline["f_bus"])

            JuMP.set_lower_bound(q_dc_sec[f_idx], dcline["qminf"])
            JuMP.set_upper_bound(q_dc_sec[f_idx], dcline["qmaxf"])

            JuMP.set_lower_bound(q_dc_sec[t_idx], dcline["qmint"])
            JuMP.set_upper_bound(q_dc_sec[t_idx], dcline["qmaxt"])
        end
    end

    for (l,dcline) in ref(pm, nw, :dcline)
        if haskey(dcline, "qf")
            f_idx = (l, dcline["f_bus"], dcline["t_bus"])
            JuMP.set_start_value(q_dc_sec[f_idx], dcline["qf"])
        end

        if haskey(dcline, "qt")
            t_idx = (l, dcline["t_bus"], dcline["f_bus"])
            JuMP.set_start_value(q_dc_sec[t_idx], dcline["qt"])
        end
    end

    arcs_from_sec = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs_from_dc) for (c, contingency) in enumerate(contingencies) if l ∉ get(contingency, :dcline, [])]
    arcs_to_sec = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs_to_dc) for (c, contingency) in enumerate(contingencies) if l ∉ get(contingency, :dcline, [])]

    report && sol_component_value_edge(pm, nw, :dcline, :qf, :qt, arcs_from_sec, arcs_to_sec, q_dc_sec)
end

function variable_dcline_power_sec_complete(pm::AbstractPowerModel, contingencies::Vector{Contingency}; kwargs...)
    variable_dcline_power_real_sec_complete(pm, contingencies::Vector{Contingency}; kwargs...)
    #variable_dcline_power_imaginary_complete(pm; kwargs...)
end

"variable: `p_dc[l,i,j]` for `(l,i,j)` in `arcs_dc`"
function variable_dcline_power_real_sec_complete(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    # one variable for each arc in each contingency where it is not disabled
    sec_arcs = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs_dc) for (c, contingency) in enumerate(contingencies)]
    
    p_dc_sec = var(pm, nw)[:p_dc_sec] = JuMP.@variable(pm.model,
        [sec_arc in sec_arcs], base_name="$(nw)_p_dc_sec",
    )

    for (c, contingency) in enumerate(contingencies)
        if bounded
            for (l,dcline) in ref(pm, nw, :dcline)
                f_idx = ((l,c), dcline["f_bus"], dcline["t_bus"])
                t_idx = ((l,c), dcline["t_bus"], dcline["f_bus"])

                JuMP.set_lower_bound(p_dc_sec[f_idx], dcline["pminf"])
                JuMP.set_upper_bound(p_dc_sec[f_idx], dcline["pmaxf"])

                JuMP.set_lower_bound(p_dc_sec[t_idx], dcline["pmint"])
                JuMP.set_upper_bound(p_dc_sec[t_idx], dcline["pmaxt"])
            end
        end

        for (l,dcline) in ref(pm, nw, :dcline)
            if haskey(dcline, "pf")
                f_idx = ((l,c), dcline["f_bus"], dcline["t_bus"])
                JuMP.set_start_value(p_dc_sec[f_idx], dcline["pf"])
            end

            if haskey(dcline, "pt")
                t_idx = ((l,c), dcline["t_bus"], dcline["f_bus"])
                JuMP.set_start_value(p_dc_sec[t_idx], dcline["pt"])
            end
        end
    end

    arcs_from_sec = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs_from_dc) for (c, contingency) in enumerate(contingencies)]
    arcs_to_sec = [((l,c),i,j) for (l,i,j) in ref(pm, nw, :arcs_to_dc) for (c, contingency) in enumerate(contingencies)]

    report && sol_component_value_edge(pm, nw, :dcline, :pf, :pt, arcs_from_sec, arcs_to_sec, p_dc_sec)
end

"wrapper function for gen, branch and dcline power as well as voltage angle difference slack variables"
function variable_component_slack(pm::AbstractPowerModel, slack_types::Vector{Symbol}; kwargs...)
    if :pg       ∈ slack_types; variable_gen_power_real_slack(pm; kwargs...);        end
    if :qg       ∈ slack_types; variable_gen_power_imaginary_slack(pm; kwargs...);   end
    if :ratea    ∈ slack_types; variable_rate_a_slack(pm; kwargs...);                end
    if :dclinepf ∈ slack_types; variable_dcline_power_slack(pm; kwargs...);          end
    if :vad      ∈ slack_types; variable_vad_slack(pm; kwargs...);                   end
    if :vm       ∈ slack_types; variable_bus_voltage_magnitude_slack(pm; kwargs...); end
    if :pd       ∈ slack_types; variable_pd_slack(pm; kwargs...);                    end
    if :qd       ∈ slack_types; variable_qd_slack(pm; kwargs...);                    end
end

"wrapper function for gen, branch and dcline power as well as voltage angle difference slack binary variables"
function variable_component_slack_on_off(pm::AbstractPowerModel, slack_types::Vector{Symbol}; kwargs...)
    if :pg       ∈ slack_types; variable_gen_power_real_slack_on_off(pm; kwargs...);        end
    if :qg       ∈ slack_types; variable_gen_power_imaginary_slack_on_off(pm; kwargs...);   end
    if :ratea    ∈ slack_types; variable_rate_a_slack_on_off(pm; kwargs...);                end
    if :dclinepf ∈ slack_types; variable_dcline_power_slack_on_off(pm; kwargs...);          end
    if :vad      ∈ slack_types; variable_vad_slack_on_off(pm; kwargs...);                   end
    if :vm       ∈ slack_types; variable_bus_voltage_magnitude_slack_on_off(pm; kwargs...); end
    if :pd       ∈ slack_types; variable_pd_slack_on_off(pm; kwargs...);                    end
end

"variable: `0 <= z_branch[l] <= 1` for `l` in `branch`es"
function variable_branch_indicator(pm::AbstractPowerModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true, fix::Bool=true)
    sw_br_ids = [br_id for (br_id, branch) in ref(pm, nw, :branch) if is_branch_switchable(branch)]
    if !relax
        z_branch = var(pm, nw)[:z_branch] = JuMP.@variable(pm.model,
            [l in sw_br_ids], base_name="$(nw)_z_branch",
            binary = true,
            start = comp_start_value(ref(pm, nw, :branch, l), "z_branch_start", ref(pm, nw, :branch, l, "status"))
        )
    else
        z_branch = var(pm, nw)[:z_branch] = JuMP.@variable(pm.model,
            [l in sw_br_ids], base_name="$(nw)_z_branch",
            lower_bound = 0.0,
            upper_bound = 1.0,
            start = comp_start_value(ref(pm, nw, :branch, l), "z_branch_start", ref(pm, nw, :branch, l, "status"))
        )
    end

    if get(pm.data, "FIX_TS", false) && ! get(pm.data, "SUBPROBLEM", false)
        for l ∈ sw_br_ids
            JuMP.fix(z_branch[l], JuMP.start_value(z_branch[l]), force=true)
        end
    end

    if PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP] == 0
        for l ∈ sw_br_ids
            branch = ref(pm, nw, :branch, l)
            JuMP.fix(z_branch[l], get(branch, "status", branch["br_status"]), force=true)
        end
    end

    report && sol_component_value(pm, nw, :branch, :br_status, sw_br_ids, z_branch)
end

"variable: `0 <= z_branch_change[l] <= 1` for `l` in `branch`es"
function variable_branch_change_indicator(pm::AbstractPowerModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)
    sw_br_ids = [br_id for (br_id, branch) in ref(pm, nw, :branch) if is_branch_switchable(branch)]
    if !relax
        z_branch_change = var(pm, nw)[:z_branch_change] = JuMP.@variable(pm.model,
            [l in sw_br_ids], base_name="$(nw)_z_branch_change",
            binary = true,
            start = comp_start_value(ref(pm, nw, :branch, l), "z_branch_change_start", 1.0)
        )
    else
        z_branch_change = var(pm, nw)[:z_branch_change] = JuMP.@variable(pm.model,
            [l in sw_br_ids], base_name="$(nw)_z_branch_change",
            lower_bound = 0.0,
            upper_bound = 1.0,
            start = comp_start_value(ref(pm, nw, :branch, l), "z_branch_change_start", 1.0)
        )
    end

    if get(pm.data, "FIX_TS", false) && ! get(pm.data, "SUBPROBLEM", false)
        for l ∈ sw_br_ids
            JuMP.fix(z_branch_change[l], JuMP.start_value(z_branch_change[l]), force=true)
        end
    end

    if PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP]  == 0
        for l ∈ sw_br_ids
            JuMP.fix(z_branch_change[l], 0, force=true)
        end
    end

    # if branch status = 0, branch change means z_branch_change = 1 = z_branch and if status = 1, branch change means z_branch_change = 1 = 1 - z_branch
    # (only if single timestep)
    z_branch = var(pm, nw, :z_branch)
    for l in sw_br_ids
        branch = ref(pm, nw, :branch, l)
        if branch["status"] == 0
            @constraint(pm.model, z_branch_change[l] == z_branch[l])
        elseif branch["status"] == 1
            @constraint(pm.model, z_branch_change[l] == 1 - z_branch[l])
        else
            @error("branch $l has invalid status $(branch["status"]) for branch change indicator")
        end
    end

    report && sol_component_value(pm, nw, :branch, :br_change, sw_br_ids, z_branch_change)
end

"variable: `-1 <= tr_shift[l] <= 1` for `l` in `branch`es"
function variable_transformer_shift(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, discrete::Bool=PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_DISCRETE])
    shiftable_trafo_ids = [tr_id for (tr_id, tr) in ref(pm, nw, :branch) if is_branch_shiftable(pm, nw, tr_id)]
    tr_shift = var(pm, nw)[:tr_shift] = JuMP.@variable(pm.model,
        [l in shiftable_trafo_ids], base_name="$(nw)_tr_shift",
        start = comp_start_value(ref(pm, nw, :branch, l), "tr_shift_start", ref(pm, nw, :branch, l, "cur_shift_step")),
        lower_bound = ref(pm, nw, :branch, l)["min_step"],
        upper_bound = ref(pm, nw, :branch, l)["max_step"],
        integer = discrete
    )
    
    if ! PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_SHIFTING] 
        for l in shiftable_trafo_ids
            JuMP.fix(tr_shift[l], get(ref(pm, nw, :branch, i), "previous_shift_step", 0), force=true)
        end
    end

    if get(pm.data, "FIX_TR", false)
        for i in shiftable_trafo_ids
            JuMP.fix(tr_shift[i], JuMP.start_value(tr_shift[i]), force=true)
        end
    end

    report && sol_component_value(pm, nw, :branch, :tr_shift, shiftable_trafo_ids, tr_shift)
end

"variable: `0 <= tr_shift_change[l] <= 1` for `l` in `branch`es"
function variable_transformer_shift_change_indicator(pm::AbstractPowerModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)  
    if PARAMS[:OBJECTIVE][:COST_COEFFS][:TR_SW] == 0.0
        return
    end
    
    shiftable_trafo_ids = [tr_id for (tr_id, tr) in ref(pm, nw, :branch) if is_branch_shiftable(pm, nw, tr_id)]
    if !relax
        tr_shift_change = var(pm, nw)[:tr_shift_change] = JuMP.@variable(pm.model,
            [l in shiftable_trafo_ids], base_name="$(nw)_tr_shift_change",
            binary = true,
            start = comp_start_value(ref(pm, nw, :branch, l), "tr_shift_change_start", 0.0)
        )
    else
        tr_shift_change = var(pm, nw)[:tr_shift_change] = JuMP.@variable(pm.model,
            [l in shiftable_trafo_ids], base_name="$(nw)_tr_shift_change",
            lower_bound = 0.0,
            upper_bound = 1.0,
            start = comp_start_value(ref(pm, nw, :branch, l), "tr_shift_change_start", 0.0)
        )
    end

    if ! PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_SHIFTING] 
        for l in shiftable_trafo_ids
            JuMP.fix(tr_shift_change[l], 0, force=true)
        end
    end

     if get(pm.data, "FIX_TR", false)
        for l in shiftable_trafo_ids
            JuMP.fix(tr_shift_change[l], JuMP.start_value(tr_shift_change[l]), force=true)
        end
    end

    report && sol_component_value(pm, nw, :branch, :tr_shift_change, shiftable_trafo_ids, tr_shift_change)
end

"variables: `0 <= z_gen_startup[l] <= 1` and `0 <= z_gen_shutdown[l] <= 1` for `l` in `gen`"
function variable_gen_change_indicator(pm::AbstractPowerModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)
    switchable_gen_ids_one_per_gsk = [gen_id for (gen_id, gen) in ref(pm, nw, :gen) if is_gen_switchable(gen) && !haskey(gen, "gsk_name")]
    if haskey(ref(pm, nw), :gsk_factors)
        # all gens in a gsk share the same switching state
        for (_, gsk) in ref(pm, nw, :gsk_factors)
            push!(switchable_gen_ids_one_per_gsk, minimum(keys(gsk)))
        end
    end
    if !relax
        z_gen_startup = var(pm, nw)[:z_gen_startup] = JuMP.@variable(pm.model,
            [l in switchable_gen_ids_one_per_gsk], base_name="$(nw)_z_gen_startup",
            binary = true,
            start = comp_start_value(ref(pm, nw, :gen, l), "z_gen_startup_start", 0.0)
        )
        z_gen_shutdown = var(pm, nw)[:z_gen_shutdown] = JuMP.@variable(pm.model,
            [l in switchable_gen_ids_one_per_gsk], base_name="$(nw)_z_gen_shutdown",
            binary = true,
            start = comp_start_value(ref(pm, nw, :gen, l), "z_gen_shutdown_start", 0.0)
        )
    else
        z_gen_startup = var(pm, nw)[:z_gen_startup] = JuMP.@variable(pm.model,
            [l in switchable_gen_ids_one_per_gsk], base_name="$(nw)_z_gen_startup",
            lower_bound = 0.0,
            upper_bound = 1.0,
            start = comp_start_value(ref(pm, nw, :gen, l), "z_gen_startup_start", 0.0)
        )
        z_gen_shutdown = var(pm, nw)[:z_gen_shutdown] = JuMP.@variable(pm.model,
            [l in switchable_gen_ids_one_per_gsk], base_name="$(nw)_z_gen_shutdown",
            lower_bound = 0.0,
            upper_bound = 1.0,
            start = comp_start_value(ref(pm, nw, :gen, l), "z_gen_shutdown_start", 0.0)
        )
    end

    if haskey(ref(pm, nw), :gsk_factors)
        z_gen_startup_axis  = copy(z_gen_startup.axes[1])
        z_gen_startup_data  = copy(z_gen_startup.data)
        z_gen_shutdown_axis = copy(z_gen_shutdown.axes[1])
        z_gen_shutdown_data = copy(z_gen_shutdown.data)

        # all gens in a gsk share the same switching state
        for (_, gsk) in ref(pm, nw, :gsk_factors)
            gsk_i = minimum(keys(gsk))

            for (gen_id, _) in gsk
                if gen_id ≠ gsk_i
                    push!(z_gen_startup_axis, gen_id)
                    push!(z_gen_startup_data, z_gen_startup[gsk_i])

                    push!(z_gen_shutdown_axis, gen_id)
                    push!(z_gen_shutdown_data, z_gen_shutdown[gsk_i])
                end
            end
        end

        z_gen_startup  = var(pm, nw)[:z_gen_startup]  = JuMP.Containers.DenseAxisArray(z_gen_startup_data,  z_gen_startup_axis)
        z_gen_shutdown = var(pm, nw)[:z_gen_shutdown] = JuMP.Containers.DenseAxisArray(z_gen_shutdown_data, z_gen_shutdown_axis)
    end

    if get(pm.data, "FIX_UC", false) && ! get(pm.data, "SUBPROBLEM", false)
        for i in switchable_gen_ids_one_per_gsk
            JuMP.fix(z_gen_startup[i],  JuMP.start_value(z_gen_startup[i]),  force=true)
            JuMP.fix(z_gen_shutdown[i], JuMP.start_value(z_gen_shutdown[i]), force=true)

        end
    end

    report && sol_component_value(pm, nw, :gen, :gen_startup, switchable_gen_ids_one_per_gsk, z_gen_startup)
    report && sol_component_value(pm, nw, :gen, :gen_shutdown, switchable_gen_ids_one_per_gsk, z_gen_shutdown)
end

"generates variables for both `active` and `reactive` generation"
function variable_gen_power_change(pm::AbstractPowerModel; kwargs...)
    variable_gen_power_change_real(pm; kwargs...)
    #variable_gen_power_change_imaginary(pm; kwargs...)
end

"variable: `pg_delta[j]` for `j` in `gen`"
function variable_gen_power_change_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    pg_delta_plus = var(pm, nw)[:pg_delta_plus] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :gen)], base_name="$(nw)_pg_delta_plus",
         #lower_bound = 0.0,
         #upper_bound = ref(pm, nw, :gen, i)["RDP+"],
        start = comp_start_value(ref(pm, nw, :gen, i), "pg_delta_plus_start", 0.0)
    )

    
    pg_delta_minus = var(pm, nw)[:pg_delta_minus] = JuMP.@variable(pm.model,
    [i in ids(pm, nw, :gen)], base_name="$(nw)_pg_delta_minus",
         #lower_bound = 0.0,
         #upper_bound = ref(pm, nw, :gen, i)["RDP-"],
    start = comp_start_value(ref(pm, nw, :gen, i), "pg_delta_minus_start", 0.0)
    )

    for i in ids(pm, nw, :gen)
        JuMP.@constraint(pm.model, 0.0 <= pg_delta_plus[i], base_name="$(nw)_pg_delta_plus_lb_$(i)")
        JuMP.@constraint(pm.model, pg_delta_plus[i] <= ref(pm, nw, :gen, i)["RDP+"], base_name="$(nw)_pg_delta_plus_ub_$(i)")

        JuMP.@constraint(pm.model, 0.0 <= pg_delta_minus[i], base_name="$(nw)_pg_delta_minus_lb_$(i)")
        JuMP.@constraint(pm.model, pg_delta_minus[i] <= ref(pm, nw, :gen, i)["RDP-"], base_name="$(nw)_pg_delta_minus_ub_$(i)")
    end

    if get(pm.data, "FIX_UC", false) && ! get(pm.data, "SUBPROBLEM", false)
        for (i, gen) in ref(pm, nw, :gen)
            JuMP.fix(pg_delta_plus[i],  JuMP.start_value(pg_delta_plus[i]),  force=true)
            JuMP.fix(pg_delta_minus[i], JuMP.start_value(pg_delta_minus[i]), force=true)
        end
    end

    if PARAMS[:LIMITS][:REDISPATCH][:PG_DELTA_NUM_DAY]  == 0
        for (i, gen) in ref(pm, nw, :gen)
            JuMP.fix(pg_delta_plus[i], 0, force=true)
            JuMP.fix(pg_delta_minus[i], 0, force=true)
        end
    end

    report && sol_component_value(pm, nw, :gen, :pg_delta_plus, ids(pm, nw, :gen), pg_delta_plus)
    report && sol_component_value(pm, nw, :gen, :pg_delta_minus, ids(pm, nw, :gen), pg_delta_minus)
end

"variable: `qg_delta[j]` for `j` in `gen`"
function variable_gen_power_change_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    qg_delta = var(pm, nw)[:qg_delta] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :gen)], base_name="$(nw)_qg_delta",
        start = 0
    )

    if bounded
        for (i, gen) in ref(pm, nw, :gen)
            # TODO: handle qg negative case
            JuMP.set_lower_bound(qg_delta[i], gen["qg_0"] * (1-gen["qg_delta_percentage"]))
            JuMP.set_upper_bound(qg_delta[i], gen["qg_0"] * (1+gen["qg_delta_percentage"]))
        end
    end

    report && sol_component_value(pm, nw, :gen, :qg_delta, ids(pm, nw, :gen), qg_delta)
end

"variable: `pg_delta[j]` for `j` in `gen`"
function variable_gen_power_change_on_off(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    pg_delta_on_off = var(pm, nw)[:pg_delta_on_off] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :gen)], base_name="$(nw)_pg_delta_on_off",
        start = comp_start_value(ref(pm, nw, :gen, i), "pg_delta_on_off_start", 0),
        binary = true
    )

    if get(pm.data, "FIX_UC", false) && ! get(pm.data, "SUBPROBLEM", false)
        for (i, gen) in ref(pm, nw, :gen)
            JuMP.fix(pg_delta_on_off[i],  JuMP.start_value(pg_delta_on_off[i]), force=true)
        end
    end

    if PARAMS[:LIMITS][:REDISPATCH][:PG_DELTA_NUM_DAY]  == 0
        for (i, gen) in ref(pm, nw, :gen)
            JuMP.fix(pg_delta_on_off[i], 0, force=true)
        end
    end

    report && sol_component_value(pm, nw, :gen, :pg_delta_on_off, ids(pm, nw, :gen), pg_delta_on_off)
end

function variable_gen_power_on_off(pm::AbstractPowerModel; kwargs...)
    variable_gen_power_real_on_off(pm; kwargs...)
    variable_gen_power_imaginary_on_off(pm; kwargs...)
end

"generates variables for both `active` and `reactive` generation"
function variable_gen_power(pm::AbstractPowerModel; bounded::Bool=true, kwargs...)
    slack_types = get(ismultinetwork(pm) ? pm.data["nw"]["1"] : pm.data, "slack_types", [])
    variable_gen_power_real(     pm, bounded=:pg ∉ slack_types && bounded; kwargs...)
    variable_gen_power_imaginary(pm, bounded=:qg ∉ slack_types && bounded; kwargs...)
end

"variable: `pg[j]` for `j` in `gen`"
function variable_gen_power_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    pg = var(pm, nw)[:pg] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :gen)], base_name="$(nw)_pg",
        start = comp_start_value(ref(pm, nw, :gen, i), "pg_start")
    )

    # all gens are on, and no limit on pmin
    if PARAMS[:LIMITS][:GENERATOR][:GEN_REDISPATCH_OVERRIDE] 
        for (l, gen) in ref(pm, nw, :gen)
            gen["pmin"] = min(gen["pmin"], 0)
        end
    end

    if bounded
        for (i, gen) in ref(pm, nw, :gen)
            JuMP.set_upper_bound(pg[i], gen["pmax"])
            if PARAMS[:LIMITS][:GENERATOR][:GEN_SWITCHING] 
                JuMP.set_lower_bound(pg[i], min(0, gen["pmin"]))
            else
                status = haskey(gen, "status") ? gen["status"] : gen["gen_status"]
                JuMP.set_lower_bound(pg[i], gen["pmin"] * status)
            end
        end
    else
        for (i, gen) in ref(pm, nw, :gen)
            JuMP.set_upper_bound(pg[i], gen["pmax"] * PARAMS[:BIG_M_MULT][:PG]  )
            JuMP.set_lower_bound(pg[i], min(0, gen["pmin"]))
        end
    end

    if get(pm.data, "FIX_UC", false) && ! get(pm.data, "SUBPROBLEM", false)
        for (i, gen) in ref(pm, nw, :gen)
            JuMP.fix(pg[i], JuMP.start_value(pg[i]), force=true)
        end
    end

    if PARAMS[:LIMITS][:REDISPATCH][:PG_DELTA_NUM_DAY]  == 0
        for (i, gen) in ref(pm, nw, :gen)
            JuMP.fix(pg[i], get(gen, "pg_0", gen["pg"]), force=true)
        end
    end

    report && sol_component_value(pm, nw, :gen, :pg, ids(pm, nw, :gen), pg)
end

"variable: `0 <= z_gen[j] <= 1` for `j` in `gen`"
function variable_gen_indicator(pm::AbstractPowerModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)
    switchable_gen_ids_one_per_gsk = [gen_id for (gen_id, gen) in ref(pm, nw, :gen) if is_gen_switchable(gen) && !haskey(gen, "gsk_name")]
    if haskey(ref(pm, nw), :gsk_factors)
        # all gens in a gsk share the same switching state
        for (_, gsk) in ref(pm, nw, :gsk_factors)
            gsk_i = minimum(keys(gsk))
            if is_gen_switchable(ref(pm, nw, :gen, gsk_i))
                push!(switchable_gen_ids_one_per_gsk, gsk_i)
            end
        end
    end

    if !relax
        z_gen = var(pm, nw)[:z_gen] = JuMP.@variable(pm.model,
            [j in switchable_gen_ids_one_per_gsk], base_name="$(nw)_z_gen",
            binary = true,
            start = comp_start_value(ref(pm, nw, :gen, j), "z_gen_start", 1.0)
        )
    else
        z_gen = var(pm, nw)[:z_gen] = JuMP.@variable(pm.model,
            [j in switchable_gen_ids_one_per_gsk], base_name="$(nw)_z_gen",
            lower_bound = 0.0,
            upper_bound = 1.0,
            start = comp_start_value(ref(pm, nw, :gen, j), "z_gen_start", 1.0)
        )
    end

    if haskey(ref(pm, nw), :gsk_factors)
        z_gen_axis  = copy(z_gen.axes[1])
        z_gen_data  = copy(z_gen.data)

        # all gens in a gsk share the same switching state
        for (_, gsk) in ref(pm, nw, :gsk_factors)
            gsk_i = minimum(keys(gsk))
            if is_gen_switchable(ref(pm, nw, :gen, gsk_i))
                for (gen_id, _) in gsk
                    @assert is_gen_switchable(ref(pm, nw, :gen, gen_id))
                    if gen_id ≠ gsk_i
                        push!(z_gen_axis, gen_id)
                        push!(z_gen_data, z_gen[gsk_i])
                    end
                end
            end
        end

        z_gen = var(pm, nw)[:z_gen] = JuMP.Containers.DenseAxisArray(z_gen_data, z_gen_axis)
    end

    for l in switchable_gen_ids_one_per_gsk
        if get(pm.data, "FIX_UC", false) && ! get(pm.data, "SUBPROBLEM", false)
            JuMP.fix(z_gen[l], JuMP.start_value(z_gen[l]), force=true)
        end
    end
        
    report && sol_component_value(pm, nw, :gen, :gen_status, switchable_gen_ids_one_per_gsk, z_gen)
end