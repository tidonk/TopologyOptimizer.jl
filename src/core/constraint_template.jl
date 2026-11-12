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

#### "normal" constraints
"""
This constraint captures problem agnostic constraints that are used to link
the model's voltage variables together, in addition to the standard problem
formulation constraints.  The on/off name indicates that the voltages in this
constraint can be set to zero via an indicator variable

Notable examples include the constraints linking the voltages in the
ACTPowerModel, constraints linking convex relaxations of voltage variables.
"""

""
function constraint_theta_ref_sec(pm::AbstractPowerModel, i::Int, c::ContingencyID, contingency::Dict; nw::Int=nw_id_default)
    constraint_theta_ref_sec(pm, nw, i, c, contingency)
end

""
function constraint_power_balance_sec(pm::AbstractPowerModel, i::Int, c::ContingencyID, contingency::Dict; nw::Int=nw_id_default)
    bus = ref(pm, nw, :bus, i)
    bus_arcs = ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = ref(pm, nw, :bus_arcs_dc, i)
    bus_arcs_sw = ref(pm, nw, :bus_arcs_sw, i)
    bus_gens = ref(pm, nw, :bus_gens, i)
    bus_loads = ref(pm, nw, :bus_loads, i)
    bus_shunts = ref(pm, nw, :bus_shunts, i)
    bus_storage = ref(pm, nw, :bus_storage, i)

    bus_pd = Dict(k => ref(pm, nw, :load, k, "pd") for k in bus_loads)
    bus_qd = Dict(k => ref(pm, nw, :load, k, "qd") for k in bus_loads)

    bus_gs = Dict(k => ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_power_balance_sec(pm, nw, i, c, contingency, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
end

""
function constraint_power_balance_all(pm::AbstractPowerModel; nw::Int=nw_id_default)
    bus_arcs = Vector()
    bus_arcs_dc = Vector()
    bus_arcs_sw = Vector()
    bus_gens = reduce(vcat, [ref(pm, nw, :bus_gens, i) for i in ids(pm, nw, :bus)])
    bus_loads = reduce(vcat, [ref(pm, nw, :bus_loads, i) for i in ids(pm, nw, :bus)])
    bus_shunts = reduce(vcat, [ref(pm, nw, :bus_shunts, i) for i in ids(pm, nw, :bus)])
    bus_storage = reduce(vcat, [ref(pm, nw, :bus_storage, i) for i in ids(pm, nw, :bus)])

    bus_pd = Dict(k => ref(pm, nw, :load, k, "pd") for k in bus_loads)
    bus_qd = Dict(k => ref(pm, nw, :load, k, "qd") for k in bus_loads)

    bus_gs = Dict(k => ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_power_balance(pm, nw, 0, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
end

""
function constraint_ohms_yt_from_on_off(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = calc_branch_y(branch)
    tr, ti = calc_branch_t(branch)
    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]
    tm = branch["tap"]

    if PARAMS[:BIG_M_MULT][:SET_MAX_ANGLE_CLOSE] && is_branch_switchable(branch) && haskey(branch, "angmin_close") && haskey(branch, "angmax_close")
        vad_max = branch["angmax_close"]
        vad_min = branch["angmin_close"]
    else
        vad_max = ref(pm, nw, :off_angmax)
        vad_min = ref(pm, nw, :off_angmin)
    end

    constraint_ohms_yt_from_on_off(pm, nw, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max)
end

""
function constraint_ohms_yt_to_on_off(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = calc_branch_y(branch)
    tr, ti = calc_branch_t(branch)
    g_to = branch["g_to"]
    b_to = branch["b_to"]
    tm = branch["tap"]

    if PARAMS[:BIG_M_MULT][:SET_MAX_ANGLE_CLOSE] && is_branch_switchable(branch) && haskey(branch, "angmin_close") && haskey(branch, "angmax_close")
        vad_max = branch["angmax_close"]
        vad_min = branch["angmin_close"]
    else
        vad_max = ref(pm, nw, :off_angmax)
        vad_min = ref(pm, nw, :off_angmin)
    end

    constraint_ohms_yt_to_on_off(pm, nw, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max)
end

""
function constraint_ohms_yt_from_sec(pm::AbstractPowerModel, i::Int, c::ContingencyID, contingency::Dict; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = ((i,c), f_bus, t_bus)
    t_idx = ((i,c), t_bus, f_bus)

    g, b = calc_branch_y(branch)
    tr, ti = calc_branch_t(branch)
    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]
    tm = branch["tap"]

    constraint_ohms_yt_from_sec(pm, nw, i, c, contingency, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
end

""
function constraint_ohms_yt_from_sec_on_off(pm::AbstractPowerModel, i::Int, c::ContingencyID, contingency::Dict; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = ((i,c), f_bus, t_bus)
    t_idx = ((i,c), t_bus, f_bus)

    g, b = calc_branch_y(branch)
    tr, ti = calc_branch_t(branch)
    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]
    tm = branch["tap"]

    if PARAMS[:BIG_M_MULT][:SET_MAX_ANGLE_CLOSE] && is_branch_switchable(branch) && haskey(branch, "angmin_close") && haskey(branch, "angmax_close")
        vad_max = branch["angmax_close"]
        vad_min = branch["angmin_close"]
    else
        vad_max = ref(pm, nw, :off_angmax)
        vad_min = ref(pm, nw, :off_angmin)
    end

    constraint_ohms_yt_from_sec_on_off(pm, nw, i, c, contingency, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max)
end

""
function constraint_ohms_yt_to_sec(pm::AbstractPowerModel, i::Int, c::ContingencyID, contingency::Dict; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = ((i,c), f_bus, t_bus)
    t_idx = ((i,c), t_bus, f_bus)

    g, b = calc_branch_y(branch)
    tr, ti = calc_branch_t(branch)
    g_to = branch["g_to"]
    b_to = branch["b_to"]
    tm = branch["tap"]

    constraint_ohms_yt_to_sec(pm, nw, i, c, contingency, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
end

""
function constraint_ohms_yt_to_sec_on_off(pm::AbstractPowerModel, i::Int, c::ContingencyID, contingency::Dict; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = ((i,c), f_bus, t_bus)
    t_idx = ((i,c), t_bus, f_bus)

    g, b = calc_branch_y(branch)
    tr, ti = calc_branch_t(branch)
    g_to = branch["g_to"]
    b_to = branch["b_to"]
    tm = branch["tap"]

    if PARAMS[:BIG_M_MULT][:SET_MAX_ANGLE_CLOSE] && is_branch_switchable(branch) && haskey(branch, "angmin_close") && haskey(branch, "angmax_close")
        vad_max = branch["angmax_close"]
        vad_min = branch["angmin_close"]
    else
        vad_max = ref(pm, nw, :off_angmax)
        vad_min = ref(pm, nw, :off_angmin)
    end

    constraint_ohms_yt_to_sec_on_off(pm, nw, i, c, contingency, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max)
end

""
function constraint_power_magnitude_sqr_sec_on_off(pm::AbstractPowerModel, i::Int, c::Int, contingency::Contingency; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    arc_from = ((i,c), f_bus, t_bus)

    tm = branch["tap"]

    constraint_power_magnitude_sqr_sec_on_off(pm, nw, i, c, f_bus, arc_from, tm)
end

""
function constraint_power_magnitude_link_sec_on_off(pm::AbstractPowerModel, i::Int, c::Int, contingency::Contingency; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    arc_from = ((i,c), f_bus, t_bus)

    g, b = calc_branch_y(branch)
    tr, ti = calc_branch_t(branch)
    # CHECK: Do I need all these variables?
    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]
    g_to = branch["g_to"]
    b_to = branch["b_to"]
    tm = branch["tap"]

    constraint_power_magnitude_link_sec_on_off(pm, nw, i, c, arc_from, g, b, g_fr, b_fr, g_to, b_to, tr, ti, tm)
end

#### bound constraints
### thermal limit
## from

""
function constraint_thermal_limit_from_slack(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    if :ratea ∉ pm.data["slack_types"]
        constraint_thermal_limit_from(pm, i, nw=nw)
        return
    end
        
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    if !haskey(branch, "rate_a")
        Memento.error(_LOGGER, "constraint_thermal_limit_from_slack requires a rate_a value on all branches, calc_thermal_limits! can be used to generate reasonable values")
    end

    constraint_thermal_limit_from_slack(pm, nw, i, f_idx, get_adjusted_rate_a(pm, branch))
end

""
function constraint_thermal_limit_from_slack_on_off(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    if :ratea ∉ pm.data["slack_types"]
        constraint_thermal_limit_from_on_off(pm, i, nw=nw)
        return
    end
        
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    if !haskey(branch, "rate_a")
        Memento.error(_LOGGER, "constraint_thermal_limit_from_on_off requires a rate_a value on all branches, calc_thermal_limits! can be used to generate reasonable values")
    end

    constraint_thermal_limit_from_slack_on_off(pm, nw, i, f_idx, get_adjusted_rate_a(pm, branch))
end

""
function constraint_thermal_limit_from_sec(pm::AbstractPowerModel, i::Int, c::ContingencyID, contingency::Dict; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = ((i, c), f_bus, t_bus)

    if !haskey(branch, "rate_a")
        Memento.error(_LOGGER, "constraint_thermal_limit_from requires a rate_a value on all branches, calc_thermal_limits! can be used to generate reasonable values")
    end

    constraint_thermal_limit_from_sec(pm, nw, i, f_idx, get_adjusted_rate_a(pm, branch))
end

""
function constraint_thermal_limit_from_sec_on_off(pm::AbstractPowerModel, i::Int, c::ContingencyID, contingency::Dict; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = ((i, c), f_bus, t_bus)

    if !haskey(branch, "rate_a")
        Memento.error(_LOGGER, "constraint_thermal_limit_from_on_off requires a rate_a value on all branches, calc_thermal_limits! can be used to generate reasonable values")
    end

    constraint_thermal_limit_from_sec_on_off(pm, nw, i, f_idx, get_adjusted_rate_a(pm, branch))
end

""
function constraint_thermal_limit_from_slack_sec(pm::AbstractPowerModel, i::Int, c::ContingencyID, contingency::Dict; nw::Int=nw_id_default)
    if :ratea ∉ pm.data["slack_types"]
        constraint_thermal_limit_from_sec(pm, i, c, contingency, nw=nw)
        return
    end

    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = ((i, c), f_bus, t_bus)

    if !haskey(branch, "rate_a")
        Memento.error(_LOGGER, "constraint_thermal_limit_from_slack_sec requires a rate_a value on all branches, calc_thermal_limits! can be used to generate reasonable values")
    end

    constraint_thermal_limit_from_slack_sec(pm, nw, i, f_idx, get_adjusted_rate_a(pm, branch))
end

""
function constraint_thermal_limit_from_slack_sec_on_off(pm::AbstractPowerModel, i::Int, c::ContingencyID, contingency::Dict; nw::Int=nw_id_default)
    if :ratea ∉ pm.data["slack_types"]
        constraint_thermal_limit_from_sec_on_off(pm, i, c, contingency, nw=nw)
        return
    end

    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = ((i, c), f_bus, t_bus)

    if !haskey(branch, "rate_a")
        Memento.error(_LOGGER, "constraint_thermal_limit_from_on_off requires a rate_a value on all branches, calc_thermal_limits! can be used to generate reasonable values")
    end

    constraint_thermal_limit_from_slack_sec_on_off(pm, nw, i, f_idx, get_adjusted_rate_a(pm, branch))
end

## to

""
function constraint_thermal_limit_to_slack(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    if :ratea ∉ pm.data["slack_types"]
        constraint_thermal_limit_to(pm, i, nw=nw)
        return
    end

    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    t_idx = (i, t_bus, f_bus)

    if !haskey(branch, "rate_a")
        Memento.error(_LOGGER, "constraint_thermal_limit_to_slack requires a rate_a value on all branches, calc_thermal_limits! can be used to generate reasonable values")
    end

    constraint_thermal_limit_to_slack(pm, nw, i, t_idx, get_adjusted_rate_a(pm, branch))
end

""
function constraint_thermal_limit_to_slack_on_off(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    if :ratea ∉ pm.data["slack_types"]
        constraint_thermal_limit_to_on_off(pm, i, nw=nw)
        return
    end

    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    t_idx = (i, t_bus, f_bus)

    if !haskey(branch, "rate_a")
        Memento.error(_LOGGER, "constraint_thermal_limit_to_on_off requires a rate_a value on all branches, calc_thermal_limits! can be used to generate reasonable values")
    end

    constraint_thermal_limit_to_slack_on_off(pm, nw, i, t_idx, get_adjusted_rate_a(pm, branch))
end

""
function constraint_thermal_limit_to_sec(pm::AbstractPowerModel, i::Int, c::ContingencyID, contingency::Dict; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    t_idx = ((i, c), t_bus, f_bus)

    if !haskey(branch, "rate_a")
        Memento.error(_LOGGER, "constraint_thermal_limit_to requires a rate_a value on all branches, calc_thermal_limits! can be used to generate reasonable values")
    end

    constraint_thermal_limit_to_sec(pm, nw, i, t_idx, get_adjusted_rate_a(pm, branch))
end

""
function constraint_thermal_limit_to_sec_on_off(pm::AbstractPowerModel, i::Int, c::ContingencyID, contingency::Dict; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    t_idx = ((i, c), t_bus, f_bus)

    if !haskey(branch, "rate_a")
        Memento.error(_LOGGER, "constraint_thermal_limit_to_on_off requires a rate_a value on all branches, calc_thermal_limits! can be used to generate reasonable values")
    end

    constraint_thermal_limit_to_sec_on_off(pm, nw, i, t_idx, get_adjusted_rate_a(pm, branch))
end

""
function constraint_thermal_limit_to_slack_sec(pm::AbstractPowerModel, i::Int, c::ContingencyID, contingency::Dict; nw::Int=nw_id_default)
    if :ratea ∉ pm.data["slack_types"]
        constraint_thermal_limit_to_sec(pm, i, c, contingency, nw=nw)
        return
    end

    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    t_idx = ((i, c), t_bus, f_bus)

    if !haskey(branch, "rate_a")
        Memento.error(_LOGGER, "constraint_thermal_limit_to_slack_sec requires a rate_a value on all branches, calc_thermal_limits! can be used to generate reasonable values")
    end

    constraint_thermal_limit_to_slack_sec(pm, nw, i, t_idx, get_adjusted_rate_a(pm, branch))
end

""
function constraint_thermal_limit_to_slack_sec_on_off(pm::AbstractPowerModel, i::Int, c::ContingencyID, contingency::Dict; nw::Int=nw_id_default)
    if :ratea ∉ pm.data["slack_types"]
        constraint_thermal_limit_to_sec_on_off(pm, i, c, contingency, nw=nw)
        return
    end

    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    t_idx = ((i, c), t_bus, f_bus)

    if !haskey(branch, "rate_a")
        Memento.error(_LOGGER, "constraint_thermal_limit_to_on_off requires a rate_a value on all branches, calc_thermal_limits! can be used to generate reasonable values")
    end

    constraint_thermal_limit_to_slack_sec_on_off(pm, nw, i, t_idx, get_adjusted_rate_a(pm, branch))
end

### vad

""
function constraint_voltage_angle_difference(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    pair = (f_bus, t_bus)
    buspair = ref(pm, nw, :buspairs, pair)

    #if buspair["branch"] == i
        constraint_voltage_angle_difference(pm, nw, f_idx, buspair["angmin"], buspair["angmax"])
    #end
end

""
function constraint_voltage_angle_difference_on_off(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_idx = (i, branch["f_bus"], branch["t_bus"])

    if PARAMS[:BIG_M_MULT][:SET_MAX_ANGLE_CLOSE] && is_branch_switchable(branch) && haskey(branch, "angmin_close") && haskey(branch, "angmax_close")
        vad_max = branch["angmax_close"]
        vad_min = branch["angmin_close"]
    else
        vad_max = ref(pm, nw, :off_angmax)
        vad_min = ref(pm, nw, :off_angmin)
    end

    constraint_voltage_angle_difference_on_off(pm, nw, f_idx, branch["angmin"], branch["angmax"], vad_min, vad_max)
end

""
function constraint_voltage_angle_difference_slack(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    if :vad ∉ pm.data["slack_types"]
        constraint_voltage_angle_difference(pm, i, nw=nw)
        return
    end

    branch = ref(pm, nw, :branch, i)
    f_idx = (i, branch["f_bus"], branch["t_bus"])

    if PARAMS[:BIG_M_MULT][:SET_MAX_ANGLE_CLOSE] && is_branch_switchable(branch) && haskey(branch, "angmin_close") && haskey(branch, "angmax_close")
        vad_max = branch["angmax_close"]
        vad_min = branch["angmin_close"]
    else
        vad_max = ref(pm, nw, :off_angmax)
        vad_min = ref(pm, nw, :off_angmin)
    end

    constraint_voltage_angle_difference_slack(pm, nw, f_idx, branch["angmin"], branch["angmax"], vad_min, vad_max)
end

""
function constraint_voltage_angle_difference_slack_on_off(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    if :vad ∉ pm.data["slack_types"]
        constraint_voltage_angle_difference_on_off(pm, i, nw=nw)
        return
    end

    branch = ref(pm, nw, :branch, i)
    f_idx = (i, branch["f_bus"], branch["t_bus"])

    if PARAMS[:BIG_M_MULT][:SET_MAX_ANGLE_CLOSE] && is_branch_switchable(branch) && haskey(branch, "angmin_close") && haskey(branch, "angmax_close")
        vad_max = branch["angmax_close"]
        vad_min = branch["angmin_close"]
    else
        vad_max = ref(pm, nw, :off_angmax)
        vad_min = ref(pm, nw, :off_angmin)
    end

    constraint_voltage_angle_difference_slack_on_off(pm, nw, f_idx, branch["angmin"], branch["angmax"], vad_min, vad_max)
end

""
function constraint_voltage_angle_difference_sec(pm::AbstractPowerModel, i::Int, c::ContingencyID, contingency::Dict; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = ((i,c), f_bus, t_bus)
    pair = (f_bus, t_bus)
    buspair = ref(pm, nw, :buspairs, pair)

    #if buspair["branch"] == i
        constraint_voltage_angle_difference_sec(pm, nw, c, contingency, f_idx, buspair["angmin"], buspair["angmax"])
    #end
end

""
function constraint_voltage_angle_difference_sec_on_off(pm::AbstractPowerModel, i::Int, c::ContingencyID, contingency::Dict; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_idx = ((i,c), branch["f_bus"], branch["t_bus"])

    if PARAMS[:BIG_M_MULT][:SET_MAX_ANGLE_CLOSE] && is_branch_switchable(branch) && haskey(branch, "angmin_close") && haskey(branch, "angmax_close")
        vad_max = branch["angmax_close"]
        vad_min = branch["angmin_close"]
    else
        vad_max = ref(pm, nw, :off_angmax)
        vad_min = ref(pm, nw, :off_angmin)
    end

    constraint_voltage_angle_difference_sec_on_off(pm, nw, c, contingency, f_idx, branch["angmin"], branch["angmax"], vad_min, vad_max)
end

""
function constraint_voltage_angle_difference_slack_sec(pm::AbstractPowerModel, i::Int, c::ContingencyID, contingency::Dict; nw::Int=nw_id_default)
    if :vad ∉ pm.data["slack_types"]
        constraint_voltage_angle_difference_sec(pm, i, c, contingency, nw=nw)
        return
    end

    branch = ref(pm, nw, :branch, i)
    f_idx = ((i,c), branch["f_bus"], branch["t_bus"])

    if PARAMS[:BIG_M_MULT][:SET_MAX_ANGLE_CLOSE] && is_branch_switchable(branch) && haskey(branch, "angmin_close") && haskey(branch, "angmax_close")
        vad_max = branch["angmax_close"]
        vad_min = branch["angmin_close"]
    else
        vad_max = ref(pm, nw, :off_angmax)
        vad_min = ref(pm, nw, :off_angmin)
    end

    constraint_voltage_angle_difference_slack_sec(pm, nw, c, contingency, f_idx, branch["angmin"], branch["angmax"], vad_min, vad_max)
end

""
function constraint_voltage_angle_difference_slack_sec_on_off(pm::AbstractPowerModel, i::Int, c::ContingencyID, contingency::Dict; nw::Int=nw_id_default)
    if :vad ∉ pm.data["slack_types"]
        constraint_voltage_angle_difference_sec_on_off(pm, i, c, contingency, nw=nw)
        return
    end
    
    branch = ref(pm, nw, :branch, i)
    
    f_idx = ((i,c), branch["f_bus"], branch["t_bus"])
    
    if PARAMS[:BIG_M_MULT][:SET_MAX_ANGLE_CLOSE] && is_branch_switchable(branch) && haskey(branch, "angmin_close") && haskey(branch, "angmax_close")
        vad_max = branch["angmax_close"]
        vad_min = branch["angmin_close"]
    else
        vad_max = ref(pm, nw, :off_angmax)
        vad_min = ref(pm, nw, :off_angmin)
    end
    
    if is_branch_switchable(branch)
        # for breakers and couplers: angmin/angmax = 0, i.e., the angle will simply be conserved
        constraint_voltage_angle_difference_slack_sec_on_off(pm, nw, c, contingency, f_idx, branch["angmin"], branch["angmax"], vad_min, vad_max)
    else
        constraint_voltage_angle_difference_slack_sec(pm, nw, c, contingency, f_idx, branch["angmin"], branch["angmax"], vad_min, vad_max)
    end
end

### dcline

""
function constraint_dcline_power_losses_sec(pm::AbstractPowerModel, i::Int, c::ContingencyID, contingency::Dict; nw::Int=nw_id_default)
    dcline = ref(pm, nw, :dcline, i)
    f_bus = dcline["f_bus"]
    t_bus = dcline["t_bus"]
    f_idx = ((i, c), f_bus, t_bus)
    t_idx = ((i, c), t_bus, f_bus)
    loss0 = dcline["loss0"]
    loss1 = dcline["loss1"]

    constraint_dcline_power_losses_sec(pm, nw, f_bus, t_bus, f_idx, t_idx, loss0, loss1)
end

"limit a gen's power to 0 if it is off, and to pmax otherwise"
function constraint_gen_power_on_off(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    gen = ref(pm, nw, :gen, i)

    if !is_gen_switchable(gen)
        return
    end

    if :pg ∉ ref(pm, nw, :slack_types)
        constraint_gen_power_on_off(pm, nw, i, gen["pmin"], gen["pmax"], gen["qmin"], gen["qmax"])
    else
        constraint_gen_power_slack_on_off(pm, nw, i, gen["pmin"], gen["pmax"], gen["qmin"], gen["qmax"])
    end

    #constraint_gen_power_real_on_off(pm, nw, i, gen["pmin"], gen["pmax"])
    #constraint_gen_power_imaginary_on_off(pm, nw, i, gen["qmin"], gen["qmax"])
end

### wrapper
"helper function to get component lower and upper bounds for slacks"
function _get_slack_bounds(pm::AbstractPowerModel, nw::Int, i::Int, component::Dict, slack_type::Symbol)
    lb, ub, lb2, ub2 = 0, 0, 0, 0
    if slack_type == :branch
        @warn("use :ratea instead of :branch slack type")
        lbs, ubs = ref_calc_branch_flow_bounds(ref(pm, nw, :branch), ref(pm, nw, :bus))
        lb, ub = (lbs[i], ubs[i])
    elseif slack_type == :pg
        lb, ub = (component["pmin"], component["pmax"])
    elseif slack_type == :qg
        lb, ub = (component["qmin"], component["qmax"])
    elseif slack_type == :dclinepf
        lb,  ub  = (component["pminf"], component["pmaxf"])
        lb2, ub2 = (component["pmint"], component["pmaxt"])
    elseif slack_type == :vm
        lb = component["vmin"]
        ub = component["vmax"]
    elseif slack_type == :ratea
        ub = get_adjusted_rate_a(pm, component)
    elseif slack_type == :vad
        lb, ub = (component["angmin"], component["angmax"])
    elseif slack_type == :pd
        lb, ub = (-abs(component["pd"]), abs(component["pd"]))
    elseif slack_type == :qd
        lb, ub = (-abs(component["qd"]), abs(component["qd"]))
    else
        @warn("bounds for slack type unknown, using 0")
    end

    return lb, ub, lb2, ub2
end

"wrapper function for linking slack variables to their corresponding original variables"
function constraint_component_slack(pm::AbstractPowerModel, i::Int, slack_type::Symbol; nw::Int=nw_id_default)
    slack_component_type = CONSTS[:INTERNAL_MAPPINGS][:SLACK_VARIABLES][:SLACK_COMPONENT_TYPE][slack_type]
    component = ref(pm, nw, slack_component_type)[i]

    lb, ub, lb2, ub2 = _get_slack_bounds(pm, nw, i, component, slack_type)

    constraint_component_slack(pm, nw, i, slack_type, component, lb, ub, lb2, ub2)
end

"wrapper function for linking slack variables to their corresponding original variables for security constraints"
function constraint_component_slack_sec(pm::AbstractPowerModel, i::Int, c::ContingencyID, contingency::Contingency, slack_type::Symbol; nw::Int=nw_id_default)
    slack_component_type = CONSTS[:INTERNAL_MAPPINGS][:SLACK_VARIABLES][:SLACK_COMPONENT_TYPE][slack_type]
    component = ref(pm, nw, slack_component_type)[i]

    lb, ub, lb2, ub2 = _get_slack_bounds(pm, nw, i, component, slack_type)

    constraint_component_slack_sec(pm, nw, i, c, contingency, slack_type, component, lb, ub, lb2, ub2)
end

"wrapper function for linking binary slack variables to their corresponding normal slack variables"
function constraint_component_slack_on_off(pm::AbstractPowerModel, i::Int, slack_type::Symbol; nw::Int=nw_id_default)
    slack_component_type = CONSTS[:INTERNAL_MAPPINGS][:SLACK_VARIABLES][:SLACK_COMPONENT_TYPE][slack_type]    
    component = ref(pm, nw, slack_component_type)[i]

    # get bounds based on the data
    lb, ub, lb2, ub2 = _get_slack_bounds(pm, nw, i, component, slack_type)
    
    # fix index, since branches are adressed by their (f_bus, t_bus), not idx
    if slack_component_type == :branch
        f_bus = component["f_bus"]
        t_bus = component["t_bus"]
        f_idx = (i, f_bus, t_bus)
        t_idx = (i, t_bus, f_bus)
        constraint_component_slack_on_off(pm, nw, i, slack_type, f_idx, lb, ub, lb2, ub2)
        constraint_component_slack_on_off(pm, nw, i, slack_type, t_idx, lb, ub, lb2, ub2)
    else
        f_idx = i
        constraint_component_slack_on_off(pm, nw, i, slack_type, f_idx, lb, ub, lb2, ub2)
    end
end

""
function constraint_gen_power_delta_linking(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    gen = ref(pm, nw, :gen, i)
 
    constraint_gen_power_delta_linking(pm, nw, i, gen["pg_0"])
    #constraint_gen_power_imaginary_delta_linking(pm, nw, i, gen["qg_0"])
end

"link gens that should produce the same power with regard to a specific GSK factor"
function constraint_gen_power_gsk_linking(pm::AbstractPowerModel; nw::Int=nw_id_default)
    if haskey(ref(pm, nw), :gsk_factors)

        for (gsk_name, gsk) in ref(pm, nw, :gsk_factors)
            constraint_gen_power_gsk_linking_real(pm, nw, string(gsk_name), gsk)
            #constraint_gen_power_gsk_linking_imaginary(pm, nw, i, gen["qg_0"])
        end

    end
end

""
function constraint_theta_ref(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    constraint_theta_ref(pm, nw, i)
end

""
function constraint_power_balance(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    bus = ref(pm, nw, :bus, i)
    bus_arcs = ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = ref(pm, nw, :bus_arcs_dc, i)
    bus_arcs_sw = ref(pm, nw, :bus_arcs_sw, i)
    bus_gens = ref(pm, nw, :bus_gens, i)
    bus_loads = ref(pm, nw, :bus_loads, i)
    bus_shunts = ref(pm, nw, :bus_shunts, i)
    bus_storage = ref(pm, nw, :bus_storage, i)

    bus_pd = Dict(k => ref(pm, nw, :load, k, "pd") for k in bus_loads)
    bus_qd = Dict(k => ref(pm, nw, :load, k, "qd") for k in bus_loads)

    bus_gs = Dict(k => ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_power_balance(pm, nw, i, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
end

""
function constraint_ohms_yt_from(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = calc_branch_y(branch)
    tr, ti = calc_branch_t(branch)
    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]
    tm = branch["tap"]

    constraint_ohms_yt_from(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
end


""
function constraint_ohms_yt_to(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = calc_branch_y(branch)
    tr, ti = calc_branch_t(branch)
    g_to = branch["g_to"]
    b_to = branch["b_to"]
    tm = branch["tap"]

    constraint_ohms_yt_to(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
end