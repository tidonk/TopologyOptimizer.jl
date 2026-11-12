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

"apo models ignore reactive power flows"
function variable_branch_power_imaginary_slack(pm::AbstractActivePowerModel; nw::Int=nw_id_default, report::Bool=true, kwargs...)
    report && sol_component_fixed(pm, nw, :branch, :qf_slack, ids(pm, nw, :branch), NaN)
    report && sol_component_fixed(pm, nw, :branch, :qt_slack, ids(pm, nw, :branch), NaN)
end

"apo models ignore reactive power flows"
function variable_branch_power_imaginary_slack_on_off(pm::AbstractActivePowerModel; nw::Int=nw_id_default, report::Bool=true, kwargs...)
    report && sol_component_fixed(pm, nw, :branch, :qf_slack_on_off, ids(pm, nw, :branch), NaN)
    report && sol_component_fixed(pm, nw, :branch, :qt_slack_on_off, ids(pm, nw, :branch), NaN)
end

"apo models ignore reactive power flows"
function variable_branch_power_imaginary_sec(pm::AbstractActivePowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, report::Bool=true, kwargs...)
    report && sol_component_fixed(pm, nw, :branch, :qf_sec, ids(pm, nw, :branch), NaN)
    report && sol_component_fixed(pm, nw, :branch, :qt_sec, ids(pm, nw, :branch), NaN)
end

"apo models ignore reactive power flows"
function variable_branch_power_imaginary_sec_complete(pm::AbstractActivePowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, report::Bool=true, kwargs...)
    report && sol_component_fixed(pm, nw, :branch, :qf_sec, ids(pm, nw, :branch), NaN)
    report && sol_component_fixed(pm, nw, :branch, :qt_sec, ids(pm, nw, :branch), NaN)
end

"ensure that all slack variables for a component for a parameter are equal in all timesteps"
function constraint_slack_linking(pm, n_timesteps::Int)
    slack_types = pm.data["slack_types"]
    for (slack_type, component_slack_variables) in CONSTS[:INTERNAL_MAPPINGS][:SLACK_VARIABLES][:SLACK_VARIABLE_INTERNAL_NAME_MAP]
        if slack_type in slack_types
            slack_component_type = CONSTS[:INTERNAL_MAPPINGS][:SLACK_VARIABLES][:SLACK_COMPONENT_TYPE][slack_type]
            for component_id in ids(pm, 1, slack_component_type)
                for slack_variable in component_slack_variables
                    if slack_type == :qg && !(typeof(pm) <: AbstractACPModel)
                        continue
                    end
                    
                    if slack_component_type == :branch
                        branch = ref(pm, 1, :branch, component_id)
                        f_bus = branch["f_bus"]
                        t_bus = branch["t_bus"]
                        f_idx = (component_id, f_bus, t_bus)
                    else 
                        f_idx = component_id
                    end
                    for n in 2:n_timesteps
                        slack_var_1 = var(pm, n-1, slack_variable, f_idx)
                        slack_var_2 = var(pm, n, slack_variable, f_idx)
                        JuMP.@constraint(pm.model, slack_var_1 == slack_var_2)
                    end
                end
            end
        end
    end
end

""
function constraint_pg_equal_pgsec(pm::AbstractActivePowerModel; nw::Int=nw_id_default)
    pg      = get(var(pm, nw),   :pg, Dict())
    pg_sec  = get(var(pm, nw),   :pg_sec, Dict())

    for i in ids(pm, nw, :gen)
        for (c, contingency) in enumerate(ref(pm, nw, :contingencies))
            if i ∉ get(contingency, :gen, [])
                JuMP.@constraint(pm.model, pg[i] == pg_sec[(i, c)], base_name = "$(nw)_pg_eq_pgsec_$(i)")
            end
        end
    end
end

""
function constraint_thermal_limit_from_slack(pm::AbstractActivePowerModel, n::Int, i, f_idx, rate_a)
    @assert(:ratea ∈ ref(pm, n, :slack_types))

    p_fr = var(pm, n, :p, f_idx)

    if is_branch_monitored(pm, n, i)
        rate_a_slack = var(pm, n, :rate_a_slack, f_idx)
        JuMP.@constraint(pm.model, p_fr <=  rate_a + rate_a_slack, base_name = "$(n)_rating_from_slack_ub_$(i)")
        JuMP.@constraint(pm.model, p_fr >= -rate_a - rate_a_slack, base_name = "$(n)_rating_from_slack_lb_$(i)")
    end
end

""
function constraint_thermal_limit_from_slack_on_off(pm::AbstractActivePowerModel, n::Int, i, f_idx, rate_a)
    @assert(:ratea ∈ ref(pm, n, :slack_types))

    p_fr = var(pm, n, :p, f_idx)

    p_max = ref(pm, n, :p_max)

    if is_branch_monitored(pm, n, i)
        # bound
        rate_a_slack = var(pm, n, :rate_a_slack, f_idx)
        JuMP.@constraint(pm.model, p_fr <=  rate_a + rate_a_slack, base_name = "$(n)_rating_from_slack_ub_$(i)")
        JuMP.@constraint(pm.model, p_fr >= -rate_a - rate_a_slack, base_name = "$(n)_rating_from_slack_lb_$(i)")
    end
    if is_branch_switchable(pm, n, i)
        # switching
        # TODO: make this an indicator constraint
        z = var(pm, n, :z_branch, i)
        JuMP.@constraint(pm.model, p_fr <=   p_max * z, base_name = "$(n)_rating_from_slack_on_off_ub_$(i)")
        JuMP.@constraint(pm.model, p_fr >= - p_max * z, base_name = "$(n)_rating_from_slack_on_off_lb_$(i)")
    end
end

""
function constraint_thermal_limit_from_sec(pm::AbstractActivePowerModel, n::Int, i::Int, f_idx, rate_a)
    p_fr_sec = var(pm, n, :p_sec, f_idx)

    (l,c),i,j = f_idx

    JuMP.@constraint(pm.model, p_fr_sec <=  rate_a, base_name = "$(n)_C$(c)_rating_from_ub_$(l)")
    JuMP.@constraint(pm.model, p_fr_sec >= -rate_a, base_name = "$(n)_C$(c)_rating_from_lb_$(l)")
end

""
function constraint_thermal_limit_from_sec_on_off(pm::AbstractActivePowerModel, n::Int, i, f_idx, rate_a)
    p_fr_sec = var(pm, n, :p_sec, f_idx)
    z = var(pm, n, :z_branch, i)

    (l,c),i,j = f_idx

    JuMP.@constraint(pm.model, p_fr_sec <=  rate_a * z, base_name = "$(n)_C$(c)_rating_from_on_off_ub_$(l)")
    JuMP.@constraint(pm.model, p_fr_sec >= -rate_a * z, base_name = "$(n)_C$(c)_rating_from_on_off_lb_$(l)")
end

""
function constraint_thermal_limit_from_slack_sec(pm::AbstractActivePowerModel, n::Int, i, f_idx, rate_a)
    @assert(:ratea ∈ ref(pm, n, :slack_types))
    
    p_fr_sec = var(pm, n, :p_sec, f_idx)

    (l,c),i,j = f_idx
    
    if is_branch_monitored(pm, n, l)
        rate_a_slack = var(pm, n, :rate_a_slack_sec, f_idx)
        JuMP.@constraint(pm.model, p_fr_sec <=  rate_a + rate_a_slack, base_name = "$(n)_C$(c)_rating_from_slack_ub_$(l)")
        JuMP.@constraint(pm.model, p_fr_sec >= -rate_a - rate_a_slack, base_name = "$(n)_C$(c)_rating_from_slack_lb_$(l)")
    end
end

""
function constraint_thermal_limit_from_slack_sec_on_off(pm::AbstractActivePowerModel, n::Int, i, f_idx, rate_a)
    @assert(:ratea ∈ ref(pm, n, :slack_types))
    
    p_fr_sec = var(pm, n, :p_sec, f_idx)

    p_max = ref(pm, n, :p_max)
    
    (l,c),i,j = f_idx
    
    if is_branch_monitored(pm, n, l)
        # bound
        rate_a_slack = var(pm, n, :rate_a_slack_sec, f_idx)
        JuMP.@constraint(pm.model, p_fr_sec <=  rate_a + rate_a_slack, base_name = "$(n)_C$(c)_rating_from_ub_$(l)")
        JuMP.@constraint(pm.model, p_fr_sec >= -rate_a - rate_a_slack, base_name = "$(n)_C$(c)_rating_from_lb_$(l)")
    end
    if is_branch_switchable(pm, n, l)
        # switching
        z = var(pm, n, :z_branch, l)
        JuMP.@constraint(pm.model, p_fr_sec <=  p_max * z, base_name = "$(n)_C$(c)_rating_from_on_off_$(l)")
        JuMP.@constraint(pm.model, p_fr_sec >= -p_max * z, base_name = "$(n)_C$(c)_rating_from_on_off_$(l)")
    end
end

""
function constraint_thermal_limit_to_slack(pm::AbstractActivePowerModel, n::Int, i, t_idx, rate_a)
    @assert(:ratea ∈ ref(pm, n, :slack_types))
    
    p_to = var(pm, n, :p, t_idx)
    
    if is_branch_monitored(pm, n, i)
        rate_a_slack = var(pm, n, :rate_a_slack, t_idx)
        JuMP.@constraint(pm.model, p_to <=  rate_a + rate_a_slack, base_name = "$(n)_rating_to_slack_ub_$(i)")
        JuMP.@constraint(pm.model, p_to >= -rate_a - rate_a_slack, base_name = "$(n)_rating_to_slack_lb_$(i)")
    end
end

""
function constraint_thermal_limit_to_slack_on_off(pm::AbstractActivePowerModel, n::Int, i, t_idx, rate_a)
    @assert(:ratea ∈ ref(pm, n, :slack_types))

    p_to = var(pm, n, :p, t_idx)

    p_max = ref(pm, n, :p_max)

    if is_branch_monitored(pm, n, i)
        # bound
        rate_a_slack = var(pm, n, :rate_a_slack, t_idx)
        JuMP.@constraint(pm.model, p_to <=  rate_a + rate_a_slack, base_name = "$(n)_rating_to_slack_ub_$(i)")
        JuMP.@constraint(pm.model, p_to >= -rate_a - rate_a_slack, base_name = "$(n)_rating_to_slack_lb_$(i)")
    end
    if is_branch_switchable(pm, n, i)
        # switching
        z = var(pm, n, :z_branch, i)
        JuMP.@constraint(pm.model, p_to <=  p_max * z, base_name = "$(n)_rating_to_slack_on_off_ub_$(i)")
        JuMP.@constraint(pm.model, p_to >= -p_max * z, base_name = "$(n)_rating_to_slack_on_off_lb_$(i)")
    end
end

""
function constraint_thermal_limit_to_sec(pm::AbstractActivePowerModel, n::Int, i::Int, t_idx, rate_a)
    p_to_sec = var(pm, n, :p_sec, t_idx)
    if isa(p_to_sec, JuMP.VariableRef) && JuMP.has_lower_bound(p_to_sec)
        cstr = JuMP.LowerBoundRef(p_to_sec)
        JuMP.lower_bound(p_to_sec) < -rate_a && JuMP.set_lower_bound(p_to_sec, -rate_a)
        if JuMP.has_upper_bound(p_to_sec)
            JuMP.upper_bound(p_to_sec) >  rate_a && JuMP.set_upper_bound(p_to_sec,  rate_a)
        end
    else
        cstr = JuMP.@constraint(pm.model, p_to_sec <= rate_a)
    end

    if _IM.report_duals(pm)
        sol(pm, n, :branch, t_idx[1])[:mu_sm_to] = cstr
    end
end

""
function constraint_thermal_limit_to_sec_on_off(pm::AbstractActivePowerModel, n::Int, i, t_idx, rate_a)
    p_to_sec = var(pm, n, :p_sec, t_idx)
    z = var(pm, n, :z_branch, i)

    (l,c),i,j = t_idx

    JuMP.@constraint(pm.model, p_to_sec <=  rate_a * z, base_name = "$(n)_C$(c)_rating_to_on_off_ub_$(l)")
    JuMP.@constraint(pm.model, p_to_sec >= -rate_a * z, base_name = "$(n)_C$(c)_rating_to_on_off_lb_$(l)")
end

""
function constraint_thermal_limit_to_slack_sec(pm::AbstractActivePowerModel, n::Int, i, t_idx, rate_a)
    @assert(:ratea ∈ ref(pm, n, :slack_types))

    p_to_sec = var(pm, n, :p_sec, t_idx)

    (l,c),i,j = t_idx

    if is_branch_monitored(pm, n, l)
        rate_a_slack = var(pm, n, :rate_a_slack_sec, t_idx)
        JuMP.@constraint(pm.model, p_to_sec <=  rate_a + rate_a_slack, base_name = "$(n)_C$(c)_rating_to_slack_ub_$(l)")
        JuMP.@constraint(pm.model, p_to_sec >= -rate_a - rate_a_slack, base_name = "$(n)_C$(c)_rating_to_slack_lb_$(l)")
    end
end

""
function constraint_thermal_limit_to_slack_sec_on_off(pm::AbstractActivePowerModel, n::Int, i, t_idx, rate_a)
    @assert(:ratea ∈ ref(pm, n, :slack_types))
    
    p_to_sec = var(pm, n, :p_sec, t_idx)
    
    p_max = ref(pm, n, :p_max)

    (l,c),i,j = t_idx
    
    if is_branch_monitored(pm, n, l)
        # bound
        rate_a_slack = var(pm, n, :rate_a_slack_sec, t_idx)
        JuMP.@constraint(pm.model, p_to_sec <=  rate_a + rate_a_slack, base_name = "$(n)_C$(c)_rating_to_slack_ub_$(l)")
        JuMP.@constraint(pm.model, p_to_sec >= -rate_a - rate_a_slack, base_name = "$(n)_C$(c)_rating_to_slack_lb_$(l)")
    end
    if is_branch_switchable(pm, n, l)
        # switching
        z = var(pm, n, :z_branch, l)
        JuMP.@constraint(pm.model, p_to_sec <=  p_max * z, base_name = "$(n)_C$(c)_rating_to_slack_on_off_$(l)")
        JuMP.@constraint(pm.model, p_to_sec >= -p_max * z, base_name = "$(n)_C$(c)_rating_to_slack_on_off_$(l)")
    end
end

"apo models ignore reactive power flows"
function variable_gen_power_imaginary_sec(pm::AbstractActivePowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, report::Bool=true, kwargs...)
    # already reported in `variable_gen_power_imaginary`
end

"apo models ignore reactive power flows"
function variable_gen_power_imaginary_sec_on_off(pm::AbstractActivePowerModel, contingencies::Vector{Contingency}; nw::Int=nw_id_default, report::Bool=true, kwargs...)
    # already reported in `variable_gen_power_imaginary`
end

""
function constraint_power_balance_sec(pm::AbstractActivePowerModel, n::Int, i::Int, c::ContingencyID, contingency::Contingency, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    p_sec    = get(var(pm, n),    :p_sec, Dict()); _check_var_keys(p_sec, bus_arcs, "active power", "branch")
    pg_sec   = get(var(pm, n),   :pg_sec, Dict()); _check_var_keys(pg_sec, bus_gens, "active power", "generator_sec")
    ps   = get(var(pm, n),   :ps, Dict()); _check_var_keys(ps, bus_storage, "active power", "storage")
    psw  = get(var(pm, n),  :psw, Dict()); _check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    p_dc_sec = get(var(pm, n), :p_dc_sec, Dict()); _check_var_keys(p_dc_sec, bus_arcs_dc, "active power", "dcline")
    if :pd in get(pm.data, "slack_types", [])
        pd_slack_lower = get(var(pm, n), :pd_slack_lower, Dict()); _check_var_keys(pd_slack_lower, bus_arcs_dc, "pd slack lower", "load")
        pd_slack_upper = get(var(pm, n), :pd_slack_upper, Dict()); _check_var_keys(pd_slack_upper, bus_arcs_dc, "pd slack upper", "load")
        bus_loads = [l_id for l_id in ids(pm, n, :load) if ref(pm, n, :load, l_id)["load_bus"] == i]
    else
        pd_slack_lower = Dict()
        pd_slack_upper = Dict()
        bus_loads = []
    end

    # TODO these contingecies should be filtered out in advance
    if length([p_sec[((l,c),i,j)] for (l,i,j) in bus_arcs if l ∉ get(contingency, :branch, [])]) == 0
        return
    end

    # not implemented
    @assert length(ps) == length(psw) == 0

    cstr = JuMP.@constraint(pm.model,
        sum(p_sec[((l,c),i,j)] for (l,i,j) in bus_arcs if l ∉ get(contingency, :branch, []))
        + sum(p_dc_sec[((l_dc,c),i_dc,j_dc)] for (l_dc,i_dc,j_dc) in bus_arcs_dc if l_dc ∉ get(contingency, :dcline, []))
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        ==
        sum(pg_sec[(g, c)] for g in bus_gens if g ∉ get(contingency, :gen, []))
        - sum(ps_sec[(s,c)] for s in bus_storage if s ∉ get(contingency, :storage, []))
        - sum(pd for pd in values(bus_pd))
        - sum(pd_slack_upper[load] for load in bus_loads)
        + sum(pd_slack_lower[load] for load in bus_loads)
        - sum(gs for gs in values(bus_gs))*1.0^2,
        base_name = "$(n)_C$(c)_power_balance_$(i)"
    )

    if _IM.report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr
        sol(pm, n, :bus, i)[:lam_kcl_i] = NaN
    end
end

"ensure that pg is only != 0 if the gen is on and that pg is exactly pg_0 + pg_delta for all gens"
function constraint_gen_power_real_on_off(pm::AbstractActivePowerModel, n::Int, i::Int, pmin, pmax)
    pg = var(pm, n, :pg, i)
    pmax_slack = get(var(pm, n), :pmax_slack, Dict(i => 0))[i]
    pmin_slack = get(var(pm, n), :pmin_slack, Dict(i => 0))[i]
    if pmax_slack == 0
        pg_big_m_mult = 1
    else
        pg_big_m_mult = PARAMS[:BIG_M_MULT][:PG]  
    end
    z_gen = var(pm, n, :z_gen, i)
    
    # reconfigure bounds to allow for switching off generators
    JuMP.set_lower_bound(pg, min(0, ref(pm, n, :gen, i)["pmin"]))
    JuMP.set_upper_bound(pg, pmax * pg_big_m_mult)

    # generation has to be between pmax and pmin if the gen is on, and else pg=0
    JuMP.@constraint(pm.model, pg <= z_gen * (pmax + pmax_slack), base_name = "$(n)_pg_on_off_ub_$(i)")
    JuMP.@constraint(pm.model, pg >= z_gen * (pmin - pmin_slack), base_name = "$(n)_pg_on_off_lb_$(i)")
end

""
function constraint_gen_power_min_switching(pm::AbstractActivePowerModel, n::Int, i::Int, pg_0::Float64, RDP_minus::Float64)
    z_gen = var(pm, n, :z_gen, i)

    if pg_0 > RDP_minus
        # gen can only be switched off if its pg_delta_minus can got to 0
        JuMP.@constraint(pm.model, z_gen == 1, base_name = "$(n)_pg_delta_min_switching_$(i)")
    end

end

"each virtual gen i with factor factor[i] of this GSK has to produce exactly as much as pg[j] / factor[j] * factor[i] for each j ≠ i with factor[j] of this GSK"
function constraint_gen_power_gsk_linking_real(pm::AbstractActivePowerModel, n::Int, gsk_name::String, gsk_factors::Dict{Int, Float64})
    for (i, gsk_factor) in gsk_factors
        other_gens = [j for (j, factor) in gsk_factors if j != i]
        pg_delta_plus = var(pm, n, :pg_delta_plus, i)
        pg_delta_minus = var(pm, n, :pg_delta_minus, i)
        JuMP.@constraint(pm.model, pg_delta_plus == (pg_delta_plus + sum(var(pm, n, :pg_delta_plus, j) for j in other_gens)) * gsk_factor, base_name = "$(n)_pg_gsk_linking_plus_$(gsk_name)_$(i)")
        JuMP.@constraint(pm.model, pg_delta_minus == (pg_delta_minus + sum(var(pm, n, :pg_delta_minus, j) for j in other_gens)) * gsk_factor, base_name = "$(n)_pg_gsk_linking_minus_$(gsk_name)_$(i)")
    end 
end

# Modified from PowerModels.jl (https://github.com/lanl-ansi/PowerModels.jl), see LICENSE
""
function constraint_power_balance(pm::AbstractActivePowerModel, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    p    = get(var(pm, n),    :p, Dict()); _check_var_keys(p, bus_arcs, "active power", "branch")
    pg   = get(var(pm, n),   :pg, Dict()); _check_var_keys(pg, bus_gens, "active power", "generator")
    ps   = get(var(pm, n),   :ps, Dict()); _check_var_keys(ps, bus_storage, "active power", "storage")
    psw  = get(var(pm, n),  :psw, Dict()); _check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    p_dc = get(var(pm, n), :p_dc, Dict()); _check_var_keys(p_dc, bus_arcs_dc, "active power", "dcline")
    if :pd in get(pm.data, "slack_types", [])
        pd_slack_lower = get(var(pm, n), :pd_slack_lower, Dict()); _check_var_keys(pd_slack_lower, bus_arcs_dc, "pd slack lower", "load")
        pd_slack_upper = get(var(pm, n), :pd_slack_upper, Dict()); _check_var_keys(pd_slack_upper, bus_arcs_dc, "pd slack upper", "load")
        bus_loads = [l_id for l_id in ids(pm, n, :load) if ref(pm, n, :load, l_id)["load_bus"] == i]
    else
        pd_slack_lower = Dict()
        pd_slack_upper = Dict()
        bus_loads = []
    end

    cstr = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(p_dc[a_dc] for a_dc in bus_arcs_dc)
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        >=
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pd for pd in values(bus_pd))
        - sum(pd_slack_upper[load] for load in bus_loads)
        + sum(pd_slack_lower[load] for load in bus_loads)
        - sum(gs for gs in values(bus_gs))*1.0^2,
        base_name = "$(n)_power_balance_ub_$(i)"
    )

    cstr = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(p_dc[a_dc] for a_dc in bus_arcs_dc)
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        <=
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pd for pd in values(bus_pd))
        - sum(pd_slack_upper[load] for load in bus_loads)
        + sum(pd_slack_lower[load] for load in bus_loads)
        - sum(gs for gs in values(bus_gs))*1.0^2,
        base_name = "$(n)_power_balance_lb_$(i)"
    )

    if _IM.report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr
        sol(pm, n, :bus, i)[:lam_kcl_i] = NaN
    end
end

# Modified from PowerModels.jl (https://github.com/lanl-ansi/PowerModels.jl), see LICENSE
""
function constraint_thermal_limit_from_on_off(pm::AbstractActivePowerModel, n::Int, i, f_idx, rate_a)
    p_fr = var(pm, n, :p, f_idx)
    z = var(pm, n, :z_branch, i)

    JuMP.@constraint(pm.model, p_fr <=  rate_a*z, base_name = "$(n)_rating_from_on_off_ub_$(i)")
    JuMP.@constraint(pm.model, p_fr >= -rate_a*z, base_name = "$(n)_rating_from_on_off_lb_$(i)")
end

""
function constraint_branch_from_on_off(pm::AbstractActivePowerModel, n::Int, i, f_idx, rate_a)
    p_fr = var(pm, n, :p, f_idx)

    p_max = ref(pm, n, :p_max)

    if is_branch_switchable(pm, n, i)
        # switching
        z = var(pm, n, :z_branch, i)
        JuMP.@constraint(pm.model, p_fr <= p_max * z, base_name = "$(n)_branch_from_on_off_ub_$(i)")
        JuMP.@constraint(pm.model, p_fr >= -p_max * z, base_name = "$(n)_branch_from_on_off_lb_$(i)")
    end
end

# Modified from PowerModels.jl (https://github.com/lanl-ansi/PowerModels.jl), see LICENSE
""
function constraint_thermal_limit_to_on_off(pm::AbstractActivePowerModel, n::Int, i, t_idx, rate_a)
    p_to = var(pm, n, :p, t_idx)
    z = var(pm, n, :z_branch, i)

    JuMP.@constraint(pm.model, p_to <=  rate_a*z, base_name = "$(n)_rating_to_ub_$(i)")
    JuMP.@constraint(pm.model, p_to >= -rate_a*z, base_name = "$(n)_rating_to_lb_$(i)")
end

""
function constraint_branch_to_on_off(pm::AbstractActivePowerModel, n::Int, i, t_idx, rate_a)
    p_to = var(pm, n, :p, t_idx)

    p_max = ref(pm, n, :p_max)

    if is_branch_switchable(pm, n, i)
        # switching
        z = var(pm, n, :z_branch, i)
        JuMP.@constraint(pm.model, p_to <=  p_max * z, base_name = "$(n)_branch_to_on_off_ub_$(i)")
        JuMP.@constraint(pm.model, p_to >= -p_max * z, base_name = "$(n)_branch_to_on_off_lb_$(i)")
    end
end

"ensure that between two timesteps the delta between the pg's is not too large"
function constraint_gen_power_change_linking(pm, nw::Int, i::Int)
    gen         = ref(pm, nw, :gen, i)
    pg_current  = var(pm, nw, :pg, i)
    pg_previous = nw > 1 ? var(pm, nw-1, :pg, i) : return #get(gen, "pg_previous", gen["pg_0"])

    if ! haskey(gen, "ramp_60_pos")
        @warn("Generator $i has no value for `ramp_60_pos`. Creating no constraints for its ramp.")
        return
    else
        ramp_positive = gen["ramp_60_pos"]
        ramp_negative = gen["ramp_60_neg"]
        #@info("Generator $(gen["index"]) has ramp $ramp_positive and - $ramp_negative")
    end

    # ensure that between timesteps, ramps are fulfilled
    JuMP.@constraint(pm.model, pg_previous - pg_current  <= ramp_negative, base_name = "$(nw)_pg_ramp_neg_linking_$(i)")
    JuMP.@constraint(pm.model, pg_current  - pg_previous <= ramp_positive, base_name = "$(nw)_pg_ramp_pos_linking_$(i)")
end


"ensure that tr_shift[i] is not exceeding limits"
function constraint_transformer_shift_change(pm, i::Int; nw::Int=nw_id_default, tral::Real=PARAMS[:LIMITS][:TRANSFORMER][:TRAL_COMPONENT_AND_TIMESTEP])
    if ! is_branch_shiftable(pm, nw, i) || ! (tral < Inf)  
        return
    end
    tr_shift          = var(pm, nw, :tr_shift, i) # in "steps"
    tr_shift_previous = nw > 1 ? var(pm, nw-1, :tr_shift, i) : get(ref(pm, nw, :branch, i), "previous_shift_step", 0) # in "steps"

    JuMP.@constraint(pm.model, tr_shift <= tr_shift_previous + tral, base_name = "$(nw)_tr_shift_ub_$(i)")
    JuMP.@constraint(pm.model, tr_shift >= tr_shift_previous - tral, base_name = "$(nw)_tr_shift_lb_$(i)")
end