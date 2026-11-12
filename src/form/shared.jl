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

"`t[ref_bus] == 0`"
function constraint_theta_ref_sec(pm::AbstractPolarModels, n::Int, i::Int, c::ContingencyID, contingency::Dict)
    @assert(get(contingency, :bus, "") ≠ i, "The reference bus cannot be a contingency")
    JuMP.@constraint(pm.model, var(pm, n, :va_sec)[(i,c)] == 0, base_name = "$(n)_$(c)_theta_ref_$(i)")
end

"Do nothing, no way to represent this in these variables"
function constraint_theta_ref_sec(pm::AbstractWModels, n::Int, ref_bus::Int, c::ContingencyID, contingency::Dict)
end

""
function constraint_power_balance_sec(pm::AbstractWModels, n::Int, i, c::ContingencyID, contingency::Contingency, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    w_sec    = var(pm, n, :w_sec, (i,c))
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
        - sum(gs for gs in values(bus_gs))*w_sec
    )
    cstr_q = JuMP.@constraint(pm.model,
        sum(q_sec[((l,c),i,j)] for (l,i,j) in bus_arcs if l ∉ get(contingency, :branch, []))
        + sum(q_dc_sec[((l_dc,c),i_dc,j_dc)] for (l_dc,i_dc,j_dc) in bus_arcs_dc if l_dc ∉ get(contingency, :dcline, []))
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        ==
        sum(qg_sec[(g, c)] for g in bus_gens if g ∉ get(contingency, :gen, []))
        - sum(qs[(s,c)] for s in bus_storage if s ∉ get(contingency, :storage, []))
        - sum(qd for qd in values(bus_qd))
        + sum(bs for bs in values(bus_bs))*w_sec
    )

    if _IM.report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end

"the total load needs to be between the minimum and maximum amount of power generated by all switched on gens, ie sum(pmin * zgen) ≤ sum(pd) ≤ sum(pmax * zgen)"
function constraint_power_balance_switching(pm::AbstractPowerModel; nw::Int=nw_id_default)
    z_gen = var(pm, nw, :z_gen)
    gens   = ref(pm, nw,   :gen)
    pd_sum = sum([load["pd"] for (_, load) in ref(pm, nw, :load)])

    JuMP.@constraint(pm.model, sum((gen["pg_0"] - gen["RDP-"]) * (is_gen_switchable(gen) ? z_gen[i] : gen["status"]) for (i, gen) in gens) ≤ pd_sum, base_name="$(nw)_constraint_power_balance_switching_min")
    JuMP.@constraint(pm.model, sum((gen["pg_0"] + gen["RDP+"]) * (is_gen_switchable(gen) ? z_gen[i] : gen["status"]) for (i, gen) in gens) ≥ pd_sum, base_name="$(nw)_constraint_power_balance_switching_max")
end

"""
```
t[f_bus] - t[t_bus] <= angmax
t[f_bus] - t[t_bus] >= angmin
```
"""
### Only adds constraint name
function constraint_voltage_angle_difference(pm::AbstractPolarModels, n::Int, f_idx, angmin, angmax)
    i, f_bus, t_bus = f_idx

    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)

    JuMP.@constraint(pm.model, va_fr - va_to <= angmax, base_name = "$(n)_vad_ub_$(i)")
    JuMP.@constraint(pm.model, va_fr - va_to >= angmin, base_name = "$(n)_vad_lb_$(i)")
end

"Link the gen `i`s power to the `contingency`'s activation."
function constraint_gen_power_real_sec_complete(pm, i::Int64, c::ContingencyID, contingency::Contingency, pmin, pmax; n::Int=nw_id_default)
    pg_sec = get(var(pm, n), :pg_sec, Dict())#; _check_var_keys(pg_sec, bus_gens, "active power", "generator_sec")
    #m_contingency = get(var(pm, n), :m_contingency, Dict())
    m = var(pm, n, :m_contingency, c)

    if i in get(contingency, :gen, [])
        # if the contingency is active, the gen produces no power, otherwise it needs to follow its usual bounds
        JuMP.@constraint(pm.model, pg_sec[(i,c)] <= (1-m)*pmax)
        JuMP.@constraint(pm.model, pg_sec[(i,c)] >= (1-m)*pmin)
    end
end

"Link the gen `i`s power to the `contingency`'s activation."
function constraint_gen_power_imaginary_sec_complete(pm, i::Int64, c::ContingencyID, contingency::Contingency, qmin, qmax; n::Int=nw_id_default)
    qg_sec = get(var(pm, n), :qg_sec, Dict())#; _check_var_keys(pg_sec, bus_gens, "active power", "generator_sec")
    #m_contingency = get(var(pm, n), :m_contingency, Dict())
    m = var(pm, n, :m_contingency, c)

    if i in get(contingency, :gen, [])
        # if the contingency is active, the gen produces no power, otherwise it needs to follow its usual bounds
        JuMP.@constraint(pm.model, qg_sec[(i,c)] <= (1-m)*qmax)
        JuMP.@constraint(pm.model, qg_sec[(i,c)] >= (1-m)*qmin)
    end
end

# slacks
"wrapper function for linking slack variables to original variables"
function constraint_component_slack(pm::AbstractPowerModel, nw::Int, i::Int, slack_type::Symbol, component::Dict, lb, ub, lb2, ub2)
    if     slack_type == :pg; constraint_gen_slack_real(             pm, nw, i, lb, ub);
    elseif slack_type == :qg; constraint_gen_slack_imaginary(        pm, nw, i, lb, ub);
    elseif slack_type == :dclinepf; constraint_dcline_slack(         pm, nw, i, component, lb, ub, lb2, ub2);
    elseif slack_type == :vm; constraint_bus_voltage_magnitude_slack(pm, nw, i, component, lb, ub);
    elseif slack_type == :ratea; # do nothing, this is handled in `constraint_thermal_limit*`
    elseif slack_type == :vad;   # do nothing, this is handled in `constraint_ohms_yt*`
    elseif slack_type == :pd; # do nothing, this is handled in `constraint_power_balance*`
    elseif slack_type == :qd;   # do nothing, this is handled in `constraint_power_balance*`
    else
        @warn("Slack type $slack_type not considered")
    end
end

"links slack variables to original variables for gens"
function constraint_gen_slack_real(pm, n::Int, i::Int, pmin, pmax)
    pg = get(var(pm, n), :pg, Dict())
    pmax_slack = get(var(pm, n), :pmax_slack, Dict())
    pmin_slack = get(var(pm, n), :pmin_slack, Dict())
    
    if PARAMS[:LIMITS][:GENERATOR][:GEN_SWITCHING] && is_gen_switchable(ref(pm, n, :gen, i))
        z_gen = get(var(pm, n), :z_gen, Dict())[i]
    else
        z_gen = ref(pm, n, :gen, i)["status"]
    end

    # if on, it can produce nearly as much as it wants, else 0
    JuMP.@constraint(pm.model, pg[i] <= z_gen * (pmax * PARAMS[:BIG_M_MULT][:PG]  ))
    JuMP.@constraint(pm.model, pg[i] >= 0)

    # if on, a slack is required if pmax/pmin is exceeded
    JuMP.@constraint(pm.model, pg[i] <= pmax + pmax_slack[i])
    JuMP.@constraint(pm.model, pg[i] >= z_gen * pmin - pmin_slack[i])
end

"links slack variables to original variables for dclines"
function constraint_dcline_slack(pm, n::Int, i::Int, dcline, pminf, pmaxf, pmint, pmaxt)
    p_dc = get(var(pm, n), :p_dc, Dict())
    p_dc_slack_upper = get(var(pm, n), :p_dc_slack_upper, Dict())
    p_dc_slack_lower = get(var(pm, n), :p_dc_slack_lower, Dict())

    JuMP.@constraint(pm.model, p_dc[dcline] <= pmaxf + p_dc_slack_upper[dcline])
    JuMP.@constraint(pm.model, p_dc[dcline] >= pminf - p_dc_slack_lower[dcline])

    JuMP.@constraint(pm.model, p_dc[dcline] <= pmaxt + p_dc_slack_upper[dcline])
    JuMP.@constraint(pm.model, p_dc[dcline] >= pmint - p_dc_slack_lower[dcline])
end

"links binary slack variables to normal slack variables"
function constraint_component_slack_on_off(pm::AbstractPowerModel, nw::Int, i::Int, slack_type::Symbol, f_idx, lb, ub, lb2, ub2; slack_mults = SLACK_MULTS)
    var_symbols        = get_slack_var_symbols(       pm, [slack_type])
    var_symbols_on_off = get_slack_var_symbols_on_off(pm, [slack_type])

    for (s, var_symbol) in enumerate(var_symbols)
        var_on_off_symbol = var_symbols_on_off[s]

        slack        = get(var(pm, nw), var_symbol,        Dict())
        slack_on_off = get(var(pm, nw), var_on_off_symbol, Dict())

        # no variable of this type
        if isempty(slack)
            @assert(isempty(slack_on_off))
            continue
        end

        # Choose bound by using the name of the slack to determine whether this is upper or lower.
        # Note that we usually would have to add two constraints for dcline, however this is just an approximation.
        if slack_type ∈ [:pg, :qg, :ratea, :dclinepf, :vm, :pd, :vad]
            bound = occursin("lower", "$var_symbol") ? max(lb, lb2) : max(ub, ub2)
            max_slack_mult = slack_mults[slack_type]
        else
            # original bound unclear
            @warn("Original bound unclear!")
            max_slack_mult = 100
            bound = 1
        end
        
        JuMP.@constraint(pm.model, slack[f_idx] <= slack_on_off[f_idx] * max_slack_mult * (1+bound))
    end
end

"wrapper function for linking slack variables to original variables for security constraints"
function constraint_component_slack_sec(pm::AbstractPowerModel, nw::Int, i::Int, c::ContingencyID, contingency::Contingency, slack_type::Symbol, component::Dict, lb, ub, lb2, ub2)
    if     slack_type == :pg; constraint_gen_slack_sec_real(             pm, nw, i, c, contingency, lb, ub); 
    elseif slack_type == :qg; constraint_gen_slack_sec_imaginary(        pm, nw, i, c, contingency, lb, ub); 
    elseif slack_type == :dclinepf; constraint_dcline_slack_sec(         pm, nw, i, c, contingency, component, lb, ub, lb2, ub2);     
    elseif slack_type == :vm; constraint_bus_voltage_magnitude_slack_sec(pm, nw, i, c, contingency, component, lb, ub);
    elseif slack_type == :ratea; # do nothing, this is handled in `constraint_thermal_limit*`
    elseif slack_type == :vad;   # do nothing, this is handled in `constraint_ohms_yt*`
    else
        @warn("Slack type $slack_type not considered")
    end
end

"links slack variables to original variables for gens for security constraints"
function constraint_gen_slack_sec_real(pm, n::Int, i::Int, c::ContingencyID, contingency::Contingency, pmin, pmax)
    pg_sec = get(var(pm, n), :pg_sec, Dict())
    pmax_slack = get(var(pm, n), :pmax_slack, Dict())
    pmin_slack = get(var(pm, n), :pmin_slack, Dict())

    if PARAMS[:LIMITS][:GENERATOR][:GEN_SWITCHING] && is_gen_switchable(ref(pm, n, :gen, i))
        z_gen = get(var(pm, n), :z_gen, Dict())[i]
    else
        z_gen = 1.0
    end

    # if on, it can produce nearly as much as it wants, else 0
    JuMP.@constraint(pm.model, pg_sec[(i,c)] <= z_gen * (pmax * PARAMS[:BIG_M_MULT][:PG]  ))
    JuMP.@constraint(pm.model, pg_sec[(i,c)] >= 0)

    # if on, a slack is required if pmax/pmin is exceeded
    JuMP.@constraint(pm.model, pg_sec[(i,c)] <= pmax + pmax_slack[i])
    JuMP.@constraint(pm.model, pg_sec[(i,c)] >= z_gen * pmin - pmin_slack[i])
end

"links slack variables to original variables for branches for security constraints"
function constraint_branch_slack_sec(pm::AbstractACPModel, n::Int, i::Int, c::ContingencyID, contingency::Contingency, f_idx, flow_lb, flow_ub)
    # handled via thermal limit
end

"links slack variables to original variables for dclines for security constraints"
function constraint_dcline_slack_sec(pm, n::Int, i::Int, arc, c::ContingencyID, contingency::Contingency, pminf, pmaxf, pmint, pmaxt)
    p_dc_sec = get(var(pm, n), :p_dc_sec, Dict())
    p_dc_slack_upper = get(var(pm, n), :p_dc_slack_upper, Dict())
    p_dc_slack_lower = get(var(pm, n), :p_dc_slack_lower, Dict())

    l,i,j = arc

    JuMP.@constraint(pm.model, p_dc_sec[((l,c),i,j)] <= pmaxf + p_dc_slack_upper[((l,c),i,j)])
    JuMP.@constraint(pm.model, p_dc_sec[((l,c),i,j)] >= pminf - p_dc_slack_lower[((l,c),i,j)])

    JuMP.@constraint(pm.model, p_dc_sec[((l,c),i,j)] <= pmaxt + p_dc_slack_upper[((l,c),i,j)])
    JuMP.@constraint(pm.model, p_dc_sec[((l,c),i,j)] >= pmint - p_dc_slack_lower[((l,c),i,j)])
end

"`t[ref_bus] == 0`"
function constraint_theta_ref(pm::AbstractPolarModels, n::Int, i::Int)
    JuMP.@constraint(pm.model, var(pm, n, :va)[i] == 0, base_name = "$(n)_theta_ref_$(i)")
end

"links power gen variable to pg_delta variables"
function constraint_gen_power_delta_linking(pm::AbstractPowerModel, n::Int, i::Int, pg_0)
    pg = var(pm, n, :pg, i)
    pg_delta_plus   = var(pm, n, :pg_delta_plus, i)
    pg_delta_minus  = var(pm, n, :pg_delta_minus, i)

    # linking that pg is exactly pg_0 + pg_delta for all gens
    JuMP.@constraint(pm.model, pg == pg_0 + pg_delta_plus - pg_delta_minus, base_name = "$(n)_pg_delta_linking_$(i)")
end

"lower bound on nested benders IP cost based on switching states: if an active gen should be switched off, we know the exact redispatch costs already, if it is switched on, we know a lower bound on the redispatch cost already"
function constraint_benders_min_gen_cost(pm, nw::Int)
    if PARAMS[:LIMITS][:GENERATOR][:GEN_SWITCHING]
        gens = ref(pm, nw, :gen)

        JuMP.@constraint(pm.model,
                          sum((1 - var(pm, nw, :z_gen, gen_id)) * (                        gens[gen_id]["pg_0"]  * (gens[gen_id]["cost_change"] - gens[gen_id]["cost_produce"]) * 0.01) for gen_id in axes(var(pm, nw, :z_gen),1) if  Bool(gens[gen_id]["status"]))
                        + sum(     var(pm, nw, :z_gen, gen_id)  * ((gens[gen_id]["pmin"] - gens[gen_id]["pg_0"]) * (gens[gen_id]["cost_change"] + gens[gen_id]["cost_produce"]) * 0.01) for gen_id in axes(var(pm, nw, :z_gen),1) if !Bool(gens[gen_id]["status"]))
                        <= var(pm, nw, :η_rd),
                        base_name = "benders_min_gen_cost"
                        )
    end
end

"per breaker, at most one of its breaker branches may be active at a given time"
function constraint_breaker_branch_limit(pm, i::Int; nw=nw_id_default)
    bus = ref(pm, nw, :bus, i)
    if bus["type"] ≠ "breaker"
        return
    end

    breaker_branches = bus["breaker_branches"]
    @assert(length(breaker_branches) > 0, "Breaker bus $(i) has no breaker branches - what is its purpose?")

    z_branch = var(pm, nw, :z_branch)

    JuMP.@constraint(pm.model, sum(z_branch[j] for j in breaker_branches) <= 1, base_name = "$(nw)_breaker_branch_limit_$(i)")
end

"limits the number of simultaneous topological actions"
function constraint_top_actions_limit(pm; n_top_actions::Int64=PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP], nw=nw_id_default)
    # introduce branch change indicator and link it
    variable_branch_change_indicator(pm, nw=nw)
    sw_br_ids = [i for (i, branch) in ref(pm, nw, :branch) if is_branch_switchable(branch)]
    for i in sw_br_ids
        constraint_branch_change_linking(pm, i, nw=nw)
    end

    # create constraint to limit number of branch changes
    if length(sw_br_ids) > 0
        z_branch_change = var(pm, :z_branch_change)
        JuMP.@constraint(pm.model, sum(z_branch_change) <= n_top_actions)
    end
end

"ensure that z_branch_change[i] is 1 iff the value of z_branch[i] changed in the current time step"
function constraint_branch_change_linking(pm, i::Int; nw::Int=nw_id_default)
    z_branch_change     = var(pm, nw, :z_branch_change, i)
    z_branch_current    = var(pm, nw, :z_branch, i)
    z_branch_previous   = nw > 1 ? var(pm, nw-1, :z_branch, i) : get(ref(pm, nw, :branch, i), "status_previous", get(ref(pm, nw, :branch, i), "status", ref(pm, nw, :branch, i)["br_status"]))

    JuMP.@constraint(pm.model, -z_branch_change <= z_branch_previous - z_branch_current                   , base_name = "$(nw)_br_change_linking_ub_$(i)")
    JuMP.@constraint(pm.model,                     z_branch_previous - z_branch_current <= z_branch_change, base_name = "$(nw)_br_change_linking_lb_$(i)")
end

"ensure that tr_shift_change[i] is 1 iff the value of tr_shift[i] changed in the current time step"
function constraint_transformer_shift_change_linking(pm, i::Int; nw::Int=nw_id_default, tral::Real=PARAMS[:LIMITS][:TRANSFORMER][:TRAL_COMPONENT_AND_TIMESTEP])
    tr_shift_change     = get(var(pm, nw), :tr_shift_change, Dict(i => 1.0))[i]
    tr_shift_current    = var(pm, nw, :tr_shift, i)
    tr_shift_previous   = nw > 1 ? var(pm, nw-1, :tr_shift, i) : get(ref(pm, nw, :branch, i), "previous_shift_step", 0)

    operational_gradient = min(tral, 1 + 2 * ref(pm, nw, :branch, i)["n_shift_steps"]) # shift at max tral, or the bound of tr_shift
    
    JuMP.@constraint(pm.model, tr_shift_current  <= tr_shift_previous + tr_shift_change * operational_gradient, base_name = "$(nw)_tr_shift_linking_ub_$(i)")
    JuMP.@constraint(pm.model, tr_shift_current  >= tr_shift_previous - tr_shift_change * operational_gradient, base_name = "$(nw)_tr_shift_linking_lb_$(i)")
end

"ensure that z_gen_shutdown[i]/z_gen_startup[i] is 1 iff the gen is shutdown/started up in the current time step"
function constraint_gen_change_linking(pm, nw::Int, i::Int)
    gen_current = ref(pm, nw, :gen, i)

    # don't create the constraint if everything is fixed
    if !is_gen_switchable(gen_current) && (nw <= 1 || !is_switchable(ref(pm, nw-1, :gen, i)))
        return
    end

    if is_gen_switchable(gen_current)
        z_gen_current  = var(pm, nw, :z_gen, i)
        z_gen_shutdown = var(pm, nw, :z_gen_shutdown, i)
        z_gen_startup  = var(pm, nw, :z_gen_startup, i)
    else
        @assert gen_current["status"] == 1

        z_gen_current  = 1
        z_gen_shutdown = 0
        z_gen_startup  = 0
    end

    if nw <= 1
        z_gen_previous = gen_current["status_previous"]
    elseif is_switchable(ref(pm, nw-1, :gen, i)) # FIXME: if gen was off and not switchable in the previous timestep, there will not be a ref(pm, nw-1, :gen, i)
        z_gen_previous = var(pm, nw-1, :z_gen, i)
    else
        z_gen_previous = gen["status"]
    end

    # a gen's status is defined by the previous status and shutdown or startup
    JuMP.@constraint(pm.model, z_gen_current == z_gen_previous + z_gen_startup - z_gen_shutdown, base_name = "$(nw)_gen_status_change_linking_$(i)")
end

"limit the number of total topological changes per time step"
function constraint_branch_max_changes_timestep(pm, nw::Int; tal_timestep::Real=PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP] )
    if tal_timestep < length(ref(pm, nw, :branch))
        JuMP.@constraint(pm.model, sum(z_branch_change for z_branch_change in var(pm, nw, :z_branch_change)) <= tal_timestep, base_name = "$(nw)_br_change_max_hour")
        #JuMP.@constraint(pm.model, sum(z_branch_change for z_branch_change in var(pm, nw, :z_branch_change)) >= 1, base_name = "$(nw)_br_change_min_hour") # to debug whether switching works
    end
end

"limit the number of total topological changes over all time steps and components"
function constraint_branch_max_changes_total(pm, number_nws::Int; tal_day::Real=PARAMS[:LIMITS][:BRANCH][:TAL_DAY] )
    if tal_day < Inf && number_nws > 1
        JuMP.@constraint(pm.model, sum(sum(z_branch_change for z_branch_change in var(pm, nw, :z_branch_change)) for nw in 1:number_nws) <= tal_day, base_name = "br_change_max_day")
    end
end

"limit the number of total topological changes over all time steps for one component"
function constraint_branch_max_changes_component(pm, number_nws::Int, i::Int; tal_component_and_day::Real=PARAMS[:LIMITS][:BRANCH][:TAL_COMPONENT_AND_DAY] )
    if tal_component_and_day < Inf && number_nws > 1
        JuMP.@constraint(pm.model, sum(sum(z_branch_change for z_branch_change in var(pm, nw, :z_branch_change, i)) for nw in 1:number_nws) <= tal_component_and_day, base_name = "br_change_max_day_comp_$(i)")
    end
end

"force topological actions to last for a given duration"
function constraint_switch_duration(pm, number_nws::Int, i::Int; remain_switched_for::Int=PARAMS[:LIMITS][:BRANCH][:REMAIN_SWITCHED_FOR] )
    if remain_switched_for > 0
        # TODO: consider time steps before the first
        # in each interval of length `remain_switched_for``, at most one topological actions may happen on this component
        remain_switched_for = min(remain_switched_for, number_nws) # otherwise the list might be empty

        for nw_start in 1:number_nws - remain_switched_for + 1
            JuMP.@constraint(pm.model, sum(var(pm, nw, :z_branch_change, i) for nw in nw_start:(nw_start + remain_switched_for - 1)) <= 1, base_name = "br_change_duration_$(i)_$(nw_start)+")
        end
    end
end

"limit the number of total topological changes per time step"
function constraint_transformer_max_changes_timestep(pm, nw::Int; tral_timestep::Real=PARAMS[:LIMITS][:TRANSFORMER][:TRAL_TIMESTEP] )
    if tral_timestep < Inf
        JuMP.@constraint(pm.model, sum(tr_shift_change for tr_shift_change in var(pm,nw,:tr_shift_change)) <= tral_timestep, base_name = "$(nw)_tr_shift_max_hour")
    end
end

"limit the number of total topological changes over all time steps and components"
function constraint_transformer_max_changes_total(pm, number_nws::Int; tral_day::Real=PARAMS[:LIMITS][:TRANSFORMER][:TRAL_DAY] )
    if tral_day < Inf
        JuMP.@constraint(pm.model, sum(sum(tr_shift_change for tr_shift_change in var(pm,nw,:tr_shift_change)) for nw in 1:number_nws) <= tral_day, base_name = "tr_shift_change_max_day")
    end
end

"limit the number of total topological changes over all time steps for one component"
function constraint_transformer_max_changes_component(pm, number_nws::Int, i::Int; tral_component_and_day::Real=PARAMS[:LIMITS][:TRANSFORMER][:TRAL_COMPONENT_AND_DAY] )
    if tral_component_and_day < Inf
        JuMP.@constraint(pm.model, sum(var(pm,nw,:tr_shift_change,i) for nw in 1:number_nws) <= tral_component_and_day, base_name = "tr_shift_max_day_comp_$(i)")
    end
end

"limit the number of total generator dispatch changes over all time steps and generators"
function constraint_gen_max_changes_total(pm, number_nws::Int; pg_delta_num_day::Real=PARAMS[:LIMITS][:REDISPATCH][:PG_DELTA_NUM_DAY] )
    if pg_delta_num_day < Inf
        JuMP.@constraint(pm.model, sum(sum(pg_delta_on_off for pg_delta_on_off in var(pm, nw, :pg_delta_on_off)) for nw in 1:number_nws) <= pg_delta_num_day, base_name = "pg_delta_max_day")
    end
end

"limit the number of generator shutdowns over all time steps and generators"
function constraint_gen_max_shutdown_total(pm, number_nws::Int; gen_shutdown_num_day::Real=PARAMS[:LIMITS][:GENERATOR][:GEN_SHUTDOWN_NUM_DAY] )
    if gen_shutdown_num_day < Inf
        JuMP.@constraint(pm.model, sum(sum(z_gen_shutdown for z_gen_shutdown in var(pm, nw, :z_gen_shutdown)) for nw in 1:number_nws) <= gen_shutdown_num_day, base_name = "gen_shutdown_max_day")
    end
end

"limit the number of total topological changes over all time steps and components"
function constraint_gen_max_startup_total(pm, number_nws::Int; gen_startup_num_day::Real=PARAMS[:LIMITS][:GENERATOR][:GEN_STARTUP_NUM_DAY] )
    if gen_startup_num_day < Inf
        JuMP.@constraint(pm.model, sum(sum(z_gen_startup for z_gen_startup in var(pm, nw, :z_gen_startup)) for nw in 1:number_nws) <= gen_startup_num_day, base_name = "gen_startup_max_day")
    end
end

"force startup to last for a given duration: sum_{t+1 ≤ t' ≤ t+min_up} z_shutdown[t'] ≤ z_gen[t]"
function constraint_startup_minup(pm, T::Int, i::Int; min_up::Int=PARAMS[:LIMITS][:GENERATOR][:GEN_MIN_UP] )
    # in each interval of length `min_up`, the gen may not be shut down if it was started up just before
    if min_up > 0
        min_up = min(min_up, T) # otherwise the list might be empty

        for t in min_up:T
            JuMP.@constraint(pm.model, sum(var(pm, q, :z_gen_startup, i) for q in (t - min_up + 1):t) <= var(pm, t, :z_gen, i), base_name = "gen_minup_$(i)")
        end
    end
end

"force shutdown to last for a given duration: sum_{t+1 ≤ t' ≤ t+min_down} z_startup[t'] ≤ z_shutdown[t-1]"
function constraint_shutdown_mindown(pm, T::Int, i::Int; min_down::Int=PARAMS[:LIMITS][:GENERATOR][:GEN_MIN_DOWN] )
    # in each interval of length `min_down`, the gen may not be started up if it was shut down just before
    if min_down > 0
        min_down = min(min_down, T) # otherwise the list might be empty

        for t in min_down:T
            JuMP.@constraint(pm.model, sum(var(pm, q, :z_gen_shutdown, i) for q in (t - min_down + 1):t) <= 1 - var(pm, t, :z_gen, i), base_name = "gen_mindown_$(i)")
        end
    end
end

"constraint ensuring that at least one disrespected_contingency is now respected"
function constraint_contingency_coverage(pm, disrespected_contingencies::Vector{Contingency}, contingencies::Vector{Contingency})
    disrespected_contingencies_ids::Vector{Int64} = [contingency_id 
            for (contingency_id, contingency) in enumerate(contingencies) 
            if contingency in disrespected_contingencies
        ]

    ms = var(pm, :m_contingency)
    #print("Adding constraint $(sum(ms[i] for i in disrespected_contingencies_ids)) >= 1")
    JuMP.@constraint(pm.model, sum(ms[i] for i in disrespected_contingencies_ids) >= 1)
end

"links binary gen delta variables to normal gen delta"
function constraint_gen_power_change_on_off(pm::AbstractPowerModel, nw::Int, i::Int)
    pg_delta_plus   = var(pm, nw, :pg_delta_plus, i)
    pg_delta_minus  = var(pm, nw, :pg_delta_minus, i)
    pg_delta_on_off = var(pm, nw, :pg_delta_on_off, i)
    pg_delta_minus_limit = ref(pm, nw, :gen, i)["RDP-"]
    pg_delta_plus_limit  = ref(pm, nw, :gen, i)["RDP+"]

    JuMP.@constraint(pm.model, pg_delta_minus <= pg_delta_minus_limit * pg_delta_on_off, base_name = "$(nw)_pg_delta_minus_on_off_linking_$(i)")
    JuMP.@constraint(pm.model, pg_delta_plus  <= pg_delta_plus_limit  * pg_delta_on_off, base_name = "$(nw)_pg_delta_plus_on_off_linking_$(i)")
end

""
function constraint_power_balance(pm::AbstractWModels, n::Int, i, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    w    = var(pm, n, :w, i)
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
        - sum(gs for gs in values(bus_gs))*w
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
        + sum(bs for bs in values(bus_bs))*w
    )

    if _IM.report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end