"Return a completed standalone restricted solve through the public entry point."
struct RestrictedWarmstartFinished <: Exception
    result::Dict{String,Any}
end

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

"run SCOPF and insert values as warmstart"
function warmstart_scopf!(data, solver_convex, slack_types; fix_uc=false, disable_sc=false, check_contingencies=false)
    if PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP] == 0
        return
    end
    
    #@assert (disable_sc && length(data["nw"]["1"]["contingencies"]) > 0)
    @info("Warmstart: $(disable_sc ? "OPF" : "SCOPF")")
    collect_admittance = !disable_sc && get(PARAMS[:WARMSTARTING][:RSCOTS][:WEIGHTS], :ADMITTANCE_SECANT, 0.0) > 0
    captured = Ref{Any}()
    capture_model(pm) = (captured[] = pm; pm)
    old_time_limit = get(data, "_TimeLimit", nothing)
    old_controls = Dict(k => data[k] for k in ("_MIPGap", "_MIPFocus") if haskey(data, k))
    data["_TimeLimit"] = get(PARAMS[:WARMSTARTING], :SCOPF_TIME_LIMIT, 1800)
    data["_collect_admittance_scores"] = collect_admittance
    # save original params
    TAL_copy = deepcopy(PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP])
    pd_enabled_copy             = deepcopy(PARAMS[:BENDERS][:PARTIAL_DECOMP][:ENABLED])
    pd_include_base_pf_copy     = deepcopy(PARAMS[:BENDERS][:PARTIAL_DECOMP][:INCLUDE_BASE_PF])
    pd_n_max_copy               = deepcopy(PARAMS[:BENDERS][:PARTIAL_DECOMP][:N_MAX_CONTINGENCIES])
    pd_n_max_pct_copy           = deepcopy(PARAMS[:BENDERS][:PARTIAL_DECOMP][:N_MAX_CONTINGENCIES_PCT])
    #TRANSFORMER_SHIFTING_copy = deepcopy(PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_SHIFTING])
    n_threads_benders_copy = deepcopy(PARAMS[:SOLVER][:MULTITHREADING][:N_THREADS_BENDERS])
    disable_sc && (contingency_copy = deepcopy(data["nw"]["1"]["contingencies"]))

    result = nothing
    try
        # set params
        disable_sc && (data["nw"]["1"]["contingencies"] = Vector{Contingency}()) # no SC if disabled
        PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP] = 0 # no TS
        PARAMS[:SOLVER][:MULTITHREADING][:N_THREADS_BENDERS] = 1 # no parallelization
        #PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_SHIFTING] = false # no transformer shifting
        PARAMS[:BENDERS][:PARTIAL_DECOMP][:ENABLED]                 = true  # build PFs in master
        PARAMS[:BENDERS][:PARTIAL_DECOMP][:INCLUDE_BASE_PF]         = true  # calculate base power flow in master
        PARAMS[:BENDERS][:PARTIAL_DECOMP][:N_MAX_CONTINGENCIES]     = 0     # disable contingency-PF / partial decomposition during warmstart
        PARAMS[:BENDERS][:PARTIAL_DECOMP][:N_MAX_CONTINGENCIES_PCT] = 0.0
        refresh_derived_params!()
        data["_MIPGap"] = 0.05 # no need for optimal solution
        data["_MIPFocus"] = 1 # find a good primal solution fast
        #data["_NumericFocus"] = 0 # find a good primal solution fast

        if check_contingencies
            # set rate a for N-1 to `PARAMS[:WARMSTARTING][:RSCOTS][:CHECK_CONTINGENCIES_RATEA_MULT]` to get more overloads
            reduce_max_loading!(data["nw"]["1"], ub_mult=PARAMS[:WARMSTARTING][:RSCOTS][:CHECK_CONTINGENCIES_RATEA_MULT])
        end

        # solve SCOPF and set values as warmstart
        result = solve_secure_ots_time_benders(data, DCPPowerModel, solver_convex, weighted_obj_func_factory([objective_min_num_gen_switch, objective_min_num_tr_switch, objective_min_num_branch_switch, objective_min_redispatch_cost, objective_min_num_redispatch, objective_min_slack_relative_linear, objective_min_num_slacks], [data["cost_coeffs"]["gen_sw"], data["cost_coeffs"]["tr_sw"], data["cost_coeffs"]["br_sw"], data["cost_coeffs"]["redispatch"], data["cost_coeffs"]["num_redispatch"], data["cost_coeffs"]["overload"], data["cost_coeffs"]["num_overload"]]), slack_types=slack_types, solver_specific_build_function=capture_model)
       

        if ! PARAMS[:BENDERS][:BASE_PF_IN_MASTER] && ! check_contingencies
            # calculate power flow for base case
            calculate_power_flow!(data, result, solver_convex)
        elseif check_contingencies
            # calculate power flow for contingencies
            _, result_with_outages, _ = check_relevant_overloads(deepcopy(data), result, solver_convex, for_validation=false)
            result["solution"]["nw"]["1"]["branch"] = result_with_outages["solution"]["nw"]["1"]["branch"]

            # reset rate a to original (true) value
            reset_max_loading!(data["nw"]["1"], ub_mult=PARAMS[:WARMSTARTING][:RSCOTS][:CHECK_CONTINGENCIES_RATEA_MULT])
        end

        if result["termination_status"] == JuMP.MOI.OPTIMAL
            data["result"] = result
            data["WARMSTART_OBJ"] = result["objective"]
            data["WARMSTART_PROVIDED"] = true
            set_dc_pf_start_values!(data, result["solution"], fix_uc=fix_uc, fix_ts=false, fix_benders=false, fix_pf=false, fix_slack=false)

            if collect_admittance
                diagnostic = collect_admittance_scores!(captured[])
                diagnostic["scopf"] = Dict("usable" => true, "status" => string(result["termination_status"]), "objective" => result["objective"])
                diagnostic["scopf_solve_seconds"] = result["solve_time"]
                diagnostic["scopf_formulation"] = "Benders"
                data["admittance_warmstart"] = diagnostic
            end

            @info("Adding sw_heur_scores to data")
            # solve RSCOTS and set values as warmstart
            calc_and_set_br_sw_heur_scores!(data["nw"]["1"], 
                                            set_switchability=false,
                                            nw_sol=result["solution"]["nw"]["1"], 
                                            params=PARAMS[:WARMSTARTING][:RSCOTS], 
                                            heur_params=PARAMS[:WARMSTARTING][:RSCOTS][:HEUR])
        else
            data["WARMSTART_PROVIDED"] = false
        end

    finally
        # revert to original params
        PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP] = TAL_copy
        PARAMS[:BENDERS][:PARTIAL_DECOMP][:ENABLED]                 = pd_enabled_copy
        PARAMS[:BENDERS][:PARTIAL_DECOMP][:INCLUDE_BASE_PF]         = pd_include_base_pf_copy
        PARAMS[:BENDERS][:PARTIAL_DECOMP][:N_MAX_CONTINGENCIES]     = pd_n_max_copy
        PARAMS[:BENDERS][:PARTIAL_DECOMP][:N_MAX_CONTINGENCIES_PCT] = pd_n_max_pct_copy
        refresh_derived_params!()
        #PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_SHIFTING] = TRANSFORMER_SHIFTING_copy
        PARAMS[:SOLVER][:MULTITHREADING][:N_THREADS_BENDERS] = n_threads_benders_copy
        delete!(data, "_MIPGap")
        delete!(data, "_MIPFocus")
        merge!(data, old_controls)
        #delete!(data, "_NumericFocus")
        disable_sc && (data["nw"]["1"]["contingencies"] = contingency_copy)
        delete!(data, "_collect_admittance_scores")
        old_time_limit === nothing ? delete!(data, "_TimeLimit") : (data["_TimeLimit"] = old_time_limit)
    end
    score_time = collect_admittance && haskey(data, "admittance_warmstart") ? data["admittance_warmstart"]["score_solve_seconds"] : 0.0
    data["warmstart_time"] = get(data, "warmstart_time", 0) + result["solve_time"] + get(result, "power_flow_solve_time", 0.0) + score_time
    if collect_admittance && !data["WARMSTART_PROVIDED"]
        data["admittance_warmstart"] = Dict("scopf" => Dict("usable" => false, "status" => string(result["termination_status"])),
            "scopf_formulation" => "Benders", "scopf_solve_seconds" => result["solve_time"], "score_solve_seconds" => 0.0,
            "selected" => Int[])
    end

    if !data["WARMSTART_PROVIDED"]
        @warn("Warmstart SCOPF failed with termination status $(result["termination_status"])")
    end

    return result
end

"run OTS and insert values as warmstart"
function warmstart_ots!(data, solver_convex, slack_types; restricted_TS=false)
    @assert PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP] > 0

    # run OTS: filter set of switchable branches and calculate a good primal solution
    t = @elapsed begin
        # save original params
        TAL_copy = deepcopy(PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP])
        contingency_copy = deepcopy(data["nw"]["1"]["contingencies"])

        # set params
        data["nw"]["1"]["contingencies"] = Vector{Contingency}() # no SC
        PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP] = 1 # limited TS
        data["_MIPGap"] = 0.05 # no need for optimal solution

        # solve OTS and set values as warmstart
        if restricted_TS
            calc_and_set_br_sw_heur_scores!(data["nw"]["1"],
                                            nw_sol=data["result"]["solution"]["nw"]["1"],
                                            params=PARAMS[:WARMSTARTING][:RSCOTS],
                                            heur_params=PARAMS[:WARMSTARTING][:RSCOTS][:HEUR])
        end

        result = solve_secure_ots_time_benders(data, DCPPowerModel, solver_convex, weighted_obj_func_factory([objective_min_num_gen_switch, objective_min_num_tr_switch, objective_min_num_branch_switch, objective_min_redispatch_cost, objective_min_num_redispatch, objective_min_slack_relative_linear, objective_min_num_slacks], [data["cost_coeffs"]["gen_sw"], data["cost_coeffs"]["tr_sw"], data["cost_coeffs"]["br_sw"], data["cost_coeffs"]["redispatch"], data["cost_coeffs"]["num_redispatch"], data["cost_coeffs"]["overload"], data["cost_coeffs"]["num_overload"]]), slack_types=slack_types)
       
        if ! PARAMS[:BENDERS][:BASE_PF_IN_MASTER]
            # calculate power flow for base case
            calculate_power_flow!(data, result, solver_convex)
        end
        if result["termination_status"] == JuMP.MOI.OPTIMAL
            data["WARMSTART_OBJ"] = result["objective"]
            data["WARMSTART_PROVIDED"] = true
            set_dc_pf_start_values!(data, result["solution"], fix_uc=false, fix_ts=false, fix_benders=false, fix_pf=false, fix_slack=false)
        else
            data["WARMSTART_PROVIDED"] = false
            @error("Warmstart OTS failed with termination status $(result["termination_status"])")
        end

        if restricted_TS
            reset_switchability!(data["nw"]["1"])
        end

        # revert to original params
        PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP] = TAL_copy
        delete!(data, "_BestObjStop")
        delete!(data, "_MIPGap")
        data["nw"]["1"]["contingencies"] = contingency_copy
    end
    data["warmstart_time"] = get(data, "warmstart_time", 0) + result["solve_time"] + get(result, "power_flow_solve_time", 0.0)
    return result
end

"Use the existing RSCOTS selection while keeping the full model for score recovery and exchange."
function rscots_stage(data, result, tal)
    nw = data["nw"]["1"]
    config = PARAMS[:WARMSTARTING][:RSCOTS]
    try
        calc_and_set_br_sw_heur_scores!(nw; nw_sol=result["solution"]["nw"]["1"],
            params=config, heur_params=config[:HEUR])
    finally
        reset_switchability!(nw)
    end
    return Dict{String,Any}("tal" => tal,
        "selected" => sort([b["index"] for b in values(nw["branch"]) if b["heur_switchable"]]),
        "retained" => sort([b["index"] for b in values(nw["branch"]) if get(b, "heur_switched", false)]),
        "scores" => Dict(b["index"] => b["sw_heur_score"] for b in values(nw["branch"]) if is_branch_switchable(b)))
end

"Apply stage restrictions to a fresh model without changing the input topology."
function restrict_rscots_stage!(pm, stage)
    selected = Set(stage["selected"])
    for (i, branch) in ref(pm, 1, :branch)
        is_branch_switchable(branch) || continue
        i in selected && continue
        status = branch["status"]
        z = var(pm, 1, :z_branch)[i]
        @assert !is_fixed(z) || fix_value(z) == status "Stage conflicts with an existing branch fixing"
        fix(z, status; force=true)
    end
    return pm
end

"run RSCOTS and insert values as warmstart"
function warmstart_rscots!(data, solver_convex, slack_types, result; TAL=PARAMS[:WARMSTARTING][:RSCOTS][:TAL], check_stall=false,
                          stage=nothing, refresh_scores=false)
     if ! data["WARMSTART_PROVIDED"] || TAL == 0 || TAL > PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP]
        return result # no warmstart provided means RSCOTS cannot be run for any sensible switching heuristic
    end

    if check_stall && haskey(data, "OLD_WARMSTART_OBJ") && data["WARMSTART_OBJ"] >= data["OLD_WARMSTART_OBJ"]
        @info("Skipping RSCOTS as no improvement in WARMSTART_OBJ")
        return result
    end

    previous_result = result
    # keep copy of original params
    TAL_copy = deepcopy(PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP])

    old_controls = Dict(k => data[k] for k in ("_BestObjStop", "_MIPGap", "_TimeLimit") if haskey(data, k))
    captured = Ref{Any}()
    old_collect = get(data, "_collect_admittance_scores", nothing)
    try
        # set params
        PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP] = TAL # set TAL for RSCOTS
        data["OLD_WARMSTART_OBJ"] = data["WARMSTART_OBJ"] # keep old warmstart obj
        # _BestObjStop: terminate the warmstart MIP once the incumbent reaches this fraction of the previous warmstart objective.
        config = PARAMS[:WARMSTARTING][:RSCOTS]
        fraction = config[:BEST_OBJ_STOP_FRACTION]
        data["_BestObjStop"] = fraction > 0 ? fraction * data["WARMSTART_OBJ"] : -1e100
        data["_MIPGap"] = config[:MIP_GAP]
        #data["_MIPFocus"] = 1 # find a good primal solution fast
        #data["_NumericFocus"] = 0 # find a good primal solution fast
        data["_TimeLimit"] = PARAMS[:WARMSTARTING][:RSCOTS][:TIME_LIMIT] # no need for optimal solution

        # solve RSCOTS and set values as warmstart
        if stage === nothing
            calc_and_set_br_sw_heur_scores!(data["nw"]["1"],
                                        nw_sol=result["solution"]["nw"]["1"],
                                        params=PARAMS[:WARMSTARTING][:RSCOTS],
                                        heur_params=PARAMS[:WARMSTARTING][:RSCOTS][:HEUR])
        else
            data["_rscots_stage"] = stage
            data["_collect_admittance_scores"] = refresh_scores
        end

        if haskey(data, "admittance_warmstart")
            diagnostic = data["admittance_warmstart"]
            diagnostic["selected"] = stage === nothing ?
                sort([b["index"] for b in values(data["nw"]["1"]["branch"]) if is_branch_switchable(b)]) : copy(stage["selected"])
            diagnostic["restricted_switchable"] = copy(diagnostic["selected"])
        end
        function record_rscots_model(pm)
            captured[] = pm
            stage === nothing || restrict_rscots_stage!(pm, stage)
            if haskey(data, "admittance_warmstart")
                data["admittance_warmstart"]["rscots_actual_time_limit"] = time_limit_sec(pm.model)
                data["admittance_warmstart"]["rscots_actual_threads"] = get_optimizer_attribute(pm.model, "Threads")
                data["admittance_warmstart"]["rscots_actual_mipgap"] = get_optimizer_attribute(pm.model, "MIPGap")
                data["admittance_warmstart"]["rscots_actual_bestobjstop"] = get_optimizer_attribute(pm.model, "BestObjStop")
            end
            if stage !== nothing
                stage["actual_mipgap"] = get_optimizer_attribute(pm.model, "MIPGap")
                stage["actual_bestobjstop"] = get_optimizer_attribute(pm.model, "BestObjStop")
            end
            return pm
        end

        # initialize cut pool for carry-forward if enabled
        if PARAMS[:WARMSTARTING][:RSCOTS][:CARRY_FORWARD_CUTS]
            data["_benders_cuts_pool"] = SerializedCut[]
        end

        result = solve_secure_ots_time_benders(data, DCPPowerModel, solver_convex, weighted_obj_func_factory([objective_min_num_gen_switch, objective_min_num_tr_switch, objective_min_num_branch_switch, objective_min_redispatch_cost, objective_min_num_redispatch, objective_min_slack_relative_linear, objective_min_num_slacks], [data["cost_coeffs"]["gen_sw"], data["cost_coeffs"]["tr_sw"], data["cost_coeffs"]["br_sw"], data["cost_coeffs"]["redispatch"], data["cost_coeffs"]["num_redispatch"], data["cost_coeffs"]["overload"], data["cost_coeffs"]["num_overload"]]), slack_types=slack_types, solver_specific_build_function=record_rscots_model)

        data["warmstart_time"] = get(data, "warmstart_time", 0.0) + result["solve_time"]
        accepted = result["primal_status"] == MOI.FEASIBLE_POINT && isfinite(result["objective"]) && result["objective"] <= previous_result["objective"]
        if haskey(data, "admittance_warmstart")
            merge!(data["admittance_warmstart"], Dict("rscots_solve_seconds" => result["solve_time"],
                "rscots_objective" => result["objective"], "rscots_status" => string(result["termination_status"]),
                "accepted_rscots" => accepted))
        end
        if !accepted
            @warn("RSCOTS did not improve the feasible warmstart; retaining the preceding solution.")
            result = previous_result
        end

        if ! PARAMS[:BENDERS][:BASE_PF_IN_MASTER]
            # Recover the base power flow and include this additional solver call.
            old_pf_time = get(result, "power_flow_solve_time", 0.0)
            calculate_power_flow!(data, result, solver_convex)
            pf_time = result["power_flow_solve_time"] - old_pf_time
            data["warmstart_time"] += pf_time
            if haskey(data, "admittance_warmstart")
                data["admittance_warmstart"]["power_flow_solve_seconds"] = pf_time
            end
        end

        # filter active cuts and accumulate for carry-forward
        if PARAMS[:WARMSTARTING][:RSCOTS][:CARRY_FORWARD_CUTS]
            active_cuts = filter_active_cuts!(data, result, PARAMS[:WARMSTARTING][:RSCOTS][:CARRY_FORWARD_ACTIVITY_TOL])
            mode = PARAMS[:WARMSTARTING][:RSCOTS][:CARRY_FORWARD_MODE]
            accumulate_carried_cuts!(data, mode)
            n_carried_total = length(get(data, "_carried_benders_cuts", SerializedCut[]))
            @info("RSCOTS cut carry-forward: $(length(active_cuts)) active cuts from this stage, $(n_carried_total) total carried cuts.")
        end

        data["WARMSTART_OBJ"] = result["objective"]
        data["WARMSTART_PROVIDED"] = true
        set_dc_pf_start_values!(data, result["solution"], fix_uc=false, fix_ts=false, fix_tr=false, fix_benders=false, fix_pf=false, fix_slack=false)
        if stage !== nothing
            topology = Dict(b["index"] => round(Int, get(result["solution"]["nw"]["1"]["branch"][id], "br_status", b["status"]))
                for (id, b) in data["nw"]["1"]["branch"])
            merge!(stage, Dict("accepted" => accepted, "objective" => result["objective"],
                "result_topology" => topology,
                "solve_seconds" => solve_time(captured[].model),
                "status" => string(termination_status(captured[].model))))
            if refresh_scores && accepted
                stage["score_recovery"] = collect_admittance_scores!(captured[]; topology=topology)
                data["warmstart_time"] += stage["score_recovery"]["score_solve_seconds"]
            end
        end
        #PARAMS[:PRINTING][:BENDERS][:LOG_LINE] = true # print log line for RSCOTS
        #PARAMS[:PRINTING][:BENDERS][:LOG_LINE_HEADER] = true # print log line for RSCOTS

    finally
        reset_switchability!(data["nw"]["1"])
        if stage !== nothing
            for b in values(data["nw"]["1"]["branch"])
                b["heur_switchable"] = b["index"] in stage["selected"]
            end
            delete!(data, "_rscots_stage")
            old_collect === nothing ? delete!(data, "_collect_admittance_scores") : (data["_collect_admittance_scores"] = old_collect)
        end
        if haskey(data, "admittance_warmstart")
            data["admittance_warmstart"]["restored_switchable"] = sort([b["index"] for b in values(data["nw"]["1"]["branch"]) if is_branch_switchable(b)])
        end
        PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP] = TAL_copy
        for key in ("_BestObjStop", "_MIPGap", "_TimeLimit")
            delete!(data, key)
        end
        merge!(data, old_controls)
    end
    return result
end

"Increase TAL in stages, refreshing scores around each accepted solution."
function warmstart_rscots_sequence!(data, solver_convex, slack_types, result; helper=false)
    data["WARMSTART_PROVIDED"] || return result
    target = PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP]
    config = PARAMS[:WARMSTARTING][:RSCOTS]
    mode = config[:STAGE_MODE]
    mode in (:single, :sequential) || error("Unknown restricted stage mode: $mode")
    helper && mode == :single && error("The parallel helper requires sequential stages")
    # The TAL=1 sequential experiment follows the original path exactly.
    target == 1 && return warmstart_rscots!(data, solver_convex, slack_types, result; TAL=1)
    @assert !data["BUSBARS"] && length(data["nw"]) == 1 "Staged RSCOTS currently supports one timestep without busbar reconfiguration"
    pd = PARAMS[:BENDERS][:PARTIAL_DECOMP]
    original_pd = deepcopy(pd)
    original_limit = config[:TIME_LIMIT]
    initial_time = data["warmstart_time"]
    stages = Any[]
    data["rscots_stages"] = stages
    refresh = get(config[:WEIGHTS], :ADMITTANCE_SECANT, 0.0) > 0
    try
        # The accepted model supplies base flows and cut duals for the next stage.
        if refresh
            merge!(pd, Dict(:ENABLED => true, :INCLUDE_BASE_PF => true,
                :N_MAX_CONTINGENCIES => 0, :N_MAX_CONTINGENCIES_PCT => 0.0))
            refresh_derived_params!()
        end
        for tal in (mode == :single ? (target:target) : (1:(helper ? target - 1 : target)))
            if helper
                remaining = original_limit - (data["warmstart_time"] - initial_time)
                remaining <= 0 && break
                config[:TIME_LIMIT] = min(config[:STAGE_TIME_LIMIT], remaining)
            end
            stage = rscots_stage(data, result, tal)
            push!(stages, stage)
            result = warmstart_rscots!(data, solver_convex, slack_types, result;
                TAL=tal, stage=stage, refresh_scores=refresh && tal < target)
        end
        if helper
            data["_rscots_final_stage"] = rscots_stage(data, result, target)
            data["_rscots_final_stage"]["time_limit"] = max(0.0, original_limit - (data["warmstart_time"] - initial_time))
        end
    finally
        config[:TIME_LIMIT] = original_limit
        PARAMS[:BENDERS][:PARTIAL_DECOMP] = original_pd
        refresh_derived_params!()
    end
    return result
end

"perform different types of warmstart, and adjust the parameters accordingly (temporarily)"
function warm_start!(data, solver_convex, slack_types)
    # warm start
    BENDERS_CUT_STRENGTHENING_copy = deepcopy(PARAMS[:BENDERS][:CUT_STRENGTHENING][:ENABLED])
    FRACTIONAL_SEPA_NODELIM_copy = deepcopy(PARAMS[:BENDERS][:CUTS][:FRACTIONAL_SEPA_NODELIM])
    PARTIAL_DECOMP_copy = deepcopy(PARAMS[:BENDERS][:PARTIAL_DECOMP][:ENABLED])
    WARMSTART_METHOD_copy = deepcopy(PARAMS[:WARMSTARTING][:WARMSTART_METHOD])
    PARAMS[:BENDERS][:CUT_STRENGTHENING][:ENABLED] = false # disable cut strengthening for warm starting
    PARAMS[:BENDERS][:CUTS][:FRACTIONAL_SEPA_NODELIM] = -2 # disable fractional separation for warm starting
    PARAMS[:BENDERS][:PARTIAL_DECOMP][:ENABLED] = false
    refresh_derived_params!()

    # ensure no stale state from a previous run leaks into this warmstart
    delete!(data, "_carried_benders_cuts")
    delete!(data, "_benders_cuts_pool")
    delete!(data, "_active_benders_cuts_stage")
    delete!(data, "OLD_WARMSTART_OBJ")
    delete!(data, "WARMSTART_OBJ")
    delete!(data, "admittance_warmstart")
    delete!(data, "rscots_stages")
    delete!(data, "_rscots_final_stage")
    data["WARMSTART_PROVIDED"] = false  # reset, keep key (read by create_master_problem!)

    main_TAL = PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP]
    rscots_TAL = PARAMS[:WARMSTARTING][:RSCOTS][:TAL]

    if PARAMS[:WARMSTARTING][:WARMSTART_METHOD] == :rscots2 && main_TAL < 2
        @info("Warmstart mode RSCOTS2 only sensible for main TAL > 1. Falling back to :rscots.")
        PARAMS[:WARMSTARTING][:WARMSTART_METHOD] = :rscots
    end

    if PARAMS[:WARMSTARTING][:WARMSTART_METHOD] == :opf # warm start with first solving OPF
        result = warmstart_scopf!(data, solver_convex, slack_types, disable_sc=true)
    elseif PARAMS[:WARMSTARTING][:WARMSTART_METHOD] == :ots # warm start with first solving OTS
        result = warmstart_ots!(data, solver_convex, slack_types)
    elseif PARAMS[:WARMSTARTING][:WARMSTART_METHOD] == :scopf # warm start with first solving SCOPF
        result = warmstart_scopf!(data, solver_convex, slack_types)
    elseif PARAMS[:WARMSTARTING][:WARMSTART_METHOD] == :rscots # warm start with first solving SCOTS with only promising switching variables; TAL from settings
        result0 = warmstart_scopf!(data, solver_convex, slack_types)
        result = warmstart_rscots!(data, solver_convex, slack_types, result0, TAL=rscots_TAL)
    elseif PARAMS[:WARMSTARTING][:WARMSTART_METHOD] == :rscots_sequential
        result0 = warmstart_scopf!(data, solver_convex, slack_types)
        helper = get(ENV, "INNOPTEM_EXCHANGE_HELPER", "") == "rscots-sequential" && get(ENV, "INNOPTEM_EXCHANGE_ROLE", "") == "donor"
        result = warmstart_rscots_sequence!(data, solver_convex, slack_types, result0; helper=helper)
    elseif PARAMS[:WARMSTARTING][:WARMSTART_METHOD] == :rscots2 # sequential RSCOTS warmstart at TAL = 1, 2
        result0 = warmstart_scopf!(data, solver_convex, slack_types)
        result1 = warmstart_rscots!(data, solver_convex, slack_types, result0, TAL=1)
        result = warmstart_rscots!(data, solver_convex, slack_types, result1, TAL=2, check_stall=true)
    elseif PARAMS[:WARMSTARTING][:WARMSTART_METHOD] == :rscots3 # sequential RSCOTS warmstart at TAL = 1, 2, 3
        result0 = warmstart_scopf!(data, solver_convex, slack_types)
        result1 = warmstart_rscots!(data, solver_convex, slack_types, result0, TAL=1)
        result2 = warmstart_rscots!(data, solver_convex, slack_types, result1, TAL=2, check_stall=true)
        result = warmstart_rscots!(data, solver_convex, slack_types, result2, TAL=3, check_stall=true)
    else
        error("Unknown WARMSTART_METHOD: $(PARAMS[:WARMSTARTING][:WARMSTART_METHOD])")
    end

    # record number of carried-forward Benders cuts (set by inject_carried_cuts! during the main solve)
    data["n_carried_benders_cuts"] = length(get(data, "_carried_benders_cuts", SerializedCut[]))

    PARAMS[:BENDERS][:CUT_STRENGTHENING][:ENABLED] = BENDERS_CUT_STRENGTHENING_copy
    PARAMS[:BENDERS][:CUTS][:FRACTIONAL_SEPA_NODELIM] = FRACTIONAL_SEPA_NODELIM_copy
    PARAMS[:BENDERS][:PARTIAL_DECOMP][:ENABLED] = PARTIAL_DECOMP_copy
    refresh_derived_params!()
    PARAMS[:WARMSTARTING][:WARMSTART_METHOD] = WARMSTART_METHOD_copy

    if PARAMS[:WARMSTARTING][:RESTRICTED_ONLY]
        data["WARMSTART_PROVIDED"] || error("No feasible restricted solution is available")
        result["primal_status"] == MOI.FEASIBLE_POINT || error("The restricted solution is infeasible")
        diagnostic = get(data, "admittance_warmstart", Dict{String,Any}())
        result["admittance_warmstart"] = diagnostic
        result["rscots_stages"] = get(data, "rscots_stages", Any[])
        result["time"] = Dict("warmstart" => data["warmstart_time"], "solve" => 0.0,
            "restricted" => get(diagnostic, "rscots_solve_seconds", 0.0))
        result["info"] = Dict("status" => get(diagnostic, "rscots_status", string(result["termination_status"])))
        result["unrestricted_solve_executed"] = false
        delete!(data, "OLD_WARMSTART_OBJ")
        throw(RestrictedWarmstartFinished(result))
    end

    # after the restore, so the selection sees the params of the actual run
    prepare_contingency_selection!(data, result=result)

    # clear stale warmstart bookkeeping so repeated runs don't leak state
    delete!(data, "OLD_WARMSTART_OBJ")
end
