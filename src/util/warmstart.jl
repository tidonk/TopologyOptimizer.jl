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
    t = @elapsed begin
        # save original params
        TAL_copy = deepcopy(PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP])
        pd_enabled_copy             = deepcopy(PARAMS[:BENDERS][:PARTIAL_DECOMP][:ENABLED])
        pd_include_base_pf_copy     = deepcopy(PARAMS[:BENDERS][:PARTIAL_DECOMP][:INCLUDE_BASE_PF])
        pd_n_max_copy               = deepcopy(PARAMS[:BENDERS][:PARTIAL_DECOMP][:N_MAX_CONTINGENCIES])
        pd_n_max_pct_copy           = deepcopy(PARAMS[:BENDERS][:PARTIAL_DECOMP][:N_MAX_CONTINGENCIES_PCT])
        #TRANSFORMER_SHIFTING_copy = deepcopy(PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_SHIFTING])
        n_threads_benders_copy = deepcopy(PARAMS[:SOLVER][:MULTITHREADING][:N_THREADS_BENDERS])
        disable_sc && (contingency_copy = deepcopy(data["nw"]["1"]["contingencies"]))

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
        result = solve_secure_ots_time_benders(data, DCPPowerModel, solver_convex, weighted_obj_func_factory([objective_min_num_gen_switch, objective_min_num_tr_switch, objective_min_num_branch_switch, objective_min_redispatch_cost, objective_min_num_redispatch, objective_min_slack_relative_linear, objective_min_num_slacks], [data["cost_coeffs"]["gen_sw"], data["cost_coeffs"]["tr_sw"], data["cost_coeffs"]["br_sw"], data["cost_coeffs"]["redispatch"], data["cost_coeffs"]["num_redispatch"], data["cost_coeffs"]["overload"], data["cost_coeffs"]["num_overload"]]), slack_types=slack_types)
       

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
        #delete!(data, "_NumericFocus")
        disable_sc && (data["nw"]["1"]["contingencies"] = contingency_copy)
    end
    
    data["warmstart_time"] = get(data, "warmstart_time", 0) + t

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
    data["warmstart_time"] = get(data, "warmstart_time", 0) + t 
    return result
end

"run RSCOTS and insert values as warmstart"
function warmstart_rscots!(data, solver_convex, slack_types, result; TAL=PARAMS[:WARMSTARTING][:RSCOTS][:TAL], check_stall=false)
     if ! data["WARMSTART_PROVIDED"] || TAL == 0 || TAL > PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP]
        return result # no warmstart provided means RSCOTS cannot be run for any sensible switching heuristic
    end

    if check_stall && haskey(data, "OLD_WARMSTART_OBJ") && data["WARMSTART_OBJ"] >= data["OLD_WARMSTART_OBJ"]
        @info("Skipping RSCOTS as no improvement in WARMSTART_OBJ")
        return result
    end

    # keep copy of original params
    TAL_copy = deepcopy(PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP])

    # set params
    PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP] = TAL # set TAL for RSCOTS
    data["OLD_WARMSTART_OBJ"] = data["WARMSTART_OBJ"] # keep old warmstart obj
    # _BestObjStop: terminate the warmstart MIP once the incumbent reaches this fraction of the previous warmstart objective.
    data["_BestObjStop"] = PARAMS[:WARMSTARTING][:RSCOTS][:BEST_OBJ_STOP_FRACTION] * data["WARMSTART_OBJ"]
    data["_MIPGap"] = 0.01 # no need for optimal solution
    #data["_MIPFocus"] = 1 # find a good primal solution fast
    #data["_NumericFocus"] = 0 # find a good primal solution fast
    data["_TimeLimit"] = PARAMS[:WARMSTARTING][:RSCOTS][:TIME_LIMIT] # no need for optimal solution

    # solve RSCOTS and set values as warmstart
    calc_and_set_br_sw_heur_scores!(data["nw"]["1"], 
                                    nw_sol=result["solution"]["nw"]["1"], 
                                    params=PARAMS[:WARMSTARTING][:RSCOTS], 
                                    heur_params=PARAMS[:WARMSTARTING][:RSCOTS][:HEUR])
    
    # initialize cut pool for carry-forward if enabled
    if PARAMS[:WARMSTARTING][:RSCOTS][:CARRY_FORWARD_CUTS]
        data["_benders_cuts_pool"] = SerializedCut[]
    end

    result = solve_secure_ots_time_benders(data, DCPPowerModel, solver_convex, weighted_obj_func_factory([objective_min_num_gen_switch, objective_min_num_tr_switch, objective_min_num_branch_switch, objective_min_redispatch_cost, objective_min_num_redispatch, objective_min_slack_relative_linear, objective_min_num_slacks], [data["cost_coeffs"]["gen_sw"], data["cost_coeffs"]["tr_sw"], data["cost_coeffs"]["br_sw"], data["cost_coeffs"]["redispatch"], data["cost_coeffs"]["num_redispatch"], data["cost_coeffs"]["overload"], data["cost_coeffs"]["num_overload"]]), slack_types=slack_types)
    
    if ! PARAMS[:BENDERS][:BASE_PF_IN_MASTER]
        # calculate power flow for base case
        calculate_power_flow!(data, result, solver_convex)
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
    #PARAMS[:PRINTING][:BENDERS][:LOG_LINE] = true # print log line for RSCOTS
    #PARAMS[:PRINTING][:BENDERS][:LOG_LINE_HEADER] = true # print log line for RSCOTS
    
    reset_switchability!(data["nw"]["1"])
    
    # revert to original data
    PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP] = TAL_copy
    delete!(data, "_BestObjStop")
    delete!(data, "_MIPGap")
    #delete!(data, "_MIPFocus")
    #delete!(data, "_NumericFocus")
    delete!(data, "_TimeLimit")
    
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

    # after the restore, so the selection sees the params of the actual run
    prepare_contingency_selection!(data, result=result)

    # clear stale warmstart bookkeeping so repeated runs don't leak state
    delete!(data, "OLD_WARMSTART_OBJ")
end