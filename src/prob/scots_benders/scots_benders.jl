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

import JuMP
import Gurobi
using Printf
using Dualization
using Graphs
using MetaGraphsNext
using Dates
import StatsBase: countmap
using Base.Threads
PAD = " "^82
PAD2 = " "^147

struct WarmstartError <: Exception
    message::String  # A field to hold a descriptive message
end

#-----creation----------------------------------------------------------------------------------------------------------#

""
function solve_secure_ots_time_benders(input::Any, model_type::Type, optimizer::Any, obj_funcs::Union{Function, Tuple{Vector{<:Function}, Vector{<:Real}}}; solver_specific_build_function=ssbf(pm)=pm, multinetwork=true, slack_types=[:pg, :ratea, :dclinepf, :vad], kwargs...)
    return solve_secure_model(input, model_type, optimizer, solver_specific_build_function ∘ build_secure_ots_time_benders, slack_types; ref_extensions=[ref_add_on_off_va_bounds!], multinetwork=true, kwargs...)
end

"create Benders master problem"
function create_master_problem!(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int64=nw_id_default, report::Bool=true)::AbstractPowerModel
    # set up gurobi
    gurobi_master_env = Gurobi.Env()
    standard_params = Dict("LogToConsole"       => PARAMS[:PRINTING][:SOLVER][:LOG_TO_CONSOLE_MP], 
                           "LazyConstraints"   => 1, 
                           "Threads"           => PARAMS[:SOLVER][:MULTITHREADING][:N_THREADS_MIP_MP], 
                           "LogFile"           => get_log_file_name(pm.data),#joinpath(CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:SOLVER_LOG_DIR], "GRB_MP_$(pm.data["name"])-$(replace(pm.data["start_time"], " " => "_")).out"),
                           "TimeLimit"         => PARAMS[:SOLVER][:SOLVING_TIME_LIMIT],
                           "NumericFocus"      => 1,
                           "StartNodeLimit"    => pm.data["WARMSTART_PROVIDED"] ? -1 : -3
                          )
    
    if get(pm.data, "WARMSTART_PROVIDED", false)
        standard_params["Presolve"] = -1
    end

    if haskey(pm.data, "_MIPGap")
        standard_params["MIPGap"] = pm.data["_MIPGap"]
    end
    if haskey(pm.data, "_MIPFocus")
        standard_params["MIPFocus"] = pm.data["_MIPFocus"]
    end
    if haskey(pm.data, "_TimeLimit")
        standard_params["TimeLimit"] = pm.data["_TimeLimit"]
    elseif get(PARAMS[:SOLVER], :INCLUDE_WARMSTART_TIME_LIMIT, false)
        standard_params["TimeLimit"] = max(0.0, standard_params["TimeLimit"] - get(pm.data, "warmstart_time", 0.0))
        pm.data["main_time_limit"] = standard_params["TimeLimit"]
    end
    if haskey(pm.data, "_BestObjStop")
        standard_params["BestObjStop"] = pm.data["_BestObjStop"]
    end
    if haskey(pm.data, "_NumericFocus")
        standard_params["NumericFocus"] = pm.data["_NumericFocus"]
    end
    
    JuMP.set_optimizer(pm.model, optimizer_with_attributes(() -> Gurobi.Optimizer(gurobi_master_env), standard_params..., ))

    # build master problem (contingency PFs come from ref(pm, t, :contingencies); empty when partial decomp disabled or caps = 0)
    build_secure_ots_time(pm, with_base_pf=PARAMS[:BENDERS][:BASE_PF_IN_MASTER])

    # add variable for Benders
    eta_start = haskey(ref(pm, nw), :η) ? ref(pm, nw, :η) : Dict(c => Dict("contingency" => contingency, "η_start" => 0.0) for (c, contingency) in enumerate(contingencies))
    η = var(pm, nw)[:η] = JuMP.@variable(pm.model, 
                                        [c in eachindex(contingencies)], 
                                        base_name = "benders_eta", 
                                        lower_bound = 0, 
                                        start = comp_start_value(eta_start[c], "η_start"))
    
    if get(pm.data, "FIX_BENDERS", false)
        @assert(! get(pm.data, "SUBPROBLEM", false))
        @assert(length(eta_start) == length(contingencies))
        for (c, contingency) in enumerate(contingencies)
            @assert(eta_start[c]["contingency"] == contingency)
            JuMP.fix(η[c], JuMP.start_value(η[c]), force=true)
        end
    end
    report && sol_component_value(pm, nw, :eta, :η, eachindex(contingencies), η)
    
    # define objective function (slack/overload cost terms refer to base-PF slack vars and are only valid when base PF is in master)
    if PARAMS[:BENDERS][:BASE_PF_IN_MASTER]
        weighted_obj_func_factory([objective_min_num_gen_switch, objective_min_num_tr_switch, objective_min_num_branch_switch, 
            objective_min_redispatch_cost, objective_min_num_redispatch, objective_min_slack_relative_linear, objective_min_num_slacks], 
            [pm.data["cost_coeffs"]["gen_sw"], pm.data["cost_coeffs"]["tr_sw"], pm.data["cost_coeffs"]["br_sw"], pm.data["cost_coeffs"]["redispatch"], 
            pm.data["cost_coeffs"]["num_redispatch"], pm.data["cost_coeffs"]["overload"], pm.data["cost_coeffs"]["num_overload"]])(pm)
    else
        weighted_obj_func_factory([objective_min_num_gen_switch, objective_min_num_tr_switch, objective_min_num_branch_switch, 
            objective_min_redispatch_cost, objective_min_num_redispatch], 
            [pm.data["cost_coeffs"]["gen_sw"], pm.data["cost_coeffs"]["tr_sw"], pm.data["cost_coeffs"]["br_sw"], pm.data["cost_coeffs"]["redispatch"], 
            pm.data["cost_coeffs"]["num_redispatch"]])(pm)
    end
    
    # if a custom objective function is provided, set MP to it. Note: Subproblem objectives will remain to be "min overloads"
    if haskey(pm.data, "objective_function")
        if pm.data["objective_function"] == "objective_min_line_load"
            @info("Custom objective function given. Setting Benders subproblem objective function to 'objective_min_line_load'.")
            objective_min_line_load(pm)
        end
    end

    # add Benders variable to objective function
    @assert JuMP.objective_sense(pm.model) == JuMP.MIN_SENSE
    JuMP.@objective(pm.model, Min,
        JuMP.objective_function(pm.model)
        + sum(
            sum( var(pm, n, :η), init=0.0)
        for (n, nw_ref) in nws(pm))
    )
    if PARAMS[:PRECISION][:SOLVER_PREC_DIGITS] >= 0 && PARAMS[:ROUNDING][:MODEL_COEFFS]
        pm.model = round_coeffs!(pm.model)
    end
    #@info("Writing file to mp_$(Dates.format(now(), "HH:MM:SS")).mps and mp_$(Dates.format(now(), "HH:MM:SS")).lp")
    #JuMP.write_to_file(pm.model, "mp_$(Dates.format(now(), "HH:MM:SS")).mps")
    #JuMP.write_to_file(pm.model, "mp_$(Dates.format(now(), "HH:MM:SS")).lp")
    
    # inject carried-forward Benders cuts from previous warmstarting stages
    # Active during RSCOTS stages (pool exists) AND during the main solve (pool was deleted by
    # filter_active_cuts! at the end of the last warmstart stage but _carried_benders_cuts persists).
    if PARAMS[:WARMSTARTING][:RSCOTS][:CARRY_FORWARD_CUTS] && (haskey(pm.data, "_benders_cuts_pool") || !isempty(get(pm.data, "_carried_benders_cuts", [])))
        inject_carried_cuts!(pm, pm.data, contingencies, nw=nw)
    end
    
    return pm
end

"create Benders subproblems"
function create_subproblem!(pm::AbstractPowerModel; nw::Int64=nw_id_default)::AbstractPowerModel
    pm.data["SUBPROBLEM"] = true
    pm_sub = build_ots_slack_benders(pm, n=nw)

    print_cons_with_large_rhs(pm_sub, rhs_lim=1e+7)

    gurobi_sub_env = Gurobi.Env()
    JuMP.set_optimizer(pm_sub.model, optimizer_with_attributes(() -> Gurobi.Optimizer(gurobi_sub_env), 
                    "LogToConsole"      => PARAMS[:PRINTING][:SOLVER][:LOG_TO_CONSOLE_SP], 
                    #"LogFile"           => get_log_file_name(pm.data),#joinpath(CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:SOLVER_LOG_DIR], "GRB_SP_$(pm.data["name"])-$(replace(pm.data["start_time"], " " => "_")).out"),
                    "DualReductions"    => 0,       # required to obtain InfUnbdInfo
                    "InfUnbdInfo"       => 1,       # required to know to add cut
                    "LPWarmStart"       => 0,       # do not use LP warm start
                    "NumericFocus"      => 1,       # some SPs are numerically challenging
                    "Threads"           => PARAMS[:SOLVER][:MULTITHREADING][:N_THREADS_MIP_SP]))
                    
    _weighted_obj_func_factory(([objective_min_slack_relative_linear], [pm.data["cost_coeffs"]["overload"]]))(pm_sub)
    if PARAMS[:PRECISION][:SOLVER_PREC_DIGITS] >= 0 && PARAMS[:ROUNDING][:MODEL_COEFFS]
        pm_sub.model = round_coeffs!(pm_sub.model)
    end
    return pm_sub
end

"remove z_branch variable from cut"
function remove_var_from_expr!(lhs::JuMP.AffExpr, master_vars::Dict{Symbol, AbstractVector}, outage_branches::Vector{Int64})
    # remove outage branches from cut
    for b_id in outage_branches
        # remove the switching variable from the cut
        if PARAMS[:LIMITS][:BRANCH][:BRANCH_SWITCHING]
            if b_id in axes(master_vars[:z_branch], 1)
                z_branch = master_vars[:z_branch][b_id]
                if haskey(lhs.terms, z_branch)
                    lhs.terms[z_branch] = 0.0
                end
            end
        end

        # remove the outaged branches' transformer shifting variable from the cut
        if PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_SHIFTING] && :tr_shift ∈ keys(master_vars)
            if b_id in axes(master_vars[:tr_shift], 1)
                tr_shift = master_vars[:tr_shift][b_id]
                if haskey(lhs.terms, tr_shift)
                    lhs.terms[tr_shift] = 0.0
                end
            end
        end
    end
end

"create Benders (reduced-size) subproblem linear programs"
function create_lp!(pm::AbstractPowerModel; nw::Int64=nw_id_default)::AbstractPowerModel
    pm_lp = build_ots_benders(pm, relax=false, nw=nw)
    
    gurobi_sub_env = Gurobi.Env()
    JuMP.set_optimizer(pm_lp.model, optimizer_with_attributes(() -> Gurobi.Optimizer(gurobi_sub_env), 
                    "LogToConsole"      => PARAMS[:PRINTING][:SOLVER][:LOG_TO_CONSOLE_SP], 
                    "DualReductions"    => 0,       # required to obtain InfUnbdInfo
                    "InfUnbdInfo"       => 1,       # required to know to add cut
                    "LPWarmStart"       => 0,       # do not use LP warm start
                    "Threads"           => PARAMS[:SOLVER][:MULTITHREADING][:N_THREADS_MIP_SP]))

    if PARAMS[:PRECISION][:SOLVER_PREC_DIGITS] >= 0 && PARAMS[:ROUNDING][:MODEL_COEFFS]
        pm_lp.model = round_coeffs(pm_lp.model)
    end
    return pm_lp
end

# -----helper functions-------------------------------------------------------------------------------------------------------#
"initialize the callback stats tracking dictionary"
function init_callback_stats(;nested::Bool=false)::Dict{Symbol,Int64}
    base_stats = [:n_callbacks, 
                  :n_solved_subs_total, 
                  :n_benders_cuts_added, 
                  :n_benders_opt_cuts_found, 
                  :n_benders_feas_cuts_found, 
                  :n_benders_cuts_found, 
                  :n_connect_cuts_added, 
                  :n_leq_fails, 
                  :n_leq_residual_failures,
                  :n_questionable_sol, 
                  :n_primal_lp,
                  :n_fallback_primal_lp,
                  :n_dual_lp, 
                  :n_dual_solver, 
                  :n_cut_value_fail, 
                  :n_cut_value_disagree, 
                  :n_cut_value_agree,
                  :n_rf_new,
                  :n_rf_update,
                  :n_frac_sepas
                 ]

    if nested
        stats = [base_stats; :n_benders_cuts_added_mp_gens]
    else
        stats = base_stats
    end

    return Dict{Symbol,Int64}(s => 0 for s in stats)
end

# Subproblem workers can enter the LP path concurrently.
const _LP_CALL_STATS_LOCK = ReentrantLock()

"initialize the callback times tracking dictionary"
function init_callback_times(;nested::Bool=false)::Dict{Symbol,Float64}
    base_times = [:pre, 
                  :connect, 
                  :contingency, 
                  :contingency_primal, 
                  :contingency_primal_check, 
                  :contingency_dual, 
                  :contingency_dual_check,
                  :contingency_fallback,  
                  :contingency_cut_build, 
                  :contingency_cut_check, 
                  :contingency_cut_opt, 
                  :contingency_cut_feas, 
                  :contingency_cut_none,
                  :refactorization
                 ]

    if nested
        times = base_times
    else
        times = base_times
    end

    return Dict{Symbol,Float64}(s => 0.0 for s in times)
end

"plug in each master variable value into the respective variable in the expression"
function plug_in_values(lhs::JuMP.AffExpr, master_var_vals::Dict{Symbol, AbstractVector})
    # for each variable, plug in the value of the variable
    val_evaluated = lhs.constant
    for (var, coeff) in lhs.terms
        # plug in the value of the variable
        if coeff != 0.0
            var_name, var_idx = get_var_basename_idx(var)

            val_evaluated += coeff * master_var_vals[var_name][var_idx]::Float64
        end
    end
    return val_evaluated
end
function plug_in_values(lhs::JuMP.AffExpr, master_var_vals_var::Dict{VariableRef, Float64})
    # for each variable, plug in the value of the variable
    val_evaluated = lhs.constant
    for (var, coeff) in lhs.terms
        # plug in the value of the variable
        if coeff != 0.0
            val_evaluated += coeff * master_var_vals_var[var]
        end
    end
    return val_evaluated
end

"plug in the values of the master variables into the cut and compare against new_η_val. Return evaluated LHS value and whether they agree"
function check_lhs_value(lhs::JuMP.AffExpr, cb_data::Gurobi.CallbackData, val::Float64, outage_branches::Vector{Int64}; is_opt_cut::Bool=true)::Tuple{Float64, Bool}
    val_evaluated = JuMP.value(x ->JuMP.callback_value.(cb_data, x), lhs)

    return check_lhs_value(val_evaluated, val, outage_branches, is_opt_cut=is_opt_cut)
end
function check_lhs_value(lhs::JuMP.AffExpr, master_var_vals::Dict{Symbol, AbstractVector}, val::Float64, outage_branches::Vector{Int64}; is_opt_cut::Bool=true)::Tuple{Float64, Bool}
    val_evaluated = plug_in_values(lhs, master_var_vals)

    return check_lhs_value(val_evaluated, val, outage_branches, is_opt_cut=is_opt_cut)
end
function check_lhs_value(lhs::JuMP.AffExpr, master_var_vals::Dict{VariableRef, Float64}, val::Float64, outage_branches::Vector{Int64}; is_opt_cut::Bool=true)::Tuple{Float64, Bool}
    val_evaluated = plug_in_values(lhs, master_var_vals)

    return check_lhs_value(val_evaluated, val, outage_branches, is_opt_cut=is_opt_cut)
end
function check_lhs_value(val_evaluated::Float64, val::Float64, outage_branches::Vector{Int64}; is_opt_cut::Bool=true)::Tuple{Float64, Bool}
    if is_opt_cut
        # check whether val evaluated is the same as the val
        if ! isapprox(val_evaluated, val, atol=PARAMS[:BENDERS][:CHECKS][:CHECK_CUT_VALUE_TOL])
            str = @sprintf("Cut LHS = %12.4f; SP obj. = %12.4f; relative diff = %6.2f %%", val_evaluated, val, abs(val_evaluated - val) / max(abs(val_evaluated), abs(val)) * 100) 
            @debug("[Check cut val] LHS value is different from η value (relative difference. $str. outage branches = $outage_branches")
            callback_stats[:n_cut_value_disagree] += 1

            return val_evaluated, false
        end
    else
        # check whether val_evaluated is > 0
        if val_evaluated ≤ PARAMS[:PRECISION][:SOLVER_FEASTOL]
            callback_stats[:n_cut_value_disagree] += 1

            return val_evaluated, false
        end
    end

    callback_stats[:n_cut_value_agree] += 1
    return val_evaluated, true
end

"for a given (solved) model, check whether the reduced cost of binary variables is zero"
function check_redcost(model)
    rc_sum = 0.0
    for var in JuMP.all_variables(model)
        if JuMP.reduced_cost(var) != 0.0 && JuMP.is_binary(var)
            rc_sum += JuMP.reduced_cost(var) * JuMP.value(var)
            #@show(var, JuMP.reduced_cost(var), JuMP.value(var))
        end
    end
    if rc_sum != 0.0
        @warn("Reduced cost of binary variable is not zero: $(rc_sum)")
    end
end

#-----solving----------------------------------------------------------------------------------------------------------#
"adds the time on the clock to 'time' and updates the clock to the current time"
function count_time!(clocks::Vector{Float64}, callback_times, type::Symbol, counter_fallback::Bool)
    # track fallback time separately
    if counter_fallback && type ∈ [:contingency_primal, :contingency_dual, :contingency_primal_check, :contingency_dual_check, :contingency_cut_build, :contingency_cut_check]
        type = :contingency_fallback
    end

    if type ∈ [:contingency_primal, :contingency_dual, :contingency_primal_check, :contingency_dual_check, :contingency_cut_build, :contingency_cut_check]
        callback_times[type] += time() - clocks[1]
        clocks[1] = time()
    elseif type ∈ [:contingency_cut_opt, :contingency_cut_feas, :contingency_cut_none]
        callback_times[type] += time() - clocks[2]
        clocks[2] = time()
    elseif type == :contingency_fallback
        callback_times[type] += time() - clocks[1]
        clocks[1] = time()
    else
        @error("Unknown type $(type) for callback time")
    end
end

"solve a single Benders subproblem for a given contingency"
function solve_subproblem(pm_sub::AbstractPowerModel,
                          c::Int64, 
                          cb_data,
                          master_varss::Vector{Dict{Symbol, AbstractVector}},                  # one per master (nested)
                          master_var_vals::Dict{Symbol, AbstractVector},
                          master_var_vals_var::Dict{VariableRef, Float64},                     # same values as master_var_vals
                          lhsss::Vector{Dict{Symbol, Dict{Int, JuMP.GenericAffExpr}}},         # one per master (nested)
                          dualss::Dict{Symbol, Dict{Int, Float64}},                            # allocated with dummy values
                          outage_branches::Vector{Int64},
                          solver::AbstractSubproblemSolver,
                          cut_exprs::Vector{AffExpr};                                          # (partially) allocated for `cut_expr ≤ η`, one per master (nested)
                          nw::Int64=nw_id_default,
                          int_sol::Bool=true,
                          hint_feascut::Bool=false,
                          counter_fallback::Bool=false,
                          prec_mult=1.0
    )::Tuple{<:Union{Vector{JuMP.ScalarConstraint}, Nothing}, Float64, Symbol, Type}

    clocks = [time(), time()] # timestamps to track time, categorized by functionality (e.g., cut check etc.) and generated cut type

    ### Solve Primal Problem ###
    # Count actual primal LP attempts, including those that produce no cut.
    # Fallback recursion passes the underlying LP solver, not FallbackSolver.
    if !hint_feascut && solver isa Union{FullLPSolver, LPSolver}
        lock(_LP_CALL_STATS_LOCK) do
            callback_stats[:n_primal_lp] += 1
            callback_stats[:n_fallback_primal_lp] += Int(counter_fallback)
        end
    end
    success_primal = ! hint_feascut && solve_primal!(solver)
    if success_primal && solver isa AbstractMatrixSolver
        success_primal = align_component_angles!(solver, master_var_vals, outage_branches)
    end
    count_time!(clocks, callback_times, :contingency_primal, counter_fallback)
    ### -- ###

    ### Check Primal Solution & Calculate Dual Solution ###
    if success_primal
        # examine the given primal solution (saved in `solver`)

        ### Check Primal Solution ###
        has_new_cut, is_opt_cut, obj_func_value, betas = check_primal_solution(solver, pm_sub, c, master_var_vals, outage_branches, nw=nw)
        
        if PARAMS[:BENDERS][:CHECKS][:CHECK_WITH_FALLBACK] && isa(solver, AbstractMatrixSolver)
            fallback_start = time()
            (_, eta, type, solver_type) = __check_primal_solution_with_fallback(solver, pm_sub, c, cb_data, master_var_vals, master_var_vals_var, outage_branches, master_varss[1], lhsss[1], dualss, cut_exprs, nw=nw)
            clocks[1] += time() - fallback_start
            type_orig = has_new_cut ? (is_opt_cut ? :Optimality : :Feasibility) : :Nothing

            if !isapprox(eta, obj_func_value, atol=0.1) || type_orig != type
                solver_orig = PARAMS[:BENDERS][:SUBPROBLEM_SOLVER_TYPE]

                @printf("%s ├ %4s │ LEQ Solution check: Disagreement between subproblem solvers! Contingency: %s\n", PAD2, "FAIL", pm_sub.data["nw"]["1"]["contingencies"][c][:branch])
                @printf("%s ├ %4s │ %1s │ %9.2f\n", PAD2, CONSTS[:INTERNAL_MAPPINGS][:SLACK_VARIABLES][:SOLVER_NAME_MAP][Symbol(nameof(solver_orig))], string(type_orig)[1], obj_func_value)
                @printf("%s └ %4s │ %1s │ %9.2f\n", PAD2, CONSTS[:INTERNAL_MAPPINGS][:SLACK_VARIABLES][:SOLVER_NAME_MAP][Symbol(nameof(solver_type))], string(type)[1], eta)

                # for debugging - infeasible subproblems might hint at data errors
                if JuMP.termination_status(pm_sub.model) == JuMP.MOI.INFEASIBLE && type_orig != :Feasibility
                    @info("Optimizing IIS of Subproblem:")
                    JuMP.compute_conflict!(pm_sub.model)
                    iis_model, _ = JuMP.copy_conflict(pm_sub.model)
                    print(iis_model)
                    # from iis_model, after "Subject to", regex b_ids out of all p[(b_id, ...)]
                    matches = parse.(Int64, unique!(collect(m.captures[1] for m in eachmatch(r"p\[\((\d+),", string(iis_model)))))
                    
                    # for the LEQ solution, print power flows and angles of all branches contained in the IIS
                    for (branch_id::Int64, branch::Dict{Symbol, Any}) in solver.branch_data
                        f_bus = branch[:f_bus]
                        t_bus = branch[:t_bus]
                        vad1 = get_vad(solver, f_bus, t_bus)
                        pf1 = get_pf(solver, branch, vad1)
                        if branch_id ∈ matches
                            @printf("Branch %4d: pf  = %10.6f, rate_a = %10.6f, overload factor = %4.2f\n", branch_id, pf1, branch[:rate_a], pf1/branch[:rate_a])
                            @printf("Branch %4d: vad = %10.6f\n", branch_id, vad1)
                        end
                    end
                    throw(InvalidStateException("Fallback solver says infeasible, but LEQ solution is feasible", :SolutionCheckDisagreement))
                end

                # check cost components during their collection
                if PARAMS[:BENDERS][:CHECKS][:CHECK_COST_COMPONENTS] && isa(solver, AbstractMatrixSolver)
                    check_primal_solution(solver, pm_sub, c, master_var_vals, outage_branches, nw=nw, print_costs=true)
                end
            end
        end

        count_time!(clocks, callback_times, :contingency_primal_check, counter_fallback)
        ### -- ###

        ### Calculate Dual Solution ###
        if ! has_new_cut
            if typeof(solver) in [FullLPSolver, LPSolver] && JuMP.termination_status(solver.pm.model) != JuMP.MOI.OPTIMAL
                return nothing, 0, :NumericalError, FullLPSolver
            end
            return nothing, 0, :Nothing, typeof(solver)
        end
            
        ## we need to add a cut, so calculate the duals
        # get dual values
        # solve opt/feas dual problem, set `success_dual` depending on whether the solver was able to handle the problem
        if typeof(solver) == FullLPSolver || (PARAMS[:BENDERS][:SUBPROBLEM_SOLVER_DUAL_SOLVE] && !(typeof(solver) <: AbstractMatrixSolver && !is_opt_cut))
            # compute dual var values
            update_dual_opt!(solver, pm_sub, master_var_vals, betas, nw=nw) # update dual problem with primal solution
            success_dual = solve_dual_opt!(solver) && check_dual_solution(solver)

        # feasibility cut
        elseif PARAMS[:BENDERS][:SUBPROBLEM_SOLVER_DUAL_SOLVE] && typeof(solver) <: AbstractMatrixSolver && !is_opt_cut
            update_dual_feas!(solver, pm_sub, master_var_vals, outage_branches, betas, nw=nw) # update dual problem with primal solution
            success_dual = solve_dual_feas!(solver) && check_dual_solution(solver) &&
                           check_feasibility_ray(solver, master_var_vals, outage_branches, betas)
        else
            success_dual = false
        end
        
        callback_stats[:n_dual_solver]    += Int(success_dual)
        callback_stats[:n_dual_lp]        += Int(success_dual && typeof(solver) == FullLPSolver)

        count_time!(clocks, callback_times, :contingency_dual, counter_fallback)
        ### -- ###
    else
        success_dual = false
    end
    ### -- ###

    ### Get Dual Solution ###
    # compute dual values for optimality/feasibility cut
    # if we have calculated a primal and dual solution (success_dual == true), fill dualss with the necessary dual values
    if success_dual
        # reset dualss (since most dual values are 0, it is simpler to only set the non-zeroes in update_dualss_*!)
        for con_symbol in keys(dualss)
            for comp_id in keys(dualss[con_symbol])
                dualss[con_symbol][comp_id] = 0.0
            end
        end

        # compute dual values for optimality/feasibility cut
        if is_opt_cut
            update_dualss_opt!(dualss, pm_sub, master_var_vals, outage_branches, solver, betas, nw=nw)
        else
            update_dualss_feas!(dualss, pm_sub, master_var_vals, outage_branches, solver, betas, nw=nw)
        end

    # if we couldn't solve using `solver` for some reason, solve using its fallback solver, unless it doesn't have one
    elseif hasproperty(solver, :fallback_solver)
        # apply the buffered update, then solve, then reset the buffered update, then return the computed cut
        old_rhss = apply_fallback_update!(solver.fallback_solver)
        fallback_start = time()
        (new_con, eta, type, solver_type) = solve_subproblem(pm_sub, c, cb_data, master_varss, master_var_vals, master_var_vals_var, lhsss, dualss, outage_branches, solver.fallback_solver.solver, cut_exprs, nw=nw, int_sol=int_sol, counter_fallback=true)
        clocks[1] += time() - fallback_start
        reset_fallback_update!(solver.fallback_solver, old_rhss)
        count_time!(clocks, callback_times, :contingency_dual_check, counter_fallback)
        return (new_con, eta, type, solver_type)

    # this really shouldn't happen, unless perhaps a FullLPSolver is not at the end of a fallback chain
    else
        @warn("Solver $(typeof(solver)) failed, and there is no fallback solver remaining. Unless another cut is found in this iteration, the solution might no longer be optimal.")
        count_time!(clocks, callback_times, :contingency_dual_check, counter_fallback)
        return nothing, 0, :Nothing, typeof(solver)
    end

    count_time!(clocks, callback_times, :contingency_dual_check, counter_fallback)
    ### -- ###

    ### Build Cut(s) ###
    cuts = Vector{JuMP.ScalarConstraint}() # one for each version of master vars provided by lhsss

    constant = is_opt_cut ? JuMP.objective_function(pm_sub.model).constant : 0.0
    
    new_η_val = obj_func_value
    is_cut_val_correct = true

    # go over lhsss from different levels of nested Benders, if applicable 
    for (i, lhss) in enumerate(lhsss)
        master_vars = master_varss[i]
        cut_expr = cut_exprs[i]
        
        cut_expr = update_benders_expression!(cut_expr, constant, lhss, dualss, master_vars, outage_branches)
        
        # check if a cut should be added
        if (success_dual && is_opt_cut) || (typeof(solver) == FullLPSolver && JuMP.termination_status(pm_sub.model) == JuMP.MOI.OPTIMAL)
            # get the master η variable object for the cut
            η_val = master_var_vals[:η][c]
            
            if i == 1 && PARAMS[:BENDERS][:CHECKS][:CHECK_CUT_VALUE] && hasproperty(solver, :fallback_solver) # i == 1 means we're in the usual cuts, not those from the intermediate problem
                # check if expression LHS value is same as obj_func_value
                new_η_val, is_cut_val_correct = check_lhs_value(cut_expr, master_var_vals_var, obj_func_value, outage_branches)
            end

            if _isapproxless(η_val, new_η_val, atol=prec_mult*SOLVER_PRECISION)
                η = master_vars[:η][c]
                type = :Optimality

                callback_stats[:n_benders_opt_cuts_found] += 1
                count_time!(clocks, callback_times, :contingency_cut_opt, counter_fallback)
            else
                # check inconsistent cuts with the LP before returning a miss
                if !is_cut_val_correct && hasproperty(solver, :fallback_solver)
                    old_rhss = apply_fallback_update!(solver.fallback_solver)
                    try
                        fallback_result = solve_subproblem(pm_sub, c, cb_data,
                            master_varss, master_var_vals, master_var_vals_var,
                            lhsss, dualss, outage_branches, solver.fallback_solver.solver,
                            cut_exprs; nw=nw, int_sol=int_sol, counter_fallback=true,
                            prec_mult=prec_mult)
                        callback_stats[:n_cut_value_fail] += 1
                        return fallback_result
                    finally
                        reset_fallback_update!(solver.fallback_solver, old_rhss)
                    end
                end
                if ! PARAMS[:BENDERS][:CHECKS][:CHECK_WITH_FALLBACK] && PARAMS[:BENDERS][:SUBPROBLEM_SOLVER_TYPE] != typeof(solver) && new_η_val - η_val > SOLVER_PRECISION
                    @warn("The last fallback solver in the chain of type $(typeof(solver)) did not find any cuts for contingency $c, although the previous solver $(PARAMS[:BENDERS][:SUBPROBLEM_SOLVER_TYPE]) detected that there might have to be a cut with violation $(new_η_val) - $(η_val) = $(new_η_val - η_val). This might hint at a poorly scaled cut.")
                end

                count_time!(clocks, callback_times, :contingency_cut_none, counter_fallback)
                return nothing, new_η_val, :Nothing, typeof(solver)
            end
        # feasibility cut
        elseif (success_dual && !is_opt_cut) || (typeof(solver) == FullLPSolver && JuMP.termination_status(pm_sub.model) == JuMP.MOI.INFEASIBLE)
            η = JuMP.AffExpr(0)
            type = :Feasibility

            if i == 1 # i == 1 means we're in the usual cuts, not those from the intermediate problem
                if PARAMS[:BENDERS][:CHECKS][:CHECK_FEAS_CUT_VALUE]
                    new_η_val, is_cut_val_correct = check_lhs_value(cut_expr, master_var_vals_var, PARAMS[:PRECISION][:SOLVER_FEASTOL], outage_branches, is_opt_cut=false)
                else
                    new_η_val = Inf
                end
            end

            callback_stats[:n_benders_feas_cuts_found] += 1
            count_time!(clocks, callback_times, :contingency_cut_feas, counter_fallback)
        else
            @error("Unhandled termination status: $(JuMP.termination_status(pm_sub.model)) (primal status: $(JuMP.primal_status(pm_sub.model)), dual status: $(JuMP.dual_status(pm_sub.model))")
        end

        if ! is_cut_val_correct && PARAMS[:BENDERS][:CHECKS][:CHECK_COST_COMPONENTS] && isa(solver, AbstractMatrixSolver)
            check_primal_solution(solver, pm_sub, c, master_var_vals, outage_branches, nw=nw, print_costs=true)
        end

        cut = JuMP.@build_constraint(cut_expr <= η)

        count_time!(clocks, callback_times, :contingency_cut_build, counter_fallback)
        ### -- ###

        ### Check Cut ###
        # if the SP obj was not the same as the LHS val, something might be wrong with the cut, so run the fallback solver on it to be sure
        if !is_cut_val_correct && (!is_opt_cut || PARAMS[:BENDERS][:CHECKS][:CHECK_FALLBACK_ON_CUT_VAL_DISAG])
            if hasproperty(solver, :fallback_solver)
                @info("Checking with fallback solver, because cut value = $(round(new_η_val, digits=3)) ≠ $(round(obj_func_value, digits=3)) = LHS value and type is $type.")
                old_rhss = apply_fallback_update!(solver.fallback_solver)
                fallback_start = time()
                local new_con_fb, eta_fb, type_fb, solver_type
                try
                    (new_con_fb, eta_fb, type_fb, solver_type) = solve_subproblem(pm_sub, c, cb_data, master_varss, master_var_vals, master_var_vals_var, lhsss, dualss, outage_branches, solver.fallback_solver.solver, cut_exprs, nw=nw, int_sol=int_sol, counter_fallback=true, prec_mult=prec_mult)
                    # Read the LP result before resetting its constraints invalidates it.
                    fb_solver = solver.fallback_solver.solver
                    if typeof(fb_solver) in [FullLPSolver, LPSolver]
                        fb_model = fb_solver.pm.model
                        @info "[FB on cut val disagree] LP result" contingency=c outage_branches=outage_branches status=JuMP.termination_status(fb_model) primal_status=JuMP.primal_status(fb_model) lp_objective=(JuMP.has_values(fb_model) ? JuMP.objective_value(fb_model) : nothing) master_eta=master_var_vals[:η][c] leq_objective=obj_func_value leq_cut_value=new_η_val fallback_type=type_fb
                    end
                finally
                    clocks[1] += time() - fallback_start
                    reset_fallback_update!(solver.fallback_solver, old_rhss)
                end

                # Infeasibility agreement does not validate the rejected ray.
                # Return all fallback cuts, including nested levels and solver provenance.
                callback_stats[:n_cut_value_fail] += 1
                count_time!(clocks, callback_times, :contingency_cut_check, counter_fallback)
                return new_con_fb, eta_fb, type_fb, solver_type
            else
                @warn("Last fallback solver $(typeof(solver)) in the solver chain added a cut that it doesn't agree with itself. This could hint at a poorly scaled cut. \n lhs val = $(round(new_η_val, digits=3)), SP obj func = $(round(obj_func_value, digits=3)), difference: $(round(new_η_val - obj_func_value, digits=3))")
            end
        end
        
        push!(cuts, cut)
        count_time!(clocks, callback_times, :contingency_cut_check, counter_fallback)
    end
    ### -- ###

    return cuts, new_η_val, type, typeof(solver)
end


"for all monitored branches, apply the adjusted rate a in the constraint `constraint_thermal_limit_from{_on_off}"
function apply_adjusted_rate_a!(pm_sub::AbstractPowerModel, 
                                sub_cons_for_duals::Dict{Symbol, Dict{Int, JuMP.ConstraintRef}},
                                lhsss::Vector{Dict{Symbol, Dict{Int, JuMP.GenericAffExpr}}};
                                for_n1::Bool=true,
                                nw::Int64=nw_id_default)
    for (b_id, branch) in ref(pm_sub, nw, :branch)
        if ! is_branch_electrical(branch)
            continue
        end

        adjusted_rate_a_n1 = get_adjusted_rate_a(branch, for_n1)

        dif = get_adjusted_rate_a(branch, !for_n1) - adjusted_rate_a_n1

        if is_branch_monitored(branch)
            index = branch["index"]::Int64

            JuMP.set_normalized_rhs(sub_cons_for_duals[:rating_from_slack_lb][index], -adjusted_rate_a_n1)
            JuMP.set_normalized_rhs(sub_cons_for_duals[:rating_from_slack_ub][index],  adjusted_rate_a_n1)

            for lhss in lhsss
                lhss[:rating_from_slack_lb][index].constant += dif
                lhss[:rating_from_slack_ub][index].constant -= dif
            end
        end

        #if is_branch_switchable(branch) # this shouldn't be necessary, since it is a big M constraint, only forcing 0/any, and no limit
        #end
    end
end

"get the right-hand side of a constraint"
function get_con_rhs(con::JuMP.ScalarConstraint)
    return con.set.upper
end

"Solve Benders subproblems in parallel. Each thread executes this function and continually tries to grab the next subproblem."
function solve_subproblem_multi(thread_id::Int64, 
                                pm_sub::AbstractPowerModel,
                                cb_data,
                                subproblem_prio_list::Vector{Dict{Symbol, Any}},
                                contingencies::Vector{Contingency}, 
                                sub_cons_for_duals::Dict{Symbol, Dict{Int, JuMP.ConstraintRef}},
                                master_varss::Vector{Dict{Symbol, AbstractVector}},              # one per master (nested)
                                master_var_vals::Dict{Symbol, AbstractVector},
                                master_var_vals_var::Dict{VariableRef, Float64},                 # same values as master_var_vals
                                lhsss::Vector{Dict{Symbol, Dict{Int, JuMP.GenericAffExpr}}},     # one per master (nested)
                                dualss::Dict{Symbol, Dict{Int, Float64}},                        # allocated with dummy values
                                cur_pos::Vector{Int64},
                                cur_pos_sem::Base.Semaphore, 
                                counters::Dict{Symbol, Real},
                                counters_sem::Base.Semaphore,
                                solver::AbstractSubproblemSolver,
                                cut_exprs::Vector{AffExpr};                                      # (partially) allocated for `cut_expr ≤ η`, one per master (nested)
                                nw::Int64=nw_id_default,
                                _post_func::Function=_prio_post_set_call,
                                int_sol::Bool=true,
                                prec_mult=1.0
                               )::Tuple{Vector{Vector{JuMP.ScalarConstraint}}, Vector{Float64}, Vector{Int64}, Vector{Symbol}}
    
    cons = Vector{Vector{JuMP.ScalarConstraint}}()
    violations = Vector{Float64}() # how much more expensive is the solution than estimated by the master?
    subproblem_ids = Vector{Int64}()
    cut_types = Vector{Symbol}()
    start_time = time()

    apply_master!(solver, pm_sub, master_var_vals, nw=nw)
    while true
        # get next subproblem to work on
        Base.acquire(cur_pos_sem)
        pos         = cur_pos[1]
        cur_pos[1] += 1
        Base.release(cur_pos_sem)
        
        # get contingency and outage branches
        if pos > length(subproblem_prio_list)
            break
        end
        prio            = subproblem_prio_list[pos][:priority]
        c_id            = subproblem_prio_list[pos][:contingency_id]
        hint_feascut    = get(subproblem_prio_list[pos], :hint_feascut, false)
        outage_branches = get_outage_branches(pm_sub, contingencies[c_id], master_var_vals, filter=true, nw=nw)
        # Already-open failed lines are also exempt from off-angle bounds, but
        # must not be removed from the PF matrix or switched back on a second time.
        all_outage_branches = get_outage_branches(pm_sub, contingencies[c_id], master_var_vals, nw=nw)

        # if no outage branch, apply adjusted rate a for master problem
        if isempty(all_outage_branches)
            pm_sub.data["SUBPROBLEM"] = false
            apply_adjusted_rate_a!(pm_sub, sub_cons_for_duals, lhsss, for_n1=false, nw=nw)
        end
        # turn off z variable of outage branches
        switch_off_outage_branches!(master_var_vals, outage_branches)
        # apply additional, solver specific, modifications
        solver_outages = solver isa AbstractMatrixSolver ? outage_branches : all_outage_branches
        old_rhss = apply_contingency!(solver, pm_sub, solver_outages, master_var_vals, int_sol, nw=nw)
        if solver isa AbstractMatrixSolver
            apply_contingency!(solver.fallback_solver, pm_sub, all_outage_branches, master_var_vals, int_sol, nw=nw)
        end

        # solve subproblem
        t = @elapsed begin 
            (new_con, eta, type, solver_type) = solve_subproblem(pm_sub, 
                                                                 c_id, 
                                                                 cb_data, 
                                                                 master_varss, 
                                                                 master_var_vals,
                                                                 master_var_vals_var,
                                                                 lhsss, 
                                                                 dualss, 
                                                                 all_outage_branches, # these vars are excluded from the cut + used to check primal and dual feasibility
                                                                 solver,
                                                                 cut_exprs,
                                                                 nw=nw, 
                                                                 int_sol=int_sol,
                                                                 hint_feascut=hint_feascut,
                                                                 prec_mult=prec_mult)
        end

        if haskey(pm_sub.ext, :admittance_cuts) && type in (:Optimality, :Feasibility)
            record_admittance_cut!(pm_sub, c_id, first(new_con), dualss)
        end

        # reset additional modifications
        reset_contingency!(solver, pm_sub, solver_outages, old_rhss, nw=nw)
        # turn contingency branches back on if they were on before
        switch_on_outage_branches!(master_var_vals, outage_branches)
        # if no outage branch, apply adjusted rate a for master problem
        if isempty(all_outage_branches)
            pm_sub.data["SUBPROBLEM"] = true
            apply_adjusted_rate_a!(pm_sub, sub_cons_for_duals, lhsss, for_n1=true, nw=nw)
        end

        # a cut was generated
        violation = 0.0
        if type == :Feasibility && ! PARAMS[:DEBUG][:ENFORCE_MAX_ANG_DIFF]  
            # skip this cut
            count_sym = :n_unrobust_hits
            PARAMS[:PRINTING][:BENDERS][:LOG_LINE] && @printf("%s%5d | %5d | %3d | %+.2f | %.2f | %*s | %9.2f | %9.2f | %4s\n", PAD, c_id, pos, thread_id, round(prio, digits=3), round(t, digits=2), ndigits(length(contingencies)), string(type)[1], eta, master_var_vals[:η][c_id], CONSTS[:INTERNAL_MAPPINGS][:SLACK_VARIABLES][:SOLVER_NAME_MAP][Symbol(nameof(solver_type))])
        elseif type ∈ [:Feasibility, :Optimality]
            # check if the cut is unrobust, i.e., has a very small violation
            rhs = get_con_rhs(new_con[1])

            if type == :Optimality && eta - master_var_vals[:η][c_id] < PARAMS[:BENDERS][:PARTIAL][:ETA_ROBUSTNESS]
                count_sym = :n_unrobust_hits
            # here, eta is the LHS value, so check if lhs/rhs > PARAMS[:FEASTOL]
            elseif type == :Feasibility && eta <= PARAMS[:BENDERS][:PARTIAL][:ETA_ROBUSTNESS_FEAS]
                #@info("Unrobust feasibility cut due to eta = $(round(eta, digits=8)) < $(PARAMS[:BENDERS][:PARTIAL][:ETA_ROBUSTNESS_FEAS]) = eta robustness feas")
                count_sym = :n_unrobust_hits
            elseif type == :Feasibility && eta/rhs <= PARAMS[:PRECISION][:SOLVER_FEASTOL]
                #@info("Unrobust feasibility cut due to (lhs+rhs)/rhs = $(round(eta/rhs, digits=8)) < $(PARAMS[:PRECISION][:SOLVER_FEASTOL]) = solver feasibility tolerance where lhs = $(round(eta, digits=8)) and rhs = $(round(rhs, digits=8))")
                count_sym = :n_unrobust_hits
            else
                count_sym = :n_hits
                subproblem_prio_list[pos][:n_cuts] += 1
            end
            push!(cons, new_con)
            if type == :Optimality
                violation = eta - master_var_vals[:η][c_id]
            else # Feasibility
                violation = eta # workaround: eta here means the LHS value
            end
            push!(violations, violation)
            push!(subproblem_ids, c_id)
            push!(cut_types, type)
            PARAMS[:PRINTING][:BENDERS][:LOG_LINE] && @printf("%s%5d | %5d | %3d | %+.2f | %.2f | %*s | %9.2f | %9.2f | %4s\n", PAD, c_id, pos, thread_id, round(prio, digits=3), round(t, digits=2), ndigits(length(contingencies)), string(type)[1], eta, master_var_vals[:η][c_id], CONSTS[:INTERNAL_MAPPINGS][:SLACK_VARIABLES][:SOLVER_NAME_MAP][Symbol(nameof(solver_type))])
        elseif type == :NumericalError
            count_sym   = :num_errors
            PARAMS[:PRINTING][:BENDERS][:LOG_LINE] && PARAMS[:PRINTING][:BENDERS][:LOG_LINE_UNSUCC] && @printf("%s%5d | %5d | %3d | %+.2f | %.2f | %*s | %s | %9.2f | %4s |\n", PAD, c_id, pos, thread_id, round(prio, digits=3), round(t, digits=2), ndigits(length(contingencies)), "N", "         ", master_var_vals[:η][c_id], CONSTS[:INTERNAL_MAPPINGS][:SLACK_VARIABLES][:SOLVER_NAME_MAP][Symbol(nameof(solver_type))])
        else
            @assert(type == :Nothing)
            count_sym   = :n_misses
            PARAMS[:PRINTING][:BENDERS][:LOG_LINE] && PARAMS[:PRINTING][:BENDERS][:LOG_LINE_UNSUCC] && @printf("%s%5d | %5d | %3d | %+.2f | %.2f | %*s | %s | %9.2f | %4s |\n", PAD, c_id, pos, thread_id, round(prio, digits=3), round(t, digits=2), ndigits(length(contingencies)), " ", "         ", master_var_vals[:η][c_id], CONSTS[:INTERNAL_MAPPINGS][:SLACK_VARIABLES][:SOLVER_NAME_MAP][Symbol(nameof(solver_type))])
        end
        subproblem_prio_list[pos][:n_calls] += 1

        # print :load, :sw_score, :reliability, :total_share, :n_since_call of the subproblem subproblem_prio_list[pos]
        #@printf("%s%5d | load = %5.3f | sw_score = %5.3f | rel = %5.3f | tot_share = %5.3f | n_since_call = %3d | total_viol = %12d\n", PAD, c_id, subproblem_prio_list[pos][:load], subproblem_prio_list[pos][:sw_score], subproblem_prio_list[pos][:reliability], subproblem_prio_list[pos][:total_share], subproblem_prio_list[pos][:n_calls], subproblem_prio_list[pos][:sum_violations])

        # update prio based on whether a cut was found
        _post_func(subproblem_prio_list[pos])

        # check whether this subproblem missed, hit only within μ_robustness, or hit fully, and add to global counters
        Base.acquire(counters_sem)
        counters[count_sym] += 1
        counters[:n_feascuts] += Int(type == :Feasibility && count_sym == :n_hits) # only count feasibility cuts if they are actually hits
        counters[:obj_value] += violation

        n_hits   = counters[:n_hits]
        n_misses = counters[:n_misses]
        n_feascuts = counters[:n_feascuts]
        n_interrupt_requests = counters[:n_interrupt_requests]
        obj_value = counters[:obj_value]
        Base.release(counters_sem)
        
        # check if we can interrupt forcefully
        if n_interrupt_requests >= 1
            throw(InterruptException())
        # check if we should stop early (partial solving)
        elseif PARAMS[:BENDERS][:PARTIAL][:ENABLED] && n_hits > 0

            # check if there is a stabilization round (only applicable if no feasibility cuts were found)
            if (PARAMS[:BENDERS][:PARTIAL][:STAB_ROUND_FREQ] < Inf) && # a limit should be applied
                (mod(callback_stats[:n_callbacks] - 1, PARAMS[:BENDERS][:PARTIAL][:STAB_ROUND_FREQ]) == 0) && # we are in a stabilization round
                n_feascuts == 0 # we have not found any feasibility cuts
                continue
            end
            
            if PARAMS[:BENDERS][:PARTIAL][:LIMITS][:STOP_ON_PB] && ! isnothing(cb_data)
                sol_count = Ref{Cint}()
                GRBcbget(cb_data, GRB_CB_MIPSOL, GRB_CB_MIPSOL_SOLCNT, sol_count)
                primal_bound = Ref{Cdouble}()
                GRBcbget(cb_data, GRB_CB_MIPSOL, GRB_CB_MIPSOL_OBJBST, primal_bound)
            end

            if      n_hits >= PARAMS[:BENDERS][:PARTIAL][:LIMITS][:STOP_AFTER_N_HITS]
                #println(" hits = $(counters[:n_hits]) >= $(PARAMS[:BENDERS][:PARTIAL][:LIMITS][:STOP_AFTER_N_HITS]) = stop_after_n_hits (with n_unrobust_hits = $(counters[:n_unrobust_hits]))")
                break
            elseif  n_hits >= PARAMS[:BENDERS][:PARTIAL][:LIMITS][:STOP_AFTER_N_HITS_REL] * length(contingencies)
                #println(" hits = $(counters[:n_hits]) >= $(PARAMS[:BENDERS][:PARTIAL][:LIMITS][:STOP_AFTER_N_HITS]) * length(contingencies) = stop_after_n_hits_rel (with n_unrobust_hits = $(counters[:n_unrobust_hits]))")
                break
            elseif  n_misses >= PARAMS[:BENDERS][:PARTIAL][:LIMITS][:STOP_AFTER_N_MISSES]
                #println(" misses = $(counters[:n_misses]) >= $(PARAMS[:BENDERS][:PARTIAL][:LIMITS][:STOP_AFTER_N_MISSES]) = stop_after_n_misses")
                break
            elseif  prio < PARAMS[:BENDERS][:PARTIAL][:LIMITS][:MIN_PRIO]
                #println(" prio = $(prio) < $min_prio = min_prio")
                break
            elseif time() - start_time >= PARAMS[:BENDERS][:PARTIAL][:LIMITS][:STOP_AFTER_TIME]
                #println(" time = $(time() - start_time) > $(PARAMS[:BENDERS][:PARTIAL][:TIME_LIMIT]) = time_limit")
                break 
            elseif n_feascuts >= PARAMS[:BENDERS][:PARTIAL][:LIMITS][:STOP_AFTER_N_FEAS_CUTS]
                #println(" n_feas_cuts = $(n_feas_cuts) >= $(PARAMS[:BENDERS][:PARTIAL][:LIMITS][:STOP_AFTER_N_FEAS_CUTS]) = stop_after_n_feas_cuts")
                break
            elseif PARAMS[:BENDERS][:PARTIAL][:LIMITS][:STOP_ON_PB] && ! isnothing(cb_data) && (sol_count[] > 0 && obj_value * (1 - PARAMS[:BENDERS][:PARTIAL][:LIMITS][:STOP_ON_PB_REL_BUFFER]) > primal_bound[] + SOLVER_PRECISION)
                break
            end
        end
    end

    reset_master!(solver, pm_sub, master_var_vals, nw=nw)
    return cons, violations, subproblem_ids, cut_types
end

"wait for all threads to finish, check for interrupts. Return true if hard interrupt was requested."
function wait_for_threads(threads, counters::Dict{Symbol, Real}, counters_sem::Base.Semaphore)
    while true
        try
            # if all threads are done, break the loop
            for t in threads
                Threads.wait(t)
            end
            break
        catch e
            # if an interrupt was signaled, add an interrupt request
            if isa(e, InterruptException)
                Base.acquire(counters_sem)
                counters[:n_interrupt_requests] += 1
                Base.release(counters_sem)
            else
                rethrow(e)
            end

            # if 1 request was made (0 not possible, since we're here), the loop will terminate after at least one cut was found
            if counters[:n_interrupt_requests] >= 1
                @warn("Interrupt requested. Benders loop will be interrupted after the currently running $(length(threads)) subproblems are finished.")
                for t in threads
                    Threads.wait(t)
                end
                rethrow(e)
            end

            break
        end
    end
end

"builds a 'new' heuristic solution based on the current master solution: x->x^*, η->new value"
function build_benders_solution(master_vars::Dict, η_val_news::Vector, cb_data, η_key::Symbol; master_var_vals=Dict{Symbol, AbstractVector}())::Tuple{Vector{JuMP.VariableRef}, Vector{Float64}}
    variables = Vector{JuMP.VariableRef}()
    values = Vector{Float64}()

    # keep the value of all normal master vars
    for (key, vars) in master_vars
        
        for var in vars
            if key ≠ η_key
                push!(variables, var)
                
                if isempty(master_var_vals) # no master_var_vals given, so we can use the callback data
                    value = JuMP.callback_value(cb_data, var)
                else
                    key = get_var_basename(var)
                    idx = get_var_idx(var)
                    value = master_var_vals[key][idx]
                end

                push!(values, value)
            end
        end
    end

    # update all eta vars to their new value
    for (c_id, value) in enumerate(η_val_news)
        push!(variables, master_vars[η_key][c_id])
        push!(values,    value)
    end
    
    return variables, values
end

"Creates a thread-local copy of lhsss, sharing VariableRef keys but copying AffExpr objects (and their constants) to avoid race conditions when apply_adjusted_rate_a! mutates constants."
function copy_lhsss(lhsss::Vector{Dict{Symbol, Dict{Int, JuMP.GenericAffExpr}}})::Vector{Dict{Symbol, Dict{Int, JuMP.GenericAffExpr}}}
    return [Dict(
        con_symbol => Dict(
            comp_id => copy(lhs_val) for (comp_id, lhs_val) in lhs_dict
        )
        for (con_symbol, lhs_dict) in lhss
    ) for lhss in lhsss]
end

"Spawns threads to solve subproblems in parallel. Called once per callback."
function solve_subproblems_multi(n_threads::Int64, 
                                 pm_subs,
                                 cb_data::Union{Gurobi.CallbackData, Nothing},
                                 contingencies::Vector{Contingency}, 
                                 sub_cons_for_dualss::Vector{Dict{Symbol, Dict{Int, JuMP.ConstraintRef}}}, 
                                 master_varss::Vector{Dict{Symbol, AbstractVector}},                         # one per master (nested)
                                 master_var_valss::Vector{Dict{Symbol, AbstractVector}},
                                 master_var_valss_var::Vector{Dict{VariableRef, Float64}},                   # same values as master_var_valss
                                 lhsss::Vector{Dict{Symbol, Dict{Int, JuMP.GenericAffExpr}}},                # one per master (nested)
                                 dualsss::Vector{Dict{Symbol, Dict{Int, Float64}}},                          # allocated with dummy values
                                 subproblem_prio_list::Vector{Dict{Symbol, Any}},
                                 solvers::Vector{<:AbstractSubproblemSolver},
                                 cut_exprss::Vector{Vector{AffExpr}};                                        # (partially) allocated for `cut_expr ≤ η`, one per master (nested)
                                 nw::Int64=nw_id_default,
                                 _pre_func!::Function=_prio_pre_weighted_sum!,
                                 int_sol::Bool=true,
                                 prec_mult=1.0,
                                 validation_status=nothing,
                                 validated_repair=nothing
                                )::Tuple{Vector{Vector{JuMP.ScalarConstraint}}, Vector{Float64}, Vector{Int64}, Vector{Tuple{Vector{VariableRef}, Vector{Float64}}}, Int, Int, Float64}
    
    callback_time = time()
    sol_obj_value = Ref{Cdouble}() # objective value of current MIPSOL solution
    isnothing(cb_data) ? sol_obj_value[] = 0.0 : GRBcbget(cb_data, GRB_CB_MIPSOL, GRB_CB_MIPSOL_OBJ, sol_obj_value)

    threads         = []  # list of threads that each get one copy of the subproblem, sub_cons_for_dualss, and master_var_valss
    cur_pos         = [1] # current position in the priority queue of subproblems (Vector since an Int would be copied, a vector is only referenced)
    counters        = Dict(:n_hits => 0, :n_misses => 0, :n_feascuts => 0, :n_unrobust_hits => 0, :n_misses_after_first_hit => 0, :num_errors => 0, :n_interrupt_requests => 0, :obj_value => sol_obj_value[]) # counters for hits and misses to terminate correctly
    cur_pos_sem     = Base.Semaphore(1) # semaphore to write to cur_pos
    counters_sem    = Base.Semaphore(1) # semaphore to write to counters dict

    # set up subproblem priority list for current master solution
    if ! isempty(subproblem_prio_list)
        _pre_func!(subproblem_prio_list, counters, data=ref(pm_subs[1], nw))
        sort!(subproblem_prio_list, by=(entry->entry[:priority]), rev=true)
    end

    # print subproblem statistics table header
    PARAMS[:PRINTING][:BENDERS][:LOG_LINE_HEADER] && println("$PAD  #SP |  Pos  | Thr | Score | Time | T |         η |  η_master | Slvr | C")
    
    # create per-thread copies of lhsss to avoid race conditions when apply_adjusted_rate_a! mutates shared constants
    n_actual_threads = min(n_threads, length(pm_subs))
    lhsss_per_thread = [thread_id == 1 ? lhsss : copy_lhsss(lhsss) for thread_id in 1:n_actual_threads]

    # initialize all threads that will each solve multiple different subproblems over the course of solving
    for thread_id in 1:n_actual_threads
        t = Threads.@spawn solve_subproblem_multi(thread_id, pm_subs[thread_id], cb_data, subproblem_prio_list, 
                contingencies, sub_cons_for_dualss[thread_id], master_varss, master_var_valss[thread_id], master_var_valss_var[thread_id], lhsss_per_thread[thread_id], dualsss[thread_id], cur_pos, cur_pos_sem, counters, counters_sem, solvers[thread_id], cut_exprss[thread_id], nw=nw, int_sol=int_sol, prec_mult=prec_mult)
        push!(threads, t)
    end
    # wait for all threads to finish, check for interrupts
    wait_for_threads(threads, counters, counters_sem)
    
    # get cuts from all subproblems
    data::Vector{Tuple{Vector{Vector{JuMP.ScalarConstraint}}, Vector{Float64}, Vector{Int64}, Vector{Symbol}}} = fetch.(threads)
    cons::Vector{Vector{JuMP.ScalarConstraint}} = vcat([t[1] for t in data]...)
    violations::Vector{Float64}                 = vcat([t[2] for t in data]...)
    subproblem_ids::Vector{Int64}               = vcat([t[3] for t in data]...)
    cut_types::Vector{Symbol}                   = vcat([t[4] for t in data]...)
    
    # print subproblem statistics summary
    n_solved_subs       = counters[:n_hits] + counters[:n_unrobust_hits] + counters[:n_misses]
    n_generated_cuts    = counters[:n_hits] + counters[:n_unrobust_hits]
    if !isnothing(validation_status)
        validation_status[] = int_sol && n_solved_subs == length(contingencies) &&
                              counters[:num_errors] == 0 && n_generated_cuts == 0
    end
    if counters[:num_errors] > 0 
        if counters[:n_hits] == 0
            @error("Could not solve $(counters[:num_errors]) subproblems due to numerical troubles, and no other cut found.") 
        else
            @warn("Could not solve $(counters[:num_errors]) subproblems due to numerical troubles, but we were able to recover with another cut.") 
        end
    end

    # if we have cuts, cut selection is enabled, and we are not in a cut stabilization round, select cuts
    (cuts, n_feas_selected, n_opt_selected, sum_selected_violations) = select_cuts!(ref(pm_subs[1], nw), subproblem_prio_list, subproblem_ids, cons, violations, cb_data, cut_types, n_callbacks=callback_stats[:n_callbacks])
    #@info("Selected $(length(cuts)) cuts from $(length(cons)) generated cuts with $(counters[:n_hits]) hits, $(counters[:n_unrobust_hits]) unrobust hits, and $(counters[:n_misses]) misses.")

    callback_stats[:n_solved_subs_total]  += n_solved_subs
    callback_stats[:n_benders_cuts_found] += n_generated_cuts
    callback_stats[:n_benders_cuts_added] += length(cuts)
    
    # A donor can stream the fully checked recourse repair even when Gurobi's
    # own incumbent heuristic does not submit it or accepts it without another
    # no-cut MIPSOL callback. Never repair feasibility or incomplete separation.
    if !isnothing(validated_repair)
        validated_repair[] = nothing
        if int_sol && n_solved_subs == length(contingencies) && counters[:num_errors] == 0 &&
                !(:Feasibility in cut_types)
            eta_values = [master_var_valss[1][:η][i] for i in eachindex(contingencies)]
            for (index, s_id) in enumerate(subproblem_ids)
                eta_values[s_id] += violations[index]
            end
            validated_repair[] = build_benders_solution(master_varss[1], eta_values, cb_data, :η;
                                                        master_var_vals=master_var_valss[1])
        end
    end

    # we might have a new solution to provide to the solver
    sols = Vector{Tuple{Vector{VariableRef}, Vector{Float64}}}()
    if PARAMS[:BENDERS][:SUBMIT_HEUR_SOLUTION] && n_solved_subs == length(contingencies)
        if ! isnothing(cb_data) 
            sol_count = Ref{Cint}()
            GRBcbget(cb_data, GRB_CB_MIPSOL, GRB_CB_MIPSOL_SOLCNT, sol_count)
            primal_bound = Ref{Cdouble}()
            GRBcbget(cb_data, GRB_CB_MIPSOL, GRB_CB_MIPSOL_OBJBST, primal_bound)
        end

        # the first solution or a better solution
        if isnothing(cb_data) || sol_count[] == 0 || counters[:obj_value] < sol_obj_value[]
            η_values = [master_var_valss[1][:η][i] for i in 1:length(contingencies)]
            for (index, s_id) in enumerate(subproblem_ids)
                η_values[s_id] += violations[index]
            end
            
            variables, values = build_benders_solution(master_varss[1], η_values, cb_data, :η, master_var_vals = master_var_valss[1])
            
            sols = [(variables, values)]
        end
    end
    
    (PARAMS[:PRINTING][:BENDERS][:LOG_LINE] || PARAMS[:PRINTING][:BENDERS][:LOG_LINE_SUMMARY]) && @printf("%sΣ%4d |       | %3d |       | %.2f | %*d | %9.2f | %9.2f |      | 🎯 %6.2f%%\n", PAD, n_solved_subs, n_threads, time() - callback_time, ndigits(length(contingencies)), n_generated_cuts, counters[:obj_value], sol_obj_value[], n_generated_cuts/n_solved_subs*100)

    return (cuts, violations, subproblem_ids, sols, n_feas_selected, n_opt_selected, sum_selected_violations)
end
function solve_subproblems_multi(n_threads::Int64, 
                                 pm_subs,
                                 cb_data,
                                 contingencies::Vector{Contingency}, 
                                 sub_cons_for_dualss::Vector{Dict{Symbol, Dict{Int, JuMP.ConstraintRef}}}, 
                                 master_vars::Dict{Symbol, AbstractVector}, 
                                 master_var_valss::Vector{Dict{Symbol, AbstractVector}},
                                 master_var_valss_var::Vector{Dict{VariableRef, Float64}},                   # same values as master_var_valss
                                 lhss::Dict{Symbol, Dict{Int, JuMP.GenericAffExpr}},
                                 dualsss::Vector{Dict{Symbol, Dict{Int, Float64}}},                          # allocated with dummy values
                                 subproblem_prio_list::Vector{Dict{Symbol, Any}},
                                 solvers::Vector{<:AbstractSubproblemSolver},
                                 cut_exprs::Vector{AffExpr};                                                 # (partially) allocated for `cut_expr ≤ η`
                                 nw::Int64=nw_id_default,
                                 _pre_func!::Function=_prio_pre_weighted_sum!,
                                 int_sol::Bool=true,
                                 prec_mult=1.0,
                                 validation_status=nothing,
                                 validated_repair=nothing
                                )::Tuple{Vector{JuMP.ScalarConstraint}, Vector{Float64}, Vector{Int64}, Vector{Tuple{Vector{VariableRef}, Vector{Float64}}}, Int, Int, Float64}
    (cuts, violations, subproblem_ids, sols, n_feas_selected, n_opt_selected, sum_selected_violations) = solve_subproblems_multi(n_threads, pm_subs, cb_data, contingencies, sub_cons_for_dualss, [master_vars], master_var_valss, master_var_valss_var, [lhss], dualsss, subproblem_prio_list, solvers, [[cut_expr] for cut_expr in cut_exprs], nw=nw, _pre_func! =_pre_func!, int_sol=int_sol, prec_mult=prec_mult, validation_status=validation_status, validated_repair=validated_repair)

    return ([cut[1] for cut in cuts], violations, subproblem_ids, sols, n_feas_selected, n_opt_selected, sum_selected_violations)
end


"Return constraints that prohibit the current solution if the graph is not connected"
function generate_connectedness_con(pm_subs::Vector, graphs::Vector, subproblem_prio_list, master_vars, master_var_vals::Dict, switched_off_by_master::Vector{Int64}, switched_on_by_master::Vector{Int64}; nw::Int64=nw_id_default)#::Tuple{Vector{JuMP.ScalarConstraint}, Vector{Contingency}}
    # TODO: make use of the multiple pm_sub and g for multithreading
    pm_sub          = pm_subs[1]
    g               = graphs[1]

    if PARAMS[:LIMITS][:BRANCH][:BRANCH_SWITCHING]
        off_branches = [branch_id for branch_id in axes(master_var_vals[:z_branch], 1) if isapprox(master_var_vals[:z_branch][branch_id], 0.0, atol=SOLVER_PRECISION)]
    else
        off_branches = [branch["index"] for (_, branch) in ref(pm_sub, nw, :branch) if branch["status"] == 0]
    end
    busbars = pm_sub.data["BUSBARS"]
    # `graph` includes ALL components (also isolated/unpowered ones)
    graph = _create_underlying_graph(pm_sub.data["nw"]["$nw"], merge_substations=busbars, off_branches=Vector{String}(string.(off_branches)), only_largest_component=false)

    base_powered_disconnected = count_powered_components(
        graph,
        pm_sub.data["nw"]["$nw"]["bus"],
        pm_sub.data["nw"]["$nw"]["load"],
        pm_sub.data["nw"]["$nw"]["gen"],
        early_stop=2,
    ) > 1

    splitting_contingencies_powered = Vector{Contingency}()
    splitting_contingencies = Vector{Contingency}()
    # base case
    if base_powered_disconnected
        if !isempty(subproblem_prio_list)
            splitting_contingencies = [subproblem_prio_list[1][:contingency]] # (plus an arbitrary contingency to strengthen the cut)
            splitting_contingencies_powered = [subproblem_prio_list[1][:contingency]]
        end
    else
        # graph may still contain benign unpowered isolated components
        graph_main = _keep_largest_component!(deepcopy(graph))
        cont = [e[:contingency] for e in subproblem_prio_list]
        splitting_contingencies_powered_new, splitting_contingencies_new = get_splitting_contingencies(graph_main, cont, ref(pm_sub, nw), only_powered=true)
        append!(splitting_contingencies_powered, splitting_contingencies_powered_new)
        # connected in all contingencies & base case
        if isempty(splitting_contingencies_powered)
            # check if it is also connected if we consider the topology itself, even if split up components are unpowered
            # we need that information to later on run LEQ directly for these contingencies
            append!(splitting_contingencies, splitting_contingencies_new)

            return Vector{JuMP.ScalarConstraint}(), splitting_contingencies
        end
    end

    # data for min cut cuts
    min_cut, capacitated = PARAMS[:BENDERS][:CONN_CUT_MIN_CUT], PARAMS[:BENDERS][:CONN_CUT_CAPACITATED]
    if min_cut || capacitated

        branch_data = ref(pm_sub, nw, :branch)

        off_branch_buses = Dict{Int64, Tuple{Int64, Int64}}(b_id => (branch_data[b_id]["f_bus"], branch_data[b_id]["t_bus"]) for b_id in off_branches)
        powered_buses = _powered_buses(_bus_load_gen(ref(pm_sub, nw))[2:3]...)
        single_contingency_branches = Set{Int64}(only(e[:contingency][:branch]) for e in subproblem_prio_list if keys(e[:contingency]) == Set([:branch]) && length(e[:contingency][:branch]) == 1)

    end

    # check base case (no contingency), then each contingency case
    cons = Vector{JuMP.ScalarConstraint}()
    for contingency in splitting_contingencies_powered
        disabled_branches_wo_contingency = [b_id for b_id in off_branches if b_id::Int64 ∉ contingency[:branch]::Vector{Int64}]
        disabled_branches_contingency = ["$b_id" for b_id in contingency[:branch]]

        for b_id in contingency[:branch]
            if haskey(ref(pm_sub, nw, :branch, b_id), "orig_index")
                push!(disabled_branches_contingency, "$(ref(pm_sub, nw, :branch, b_id)["orig_index"])")
            end
        end
        disabled_branches_bid_fbus_tbus  = [["$b_id"; string.(_get_branch_fbus_tbus(b_id, pm_sub, nw=nw))] for b_id in disabled_branches_wo_contingency]
        st_size = _underlying_graph_steiner_tree_approx(pm_sub, master_var_vals, graph, disabled_branches_contingency, disabled_branches_bid_fbus_tbus, nw=nw, only_powered=true)

        if st_size > 0
            if !base_powered_disconnected && isempty(switched_off_by_master) && !any(b_id ∈ switched_on_by_master for b_id in contingency[:branch])
                error("Contingency $contingency (without any additional switching) leads to unconnectedness. This contingency should have been filtered out in the preprocessing!")
            else
                ## create connectedness cuts
                # get the buses on this side of the min cut
                side = (min_cut || capacitated) ? get_min_cut_side(graph, contingency, off_branch_buses, powered_buses) : nothing

                if isnothing(side)
                    # at least one off branch needs to be switched on
                    _push_new!(cons, JuMP.@build_constraint(sum(master_vars[:z_branch][b_id] for b_id in disabled_branches_wo_contingency) >= st_size))

                else
                    # crosses are all branches crossing from "this" side over, cut_branches excludes contingencies
                    crosses = b_id -> (branch_data[b_id]["f_bus"] ∈ side) ≠ (branch_data[b_id]["t_bus"] ∈ side)
                    cut_branches = [b_id for b_id in disabled_branches_wo_contingency if crosses(b_id)]

                    # capacitated cut (non-contingency has capacity 2, contingency has 1, need at least 2)
                    cap_con = capacitated ? _capacitated_connectedness_cut([cut_branches; [b_id for b_id in contingency[:branch] if haskey(branch_data, b_id) && crosses(b_id)]], single_contingency_branches, master_vars, master_var_vals) : nothing

                    isnothing(cap_con) || _push_new!(cons, cap_con)

                    # min cut (at least one branch in the min cut needs to be turned on)
                    (min_cut || isnothing(cap_con)) && _push_new!(cons, JuMP.@build_constraint(sum(master_vars[:z_branch][b_id] for b_id in cut_branches) >= 1))
                end
                PARAMS[:PRINTING][:BENDERS][:LOG_LINE] && println("Contingency $contingency: Graph is not connected. Adding cut.")
                PARAMS[:BENDERS][:CONN_CUT_ALL_CONTINGENCIES] || return cons, splitting_contingencies
            end
        end
    end
    return cons, splitting_contingencies
end

"push `con` to `cons` if no equal constraint is contained"
function _push_new!(cons::Vector{JuMP.ScalarConstraint}, con::JuMP.ScalarConstraint)
    any(c -> c.func == con.func && c.set == con.set, cons) || push!(cons, con)
    return cons
end

"capacitated connectedness cut over a cut separating powered buses"
function _capacitated_connectedness_cut(cut_branches::Vector{Int64}, unsafe_branches::Set{Int64}, master_vars, master_var_vals::Dict)::Union{JuMP.ScalarConstraint, Nothing}
    z_branch, z_branch_val = master_vars[:z_branch], master_var_vals[:z_branch]
    terms = Vector{Tuple{Float64, JuMP.VariableRef}}()
    rhs, value = 2.0, 0.0

    for b_id in cut_branches
        coeff = b_id ∈ unsafe_branches ? 1.0 : 2.0

        if b_id ∈ axes(z_branch, 1)
            push!(terms, (coeff, z_branch[b_id]))
            value += coeff * z_branch_val[b_id]
        else
            rhs -= coeff # not switchable, i.e., on
        end

    end

    value < rhs - SOLVER_PRECISION || return nothing
    
    return JuMP.@build_constraint(sum(coeff * z for (coeff, z) in terms) >= rhs)
end

"fix switching variables based on connectedness"
function fix_forced_on_branches!(master_vars::Dict, data::Dict, contingencies::Vector{Contingency})::Int64
    (!PARAMS[:BENDERS][:FIX_FORCED_ON_BRANCHES] || !haskey(master_vars, :z_branch) || isempty(contingencies)) && return 0

    z_branch = master_vars[:z_branch]
    n_fixed = 0
    for branch_id in get_forced_on_branches(data, contingencies)
        
        if branch_id ∈ axes(z_branch, 1) && !JuMP.is_fixed(z_branch[branch_id])
            JuMP.fix(z_branch[branch_id], 1.0; force=true)
            n_fixed += 1
        end
    end
    
    return n_fixed
end

"add hint to subproblem priority list to use fallback solver for certain contingencies"
function add_feascut_hint_to_subproblem_prio_list!(subproblem_prio_list, contingencies_hint_feascut::Vector{Contingency}, subproblem_solver_type::Type)
    # check if subproblem_solver_type is a solver that is not able to generate feasibility cuts (i.e., an AbstractMatrixSolver)
    if subproblem_solver_type <: AbstractMatrixSolver
        for i in eachindex(subproblem_prio_list)
            subproblem_prio_list[i][:hint_feascut] = subproblem_prio_list[i][:contingency] ∈ contingencies_hint_feascut
        end
    end
end

"add b to all electrical branches if busbars, else add b to all branches"
function add_branch_b!(pm::AbstractPowerModel; nw::Int64=nw)
    # add `b`s to all branches
    for (_, branch) in ref(pm, nw, :branch)
        if pm.data["BUSBARS"]
            if is_branch_electrical(branch)
                _, branch["b"] = calc_branch_y(branch)
            end
        else
            _, branch["b"] = calc_branch_y(branch)
        end
    end
end

"check whether the subproblem solver type is supported, and modify if required (necessary for busbars)"
function check_subproblem_solver_type(pm, subproblem_solver_type::Type)
    # currently not supported
    if pm.data["BUSBARS"] && subproblem_solver_type ∉ [FullLPSolver, LinEqSysBBSolver, LinEqSysBBRfSolver]
        if subproblem_solver_type == LinEqSysSolver
            new_solver_type = LinEqSysBBSolver
        elseif subproblem_solver_type == LinEqSysRfSolver
            new_solver_type = LinEqSysBBRfSolver
        else
            new_solver_type = FullLPSolver
        end
        
        @info("Busbars not supported with solver type $subproblem_solver_type. Setting to $new_solver_type.")
        subproblem_solver_type = new_solver_type
        PARAMS[:BENDERS][:SUBPROBLEM_SOLVER_TYPE] = new_solver_type
    end

    return subproblem_solver_type
end

"initialize fractional solvers (only FullLP and LP supported)"
function init_frac_solvers(subproblem_solver_type, pm, pm_subs, n_threads, sub_cons_for_dualss; nw=nw) 
    if subproblem_solver_type == LPSolver
        frac_solvers = init_lp_solvers(pm, n_threads, nw=nw)
    else
        frac_solvers = init_fulllp_solvers(pm_subs, n_threads, sub_cons_for_dualss)
    end
end

"initialize the integer solvers"
function init_int_solvers(subproblem_solver_type, pm, pm_subs, n_threads, sub_cons_for_dualss; nw=nw)
    if subproblem_solver_type == LPSolver
        int_solvers = init_lp_solvers(pm, n_threads, nw=nw)
    elseif subproblem_solver_type == LinEqSysSolver
        int_solvers = init_lineqsys_solvers(pm, pm_subs, n_threads, sub_cons_for_dualss, nw=nw)
    elseif subproblem_solver_type == LinEqSysBBSolver
        int_solvers = init_lineqsysbb_solvers(pm, pm_subs, n_threads, sub_cons_for_dualss, nw=nw)
    #elseif subproblem_solver_type == CholeskySolver
        #int_solvers = init_cholesky_solver(pm, n_threads, nw=nw)
    elseif subproblem_solver_type == FullLPSolver
        int_solvers = init_fulllp_solvers(pm_subs, n_threads, sub_cons_for_dualss)
    elseif subproblem_solver_type == LinEqSysRfSolver
        int_solvers = init_lineqsysrf_solvers(pm, pm_subs, n_threads, sub_cons_for_dualss, nw=nw)
    elseif subproblem_solver_type == LinEqSysBBRfSolver
        int_solvers = init_lineqsysbbrf_solvers(pm, pm_subs, n_threads, sub_cons_for_dualss, nw=nw)
    else
        throw(ArgumentError("Unknown subproblem solver type $(subproblem_solver_type)!"))
    end
end

"get lhs of the constraints to be multiplied with duals"
function get_empty_lhss_dualss(sub_cons_for_duals, master_vars)::Tuple{Dict{Symbol, Dict{Int, JuMP.GenericAffExpr}}, Dict{Symbol, Dict{Int, Float64}}}
    lhss   = Dict{Symbol, Dict{Int, JuMP.GenericAffExpr}}()
    dualss = Dict{Symbol, Dict{Int, Float64}}()
    for (con_symbol, cons_dict) in sub_cons_for_duals
        lhss[con_symbol] = get_lhs(cons_dict, master_vars)
        dualss[con_symbol] = Dict(comp_id=>0.0 for comp_id in keys(lhss[con_symbol]))
    end
    return lhss, dualss
end

"get lhs to be multiplied with duals, pre-allocate duals, for multiple sub_cons_for_duals"
function get_empty_lhss_dualsss(sub_cons_for_dualss, master_vars)::Tuple{Dict{Symbol, Dict{Int, JuMP.GenericAffExpr}}, Vector{Dict{Symbol, Dict{Int, Float64}}}}
    dualsss = Vector{Dict{Symbol, Dict{Int, Float64}}}()

    lhss, dualss = get_empty_lhss_dualss(sub_cons_for_dualss[1], master_vars)
    push!(dualsss, dualss)

    for _ in 1:length(sub_cons_for_dualss)-1
        push!(dualsss, deepcopy(dualss))
    end
    return lhss, dualsss
end

"get lhs to be multiplied with duals, pre-allocate duals, for multiple sub_cons_for_duals"
function get_empty_lhsss_dualsss(sub_cons_for_dualss, master_vars)::Tuple{Vector{Dict{Symbol, Dict{Int, JuMP.GenericAffExpr}}}, Vector{Dict{Symbol, Dict{Int, Float64}}}}
    lhsss =   Vector{Dict{Symbol, Dict{Int, JuMP.GenericAffExpr}}}()
    dualsss = Vector{Dict{Symbol, Dict{Int, Float64}}}()
    for sub_cons_for_duals in sub_cons_for_dualss
        lhss, dualss = get_empty_lhss_dualss(sub_cons_for_duals, master_vars)
        push!(lhsss, lhss)
        push!(dualsss, dualss)
    end
    return lhsss, dualsss
end

"using in-out strengthening, we relax the integrality of the master problem variables, solve the LP relaxation, and add cuts to the master problem"
function strengthen_master_lp_formulation!(pm, pm_subs, contingencies, n_threads, sub_cons_for_dualss, master_vars, lhss, dualsss, subproblem_prio_list, solvers, cut_exprs::Vector{AffExpr}; nw=1, retain_integrality=true)
    "relax integrality of the master problem variables (assuming that all variables of one type have the same domain), and return the map of their type to their domain"
    function relax_integrality!(model::JuMP.Model, master_vars)
        # fix variables and unset their binary status
        var_type_map = Dict{Symbol, Symbol}()
        for (var_name, vars) in master_vars
            if any(JuMP.is_binary.(vars))
                JuMP.unset_binary.(vars)
                var_type_map[var_name] = :binary
            elseif any(JuMP.is_integer.(vars))
                JuMP.unset_integer.(vars)
                var_type_map[var_name] = :integer
            else
                var_type_map[var_name] = :continuous
            end
        end
        return var_type_map
    end

    "using the map of variable types, set the variables back to their original type"
    function enforce_integrality!(model::JuMP.Model, var_type_map)
        # set variables back to their original type
        for (var_name, var_type) in var_type_map
            if var_type == :binary
                JuMP.set_binary.(var(pm, nw, var_name))
            elseif var_type == :integer
                JuMP.set_integer.(var(pm, nw, var_name))
            end
        end
    end

    "get master variable values"
    function get_master_var_valss(master_vars, n_threads)::Vector{Dict{Symbol, AbstractVector}}
        master_var_valss = Dict{Symbol, AbstractVector}()
        # get the values of the master variables for each thread
        for (key, vars) in master_vars
            master_var_valss[key] = JuMP.value.(vars)
        end
        return [copy(master_var_valss) for _ in 1:n_threads]
    end

    "optimize relaxed master problem LP relaxation without switching"
    function get_master_var_vals_no_rd_no_sw(pm, n_threads, master_vars; nw=nw_id_default)
        no_switch = JuMP.@constraint(pm.model, sum(master_vars[:z_branch_change]) <= 0)
        no_rd = JuMP.@constraint(pm.model, sum(var(pm, nw, :pg_delta_on_off)) <= 0)

        JuMP.optimize!(pm.model)

        master_var_valss_no_switch_no_rd = get_master_var_valss(master_vars, n_threads)

        delete(pm.model, no_switch)
        delete(pm.model, no_rd)

        return master_var_valss_no_switch_no_rd
    end

    "get master variable values initialized to ones"
    function get_master_var_vals_ones(master_vars)
        # get the values of the master variables for each thread, initialized to ones
        master_var_valss_ones = Dict{Symbol, AbstractVector}()
        for (key, vars) in master_vars
            master_var_valss_ones[key] = ones(length(vars))
        end
        return master_var_valss_ones
    end

    function print_iter_and_update_stall_counter(pm, n_iters, strengthened_lp_obj_value, n_added_in_out_cuts, t, dual_stall_phase, dual_stall_counter)
        # print and check for required phase change
        if n_iters == 0
            unstrengthened_lp_obj_value = JuMP.objective_value(pm.model)
        else
            prev_strengthened_lp_obj_value = strengthened_lp_obj_value
            strengthened_lp_obj_value = JuMP.objective_value(pm.model)
            @printf(" %2d | %8.2f | %15.8f | %5d |\n", n_iters, time() - t, strengthened_lp_obj_value, n_added_in_out_cuts)

            if isapprox(strengthened_lp_obj_value, prev_strengthened_lp_obj_value, atol=0.001)
                dual_stall_counter += 1

                if dual_stall_counter >= PARAMS[:BENDERS][:IN_OUT_STRENGTHENING][:STALL_PHASE_SWITCH_LIM]
                    println("Dual stall detected, switching to dual stall phase $(dual_stall_phase + 1).")
                    dual_stall_phase += 1
                    dual_stall_counter = 0
                end
            end
        end
        return dual_stall_phase, dual_stall_counter, strengthened_lp_obj_value, unstrengthened_lp_obj_value
    end

    "update master_var_vals_stable to 0.5 * master_var_vals_stable + 0.5 * master_var_vals and then set master_var_vals to lambda * master_var_vals + (1 - lambda) * master_var_vals_stable + delta * 1 (the vector where all vars are set to one)"
    function set_in_out_master_var_valss!(master_var_valss::Vector{Dict{Symbol, AbstractVector}},
                                          master_var_vals_stable::Dict{Symbol, AbstractVector},
                                          master_var_vals_initial::Dict{Symbol, AbstractVector},
                                          var_type_map::Dict{Symbol, Symbol},
                                          dual_stall_phase::Int64;
                                          lambda::Float64=PARAMS[:BENDERS][:IN_OUT_STRENGTHENING][:LAMBDA], 
                                          delta::Float64=PARAMS[:BENDERS][:IN_OUT_STRENGTHENING][:DELTA],
                                          round_int_vars::Bool=true)

        # dual_stall_phase key:
        # 0: set lambda to setting, set delta to 0.0002
        # 1: set lambda to 1, set delta to 0.0002
        # 2: set lambda to 1, set delta to 0.0
        # 3: only for one iteration until the stopping criterion is met

        if dual_stall_phase == 0
            λ = lambda
            δ = delta
        elseif dual_stall_phase == 1
            λ = 1.0
            δ = delta
        elseif dual_stall_phase >= 2
            λ = 1.0
            δ = 0.0
        end

        # move master_var_vals_stable halfway towards master_var_valss
        for (key, var_vals) in master_var_valss[1]
            master_var_vals_stable[key] = 0.5 .* var_vals .+ 0.5 .* master_var_vals_stable[key]
        end

        # set master_var_valss to the "intermediate point" (linear combination of the stabilizing point and initial point)
        for master_var_vals in master_var_valss
            for (key, var_vals) in master_var_vals
                # for all variables, we combine the values of the first thread with the current values
                var_vals = λ .* var_vals .+ (1 - λ) .* master_var_vals_stable[key] .+ δ .* 1#master_var_vals_initial[key]

                # floor all integer vars according to var_type_map
                if round_int_vars
                    if var_type_map[key] == :integer || var_type_map[key] == :binary
                        var_vals = floor.(var_vals)
                    end
                end
            end
        end
    end

    # print header line for the in-out strengthening
    @info("Strengthening master problem LP formulation with in-out cuts.")
    @printf(" %5s | %8s | %15s | %5s |\n", "it", "time", "LP obj.", "ncuts")
    t = time()
    n_iters = 0
    n_added_in_out_cuts = 0
    unstrengthened_lp_obj_value = 0.0
    strengthened_lp_obj_value = 0.0
    dual_stall_counter = 0
    dual_stall_phase = 0
    
    # disable logging of solver attached to pm.model
    JuMP.set_optimizer_attribute(pm.model, "LogToConsole", 0)

    # record domain of variables and relax them
    var_type_map = relax_integrality!(pm.model, master_vars)

    # optimize master problem LP relaxation with switching
    JuMP.optimize!(pm.model)
    master_var_valss_initial = get_master_var_valss(master_vars, n_threads)
    master_var_valss = get_master_var_valss(master_vars, n_threads)
    master_var_vals_ones = get_master_var_vals_ones(master_vars)
    #master_var_valss_no_sw_no_rd = get_master_var_vals_no_rd_no_sw(pm, n_threads, master_vars, nw=nw)

    # set params for in-out strengthening
    STOP_ON_PB_copy = copy(PARAMS[:BENDERS][:PARTIAL][:LIMITS][:STOP_ON_PB])
    STOP_AFTER_N_HITS_copy = copy(PARAMS[:BENDERS][:PARTIAL][:LIMITS][:STOP_AFTER_N_HITS])
    STOP_AFTER_TIME_copy = copy(PARAMS[:BENDERS][:PARTIAL][:LIMITS][:STOP_AFTER_TIME])
    #CUT_SELECTION_copy = copy(PARAMS[:BENDERS][:CUTS][:SELECTION][:ENABLED])
    TAL_TIMESTEP_copy = copy(PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP])
    PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP] = 100000
    PARAMS[:BENDERS][:PARTIAL][:LIMITS][:STOP_ON_PB] = true
    #PARAMS[:BENDERS][:CUTS][:SELECTION][:ENABLED] = false
    PARAMS[:BENDERS][:PARTIAL][:LIMITS][:STOP_AFTER_N_HITS] = PARAMS[:BENDERS][:IN_OUT_STRENGTHENING][:N_HITS_LIM]
    PARAMS[:BENDERS][:PARTIAL][:LIMITS][:STOP_AFTER_TIME] = PARAMS[:BENDERS][:IN_OUT_STRENGTHENING][:ITER_TIME_LIM]

    while n_iters < PARAMS[:BENDERS][:IN_OUT_STRENGTHENING][:N_ITER_LIM] && dual_stall_phase < 3
        dual_stall_phase, dual_stall_counter, strengthened_lp_obj_value, unstrengthened_lp_obj_value = print_iter_and_update_stall_counter(pm, n_iters, strengthened_lp_obj_value, n_added_in_out_cuts, t, dual_stall_phase, dual_stall_counter)

        if ! isempty(pm_subs)
            set_in_out_master_var_valss!(master_var_valss, master_var_vals_ones, master_var_valss_initial[1], var_type_map, dual_stall_phase, round_int_vars=retain_integrality)

            master_var_valss_var = get_master_var_vals_vars(master_vars, master_var_valss[1], n_threads)

            # add line loads of current master solution to the subproblem priority list
            _add_data_to_subproblem_prio_list!.(subproblem_prio_list, [pm], [get(var(pm, nw), :p, Dict())], nw=nw)

            # solve powerflow subproblems for all contingencies
            (powerflow_cons, violations, _, sols, _, _, _) = solve_subproblems_multi(n_threads, 
                                                    pm_subs,
                                                    nothing,
                                                    contingencies, 
                                                    sub_cons_for_dualss, 
                                                    master_vars, 
                                                    master_var_valss,
                                                    master_var_valss_var,
                                                    lhss, 
                                                    dualsss, 
                                                    subproblem_prio_list, 
                                                    solvers,
                                                    cut_exprs,
                                                    nw=nw, 
                                                    int_sol=retain_integrality)

            # add the cuts to the master problem
            if ! isempty(powerflow_cons)
                for con in powerflow_cons
                    n_added_in_out_cuts += 1
                    JuMP.add_constraint(pm.model, con, "benders_in_out_cut_$(n_added_in_out_cuts)")
                end
            end
        end
        n_iters += 1

        # get the LP solution for printing and warmstart and the next iteration
        JuMP.optimize!(pm.model)
        master_var_valss = get_master_var_valss(master_vars, n_threads)
    end
    @info("Added $n_added_in_out_cuts in-out cuts to the master problem LP formulation, strengthening it from $unstrengthened_lp_obj_value to $(strengthened_lp_obj_value).")

    # reset params for in-out strengthening
    PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP] = TAL_TIMESTEP_copy
    #PARAMS[:BENDERS][:CUTS][:SELECTION][:ENABLED] = CUT_SELECTION_copy
    PARAMS[:BENDERS][:PARTIAL][:LIMITS][:STOP_ON_PB] = STOP_ON_PB_copy
    PARAMS[:BENDERS][:PARTIAL][:LIMITS][:STOP_AFTER_N_HITS] = STOP_AFTER_N_HITS_copy
    PARAMS[:BENDERS][:PARTIAL][:LIMITS][:STOP_AFTER_TIME] = STOP_AFTER_TIME_copy

    # enforce integrality of the master problem variables for full solving afterwards
    enforce_integrality!(pm.model, var_type_map)        
    # re-enable logging of solver attached to pm.model
    JuMP.set_optimizer_attribute(pm.model, "LogToConsole", PARAMS[:PRINTING][:SOLVER][:LOG_TO_CONSOLE_MP])
end

"add master variable values to the pool, remove oldest entry if pool size exceeded"
function _add_master_var_valss_to_pool!(master_var_vals_pool::Vector{Dict{Symbol, AbstractVector}}, master_var_vals::Dict{Symbol, AbstractVector}; max_pool_size::Int64=PARAMS[:BENDERS][:CUT_STRENGTHENING][:HISTORY][:POOL_SIZE])
    # pop the last entry if pool would overflow
    length(master_var_vals_pool) >= max_pool_size && pop!(master_var_vals_pool)

    # push the new master variable values to the front of the pool
    pushfirst!(master_var_vals_pool, master_var_vals)
end

"linearly combine master_var_vals with the pool, using the depreciation rate"
function combine_master_var_valss(master_var_vals::Dict{Symbol, AbstractVector}, master_var_vals_pool::Vector{Dict{Symbol, AbstractVector}}; round_ints::Bool=PARAMS[:BENDERS][:CUT_STRENGTHENING][:ROUND_INTS])::Dict{Symbol, AbstractVector}
    combined_master_var_vals = copy(master_var_vals)

    if PARAMS[:BENDERS][:CUT_STRENGTHENING][:HISTORY][:WEIGHTING_MODE] == :EQUAL # equal for each solution (including the first one)
        sol_weights = [1/length(master_var_vals_pool) for _ in 1:length(master_var_vals_pool)]
    elseif PARAMS[:BENDERS][:CUT_STRENGTHENING][:HISTORY][:WEIGHTING_MODE] == :EXPONENTIAL # 2^-i = 0.5, 0.25, 0.125, ...
        sol_weights = [(2.0^-i) for i in 1:length(master_var_vals_pool)]
    else
        throw("Unknown master variable values weighting mode $(PARAMS[:BENDERS][:CUT_STRENGTHENING][:HISTORY][:WEIGHTING_MODE]).")
    end

    # scale up sol_weights such that they sum up to 1
    sol_weights = sol_weights ./ sum(sol_weights)

    # combine the pool of master variable values with the current master variable values
    for i in eachindex(master_var_vals_pool)
        for (key, var_vals) in master_var_vals_pool[i]
            for index in eachindex(var_vals)
                if i == 1
                    # combined_master_var_vals already contains the first solution, so we just scale it, not add it
                    combined_master_var_vals[key][index] *= sol_weights[i]
                else
                    combined_master_var_vals[key][index] += sol_weights[i] * var_vals[index]
                end
            end 
        end
    end

    # if round_ints, round the values to the nearest integer
    if round_ints
        # round the values to the nearest integer
        for (key, var_vals) in combined_master_var_vals
            if key in [:z_branch, :z_gen]
                # round the values to the nearest integer
                combined_master_var_vals[key] = round.(Int, var_vals) # round the values to the nearest integer
            end
        end
    end

    return combined_master_var_vals
end

"main function to create master and sub problem and define Benders callback"
function build_secure_ots_time_benders(pm::AbstractPowerModel; nw::Int64=1, n_threads::Int64=PARAMS[:SOLVER][:MULTITHREADING][:N_THREADS_BENDERS], subproblem_solver_type::Type = PARAMS[:BENDERS][:SUBPROBLEM_SOLVER_TYPE])::AbstractPowerModel
    # save original contingencies and build priority list
    contingencies = ref(pm, nw, :contingencies) # list of contingencies that represent the subproblems
    contingencies_full = ref(pm, nw, :contingencies) # to restore later
    n_contingencies = length(contingencies)
    n_threads = min(n_threads, length(contingencies))
    ref(pm, nw)[:contingencies] = Vector{Contingency}()

    # # Partial Benders Decomposition: leave some contingencies in the master problem
    if PARAMS[:BENDERS][:PARTIAL_DECOMP][:ENABLED]
        n_keep_req = n_master_contingencies(n_contingencies)
        if n_keep_req > 0 # if any contingencies (and not just possibly the base case) should be kept, do this
            contingencies_relevant = get(ref(pm, nw), :contingencies_relevant, Vector())
            if !isempty(contingencies_relevant)
                n_keep = min(n_keep_req, length(contingencies_relevant))
                contingencies_keep = contingencies_relevant[1:n_keep]
                ref(pm, nw)[:contingencies] = contingencies_keep # move kept contingencies back into the master
                contingencies = setdiff(contingencies, contingencies_keep) # reduce subproblem contingency list by the kept contingencies
                @info("Partial Benders Decomposition: Leaving $(length(contingencies_keep)) contingencies in the master problem, leaving the rest in the original contingency list for decomposition.")
                @assert(length(contingencies) + length(ref(pm, nw)[:contingencies]) == n_contingencies, "Not all contingencies were accounted for")
            else
                @warn("Partial Benders Decomposition is enabled but no relevant contingencies were detected. Partial decomposition will be skipped.")
            end
        end
    end

    # track scores and order subproblems accordingly
    subproblem_prio_list::Vector{Dict{Symbol, Any}} = [Dict(:contingency_id=>c, 
                                                            :contingency=>cont, 
                                                            :priority=>0.0, 
                                                            :load=>0.0, 
                                                            :reliability => 0.0, 
                                                            :n_calls => 0.0,
                                                            :diversity => 0.0,
                                                            :n_cuts => 0.0,
                                                            :n_cuts_used => 0.0, # in case we want to limit the number of cuts used
                                                            :sum_violations => 0.0,
                                                            :round_called => -1) for (c, cont) in enumerate(contingencies)]

    # check subproblem solver compatibility for current settings
    subproblem_solver_type = check_subproblem_solver_type(pm, subproblem_solver_type)
    
    # add required subproblem solver data
    if subproblem_solver_type <: AbstractMatrixSolver # TODO: or if any fallback solver is an AbstractMatrixSolver
        add_branch_b!(pm; nw=nw)
    end
    printst("Benders", "Preprocessing done.")

    # build problems as full LP models (TODO: only do that if FullLP is included as solver or any fallback solver)
    pm_subs             = Vector{AbstractPowerModel}(undef, n_threads)
    Threads.@threads for i in 1:n_threads
        pm_subs[i] = create_subproblem!(deepcopy(pm), nw=nw)
    end
    sub_cons_for_dualss = get_sub_cons_for_duals.(pm_subs)
    if get(pm.data, "_collect_admittance_scores", false) && !isempty(pm_subs)
        pm.ext[:admittance_subproblem] = first(pm_subs)
        first(pm_subs).ext[:admittance_cuts] = Any[]
    end
    printst("Benders", "Built $n_threads subproblem$(n_threads > 1 ? "s" : "").")

    # if contingencies exist, initialize subproblem solvers
    if !isempty(contingencies)
        int_solvers = init_int_solvers(subproblem_solver_type, pm, pm_subs, n_threads, sub_cons_for_dualss, nw=nw)
        if get(pm.data, "_collect_admittance_scores", false)
            pm.ext[:admittance_solver] = first(int_solvers)
        end

        # if cuts on fractional solutions should be generated, initialize the solvers
        if PARAMS[:BENDERS][:CUTS][:FRACTIONAL_SEPA_NODELIM] >= -1 || PARAMS[:BENDERS][:IN_OUT_STRENGTHENING][:ENABLED] || PARAMS[:BENDERS][:CUT_STRENGTHENING][:ENABLED]
           frac_solvers = init_frac_solvers(subproblem_solver_type, pm, pm_subs, n_threads, sub_cons_for_dualss, nw=nw)
        end
    end

    pm = create_master_problem!(pm, contingencies, nw=nw)
    printst("Benders", "Built master problem.")

    if PARAMS[:BRANCHING][:BRANCHING_PRIORITY_MODE] == :SW_HEUR_PRIO
        add_branching_priorities(pm, nw=nw)
    end

    master_vars = get_master_vars(pm, nw=nw)
    fix_forced_on_branches!(master_vars, ref(pm, nw), contingencies)

    graphs = [_create_underlying_graph(ref(pm, nw)) for _ in 1:n_threads]

    # get lhs to be multiplied with duals, pre-allocate duals, pre-allocate cut expressions
    if ! isempty(contingencies)
        lhss, dualsss = get_empty_lhss_dualsss(sub_cons_for_dualss, master_vars)
    end

    cut_exprs = [JuMP.@expression(pm.model, JuMP.AffExpr(0)) for _ in 1:n_threads]::Vector{AffExpr}
    reset_cut_exprs!(cut_exprs, master_vars)

    global callback_stats = init_callback_stats(nested=false)
    global callback_times = init_callback_times(nested=false)

    # add some cuts to the master problem
    PARAMS[:BENDERS][:IN_OUT_STRENGTHENING][:ENABLED] && strengthen_master_lp_formulation!(pm, pm_subs, contingencies, n_threads, sub_cons_for_dualss, master_vars, lhss, dualsss, subproblem_prio_list, frac_solvers, cut_exprs, nw=nw)

    master_var_vals_pool = Vector{Dict{Symbol, AbstractVector}}()#undef, PARAMS[:BENDERS][:CUT_STRENGTHENING][:HISTORY][:POOL_SIZE])
    
    grb_mip_log = false
    grb_warmstart_success = false
    powerflow_cons = Vector{JuMP.ScalarConstraint}()
    violations = Vector{Float64}()
    n_feas_selected = 0
    n_opt_selected = 0
    sum_selected_violations = 0.0

    master_var_vals = Dict{Symbol, AbstractVector}()
    master_var_vals_var = Dict{VariableRef, Float64}()
    incumbent_exchange = init_incumbent_exchange(pm; contingencies)
    if !isnothing(incumbent_exchange)
        pm.data["_parallel_incumbent_exchange"] = incumbent_exchange
    end

    "update master_var_vals_var and master_var_vals with the values from the callback data"
    function update_master_var_vals_var!(cb_data, int_sol::Bool)
        if isempty(master_var_vals)
            for (key, vars) in master_vars
                master_var_vals[key] = JuMP.callback_value.(cb_data, vars)
            end
            for (_, vars) in master_vars
                for var in vars
                    master_var_vals_var[var] = JuMP.callback_value.(cb_data, var)
                end
            end
        else
            for (key, vars) in master_vars
                _master_var_vals = master_var_vals[key]

                for id in axes(vars, 1)
                    var = vars[id]
                    
                    val = JuMP.callback_value(cb_data, var)

                    _master_var_vals[id] = val
                    master_var_vals_var[var] = val
                end
            end
        end

        # round discrete separation values; keep raw values for cut evaluation
        for (key, vars) in master_vars
            for id in axes(vars, 1)
                var = vars[id]
                if int_sol && (JuMP.is_binary(var) || JuMP.is_integer(var))
                    val = master_var_vals[key][id]
                    @assert(isfinite(val) && abs(val - round(val)) <= 1e-5, "Nonintegral discrete value in integer separation: $val")
                    master_var_vals[key][id] = round(val)
                end
            end
        end
    end


    "callback that adds Benders cuts to the master problem"
    function benders_callback(cb_data, cb_where::Cint)
        exchange_callback!(incumbent_exchange, pm, cb_data, cb_where) === :stop && return
        # Check that we are in a valid fractional or integer node
        if cb_where != GRB_CB_MIPSOL && cb_where != GRB_CB_MIPNODE && cb_where != GRB_CB_MESSAGE
            return
        end

        # Check whether the cuts we added prevented a new incumbent solution from being accepted
        if cb_where == GRB_CB_MESSAGE
            # get message text
            resultP = Ref{Cstring}()
            GRBcbget(cb_data, cb_where, GRB_CB_MSG_STRING, resultP)
            GRB_MSG = unsafe_string(resultP[])
            # if we are not yet in the MIP log section, check if we are now
            if ! grb_mip_log
                if startswith(GRB_MSG, "    Nodes    |    Current Node    ")
                    grb_mip_log = true
                elseif startswith(GRB_MSG, "Loaded user MIP start with objective") || ! pm.data["WARMSTART_PROVIDED"] || PARAMS[:BENDERS][:PARTIAL_DECOMP][:ENABLED] # TODO check partial decomp compatibility
                    grb_warmstart_success = true
                end
                return
            # if we are in the MIP log section, check if the line starts with "H" or "*" (indicating a newly found incumbent solution)
            else
                if ! grb_warmstart_success
                    throw(WarmstartError("Provided warmstart was not accepted by Gurobi."))
                    grb_warmstart_success = true # set this to true to not print this line after it was printed once
                end
                if GRB_MSG[1] == 'H' || GRB_MSG[1] == '*'
                    @debug("Incumbent found. Previous iteration: $(length(powerflow_cons)) cuts ($(n_feas_selected) feas, $(n_opt_selected) opt), sum_selected_violations=$(sum_selected_violations)")
                    if sum_selected_violations > PARAMS[:BENDERS][:CHECKS][:CHECK_SOL_ACCEPT_TOL]
                        @error("New incumbent solution accepted even though $(length(powerflow_cons)) new cuts ($(n_feas_selected) feas, $(n_opt_selected) opt) were added with a cumulated selected violation of $(sum_selected_violations) in the previous iteration.\n Solution might be incorrect.")
                        for cut in powerflow_cons
                            @error("Cut: $(cut)")
                            # write cut into file cut.lp
                            if false
                                open("cut.lp", "w") do io
                                    func = cut.func
                                    for (var, coeff) in func.terms
                                        if coeff != 0.0
                                            println(io, "$(coeff) * $(JuMP.name(var)) + ")
                                        end
                                    end
                                    rhs = func.constant
                                    println(io, " <= $(rhs)")
                                end
                            end
                        end
                    end
                end
            end
            return
        end

        # If we're in a fractional node, check that it is optimal and that we still want to separate on fractional solutions
        if cb_where == GRB_CB_MIPNODE
            resultP = Ref{Cint}()
            GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_STATUS, resultP)
            if resultP[] != GRB_OPTIMAL
                return  # Solution is something other than optimal.
            end
            
            # run rounding heuristic
            # Gurobi.load_callback_variable_primal(cb_data, cb_where)
            # master_var_valss        = Vector{Dict{Symbol, AbstractVector}}([Dict(key => JuMP.callback_value.(cb_data, vars) for (key, vars) in master_vars) for thread in 1:n_threads])
            # sols = scots_rounding_heuristic_all_frac(pm, master_vars, master_var_valss[1], cb_data)
            # for (variables, values) in sols
                # JuMP.MOI.submit(pm.model, JuMP.MOI.HeuristicSolution(cb_data), variables, values)
            # end

            # check if number of solved nodes is > PARAMS[:BENDERS][:CUTS][:FRACTIONAL_SEPA_NODELIM] , then return
            resultP = Ref{Cdouble}()
            GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_NODCNT, resultP)
            if (PARAMS[:BENDERS][:CUTS][:FRACTIONAL_SEPA_NODELIM] == -1 || resultP[] <= PARAMS[:BENDERS][:CUTS][:FRACTIONAL_SEPA_NODELIM]) && callback_stats[:n_frac_sepas] < PARAMS[:BENDERS][:CUTS][:FRACTIONAL_SEPA_ITERLIM]
                @debug("Node count: $(resultP[]) <= $(PARAMS[:BENDERS][:CUTS][:FRACTIONAL_SEPA_NODELIM]), separating on fractional solution.")
                solvers = frac_solvers
                int_sol = false
                callback_stats[:n_frac_sepas] += 1
            else
                return
            end
        elseif cb_where == GRB_CB_MIPSOL
            if !isempty(pm_subs)
                solvers = int_solvers
            end

            int_sol = true
        end

        callback_stats[:n_callbacks] += 1

        # do preprocessing required for both, connectedness and powerflow subproblems
        pre_start_time = time()
        Gurobi.load_callback_variable_primal(cb_data, cb_where)
        
        # compile or update dictionary that contains the variables and the values of the variables
        update_master_var_vals_var!(cb_data, int_sol)

        # Use deepcopy to avoid shared references between threads (for the vals, copy suffices, since nothing is mutated)
        master_var_valss::Vector{Dict{Symbol, AbstractVector}} = [deepcopy(master_var_vals) for _ in 1:n_threads]
        master_var_valss_var = [copy(master_var_vals_var) for _ in 1:n_threads]::Vector{Dict{VariableRef, Float64}}
        
        _add_master_var_valss_to_pool!(master_var_vals_pool, master_var_vals)

        (switched_off_by_master, switched_on_by_master)  = get_master_decisions(pm, cb_data, nw=nw)
        callback_times[:pre] += time() - pre_start_time

        # solve connectedness subproblems for all contingencies
        if int_sol && ! isempty(pm_subs) #&& PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP] > 0
            connect_start_time = time()
            connectedness_cons, contingencies_hint_feascut = generate_connectedness_con(pm_subs, graphs, subproblem_prio_list, master_vars, master_var_valss[1], switched_off_by_master, switched_on_by_master, nw=nw)
            callback_times[:connect] += time() - connect_start_time

            if ! isempty(connectedness_cons)
                if PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP] == 0
                    @warn("Connectedness cut added for network without any switching.")
                end
                for connectedness_con in connectedness_cons
                    JuMP.MOI.submit(pm.model, JuMP.MOI.LazyConstraint(cb_data), connectedness_con)
                end
                callback_stats[:n_connect_cuts_added] += length(connectedness_cons)
                return
            end
        else
            contingencies_hint_feascut = Vector{Contingency}()
        end

        t = time()
        validated_incumbent = Ref(int_sol && isempty(pm_subs))
        validated_repair = !isnothing(incumbent_exchange) && incumbent_exchange.role == :donor ? Ref{Any}(nothing) : nothing
        if ! isempty(pm_subs)
            # solve contingency subproblems
            _add_data_to_subproblem_prio_list!.(subproblem_prio_list, [pm], [cb_data], nw=nw)

            # add feascut hints to subproblem_prio_list
            add_feascut_hint_to_subproblem_prio_list!(subproblem_prio_list, contingencies_hint_feascut, subproblem_solver_type)

            # solve powerflow subproblems for all contingencies
            (powerflow_cons, violations, subproblem_ids, sols, n_feas_selected, n_opt_selected, sum_selected_violations) = solve_subproblems_multi(n_threads, 
                                                     pm_subs, 
                                                     cb_data, 
                                                     contingencies, 
                                                     sub_cons_for_dualss, 
                                                     master_vars, 
                                                     master_var_valss,
                                                     master_var_valss_var,
                                                     lhss, 
                                                     dualsss, 
                                                     subproblem_prio_list, 
                                                     solvers,
                                                     cut_exprs,
                                                     nw=1, 
                                                     int_sol=int_sol,
                                                     validation_status=validated_incumbent,
                                                     validated_repair=validated_repair)

            for con in powerflow_cons
                JuMP.MOI.submit(pm.model, JuMP.MOI.LazyConstraint(cb_data), con)
            end
            audit_exchange_cuts!(incumbent_exchange, powerflow_cons)
            # serialize cuts for carry-forward if enabled
            if PARAMS[:WARMSTARTING][:RSCOTS][:CARRY_FORWARD_CUTS] && haskey(pm.data, "_benders_cuts_pool")
                serialize_constraints!(pm.data["_benders_cuts_pool"], powerflow_cons, contingencies, subproblem_ids)
            end
            for (variables, values) in sols
                accept = JuMP.MOI.submit(pm.model, JuMP.MOI.HeuristicSolution(cb_data), variables, values)
            end
        end
        callback_times[:contingency] += time() - t
        validated_incumbent[] && publish_incumbent!(incumbent_exchange, pm, cb_data)
        !isnothing(validated_repair) && publish_repaired_incumbent!(incumbent_exchange, pm, cb_data, validated_repair[])

        # run Benders cut strengthening
        if PARAMS[:BENDERS][:CUT_STRENGTHENING][:ENABLED] && ! isempty(powerflow_cons) && length(master_var_vals_pool) > 1
            combined_master_var_vals = combine_master_var_valss(master_var_valss[1], master_var_vals_pool, round_ints = PARAMS[:BENDERS][:CUT_STRENGTHENING][:ROUND_INTS])
            combined_master_var_valss = [combined_master_var_vals for _ in 1:n_threads]

            n_hits_copy = copy(PARAMS[:BENDERS][:PARTIAL][:LIMITS][:STOP_AFTER_N_HITS])
            PARAMS[:BENDERS][:PARTIAL][:LIMITS][:STOP_AFTER_N_HITS] = 1

            # solve powerflow subproblems for all contingencies
            (powerflow_cons_strengthening, violations_strengthening, subproblem_ids_strengthening, _, _, _, _) = solve_subproblems_multi(n_threads, 
                                                     pm_subs, 
                                                     cb_data, 
                                                     contingencies, 
                                                     sub_cons_for_dualss, 
                                                     master_vars, 
                                                     combined_master_var_valss,
                                                     master_var_valss_var, 
                                                     lhss, 
                                                     dualsss, 
                                                     subproblem_prio_list, 
                                                     frac_solvers,
                                                     cut_exprs,
                                                     nw=1, 
                                                     int_sol=PARAMS[:BENDERS][:CUT_STRENGTHENING][:ROUND_INTS])

            # reset the number of hits
            PARAMS[:BENDERS][:PARTIAL][:LIMITS][:STOP_AFTER_N_HITS] = n_hits_copy

            for (i, con) in enumerate(powerflow_cons_strengthening)
                JuMP.MOI.submit(pm.model, JuMP.MOI.LazyConstraint(cb_data), con)
                @debug("[CS] Submitting valid inequality with violation $(violations_strengthening[i]).")
            end
            # serialize strengthening cuts for carry-forward if enabled
            if PARAMS[:WARMSTARTING][:RSCOTS][:CARRY_FORWARD_CUTS] && haskey(pm.data, "_benders_cuts_pool")
                serialize_constraints!(pm.data["_benders_cuts_pool"], powerflow_cons_strengthening, contingencies, subproblem_ids_strengthening)
            end
        
        end

        # run heuristic
        #resultP2 = Ref{Cdouble}()
        #GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_NODCNT, resultP2)
        #println("GRB_CB_MIPNODE_NODCNT = $(resultP2[])")
        if pm.data["WARMSTART_PROVIDED"] && PARAMS[:PRIMAL_HEURISTICS][:ENABLED] #&& resultP2[] > 0
            _, master_var_vals_modifieds = scots_rounding_heuristic_sw_score(pm, master_vars, master_var_valss[1], cb_data)

            for (i, master_var_vals_modified) in enumerate(master_var_vals_modifieds)
                # @info("Fixing heuristic solution $(i)/$(length(master_var_vals_modifieds)) with modified master variable values.")
                # # fix z_gen/z_branch/tr_shift(if TRANSFORMER_DISCRETE) variables in the master problem
                # for (var_name, vars) in master_vars
                #     if haskey(master_var_vals_modified, var_name)
                #         if var_name == :z_gen || var_name == :z_branch
                #             JuMP.unset_binary.(vars)
                #             JuMP.fix.(vars, master_var_vals_modified[var_name])
                #         elseif var_name == :tr_shift && PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_DISCRETE]
                #             JuMP.unset_integer.(vars)
                #             JuMP.fix.(vars, master_var_vals_modified[var_name])
                #         end
                #     end
                # end

                @info("Solving subproblems for heuristic solution $(i)/$(length(master_var_vals_modifieds)).")
                master_var_valss_modified = [master_var_vals_modifieds[i] for _ in 1:n_threads]

                PARAMS[:BENDERS][:PARTIAL][:ENABLED] = false
                (cuts, violations, subproblem_ids, sols, _, _, _) = solve_subproblems_multi(n_threads, 
                                            pm_subs, 
                                            nothing, 
                                            contingencies, 
                                            sub_cons_for_dualss, 
                                            master_vars, 
                                            master_var_valss_modified,
                                            master_var_valss_var, # TODO modify
                                            lhss, 
                                            dualsss, 
                                            subproblem_prio_list, 
                                            solvers,
                                            cut_exprs,
                                            nw=1, 
                                            int_sol=int_sol,
                                            _pre_func! = _prio_pre_weighted_sum!,
                                            prec_mult = 1.0)
                PARAMS[:BENDERS][:PARTIAL][:ENABLED] = true

                # add cuts that were separated on the heuristic solutions to the model
                @info("Adding $(length(cuts)) cuts (based on heur. sol.) to the master problem.")
                for (i, cut) in enumerate(cuts)
                    JuMP.MOI.submit(pm.model, JuMP.MOI.LazyConstraint(cb_data), cut)
                    sp_idx = subproblem_ids[i]
                    master_var_valss_modified[1][:η][sp_idx] = ceil(violations[i], digits=PARAMS[:PRECISION][:SOLVER_PREC_DIGITS]-1)
                end
                    
                # @info("Re-optimizing master problem to fix heuristic solution after adding cuts.")
                # JuMP.optimize!(pm.model) # re-optimize the master problem to fix the solution

                #@info("Updating master variable values after heuristic solution.")
                #master_var_valss = get_master_var_valss(master_vars, n_threads)
                variables, values = build_benders_solution(master_vars, 
                                                master_var_valss_modified[1][:η], 
                                                nothing, 
                                                :η, 
                                                master_var_vals = master_var_valss_modified[1])

                #@show(variables, values)
                @info("Adding heuristic solution (incl. eta values) to the master problem.")
                accepted = JuMP.MOI.submit(pm.model, JuMP.MOI.HeuristicSolution(cb_data), variables, values)

                # @info("Unfixing heuristic solution $(i)/$(length(master_var_vals_modifieds)).")
                # # fix z_gen/z_branch/tr_shift(if TRANSFORMER_DISCRETE) variables in the master problem
                # for (var_name, vars) in master_vars
                #     if haskey(master_var_vals_modified, var_name)
                #         if var_name == :z_gen || var_name == :z_branch
                #             JuMP.set_binary.(vars)
                #             JuMP.unfix.(vars)
                #         elseif var_name == :tr_shift && PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_DISCRETE]
                #             JuMP.set_integer.(vars)
                #             JuMP.unfix.(vars)
                #         end
                #     end
                # end
            end
        end
    end

    # add callback to solver
    JuMP.MOI.set(pm.model, Gurobi.CallbackFunction(), benders_callback)

    ref(pm, nw)[:contingencies] = contingencies_full
    return pm
end
