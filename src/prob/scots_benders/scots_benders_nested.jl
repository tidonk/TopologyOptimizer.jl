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
function solve_secure_ots_time_benders_nested(input::Any, model_type::Type, optimizer::Any; solver_specific_build_function=ssbf(pm)=pm, multinetwork=true, slack_types=[:pg, :ratea, :dclinepf, :vad], kwargs...)
    return solve_secure_model(input, model_type, optimizer, solver_specific_build_function ∘ build_secure_ots_time_benders_nested, slack_types; ref_extensions=[ref_add_on_off_va_bounds!], multinetwork=true, kwargs...)
end

"create Benders master problem for switching (top level problem)"
function create_master_switching!(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int64=nw_id_default, report::Bool=true)::AbstractPowerModel
    # set up gurobi
    gurobi_master_sw_env = Gurobi.Env()
    standard_params = Dict("LogToConsole"       => PARAMS[:PRINTING][:SOLVER][:LOG_TO_CONSOLE_MP], 
                           "LazyConstraints"   => 1, 
                           "Threads"           => PARAMS[:SOLVER][:MULTITHREADING][:N_THREADS_MIP_MP], 
                           "LogFile"           => get_log_file_name(pm.data),#joinpath(CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:SOLVER_LOG_DIR], "GRB_MP_$(pm.data["name"])-$(replace(pm.data["start_time"], " " => "_")).out"),
                           "TimeLimit"         => PARAMS[:SOLVER][:SOLVING_TIME_LIMIT],
                           "NumericFocus"      => 1,
                           "Heuristics"        => PARAMS[:BENDERS][:NESTED][:MP_PRIMAL_HEURISTICS],
                           "StartNodeLimit"    => pm.data["WARMSTART_PROVIDED"] ? -1 : -3
                          )
    
    if get(pm.data, "WARMSTART_PROVIDED", false)
        standard_params["Presolve"] = -1
    end
    if haskey(pm.data, "_MIPGap")
        standard_params["MIPGap"] = pm.data["_MIPGap"]
    end
    if haskey(pm.data, "_TimeLimit")
        standard_params["TimeLimit"] = pm.data["_TimeLimit"]
    end
    if haskey(pm.data, "_BestObjStop")
        standard_params["BestObjStop"] = pm.data["_BestObjStop"]
    end
    
    JuMP.set_optimizer(pm.model, optimizer_with_attributes(() -> Gurobi.Optimizer(gurobi_master_sw_env), standard_params..., ))

    # build switching master problem
    build_switching_model_time(pm, with_pf=PARAMS[:BENDERS][:BASE_PF_IN_MASTER])

    # define objective function
    weighted_obj_func_factory([objective_min_num_tr_switch, objective_min_num_branch_switch, objective_min_num_gen_switch],
        [pm.data["cost_coeffs"]["tr_sw"], pm.data["cost_coeffs"]["br_sw"], pm.data["cost_coeffs"]["gen_sw"]])(pm)
    @assert JuMP.objective_sense(pm.model) == JuMP.MIN_SENSE

    # add variable for Benders
    @assert PARAMS[:BENDERS][:NESTED][:ADD_IP_CUTS] || PARAMS[:BENDERS][:NESTED][:ADD_SP_CUTS]
    η_rd = var(pm, nw)[:η_rd] = JuMP.@variable(pm.model, 
                                        [1:1],
                                        base_name = "benders_eta_rd", 
                                        lower_bound = 0.0, 
                                        start = 0.0)

    PARAMS[:LIMITS][:GENERATOR][:GEN_SWITCHING] && constraint_benders_min_gen_cost(pm, nw)
    
    # Note: η_rd is a scalar variable representing total redispatch cost across all contingencies,
    # so we don't store it per-contingency in the solution

    lhs = JuMP.@expression(pm.model, 0.0)

    # add Benders variable to objective function
    JuMP.@objective(pm.model, Min,
        JuMP.objective_function(pm.model)
        + sum( var(pm, n, :η_rd)[1] for (n, nw_ref) in nws(pm) )
        )

    if PARAMS[:BENDERS][:NESTED][:ADD_SP_CUTS]
        @assert PARAMS[:BENDERS][:NESTED][:GENERATION_IN_MP]
        eta_start = haskey(ref(pm, nw), :η) ? ref(pm, nw, :η) : Dict(c => Dict("contingency" => contingency, "η_start" => 0.0) for (c, contingency) in enumerate(contingencies))
        η = var(pm, nw)[:η] = JuMP.@variable(pm.model, 
                                        [c in eachindex(contingencies)], 
                                        base_name = "benders_eta", 
                                        lower_bound = 0.0, 
                                        start = comp_start_value(eta_start[c], "η_start"))

        report && sol_component_value(pm, nw, :eta, :η, eachindex(contingencies), η)

        # add Benders variable to objective function
        #JuMP.@objective(pm.model, Min,
        #    JuMP.objective_function(pm.model)
        #    + sum(
        #        sum( var(pm, n, :η), init=0.0)
        #    for (n, nw_ref) in nws(pm))
        #  )

        lhs += sum(var(pm, nw, :η))
    end


    # gens are lb on eta if included
    if PARAMS[:BENDERS][:NESTED][:GENERATION_IN_MP]
        lhs += _expr_min_redispatch_cost(pm) * pm.data["cost_coeffs"]["redispatch"]
    else
        constraint_benders_min_gen_cost(pm, nw)
    end

    JuMP.@constraint(pm.model, lhs <= η_rd[1], base_name="benders_gen_cost")
    

    if PARAMS[:PRECISION][:SOLVER_PREC_DIGITS] >= 0 && PARAMS[:ROUNDING][:MODEL_COEFFS]
        pm.model = round_coeffs!(pm.model)
    end
    #@info("Writing file to mp_$(Dates.format(now(), "HH:MM:SS")).mps and mp_$(Dates.format(now(), "HH:MM:SS")).lp")
    #JuMP.write_to_file(pm.model, "mp_$(Dates.format(now(), "HH:MM:SS")).mps")
    #JuMP.write_to_file(pm.model, "mp_$(Dates.format(now(), "HH:MM:SS")).lp")
    return pm
end

"create Benders master problem for gens (second level problem)"
function create_master_gens!(pm::AbstractPowerModel, contingencies::Vector{Contingency}; nw::Int64=nw_id_default, report::Bool=true)::AbstractPowerModel
    # set up gurobi
    gurobi_master_gen_env = Gurobi.Env()
    standard_params = Dict("LogToConsole"       => PARAMS[:PRINTING][:SOLVER][:LOG_TO_CONSOLE_IP], 
                           "LazyConstraints"   => 1, 
                           "MIPFocus"          => 1,  # IMPORTANT!! with MIPSepCuts
                           "MIPSepCuts"        => 2,  # IMPORTANT!! with MIPFocus
                           "Threads"           => PARAMS[:SOLVER][:MULTITHREADING][:N_THREADS_MIP_MP], 
                           "DualReductions"    => 0,       # required to obtain InfUnbdInfo
                           "InfUnbdInfo"       => 1,       # required to know to add cut
                           "LogFile"           => get_log_file_name(pm.data),#joinpath(CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:SOLVER_LOG_DIR], "GRB_IP_$(pm.data["name"])-$(replace(pm.data["start_time"], " " => "_")).out"),
                           "TimeLimit"         => PARAMS[:SOLVER][:SOLVING_TIME_LIMIT],
                           "NumericFocus"      => 1
                          )

    JuMP.set_optimizer(pm.model, optimizer_with_attributes(() -> Gurobi.Optimizer(gurobi_master_gen_env), standard_params..., ))

    # build switching master problem
    build_gen_model_time(pm, with_pf=PARAMS[:BENDERS][:NESTED][:BASE_PF_IN_INTERMEDIATE])

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
    
    # define objective function
    if PARAMS[:BENDERS][:NESTED][:BASE_PF_IN_INTERMEDIATE]
        weighted_obj_func_factory([objective_min_redispatch_cost, objective_min_slack_relative_linear, objective_min_num_slacks], 
            [pm.data["cost_coeffs"]["redispatch"], pm.data["cost_coeffs"]["overload"], pm.data["cost_coeffs"]["num_overload"]])(pm)
    else
        weighted_obj_func_factory([objective_min_redispatch_cost], 
            [pm.data["cost_coeffs"]["redispatch"]])(pm)
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
    return pm

end

"main function to create master and sub problem and define Benders callback"
function build_secure_ots_time_benders_nested(pm::AbstractPowerModel; nw::Int64=1, n_threads::Int64=PARAMS[:SOLVER][:MULTITHREADING][:N_THREADS_BENDERS], subproblem_solver_type::Type = PARAMS[:BENDERS][:SUBPROBLEM_SOLVER_TYPE])::AbstractPowerModel
    # save original contingencies and build priority list
    contingencies = ref(pm, nw, :contingencies)
    n_threads     = min(n_threads, length(contingencies))
    ref(pm, nw)[:contingencies] = Vector{Contingency}()

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

    # build problems as full LP models (TODO: only do that if FullLP is included as solver or any fallback solver)
    pm_subs             = Vector{AbstractPowerModel}(undef, n_threads)
    Threads.@threads for i in 1:n_threads
        pm_subs[i] = create_subproblem!(deepcopy(pm), nw=nw)
    end
    sub_cons_for_dualss_gen = get_sub_cons_for_duals.(pm_subs)

    # if contingencies exist, initialize subproblem solvers
    if !isempty(contingencies)
        int_solvers = init_int_solvers(subproblem_solver_type, pm, pm_subs, n_threads, sub_cons_for_dualss_gen, nw=nw)

        # if cuts on fractional solutions should be generated, initialize the solvers
        if PARAMS[:BENDERS][:CUTS][:FRACTIONAL_SEPA_NODELIM] >= -1
            frac_solvers = init_frac_solvers(subproblem_solver_type, pm, pm_subs, n_threads, sub_cons_for_dualss_gen, nw=nw)
        end
    end
    
    pm_master_gens          = create_master_gens!(deepcopy(pm), contingencies, nw=nw)
    pm_master_switching     = create_master_switching!(pm, contingencies, nw=nw)

    switching_model = pm_master_switching.model
    gen_model       = pm_master_gens.model

    # warm start η_rd, otherwise the MIP start (η_rd = 0) can be rejected
    if get(pm.data, "WARMSTART_PROVIDED", false)
        η_rd_start = JuMP.value(v -> something(JuMP.start_value(v), 0.0), JuMP.objective_function(gen_model))
        JuMP.set_start_value(var(pm_master_switching, nw, :η_rd)[1], ceil(η_rd_start, digits=PARAMS[:PRECISION][:SOLVER_PREC_DIGITS]-2))
    end

    if PARAMS[:BRANCHING][:BRANCHING_PRIORITY_MODE] == :SW_HEUR_PRIO
        add_branching_priorities(pm_master_switching, nw=nw)
    end

    master_vars_switching, master_vars_gen = get_master_vars_nested(pm_master_switching, pm_master_gens, nw=nw)
    fix_forced_on_branches!(master_vars_switching, ref(pm_master_switching, nw), contingencies)
    gen_var_to_master_var_switching = Dict{VariableRef, VariableRef}(var=>get_master_var_from_var(var, master_vars_switching) for var in JuMP.all_variables(gen_model) if get_var_basename(var) ∈ keys(master_vars_switching))

    graphs = [_create_underlying_graph(ref(pm_master_switching, nw)) for _ in 1:n_threads]

    # track number of callback calls
    global callback_stats = init_callback_stats(nested=true)
    global callback_times = init_callback_times(nested=true)
        
    grb_mip_log = false
    grb_warmstart_success = false
    violations = Vector{Float64}()

    # extract all relevant constraints
    sub_cons_for_duals_switching = get_sub_cons_for_duals_gen_model(pm_master_gens, with_pf=PARAMS[:BENDERS][:NESTED][:BASE_PF_IN_INTERMEDIATE])
    lhss_gen, dualsss_gen = get_empty_lhss_dualsss(sub_cons_for_dualss_gen, master_vars_gen) # get lhss and dualss to use for cut to add to gen master (i.e., cuts from contingency subproblems)

    lhsss_sp = [lhss_gen]
    master_varss_sp = [master_vars_gen]
    if PARAMS[:BENDERS][:NESTED][:ADD_SP_CUTS]
        lhss_switching_sp, _  = get_empty_lhss_dualsss(sub_cons_for_dualss_gen, master_vars_switching) # get lhss and dualss to use for cut to add to switching master (i.e., cuts from contingency subproblems)

        push!(lhsss_sp,        lhss_switching_sp)
        push!(master_varss_sp, master_vars_switching)
    end
    lhss_switching, dualss_switching = get_empty_lhss_dualss(sub_cons_for_duals_switching, master_vars_switching) # get lhss and dualss to use for cut to add to switching master

    cut_expr_switching = JuMP.@expression(switching_model, JuMP.AffExpr(0))
    reset_cut_expr!(cut_expr_switching, master_vars_switching)

    cut_exprs_gen = [JuMP.@expression(gen_model, JuMP.AffExpr(0)) for _ in 1:n_threads]::Vector{AffExpr}
    reset_cut_exprs!(cut_exprs_gen, master_vars_gen)
    
    cut_exprss_sp = Vector{Vector{AffExpr}}()
    if PARAMS[:BENDERS][:NESTED][:ADD_SP_CUTS]
        for i in 1:n_threads
            cut_expr_switching_sp = JuMP.@expression(switching_model, JuMP.AffExpr(0))
            reset_cut_expr!(cut_expr_switching_sp, master_vars_switching)

            push!(cut_exprss_sp, [cut_exprs_gen[i], cut_expr_switching_sp])
        end
    else
        cut_exprss_sp = [[cut_exprs_gen[i]] for i in 1:n_threads]
    end

    best_η_val = Inf
    n_pool_cuts = 0  # number of benders cuts kept in the intermediate problem across callbacks


    cons = []
    presolved = false

    "callback that adds Benders cuts to the master problem"
    function nested_benders_callback(cb_data, cb_where::Cint)
        # Check that we are in a valid fractional or integer node
        if cb_where == GRB_CB_PRESOLVE
            presolved = true
        end

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
                elseif startswith(GRB_MSG, "Loaded user MIP start with objective") || ! pm.data["WARMSTART_PROVIDED"]
                    grb_warmstart_success = true
                end
                return
            # if we are in the MIP log section, check if the line starts with "H" or "*" (indicating a newly found incumbent solution)
            else
                if ! grb_warmstart_success
                    @warn("Provided warmstart was not accepted by Gurobi.")
                    grb_warmstart_success = true # set this to true to not print this line after it was printed once
                end
                if GRB_MSG[1] == 'H' || GRB_MSG[1] == '*'
                    if ! isempty(violations)
                        @warn("New incumbent solution accepted even though $(length(violations)) new cuts were added with a cumulated violation of $(sum(violations)) in the previous iteration.\n Solution might be incorrect, or the cut costs might be integrated in the solution cost.")
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
            
            # check if number of solved nodes is > FRACTIONAL_SEPA_NODELIM, then return
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
        # compile dictionary that contains the variables and a dictionary that contains a flag indicating whether they are formulated as integer, and the values of the variables
        master_var_vals_switching::Dict{Symbol, AbstractVector} = Dict(key => JuMP.callback_value.(cb_data, vars) for (key, vars) in master_vars_switching )
        callback_times[:pre] += time() - pre_start_time
        
        # solve connectedness subproblems for all contingencies
        if int_sol && ! isempty(pm_subs) #&& PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP] > 0
            connect_start_time = time()
            
            (switched_off_by_master, switched_on_by_master)  = get_master_decisions(pm_master_switching, cb_data, nw=nw)
            connectedness_cons, contingencies_hint_feascut = generate_connectedness_con(pm_subs, graphs, subproblem_prio_list, master_vars_switching, master_var_vals_switching, switched_off_by_master, switched_on_by_master, nw=nw)

            callback_times[:connect] += time() - connect_start_time
            if ! isempty(connectedness_cons)
                if PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP] == 0
                    @warn("Connectedness cut added for network without any switching.")
                end
                for connectedness_con in connectedness_cons
                    JuMP.MOI.submit(switching_model, JuMP.MOI.LazyConstraint(cb_data), connectedness_con)
                end
                callback_stats[:n_connect_cuts_added] += length(connectedness_cons)
                return
            end
        end

        ##############################################################################################################
        # We now solve the gen master problem (with the fixed master var vals from the switching master), and then solve the powerflow subproblems for all contingencies.
        # With the obtained cuts from the powerflow subproblems, we add them to the gen master problem and resolve.
        # We can resolve after each obtained cut, if we want.
        # Once we have resolved the gen master problem to optimality, and it yielded an objective function value that is larger than the current switching master η, we can stop.
        # We then build the cut using the dual variable values from the gen master problem and add it to the switching master problem.
        ##############################################################################################################  
        
        start_time = time()
        
        @assert(length(master_vars_gen[:z_branch]) == length(master_var_vals_switching[:z_branch]), "Some branches have z_branch in switching model, but not in gen model")
        
        # fix switching master variables in the gen master
        if PARAMS[:LIMITS][:BRANCH][:BRANCH_SWITCHING]
            for branch_id::Int in axes(master_var_vals_switching[:z_branch],1)
                JuMP.fix(master_vars_gen[:z_branch][branch_id], master_var_vals_switching[:z_branch][branch_id], force=true)
            end
        end
        if PARAMS[:LIMITS][:GENERATOR][:GEN_SWITCHING] 
            for gen_id::Int in axes(master_var_vals_switching[:z_gen],1)
                JuMP.fix(master_vars_gen[:z_gen][gen_id], master_var_vals_switching[:z_gen][gen_id], force=true)
            end
        end
        if PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_SHIFTING]
            for tr_id::Int in axes(master_var_vals_switching[:tr_shift],1)
                JuMP.fix(master_vars_gen[:tr_shift][tr_id], master_var_vals_switching[:tr_shift][tr_id], force=true)
            end
        end

        # show which branch was switched
        if PARAMS[:PRINTING][:BENDERS][:LOG_LINE_INTERPROBLEM] && PARAMS[:PRINTING][:BENDERS][:LOG_LINE_INTERPROBLEM_SW]
            sw_str = "🔌"
            first = true
            for b_id::Int in axes(master_vars_switching[:z_branch], 1)
                if master_var_vals_switching[:z_branch_change][b_id] > SOLVER_PRECISION
                    sw_str = "$sw_str$(first ? "" : ",")$b_id"
                    first = false
                end
            end
            if PARAMS[:LIMITS][:GENERATOR][:GEN_SWITCHING]
                sw_str = "$(sw_str)🔋"
                first = true
                for gen_id::Int in axes(master_var_vals_switching[:z_gen],1)
                    gen = ref(pm, nw, :gen, gen_id)
                    if master_var_vals_switching[:z_gen][gen_id] > SOLVER_PRECISION && gen["status"] == 0
                        sw_str = "$sw_str$(first ? "" : ",")$gen_id"
                        first = false
                    end
                end
                sw_str = "$(sw_str)🪫 "
                first = true
                for gen_id::Int in axes(master_var_vals_switching[:z_gen],1)
                    gen = ref(pm, nw, :gen, gen_id)
                    if master_var_vals_switching[:z_gen][gen_id] < 1 - SOLVER_PRECISION && gen["status"] == 1
                        sw_str = "$sw_str$(first ? "" : ",")$gen_id"
                        first = false
                    end
                end
            end
        else
            sw_str = ""
        end
        
        sol_obj_value = Ref{Cdouble}() # objective value of current MIPSOL solution
        sol_count = Ref{Cint}()
        primal_bound = Ref{Cdouble}()
        if cb_where == GRB_CB_MIPSOL
            GRBcbget(cb_data, cb_where, GRB_CB_MIPSOL_OBJ, sol_obj_value)
            GRBcbget(cb_data, cb_where, GRB_CB_MIPSOL_SOLCNT, sol_count)
            GRBcbget(cb_data, cb_where, GRB_CB_MIPSOL_OBJBST, primal_bound)
        else
            GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_OBJ, sol_obj_value)
            GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_SOLCNT, sol_count)
            GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_OBJBST, primal_bound)
        end
            
        if sol_count[] > 0
            best_η_val = min(best_η_val, primal_bound[])
        end

        # get Benders variable value from master ("switching master")
        if PARAMS[:BENDERS][:NESTED][:ADD_IP_CUTS]
            η_val_master = master_var_vals_switching[:η_rd][1]
        else
            η_val_master = JuMP.value(x ->JuMP.callback_value.(cb_data, x), JuMP.objective_function(switching_model))
        end

        # pre-presolve heuristics sometimes ignore previous lazy constraints -> re-submit the violated ones
        if PARAMS[:BENDERS][:NESTED][:ADD_IP_CUTS] && !presolved && !isempty(cons)
            resubmitted = false
            for con in cons
                if JuMP.value(x -> JuMP.callback_value.(cb_data, x), con.func) - con.set.upper > SOLVER_PRECISION
                    JuMP.MOI.submit(switching_model, JuMP.MOI.LazyConstraint(cb_data), con)
                    resubmitted = true
                end
            end
            resubmitted && return
        end
        
        # solve until converged, i.e., no more cuts were found (could also use an η tolerance)
        η_val_new = -Inf
        n_iter = 0
        master_var_valss_gens = Vector{Dict{Symbol, AbstractVector}}()

        prec_mult = 1.0
        n_iter_same = 0
        
        stop_sym = ' '
        temp_cuts = Vector{JuMP.ConstraintRef}() # save cuts that are not in the pool, and thus have to be removed from the intermediate problem
        temp_cuts_ids = Vector{Int}() # save indices of cuts that are not in the pool, and thus have to be removed from the intermediate problem
        cut_count = n_pool_cuts
        while true
            # solve pm_master_gens after fixing the master variables from the switching master    
            solve_time = time()
            JuMP.optimize!(gen_model)

            if JuMP.termination_status(gen_model) == JuMP.MOI.INFEASIBLE
                PARAMS[:PRINTING][:BENDERS][:LOG_LINE_INTERPROBLEM] && @printf("%sf  %2d /                       %4.2f /               /           /        %s \n", PAD, n_iter, time() - solve_time, sw_str)
                break
            end
            
            # get objective function value
            η_val_before = η_val_new
            η_val_new = JuMP.objective_value(gen_model)

            # if the solution does not change for a while, try to increase the optimality cut acceptance tolerance
            if isapprox(η_val_before, η_val_new, atol=SOLVER_PRECISION)
                n_iter_same += 1

                if n_iter_same ≥ PARAMS[:BENDERS][:NESTED][:SAME_SOL_THRESH]
                    prec_mult *= 10.0
                    @warn("Same solution for $(PARAMS[:BENDERS][:NESTED][:SAME_SOL_THRESH]) iterations. Increasing precision to check cut violation by x$prec_mult")
                end

            else
                n_iter_same = 0
                prec_mult = 1.0
            end

            PARAMS[:PRINTING][:BENDERS][:LOG_LINE_INTERPROBLEM] && @printf("%s   %2d /                       %4.2f /     %9.2f / %9.2f /        %s \n", PAD, n_iter, time() - solve_time, η_val_new, η_val_master, sw_str)

            master_var_vals_gens::Dict{Symbol, AbstractVector} = Dict(key => JuMP.value.(vars) for (key, vars) in master_vars_gen)
            master_var_vals_var_gens =Dict{VariableRef, Float64}()
            for (_, vars) in master_vars_gen
                for var in vars
                    master_var_vals_var_gens[var] = JuMP.value(var)
                end
            end
            master_var_valss_gens::Vector{Dict{Symbol, AbstractVector}} = [copy(master_var_vals_gens) for _ in 1:n_threads]
            master_var_valss_var_gens::Vector{Dict{VariableRef, Float64}} = [copy(master_var_vals_var_gens) for _ in 1:n_threads]

            if PARAMS[:BENDERS][:NESTED][:PRIMAL_BOUND_STOP] # already over primal bound

                # if a solution has been found by Gurobi and us, check if we are already worse (and can stop generating more cuts to get even worser)
                if sol_count[] > 0 || best_η_val < Inf
                    # `sol_obj_value - η_val_master` is equal to the objective function without benders cuts, so this gives us the obj func value after """replacing""" the tightest old cut by this new one
                    if η_val_new + (sol_obj_value[] - η_val_master) > primal_bound[] || η_val_new + (sol_obj_value[] - η_val_master) > best_η_val
                        stop_sym = 'b'
                        break
                    end
                end
            end

            if ! isempty(pm_subs)
                # solve contingency subproblems
                _add_data_to_subproblem_prio_list!.(subproblem_prio_list, [pm_master_gens], [get(var(pm_master_gens, nw), :p, Dict())], nw=nw)
                
                # add feascut hints to subproblem_prio_list
                add_feascut_hint_to_subproblem_prio_list!(subproblem_prio_list, contingencies_hint_feascut, subproblem_solver_type)

                # solve powerflow subproblems for all contingencies, only use master_vars_gen, since any others cannot be in the cut (the variable does not exist in the gen master)
                (powerflow_conss, violations, _, _) = solve_subproblems_multi(n_threads, pm_subs, cb_data, contingencies, sub_cons_for_dualss_gen, master_varss_sp, master_var_valss_gens, master_var_valss_var_gens, lhsss_sp, dualsss_gen, subproblem_prio_list, solvers, cut_exprss_sp, nw=1, int_sol=int_sol, prec_mult=prec_mult)
                
                # convergence
                if isempty(powerflow_conss) # no new cuts
                    stop_sym = '∅'

                    # sometimes the cuts are not respected and we need to re-submit them :c
                    if n_iter == 0 && PARAMS[:BENDERS][:NESTED][:ADD_SP_CUTS]
                        for con in cons
                            val = JuMP.value(x ->JuMP.callback_value.(cb_data, x), con.func) - con.set.upper
                            if val > 0.0
                                JuMP.MOI.submit(switching_model, JuMP.MOI.LazyConstraint(cb_data), con)
                            end
                        end
                    end

                    break
                elseif PARAMS[:BENDERS][:NESTED][:NO_CHANGE_STOP] && (η_val_new < η_val_before + SOLVER_PRECISION) # no improvement through cut
                    @warn("solution might no longer be optimal")
                    stop_sym = '='
                    break
                end
                
                for (i, powerflow_cons) in enumerate(powerflow_conss)
                    # add constraint to gen_model
                    cut_id = callback_stats[:n_benders_cuts_added_mp_gens]
                    c_ref = JuMP.add_constraint(gen_model, powerflow_cons[1], "benders_cut_$cut_id")
                    
                    # add to sub_cons_for_duals_switching, lhss, dualss
                    if PARAMS[:BENDERS][:NESTED][:ADD_IP_CUTS]
                        sub_cons_for_duals_switching[:benders_cut][cut_id] = c_ref
                        lhss_switching[:benders_cut][cut_id] = get_lhs_benders_cut(c_ref, gen_var_to_master_var_switching)
                        dualss_switching[:benders_cut][cut_id] = 0.0
                    end

                    # if cut pool is already full, save index to delete later on
                    if cut_count >= PARAMS[:BENDERS][:NESTED][:CUT_POOL_SIZE] || violations[i] == Inf
                        push!(temp_cuts, c_ref)
                        push!(temp_cuts_ids, cut_id)
                    else
                        cut_count += 1
                    end

                    # insert values into con
                    callback_stats[:n_benders_cuts_added_mp_gens] += 1

                    # add to switching master
                    if PARAMS[:BENDERS][:NESTED][:ADD_SP_CUTS]
                        push!(cons, powerflow_cons[2])
                        JuMP.MOI.submit(switching_model, JuMP.MOI.LazyConstraint(cb_data), powerflow_cons[2])
                    end
                end
            else
                break
            end
            n_iter += 1
        end

        PARAMS[:PRINTING][:BENDERS][:LOG_LINE_INTERPROBLEM] && @printf("%s%s  %2d /                       %4.2f /     %9.2f / %9.2f /        %s \n", PAD, stop_sym, n_iter, time() - start_time, η_val_new, η_val_master, sw_str)

        if PARAMS[:BENDERS][:NESTED][:ADD_IP_CUTS]
            update_dualss_opt!(dualss_switching, sub_cons_for_duals_switching)

            # check which cut to add (feasibility or optimality)
            if JuMP.termination_status(gen_model) == JuMP.MOI.OPTIMAL
                η_rd = var(pm_master_switching, nw, :η_rd)[1]
            elseif JuMP.termination_status(gen_model) == JuMP.MOI.INFEASIBLE
                η_rd = 0
            else
                @error("Model status is not optimal or infeasible. Status: $(JuMP.termination_status(gen_model))")
            end

            # check if the constraint will be binding in the switching master
            if !(JuMP.termination_status(gen_model) == JuMP.MOI.OPTIMAL && η_val_new <= η_val_master + SOLVER_PRECISION)
                # build the cut using the dual variable values from the gen master problem and add it to the switching master problem
                constant = JuMP.termination_status(gen_model) == JuMP.MOI.OPTIMAL ? JuMP.objective_function(gen_model).constant : 0.0
                update_benders_expression!(cut_expr_switching, constant, lhss_switching, dualss_switching)

                # finish the cut
                con = JuMP.@build_constraint(cut_expr_switching <= η_rd)

                # remember cuts added before presolve, to re-submit them if they are ignored
                if !presolved
                    push!(cons, con)
                end

                JuMP.MOI.submit(switching_model, JuMP.MOI.LazyConstraint(cb_data), con)

                # if optimal, we might have a (better) primal solution by simply updating the value for eta, but Gurobi would disregard the entire solution
                if PARAMS[:BENDERS][:SUBMIT_HEUR_SOLUTION] && JuMP.termination_status(gen_model) == JuMP.MOI.OPTIMAL
                    # total objective value
                    obj_val_new = η_val_new + (sol_obj_value[] - η_val_master)

                    # if there is no incumbent, this solution is obviously good, otherwise check if better than incumbent
                    new_sol = (sol_count[] == 0 && best_η_val ≥ Inf) || obj_val_new < best_η_val

                    if new_sol
                        best_η_val = obj_val_new

                        variables, values = build_benders_solution(master_vars_switching, [η_val_new], cb_data, :η_rd)
                        
                        JuMP.MOI.submit(switching_model, JuMP.MOI.HeuristicSolution(cb_data), variables, values)
                    end
                end
            end
        end

        n_pool_cuts = cut_count
        JuMP.delete.(gen_model, temp_cuts)
        if PARAMS[:BENDERS][:NESTED][:ADD_IP_CUTS]
            filter!(id_con  -> id_con[1]::Int  ∉ temp_cuts_ids, sub_cons_for_duals_switching[:benders_cut])
            filter!(id_lhs  -> id_lhs[1]::Int  ∉ temp_cuts_ids, lhss_switching[:benders_cut])
            filter!(id_dual -> id_dual[1]::Int ∉ temp_cuts_ids, dualss_switching[:benders_cut])
        end
    end

    # add callback to solver
    JuMP.MOI.set(switching_model, Gurobi.CallbackFunction(), nested_benders_callback)
    return pm_master_switching
end