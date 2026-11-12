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

"returns the to var corresponding master var"
function get_master_var_from_var(var, master_vars::Dict{Symbol, AbstractVector})
    c_id = parse(Int64, match(r"\[(\d+)\]", JuMP.name(var)).captures[1])

    return master_vars[get_var_basename(var)][c_id]
end

"gets lhs of subproblem constraints to be multiplied with dual"
function get_lhs(cons_dict::Dict{Int, JuMP.ConstraintRef}, master_vars::Dict{Symbol, AbstractVector})::Dict{Int, JuMP.GenericAffExpr}
    return Dict(comp_id => get_lhs(con, master_vars) for (comp_id, con) in cons_dict)
end
function get_lhs(con::JuMP.ConstraintRef, master_vars::Dict{Symbol, AbstractVector})::JuMP.GenericAffExpr
    # get all variables that are part of the constraint if they are contained in the master variables, i.e., should be used for the cut
    vars = [(var, coeff) for (var, coeff) in JuMP.constraint_object(con).func.terms if get_var_basename(var) ∈ keys(master_vars)]

    # init expr to rhs
    rhs  = JuMP.normalized_rhs(con)

    rounding_digits = PARAMS[:PRECISION][:SOLVER_PREC_DIGITS]
    if rounding_digits >= 0 && PARAMS[:ROUNDING][:MODEL_COEFFS]
        rhs = round(rhs, digits=rounding_digits)
    end

    expr = JuMP.AffExpr(rhs)
    for (var, coeff) in vars
        JuMP.add_to_expression!(expr, -coeff, get_master_var_from_var(var, master_vars))
    end
    
    return expr
end
"gets lhs of subproblem constraints to be multiplied with dual *specific for benders cuts*"
function get_lhs_benders_cut(con::JuMP.ConstraintRef, var_to_master_vars::Dict{VariableRef, VariableRef})::JuMP.GenericAffExpr
    # init expr to rhs
    rhs  = JuMP.normalized_rhs(con)

    rounding_digits = PARAMS[:PRECISION][:SOLVER_PREC_DIGITS]
    if rounding_digits >= 0 && PARAMS[:ROUNDING][:MODEL_COEFFS]
        rhs = round(rhs, digits=rounding_digits)
    end

    expr = JuMP.AffExpr(rhs)
    for (var, coeff) in JuMP.constraint_object(con).func.terms
        if haskey(var_to_master_vars, var)
            JuMP.add_to_expression!(expr, -coeff, var_to_master_vars[var])
        end
    end
    
    return expr
end

"compares pf and vad of `solver`` with those of its fallback solver"
function __check_primal_solution_with_fallback(solver, pm::AbstractPowerModel, c::Int64, cb_data, master_var_vals::Dict{Symbol, AbstractVector}, master_var_vals_var::Dict{VariableRef, Float64}, outage_branches::Vector{Int64}, master_vars, lhss, dualss, cut_exprs::Vector{AffExpr}; nw::Int64=nw_id_default)::Tuple{<:Union{Vector{JuMP.ScalarConstraint}, Nothing}, Float64, Symbol, Type}
    if hasfield(typeof(solver), :fallback_solver)
        old_rhss = apply_fallback_update!(solver.fallback_solver)

        (new_con, eta, type, solver_type) = solve_subproblem(pm, 
                                                                c, 
                                                                cb_data, 
                                                                [master_vars], 
                                                                master_var_vals,
                                                                master_var_vals_var,
                                                                [lhss], 
                                                                dualss, 
                                                                outage_branches, # these vars are excluded from the cut + used to check primal and dual feasibility
                                                                solver.fallback_solver.solver,
                                                                cut_exprs,
                                                                nw=nw, 
                                                                int_sol=true,
                                                                counter_fallback=true)

    else 
        return
    end

    # if FullLP solver returned a feasibility cut, so there is no solution from which we could compare power flows
    if type == :Feasibility
        println("Fallback solver returned a feasibility cut, skipping primal solution check.")
        eta = -1.0
    else

        master_var_vals_branch = master_var_vals[:z_branch]::JuMP.Containers.DenseAxisArray{Float64, 1, Tuple{Vector{Int64}}, Tuple{JuMP.Containers._AxisLookup{Dict{Int64, Int64}}}}
        if PARAMS[:BENDERS][:CHECKS][:CHECK_PF_WITH_FALLBACK]
            wrong_branches = Vector{Int}()
            for (branch_id, branch::Dict{Symbol, Any}) in solver.branch_data
                # only check "normal" branches if in busbar mode. We do not use this information for the cut, so we do not have to check it either.
                if pm.data["BUSBARS"] && ! is_branch_electrical(branch)
                    continue
                end

                f_bus::Int64 = branch[:f_bus]
                t_bus::Int64 = branch[:t_bus]

                # the matrix solver will always return a flow (based on the buses' angles), even if disconnected, off, or outaged. We do not use this information for the cut, so we do not have to check it either.
                if typeof(solver) <: AbstractMatrixSolver 
                    if solver.bus_to_col[f_bus] == 0 || solver.bus_to_col[t_bus] == 0 || branch_id ∈ outage_branches
                        continue
                    end

                    # else
                    is_off = is_branch_switchable(branch) ? (master_var_vals_branch[branch_id]::Float64 < 0.5) : (branch[:status]::Int == 0)
                    if is_off
                        continue
                    end
                end

                # for the current ep solution, calculate vad and pf
                vad1 = get_vad(solver, f_bus, t_bus)
                pf1 = get_pf(solver, branch, vad1)
                vad2, pf2 = get_vad_pf(solver.fallback_solver.solver, branch, nw=nw)
                #pf2 = get_pf(solver,branch,vat)
                
                adjusted_rate_a = get_adjusted_rate_a(pm, branch)

                load1 = pf1/adjusted_rate_a
                load2 = pf2/adjusted_rate_a

                if !isapprox(load1, load2, atol=0.01) #|| !isapprox(pf1, pf2, atol=SOLVER_PRECISION*10000)
                    push!(wrong_branches, branch_id)
                    #@show ref(pm, nw, :branch, branch_id)
                    @printf("Branch %4d: pf  = %10.6f, pf_fb  = %10.6f\n", branch_id, pf1, pf2)
                    @printf("Branch %4d: vad = %10.6f, vad_fb = %10.6f\n", branch_id, vad1, vad2)
                    if ! is_branch_switchable(branch)
                        if is_branch_electrical(branch)
                            #@show(solver.fallback_solver.solver.sub_cons_for_duals[:ohm_fr_ub][branch_id])
                            #@show(solver.fallback_solver.solver.sub_cons_for_duals[:ohm_fr_lb][branch_id])
                            #@show(JuMP.normalized_rhs(solver.fallback_solver.solver.sub_cons_for_duals[:ohm_fr_ub][branch_id]))
                            #@show(JuMP.normalized_rhs(solver.fallback_solver.solver.sub_cons_for_duals[:ohm_fr_lb][branch_id]))
                        end
                        #@show(JuMP.normalized_rhs(solver.fallback_solver.solver.sub_cons_for_duals[:vad_ub][branch_id]))
                        #@show(JuMP.normalized_rhs(solver.fallback_solver.solver.sub_cons_for_duals[:vad_lb][branch_id]))
                    else
                        @show(master_var_vals[:z_branch][branch_id])
                        #@show(master_var_vals[:z_branch_change][branch_id])
                    end
                    #@show(solver.vas[solver.bus_to_col[branch[:f_bus]]], solver.vas[solver.bus_to_col[branch[:t_bus]]])
                end

                # only print all those that now incur additional costs in exactly one of the two solvers
                if (pf1 > adjusted_rate_a + SOLVER_PRECISION && pf2 <= adjusted_rate_a - SOLVER_PRECISION) || (pf1 <= adjusted_rate_a - SOLVER_PRECISION && pf2 > adjusted_rate_a + SOLVER_PRECISION)
                    @warn("[Branch Overload Mismatch] Solver: $(pf1 > adjusted_rate_a), FB solver: $(pf2 > adjusted_rate_a)")
                    @printf("                           Branch %4d: pf  = %10.6f, pf_fb  = %10.6f, rate_a = %10.6f\n", branch_id, pf1, pf2, adjusted_rate_a)
                    #@printf("Branch %4d: vad = %10.6f, vad_fb = %10.6f\n", branch_id, vad1, vad2)
                end

            end

            # if the solver is a busbar solver, check that the components match, i.e. that the fallback va of the components of the solver is the same
            if hasproperty(solver, :components)
                for component in solver.components
                    ref_bus_id = isempty(component.breaker_buses) ? component.busbar_buses[1] : component.breaker_buses[1]
                    ref_va = get_va(solver.fallback_solver.solver, ref_bus_id, nw=nw)

                    # check for differences
                    is_va_different = false
                    for bus_id in component.busbar_buses
                        if !isapprox(get_va(solver.fallback_solver.solver, bus_id, nw=nw), ref_va, atol = SOLVER_PRECISION)
                            is_va_different = true
                            break
                        end
                    end
                    if !is_va_different
                        for bus_id in component.breaker_buses
                            if !isapprox(get_va(solver.fallback_solver.solver, bus_id, nw=nw), ref_va, atol = SOLVER_PRECISION)
                                is_va_different = true
                                break
                            end
                        end
                    end

                    # print violation info
                    if is_va_different
                        @warn("Buses of solver component have differing fallback solver voltage angles:")
                        println("Busbar buses:")
                        for bus_id in component.busbar_buses
                            println("$bus_id: va $(get_va(solver.fallback_solver.solver, bus_id, nw=nw))")
                        end
                        println("Breaker buses:")
                        for bus_id in component.breaker_buses
                            println("$bus_id: va $(get_va(solver.fallback_solver.solver, bus_id, nw=nw))")
                        end
                        println("Coupler branches:")
                        for branch_id in component.coupler_branches
                            branch = ref(pm, nw, :branch, branch_id)
                            status = round(master_var_vals[:z_branch][branch_id])
                            println("$branch_id: $(branch["f_bus"]) -> $(branch["t_bus"]), status $status")
                        end
                        println("Breaker branches:")
                        for branch_id in component.breaker_branches
                            branch = ref(pm, nw, :branch, branch_id)
                            status = round(master_var_vals[:z_branch][branch_id])
                            println("$branch_id: $(branch["f_bus"]) -> $(branch["t_bus"]), status $status")
                        end
                    end
                end
            end

            if length(wrong_branches) > 10 
                throw(InvalidStateException("Primal solution of solver does not match 10 or more of the fallback solver for branches: $wrong_branches", :primal_solution_does_not_match_fallback))
            end
        end
    end
    
    reset_fallback_update!(solver.fallback_solver, old_rhss)
    return new_con, eta, type, solver_type
end

"checks the current primal solution computed by `solver`, 
 returns (true iff a new cut should exist; true iff only ratea constraints where violated; the objective function value (-1 if not only ratea violations); betas, a Dict containing violating branches for optimality/feasibility cuts; whether a dual solution can be immediately derived from betas (currently only for certain feas cuts))"
function check_primal_solution(solver::Union{AbstractMatrixSolver, LPSolver}, pm::AbstractPowerModel, c::Int64, master_var_vals::Dict{Symbol, AbstractVector}, outage_branches::Vector{Int64}; nw::Int64=nw_id_default, print_costs::Bool=false)::Tuple{Bool,Bool,Float64,Dict{Int64, Float64}}
    ## either 
    ## (1) the easy primal (ep) is infeasible         -> compute duals & (feasibility?) cuts
    ## (2) the ep violates switched off branches      -> :(
    ## (3) the ep violates only [ratea]s              -> the objective function value of the ep is more/less than the master computed for it:
    ##   (3a) more                                       -> compute duals & optimality cuts
    ##   (3b) less                                       -> no new cuts, continue
    ## (4) the ep violates at least one vad 
    ## (4.1) the ep violates at least one rate_a*big_m_mult
    ##     (or maybe other constraints)               -> compute duals & feasibility cuts
    ## (5) vad violations of off branches             -> compute duals & feasibility cuts
    ## (6) no violations                              -> no new cuts, continue
    ##

    # (not (1), otherwise the method would not have been called)

    has_new_cut             = false
    only_ratea_violations   = true
    obj_func_value          = 0.0

    betas = Dict{Int64, Float64}() # duals corresponding to rating_from_slack_lb/ub: here it is -[the objective function coefficient] of all overloaded branches, else 0
    vad_violation = Dict{Int64, Float64}() # first branch with a vad violation (4)/(5), stored in case of marginal violation

    master_var_vals_branch = master_var_vals[:z_branch]::JuMP.Containers.DenseAxisArray{Float64, 1, Tuple{Vector{Int64}}, Tuple{JuMP.Containers._AxisLookup{Dict{Int64, Int64}}}}

    apply_ub_mults=(PARAMS[:LIMIT_MULT]::Dict{Symbol, Any})[:APPLY_UB_MULTS]::Bool
    for_n1 = get(pm.data, "SUBPROBLEM", false)::Bool
    for (branch_id, branch::Dict{Symbol, Any}) in solver.branch_data
        f_bus::Int64 = branch[:f_bus]
        t_bus::Int64 = branch[:t_bus]

        is_off = is_branch_switchable(branch) ? (master_var_vals_branch[branch_id]::Float64 < 0.5) : (branch[:status]::Int == 0)

        # the matrix solver will always return a flow (based on the buses' angles), even if disconnected, off, or outaged. We do not use this information for the cut, so we do not have to check it either. However, we need to check vads for off branches
        if typeof(solver) <: AbstractMatrixSolver 
            if solver.bus_to_col[f_bus] ≤ 0 || solver.bus_to_col[t_bus] ≤ 0 || branch_id ∈ outage_branches
                continue
            end

            # else
            if is_off

                if PARAMS[:BIG_M_MULT][:SET_MAX_ANGLE_CLOSE]
                    # (5) check all off branches for vad violations - also non-electrical ones
                    vad = get_vad(solver, f_bus, t_bus)

                    angmin = branch[:angmin_close]::Float64
                    angmax = branch[:angmax_close]::Float64

                    # check (5): vad violations of off branches: check whether not angmin ≤ vad ≤ angmax
                    if vad < angmin - 5*SOLVER_PRECISION || angmax + 5*SOLVER_PRECISION < vad
                        return (true, false, -1.0, Dict{Int64, Float64}(branch_id => vad))
                    elseif (vad < angmin - SOLVER_PRECISION || angmax + SOLVER_PRECISION < vad) && isempty(vad_violation)
                        vad_violation[branch_id] = vad # marginal -> look for better cut
                    end
                end

                continue
            end
        end

        # Active busbar couplers are represented by equal angles in the matrix.
        # Off couplers still have angle bounds and were checked above.
        is_branch_electrical(branch) || continue

        # calculate vad and pf
        vad = get_vad(solver, f_bus, t_bus)
        pf = get_pf(solver, branch, vad)
        rate_a = get_adjusted_rate_a(branch, for_n1, apply_ub_mults=apply_ub_mults)::Float64

        # check (2): if the branch is off, pf needs to be 0
        if typeof(solver) == :LPSolver
            
            if is_off || branch_id ∈ outage_branches
                # check (2)
                if _isapproxgreater(abs(pf), 0, atol=SOLVER_PRECISION)
                    has_new_cut           = true
                    only_ratea_violations = false
                    throw(InvalidStateException("Easy primal solution violates limits of switched-off branches for branch $branch_id with power flow $pf.", :easy_primal_violates_switched_off_branches))
                end
                continue # there cannot be any rate a or angle violations for switched-off branches, so continue
            end
        end
        angmin = branch[:angmin]::Float64
        angmax = branch[:angmax]::Float64

        # check (4): vad violations: check whether not angmin ≤ vad ≤ angmax
        if vad < angmin - SOLVER_PRECISION || angmax + SOLVER_PRECISION < vad
            if print_costs
                @printf("Angle violation: %4d with vad %8.4f not in [%8.4f, %8.4f]\n", branch_id, vad, angmin, angmax)
            end
            if vad < angmin - 5*SOLVER_PRECISION || angmax + 5*SOLVER_PRECISION < vad
                return (true, false, -1.0, Dict{Int64, Float64}(branch_id => vad)) # we need a feasibility cut, so one violating branch is sufficient
            end
            isempty(vad_violation) && (vad_violation[branch_id] = vad) # marginal -> look for better cut
        end

        # check (4.1): rate a violations that cannot be covered with slack (see constraint_thermal_limit_*_slack_on_off)
        #if (rate_a * PARAMS[:BIG_M_MULT][:RATE_A])-abs(pf) < -SOLVER_PRECISION && branch[:monitored]
        #    @warn("Big M multiplier for rate_a is not large enough for branch $branch_id with power flow $pf for some topology.")
        #    return (true, false, 0.0, betas)
        #end
        
        # collect (3): rate_a violations
        if rate_a-abs(pf) < -SOLVER_PRECISION::Float64 && branch[:monitored]
            has_new_cut      = true
            obj_coeff = branch[:obj_coeff]
            obj_func_value  += (abs(pf) - rate_a) * obj_coeff # objective_min_slack_*_*
            #@info("Overloaded branch: $branch_id with power flow $pf > rate_a $rate_a. Cost: $((abs(pf) - rate_a) * obj_coeff)")
            betas[branch_id] = (pf > 0 ? -1 : 1) * obj_coeff # we need an optimality cut, so we start collecting the value of it

            if print_costs
                @printf("Overloaded branch: %4d with power flow %12.4f > rate_a %12.4f. Cost: %12.4f\n", branch_id, abs(pf), rate_a, (abs(pf) - rate_a) * obj_coeff)
            end
        end
    end

    # skip if no violations (6), else (3a) or (3b)
    if has_new_cut && only_ratea_violations
        # check (3b)
        if _isapproxgreatereq(master_var_vals[:η][c], obj_func_value, atol=SOLVER_PRECISION)
            has_new_cut = false
        end
    end

    # marginal (4)/(5) violation (at most 5x the tolerance) might yield numerically unstable feas cuts -> only choose them if no better cut was found
    if ! isempty(vad_violation) && ! has_new_cut
        return (true, false, -1.0, vad_violation)
    end

    return (has_new_cut, only_ratea_violations, obj_func_value, betas)
end

"checks the current primal solution computed by `solver`, 
 returns (true iff a new cut should exist; true iff only ratea constraints where violated; the objective function value (-1 if not only ratea violations); nothing; true, indicating that the dual problem has already been solved)"
function check_primal_solution(solver::FullLPSolver, pm::AbstractPowerModel, c::Int64, master_var_vals::Dict{Symbol, AbstractVector}, outage_branches::Vector{Int64}; nw::Int64=nw_id_default)::Tuple{Bool,Bool,Float64,Dict{Int64, Float64}}
    @assert JuMP.termination_status(solver.pm.model) ≠ JuMP.MOI.INFEASIBLE_OR_UNBOUNDED

    # print power flow for this subproblem solution
    #for (branch_id, branch) ∈ ref(pm, nw, :branch)
        #pf = JuMP.value(var(solver.pm, nw, :p, (branch["index"], branch["f_bus"], branch["t_bus"])))
        #@printf("Branch %4d: %12.4f / %12.4f (%6.2f %%) %s %s\n", branch_id, abs(pf), get_adjusted_rate_a(pm, branch), abs(pf) / get_adjusted_rate_a(pm, branch) * 100, is_branch_monitored(pm, nw, branch_id) ? "m" : " ", abs(pf) > get_adjusted_rate_a(pm, branch) && is_branch_monitored(pm, nw, branch_id) ? " --- OVERLOADED ---" : "")
    #    pf = JuMP.value(var(solver.pm, nw, :p, (branch["index"], branch["f_bus"], branch["t_bus"])))
    #    rate_a_slack = JuMP.value(var(solver.pm, nw, :rate_a_slack, (branch["index"], branch["f_bus"], branch["t_bus"])))
    #    if rate_a_slack > 0
    #        @printf("Branch %4d (rate_a_slack=%12.4f): %12.4f / %12.4f (%6.2f %%) %s %s\n", branch_id, rate_a_slack, abs(pf), get_adjusted_rate_a(pm, branch), abs(pf) / get_adjusted_rate_a(pm, branch) * 100, is_branch_monitored(pm, nw, branch_id) ? "m" : " ", abs(pf) > rate_a_slack && is_branch_monitored(pm, nw, branch_id) ? " --- OVERLOADED ---" : "")
    #        @show(branch)
    #    end
    #end

    return JuMP.termination_status(solver.pm.model) == JuMP.MOI.INFEASIBLE || (JuMP.termination_status(solver.pm.model) == JuMP.MOI.OPTIMAL && JuMP.objective_value(solver.pm.model) > 0),
           JuMP.termination_status(solver.pm.model) == JuMP.MOI.OPTIMAL,
           JuMP.termination_status(solver.pm.model) == JuMP.MOI.OPTIMAL ? JuMP.objective_value(solver.pm.model) : -1,
           Dict{Int64, Float64}()
end

"checks the current primal solution computed by `solver`, 
 returns (true iff a new cut should exist; true iff only ratea constraints where violated; the objective function value (-1 if not only ratea violations); betas, a Dict containing violating branches for optimality/feasibility cuts; whether a dual solution can be immediately derived from betas (currently only for certain feas cuts))"
function check_primal_solution(solver::FallbackSolver, pm::AbstractPowerModel, c::Int64, master_var_vals::Dict{Symbol, AbstractVector}, outage_branches::Vector{Int64}; nw::Int64=nw_id_default)::Tuple{Bool,Bool,Float64,Dict{Int64, Float64}}
    return check_primal_solution(solver.solver, pm, c, master_var_vals, outage_branches, nw=nw)
end

"prints the range of the coefficients of the constraint function and the right-hand side"
function print_coeff_range(c::JuMP.ConstraintRef, con::JuMP.ScalarConstraint)
    # collect coeffs
    coeffs = Vector{Float64}()

    func = con.func
    for (var, coeff) in func.terms
        push!(coeffs, abs(coeff))
    end

    @printf("Coeff range: [%1.0e, %1.0e], RHS: %3.1e\n", minimum(coeffs), maximum(coeffs), JuMP.normalized_rhs(c))
end

"prints the range of the coefficients of the constraint function and the right-hand side"
function print_coeff_range(con::JuMP.ScalarConstraint)
    # collect coeffs
    coeffs = Vector{Float64}()

    func = con.func
    for (var, coeff) in func.terms
        push!(coeffs, abs(coeff))
    end

    @printf("Coeff range: [%1.0e, %1.0e]\n", minimum(coeffs), maximum(coeffs))
end

"add branching priority to model"
function add_branching_priorities(pm::AbstractPowerModel; nw::Int64=nw_id_default)
    MOIU.attach_optimizer(pm.model)
    for (b_id, branch) in ref(pm, nw, :branch)
        if is_branch_switchable(branch) && pm.data["WARMSTART_PROVIDED"]
            # for each branch retrieve the variable from the model
            z = var(pm, nw, :z_branch, b_id)
            # and set the priority
            MOI.set(JuMP.backend(pm.model).optimizer, Gurobi.VariableAttribute("BranchPriority"), z.index, branch["sw_heur_priority"])
        end
    end
end

"compiles a dictionary that contains the variables for easier access"
function get_master_vars(pm::AbstractPowerModel; nw::Int64=nw_id_default)::Dict{Symbol, AbstractVector}
    master_vars     = Dict{Symbol, AbstractVector}(:pg         => var(pm, nw, :pg), 
                                                   :η          => var(pm, nw, :η))
    if PARAMS[:LIMITS][:GENERATOR][:GEN_SWITCHING] 
        push!(master_vars, :z_gen => var(pm, nw, :z_gen))
        # what about z_gen_startup and z_gen_shutdown?
    end
    if PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_SHIFTING]
        push!(master_vars, :tr_shift => var(pm, nw, :tr_shift))
        if haskey(var(pm, nw), :tr_shift_change)
            push!(master_vars, :tr_shift_change => var(pm, nw, :tr_shift_change))
        end
    end
    if PARAMS[:LIMITS][:BRANCH][:BRANCH_SWITCHING]
        push!(master_vars, :z_branch => var(pm, nw, :z_branch))
        push!(master_vars, :z_branch_change => var(pm, nw, :z_branch_change))
    end
    return master_vars
end

"compiles a dictionary that contains the variables for easier access, for nested Benders decomposition"
function get_master_vars_nested(pm_master_switching::AbstractPowerModel, pm_master_gens::AbstractPowerModel; nw::Int64=nw_id_default)::Tuple{Dict{Symbol, AbstractVector}, Dict{Symbol, AbstractVector}}
    master_vars_switching = Dict{Symbol, AbstractVector}(:η_rd => var(pm_master_switching, nw, :η_rd))
    master_vars_gen      = Dict{Symbol, AbstractVector}(:pg => var(pm_master_gens, nw, :pg),
                                                         :η => var(pm_master_gens, nw, :η))
    if PARAMS[:BENDERS][:NESTED][:ADD_SP_CUTS]
        push!(master_vars_switching, :pg => var(pm_master_switching, nw, :pg))
        push!(master_vars_switching, :η  => var(pm_master_switching, nw, :η))
    end

    if PARAMS[:LIMITS][:GENERATOR][:GEN_SWITCHING] 
        push!(master_vars_switching, :z_gen => var(pm_master_switching, nw, :z_gen))
        push!(master_vars_gen, :z_gen => var(pm_master_gens, nw, :z_gen))
    end
    if PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_SHIFTING]
        push!(master_vars_switching, :tr_shift => var(pm_master_switching, nw, :tr_shift))
        push!(master_vars_gen, :tr_shift => var(pm_master_gens, nw, :tr_shift))
    end
    if PARAMS[:LIMITS][:BRANCH][:BRANCH_SWITCHING]
        push!(master_vars_switching, :z_branch => var(pm_master_switching, nw, :z_branch))
        push!(master_vars_switching, :z_branch_change => var(pm_master_switching, nw, :z_branch_change))
        push!(master_vars_gen, :z_branch => var(pm_master_gens, nw, :z_branch))
    end
    return master_vars_switching, master_vars_gen
end