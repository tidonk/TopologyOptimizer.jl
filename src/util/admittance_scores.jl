"Base-network self-transfer factors using physical DC susceptance, component by component."
function admittance_self_transfer(nw)
    buses = sort(parse.(Int, collect(keys(nw["bus"]))))
    index = Dict(b => k for (k, b) in enumerate(buses))
    n = length(buses)
    lap = zeros(n, n)
    adjacency = [Int[] for _ in 1:n]
    for b in values(nw["branch"])
        b["status"] == 0 && continue
        f, t = index[b["f_bus"]], index[b["t_bus"]]
        y = b["br_x"] / (b["br_r"]^2 + b["br_x"]^2)
        isfinite(y) || error("Nonfinite DC susceptance")
        iszero(y) && continue
        lap[f, f] += y; lap[t, t] += y
        lap[f, t] -= y; lap[t, f] -= y
        push!(adjacency[f], t); push!(adjacency[t], f)
    end
    inverse = zeros(n, n)
    visited = falses(n)
    for start in 1:n
        visited[start] && continue
        component = [start]
        visited[start] = true
        cursor = 1
        while cursor <= length(component)
            for j in adjacency[component[cursor]]
                if !visited[j]
                    visited[j] = true
                    push!(component, j)
                end
            end
            cursor += 1
        end
        sort!(component)
        reduced = component[2:end]
        isempty(reduced) && continue
        inverse[reduced, reduced] = lap[reduced, reduced] \ Matrix{Float64}(LinearAlgebra.I, length(reduced), length(reduced))
    end
    return Dict(parse(Int, key) => begin
        f, t = index[b["f_bus"]], index[b["t_bus"]]
        y = b["status"] == 0 ? 0.0 : b["br_x"] / (b["br_r"]^2 + b["br_x"]^2)
        y * (inverse[f, f] + inverse[t, t] - inverse[f, t] - inverse[t, f])
    end for (key, b) in nw["branch"])
end

"Self-transfer factors from the existing physical power-flow matrix and grounding."
function admittance_self_transfer(pm, eligible; nw=1, topology=nothing)
    isempty(eligible) && return Dict{Int,Float64}()
    branches = ref(pm, nw, :branch)
    original_status = Dict(i => b["status"] for (i, b) in branches)
    # Only the matrix construction sees the accepted topology. Switching costs
    # and action limits must continue to use the original branch statuses.
    A, _, bus_to_col = try
        if topology !== nothing
            for (i, b) in branches
                b["status"] = topology[i]
            end
        end
        build_base_pf_matrix(pm; nw=nw, only_active=true)
    finally
        for (i, b) in branches
            b["status"] = original_status[i]
        end
    end
    factorization = lu(grounded_pf_matrix(A))
    rhs = zeros(size(A, 1))
    factors = Dict(i => 0.0 for i in eligible)
    for i in eligible
        branch = ref(pm, nw, :branch, i)
        (topology === nothing ? branch["status"] : topology[i]) == 0 && continue
        f, t = bus_to_col[branch["f_bus"]], bus_to_col[branch["t_bus"]]
        fill!(rhs, 0.0)
        rhs[f] += 1.0
        rhs[t] -= 1.0
        angles = factorization \ rhs
        factors[i] = -branch["b"] * (angles[f] - angles[t])
    end
    return factors
end

"Rank initially open eligible branches first, then normalized corrected admittance scores."
function admittance_branch_order(nw, scores, eligible; initially_open_first=true)
    ids = sort(collect(eligible))
    all(isfinite(scores[e]) for e in ids) || error("Nonfinite admittance score")
    scale = max(1e-100, maximum((abs(scores[e]) for e in ids); init=0.0))
    opened = [e for e in ids if initially_open_first && nw["branch"][string(e)]["status"] == 0]
    closed = [e for e in ids if !initially_open_first || nw["branch"][string(e)]["status"] != 0]
    # Julia orders -0.0 before +0.0, whereas both represent the same score.
    # Canonicalize rounded zeros so all such ties use the branch ID.
    rounded = Dict(e => begin
        q = round(scores[e] / scale; digits=9)
        iszero(q) ? 0.0 : q
    end for e in closed)
    sort!(closed; by=e -> (-rounded[e], e))
    return vcat(opened, closed)
end

"Keep Benders cuts and their Ohm-law multipliers for SCOPF dual recovery."
function record_admittance_cut!(pm_sub, contingency, cut, duals)
    multipliers = Dict(i => get(duals[:ohm_fr_on_off_lb], i, 0.0) + get(duals[:ohm_fr_on_off_ub], i, 0.0)
        for i in union(keys(duals[:ohm_fr_on_off_lb]), keys(duals[:ohm_fr_on_off_ub])))
    push!(pm_sub.ext[:admittance_cuts], (contingency, ScalarConstraint(copy(cut.func), cut.set), multipliers))
end

"Save fallback LP results before the existing contingency reset invalidates them."
function record_admittance_lp!(solver; nw=1)
    recovery = solver.pm.ext[:admittance_recovery]
    recovery["lp_seconds"] += solve_time(solver.pm.model)
    recovery["lp_solves"] += 1
    recovery["lp_status"] = termination_status(solver.pm.model)
    recovery["lp_objective"] = has_values(solver.pm.model) ? objective_value(solver.pm.model) : Inf
    recovery["lp_flows"] = has_values(solver.pm.model) ?
        Dict(i => get_vad_pf(solver, ref(solver.pm, nw, :branch, i); nw=nw)[2] for i in recovery["eligible"]) :
        Dict{Int,Float64}()
end

"Recover SCOPF sensitivities from the Benders master and its existing subproblem."
function collect_admittance_scores!(pm; nw=1, topology=nothing)
    @assert length(pm.data["nw"]) == 1 "Admittance scoring currently supports a single timestep"
    @assert PARAMS[:BENDERS][:BASE_PF_IN_MASTER]
    pm_sub = pm.ext[:admittance_subproblem]
    sub_cons = get_sub_cons_for_duals(pm_sub)
    solver = pm.ext[:admittance_solver]
    nw_data = pm.data["nw"][string(nw)]
    active(i) = (topology === nothing ? nw_data["branch"][string(i)]["status"] : topology[i]) != 0
    eligible = sort([b["index"] for b in values(nw_data["branch"]) if is_branch_switchable(b)])
    @assert all(is_branch_electrical(nw_data["branch"][string(i)]) for i in eligible)
    master_vars = get_master_vars(pm; nw=nw)
    scopf_objective = objective_value(pm.model)
    scopf_vals = Dict{Symbol,AbstractVector}(key => value.(vars) for (key, vars) in master_vars)
    scopf_flows = Dict(i => value(var(pm, nw, :p)[(i, nw_data["branch"][string(i)]["f_bus"],
        nw_data["branch"][string(i)]["t_bus"])]) for i in eligible)

    lhss, duals = get_empty_lhss_dualss(sub_cons, master_vars)
    records = pm_sub.ext[:admittance_cuts]
    cut_refs = ConstraintRef[]

    # Fix the discrete incumbent and recover duals on the small master LP.
    # Lazy cuts must be added explicitly before removing the callback.
    fix_discrete_variables(v -> round(value(v)), pm.model)
    MOI.set(pm.model, Gurobi.CallbackFunction(), nothing)
    for (_, cut, _) in records
        push!(cut_refs, add_constraint(pm.model, cut))
    end
    set_silent(pm.model)
    master_seconds = 0.0
    recovery = Dict{String,Any}("eligible" => eligible, "lp_seconds" => 0.0, "lp_solves" => 0)
    pm_sub.ext[:admittance_recovery] = recovery
    leq_solves = 0
    cut_exprs = [AffExpr(0.0)]
    reset_cut_exprs!(cut_exprs, master_vars)
    contingencies = nw_data["contingencies"]
    flows = [Dict{Int,Float64}() for _ in contingencies]
    objectives = zeros(length(contingencies))
    converged = false
    rounds = 0
    reuse_scopf = false
    reuse_scopf_allowed = true
    flow_contingencies = Int[]
    try
        for iteration in 1:100
            rounds = iteration
            optimize!(pm.model)
            master_seconds += solve_time(pm.model)
            @assert termination_status(pm.model) == MOI.OPTIMAL && has_duals(pm.model) "SCOPF master dual recovery failed"
            master_vals = Dict{Symbol, AbstractVector}(key => value.(vars) for (key, vars) in master_vars)
            tolerance = max(1e-4, 1e-7 * abs(objective_value(pm.model)))
            # Any primal optimum can be paired with any dual optimum of an LP.
            # If the original feasible SCOPF point attains this master bound,
            # retain it; a different LP dispatch need not be checked or repaired.
            reuse_scopf = reuse_scopf_allowed && abs(objective_value(pm.model) - scopf_objective) <= tolerance
            if reuse_scopf
                master_vals = deepcopy(scopf_vals)
                flow_contingencies = sort!(unique([c for (cut_ref, (c, _, multipliers)) in zip(cut_refs, records)
                    if !iszero(dual(cut_ref)) && any(i in eligible && !iszero(mu) for (i, mu) in multipliers)]))
                objectives .= [scopf_vals[:η][c] for c in eachindex(contingencies)]
            else
                # A nonzero gap needs the existing subproblem checks and cuts.
                flow_contingencies = collect(eachindex(contingencies))
            end
            previous_count = length(records)
            master_vals_var = Dict{VariableRef,Float64}(v => val for (key, vars) in master_vars
                for (v, val) in zip(vars, master_vals[key]))
            # Reuse the SCOPF solver, its contingency updates, checks, duals and
            # fallback chain. The regular LEQ path does not solve any LP here.
            apply_master!(solver, pm_sub, master_vals; nw=nw)
            try
                for c in flow_contingencies
                    contingency = contingencies[c]
                    outage = get_outage_branches(pm_sub, contingency, master_vals; filter=true, nw=nw)
                    all_outages = get_outage_branches(pm_sub, contingency, master_vals; nw=nw)
                    if isempty(all_outages)
                        pm_sub.data["SUBPROBLEM"] = false
                        apply_adjusted_rate_a!(pm_sub, sub_cons, [lhss]; for_n1=false, nw=nw)
                    end
                    switch_off_outage_branches!(master_vals, outage)
                    solver_outages = solver isa AbstractMatrixSolver ? outage : all_outages
                    old_rhs = apply_contingency!(solver, pm_sub, solver_outages, master_vals, true; nw=nw)
                    if solver isa AbstractMatrixSolver
                        apply_contingency!(solver.fallback_solver, pm_sub, all_outages, master_vals, true; nw=nw)
                    end
                    try
                        leq_solves += solver isa AbstractMatrixSolver
                        cuts, _, cut_type, used_solver = solve_subproblem(pm_sub, c, nothing,
                            [master_vars], master_vals, master_vals_var, [lhss], duals, all_outages, solver, cut_exprs; nw=nw)
                        @assert cut_type in (:Nothing, :Optimality, :Feasibility) "SCOPF recourse recovery failed"
                        if cut_type in (:Optimality, :Feasibility)
                            record_admittance_cut!(pm_sub, c, first(cuts), duals)
                        end
                        if used_solver == FullLPSolver
                            @assert recovery["lp_status"] in (MOI.OPTIMAL, MOI.INFEASIBLE) "SCOPF fallback failed"
                            objectives[c] = recovery["lp_objective"]
                            flows[c] = recovery["lp_flows"]
                        else
                            _, _, objectives[c], _ = check_primal_solution(solver, pm_sub, c, master_vals, all_outages; nw=nw)
                            flows[c] = Dict(i => get_vad_pf(solver, ref(pm_sub, nw, :branch, i); nw=nw)[2]
                                for i in eligible if active(i) && !(i in all_outages))
                        end
                        if cut_type == :Nothing
                            @assert abs(objectives[c] - master_vals[:η][c]) <= tolerance "Benders recourse does not match the master"
                        end
                    finally
                        reset_contingency!(solver, pm_sub, solver_outages, old_rhs; nw=nw)
                        switch_on_outage_branches!(master_vals, outage)
                        if isempty(all_outages)
                            pm_sub.data["SUBPROBLEM"] = true
                            apply_adjusted_rate_a!(pm_sub, sub_cons, [lhss]; for_n1=true, nw=nw)
                        end
                    end
                end
            finally
                reset_master!(solver, pm_sub, master_vals; nw=nw)
            end
            if length(records) == previous_count
                converged = true
                break
            end
            # New cuts invalidate the saved primal point, even when its
            # objective remains within tolerance of the current LP bound.
            # Continue from the reoptimized LP instead of checking the
            # same saved point and adding the same cuts repeatedly.
            reuse_scopf && (reuse_scopf_allowed = false)
            for (_, cut, _) in records[previous_count+1:end]
                push!(cut_refs, add_constraint(pm.model, cut))
            end
        end
    finally
        delete!(pm_sub.ext, :admittance_recovery)
    end
    @assert converged "SCOPF dual recovery did not converge"

    scores = Dict(i => 0.0 for i in eligible)
    for i in eligible
        branch = nw_data["branch"][string(i)]
        active(i) || continue
        pf = reuse_scopf ? scopf_flows[i] : value(var(pm, nw, :p)[(i, branch["f_bus"], branch["t_bus"])])
        for side in ("lb", "ub")
            con = constraint_by_name(pm.model, "$(nw)_ohm_fr_on_off_$(side)_$(i)")
            con === nothing || (scores[i] += pf * dual(con))
        end
    end
    for (cut_ref, (c, _, multipliers)) in zip(cut_refs, records)
        weight = -dual(cut_ref)
        for (i, pf) in flows[c]
            scores[i] += weight * pf * get(multipliers, i, 0.0)
        end
    end
    transfer = admittance_self_transfer(pm, [i for i in eligible if !iszero(scores[i])]; nw=nw, topology=topology)
    corrected = Dict(i => iszero(scores[i]) ? 0.0 :
        (abs(1 - transfer[i]) > 1e-9 ? scores[i] / (1 - transfer[i]) : 0.0) for i in eligible)
    nw_data["admittance_secant_scores"] = corrected
    return Dict{String,Any}("score_solve_seconds" => recovery["lp_seconds"] + master_seconds,
        "score_lp_solves" => recovery["lp_solves"], "score_lp_seconds" => recovery["lp_seconds"],
        "score_leq_solves" => leq_solves, "score_subproblem_solver" => string(typeof(solver)),
        "dual_recovery_master_seconds" => master_seconds, "dual_recovery_rounds" => rounds,
        "scopf_point_reused" => reuse_scopf, "score_contingencies" => flow_contingencies,
        "scopf_point_reuse_disabled" => !reuse_scopf_allowed,
        "master_dual_gap" => scopf_objective - objective_value(pm.model),
        "dual_recovery_objective" => objective_value(pm.model), "dual_recovery_cuts" => length(records),
        "score_source" => "Benders master and weighted cut duals",
        "admittance_scores" => scores, "secant_scores" => corrected,
        "subproblem_objectives" => objectives, "eligible" => eligible)
end

"Return the admittance secant scores computed from the SCOPF warmstart."
function calc_br_sw_heur_scores_admittance_secant(nw_data; only_switchable=true, kwargs...)
    @assert haskey(nw_data, "admittance_secant_scores") "ADMITTANCE_SECANT requires a successful Benders SCOPF warmstart"
    scores = nw_data["admittance_secant_scores"]
    return Dict{Int64,Float64}(b["index"] => get(scores, b["index"], 0.0)
        for b in values(nw_data["branch"]) if !only_switchable || is_branch_switchable(b))
end
