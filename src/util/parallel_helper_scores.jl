# Ported from the archived 20261008 warm-start grid selection.jl.
# Revalidated against current main; scoring uses a compact fixed-topology LP.
parallel_eligible_ids(data) = sort([b["index"] for b in values(data["nw"]["1"]["branch"]) if is_branch_switchable(b)])

function parallel_helper_objective(data)
    weighted_obj_func_factory(
        [objective_min_num_gen_switch, objective_min_num_tr_switch,
         objective_min_num_branch_switch, objective_min_redispatch_cost,
         objective_min_num_redispatch, objective_min_slack_relative_linear,
         objective_min_num_slacks],
        [data["cost_coeffs"][key] for key in
         ["gen_sw", "tr_sw", "br_sw", "redispatch", "num_redispatch", "overload", "num_overload"]])
end

function parallel_fixed_topology_duals(data, solver; perturbation_check=false,
                                      capture=nothing, time_limit=nothing)
    score_data = deepcopy(data)
    score_data["application"] = "DC_SCOTS_TIME_SLACK"
    score_data["FIX_TS"] = false
    score_data["nw"]["1"]["FIX_TS"] = false
    captured = Ref{Any}()
    fixes = Dict{Int,VariableRef}()
    function prepare(pm)
        # Full compact model, not an initial Benders master lacking cuts.
        @assert length(ref(pm, 1, :contingencies)) == length(data["nw"]["1"]["contingencies"])
        relax_integrality(pm.model)
        change = var(pm, 1, :z_branch_change)
        status = var(pm, 1, :z_branch)
        # Deleting a constraint invalidates the model's name lookup cache.
        # Resolve all references once before removing redundant constraints.
        link_refs = Dict((i, side) => constraint_by_name(pm.model,
            "1_br_change_linking_$(side)_$(i)")
            for i in parallel_eligible_ids(data) for side in ("ub", "lb"))
        for i in parallel_eligible_ids(data)
            @assert !is_fixed(change[i])
            fix(change[i], 0.0; force=true)
            fixes[i] = change[i]
            # x=0 and z=initial status imply this bound already. Remove only
            # the duplicate active bound, so it cannot absorb the fixing dual.
            # This preserves the fixed LP and permits inward x>0 perturbations.
            initial = data["nw"]["1"]["branch"][string(i)]["status"]
            if initial == 1 && has_upper_bound(status[i])
                delete_upper_bound(status[i])
            elseif initial == 0 && has_lower_bound(status[i])
                delete_lower_bound(status[i])
            end
            branch = data["nw"]["1"]["branch"][string(i)]
            @assert get(branch,"status_previous",initial) == initial
            # With a single timestep, x=|initial-z| is already encoded by the
            # exact indicator equality. These two inequalities only repeat it
            # or impose x>=0, another redundant bound at the fixing point.
            for side in ("ub","lb")
                link = link_refs[(i, side)]
                @assert link !== nothing
                delete(pm.model,link)
            end
        end
        captured[] = pm
        time_limit === nothing || set_time_limit_sec(pm.model, time_limit)
        return pm
    end
    result = solve_secure_ots_time(score_data, DCPPowerModel, solver,
        parallel_helper_objective(score_data); solver_specific_build_function=prepare, slack_types=[:ratea])
    pm = captured[]
    capture === nothing || (capture[] = (model=pm, result=result))
    time_limit === nothing || @assert time_limit_sec(pm.model) == time_limit
    if termination_status(pm.model) != MOI.OPTIMAL || !has_duals(pm.model)
        return nothing, Dict{String,Any}("usable"=>false,
            "status"=>string(termination_status(pm.model)),
            "reason"=>"Fixing duals require an optimal LP with duals")
    end
    raw = Dict(i => dual(FixRef(v)) for (i,v) in fixes)
    # For min f with x=rhs, dual is df*/d(rhs). Negative means improvement
    # when a no-change variable is allowed to increase towards switching.
    scores = Dict(i => -d for (i,d) in raw)
    diagnostic = Dict{String,Any}("usable"=>true,"status"=>string(termination_status(pm.model)),
        "objective"=>objective_value(pm.model), "fixing_duals"=>raw,
        "scores"=>scores, "n_contingencies"=>length(ref(pm, 1, :contingencies)),
        "rule"=>"rank descending -dual(z_branch_change == 0)",
        "positive_scores"=>count(>(1e-8), values(scores)))
    if perturbation_check
        f0 = objective_value(pm.model)
        checks = Any[]
        for i in sort(collect(keys(scores)); by=i -> (-scores[i], i))[1:min(3,length(scores))]
            eps = 1e-6
            fix(fixes[i], eps; force=true)
            optimize!(pm.model)
            push!(checks, Dict("branch"=>i, "epsilon"=>eps, "dual"=>raw[i],
                "status"=>string(termination_status(pm.model)),
                "finite_difference"=>has_values(pm.model) ? (objective_value(pm.model)-f0)/eps : nothing))
            fix(fixes[i], 0.0; force=true)
        end
        diagnostic["perturbation_checks"] = checks
    end
    return scores, diagnostic
end
