# Adapted from SCIP 10.0.1 src/scip/cutsel_dynamic.c, Christoph Graczyk.
# Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB).
# Licensed under the Apache License, Version 2.0.
# https://github.com/scipopt/scip/blob/v10.0.1/src/scip/cutsel_dynamic.c
# Julia port of efficacy-only selection with randnumgen=NULL and no cut pool.

"SCIP dynamic filtering with finite efficacy and absolute cosine parallelism."
function scip_dynamic_indices(efficacy::Vector{Float64}, parallelism::Matrix{Float64};
                              forced::Vector{Int}=Int[], maxselected::Int=length(efficacy),
                              filtermode::Symbol=:dynamic, mingain::Float64=0.01,
                              minortho::Float64=0.9)
    @assert filtermode in (:dynamic, :full_dynamic)
    @assert 0 <= mingain <= 1 && 0 <= minortho <= 1 && maxselected >= 0
    @assert size(parallelism) == (length(efficacy), length(efficacy))
    @assert all(isfinite, efficacy) && all(isfinite, parallelism)
    epsilon = 1e-9 # SCIP numerics/epsilon, not its feasibility tolerance.
    maxparall = 1 - minortho
    scores = copy(efficacy)
    candidates = [i for i in eachindex(efficacy) if !(i in forced)]

    function drop!(pool, position)
        pool[position] = pool[end]
        pop!(pool)
    end
    function score!(pool)
        for position in length(pool):-1:1
            if efficacy[pool[position]] <= epsilon
                drop!(pool, position)
            end
        end
    end
    function best_first!(pool)
        best = 1
        for i in 2:length(pool)
            if scores[pool[i]] > scores[pool[best]]
                best = i
            end
        end
        pool[1], pool[best] = pool[best], pool[1]
    end
    function filter!(best)
        a = efficacy[best]
        @assert a > 0
        for position in length(candidates):-1:1
            i = candidates[position]
            b = efficacy[i]
            cosine = parallelism[best, i]
            ratio = a - b >= -epsilon ? cosine * a / b : cosine * b / a
            threshold = max((a*a + b*b - (1+mingain)^2*a*a*(1-cosine*cosine))/(2*a*b), maxparall)
            if ratio - 1 >= -epsilon || cosine - threshold > epsilon
                drop!(candidates, position)
            end
        end
    end
    function project!(best)
        a = efficacy[best]
        for i in candidates
            b = efficacy[i]
            cosine = parallelism[best, i]
            scores[i] = abs(cosine-1) <= epsilon ? -1e20 :
                sqrt(a*a + b*b - 2*abs(a)*abs(b)*cosine)/sqrt(1-cosine*cosine) - a
        end
    end

    score!(candidates)
    for i in forced
        isempty(candidates) && break
        filter!(i)
    end
    selected = Int[]
    (isempty(candidates) || maxselected == 0) && return selected
    if filtermode == :full_dynamic && !isempty(forced)
        goodforced = copy(forced)
        score!(goodforced)
        if !isempty(goodforced)
            best_first!(goodforced)
            project!(goodforced[1])
        end
    end
    while !isempty(candidates)
        best_first!(candidates)
        best = popfirst!(candidates)
        push!(selected, best)
        length(selected) == maxselected && break
        filter!(best)
        filtermode == :full_dynamic && project!(best)
    end
    return selected
end

"Adapt SCIP's dynamic selector to first-level Benders cuts and existing budgets."
function select_cuts_scip_dynamic!(subproblem_prio_list, subproblem_ids, cons,
                                  physical_violations, cb_data, cut_types; params)
    isempty(cons) && return cons, subproblem_ids
    @assert params[:LIMIT_MODE] in (:abs_num, :rel_num, :rel_num_per_iter)
    first_level = [group[1] for group in cons]
    @assert all(c -> c.set isa JuMP.MOI.LessThan, first_level)
    variables = unique([v for c in first_level for v in keys(extract_cut_coefficients(c))])
    coefficients = extract_cut_coefficients_vector(first_level, variables)
    norms = sqrt.(vec(sum(abs2, coefficients; dims=2)))
    # Evaluate the actual inequalities at the callback candidate, just as SCIP
    # evaluates rows at its LP point. Outside callbacks the caller supplies values.
    violation = isnothing(cb_data) ? copy(physical_violations) :
        [JuMP.value(v -> JuMP.callback_value(cb_data, v), c.func) - c.set.upper for c in first_level]
    @assert all(isfinite, violation) "SCIP dynamic selection requires finite cut violations."
    # SCIProwGetLPEfficacy floors the norm at numerics/sumepsilon (default 1e-6).
    efficacy = [norms[i] == 0 ? 0.0 : max(0.0, violation[i])/max(norms[i], 1e-6) for i in eachindex(cons)]
    # A violated constant inequality must not disappear through zero-norm scoring.
    constants = findall(i -> norms[i] == 0 && violation[i] > 0, eachindex(cons))
    force_feasibility = get(params, :FORCE_FEASIBILITY, true)
    forced = isempty(cut_types) || !force_feasibility ? Int[] :
        findall(i -> cut_types[i] == :Feasibility && efficacy[i] > 0, eachindex(cons))
    # SCIP uses |cos(theta)|, unlike the signed cosine distance used in clustering.
    parallelism = clamp.(abs.(1 .- compute_cut_distances(coefficients, :cosine)), 0.0, 1.0)
    cap = get_num_select_cuts(params, length(cons), length(subproblem_prio_list))
    selected = scip_dynamic_indices(efficacy, parallelism; forced,
        maxselected=max(0, cap-length(forced)-length(constants)),
        filtermode=get(params, :FILTER_MODE, :dynamic),
        mingain=get(params, :MIN_GAIN, 0.01), minortho=get(params, :MIN_ORTHO, 0.9))
    retained = vcat(constants, forced, selected)
    # Lazy Benders separation must retain a candidate even when SCIP's epsilon
    # would discard all rows. This does not apply to SCIP's optional LP cuts.
    isempty(retained) && push!(retained, argmax(violation))
    return cons[retained], subproblem_ids[retained]
end
