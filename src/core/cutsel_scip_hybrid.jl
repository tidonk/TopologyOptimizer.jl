# Adapted from SCIP 10.0.1 src/scip/{cutsel_hybrid.c,misc.c,lp.c,set.c}.
# Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB), Apache-2.0.
# https://github.com/scipopt/scip/tree/v10.0.1/src/scip
# Ports the default hybrid selector, not the SCIP solver or its row processing.
# Default weights: efficacy=1, objective parallelism=0.1, integral support=0.1,
# directed cutoff distance=0; root/nonroot minortho=0.9. Newly generated rows
# have no global-pool membership, hence no 1e-4 pool bonus. The adapter uses
# original JuMP columns and an integer callback candidate, rather than SCIP's
# presolved LP rows/point. Optional forced cuts and empty-selection recovery
# are Benders caller conventions. The new campaign forces no cuts by type.
# Optional rows pass SCIP's post-selection efficacy acceptance only after
# hybrid filtering; forced rows are exempt. No ordinary refill is performed.
# Model-owned RNG state persists across batches, just as hybrid's plugin state;
# SEED is a SCIP global seed shift, not a replacement for its fixed plugin seed.

"SCIP's KISS generator, with hybrid seed 0x5EED and global seed shift."
mutable struct SCIPHybridRNG
    seed::UInt32
    xor_seed::UInt32
    mwc_seed::UInt32
    cst_seed::UInt32
end

function SCIPHybridRNG(seedshift::Integer=42; multiplier::Integer=10)
    initial = UInt32((0x5eed + (1 + 6 * multiplier) * seedshift) % UInt32)
    hash(a) = UInt32(((UInt64(a) + 0xd37e9a1ce2148403) *
                      (UInt64(initial) + 0xe5fcc163aef32782)) >> 32)
    return SCIPHybridRNG(max(hash(123456789), UInt32(1)),
        max(hash(362436000), UInt32(1)), max(hash(521288629), UInt32(1)), hash(7654321))
end

function scip_hybrid_random!(rng::SCIPHybridRNG)
    rng.seed = UInt32((UInt64(rng.seed) * 1103515245 + 12345) % UInt32)
    rng.xor_seed ⊻= rng.xor_seed << 13
    rng.xor_seed ⊻= rng.xor_seed >> 17
    rng.xor_seed ⊻= rng.xor_seed << 5
    t = UInt64(698769069) * rng.mwc_seed + rng.cst_seed
    rng.cst_seed = UInt32(t >> 32)
    rng.mwc_seed = UInt32(t % UInt32)
    return Float64(rng.seed + rng.xor_seed + rng.mwc_seed) / Float64(typemax(UInt32)) * 1e-6
end

"Draw once per row in upstream scoring order: forced, then optional rows."
function scip_hybrid_noise!(rng::SCIPHybridRNG, forced, n::Integer)
    noise = zeros(n)
    for i in forced
        noise[i] = scip_hybrid_random!(rng)
    end
    for i in 1:n
        i in forced || (noise[i] = scip_hybrid_random!(rng))
    end
    return noise
end

"Default SCIP hybrid score/filter kernel; returns ordered optional row indices."
function scip_hybrid_indices(efficacy::Vector{Float64}, parallelism::Matrix{Float64};
        objective_parallelism=zeros(length(efficacy)), integral_support=zeros(length(efficacy)),
        forced::Vector{Int}=Int[], maxselected::Int=length(efficacy), minortho::Float64=0.9,
        goodscore::Float64=0.9, badscore::Float64=0.0, score_noise=zeros(length(efficacy)))
    n = length(efficacy)
    @assert size(parallelism) == (n, n)
    @assert all(length(v) == n && all(isfinite, v) for v in
        (efficacy, objective_parallelism, integral_support, score_noise))
    @assert all(isfinite, parallelism) && 0 <= minortho <= 1 && maxselected >= 0
    @assert length(unique(forced)) == length(forced) && all(i -> 1 <= i <= n, forced)
    # Preserve C's arithmetic order and strict comparisons. All Benders rows
    # are newly generated, so the global-cut-pool bonus does not apply.
    scores = [0.1*objective_parallelism[i] + 0.1*integral_support[i] +
              efficacy[i] + score_noise[i] for i in 1:n]
    threshold = max(0.0, maximum(scores; init=0.0))
    badthreshold, goodthreshold = threshold*badscore, threshold*goodscore
    maxparall, goodmaxparall = 1.0-minortho, max(0.5, 1.0-minortho)
    pool = [i for i in 1:n if !(i in forced)]
    function filter!(row)
        for p in length(pool):-1:1
            i = pool[p]
            if parallelism[row, i] > (scores[i] >= goodthreshold ? goodmaxparall : maxparall)
                pool[p] = pool[end]
                pop!(pool)
            end
        end
    end
    for i in forced
        filter!(i)
    end
    selected = Int[]
    while !isempty(pool)
        best = 1
        for p in 2:length(pool)
            scores[pool[p]] > scores[pool[best]] && (best = p)
        end
        pool[1], pool[best] = pool[best], pool[1]
        row = pool[1]
        scores[row] < badthreshold && break
        push!(selected, row)
        length(selected) == maxselected && break
        popfirst!(pool)
        filter!(row)
    end
    return selected
end

"SCIP's separation limits, mapped to Gurobi's explored-node count."
function scip_hybrid_budget(params, n_cuts, n_scenarios; node_count=nothing)
    if params[:LIMIT_MODE] == :scip_default
        isnothing(node_count) && error("Native SCIP limits require a callback node count")
        @assert isfinite(node_count) && node_count >= 0
        cap = node_count < 1 ? get(params, :MAXCUTSROOT, 2000) : get(params, :MAXCUTS, 100)
        @assert cap >= 1
        return min(n_cuts, cap)
    end
    return get_num_select_cuts(params, n_cuts, n_scenarios)
end

function scip_hybrid_node_count(cb_data)
    isnothing(cb_data) && error("Native SCIP limits require a Gurobi callback")
    count = Ref{Cdouble}()
    rc = Gurobi.GRBcbget(cb_data, Gurobi.GRB_CB_MIPSOL, Gurobi.GRB_CB_MIPSOL_NODCNT, count)
    rc == 0 || error("Could not read Gurobi callback node count: $rc")
    return count[]
end

"Apply SCIP hybrid ranking/filtering to original affine Benders rows."
function select_cuts_scip_hybrid!(subproblem_prio_list, subproblem_ids, cons,
        physical_violations, cb_data, cut_types; params)
    isempty(cons) && return cons, subproblem_ids
    @assert params[:LIMIT_MODE] in (:abs_num, :rel_num, :rel_num_per_iter, :scip_default)
    rows = [group[1] for group in cons]
    @assert all(c -> c.set isa JuMP.MOI.LessThan, rows)
    variables = unique([v for c in rows for v in keys(extract_cut_coefficients(c))])
    coefficients = extract_cut_coefficients_vector(rows, variables)
    norms = sqrt.(vec(sum(abs2, coefficients; dims=2)))
    violation = isnothing(cb_data) ? copy(physical_violations) :
        [JuMP.value(v -> JuMP.callback_value(cb_data, v), c.func) - c.set.upper for c in rows]
    @assert all(isfinite, violation)
    # A wholly constant batch has no owning model and no geometric ranking.
    # A single violated constant contradiction suffices to reject the candidate.
    if isempty(variables)
        maximum(violation) > 0 || error("SCIP_HYBRID: generated constant cut pool has no positive raw violation; cannot separate the callback candidate.")
        retained = [argmax(violation)]
        return cons[retained], subproblem_ids[retained]
    end
    efficacy = violation ./ max.(norms, 1e-6)
    model = JuMP.owner_model(first(variables))
    acceptance_stats = get!(model.ext, :scip_hybrid_acceptance_stats) do
        Dict{Symbol,Int}(:calls => 0, :optional_rejected => 0,
            :empty_selection_rescues => 0, :all_weak_emergency_batches => 0,
            :all_weak_emergency_rows => 0)
    end
    acceptance_stats[:calls] += 1
    post_efficacy_tol = get(params, :POST_EFFICACY_TOL, 1e-6)
    @assert isfinite(post_efficacy_tol) && post_efficacy_tol >= 0
    objective = JuMP.objective_function(model)
    @assert objective isa JuMP.GenericAffExpr || objective isa JuMP.VariableRef
    objcoef(v) = objective isa JuMP.VariableRef ? Float64(v == objective) : JuMP.coefficient(objective, v)
    # SCIP's objective norm is over the whole LP, including columns absent
    # from this batch of rows. Objective constants do not contribute.
    objnorm2 = sum(abs2(objcoef(v)) for v in JuMP.all_variables(model))
    obj = objcoef.(variables)
    objective_parallelism = [norms[i]^2 * objnorm2 > 1e-9 ?
        clamp(abs(sum(coefficients[i,j]*obj[j] for j in eachindex(variables))) /
              sqrt(norms[i]^2*objnorm2), 0.0, 1.0) : 0.0 for i in eachindex(rows)]
    integral = [JuMP.is_binary(v) || JuMP.is_integer(v) for v in variables]
    integral_support = [count(j -> coefficients[i,j] != 0 && integral[j], eachindex(variables)) /
        max(1, count(!iszero, coefficients[i,:])) for i in eachindex(rows)]
    constants = findall(i -> norms[i] == 0 && violation[i] > 0, eachindex(rows))
    forced = isempty(cut_types) || !get(params, :FORCE_FEASIBILITY, false) ? Int[] :
        unique(vcat(constants, findall(i -> cut_types[i] == :Feasibility && violation[i] > 0, eachindex(rows))))
    # BLAS computes the Gram matrix; normalize afterwards so zero rows have
    # explicit parallelism zero without allocating a full denominator matrix.
    parallelism = coefficients * transpose(coefficients)
    for i in eachindex(rows), j in eachindex(rows)
        prod = norms[i]*norms[j]
        parallelism[i,j] = prod > 0 ? clamp(abs(parallelism[i,j])/prod, 0.0, 1.0) : 0.0
    end
    rng = get!(model.ext, :scip_hybrid_rng) do
        SCIPHybridRNG(get(params, :SEED, 42))
    end
    noise = scip_hybrid_noise!(rng, forced, length(rows))
    node_count = params[:LIMIT_MODE] == :scip_default ? scip_hybrid_node_count(cb_data) : nothing
    if !isnothing(node_count)
        counter = node_count < 1 ? :native_root_calls : :native_nonroot_calls
        acceptance_stats[counter] = get(acceptance_stats, counter, 0) + 1
    end
    budget = max(0, scip_hybrid_budget(params, length(rows), length(subproblem_prio_list); node_count)-length(forced))
    selected = budget == 0 ? Int[] : scip_hybrid_indices(efficacy, parallelism;
        objective_parallelism, integral_support, forced, maxselected=budget,
        minortho=get(params, :MIN_ORTHO, 0.9), score_noise=noise)
    # Match SCIP's optional-row application gate after hybrid has completed
    # scoring and similarity filtering. Retain its order without ordinary refill.
    accepted = [i for i in selected if efficacy[i] > post_efficacy_tol]
    acceptance_stats[:optional_rejected] += length(selected) - length(accepted)
    retained = vcat(forced, accepted)
    if isempty(retained)
        if any(>(post_efficacy_tol), efficacy)
            # Exceptional Benders recovery: hybrid may have filtered stronger
            # rows using a selected row that the application gate then rejects.
            push!(retained, argmax(efficacy))
            acceptance_stats[:empty_selection_rescues] += 1
        else
            # Preserve the budget and reject the candidate with one violated row.
            maximum(violation) > 0 || error("SCIP_HYBRID: generated cut pool has no positive raw violation; cannot separate the callback candidate.")
            retained = [argmax(efficacy)]
            acceptance_stats[:all_weak_emergency_batches] += 1
            acceptance_stats[:all_weak_emergency_rows] += length(retained)
        end
    end
    return cons[retained], subproblem_ids[retained]
end
