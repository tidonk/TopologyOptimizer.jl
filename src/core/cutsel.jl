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
using Clustering
using Distances
using StatsBase: sample, quantile
using Statistics: mean

"extract variable-to-coefficient mapping from a Benders cut constraint."
function extract_cut_coefficients(con)::Dict{JuMP.VariableRef, Float64}
    coeffs = Dict{JuMP.VariableRef, Float64}()
    
    # Extract coefficients from the affine expression
    expr = con.func
    if isa(expr, JuMP.GenericAffExpr)
        for (var, coeff) in expr.terms
            coeffs[var] = coeff
        end
    end
    
    return coeffs
end

"extract coefficients from multiple cuts into a matrix (rows=cuts, cols=variables)."
function extract_cut_coefficients_vector(cons::Vector, all_vars::Vector{JuMP.VariableRef})::Matrix{Float64}
    n_cuts = length(cons)
    n_vars = length(all_vars)
    
    # Create a mapping from variable to index
    var_to_idx = Dict(var => i for (i, var) in enumerate(all_vars))
    
    # Initialize the coefficient matrix
    coeff_matrix = zeros(Float64, n_cuts, n_vars)
    
    for (i, con) in enumerate(cons)
        coeffs = extract_cut_coefficients(con)
        for (var, coeff) in coeffs
            if haskey(var_to_idx, var)
                coeff_matrix[i, var_to_idx[var]] = coeff
            end
        end
    end
    
    return coeff_matrix
end

"map distance metric symbol from settings to a Distances.jl metric."
function _resolve_distance_metric(metric::Symbol)::SemiMetric
    metric == :jaccard   && return Jaccard()
    metric == :cosine    && return CosineDist()
    metric == :hamming   && return Hamming()
    metric == :euclidean && return Euclidean()
    error("Unknown distance metric: $metric. Use :jaccard, :cosine, :hamming, or :euclidean")
end

"compute pairwise distance matrix for cut coefficient rows using Distances.jl."
function compute_cut_distances(coeff_matrix::Matrix{Float64}, metric::Symbol)::Matrix{Float64}
    dist_metric = _resolve_distance_metric(metric)
    # Distances.pairwise expects observations as columns
    data = permutedims(coeff_matrix)
    # For Jaccard/Hamming, booleanize to non-zero support
    if dist_metric isa Jaccard || dist_metric isa Hamming
        data = Float64.(abs.(data) .> SOLVER_PRECISION)
    end
    return Distances.pairwise(dist_metric, data, dims=2)
end

"select a representative from each cluster (mode: :max_violation or :centroid)."
function select_cluster_representatives(clusters::Vector{Int}, violations::Vector{Float64};
                                        mode::Symbol=:max_violation,
                                        coeff_matrix::Union{Matrix{Float64}, Nothing}=nothing)::Vector{Int}
    @assert length(clusters) == length(violations) "Clusters and violations must have same length"
    if mode == :centroid
        @assert !isnothing(coeff_matrix) "coeff_matrix is required for :centroid mode"
        @assert size(coeff_matrix, 1) == length(clusters) "coeff_matrix rows must match cluster length"
    end
    
    representatives = Int[]
    unique_clusters = unique(clusters)
    
    for cluster_id in unique_clusters
        cluster_indices = findall(c -> c == cluster_id, clusters)
        @assert(! isempty(cluster_indices))
        
        if mode == :centroid && length(cluster_indices) > 1
            # Compute centroid and pick closest member
            centroid = vec(mean(coeff_matrix[cluster_indices, :], dims=1))
            best_idx = 0
            best_dist = Inf
            for ci in cluster_indices
                d = sum((coeff_matrix[ci, :] .- centroid).^2)
                if d < best_dist
                    best_dist = d
                    best_idx = ci
                end
            end
            push!(representatives, best_idx)
        else
            # Default: pick the point with maximum violation
            cluster_violations = violations[cluster_indices]
            max_idx = argmax(cluster_violations)
            push!(representatives, cluster_indices[max_idx])
        end
    end
    
    return sort(representatives)
end

"compute pairwise distance matrix using equal-frequency binning + Jaccard distance."
function compute_cut_distances_eqfreq(coeff_matrix::Matrix{Float64}; n_bins::Int=10)::Matrix{Float64}
    n_cuts, n_vars = size(coeff_matrix)
    # Build a binary indicator matrix: one-hot encode bin assignments per variable
    binned = zeros(Float64, n_cuts, n_vars * n_bins)
    for j in 1:n_vars
        col = coeff_matrix[:, j]
        nonzero_mask = abs.(col) .> SOLVER_PRECISION
        if !any(nonzero_mask)
            continue
        end
        # compute quantile edges from nonzero values
        nonzero_vals = col[nonzero_mask]
        edges = quantile(nonzero_vals, range(0.0, 1.0, length=n_bins+1))
        for i in 1:n_cuts
            if !nonzero_mask[i]
                continue  # zero coefficients get no bin assignment
            end
            # find bin: last edge ≤ value (clamp to [1, n_bins])
            bin = searchsortedlast(edges, col[i])
            bin = clamp(bin, 1, n_bins)
            binned[i, (j-1)*n_bins + bin] = 1.0
        end
    end
    return Distances.pairwise(Jaccard(), permutedims(binned), dims=2)
end


"select cuts from subproblems based on violation magnitude."
function select_cuts_by_violation!(subproblem_prio_list::Vector{Dict{Symbol, Any}}, 
                                   subproblem_ids::Vector{Int64},
                                   cons::Vector{Vector{JuMP.ScalarConstraint}},
                                   violations::Vector{Float64},
                                   cb_data::Union{Nothing, Any};
                                   params::Dict{Symbol, Any}
                                   )::Tuple{Vector{Vector{JuMP.ScalarConstraint}}, Vector{Int64}}

    @assert(length(cons) == length(violations), "Length of cuts and violations does not match.")
    @assert(length(cons) == length(subproblem_ids), "Length of cuts and subproblem_ids does not match.")

    # sort indices by descending violation
    perm = sortperm(violations, rev=true)

    # determine how many cuts to select
    selected = Int[]
    if params[:LIMIT_MODE] == :abs_num
        n = min(length(cons), params[:LIMIT_ABS_NUM])
        selected = perm[1:n]

    elseif params[:LIMIT_MODE] == :rel_num
        # constant fraction of total subproblems (instance-dependent, iteration-independent)
        alpha = params[:LIMIT_REL_ALPHA]
        n_total_sps = length(subproblem_prio_list)
        n = max(1, ceil(Int, alpha * n_total_sps))
        n = min(n, length(cons))
        selected = perm[1:n]

    elseif params[:LIMIT_MODE] == :rel_num_per_iter
        # fraction of violated cuts in this iteration
        alpha = params[:LIMIT_REL_ITER_ALPHA]
        n = max(1, ceil(Int, alpha * length(cons)))
        n = min(n, length(cons))
        selected = perm[1:n]

    elseif params[:LIMIT_MODE] == :pb
        # add cuts until cumulative violation exceeds ρ · (z_UB - z_MP)
        if ! isnothing(cb_data)
            z_ub = Ref{Cdouble}()
            z_mp = Ref{Cdouble}()
            Gurobi.GRBcbget(cb_data, Gurobi.GRB_CB_MIPSOL, Gurobi.GRB_CB_MIPSOL_OBJBST, z_ub)
            Gurobi.GRBcbget(cb_data, Gurobi.GRB_CB_MIPSOL, Gurobi.GRB_CB_MIPSOL_OBJ, z_mp)
        else
            @error("Limit mode :pb requires a callback data object to get primal bound and current objective.")
        end
        rho = params[:LIMIT_PB_RHO]
        gap = max(0.0, z_ub[] - z_mp[])
        cumulative_violation_cutoff = rho * gap

        cumulative_violation = 0.0
        for i in perm
            cumulative_violation += violations[i]
            push!(selected, i)
            if cumulative_violation > cumulative_violation_cutoff
                break
            end
        end
    else
        @error("Unknown cut selection limit mode: $(params[:LIMIT_MODE])")
    end

    cuts = cons[selected]
    c_ids_used = subproblem_ids[selected]

    return cuts, c_ids_used
end

"""
    select_cuts_by_violation_agg!(...)

Select the top-k cuts by violation as representatives, then aggregate all remaining
cuts into one combined cut.  Returns k+1 cuts total (or fewer when |cons| ≤ k).
"""
function select_cuts_by_violation_agg!(subproblem_prio_list::Vector{Dict{Symbol, Any}},
                                       subproblem_ids::Vector{Int64},
                                       cons::Vector{Vector{JuMP.ScalarConstraint}},
                                       violations::Vector{Float64},
                                       cb_data::Union{Nothing, Any};
                                       params::Dict{Symbol, Any}
                                      )::Tuple{Vector{Vector{JuMP.ScalarConstraint}}, Vector{Int64}}

    # First, select representatives using the standard violation logic
    rep_cuts, rep_ids = select_cuts_by_violation!(subproblem_prio_list, subproblem_ids, cons, violations, cb_data, params=params)

    # Aggregate all remaining (non-selected) cuts into one combined cut
    rep_id_set = Set(rep_ids)
    all_remaining = [i for i in 1:length(cons) if !(subproblem_ids[i] in rep_id_set)]
    _aggregate_remaining_cuts!(rep_cuts, rep_ids, cons, subproblem_ids, violations, all_remaining, params=params)

    return rep_cuts, rep_ids
end

"""
    _cluster_cuts(cons, violations, params; n_total_subproblems) -> (clusters, dist_matrix, coeff_matrix)

Shared clustering logic for diversity-based cut selection methods.
Computes pairwise distances on level-1 cut coefficients and returns cluster assignments.
`n_total_subproblems` is used for `:rel` cluster count mode.
"""
function _cluster_cuts(cons::Vector{Vector{JuMP.ScalarConstraint}},
                       violations::Vector{Float64};
                       params::Dict{Symbol, Any},
                       n_total_subproblems::Int=0
                      )::Tuple{Vector{Int}, Matrix{Float64}, Matrix{Float64}}

    # Use only the level-1 (switching-master) cut from each subproblem for clustering.
    first_level_cons = [group[1] for group in cons]

    # Collect all unique variables from the level-1 cuts
    all_vars = Set{JuMP.VariableRef}()
    for con in first_level_cons
        coeffs = extract_cut_coefficients(con)
        union!(all_vars, keys(coeffs))
    end
    all_vars_vec = collect(all_vars)

    # Extract coefficient matrix and pairwise distances
    coeff_matrix = extract_cut_coefficients_vector(first_level_cons, all_vars_vec)
    distance_metric = get(params, :DISTANCE_METRIC, :jaccard)
    if distance_metric == :eqfreq_bin
        n_bins = get(params, :EQFREQ_NUM_BINS, 10)
        dist_matrix = compute_cut_distances_eqfreq(coeff_matrix, n_bins=n_bins)
    else
        dist_matrix = compute_cut_distances(coeff_matrix, distance_metric)
    end

    if get(params, :USE_DBSCAN, true)
        eps = get(params, :DBSCAN_EPS, 0.5)
        min_pts = get(params, :DBSCAN_MIN_PTS, 2)
        result = Clustering.dbscan(dist_matrix, eps; min_neighbors=min_pts, metric=nothing)
        clusters = result.assignments
        next_id = maximum(clusters) + 1
        for i in eachindex(clusters)
            if clusters[i] == 0
                clusters[i] = next_id
                next_id += 1
            end
        end
    else
        n = length(cons)
        cluster_mode = get(params, :NUM_CLUSTERS_MODE, :abs)
        if cluster_mode == :rel
            alpha = get(params, :NUM_CLUSTERS_REL_ALPHA, 0.05)
            k = max(1, ceil(Int, alpha * n_total_subproblems))
            k = min(k, n)
        else  # :abs
            k = min(get(params, :NUM_CLUSTERS, 5), n)
        end
        if k >= n
            clusters = collect(1:n)
        else
            try
                clusters = Clustering.kmedoids(dist_matrix, k).assignments
            catch e
                @warn "kmedoids failed ($(typeof(e))), falling back to top-k by violation" exception=e
                perm = sortperm(violations; rev=true)
                clusters = zeros(Int, n)
                for (cid, idx) in enumerate(perm[1:k])
                    clusters[idx] = cid
                end
                selected = perm[1:k]
                for idx in perm[k+1:end]
                    _, nearest = findmin(dist_matrix[idx, j] for j in selected)
                    clusters[idx] = clusters[selected[nearest]]
                end
            end
        end
    end

    return clusters, dist_matrix, coeff_matrix
end

"""
    _compute_agg_weights(violations; normalize) -> Vector{Float64}

Compute per-cut aggregation weights from violations.  Weights equal the violation value,
with Inf violations (feasibility cuts boosted by FEASCUT_INF_VIOL) capped at the maximum
finite violation in the group.  If `normalize` is true, the weights are scaled to sum to 1.

Callers should check the `WEIGHTED_AGG` flag *before* calling this function and pass
`weights=nothing` for equal weighting.
"""
function _compute_agg_weights(violations::AbstractVector{Float64}; normalize::Bool=false)::Vector{Float64}
    max_finite = maximum((v for v in violations if isfinite(v)), init=1.0)
    raw = Float64[isfinite(v) ? v : max_finite for v in violations]
    if normalize
        s = sum(raw)
        return s > 0.0 ? raw ./ s : fill(1.0 / length(raw), length(raw))
    end
    return raw
end

"""
    aggregate_cluster_cuts(cluster_cons; weights=nothing) -> Vector{JuMP.ScalarConstraint}

Aggregate a group of cut vectors into a single cut vector by summing LHS expressions
and RHS values level-by-level.  Returns one `ScalarConstraint` per nesting level.

If `weights` is provided (a `Vector{Float64}` of the same length as `cluster_cons`),
each cut is multiplied by its weight before summing.  When `weights` is `nothing`
(the default), all cuts contribute equally (weight = 1).
"""
function aggregate_cluster_cuts(cluster_cons::Vector{Vector{JuMP.ScalarConstraint}};
                                weights::Union{Nothing, Vector{Float64}}=nothing
                               )::Vector{JuMP.ScalarConstraint}
    n_levels = maximum(length(group) for group in cluster_cons)
    agg_cuts = Vector{JuMP.ScalarConstraint}(undef, n_levels)
    for level in 1:n_levels
        agg_func = JuMP.AffExpr(0.0)
        agg_rhs  = 0.0
        for (j, group) in enumerate(cluster_cons)
            if length(group) >= level
                w = isnothing(weights) ? 1.0 : weights[j]
                JuMP.add_to_expression!(agg_func, w, group[level].func)
                agg_rhs += w * group[level].set.upper
            end
        end
        agg_cuts[level] = JuMP.ScalarConstraint(agg_func, JuMP.MOI.LessThan(agg_rhs))
    end
    return agg_cuts
end

"""
    _aggregate_remaining_cuts!(cuts, c_ids_used, cons, subproblem_ids, violations, remaining; params)

Aggregate cuts at `remaining` indices into a single combined cut, push it onto `cuts` and
`c_ids_used`.  Does nothing when `remaining` is empty.
"""
function _aggregate_remaining_cuts!(cuts::Vector{Vector{JuMP.ScalarConstraint}},
                                    c_ids_used::Vector{Int64},
                                    cons::Vector{Vector{JuMP.ScalarConstraint}},
                                    subproblem_ids::Vector{Int64},
                                    violations::Vector{Float64},
                                    remaining::Vector{Int};
                                    params::Dict{Symbol, Any})
    isempty(remaining) && return
    weighted = get(params, :WEIGHTED_AGG, false)
    normalize = get(params, :WEIGHTED_AGG_NORMALIZE, false)
    w = weighted ? _compute_agg_weights(violations[remaining], normalize=normalize) : nothing
    agg_cut = aggregate_cluster_cuts(cons[remaining], weights=w)
    push!(cuts, agg_cut)
    best = remaining[argmax(violations[remaining])]
    push!(c_ids_used, subproblem_ids[best])
end

"select cuts based on diversity using clustering (DBSCAN or k-medoids)."
function select_cuts_by_diversity!(subproblem_prio_list::Vector{Dict{Symbol, Any}},
                                   subproblem_ids::Vector{Int64},
                                   cons::Vector{Vector{JuMP.ScalarConstraint}},
                                   violations::Vector{Float64};
                                   params::Dict{Symbol, Any},
                                   plot_html::Union{String, Nothing}=nothing
                                   )::Tuple{Vector{Vector{JuMP.ScalarConstraint}}, Vector{Int64}}
    
    cuts = Vector{Vector{JuMP.ScalarConstraint}}()
    c_ids_used = Vector{Int64}()
    
    if length(cons) <= 1
        return cons, subproblem_ids
    end

    clusters, dist_matrix, coeff_matrix = _cluster_cuts(cons, violations, params=params, n_total_subproblems=length(subproblem_prio_list))

    # Select representative from each cluster.
    rep_mode = get(params, :REPRESENTATIVE_MODE, :max_violation)
    representative_indices = select_cluster_representatives(clusters, violations, mode=rep_mode, coeff_matrix=coeff_matrix)
    chosen = BitVector(i in representative_indices for i in 1:length(cons))

    # Write HTML visualization if requested
    if !isnothing(plot_html)
        write_cut_cluster_html(plot_html, dist_matrix, clusters, violations, coeff_matrix, chosen)
    end
    
    # Include the full cons[i] group for each representative subproblem
    for idx in representative_indices
        push!(cuts, cons[idx])
        push!(c_ids_used, subproblem_ids[idx])
    end
    
    return cuts, c_ids_used
end

"""
    select_cuts_by_diversity_agg!(...)

Cluster cuts, then aggregate all cuts within each cluster into a single combined cut.
Returns k aggregated cuts (one per cluster).
"""
function select_cuts_by_diversity_agg!(subproblem_prio_list::Vector{Dict{Symbol, Any}},
                                       subproblem_ids::Vector{Int64},
                                       cons::Vector{Vector{JuMP.ScalarConstraint}},
                                       violations::Vector{Float64};
                                       params::Dict{Symbol, Any},
                                       plot_html::Union{String, Nothing}=nothing
                                      )::Tuple{Vector{Vector{JuMP.ScalarConstraint}}, Vector{Int64}}

    if length(cons) <= 1
        return cons, subproblem_ids
    end

    clusters, dist_matrix, coeff_matrix = _cluster_cuts(cons, violations, params=params, n_total_subproblems=length(subproblem_prio_list))
    unique_clusters = unique(clusters)

    cuts = Vector{Vector{JuMP.ScalarConstraint}}()
    c_ids_used = Vector{Int64}()

    weighted = get(params, :WEIGHTED_AGG, false)
    normalize = get(params, :WEIGHTED_AGG_NORMALIZE, false)
    for cid in unique_clusters
        members = findall(c -> c == cid, clusters)
        w = weighted ? _compute_agg_weights(violations[members], normalize=normalize) : nothing
        agg_cut = aggregate_cluster_cuts(cons[members], weights=w)
        push!(cuts, agg_cut)
        best = members[argmax(violations[members])]
        push!(c_ids_used, subproblem_ids[best])
    end

    # Write HTML visualization if requested
    if !isnothing(plot_html)
        chosen = falses(length(cons))  # no individual representative for pure aggregation
        write_cut_cluster_html(plot_html, dist_matrix, clusters, violations, coeff_matrix, chosen)
    end

    return cuts, c_ids_used
end

"""
    select_cuts_by_diversity_agg_hyb!(...)

Cluster cuts, pick the highest-violation representative per cluster, then aggregate
the remaining (non-representative) cuts.

Two modes controlled by `params[:AGG_HYB_MODE]`:
- `:k_plus_1`  – aggregate ALL remaining cuts into 1 combined cut → k+1 cuts total
- `:two_k`     – aggregate remaining cuts per cluster → up to 2k cuts total
"""
function select_cuts_by_diversity_agg_hyb!(subproblem_prio_list::Vector{Dict{Symbol, Any}},
                                           subproblem_ids::Vector{Int64},
                                           cons::Vector{Vector{JuMP.ScalarConstraint}},
                                           violations::Vector{Float64};
                                           params::Dict{Symbol, Any},
                                           plot_html::Union{String, Nothing}=nothing
                                          )::Tuple{Vector{Vector{JuMP.ScalarConstraint}}, Vector{Int64}}

    if length(cons) <= 1
        return cons, subproblem_ids
    end

    clusters, dist_matrix, coeff_matrix = _cluster_cuts(cons, violations, params=params, n_total_subproblems=length(subproblem_prio_list))
    unique_clusters = unique(clusters)

    # Step 1: pick representatives (same as DIVERSITY)
    rep_mode = get(params, :REPRESENTATIVE_MODE, :max_violation)
    representative_indices = select_cluster_representatives(clusters, violations, mode=rep_mode, coeff_matrix=coeff_matrix)
    rep_set = Set(representative_indices)

    cuts = Vector{Vector{JuMP.ScalarConstraint}}()
    c_ids_used = Vector{Int64}()

    for idx in representative_indices
        push!(cuts, cons[idx])
        push!(c_ids_used, subproblem_ids[idx])
    end

    # Step 2: aggregate remaining cuts
    mode = get(params, :AGG_HYB_MODE, :k_plus_1)

    if mode == :two_k
        # per-cluster aggregation of non-representative members
        for cid in unique_clusters
            members = findall(c -> c == cid, clusters)
            remainder = filter(i -> !(i in rep_set), members)
            _aggregate_remaining_cuts!(cuts, c_ids_used, cons, subproblem_ids, violations, remainder, params=params)
        end
    else  # :k_plus_1
        # aggregate ALL remaining cuts into one
        all_remaining = [i for i in 1:length(cons) if !(i in rep_set)]
        _aggregate_remaining_cuts!(cuts, c_ids_used, cons, subproblem_ids, violations, all_remaining, params=params)
    end

    if !isnothing(plot_html)
        chosen = BitVector(i in representative_indices for i in 1:length(cons))
        write_cut_cluster_html(plot_html, dist_matrix, clusters, violations, coeff_matrix, chosen)
    end

    return cuts, c_ids_used
end

function select_cuts_randomly!(subproblem_prio_list::Vector{Dict{Symbol, Any}},
                                   subproblem_ids::Vector{Int64},
                                   cons::Vector{Vector{JuMP.ScalarConstraint}},
                                   violations::Vector{Float64};
                                   params::Dict{Symbol, Any}
                                   )::Tuple{Vector{Vector{JuMP.ScalarConstraint}}, Vector{Int64}}
    
    cuts = Vector{Vector{JuMP.ScalarConstraint}}()
    c_ids_used = Vector{Int64}()
    cuts_mode = get(params, :NUM_CUTS_MODE, :abs)
    if cuts_mode == :rel
        alpha = get(params, :NUM_CUTS_REL_ALPHA, 0.05)
        n_cuts = max(1, ceil(Int, alpha * length(subproblem_prio_list)))
        n_cuts = min(n_cuts, length(cons))
    else  # :abs
        n_cuts = min(length(cons), get(params, :NUM_CUTS, 5))
    end
    selected_indices = sample(1:length(cons), n_cuts; replace=false)

    for idx in selected_indices
        push!(cuts, cons[idx])
        push!(c_ids_used, subproblem_ids[idx])
    end
    
    return cuts, c_ids_used
end

"""
    select_cuts_randomly_agg!(...)

Select k random cuts as representatives, then aggregate all remaining cuts into one
combined cut.  Returns k+1 cuts total (or fewer when |cons| ≤ k).
"""
function select_cuts_randomly_agg!(subproblem_prio_list::Vector{Dict{Symbol, Any}},
                                    subproblem_ids::Vector{Int64},
                                    cons::Vector{Vector{JuMP.ScalarConstraint}},
                                    violations::Vector{Float64};
                                    params::Dict{Symbol, Any}
                                   )::Tuple{Vector{Vector{JuMP.ScalarConstraint}}, Vector{Int64}}

    # First, select representatives using the standard random logic
    rep_cuts, rep_ids = select_cuts_randomly!(subproblem_prio_list, subproblem_ids, cons, violations, params=params)

    # Aggregate all remaining (non-selected) cuts into one combined cut
    rep_id_set = Set(rep_ids)
    all_remaining = [i for i in 1:length(cons) if !(subproblem_ids[i] in rep_id_set)]
    _aggregate_remaining_cuts!(rep_cuts, rep_ids, cons, subproblem_ids, violations, all_remaining, params=params)

    return rep_cuts, rep_ids
end

"""
Hybrid cut selection: filter by mean violation, then cluster the survivors.
Keeps all cuts with violation ≥ arithmetic mean of violations (minimum `MIN_AFTER_FILTER` cuts).
Then runs diversity-based clustering on the filtered set.
"""
function select_cuts_hybrid!(subproblem_prio_list::Vector{Dict{Symbol, Any}},
                             subproblem_ids::Vector{Int64},
                             cons::Vector{Vector{JuMP.ScalarConstraint}},
                             violations::Vector{Float64};
                             params::Dict{Symbol, Any},
                             plot_html::Union{String, Nothing}=nothing
                            )::Tuple{Vector{Vector{JuMP.ScalarConstraint}}, Vector{Int64}}
    
    if length(cons) <= 1
        return cons, subproblem_ids
    end

    min_after_filter = get(params, :MIN_AFTER_FILTER, 5)

    # Step 1: filter by mean violation
    mean_viol = sum(violations) / length(violations)
    above_mean = findall(v -> v >= mean_viol, violations)
    if length(above_mean) < min_after_filter
        # keep the top min_after_filter by descending violation
        perm = sortperm(violations, rev=true)
        above_mean = perm[1:min(min_after_filter, length(cons))]
    end

    # If after filtering, we have fewer cuts than requested clusters, return all filtered
    n_clusters = get(params, :NUM_CLUSTERS, 5)
    if length(above_mean) <= n_clusters
        cuts = cons[above_mean]
        c_ids_used = subproblem_ids[above_mean]
        return cuts, c_ids_used
    end

    # Step 2: cluster the filtered cuts using diversity-based selection
    filtered_cons = cons[above_mean]
    filtered_ids = subproblem_ids[above_mean]
    filtered_violations = violations[above_mean]

    # Build diversity params for the sub-call (pass through clustering config)
    diversity_params = Dict{Symbol, Any}(
        :USE_DBSCAN => get(params, :USE_DBSCAN, false),
        :NUM_CLUSTERS => n_clusters,
        :NUM_CLUSTERS_MODE => get(params, :NUM_CLUSTERS_MODE, :abs),
        :NUM_CLUSTERS_REL_ALPHA => get(params, :NUM_CLUSTERS_REL_ALPHA, 0.05),
        :DISTANCE_METRIC => get(params, :DISTANCE_METRIC, :jaccard),
        :EQFREQ_NUM_BINS => get(params, :EQFREQ_NUM_BINS, 10),
        :DBSCAN_EPS => get(params, :DBSCAN_EPS, 0.5),
        :DBSCAN_MIN_PTS => get(params, :DBSCAN_MIN_PTS, 5),
        :REPRESENTATIVE_MODE => get(params, :REPRESENTATIVE_MODE, :max_violation),
    )

    return select_cuts_by_diversity!(subproblem_prio_list, filtered_ids, filtered_cons, 
                                     filtered_violations, params=diversity_params,
                                     plot_html=plot_html)
end

"dispatch cut selection: feasibility-first guard, then method-specific selection."
function select_cuts!(subproblem_prio_list::Vector{Dict{Symbol, Any}}, 
                     subproblem_ids::Vector{Int64},
                     cons::Vector{Vector{JuMP.ScalarConstraint}},
                     violations::Vector{Float64},
                     cb_data::Union{Gurobi.CallbackData, Nothing},
                     cut_types::Vector{Symbol}=Symbol[];
                     n_callbacks::Int=0,
                    )::Tuple{Vector{Vector{JuMP.ScalarConstraint}}, Int, Int, Float64}

    # if FEASCUT_INF_VIOL is enabled, boost feasibility cut violations to Inf
    # so they sort to the top during normal selection (also adding opt cuts).
    if PARAMS[:BENDERS][:CUTS][:SELECTION][:FEASCUT_INF_VIOL] && !isempty(cut_types)
        feas_indices = findall(t -> t == :Feasibility, cut_types)
        if !isempty(feas_indices)
            violations = copy(violations)
            for idx in feas_indices
                violations[idx] = Inf
            end
            @debug("FEASCUT_INF_VIOL: boosted $(length(feas_indices)) feasibility cuts to Inf violation, proceeding with normal selection on all $(length(cons)) cuts.")
        end
    end

    # if FEASCUT_ONLY is enabled, add all feasibility cuts and skip optimality cut selection
    if PARAMS[:BENDERS][:CUTS][:SELECTION][:FEASCUT_ONLY] && !isempty(cut_types)
        feas_indices = findall(t -> t == :Feasibility, cut_types)
        if !isempty(feas_indices)
            cuts = cons[feas_indices]
            c_ids_used = subproblem_ids[feas_indices]
            selected_viols = violations[feas_indices]
            n_feas = length(feas_indices)
            n_opt = 0
            @debug("FEASCUT_ONLY: selected $(n_feas) feasibility cuts, skipping $(length(cons) - n_feas) optimality cuts.")
            # update which cuts will actually be used
            for c_id in c_ids_used
                subproblem_prio_list[c_id][:n_cuts_used] += 1
                violation = violations[findfirst(x->x == c_id, subproblem_ids)]
                if violation < Inf
                    subproblem_prio_list[c_id][:sum_violations] += violation
                end
            end
            return (cuts, n_feas, n_opt, sum(v for v in selected_viols if v < Inf; init=0.0))
        end
    end

    # if no cons, no selection, or in a cut stabilization round, use all cuts
    if isempty(cons) || ! PARAMS[:BENDERS][:CUTS][:SELECTION][:ENABLED] || (PARAMS[:BENDERS][:CUTS][:CUT_STAB_ROUND_FREQ] > 0 && mod(n_callbacks - 1, PARAMS[:BENDERS][:CUTS][:CUT_STAB_ROUND_FREQ]) != 0)
        cuts = cons
        c_ids_used = subproblem_ids
    else
        method = PARAMS[:BENDERS][:CUTS][:SELECTION][:METHOD]
        method_params = PARAMS[:BENDERS][:CUTS][:SELECTION][method]
        # produce HTML plot on the first callback only (for DIVERSITY and HYBRID)
        plot_html = nothing
        if n_callbacks == 1 && !isempty(PARAMS[:BENDERS][:CUTS][:SELECTION][:PLOT_HTML])
            plot_html = PARAMS[:BENDERS][:CUTS][:SELECTION][:PLOT_HTML]
        end

        if method == :VIOLATION
            cuts, c_ids_used = select_cuts_by_violation!(subproblem_prio_list, subproblem_ids, cons, violations, cb_data, params=method_params)
        elseif method == :VIOLATION_AGG
            cuts, c_ids_used = select_cuts_by_violation_agg!(subproblem_prio_list, subproblem_ids, cons, violations, cb_data, params=method_params)
        elseif method == :RANDOM
            cuts, c_ids_used = select_cuts_randomly!(subproblem_prio_list, subproblem_ids, cons, violations, params=method_params)
        elseif method == :RANDOM_AGG
            cuts, c_ids_used = select_cuts_randomly_agg!(subproblem_prio_list, subproblem_ids, cons, violations, params=method_params)
        elseif method == :DIVERSITY
            cuts, c_ids_used = select_cuts_by_diversity!(subproblem_prio_list, subproblem_ids, cons, violations, 
                                                          params=method_params, plot_html=plot_html)
        elseif method == :DIVERSITY_AGG
            cuts, c_ids_used = select_cuts_by_diversity_agg!(subproblem_prio_list, subproblem_ids, cons, violations,
                                                              params=method_params, plot_html=plot_html)
        elseif method == :DIVERSITY_AGG_HYB
            cuts, c_ids_used = select_cuts_by_diversity_agg_hyb!(subproblem_prio_list, subproblem_ids, cons, violations,
                                                                   params=method_params, plot_html=plot_html)
        elseif method == :HYBRID
            cuts, c_ids_used = select_cuts_hybrid!(subproblem_prio_list, subproblem_ids, cons, violations, 
                                                    params=method_params, plot_html=plot_html)
        else
            @error("Unknown cut selection method: $method")
        end
        @assert(length(cuts) > 0)

        # ensure at least MIN_OPT_CUTS optimality cuts are in the selection (if available)
        min_opt_cuts = PARAMS[:BENDERS][:CUTS][:SELECTION][:MIN_OPT_CUTS]
        if min_opt_cuts > 0 && !isempty(cut_types)
            # count opt cuts currently selected
            n_opt_in_sel = count(c_id -> begin
                idx = findfirst(x -> x == c_id, subproblem_ids)
                cut_types[idx] != :Feasibility
            end, c_ids_used)
            # count available (unselected) opt cuts
            if n_opt_in_sel < min_opt_cuts
                selected_set = Set(c_ids_used)
                unselected_opt = [(i, subproblem_ids[i], violations[i]) for i in eachindex(subproblem_ids) 
                                  if cut_types[i] != :Feasibility && !(subproblem_ids[i] in selected_set)]
                if !isempty(unselected_opt)
                    # sort unselected opt cuts by violation descending (best first)
                    sort!(unselected_opt, by=x -> x[3], rev=true)
                    n_to_add = min(min_opt_cuts - n_opt_in_sel, length(unselected_opt))
                    for k in 1:n_to_add
                        (orig_idx, opt_c_id, _) = unselected_opt[k]
                        push!(cuts, cons[orig_idx])
                        push!(c_ids_used, opt_c_id)
                    end
                    @debug("MIN_OPT_CUTS: added $(n_to_add) optimality cut(s) on top of $(length(cuts) - n_to_add) selected cuts.")
                end
            end
        end
    end

    # count feas/opt among selected cuts and compute sum of selected violations
    sum_selected_violations = 0.0
    n_feas_selected = 0
    n_opt_selected = 0
    for c_id in c_ids_used
        subproblem_prio_list[c_id][:n_cuts_used] += 1
        idx = findfirst(x->x == c_id, subproblem_ids)
        violation = violations[idx]
        if violation < Inf
            subproblem_prio_list[c_id][:sum_violations] += violation
            sum_selected_violations += violation
        end
        if !isempty(cut_types) && cut_types[idx] == :Feasibility
            n_feas_selected += 1
        else
            n_opt_selected += 1
        end
    end
    return (cuts, n_feas_selected, n_opt_selected, sum_selected_violations)
end