"""
    pb_selected_indices(violations, is_feasibility, priority, z_ub, z_mp; rho=1.10)

Select PB cuts using unboosted physical violations and an explicit priority
permutation. Without finite bounds, retain the entire pool. With finite bounds,
the first violated feasibility cut rejects the candidate; its certificate units
are never compared with an objective gap. Otherwise retain an optimality-cut
prefix until its physical violation sum strictly exceeds `rho * max(0, UB-MP)`.
Always retain at least one cut from a nonempty optimality pool.
"""
function pb_selected_indices(violations::AbstractVector{<:Real},
                             is_feasibility::AbstractVector{Bool},
                             priority::AbstractVector{<:Integer},
                             z_ub::Real, z_mp::Real; rho::Real=1.10)::Vector{Int}
    n = length(violations)
    length(is_feasibility) == n || throw(ArgumentError("cut type mask length differs from violations"))
    length(priority) == n && sort(collect(priority)) == collect(1:n) ||
        throw(ArgumentError("priority must be a permutation of cut indices"))
    isfinite(rho) && rho >= 1 || throw(ArgumentError("rho must be finite and at least one"))
    all(v -> isfinite(v) && v >= 0, violations) ||
        throw(ArgumentError("physical violations must be finite and nonnegative"))
    if isnan(z_ub) || isnan(z_mp)
        throw(ArgumentError("bounds must not be NaN"))
    end
    isempty(priority) && return Int[]
    # Nonfinite bounds provide no usable objective-gap stopping threshold.
    (!isfinite(z_ub) || !isfinite(z_mp)) && return Int[priority...]
    for i in priority
        if is_feasibility[i] && violations[i] > 0
            return [Int(i)]
        end
    end
    cutoff = rho * max(zero(z_ub - z_mp), z_ub - z_mp)
    selected = Int[]
    cumulative_violation = 0.0
    for i in priority
        push!(selected, i)
        cumulative_violation += violations[i]
        cumulative_violation > cutoff && break
    end
    return selected
end
