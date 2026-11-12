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

"extract a canonical contingency key as sorted outaged branch IDs."
function get_contingency_branch_key(contingency::Contingency)::Vector{Int64}
    if haskey(contingency, :gen) && !isempty(contingency[:gen])
        @warn("Generator contingency identity mapping not implemented yet. Using branch key only.")
    end
    return sort(get(contingency, :branch, Int64[]))
end

"""A serialized Benders cut, stored independently of JuMP variable references.

Fields:
- `constant`: the constant term of the cut expression
- `coefficients`: maps `(var_basename, var_idx)` to the coefficient value
- `type`: `:optimality` or `:feasibility`
- `contingency_id`: the contingency index `c` at serialization time
- `contingency_key`: sorted outaged branch IDs — canonical identity for mapping across stages
"""
struct SerializedCut
    constant::Float64
    coefficients::Dict{Tuple{Symbol, Int64}, Float64}
    type::Symbol
    contingency_id::Int64
    contingency_key::Vector{Int64}
end

"serialize a JuMP cut constraint into a `SerializedCut` independent of variable references."
function serialize_cut(cut_expr::JuMP.AffExpr, type::Symbol, contingency_id::Int64, contingency_key::Vector{Int64}=Int64[])::SerializedCut
    coefficients = Dict{Tuple{Symbol, Int64}, Float64}()
    for (var_ref, coeff) in cut_expr.terms
        if coeff != 0.0
            basename, idx = get_var_basename_idx(var_ref)
            coefficients[(basename, idx)] = coeff
        end
    end
    return SerializedCut(cut_expr.constant, coefficients, type, contingency_id, contingency_key)
end

"serialize a vector of constraints with contingency identity and append to pool."
function serialize_constraints!(pool::Vector{SerializedCut}, constraints::Vector{JuMP.ScalarConstraint}, contingencies::Vector{Contingency}, subproblem_ids::Vector{Int64})
    for (i, con) in enumerate(constraints)
        # con.func is `cut_expr - η` (optimality) or `cut_expr` (feasibility)
        func = con.func::JuMP.AffExpr

        # Separate η term from the rest
        coefficients = Dict{Tuple{Symbol, Int64}, Float64}()
        type = :feasibility
        contingency_id = subproblem_ids[i]
        contingency_key = get_contingency_branch_key(contingencies[contingency_id])

        for (var_ref, coeff) in func.terms
            basename, idx = get_var_basename_idx(var_ref)
            if basename == :eta
                # This is the η variable: optimality cut (coefficient is always -1 from JuMP normalization of `expr <= η`)
                @assert isapprox(coeff, -1.0, atol=SOLVER_PRECISION) "Unexpected η coefficient: $coeff"
                type = :optimality
                # Don't include η in coefficients — it's the RHS
            elseif !isapprox(coeff, 0.0, atol=SOLVER_PRECISION)
                coefficients[(basename, idx)] = coeff
            end
        end

        push!(pool, SerializedCut(func.constant, coefficients, type, contingency_id, contingency_key))
    end
end

"filter and retain active backed cuts from the pool based on solution slack tolerance."
function filter_active_cuts!(data::Dict, result::Dict, activity_tol::Float64)::Vector{SerializedCut}
    pool = get(data, "_benders_cuts_pool", SerializedCut[])
    active_cuts = SerializedCut[]

    if !isempty(pool)
        # Extract solution values from the result
        sol = result["solution"]["nw"]["1"]
        sol_vals = Dict{Tuple{Symbol, Int64}, Float64}()

        # Map from solution component_name to a list of (json_key, var_sym) pairs we should harvest.
        # var_sym must agree with what `get_var_basename_idx` returns for the corresponding JuMP variable.
        sol_field_map = Dict{String, Vector{Tuple{String, Symbol}}}(
            "branch" => [("br_status", :z_branch), ("br_change", :z_branch_change), ("tr_shift", :tr_shift), ("tr_shift_change", :tr_shift_change)],
            "gen"    => [("pg", :pg), ("gen_status", :z_gen)],
            "eta"    => [("η", :η)],
        )

        for (comp_name, field_specs) in sol_field_map
            haskey(sol, comp_name) || continue
            for (comp_id_str, comp_sol) in sol[comp_name]
                comp_id = parse(Int64, comp_id_str)
                for (json_key, var_sym) in field_specs
                    if haskey(comp_sol, json_key)
                        sol_vals[(var_sym, comp_id)] = comp_sol[json_key]
                    end
                end
            end
        end

        n_missing = 0
        for cut in pool
            # Evaluate cut_expr at solution: constant + Σ coeff * val. Skip cuts that
            # reference master variables not present in the solution dict (conservative).
            lhs_val = cut.constant
            cut_ok = true
            for ((basename, idx), coeff) in cut.coefficients
                val = get(sol_vals, (basename, idx), nothing)
                if val === nothing
                    cut_ok = false
                    n_missing += 1
                    break
                end
                lhs_val += coeff * val
            end
            cut_ok || continue

            # RHS: η[c] for optimality, 0 for feasibility
            rhs_val = cut.type == :optimality ? get(sol_vals, (:η, cut.contingency_id), 0.0) : 0.0

            # Cut is: lhs ≤ rhs, so slack = rhs - lhs
            if rhs_val - lhs_val <= activity_tol
                push!(active_cuts, cut)
            end
        end
        if n_missing > 0
            @warn("filter_active_cuts!: $(n_missing) cut(s) referenced master variables not present in the solution dict and were skipped.")
        end
    end

    data["_active_benders_cuts_stage"] = active_cuts
    delete!(data, "_benders_cuts_pool")
    return active_cuts
end

"accumulate or replace carried Benders cuts based on the selected mode."
function accumulate_carried_cuts!(data::Dict, mode::Symbol)
    stage_cuts = pop!(data, "_active_benders_cuts_stage", SerializedCut[])

    if mode == :accumulate
        carried = get!(data, "_carried_benders_cuts", SerializedCut[])
        append!(carried, stage_cuts)
    else  # :last
        data["_carried_benders_cuts"] = copy(stage_cuts)
    end
end

"inject carried-forward Benders cuts into the master problem, folding missing vars."
function inject_carried_cuts!(pm::AbstractPowerModel, data::Dict, contingencies::Vector{Contingency}; nw::Int64=nw_id_default)
    carried = get(data, "_carried_benders_cuts", SerializedCut[])
    if isempty(carried)
        return
    end

    nw_data = data["nw"][string(nw)]

    # Build lookup from (basename, idx) → VariableRef for current master variables
    nw_vars = var(pm, nw)
    var_lookup = Dict{Tuple{Symbol, Int64}, JuMP.VariableRef}()
    for var_name in [:z_branch, :pg, :z_gen, :tr_shift, :z_branch_change, :tr_shift_change]
        if haskey(nw_vars, var_name)
            vars = nw_vars[var_name]
            for var_idx in axes(vars, 1)
                var_lookup[var_name, var_idx] = vars[var_idx]
            end
        end
    end
    # η variables
    if haskey(nw_vars, :η)
        ηs = nw_vars[:η]
        for c in axes(ηs, 1)
            var_lookup[:η, c] = ηs[c]
        end
    end

    # Build contingency identity mapping: sorted outaged branch IDs → current positional index
    # Note: generator contingency mapping is not implemented yet
    key_to_idx = Dict{Vector{Int64}, Int64}()
    for (c, cont) in enumerate(contingencies)
        key = get_contingency_branch_key(cont)
        key_to_idx[key] = c
    end

    n_injected = 0
    n_skipped_contingency = 0
    n_folded_vars = 0
    for cut in carried
        # Resolve contingency identity to current index; skip if contingency no longer exists
        new_c = get(key_to_idx, cut.contingency_key, 0)
        if new_c == 0
            n_skipped_contingency += 1
            continue
        end

        expr = JuMP.AffExpr(cut.constant)

        for ((basename, idx), coeff) in cut.coefficients
            if haskey(var_lookup, (basename, idx))
                JuMP.add_to_expression!(expr, coeff, var_lookup[basename, idx])
            else
                # Variable doesn't exist in this master — fold with its known value
                if basename == :z_branch
                    branch = get(nw_data["branch"], string(idx), nothing)
                    status = branch !== nothing ? Float64(get(branch, "status", 1)) : 1.0
                    expr.constant += coeff * status
                    n_folded_vars += 1
                elseif basename == :z_branch_change
                    # z_branch_change for non-switchable branches: fold as 0 (no change possible)
                else
                    @warn("Carried cut references variable $(basename)[$(idx)] which does not exist in the current master. Folding with value 0.")
                    # pg, z_gen, etc. — if truly missing, fold as 0
                end
            end
        end

        # Add the cut as a constraint
        if cut.type == :optimality
            η = var_lookup[:η, new_c]
            JuMP.@constraint(pm.model, expr <= η)
            n_injected += 1
        else
            JuMP.@constraint(pm.model, expr <= 0)
            n_injected += 1
        end
    end

    @info("Injected $(n_injected) carried-forward Benders cuts ($(n_folded_vars) z_branch folded, $(n_skipped_contingency) skipped due to missing contingency).")
end
