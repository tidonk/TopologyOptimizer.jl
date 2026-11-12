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

"solves a pre-(L,U,p,q,Rs)-factorized linear equation system Ax = b"
function _solve_factorized(L::SparseArrays.SparseMatrixCSC{Float64, Int64}, U::SparseArrays.SparseMatrixCSC{Float64, Int64}, row_perm::Vector{Int64}, col_perm::Vector{Int64}, row_scaling::Vector{Float64}, b::Vector{Float64})::Vector{Float64}
    z = U \ (L \ (row_scaling .* b)[row_perm])
    x = z[invperm(col_perm)]

    if !all(isfinite, x)
        throw(LinearAlgebra.SingularException(0))
    end

    return x
end

"solves a pre-(L,U,p,q,Rs)-factorized linear equation system Ax = b stored in the solver"
function _solve_factorized_opt(solver::AbstractLinEqSysRfSolver, b::Vector{Float64})::Vector{Float64}
    return _solve_factorized(solver.L, solver.U, solver.p, solver.q, solver.Rs, b)
end
"solves a pre-(L,U,p,q,Rs)-factorized linear equation system Ax = `solver.linsolve.b`"
function _solve_factorized_opt(solver::AbstractLinEqSysRfSolver)::Vector{Float64}
    return _solve_factorized_opt(solver, solver.linsolve.b)
end

"factorizes `solver.linsolve.A` into a (L,U,p,q,Rs)-factorization, then solves the equation system Ax = `solver.linsolve.b`"
function _factorize_and_solve_opt!(solver::AbstractLinEqSysRfSolver)::Vector{Float64}
    # ground the cached factors while keeping solver.A unchanged
    solver.grounding_refs = pf_component_references(solver.A)
    solver.linsolve.A = grounded_pf_matrix(solver.A, solver.grounding_refs)
    sol = LinearSolve.solve!(solver.linsolve)
            
    # couldn't solve
    if sol.retcode == LinearSolve.SciMLBase.ReturnCode.Infeasible
        throw(LinearAlgebra.SingularException(0))
    end

    F = LinearSolve.@get_cacheval(solver.linsolve, :UMFPACKFactorization)

    if F.m != 0
        solver.L  = F.L
        solver.U  = F.U
        solver.p  = F.p
        solver.q  = F.q
        solver.Rs = F.Rs

        solver.is_factorized = true
    end
    
    return sol.u
end

# shared functions for LinEqSysRfSolver and LinEqSysBBRfSolver
"solves problem in its current state, returns success flag"
function solve_primal!(solver::AbstractLinEqSysRfSolver)::Bool
    # A remove/restore sequence can damage the updated factors even if the
    # final component references are unchanged. Retry once with fresh factors.
    used_cached_factors = false
    try
        # Rank-one branch updates preserve the old grounding. Splitting or
        # reconnecting components can require a different set of references.
        if solver.is_factorized && pf_component_references(solver.A) != solver.grounding_refs
            solver.is_factorized = false
        end
        used_cached_factors = solver.is_factorized
        if !solver.is_factorized
            # new factorization
            solver.linsolve.A = solver.A
            solver.linsolve.b = solver.power_surplus
            
            solver.vas = _factorize_and_solve_opt!(solver)
        else
            solver.vas = _solve_factorized_opt(solver, solver.power_surplus)
        end
        if !valid_pf_linear_solution(solver.A, solver.vas, solver.power_surplus)
            callback_stats[:n_leq_residual_failures] += 1
            callback_stats[:n_leq_fails] += 1
            solver.is_factorized = false
            return used_cached_factors ? solve_primal!(solver) : false
        end
    catch e
        if isa(e, DimensionMismatch)
            if solver.is_factorized
                @error("$(e) while factorized, likely caused by a change in matrix size compared to the factorized version. Make sure to re-initialize the factorization after every change in matrix dimension. Sizes: L $(size(solver.L)), U $(size(solver.U)), b $(size(solver.power_surplus)), p $(size(solver.p)), q $(size(solver.q)), Rs $(size(solver.Rs)).")
            else
                @error("$(e) while not factorized. Sizes: A $(size(solver.linsolve.A)), b $(size(solver.linsolve.b)), components $(size(solver.components)).")
            end
        elseif isa(e, LinearAlgebra.SingularException)
            callback_stats[:n_leq_fails] += 1
        else
            rethrow(e)
        end
            
        solver.is_factorized = false
        return used_cached_factors ? solve_primal!(solver) : false
    end
    return true
end

"solves dual problem in its current state !after the primal was solved! to obtain a optimality cut, returns success flag"
function solve_dual_opt!(solver::AbstractLinEqSysRfSolver)::Bool
    try
        if !solver.is_factorized
            # new factorization
            solver.linsolve.b = solver.dual_rhs

            solver.duals = _factorize_and_solve_opt!(solver)
        else
            solver.duals = _solve_factorized_opt(solver, solver.dual_rhs)
        end
        if !valid_pf_linear_solution(solver.A, solver.duals, solver.dual_rhs)
            callback_stats[:n_leq_residual_failures] += 1
            callback_stats[:n_leq_fails] += 1
            solver.is_factorized = false
            return false
        end
    catch e
        if isa(e, DimensionMismatch)
            if solver.is_factorized
                @error("$(e) while factorized, likely caused by a change in matrix size compared to the factorized version. Make sure to re-initialize the factorization after every change in matrix dimension. Sizes: L $(size(solver.L)), U $(size(solver.U)), b $(size(solver.power_surplus)), p $(size(solver.p)), q $(size(solver.q)), Rs $(size(solver.Rs)).")
            else
                @error("$(e) while not factorized. Sizes: A $(size(solver.linsolve.A)), b $(size(solver.linsolve.b)), components $(size(solver.components)).")
            end
        elseif isa(e, LinearAlgebra.SingularException)
            callback_stats[:n_leq_fails] += 1
        else
            rethrow(e)
        end

        solver.is_factorized = false
        return false
    end
    return true
end

"solves dual problem in its current state !after the primal was solved! to obtain a feasibility cut, returns success flag"
function solve_dual_feas!(solver::AbstractLinEqSysRfSolver)::Bool
    try
        if solver.is_factorized

            ## We have A = LU and want L'U'=A', where
            ##
            ##                 [ A   c ]   [ L   0 ]   [ U  u ]   [  LU       Lu ]
            ##            A' = [ rᵀ  d ] = [ lᵀ  1 ] * [ 0  v ] = [ lᵀU  lᵀu + v ]
            ##
            ## First, compute new row lᵀ for L'/ col u for U':
            ##
            ##        Uᵀ l = r     L u = c
            ##
            ## Then the missing diagonal entry v for U':
            ##
            ##             v = d - lᵀu
            ##
            ## We also need to apply the row scaling factors solver.Rs as well as row and column permutations solver.q, solver.q.
            ## For this, it suffices to scale and permute the entries of the new row and column that are within the dimensions of the exisiting factorization accordingly, and then add them at the end positions, without scaling the new row.

            n, m = size(solver.A)

            new_col = (solver.feas_col .* solver.Rs)[solver.p] # c
            
            new_row = solver.feas_row[1:1, solver.q]           # rᵀ
            u = solver.L\new_col                               # L u = c
            l = transpose(solver.U)\transpose(new_row)         # Uᵀ l = r
            lᵀ = transpose(l)
            v = solver.feas_diag - (lᵀ * u)[1,1]               # v = d - lᵀu
            
            new_L = [solver.L spzeros(n); lᵀ 1]                # [ L   0 ] ; [ lᵀ  1 ]
            new_U = [solver.U u; spzeros(1, m) v]              # [ U   u ] ; [ 0   v ]
            
            Rs = append!(copy(solver.Rs), 1.0)                 # extend scaling to last row
            p = append!(copy(solver.p), n+1)
            q = append!(copy(solver.q), m+1)

            solver.duals = _solve_factorized(new_L, new_U, p, q, Rs, solver.dual_rhs)
        else
            if !_is_matrix_connected(solver.A_feas_cut)
                callback_stats[:n_leq_fails] += 1
                solver.dual_rhs=zeros(solver.n_buses)
                return false
            end
            solver.duals = solver.A_feas_cut\solver.dual_rhs
        end
    catch e
        if isa(e, LinearAlgebra.SingularException)
            callback_stats[:n_leq_fails] += 1
            solver.dual_rhs=zeros(solver.n_buses)
            return false
        else
            rethrow(e)
        end
    end
    solver.dual_rhs=zeros(solver.n_buses)
    return true
end