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

"one angle reference per component of the physical PF matrix, including isolated buses"
function pf_component_references(A::SparseArrays.SparseMatrixCSC)
    visited = falses(size(A, 1))
    references = Int64[]
    pending = Int64[]
    for root in eachindex(visited)
        visited[root] && continue
        push!(references, root)
        visited[root] = true
        push!(pending, root)
        while !isempty(pending)
            column = pop!(pending)
            for entry in nzrange(A, column)
                row = A.rowval[entry]
                # Sparse branch updates retain stored zeros; those are not edges.
                if !iszero(A.nzval[entry]) && !visited[row]
                    visited[row] = true
                    push!(pending, row)
                end
            end
        end
    end
    return references
end

"grounds every component without modifying the physical matrix used for residual checks"
function grounded_pf_matrix(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, references=pf_component_references(A))
    grounded = copy(A)
    for reference in references
        grounded[reference, reference] += max(abs(A[reference, reference]), 1.0)
    end
    return grounded
end

"checks linear solve residuals"
function valid_pf_linear_solution(A, x, rhs)
    all(isfinite, x) || return false
    # allow injection imbalance up to solver precision
    return norm(A * x - rhs, Inf) <= max(SOLVER_PRECISION, 1e-10 * norm(rhs, Inf))
end

"checks whether the matrix has an empty column"
function _is_matrix_connected(A::SparseArrays.SparseMatrixCSC)::Bool
    @inbounds for i in 1:size(A, 2)
        if A.colptr[i] == A.colptr[i+1]
            return false
        end
    end
    return true
end

# https://www.emis.de/journals/ETNA/vol.26.2007/pp161-177.dir/pp161-177.pdf
"applies the rank-1 update A' = A + u*vᵗ to the existing factorization of A L, U"
function recursive_lu_update!(L, U, u, v)
    m, n = size(U)  # assuming U is an m×n matrix

    for i in 1:m
        # Diagonal update
        U[i, i] += u[i] * v[i]
        v[i] = v[i] / U[i, i]

        # L update
        for j in (i+1):m
            u[j] -= u[i] * L[j, i]
            L[j, i] += v[i] * u[j]
        end

        # U update
        for j in (i+1):n
            U[i, j] += u[i] * v[j]
            v[j] -= v[i] * U[i, j]
        end
    end
end

"applies the rank-1 update A' = A + u*vᵗ to the existing factorization of A: `L*U == (row_scaling .* A)[row_perm, col_perm]`. u=(x, -x)ᵗ and v=(1, -1)ᵗ at positions i1,i2, otherwise zeros. This is in principle the same as `recursive_lu_update!(L, U, u, v)`, but optimized and adapted to handle row and col permutations and row scaling."
function recursive_lu_update!(L::SparseArrays.SparseMatrixCSC{Float64, Int64}, U::SparseArrays.SparseMatrixCSC{Float64, Int64}, x::Float64, i1::Int, i2::Int, row_perm::Vector{Int64}, col_perm::Vector{Int64}, row_scaling::Vector{Float64})
    # switching branch (i,j) adds x ∈ {bᵢⱼ, -bᵢⱼ} to Aᵢᵢ, Aⱼⱼ and -x to Aᵢⱼ, Aⱼᵢ
    # we thus need to update L, U by (  x  -x ) = (  x )           = u vᵗ
    #                                ( -x   x )   ( -x ) ( 1  -1 )
    m, n = size(U)  # assuming U is an m×n matrix
    
    ## build u and v
    u = zeros(Float64, m)
    v = zeros(Float64, n)
    
    # we need to find the (inversely) permutated positions of i1 and i2
    # column positions
    i1_c = findfirst(==(i1), col_perm)
    i2_c = findfirst(==(i2), col_perm)
    #row positions
    i1_r = findfirst(==(i1), row_perm)
    i2_r = findfirst(==(i2), row_perm)

    # apply row scaling (before permutation) - u corresponds to rows
    u[i1_r] =  x * row_scaling[i1]
    u[i2_r] = -x * row_scaling[i2]
    
    v[i1_c] =  1
    v[i2_c] = -1

    ## track position in sparse L matrix - this is faster than finding the position again everytime
    # L and U have vectors colptr, rowval, nzval
    #   - colptr[i] is the position of the first element of the first row of the i-th column (A[:,i]) in rowval
    #   - nzval[j] holds the value of the matrix entry and corresponds to rowval[j]
    #   - nzval[colptr[i]] = A[rowval[colptr[i]], i]
    L_id = 1

    # we need to start at the smaller row, then follow the algorithm
    for i in min(i1_r, i2_r):m
        # Diagonal update
        u_i = u[i]
        v_i = v[i]
        
        if v_i == 0.0 && u_i == 0.0
            continue
        end

        # U[i, i] += u[i] * v[i]
        if v_i != 0.0
            U_id = U.colptr[i+1] - 1 # the last row of each column of this upper triangular matrix is a diagonal entry

            if u_i != 0.0
                U_ii = U.nzval[U_id] += u_i * v_i
            else
                U_ii = U.nzval[U_id]
            end

            # v[i] = v[i] / U[i, i]
            v[i] = v_i = v_i / U_ii
        end
        
        # colptr[i] has the first position of row entries of col i, then the first entry is a diagonal entry
        L_id = L.colptr[i] + 1

        # L update
        for j in (i+1):m
            ## find L_ji
            # go until we found the row or we are out of the col range
            while L_id < L.colptr[i+1] && L.rowval[L_id] < j
                L_id += 1
            end
            rowval = L_id < L.colptr[i+1] ? L.rowval[L_id] : 0

            # if we are later in the same col, or in a later col, the entry is not 0.0 -> try to add the entry
            if rowval ≠ j
                
                if v_i * u[j] ≠ 0.0
                    # set entry to new value
                    L[j, i] = v_i * u[j]
                end
                
                continue
            end

            L_ji = L.nzval[L_id]

            # u[j] -= u[i] * L[j, i]
            u_j = u[j] -= u_i * L_ji

            # L[j, i] += v[i] * u[j]
            if v_i * u_j != 0.0
                L.nzval[L_id] = L_ji + v_i * u_j
            end
        end

        # U update
        for j in (i+1):n
            
            # find U_ij - in the range of col j, but before the diagonal entry, which is the last entry
            U_id = U.colptr[j]
            last_U_id = U.colptr[j+1] - 1
            rowval = U.rowval[U_id]
            while rowval < i && U_id < last_U_id
                U_id += 1
                rowval = U.rowval[U_id]
            end
            
            # U[i, j] += u[i] * v[j]
            u_i_x_v_j = u_i * v[j]
            if rowval ≠ i
                
                if u_i_x_v_j ≠ 0.0
                    # set entry to new value
                    U[i, j] = u_i_x_v_j
                else
                    continue
                end
                
            else # U[i, j] += u_i_x_v_j
                U.nzval[U_id] += u_i_x_v_j
            end

            # v[j] -= v[i] * U[i, j]
            val = v_i * U.nzval[U_id]
            if val ≠ 0.0
                v[j] -= val
            end
        end
    end
end

"TODO this is not quite correct, but `recursive_lu_update!` already is pretty good!
Applies the rank-1 update A' = A + u*vᵗ to the existing factorization of A L, U, where u=(x, -x)ᵗ and v=(1, -1)ᵗ at positions i<j, otherwise zeros.
This function explicitly tries to avoid handling updates w/o changes to L or U (due to zeros...)."
function recursive_lu_update_rolled_up!(L::SparseArrays.SparseMatrixCSC{Float64, Int64}, U::SparseArrays.SparseMatrixCSC{Float64, Int64}, x::Float64, i1::Int, i2::Int, row_perm::Vector{Int64}, col_perm::Vector{Int64}, row_scaling::Vector{Float64})
    # for i,j perform diagonal, L, U updates
    # Diagonal update i
    i1_c = findfirst(==(i1), col_perm)
    i2_c = findfirst(==(i2), col_perm)
    i1_r = findfirst(==(i1), row_perm)
    i2_r = findfirst(==(i2), row_perm)

    x1 =  x * row_scaling[i1]
    x2 = -x * row_scaling[i2]

    if i1_c < i2_c
        min_i_c = i1_c
        max_i_c = i2_c

        u1 = x1
        u2 = x2
    else
        min_i_c = i2_c
        max_i_c = i1_c

        u1 = x1
        u2 = x2
    end

    if i1_r < i2_r
        min_i_r = i1_r
        max_i_r = i2_r

        v1 =  1
        v2 = -1
    else
        min_i_r = i2_r
        max_i_r = i1_r

        v1 = -1
        v2 =  1
    end

    U[min_i_c, min_i_c] += u1 * v1                  # U[i,i] = U[i,i] + u[i] * v[i]    = U[i,i] + x * 1
    v1     = v1 / U[min_i_c, min_i_c]           # v[i]   =           v[i] / U[i,i] = 1 / U[i,i]

    # L update i: for j
    u2     -=  u1 * L[max_i_c,min_i_c]  # u[j]   = u[j] - u[i] * L[j,i]    = - x - x * L[j,i]       = x * (-1 - L[j,i])
    L[max_i_c,min_i_c] += v1 *  u2          # L[j,i] = L[j,i] + v[i] * u[j]

    # U update i: for j
    U[min_i_c,max_i_r] += u1 * v2                  # U[i,j] = U[i,j] + u[i] * v[j]    = U[i,j] + x * (-1)      = U[i,j] - x
    v2     -= v1 * U[min_i_c, max_i_c] # v[j]   = v[j] - v[i] * U[i,j]    = -1 - 1/U[i,i] * U[i,j] = -1 - U[i,j]/U[i,i]

    # Diagonal update j
    U[max_i_c,max_i_c] += u2 * v2          # U[j,j] = U[j,j] + u[j] * v[j]

    # i < j+1, so no more L/U updates
    return
end