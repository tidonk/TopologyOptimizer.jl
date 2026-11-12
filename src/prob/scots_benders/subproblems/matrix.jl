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

# https://stackoverflow.com/a/38029416
#function resize_sparse!(sp::SparseMatrixCSC,m,n)
#    newcolptr = sp.colptr
#    resize!(newcolptr,n+1)
#    for i=sp.n+2:n+1
#      newcolptr[i] = sp.colptr[sp.n+1] 
#    end
#    return SparseMatrixCSC(m,n,newcolptr,sp.rowval,sp.nzval)
#end
"insert empty rows & columns into `sp` such that it reaches size `m`=sp.m+1, with the empty row and column in position pos"
function resize_sparse(sp::SparseMatrixCSC,m, pos)::SparseMatrixCSC
    @assert m == sp.m+1
    sp_w_cols = [sp[1:pos-1, :]; spzeros(1,m-1); sp[pos:m-1, :]]
    
    return([sp_w_cols[:, 1:pos-1] spzeros(m,1) sp_w_cols[:, pos:m-1]])
end

"delete the row and column at position `pos` of the matrix"
function delete_matrix_row_and_col(sp::SparseMatrixCSC, pos::Int64)::SparseMatrixCSC
    return sp[[begin:pos-1;pos+1:end],[begin:pos-1;pos+1:end]]
end
  

"[unused] Returns matrix `A`, vector `power_surplus` such that A * va = power_surplus for a master solution given by `master_var_vals`,
    where vaᵤ is the voltage angle at bus u"
function build_pf_matrix(pm::AbstractPowerModel,
                         master_var_vals::Dict{Symbol, AbstractVector},
                         off_branches::Vector{Int64};
                         nw::Int64=nw_id_default,
                         ref_bus::Bool = false
                        )
    buses    = ref(pm, nw, :bus)
    branches = ref(pm, nw, :branch)
    connected_branchess = get_connected_branchess(pm, nw=nw)
    n_buses = length(buses)
    n_branches = length(branches)
    
    bus_to_col = zeros(Int64, n_buses)
    col_counter = 1
    
    # the reference bus is left out from the matrix
    ref_bus_id = 0
    if ref_bus
        @assert length(ids(pm, nw, :ref_buses)) == 1
        for id in ids(pm, nw, :ref_buses)
            ref_bus_id = id
            bus_to_col[ref_bus_id] = n_buses
            #col_to_bus[n_buses] = ref_bus_id
        end
    end
    
    # build sparse matrix
    #n_entries = 2 * (n_branches - length(off_branches) - length(connected_branchess[ref_bus_id])) + n_buses - 1
    n_entries = 2 * (n_branches - length(off_branches)) + n_buses
    if ref_bus
        n_entries -= (2*length(connected_branchess[ref_bus_id]) + 1)
    end
    col_ids = fill(0,   n_entries)
    row_ids = fill(0,   n_entries)
    values  = fill(0.0, n_entries)
    entry_counter = 1
    function _push_to_matrix(row_bus_id::Int64, col_bus_id::Int64, value::Real, pos::Int64)
        @assert(value ≠ 0)
        row_ids[pos] = bus_to_col[row_bus_id]
        col_ids[pos] = bus_to_col[col_bus_id]
        values[pos]  = value
    end

    # vector "i" for generation/demand surplus - this could also be sparse, but only if significantly fewer than half of all buses do not have any surplus
    power_surplus = fill(0.0, n_buses - (ref_bus ? 1 : 0))

    # for all buses, fill all entries in their row according to the connected branches
    for (bus_id, bus) in buses
        if bus_id == ref_bus_id; continue; end

        b_sum = 0.0
        if bus_to_col[bus_id] == 0
            bus_to_col[bus_id] = col_counter
            #col_to_bus[col_counter] = bus_id
            col_counter += 1
        end
        
        # for all branches, add their bs to the corresponding entries in the row
        connected_branches = connected_branchess[bus_id]
        for branch_id in connected_branches
            if branch_id ∈ off_branches; continue; end
            branch = branches[branch_id]
            other_bus_id = branch["f_bus"] == bus_id ? branch["t_bus"] : branch["f_bus"]
            if other_bus_id == ref_bus_id; continue; end
            
            b = branch["b"]
            
            # part of G_k,k
            b_sum -= b
            
            if bus_to_col[other_bus_id] == 0
                bus_to_col[other_bus_id] = col_counter
                #col_to_bus[col_counter] = other_bus_id
                col_counter += 1
            end
            
            @assert(b ≠ 0, "bus $bus_id has a branch with b = 0")
            # add entry G_k,l
            _push_to_matrix(bus_id, other_bus_id, b, entry_counter)
            entry_counter += 1
            
        end
        
        @assert(b_sum ≠ 0, "bus $bus_id has no connected branches, or cancelling out b values")
        # add entry G_k,k
        _push_to_matrix(bus_id, bus_id, b_sum, entry_counter)
        entry_counter += 1
    end
    
    # actually construct the matrix
    A = sparse(row_ids, col_ids, values)

    # we set the reference bus va to 0, so we can eliminate the corresponding row by substracting it from another row (eg the first row here)
    if ref_bus
        for branch_id in connected_branchess[ref_bus_id]
            if branch_id ∈ off_branches; continue; end
            branch = branches[branch_id]
            other_bus_id = branch["f_bus"] == ref_bus_id ? branch["t_bus"] : branch["f_bus"]
            A[1,other_bus_id] += branch["b"]
        end
    end

    # add load pds to correponding entries in power_surplus
    for (_, load) in ref(pm, nw, :load)
        load_bus_id = load["load_bus"]
        if load_bus_id != ref_bus_id
            power_surplus[bus_to_col[load_bus_id]] -= load["pd"]
        else
            power_surplus[1] -= load["pd"]
        end
    end

    # add gen pgs from master to correponding entries in power_surplus
    for gen_id in axes(master_var_vals[:pg],1)
        gen = ref(pm, nw, :gen, gen_id)
        gen_bus_id = gen["gen_bus"]
        if gen_bus_id != ref_bus_id
            power_surplus[bus_to_col[gen_bus_id]] += master_var_vals[:pg][gen_id]
        else
            power_surplus[1] += master_var_vals[:pg][gen_id]
        end
    end

    return A, power_surplus, bus_to_col
end

"wrapper for build_base_pf_matrix with busbar components"
function build_base_pf_matrix(pm::AbstractPowerModel, components::Vector{BBComponent}, detached_components::Vector{BBComponent}; nw::Int64=nw_id_default)
    bus_to_col::Vector{Int64}=fill(typemin(Int64), length(ref(pm, nw, :bus)))
    @assert(length(ref(pm, nw, :bus)) == maximum([bus_id for (bus_id, _) in ref(pm, nw, :bus)]), "bus ids are not consecutive")

    # for each component, map the contained buses to the components column to be
    col_counter = 0
    for component in components
        col_counter += 1

        for bus_id in component.busbar_buses
            bus_to_col[bus_id] = col_counter
        end
        for bus_id in component.breaker_buses
            bus_to_col[bus_id] = col_counter
        end
    end
    col_counter = 0
    for detached_component in detached_components
        col_counter -= 1
        for bus_id in detached_component.busbar_buses
            bus_to_col[bus_id] = col_counter
        end
    end

    return build_base_pf_matrix(pm, nw=nw, n_cols=col_counter, bus_to_col=bus_to_col, components=components)
end

"returns matrix `A`, vector `power_surplus` such that A * va = power_surplus, where vaᵤ is the voltage angle at bus u with all branches switched on/in their default position and all generators producing nothing,
 as well as vector `bus_to_col` such that the bus u corresponds to column `bus_to_col[u]` in `A`"
function build_base_pf_matrix(pm::AbstractPowerModel; nw::Int64=nw_id_default, n_cols = length(ref(pm, nw, :bus)), bus_to_col::Vector{Int64}=zeros(Int64, length(ref(pm, nw, :bus))), components::Vector{BBComponent}=Vector{BBComponent}())
    buses    = ref(pm, nw, :bus)
    branches = ref(pm, nw, :branch)
    connected_branchess = get_connected_branchess(pm, only_active=false, nw=nw)

    # get number of branches (only "true" branches, no breakers and couplers)
    n_branches = Int64(round(sum([length(c_branches) for c_branches in Base.values(connected_branchess)])/2))

    col_counter = maximum(bus_to_col) + 1
    
    # build sparse matrix
    n_entries = 2 * n_branches + n_cols
    col_ids = fill(0,   n_entries)
    row_ids = fill(0,   n_entries)
    values  = fill(0.0, n_entries)
    entry_counter = 1
    function _push_to_matrix(row_bus_id::Int64, col_bus_id::Int64, value::Real, pos::Int64)
        if pos > length(row_ids)
            @debug("increase pos to $(length(row_ids) + 1)")
            push!(row_ids, bus_to_col[row_bus_id])
            push!(col_ids, bus_to_col[col_bus_id])
            push!(values, value)
        else
            row_ids[pos] = bus_to_col[row_bus_id]
            col_ids[pos] = bus_to_col[col_bus_id]
            values[pos]  = value
        end
    end

    "handles one row of the column, ie adds all b entries into the corresponding cell in the matrix, and then returns their negative sum (to be added as a diagonal entry)"
    function _comp_bus_sus(bus_id::Int64, col_counter::Int64, entry_counter::Int64; visited_branch_ids::Set{Int64}=Set{Int64}())::Tuple{Float64, Int64, Int64}
        b_sum = 0.0

        # < 0 this is a detached breaker bus, or part of a detached component
        if bus_to_col[bus_id] < 0
            return (0.0, col_counter, entry_counter)
        end

        # 0 => this bus is uninitialized, ie it should have its own column but is not assigned one yet
        if bus_to_col[bus_id] == 0

            bus_to_col[bus_id] = col_counter
            col_counter += 1
        end
        
        # for all branches, add their bs to the corresponding entries in the row
        for branch_id in connected_branchess[bus_id]
            branch_id in visited_branch_ids && continue

            branch = branches[branch_id]
            other_bus_id = branch["f_bus"] == bus_id ? branch["t_bus"] : branch["f_bus"]
            
            if bus_to_col[other_bus_id] < 0
                continue
            end
            #@assert(branch["type"] ∉ ["K", "breaker"], "branch type is $branch[type]")
            #@assert(haskey(branch, "b"), "branch has no b: $branch")
            b = branch["b"]
            
            # part of G_k,k
            b_sum -= b
            
            if bus_to_col[other_bus_id] == 0

                bus_to_col[other_bus_id] = col_counter
                col_counter += 1
            end
            
            # add entry G_k,l
            _push_to_matrix(bus_id, other_bus_id, b, entry_counter)
            entry_counter += 1
            
            push!(visited_branch_ids, branch_id)
        end

        return (b_sum, col_counter, entry_counter)
    end
    
    # for all buses, fill all entries in their row according to the connected branches
    if !isempty(components)
        # fill the entries for each component from each of its buses. Since we cannot easily access the contents of the matrix while we are constructing it here, we need to aggregate all values that share a cell, here this should only be G_k,k
        for component in components
            component_buses = union(component.breaker_buses, component.busbar_buses) # usually, this is either a singleton component -> all lines are attached to the single bus, or a normal component -> all lines are attached to the breaker buses. However in this case, if a line is not switchable, it is directly attached to the busbar.

            if isempty(component_buses) # no buses -> ignore
                continue
            end
            # we now aggregate all columns (buses) of a single component to one singular column
            # note: in our formulation (due to plugging Kirchhoff's law into the nodal balance eq.), this means simply summing up the negative B matrix values
            b_sum = 0.0

            visited_branch_ids = Set{Int64}()
            for bus_id in component_buses
                (this_b_sum, col_counter, entry_counter) = _comp_bus_sus(bus_id, col_counter, entry_counter, visited_branch_ids=visited_branch_ids)
                
                b_sum += this_b_sum
            end

            # add entry G_k,k
            _push_to_matrix(component_buses[1], component_buses[1], b_sum, entry_counter)
            entry_counter += 1
        end
    end

    # now handle buses where no aggregation has to take place, simply construct matrix as before
    for (bus_id, _) in buses
        if bus_to_col[bus_id] == 0 || bus_to_col[bus_id] > length(components) # either not visited yet (probably not directly connected to a component), or the index is "after" the part of bus_to_col corresponding to components
            (b_sum, col_counter, entry_counter) = _comp_bus_sus(bus_id, col_counter, entry_counter)
            
            # add entry G_k,k
            _push_to_matrix(bus_id, bus_id, b_sum, entry_counter)
            entry_counter += 1
        end
    end
    # actually construct the matrix
    A = sparse(row_ids[begin:entry_counter-1], col_ids[begin:entry_counter-1], values[begin:entry_counter-1])

    # vector "i" for generation/demand surplus - this could also be sparse, but only if significantly fewer than half of all buses do not have any surplus
    power_surplus = zeros(col_counter-1) # col counter already points to the "next" element

    # add load pds to correponding entries in power_surplus
    _apply_loads!(pm, power_surplus, bus_to_col, nw=nw)

    # subtract shunt conductance power losses (P_shunt = gs * V^2, V=1.0 in DC)
    _apply_shunts!(pm, power_surplus, bus_to_col, nw=nw)

    apply_fixed_trafo_vals!(power_surplus, pm, bus_to_col, nw=nw)

    return A, power_surplus, bus_to_col
end

"Returns f and t bus ids incident to the branch"
function get_branch_bus_ids(pm::AbstractPowerModel, branch_id::Int64; nw::Int64=nw_id_default)::Tuple{Int64, Int64}
    branch = ref(pm, nw, :branch, branch_id)::Dict{String, Any}

    f_bus_id = branch["f_bus"]::Int64
    t_bus_id = branch["t_bus"]::Int64

    return f_bus_id, t_bus_id
end

"Returns f and t matrix cols incident to the branch"
function get_branch_cols(pm::AbstractPowerModel, branch_id::Int64, bus_to_col::Vector{Int64}; nw::Int64=nw_id_default)::Tuple{Int64, Int64}
    f_bus_id, t_bus_id = get_branch_bus_ids(pm, branch_id, nw=nw)

    return bus_to_col[f_bus_id], bus_to_col[t_bus_id]
end

"adds/removes a branch to/from A"
function _add_or_remove_branch!(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, pm::AbstractPowerModel, branch_id::Int64, bus_to_col::Vector{Int64}, apply_func::Function; nw::Int64=nw_id_default)::Nothing
    
    branch = ref(pm, nw, :branch, branch_id)::Dict{String, Any}
    b = branch["b"]::Float64

    f_bus = branch["f_bus"]::Int64
    t_bus = branch["t_bus"]::Int64

    @assert(f_bus != t_bus, "cannot add/remove branch $branch_id between buses $(branch["f_bus"]) -> $(branch["t_bus"]) in the same BB component")
    
    for (bus_id, other_bus_id) in [(f_bus, t_bus), (t_bus, f_bus)]
        A[bus_to_col[bus_id], bus_to_col[other_bus_id]] = apply_func(A[bus_to_col[bus_id], bus_to_col[other_bus_id]], b)
        A[bus_to_col[bus_id], bus_to_col[bus_id]]       = apply_func(A[bus_to_col[bus_id], bus_to_col[bus_id]],      -b)
    end
end

"adds/removes a branch to/from L,U"
function _add_or_remove_branch!(L::SparseArrays.SparseMatrixCSC{Float64, Int64}, U::SparseArrays.SparseMatrixCSC{Float64, Int64}, pm::AbstractPowerModel, branch_id::Int64, bus_to_col::Vector{Int64}, apply_func::Function, row_perm::Vector{Int64}, col_perm::Vector{Int64}, row_scaling::Vector{Float64}; nw::Int64=nw_id_default)::Nothing

    start_time = time()

    branch = ref(pm, nw, :branch, branch_id)::Dict{String, Any}
    b = branch["b"]::Float64

    f_bus = bus_to_col[branch["f_bus"]::Int64]
    t_bus = bus_to_col[branch["t_bus"]::Int64]

    @assert(f_bus != t_bus, "cannot add/remove branch $branch_id between buses $(branch["f_bus"]) -> $(branch["t_bus"]) in the same BB component")

    recursive_lu_update!(L, U, apply_func(0, -b), f_bus, t_bus, row_perm, col_perm, row_scaling)

    callback_times[:refactorization] += time() - start_time
    return nothing 
end

"adds a branch to A"
function add_branch!(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, pm::AbstractPowerModel, branch_id::Int64, bus_to_col::Vector{Int64}; nw::Int64=nw_id_default)
    _add_or_remove_branch!(A, pm, branch_id, bus_to_col, +, nw=nw)
end

"removes a branch from A"
function remove_branch!(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, pm::AbstractPowerModel, branch_id::Int64, bus_to_col::Vector{Int64}; nw::Int64=nw_id_default)
    _add_or_remove_branch!(A, pm, branch_id, bus_to_col, -, nw=nw)
end

"adds a branch to L,U"
function add_branch!(L::SparseArrays.SparseMatrixCSC{Float64, Int64}, U::SparseArrays.SparseMatrixCSC{Float64, Int64}, pm::AbstractPowerModel, branch_id::Int64, bus_to_col::Vector{Int64}, row_perm::Vector{Int64}, col_perm::Vector{Int64}, row_scaling::Vector{Float64}; nw::Int64=nw_id_default)
    _add_or_remove_branch!(L, U, pm, branch_id, bus_to_col, +, row_perm, col_perm, row_scaling, nw=nw)
end

"removes a branch from L,U"
function remove_branch!(L::SparseArrays.SparseMatrixCSC{Float64, Int64}, U::SparseArrays.SparseMatrixCSC{Float64, Int64}, pm::AbstractPowerModel, branch_id::Int64, bus_to_col::Vector{Int64}, row_perm::Vector{Int64}, col_perm::Vector{Int64}, row_scaling::Vector{Float64}; nw::Int64=nw_id_default)
    _add_or_remove_branch!(L, U, pm, branch_id, bus_to_col, -, row_perm, col_perm, row_scaling, nw=nw)
end

"adds an original branch to A if both ends are connected"
function add_orig_branch!(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, pm::AbstractPowerModel, branch_id::Int64, bus_to_col::Vector{Int64}; nw::Int64=nw_id_default)
    branch = ref(pm, nw, :branch, branch_id)
    orig_id = branch["orig_index"]

    orig_f_col, orig_t_col = get_branch_cols(pm, orig_id, bus_to_col, nw=nw)

    if orig_f_col > 0 && orig_t_col > 0
        add_branch!(A, pm, orig_id, bus_to_col, nw=nw)
    end
end

"removes an original branch from A if both ends are connected"
function remove_orig_branch!(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, pm::AbstractPowerModel, branch_id::Int64, bus_to_col::Vector{Int64}; nw::Int64=nw_id_default)
    branch = ref(pm, nw, :branch, branch_id)
    orig_id = branch["orig_index"]

    orig_f_col, orig_t_col = get_branch_cols(pm, orig_id, bus_to_col, nw=nw)

    if orig_f_col > 0 && orig_t_col > 0
        remove_branch!(A, pm, orig_id, bus_to_col, nw=nw)
    end
end

"adds an original branch to L,U if both ends are connected"
function add_orig_branch!(L::SparseArrays.SparseMatrixCSC{Float64, Int64}, U::SparseArrays.SparseMatrixCSC{Float64, Int64}, pm::AbstractPowerModel, branch_id::Int64, bus_to_col::Vector{Int64}, row_perm::Vector{Int64}, col_perm::Vector{Int64}, row_scaling::Vector{Float64}; nw::Int64=nw_id_default)
    branch = ref(pm, nw, :branch, branch_id)
    orig_id = branch["orig_index"]

    orig_f_col, orig_t_col = get_branch_cols(pm, orig_id, bus_to_col, nw=nw)

    if orig_f_col > 0 && orig_t_col > 0
        add_branch!(L, U, pm, orig_id, bus_to_col, row_perm, col_perm, row_scaling, nw=nw)
    end
end

"removes an original branch from L,U if both ends are connected"
function remove_orig_branch!(L::SparseArrays.SparseMatrixCSC{Float64, Int64}, U::SparseArrays.SparseMatrixCSC{Float64, Int64}, pm::AbstractPowerModel, branch_id::Int64, bus_to_col::Vector{Int64}, row_perm::Vector{Int64}, col_perm::Vector{Int64}, row_scaling::Vector{Float64}; nw::Int64=nw_id_default)
    branch = ref(pm, nw, :branch, branch_id)
    orig_id = branch["orig_index"]

    orig_f_col, orig_t_col = get_branch_cols(pm, orig_id, bus_to_col, nw=nw)

    if orig_f_col > 0 && orig_t_col > 0
        remove_branch!(L, U, pm, orig_id, bus_to_col, row_perm, col_perm, row_scaling, nw=nw)
    end
end

"applies or removes values of the current master solution to `power_surplus`"
function _apply_master_vals_gen!(pm::AbstractPowerModel, power_surplus::Vector{Float64}, bus_to_col::Vector{Int64}, master_var_vals, apply_func_pg::Function; nw::Int64=nw_id_default)
    for gen_id in axes(master_var_vals[:pg],1)
        gen = ref(pm, nw, :gen, gen_id)::Dict{String, Any}
        gen_bus_id = gen["gen_bus"]::Int64
        gen_col = bus_to_col[gen_bus_id]
        if gen_col > 0
            power_surplus[gen_col] = apply_func_pg(power_surplus[gen_col], master_var_vals[:pg][gen_id])
        end
    end
end

"substracts the (static) load demands from `power_surplus`"
function _apply_loads!(pm::AbstractPowerModel, power_surplus::Vector{Float64}, bus_to_col::Vector{Int64}; nw::Int64=nw_id_default)
    # Apply load demands
    for (_, load) in ref(pm, nw, :load)::Dict{Int64, Any}
        load_col = bus_to_col[load["load_bus"]::Int64]
        if load_col > 0
            power_surplus[load_col] -= load["pd"]::Float64
        end
    end
end

"subtracts the (static) shunt conductance power losses from `power_surplus` (P_shunt = gs * V^2, V=1.0 in DC)"
function _apply_shunts!(pm::AbstractPowerModel, power_surplus::Vector{Float64}, bus_to_col::Vector{Int64}; nw::Int64=nw_id_default)
    for (_, shunt) in ref(pm, nw, :shunt)::Dict{Int64, Any}
        shunt_bus_col = bus_to_col[shunt["shunt_bus"]::Int64]
        if shunt_bus_col > 0
            power_surplus[shunt_bus_col] -= shunt["gs"]::Float64
        end
    end
end

"applies or removes values of the current master solution to `A` & `power_surplus`"
function _apply_master_vals!(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, pm::AbstractPowerModel, power_surplus::Vector{Float64}, bus_to_col::Vector{Int64}, master_var_vals, apply_func_z_branch::Function, apply_func_pg::Function; nw::Int64=nw_id_default, gens=true)::Tuple{SparseArrays.SparseMatrixCSC{Float64, Int64}, Vector{Float64}}
    # gens
    gens && _apply_master_vals_gen!(pm, power_surplus, bus_to_col, master_var_vals, apply_func_pg::Function, nw=nw)

    # branches
    for branch_id in axes(master_var_vals[:z_branch],1)
        if master_var_vals[:z_branch][branch_id] < 0.5
            _add_or_remove_branch!(A, pm, branch_id, bus_to_col, apply_func_z_branch, nw=nw)
        end
    end

    _apply_master_trafo_vals!(power_surplus, pm, bus_to_col, master_var_vals, apply_func_z_branch, apply_func_z_branch, nw=nw)
    return A, power_surplus
end

"applies or removes values of the current master solution to `A` & `power_surplus`"
function _apply_master_vals!(L::SparseArrays.SparseMatrixCSC{Float64, Int64}, U::SparseArrays.SparseMatrixCSC{Float64, Int64}, pm::AbstractPowerModel, power_surplus::Vector{Float64}, bus_to_col::Vector{Int64}, master_var_vals, apply_func_z_branch::Function, apply_func_pg::Function, row_perm::Vector{Int64}, col_perm::Vector{Int64}, row_scaling::Vector{Float64}; nw::Int64=nw_id_default)::Tuple{SparseArrays.SparseMatrixCSC{Float64, Int64}, SparseArrays.SparseMatrixCSC{Float64, Int64}, Vector{Float64}}
    # gens
    _apply_master_vals_gen!(pm, power_surplus, bus_to_col, master_var_vals, apply_func_pg::Function, nw=nw)

    # branches
    for branch_id in axes(master_var_vals[:z_branch],1)
        if master_var_vals[:z_branch][branch_id] < 0.5
            _add_or_remove_branch!(L, U, pm, branch_id, bus_to_col, apply_func_z_branch, row_perm, col_perm, row_scaling, nw=nw)
        end
    end

    _apply_master_trafo_vals!(power_surplus, pm, bus_to_col, master_var_vals, apply_func_z_branch, apply_func_z_branch, nw=nw)
    return L, U, power_surplus
end

"add master values for trafos to solver (by adding them to solver.tr_steps) for application later on"
function add_master_trafo_vals!(solver::AbstractSubproblemSolver, master_var_vals)
    if PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_SHIFTING]
        # apply trafo shifts to solver
        for branch_id in axes(master_var_vals[:tr_shift], 1)
            solver.tr_steps[branch_id] = master_var_vals[:tr_shift][branch_id]
        end
    end
end

######################################################################################################################################################
"applies the change in trafo `shift` value of `branch` to `power_surplus`"
function _apply_trafo_val!(power_surplus::Vector{Float64}, branch::Dict{String, Any}, bus_to_col::Vector{Int64}, shift::Float64, f_func::Function, t_func::Function)
    
    f_col = bus_to_col[branch["f_bus"]]
    t_col = bus_to_col[branch["t_bus"]]
    
    f_col <= 0 && return
    t_col <= 0 && return

    b_shift = shift * branch["b"]
    
    power_surplus[f_col] = f_func(power_surplus[f_col], b_shift)
    power_surplus[t_col] = t_func(power_surplus[t_col], b_shift)
end
function _apply_trafo_val!(power_surplus::Vector{Float64}, branch::Dict{Symbol, Any}, bus_to_col::Vector{Int64}, shift::Float64, f_func::Function, t_func::Function)
    f_col = bus_to_col[branch[:f_bus]]
    t_col = bus_to_col[branch[:t_bus]]
    
    f_col <= 0 && return
    t_col <= 0 && return

    b_shift = shift * branch[:b]
    
    power_surplus[f_col] = f_func(power_surplus[f_col], b_shift)
    power_surplus[t_col] = t_func(power_surplus[t_col], b_shift)
end
"applies the trafo `shift` value of `branch` to `power_surplus`"
function apply_trafo_val!(power_surplus::Vector{Float64}, branch::Dict, bus_to_col::Vector{Int64}, shift::Float64)
    _apply_trafo_val!(power_surplus, branch, bus_to_col, shift, -, +) # -b*(va_fr - va_to - shift) = 0 => -b*va_fr + b*va_to = -b*shift    or     -b*va_to + b*va_fr = b*shift
end
"resets the trafo `shift` value of `branch` in `power_surplus`"
function reset_trafo_val!(power_surplus::Vector{Float64}, branch::Dict, bus_to_col::Vector{Int64}, shift::Float64)
    _apply_trafo_val!(power_surplus, branch, bus_to_col, shift, +, -)
end
######################################################################################################################################################

"applies the change in the master var trafo values to the solver's LEQ (by applying the shifts to its `power_surplus`)"
function _apply_master_trafo_vals!(power_surplus::Vector{Float64}, pm::AbstractPowerModel, bus_to_col::Vector{Int64}, master_var_vals, f_func::Function, t_func::Function; nw::Int64=nw_id_default)
    #if PARAMS[:LIMITS][:TRANSFORMER][:TRANSFORMER_SHIFTING]
        for (_, branch) in ref(pm, nw, :branch)::Dict{Int64, Any}
            
            # if the branch is electrical and on, get its shift
            if is_branch_electrical(branch) && (is_branch_switchable(branch) ? isapprox(master_var_vals[:z_branch][branch["index"]], 1.0, atol=SOLVER_PRECISION) : branch["status"] == 1)
                shift = get_transformer_shift_master(branch, master_var_vals)

                if shift != 0.0
                    _apply_trafo_val!(power_surplus, branch, bus_to_col, shift, f_func, t_func)
                end
            end

        end
end
"applies the master var trafo values to the solver's LEQ (by applying the shifts to its `power_surplus`)"
function apply_master_trafo_vals!(power_surplus::Vector{Float64}, pm::AbstractPowerModel, bus_to_col::Vector{Int64}, master_var_vals; nw::Int64=nw_id_default)
    _apply_master_trafo_vals!(power_surplus, pm, bus_to_col, master_var_vals, -, +, nw=nw)
end
"resets the master var trafo values to the solver's LEQ (by applying the shifts to its `power_surplus`)"
function reset_master_trafo_vals!(power_surplus::Vector{Float64}, pm::AbstractPowerModel, bus_to_col::Vector{Int64}, master_var_vals; nw::Int64=nw_id_default)
    _apply_master_trafo_vals!(power_surplus, pm, bus_to_col, master_var_vals, +, -, nw=nw)
end

"applies shifts of transformers that can not be shifted, ie either TRANSFORMER_SHIFTING false or the transformer can not be shifted at all. If `only_active`, currently disabled branches are ignored."
function apply_fixed_trafo_vals!(power_surplus::Vector{Float64}, pm::AbstractPowerModel, bus_to_col::Vector{Int64}; nw::Int64=nw_id_default, only_active::Bool=true)
    for (branch_id, branch) in ref(pm, nw, :branch)
        # if not a transformer, or shiftable and we use shift variables => we don't use a fixed value, continue
        if is_branch_shiftable(branch)
            continue
        end

        # skip disabled branches
        if only_active && get_status(branch) == 0
            continue
        end

        if is_branch_electrical(branch)
            # apply shift of branch to power_surplus
            apply_trafo_val!(power_surplus, branch, bus_to_col, get_transformer_shift(pm, nw, branch_id)::Float64)
        end
    end
end

"applies values of the current master solution to `A` & `power_surplus`"
function apply_master_vals!(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, pm::AbstractPowerModel, power_surplus::Vector{Float64}, bus_to_col::Vector{Int64}, master_var_vals; nw::Int64=nw_id_default, gens=true)::Tuple{SparseArrays.SparseMatrixCSC{Float64, Int64}, Vector{Float64}}
    _apply_master_vals!(A, pm, power_surplus, bus_to_col, master_var_vals, -, +, nw=nw, gens=gens)
    apply_master_trafo_vals!(power_surplus, pm, bus_to_col, master_var_vals, nw=nw)

    return A, power_surplus
end

"resets values of the current master solution applied to `A` & `power_surplus`"
function reset_master_vals!(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, pm::AbstractPowerModel, power_surplus::Vector{Float64}, bus_to_col::Vector{Int64}, master_var_vals; nw::Int64=nw_id_default, gens=true)::Tuple{SparseArrays.SparseMatrixCSC{Float64, Int64}, Vector{Float64}}
    _apply_master_vals!(A, pm, power_surplus, bus_to_col, master_var_vals, +, -, nw=nw, gens=gens)
    reset_master_trafo_vals!(power_surplus, pm, bus_to_col, master_var_vals, nw=nw)

    return A, power_surplus
end

"applies values of the current master solution to the factorization `LU` & `power_surplus`"
function apply_master_vals!(L::SparseArrays.SparseMatrixCSC{Float64, Int64}, U::SparseArrays.SparseMatrixCSC{Float64, Int64}, pm::AbstractPowerModel, power_surplus::Vector{Float64}, bus_to_col::Vector{Int64}, master_var_vals, row_perm::Vector{Int64}, col_perm::Vector{Int64}, row_scaling::Vector{Float64}; nw::Int64=nw_id_default, gens=true)::Tuple{SparseArrays.SparseMatrixCSC{Float64, Int64}, SparseArrays.SparseMatrixCSC{Float64, Int64}, Vector{Float64}}
    _apply_master_vals!(L, U, pm, power_surplus, bus_to_col, master_var_vals, -, +, row_perm, col_perm, row_scaling, nw=nw)
    apply_master_trafo_vals!(power_surplus, pm, bus_to_col, master_var_vals, nw=nw)

    return L, U, power_surplus
end

"resets values of the current master solution applied to to the factorization `LU` & `power_surplus`"
function reset_master_vals!(L::SparseArrays.SparseMatrixCSC{Float64, Int64}, U::SparseArrays.SparseMatrixCSC{Float64, Int64}, pm::AbstractPowerModel, power_surplus::Vector{Float64}, bus_to_col::Vector{Int64}, master_var_vals, row_perm::Vector{Int64}, col_perm::Vector{Int64}, row_scaling::Vector{Float64}; nw::Int64=nw_id_default, gens=true)::Tuple{SparseArrays.SparseMatrixCSC{Float64, Int64}, SparseArrays.SparseMatrixCSC{Float64, Int64}, Vector{Float64}}
    _apply_master_vals!(L, U, pm, power_surplus, bus_to_col, master_var_vals, +, -, row_perm, col_perm, row_scaling, nw=nw)
    reset_master_trafo_vals!(power_surplus, pm, bus_to_col, master_var_vals, nw=nw)

    return L, U, power_surplus
end

### components ###
"returns lines attached to the component"
function _get_component_lines(component::BBComponent, pm, nw::Int64=nw_id_default)::Vector{Int64}
    return [branch_id for (branch_id, branch) in ref(pm, nw, :branch) if branch["type"] ∉ ["K", "breaker"] && (branch["f_bus"] ∈ component.busbar_buses || branch["t_bus"] ∈ component.busbar_buses)]
end

"calculates connected components within each substation"
function calc_full_busbar_components(pm, off_branches::Vector{Int64}; nw::Int64=nw_id_default)::Tuple{Vector{BBComponent}, Vector{BBComponent}, Dict{Int64, BBComponent}}
    if haskey(pm.data, "nw") && haskey(pm.data["nw"]["1"], "substations")
        substations = pm.data["nw"]["1"]["substations"]
    else
        substations = pm.data["substations"]
    end
    
    components = Vector{BBComponent}() # components
    detached_components = Vector{BBComponent}()
    branch_to_comp = Dict{Int64, BBComponent}() # branch => component that contains it

    for (_, substation) in substations
        breaker_buses = [d["bus_i"] for d in substation["aux_bus"]]

        # open_branches contains all breaker branches and coupler branches (no lines) that have not been added to a component yet
        open_branches::Vector{Int64} = [parse(Int64, branch_id) for branch_id in substation["K"] if parse(Int64, branch_id) ∉ off_branches] # unvisited couplers
        for aux_bus in substation["aux_bus"] # unvisited breakers
            for branch_id in aux_bus["breaker_branches"]
                if branch_id ∉ off_branches
                    push!(open_branches, branch_id)
                end
            end
        end
        
        # each component contains buses of the substation that are connected
        if isempty(open_branches) # no branches ⇒ no connection, so singleton component
            for bus_id in parse.(Int64, substation["busbar"])
                
                component = BBComponent([bus_id],[],[],[],[])
                
                component.lines = _get_component_lines(component, pm, nw)
                
                # detached if no lines (we know there are no breaker branches)
                if isempty(component.lines)
                    push!(detached_components, component)
                else
                    push!(components, component)
                end

            end

            continue
        end

        # fill up each component until it is complete, then proceed to the next one
        all_visited_buses = Vector{Int64}()
        while !isempty(open_branches)
            # initialize the component with an arbitrary branch
            branch_id = open_branches[1]
            branch = ref(pm, nw, :branch, branch_id)

            visited_buses = [branch["f_bus"], branch["t_bus"]] # buses visited in the loop
            visited_branches = [branch_id] # branches visited in the loop
            
            open_branches = open_branches[2:end]
            
            
            # while there are updates to the current component, continue to look for branches incident to it
            new = true
            while new
                new = false
                # check if the branch is incident, and if it introduces a new vertex into the component
                for branch_id in open_branches
                    branch = ref(pm, nw, :branch, branch_id)
                    
                    if branch["f_bus"] ∈ visited_buses
                        if branch["t_bus"] ∉ visited_buses # both are contained ⇒ connected, but no update; only f_bus incident ⇒ update
                            push!(visited_buses, branch["t_bus"])
                            
                            new = true
                        end
                        
                        push!(visited_branches, branch_id)
                    elseif branch["t_bus"] ∈ visited_buses # only t_bus incident ⇒ update
                        push!(visited_buses, branch["f_bus"])
                        
                        push!(visited_branches, branch_id)
                        new = true
                    end
                end
    
                # remove visited branches from open_branches
                open_branches = [branch_id for branch_id in open_branches if branch_id ∉ visited_branches]
            end

            # construct the BBComponent: busbars, breaker buses, couplers, breaker branches
            component = BBComponent([bus_id for bus_id in visited_buses if "$bus_id" ∈ substation["busbar"]],
                                    [bus_id for bus_id in visited_buses if bus_id ∈ breaker_buses],
                                    [branch_id for branch_id in visited_branches if "$branch_id" ∈ substation["K"]],
                                    [branch_id for branch_id in visited_branches if "$branch_id" ∉ substation["K"]],
                                    [])
            component.lines = _get_component_lines(component, pm, nw)


            @assert(length(visited_buses) == length([component.breaker_buses; component.busbar_buses]), "visited buses: $visited_buses, breaker buses: $(component.breaker_buses), busbar buses: $(component.busbar_buses), branch: $branch")
            append!(all_visited_buses, visited_buses)

            # add to (detached) components
            if isempty(component.lines) && isempty(component.breaker_buses)
                push!(detached_components, component)
            else
                push!(components, component)
            end

            # fill entries in branch_to_comp
            for branch_id in [visited_branches; component.lines]
                branch_to_comp[branch_id] = component
            end
        end

        # buses without connected branches form detached singleton components
         for busbar_bus_id in parse.(Int64, substation["busbar"])
            if busbar_bus_id ∉ all_visited_buses
                push!(detached_components, BBComponent([busbar_bus_id],[],[],[],[]))
            end
        end

        # assign each detached breaker bus component idx -1
        for breaker_bus_id in breaker_buses
            if breaker_bus_id ∉ all_visited_buses
                bus_to_col[breaker_bus_id] = typemin(Int64)
            end
        end
    end

    return (components, detached_components, branch_to_comp)
end

"adds a new breaker branch to the component, where the busbar is already part of the component"
function add_breaker_branch_to_component!(component::BBComponent, branch_id::Int64, breaker_bus_id::Int64, branch_to_comp::Dict{Int64, BBComponent})
    push!(component.breaker_branches, branch_id)
    push!(component.breaker_buses,    breaker_bus_id)

    branch_to_comp[branch_id] = component
end
function add_breaker_branch_to_component!(component::BBComponent, component_col::Int64, branch_id::Int64, breaker_bus_id::Int64, bus_to_col::Vector{Int64}, branch_to_comp::Dict{Int64, BBComponent})
    add_breaker_branch_to_component!(component, branch_id, breaker_bus_id, branch_to_comp)

    bus_to_col[breaker_bus_id] = component_col
end

"removes a breaker branch from the component"
function remove_breaker_branch_from_component!(component::BBComponent, branch_id::Int64, breaker_bus_id::Int64)
    filter!(x -> x != branch_id, component.breaker_branches)
    filter!(x -> x != breaker_bus_id, component.breaker_buses)
end

"adds a new line to the component"
function add_line_to_component!(component::BBComponent, branch_id::Int64, branch_to_comp::Dict{Int64, BBComponent})
    push!(component.lines, branch_id)

    branch_to_comp[branch_id] = component
end

"removes a line from the component"
function remove_line_from_component!(component::BBComponent, branch_id::Int64)
    filter!(x -> x != branch_id, component.lines)
end

"increment the position of all following (i.e. with lower position) detached components and delete in detached_components"
function delete_from_detached_components(component_col::Int64, bus_to_col::Vector{Int64}, detached_components::Vector{BBComponent})
    for (bus_id, other_col) in enumerate(bus_to_col)
        if other_col < component_col
            bus_to_col[bus_id] = other_col + 1
        end
    end

    deleteat!(detached_components, -component_col)
end

"decrement the position of all following (i.e. with higher position) components and delete in components"
function delete_from_components(component_col::Int64, bus_to_col::Vector{Int64}, components::Vector{BBComponent})
    for (bus_id, other_col) in enumerate(bus_to_col)
        if other_col > component_col
            bus_to_col[bus_id] = other_col - 1
        end
    end

    deleteat!(components, component_col)
end

"merges componets and optionally enacts changes on the matrix"
function merge_components(scol::Int64, scomponent::BBComponent, lcol::Int64, lcomponent::BBComponent, components::Vector{BBComponent}, detached_components::Vector{BBComponent}, bus_to_col::Vector{Int64}, branch_to_comp::Dict{Int64, BBComponent}, A::SparseArrays.SparseMatrixCSC{Float64, Int64}, pm::AbstractPowerModel; nw::Int64=nw_id_default, is_scol_detached::Bool=false, is_lcol_detached::Bool=false)::SparseArrays.SparseMatrixCSC{Float64, Int64}
    
    # append all vectors of lcomponent to scomponent
    append!(scomponent.busbar_buses,     lcomponent.busbar_buses)
    append!(scomponent.breaker_buses,    lcomponent.breaker_buses)
    append!(scomponent.coupler_branches, lcomponent.coupler_branches)
    append!(scomponent.breaker_branches, lcomponent.breaker_branches)
    append!(scomponent.lines,            lcomponent.lines)

    # re-link all buses & branches associated with lcomponent to scomponent
    # we need to remove the previous connections before re-mapping bus_to_col
    # TODO the removal should already be handled by deleting row & col?
    if !is_lcol_detached
        for branch_id in lcomponent.breaker_branches
            remove_orig_branch!(A, pm, branch_id, bus_to_col, nw=nw) # remove the actual branch, not this auxiliary branch
        end
        for branch_id in lcomponent.lines
            remove_branch!(A, pm, branch_id, bus_to_col, nw=nw)
        end

        # delete lcomponent from components, and decrement all links to the following components
        delete_from_components(lcol, bus_to_col, components)

        # delete column & row in A
        #           ...lcol-1|lcol|lcol+1...
        # ...lcol-1      ⋅   |  X |   ⋅
        #      lcol      X   |  X |   X
        # lcol+1...      ⋅   |  X |   ⋅
        A = delete_matrix_row_and_col(A, lcol)
    else
        delete_from_detached_components(lcol, bus_to_col, detached_components)
    end
    
    for branch_id in lcomponent.coupler_branches
        branch_to_comp[branch_id] = scomponent
    end

    for bus_id in lcomponent.busbar_buses
        bus_to_col[bus_id] = scol
    end

    if !is_scol_detached
        for bus_id in lcomponent.breaker_buses
            bus_to_col[bus_id] = scol
        end
    
        for branch_id in lcomponent.breaker_branches
            branch_to_comp[branch_id] = scomponent

            add_orig_branch!(A, pm, branch_id, bus_to_col, nw=nw) # add the actual branch, not this auxiliary branch
        end
        for branch_id in lcomponent.lines
            branch_to_comp[branch_id] = scomponent
            
            add_branch!(A, pm, branch_id,            bus_to_col, nw=nw)
        end
    end

    return A
end

"applies the new component to all relevant data structures"
function add_component_to_matrix(component::BBComponent, components::Vector{BBComponent}, detached_components::Vector{BBComponent}, bus_to_col::Vector{Int64}, branch_to_comp::Dict{Int64, BBComponent}, A::SparseArrays.SparseMatrixCSC{Float64, Int64}, pm::AbstractPowerModel; nw::Int64=nw_id_default, is_detached::Bool=false)::SparseArrays.SparseMatrixCSC{Float64, Int64}
    if is_detached
        col = -length(detached_components)

        for bus_id in component.busbar_buses
            bus_to_col[bus_id] = col
        end
    else
        col = length(components)

        for bus_id in component.busbar_buses
            bus_to_col[bus_id] = col
        end
        for bus_id in component.breaker_buses
            bus_to_col[bus_id] = col
        end

        # add an empty row & column to A for the new component, then add the branches
        A = resize_sparse(A, maximum(bus_to_col), col)

        for branch_id in component.breaker_branches
            branch_to_comp[branch_id] = component
            
            add_orig_branch!(A, pm, branch_id, bus_to_col, nw=nw) # add the actual branch, not this auxiliary branch
        end
        for branch_id in component.lines
            branch_to_comp[branch_id] = component
            branch = ref(pm, nw, :branch, branch_id)
            
            add_branch!(A, pm, branch_id,            bus_to_col, nw=nw)
        end
    end

    for branch_id in component.coupler_branches
        branch_to_comp[branch_id] = component
    end

    return A
end

"applies values of the current master solution to `A`, `power_surplus`, `componentss` and `branch_to_comp`, set reverse=true to revert previously applied changes"
function _apply_master_vals!(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, pm::AbstractPowerModel, power_surplus::Vector{Float64}, bus_to_col::Vector{Int64}, components::Vector{BBComponent}, detached_components::Vector{BBComponent}, branch_to_comp::Dict{Int64, BBComponent}, master_var_vals, reverse::Bool; nw::Int64=nw_id_default, gens=true)::Tuple{SparseArrays.SparseMatrixCSC{Float64, Int64}, Vector{Float64}}
    # (directly attached lines are handled analogously to breaker branches in the following)
    # all buses of the same component share one column in A, we need to keep track of it (most of what we do here is there to facilitate this)
    # the ordering of operations here is important, since by altering the components, we lose information about where a given non-bb branch might have previously been attached. We need to
    # 1. remove breaker branches - otherwise we lose the information
    # 2. add couplers - doing this before removing potentially reduces the number of splits
    # 3. remove couplers
    # 4. add breaker branches - doing this last ensures that these branches don't need to be re-attached
    # whenever a "real" branch is (dis)connected, we update A, also if merging/splitting components
    # generation & load values are applied via the updated bus_to_col at the end
    master_var_vals_z_branch = master_var_vals[:z_branch]::JuMP.Containers.DenseAxisArray{Float64, 1, Tuple{Vector{Int64}}, Tuple{JuMP.Containers._AxisLookup{Dict{Int64, Int64}}}}

    for phase in 1:4
        for (_, branch) in ref(pm, nw, :branch)::Dict{Int64, Any}#axes(master_var_vals[:z_branch],1)
            if ! is_branch_switchable(branch) # if this branch is not switchable, we can not change it
                continue
            end

            branch_id = branch["index"]::Int64
            cur_status = round(master_var_vals_z_branch[branch_id]::Float64)
            #branch = ref(pm, nw, :branch, branch_id)

            status = branch["status"]::Int64
            if status != cur_status # if this changed from the starting topology...
                
                
                if reverse
                    cur_status = 1 - cur_status # add all 0 branches, remove all 1 branches
                end

                if branch["type"]  == "breaker"
                    if get(ref(pm, nw, :bus, branch["f_bus"]), "type", "") == "breaker" # no "type" => not an auxiliary bus
                        breaker_bus_id = branch["f_bus"]
                        busbar_bus_id  = branch["t_bus"]
                    else
                        breaker_bus_id = branch["t_bus"]
                        busbar_bus_id  = branch["f_bus"]
                    end

                    if phase == 1 && isapprox(cur_status, 0, atol=SOLVER_PRECISION)

                        # remove branch from A, then remove connection to branch by removing bus_to_col entry & remove from component
                        remove_orig_branch!(A, pm, branch_id, bus_to_col, nw=nw) # remove the actual branch, not this auxiliary branch, but only if its other end is connected to a component (i.e. the branch is not already removed)
                        
                        component = branch_to_comp[branch_id]

                        col = bus_to_col[busbar_bus_id]
                        bus_to_col[breaker_bus_id] = 0
                        
                        remove_breaker_branch_from_component!(component, branch_id, breaker_bus_id)
                        
                        col = bus_to_col[busbar_bus_id]
                        if isempty(component.breaker_buses)
                            @assert isempty(component.breaker_branches)
                            push!(detached_components, component)
                            
                            delete_from_components(col, bus_to_col, components)
                            A = delete_matrix_row_and_col(A, col)

                            A = add_component_to_matrix(component, components, detached_components, bus_to_col, branch_to_comp, A, pm; nw=nw, is_detached=true)
                        end


                    elseif phase == 4 && isapprox(cur_status, 1, atol=SOLVER_PRECISION)

                        # add connection to branch by mapping in bus_to_col & add to component, then add the branch to A
                        col = bus_to_col[busbar_bus_id]

                        if col > 0 # was already attached
                            component = components[col]
                        else
                            component = detached_components[-col]

                            push!(components, component)

                            delete_from_detached_components(col, bus_to_col, detached_components)

                            col = length(components)

                            # add branches back to matrix since they are not detached anymore
                            A = add_component_to_matrix(component, components, detached_components, bus_to_col, branch_to_comp, A, pm; nw=nw, is_detached=false)
                        end
                        
                        add_breaker_branch_to_component!(component, col, branch_id, breaker_bus_id, bus_to_col, branch_to_comp)

                        add_orig_branch!(A, pm, branch_id, bus_to_col, nw=nw) # add the actual branch, but only if its other end is connected to a component (i.e. the branch can actually be added)
                    end

                elseif branch["type"] == "K" # coupler
                    fbus_id = branch["f_bus"]
                    tbus_id = branch["t_bus"]

                    # _s_maller _col_umn, _l_arger _col_umn (index wise)
                    fcol, tcol = bus_to_col[fbus_id], bus_to_col[tbus_id]
                    s_id, l_id = (abs(fcol) <= abs(tcol) && !(fcol < 0 && tcol > 0)) || (fcol > 0 && tcol < 0) ? (fbus_id, tbus_id) : (tbus_id, fbus_id)
                    scol, lcol = bus_to_col[s_id], bus_to_col[l_id]
                    scomponent = scol > 0 ? components[scol] : detached_components[-scol]

                    if phase == 2 && isapprox(cur_status, 1, atol=SOLVER_PRECISION) # potential merge    

                        if scol != lcol # merge: lcomponent into scomponent (arbitrary)
                            if scol ≤ 0 && lcol ≤ 0 # both components are detached -> create a new component but don't add it to the matrix, as it is not connected to a breaker
                                @assert scol < 0 && lcol < 0

                                # already assigned: scomponent = detached_components[-scol]
                                lcomponent = detached_components[-lcol]

                                A = merge_components(scol, scomponent, lcol, lcomponent, components, detached_components, bus_to_col, branch_to_comp, A, pm; nw=nw, is_scol_detached=true, is_lcol_detached=true)

                            elseif scol ≤ 0 || lcol ≤ 0 # one of the components is detached -> add it to the other's component
                                if scol ≤ 0

                                    attached_col = lcol
                                    attached_comp = components[lcol]
                                    detached_col = scol
                                    detached_comp = detached_components[-scol]

                                elseif lcol ≤ 0

                                    attached_col = scol
                                    attached_comp = components[scol]
                                    detached_col = lcol
                                    detached_comp = detached_components[-lcol]

                                end
                                
                                scomponent = attached_comp
                                A = merge_components(attached_col, attached_comp, detached_col, detached_comp, components, detached_components, bus_to_col, branch_to_comp, A, pm; nw=nw, is_scol_detached=false, is_lcol_detached=true)
                            else
                                # already assigned: scomponent = components[scol]
                                lcomponent = components[lcol]

                                A = merge_components(scol, scomponent, lcol, lcomponent, components, detached_components, bus_to_col, branch_to_comp, A, pm; nw=nw, is_scol_detached=false, is_lcol_detached=false)
                            end

                        end

                        # add to the (now) single associated component
                        push!(scomponent.coupler_branches, branch_id)

                    elseif phase == 3 && isapprox(cur_status, 0, atol=SOLVER_PRECISION) # potential split
                        scol > 0 ? scomponent = components[scol] : scomponent = detached_components[-scol]

                        # fbus_id and tbus_id are in the same component, we need to check whether disconnecting this branch splits them
                        filter!(e->e≠branch_id, scomponent.coupler_branches)

                        # this will result in at most two (ie one additional) components
                        bbars1 = [fbus_id]
                        bbars2 = [tbus_id]
                        couplers1 = Vector{Int64}()
                        couplers2 = Vector{Int64}()
                        
                        # for each branch_id, check whether it can be attached to one of the neighboring buses
                        is_connected = false
                        open_branches = deepcopy(scomponent.coupler_branches)
                        
                        while !isempty(open_branches)

                            visited_branches = Vector{Int64}()
                            for branch_id in open_branches
                                branch = ref(pm, nw, :branch, branch_id)

                                if branch["f_bus"] ∈ bbars1
                                    if branch["t_bus"] ∈ bbars2 # connecting branch

                                        is_connected = true
                                        break
                                    else # add to component 1

                                        if branch["t_bus"] ∉ bbars1
                                            push!(bbars1, branch["t_bus"])
                                        end
                                        push!(couplers1, branch_id)
                                        push!(visited_branches, branch_id)
                                    end
                                elseif branch["f_bus"] ∈ bbars2
                                    if branch["t_bus"] ∈ bbars1 # connecting branch

                                        is_connected = true
                                        break
                                    else # add to component 2
                                        if branch["t_bus"] ∉ bbars2
                                            push!(bbars2, branch["t_bus"])
                                        end
                                        push!(couplers2, branch_id)
                                        push!(visited_branches, branch_id)
                                    end
                                else # try to attach the tbus_id
                                    if branch["t_bus"] ∈ bbars1
                                        if branch["f_bus"] ∉ bbars1
                                            push!(bbars1, branch["f_bus"])
                                        end
                                        push!(couplers1, branch_id)
                                        push!(visited_branches, branch_id)
                                    elseif branch["t_bus"] ∈ bbars2
                                        if branch["f_bus"] ∉ bbars2
                                            push!(bbars2, branch["f_bus"])
                                        end
                                        push!(couplers2, branch_id)
                                        push!(visited_branches, branch_id)
                                    end
                                end

                            end

                            # stop if the component is still connected
                            if is_connected
                                break
                            end

                            filter!(e->e∉visited_branches, open_branches)
                        end

                        # split the component if it is actually being split
                        if !is_connected

                            # shrink the remaining parts of scomponent
                            scomponent.busbar_buses     = bbars1
                            scomponent.coupler_branches = couplers1

                            # build the new component
                            new_component = BBComponent(bbars2, Vector{Int64}(), couplers2, Vector{Int64}(), [])

                            # assign breaker buses & branches
                            for branch_id in deepcopy(scomponent.breaker_branches)
                                branch = ref(pm, nw, :branch, branch_id)

                                if branch["f_bus"]     ∈ new_component.busbar_buses
                                    breaker_bus_id = branch["t_bus"]
                                elseif branch["t_bus"] ∈ new_component.busbar_buses
                                    breaker_bus_id = branch["f_bus"]
                                else # part of scomponent
                                    continue
                                end

                                add_breaker_branch_to_component!(new_component, branch_id, breaker_bus_id, branch_to_comp)
                                remove_breaker_branch_from_component!(scomponent, branch_id, breaker_bus_id)
                            end
                            # assign lines
                            for branch_id in deepcopy(scomponent.lines)
                                branch = ref(pm, nw, :branch, branch_id)

                                # do nothing if this line is not connected to one of the busbars of this new component
                                if !(branch["f_bus"] ∈ new_component.busbar_buses || branch["t_bus"] ∈ new_component.busbar_buses)
                                    continue
                                end

                                add_line_to_component!(new_component, branch_id, branch_to_comp)
                                remove_line_from_component!(scomponent, branch_id)
                            end

                            
                            # add to components
                            is_s_detached = isempty(scomponent.breaker_branches) && isempty(scomponent.lines)
                            is_new_detached = isempty(new_component.breaker_branches) && isempty(new_component.lines)

                            # if scomponent is now detached but new_component not, we swap them so we can add the new component to detached_components
                            if is_s_detached && !is_new_detached
                                tmp_component = scomponent
                                scomponent = new_component
                                new_component = tmp_component

                                components[scol] = scomponent # the entry still points to the previous scomponent

                                is_s_detached = false
                                is_new_detached = true
                            end


                            # scomponent is already in the correct vector, new_component is either detached or not
                            push!(is_new_detached ? detached_components : components, new_component)

                            ## re-link all buses & branches associated with scomponent to new_component
                            # we need to remove the previous connections before re-mapping bus_to_col
                            if scol > 0 # only if merged component was attached
                                for branch_id in new_component.breaker_branches
                                    remove_orig_branch!(A, pm, branch_id, bus_to_col, nw=nw) # remove the actual branch, not this auxiliary branch
                                end
                                for branch_id in new_component.lines
                                    remove_branch!(A, pm, branch_id,            bus_to_col, nw=nw)
                                end

                                # we need to increment the counter of the following cols to ensure that the first length(components) columns correspond to the indices of `components`
                                if !is_new_detached
                                    for (bus_id, col) in enumerate(bus_to_col)
                                        if col >= length(components)
                                            bus_to_col[bus_id] = col+1
                                        end
                                    end
                                end
                            end

                            # map buses onto new component, expand matrix if attached
                            A = add_component_to_matrix(new_component, components, detached_components, bus_to_col, branch_to_comp, A, pm; nw=nw, is_detached=is_new_detached)
                            
                        end

                    end
                else # directly attached lines - no changes in bus_to_col
                    if phase == 1 && isapprox(cur_status, 0, atol=SOLVER_PRECISION)
                        # remove branch from A & remove from component
                        remove_branch!(A, pm, branch_id, bus_to_col, nw=nw)

                    elseif phase == 4 && isapprox(cur_status, 1, atol=SOLVER_PRECISION)
                        # add the branch to A & point to in component
                        add_branch!(A, pm, branch_id, bus_to_col, nw=nw)
                    end

                end
            end
        end
        #for component ∈ components
        #    @assert(length(component.busbar_buses) == length(unique(component.busbar_buses)), "Duplicate busbars in components: $(component.busbar_buses) in phase $phase")
        #end
    end

    # replace power_surplus by load values following the current component state, then apply the master generator values
    power_surplus = zeros(maximum(bus_to_col))
    _apply_loads!(pm, power_surplus, bus_to_col, nw=nw)
    _apply_shunts!(pm, power_surplus, bus_to_col, nw=nw)
    if gens# && !reverse
        _apply_master_vals_gen!(pm, power_surplus, bus_to_col, master_var_vals, reverse ? (-) : (+), nw=nw)
    end
    _apply_master_trafo_vals!(power_surplus, pm, bus_to_col, master_var_vals, reverse ? (+) : (-), reverse ? (-) : (+), nw=nw)
    
    return A, power_surplus
end

"applies values of the current master solution to `A`, `power_surplus`, `componentss` and `branch_to_comp`. Important: `A` and `power_surplus` need to be overwritten by this funcions output arguments, otherwise they will not be altered!"
function apply_master_vals!(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, pm::AbstractPowerModel, power_surplus::Vector{Float64}, bus_to_col::Vector{Int64}, components::Vector{BBComponent}, detached_components::Vector{BBComponent}, branch_to_comp::Dict{Int64, BBComponent}, master_var_vals; nw::Int64=nw_id_default, gens=true)
    _apply_master_trafo_vals!(power_surplus, pm, bus_to_col, master_var_vals, -, +, nw=nw)
    return _apply_master_vals!(A, pm, power_surplus, bus_to_col, components, detached_components, branch_to_comp, master_var_vals, false, nw=nw, gens=gens)::Tuple{SparseArrays.SparseMatrixCSC{Float64, Int64}, Vector{Float64}}
end

"resets values of the current master solution to `A`, `power_surplus`, `componentss` and `branch_to_comp`. Important: `A` and `power_surplus` need to be overwritten by this funcions output arguments, otherwise they will not be altered!"
function reset_master_vals!(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, pm::AbstractPowerModel, power_surplus::Vector{Float64}, bus_to_col::Vector{Int64}, components::Vector{BBComponent}, detached_components::Vector{BBComponent}, branch_to_comp::Dict{Int64, BBComponent}, master_var_vals; nw::Int64=nw_id_default, gens=true)
    return _apply_master_vals!(A, pm, power_surplus, bus_to_col, components, detached_components, branch_to_comp, master_var_vals, true, nw=nw, gens=gens)::Tuple{SparseArrays.SparseMatrixCSC{Float64, Int64}, Vector{Float64}}
end