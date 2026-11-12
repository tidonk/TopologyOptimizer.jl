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

###### contains all additional types ######

##### Benders Subproblem Solver Types #####
##auxiliary types
@kwdef mutable struct BBComponent
  busbar_buses::Vector{Int64}        = Vector()
  breaker_buses::Vector{Int64}       = Vector()
  coupler_branches::Vector{Int64}    = Vector()
  breaker_branches::Vector{Int64}    = Vector()
  lines::Vector{Int64}               = Vector() # "real" branches directly attached to a busbar
end
##

abstract type AbstractSubproblemSolver end

@kwdef mutable struct FullLPSolver <: AbstractSubproblemSolver
  sub_cons_for_duals::Dict{Symbol, Dict{Int, JuMP.ConstraintRef}}
  pm                ::AbstractPowerModel
end

@kwdef mutable struct LPSolver <: AbstractSubproblemSolver
  pm         ::AbstractPowerModel
  nw         ::Int64
  branch_data::Dict{Int64, Dict{Symbol, Any}} = Dict()
end


abstract type AbstractMatrixSolver <: AbstractSubproblemSolver end # anything that solves A\b = vas
abstract type AbstractLinEqSysSolver <: AbstractMatrixSolver end # LinEqSysSolver and LinEqSysBBSolver
abstract type AbstractLinEqSysRfSolver <: AbstractMatrixSolver end # LinEqSysRfSolver and LinEqSysBBRfSolver
@kwdef mutable struct LinEqSysSolver <: AbstractLinEqSysSolver
    A              ::SparseArrays.SparseMatrixCSC{Float64, Int64}
    A_feas_cut     ::SparseArrays.SparseMatrixCSC{Float64, Int64}= deepcopy(A)
    power_surplus  ::Vector{Float64}
    bus_to_col     ::Vector{Int64}
    n_buses        ::Int64                                       = length(bus_to_col)
    vas            ::Vector{Float64}                             = Vector(undef, n_buses)
    dual_rhs       ::Vector{Float64}                             = Vector(undef, n_buses)
    duals          ::Vector{Float64}                             = Vector(undef, n_buses)
    fallback_solver::AbstractSubproblemSolver
    tr_steps       ::Dict{Float64, Float64}                      = Dict{Float64, Float64}()
    branch_data    ::Dict{Int64, Dict{Symbol, Any}}              = Dict()
end
@kwdef mutable struct LinEqSysRfSolver <: AbstractLinEqSysRfSolver # LinEqSysSolver with re-factorization
    bus_to_col     ::Vector{Int64}
    n_buses        ::Int64                                       = length(bus_to_col)
    A              ::SparseArrays.SparseMatrixCSC{Float64, Int64}
    feas_col       ::SparseArrays.SparseMatrixCSC{Float64, Int64}= deepcopy(A)
    feas_row       ::SparseArrays.SparseMatrixCSC{Float64, Int64}= feas_col
    feas_diag      ::Float64                                     = 0.0
    L              ::SparseArrays.SparseMatrixCSC{Float64, Int64}= A
    U              ::SparseArrays.SparseMatrixCSC{Float64, Int64}= A
    p              ::Vector{Int64}                               = Vector(undef, n_buses)
    q              ::Vector{Int64}                               = Vector(undef, n_buses)
    Rs             ::Vector{Float64}                             = Vector(undef, n_buses)
    is_factorized  ::Bool                                        = false
    power_surplus  ::Vector{Float64}
    linsolve       ::LinearSolve.LinearCache
    vas            ::Vector{Float64}                             = Vector(undef, n_buses)
    dual_rhs       ::Vector{Float64}                             = Vector(undef, n_buses)
    duals          ::Vector{Float64}                             = Vector(undef, n_buses)
    fallback_solver::AbstractSubproblemSolver
    tr_steps       ::Dict{Float64, Float64}                      = Dict{Float64, Float64}()
    branch_data    ::Dict{Int64, Dict{Symbol, Any}}              = Dict()
end

@kwdef mutable struct LinEqSysBBSolver <: AbstractLinEqSysSolver
    A              ::SparseArrays.SparseMatrixCSC{Float64, Int64}
    A_feas_cut     ::SparseArrays.SparseMatrixCSC{Float64, Int64}= deepcopy(A)
    power_surplus  ::Vector{Float64}
    bus_to_col     ::Vector{Int64}
    components     ::Vector{BBComponent}
    detached_components::Vector{BBComponent}
    branch_to_comp ::Dict{Int64, BBComponent}    
    n_buses        ::Int64                                       = length(bus_to_col)
    vas            ::Vector{Float64}                             = Vector(undef, n_buses)
    dual_rhs       ::Vector{Float64}                             = Vector(undef, n_buses)
    duals          ::Vector{Float64}                             = Vector(undef, n_buses)
    fallback_solver::AbstractSubproblemSolver
    tr_steps       ::Dict{Float64, Float64}                      = Dict{Float64, Float64}()
    branch_data    ::Dict{Int64, Dict{Symbol, Any}}              = Dict()
end
@kwdef mutable struct LinEqSysBBRfSolver <: AbstractLinEqSysRfSolver # LinEqSysBBSolver with re-factorization
    bus_to_col     ::Vector{Int64}
    n_buses        ::Int64                                       = length(bus_to_col)
    A              ::SparseArrays.SparseMatrixCSC{Float64, Int64}
    feas_col       ::SparseArrays.SparseMatrixCSC{Float64, Int64}= deepcopy(A)
    feas_row       ::SparseArrays.SparseMatrixCSC{Float64, Int64}= feas_col
    feas_diag      ::Float64                                     = 0.0
    L              ::SparseArrays.SparseMatrixCSC{Float64, Int64}= A
    U              ::SparseArrays.SparseMatrixCSC{Float64, Int64}= A
    p              ::Vector{Int64}                               = Vector(undef, n_buses)
    q              ::Vector{Int64}                               = Vector(undef, n_buses)
    Rs             ::Vector{Float64}                             = Vector(undef, n_buses)
    is_factorized  ::Bool                                        = false
    power_surplus  ::Vector{Float64}
    linsolve       ::LinearSolve.LinearCache
    components     ::Vector{BBComponent}
    detached_components::Vector{BBComponent}
    branch_to_comp ::Dict{Int64, BBComponent}    
    vas            ::Vector{Float64}                             = Vector(undef, n_buses)
    dual_rhs       ::Vector{Float64}                             = Vector(undef, n_buses)
    duals          ::Vector{Float64}                             = Vector(undef, n_buses)
    fallback_solver::AbstractSubproblemSolver
    tr_steps       ::Dict{Float64, Float64}                      = Dict{Float64, Float64}()
    branch_data    ::Dict{Int64, Dict{Symbol, Any}}              = Dict()
end


## Fallback Solver - buffers update info so that apply_* is called only when needed
@kwdef mutable struct FallbackSolver <: AbstractSubproblemSolver
    solver              ::AbstractSubproblemSolver
    pm_sub              ::AbstractPowerModel
    contingency_branches::Vector{Int64}                = Vector{Int64}()
    master_var_vals     ::Dict{Symbol, AbstractVector} = Dict()
    is_int              ::Bool                         = false
    nw                  ::Int64
end

##### Data Types for Contingencies #####
Contingency   = Dict{Symbol, Vector{Int64}}
ContingencyID = Int64