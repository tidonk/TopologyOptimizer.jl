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

###########################################################################
################################ Buses ####################################
###########################################################################
function create_bus_load(pd::Float64, qd::Float64, load_bus::Int, status::Int, index::Int)
    return Dict{String,Any}(
        "pd" => pd,
        "qd" => qd,
        "load_bus" => load_bus,
        "status" => status,
        "index" => index,
        "source_id" => ["bus", load_bus]
    )
end

###########################################################################
############################## Branches ###################################
###########################################################################
"return true iff branch is electrical"
function is_branch_electrical(branch::Dict{String, Any})::Bool
    # return BUSBARS || get(branch, "electrical", true)
    return branch["electrical"]::Bool
end
function is_branch_electrical(branch::Dict{Symbol, Any})::Bool
    return branch[:electrical]::Bool
end
function is_branch_electrical(pm, nw::Int64, branch_id::Int64)::Bool
    return is_branch_electrical(ref(pm, nw, :branch, branch_id))
end

"return true iff branch is monitored"
function is_branch_monitored(branch::Dict{String, Any})::Bool
    return get(branch, "monitored", true) && get(branch, "electrical", true)
end
function is_branch_monitored(branch::Dict{Symbol, Any})::Bool
    return get(branch, :monitored, true) && get(branch, :electrical, true)
end
"return true iff branch is monitored"
function is_branch_monitored(pm, nw::Int64, branch_id::Int64)
    return is_branch_monitored(ref(pm, nw, :branch, branch_id))
end

"return switchability of branches"
function is_branch_switchable(branch::Dict{String, Any})::Bool
    return branch["switchable"]::Bool
end
function is_branch_switchable(branch::Dict{Symbol, Any})::Bool
    return branch[:switchable]::Bool
end
function is_branch_switchable(pm, nw::Int64, branch_id::Int64)::Bool
    return is_branch_switchable(ref(pm, nw, :branch, branch_id))
end

"returns the `rate_a` of branch with any adjustments applied, currently only the ub mults"
function get_adjusted_rate_a(branch::Dict{String, Any}, for_n1::Bool; apply_ub_mults=(PARAMS[:LIMIT_MULT]::Dict{Symbol, Any})[:APPLY_UB_MULTS]::Bool)::Float64
    rate_a = branch["rate_a"]::Float64
    
    if apply_ub_mults
        if for_n1
            max_loading = get(branch, "max_loading_n1", 100.0)::Float64
        else
            max_loading = get(branch, "max_loading", 100.0)::Float64
        end
        
        return rate_a * max_loading/100.0
    else
        return rate_a
    end 
end

"returns the `rate_a` of branch with any adjustments applied, currently only the ub mults"
function get_adjusted_rate_a(branch::Dict{Symbol, Any}, for_n1::Bool; apply_ub_mults=(PARAMS[:LIMIT_MULT]::Dict{Symbol, Any})[:APPLY_UB_MULTS]::Bool)::Float64
    rate_a = branch[:rate_a]::Float64
    
    if apply_ub_mults
        if for_n1
            max_loading = get(branch, :max_loading_n1, 100.0)::Float64
        else
            max_loading = get(branch, :max_loading, 100.0)::Float64
        end

        return rate_a * max_loading/100.0
    else
        return rate_a
    end    
end

"returns the `rate_a` of branch with any adjustments applied, currently only the ub mults"
function get_adjusted_rate_a(data::Dict, branch::Dict; apply_ub_mults=(PARAMS[:LIMIT_MULT]::Dict{Symbol, Any})[:APPLY_UB_MULTS]::Bool, for_n1::Bool=false)::Float64
    return get_adjusted_rate_a(branch, (for_n1 || get(data, "SUBPROBLEM", false)::Bool), apply_ub_mults=apply_ub_mults)
end

"returns the `rate_a` of branch with any adjustments applied, currently only the ub mults"
function get_adjusted_rate_a(pm::AbstractPowerModel, branch::Dict; apply_ub_mults=(PARAMS[:LIMIT_MULT]::Dict{Symbol, Any})[:APPLY_UB_MULTS]::Bool, for_n1::Bool=false)::Float64
    return get_adjusted_rate_a(pm.data, branch, apply_ub_mults=apply_ub_mults, for_n1=for_n1)
end

"get the status of `branch`"
function get_status(branch::Dict)::Int
    return get(branch, "status", branch["br_status"])
end

#############################################################################
############################## Transformers #################################
#############################################################################
"return shiftability of transformers"
function is_branch_shiftable(branch::Dict{String, Any})::Bool
    return branch["shiftable"]::Bool
end
function is_branch_shiftable(branch::Dict{Symbol, Any})::Bool
    return branch[:shiftable]::Bool
end
function is_branch_shiftable(pm, nw::Int64, branch_id::Int64)::Bool
    return is_branch_shiftable(ref(pm, nw, :branch, branch_id))
end

"returns the current shift of a possible transformer at branch, where `tr_shift_var_obj` is either the :tr_shift var or its value, and `branch` is the branch dict"
function get_transformer_shift(pm, n, i)::Union{Float64, JuMP.AffExpr}
    branch  = ref(pm, n, :branch, i)
    return  is_branch_shiftable(branch) ? var(pm, n, :tr_shift, i) * branch["shift_per_step"] : branch["shift"]
end
function get_transformer_shift(solver::AbstractSubproblemSolver, branch::Dict)::Union{Float64, JuMP.AffExpr}
    return is_branch_shiftable(branch) ? solver.tr_steps[branch["index"]] * branch["shift_per_step"] : branch["shift"]
end
function get_transformer_shift(solver::AbstractSubproblemSolver, branch::Dict{Symbol, Any})::Union{Float64, JuMP.AffExpr}
    return is_branch_shiftable(branch) ? solver.tr_steps[branch[:index]] * branch[:shift_per_step] : branch[:shift]
end

"returns the shift of a transformer in the current solution or, if not shiftable, from the data"
function get_transformer_shift_master(branch::Dict{String, Any}, master_var_vals::Dict{Symbol, AbstractVector})::Union{Float64, JuMP.AffExpr}
    return is_branch_shiftable(branch) ? master_var_vals[:tr_shift][branch["index"]] * branch["shift_per_step"] : get(branch, "shift", 0.0)
end

#############################################################################
############################### Generators ##################################
#############################################################################
"return switchability of gens"
function is_gen_switchable(gen::Dict{String, Any})::Bool
    return PARAMS[:LIMITS][:GENERATOR][:GEN_SWITCHING] && get(gen, "switchable", true)
end
function is_gen_switchable(gen::Dict{Symbol, Any})::Bool
    return PARAMS[:LIMITS][:GENERATOR][:GEN_SWITCHING] && get(gen, :switchable, true)
end
function is_gen_switchable(pm, nw::Int64, gen_id::Int64)::Bool
    return is_gen_switchable(ref(pm, nw, :gen, gen_id))
end