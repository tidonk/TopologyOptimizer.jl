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

### Printing ###
"print constraints with rhs >= rhs_lim"
function print_cons_with_large_rhs(pm_sub::AbstractPowerModel; rhs_lim::Float64=5e+6)
    for con in JuMP.all_constraints(pm_sub.model, include_variable_in_set_constraints=false)
        if JuMP.normalized_rhs(con) >= rhs_lim
            @info("Constraint $(JuMP.name(con)) has rhs $(JuMP.normalized_rhs(con))")

            # extract branch ID from name
            branch_id = parse(Int, split(JuMP.name(con), "_")[end])
            @info("Branch $(branch_id): $(pm_sub.data["nw"]["1"]["branch"]["$branch_id"])")
        end
    end
end

### Getters ###
"get the maximum historic change in p0, i.e., a bound on p_delta, from data."
function get_maximum_p0_change(pm, i, pg_0_max_change_fallback=42000)
    @warn("Large pg_0_max_change_fallback :(")
    pg_0_max_change = pg_0_max_change_fallback
    n_timesteps = length(nws(pm))
    for timestep in 2:n_timesteps
        pg_0_current = ref(pm, timestep, :gen)[i]["pg_0"]
        pg_0_previous = ref(pm, timestep-1, :gen)[i]["pg_0"]
        pg_0_max_change = max(max(pg_0_previous, pg_0_current) - min(pg_0_previous, pg_0_current), pg_0_max_change)
    end
    return pg_0_max_change
end

"get a bus id from the fbus and tbus flags. will return the last occurence only."
function get_branch_id_from_fbus_tbus(data, fbus, tbus)
    # create a map from fbus, tbus and tbus, fbus to the branch id
    if !haskey(data, "branch_map")
      data["branch_map"] = Dict{Tuple{Int64, Int64}, Int64}()
      for (br_id, branch) in data["branch"]
        push!(data["branch_map"], (branch["f_bus"], branch["t_bus"]) => branch["index"])
        push!(data["branch_map"], (branch["t_bus"], branch["f_bus"]) => branch["index"])
      end
    end
    # return the branch id based on the tbus, fbus
    return data["branch_map"][(tbus, fbus)]
end

function get_pf_start_from_fbus_tbus(data, fbus::Int64, tbus::Int64)
  # create a map from fbus, tbus and tbus, fbus to the branch id
  if !haskey(data, "branch_map2")
    data["branch_map2"] = Dict{Tuple{Int64, Int64}, Int64}()
    for (br_id, branch) in data["branch"]
      #push!(data["branch_map"], (branch["f_bus"], branch["t_bus"]) => branch["index"])
      push!(data["branch_map2"], (branch["t_bus"], branch["f_bus"]) => branch["index"])
    end
  end
  # return the branch id based on the tbus, fbus
  if haskey(data["branch_map2"], (tbus, fbus))
    return data["branch"][string(data["branch_map2"][(tbus, fbus)])]["pf_start"]
  else
    return data["branch"][string(data["branch_map2"][(fbus, tbus)])]["pt_start"]
  end
end

"returns true iff x-y <= atol, ie if y is at most atol larger than x"
function _isapproxlesseq(x::Float64, y::Float64; atol::Float64=0.0)::Bool
    return x-y <= atol
end

"returns true iff -x+y <= atol, ie if y is at most atol smaller than x"
function _isapproxgreatereq(x::Float64, y::Float64; atol::Float64=0.0)::Bool
    return _isapproxlesseq(y,x,atol=atol)
end

"returns true iff -x+y < atol, ie if y is at most atol smaller than x"
function _isapproxgreater(x::Float64, y::Float64; atol::Float64=0.0)::Bool
    return !_isapproxlesseq(x,y,atol=atol)
end

"returns true iff x-y < atol, ie if y is at most atol smaller than x"
function _isapproxless(x::Float64, y::Float64; atol::Float64=0.0)::Bool
    return !_isapproxgreatereq(x,y,atol=atol)
end

"returns true iff x is more than atol different from y"
function _isapproxdifferent(x::Float64, y::Float64; atol::Float64=0.0)::Bool
    return _isapproxless(x,y,atol=atol) || _isapproxgreater(x,y,atol=atol)
end

"return the default log file path for a case name, settings, and start time in solver_log_dir."
function get_log_file_name(solver_log_dir, name, settings, start_time; prefix="")
    return joinpath(solver_log_dir, "GRB_$(prefix)$(name)_$(settings)_$(replace(replace(start_time, " " => "_"), ":" => "_")).out")
end

function get_log_file_name(data; solver_log_dir = CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:SOLVER_LOG_DIR])
    return get_log_file_name(solver_log_dir, data["name"], data["settings"], data["start_time"], prefix = get(data, "SUBPROBLEM", false) ? "SP_" : "MP_")
end