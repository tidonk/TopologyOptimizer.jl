# This file contains modified code derived from PowerModels.jl.
# Copyright (c) 2016, Los Alamos National Security, LLC. All rights reserved.
# These portions retain the license and disclaimer in LICENSES/PowerModels.txt.
# PowerModels.jl code in this file has been modified for TopologyOptimizer.jl.

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

#### General Assumptions of these OTS Models ####
#
# - if the branch status is 0 in the input, it is out of service and forced to 0 in OTS
# - the network will be maintained as one connected component (i.e. at least n-1 edges)
#

""
function solve_ots(file, model_type::Type, optimizer; kwargs...)
    return solve_model(file, model_type, optimizer, build_ots; ref_extensions=[ref_add_on_off_va_bounds!], kwargs...)
end

""
function build_ots(pm::AbstractPowerModel)
    variable_branch_indicator(pm, relax=false)
    variable_bus_voltage_on_off(pm)
    variable_gen_power(pm)
    variable_branch_power(pm)
    variable_dcline_power(pm)

    objective_min_fuel_and_flow_cost(pm)

    constraint_model_voltage_on_off(pm)

    for i in ids(pm, :ref_buses)
        constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_ohms_yt_from_on_off(pm, i)
        constraint_ohms_yt_to_on_off(pm, i)

        constraint_voltage_angle_difference_on_off(pm, i)

        constraint_thermal_limit_from_on_off(pm, i)
        constraint_thermal_limit_to_on_off(pm, i)
    end

    for i in ids(pm, :dcline)
        constraint_dcline_power_losses(pm, i)
    end

    return pm
end


""
function ref_add_on_off_va_bounds!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    apply_pm!(_ref_add_on_off_va_bounds!, ref, data)
end

"sum up all active loads pd's"
function calc_p_bound(data::Dict{String,<:Any})
    p_max = 0.0

    for i in keys(data["load"])
        if data["load"][i]["status"] == 1
            p_max += data["load"][i]["pd"]
        end
    end

    return p_max
end

""
function _ref_add_on_off_va_bounds!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    if ! haskey(ref, :off_angmin) || ! haskey(ref, :off_angmax)
        off_angmin, off_angmax = calc_theta_delta_bounds(data)
        ref[:off_angmin], ref[:off_angmax] = off_angmin, off_angmax
    end

    if ! haskey(ref, :p_max)
        p_max = calc_p_bound(data)
        ref[:p_max] = p_max
    end
end
