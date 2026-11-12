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

"helper function to add slack types to dataframe"
function add_slack_types_to_data!(data, slack_types::Vector{Symbol})
    data["slack_types"] = slack_types
    return data
end

"wrapper to initialize slack types into dataframe"
function solve_secure_model(file::String, model_type::Type, optimizer, build_method, slack_types::Vector{Symbol}; kwargs...)
    data = parse_file(file)
    return solve_secure_model(data, model_type, optimizer, build_method, slack_types; kwargs...)
end

"wrapper to initialize slack types into dataframe"
function solve_secure_model(data::Dict{String,<:Any}, model_type::Type, optimizer, build_method, slack_types::Vector{Symbol}; kwargs...)
    add_slack_types_to_data!(data, slack_types)
    return solve_secure_model(data, model_type, optimizer, build_method; kwargs...)
end

"wrapper to initialize contingencies into dataframe"
function solve_secure_model(file::String, model_type::Type, optimizer, build_method; kwargs...)
    data = parse_file(file)
    return solve_secure_model(data, model_type, optimizer, build_method; kwargs...)
end

"wrapper to initialize contingencies into dataframe"
function solve_secure_model(data::Dict{String,<:Any}, model_type::Type, optimizer, build_method; kwargs...)
    return solve_model(data, model_type, optimizer, build_method; ref_extensions=[ref_add_on_off_va_bounds!], kwargs...)
end

"rename variables so that they can be read in as LP by Gurobi"
function rename_variables!(model::JuMP.Model)
    # Get all variables in the model
    variables = JuMP.all_variables(model)

    for var in variables
        old_name = JuMP.name(var)  # Original variable name
        new_name = filter(a -> !isspace(a), old_name)  # Remove spaces
        new_name = replace(new_name, "[" => "_", "]" => "_")  # Replace brackets with underscores
        new_name = replace(new_name, "η" => "eta")  # Replace eta with eta
        #new_name = replace(new_name, )
        new_name = new_name[3:end]  # Remove the first two characters (which are the network id)
        if old_name != new_name
            JuMP.MOI.set(model, JuMP.MOI.VariableName(), var, new_name)  # Set new variable name
        end
    end
    return model
end

"round the coefficients of all variables in each constraint and the objective function to the given precision"
function round_coeffs!(model; digits=PARAMS[:PRECISION][:SOLVER_PREC_DIGITS])
    # round coeffs and rhss of constraints
    for con in JuMP.all_constraints(model, include_variable_in_set_constraints=false)     
        object = JuMP.constraint_object(con)
        # skip if the constraint is not a linear or quadratic expression
        if typeof(object.func) <: JuMP.GenericAffExpr
            # round rhs of the constraint
            JuMP.set_normalized_rhs(con, round(JuMP.normalized_rhs(con), digits=digits))
            # round coeffs of the constraint
            variables = [k for k in keys(object.func.terms)]
            for var in variables
                JuMP.set_normalized_coefficient(con, var, round(JuMP.normalized_coefficient(con, var), digits=digits))
            end
        end
    end
    object = JuMP.objective_function(model)
    # round coeffs of obj func
    if typeof(object) <: JuMP.GenericAffExpr
        variables = [k for k in keys(object.terms)]
        JuMP.@objective(model, JuMP.objective_sense(model), sum(var * round(JuMP.coefficient(object, var), digits=digits) for var in variables))
    end
    return model
end

"given a Gurobi direct_model, write presolved file"
function write_presolved_to_file(model::JuMP.Model, presolvedfile::String)
    # write presolved model
    presolvedP = Ref{Ptr{Cvoid}}(C_NULL)
    GRBpresolvemodel(JuMP.unsafe_backend(model), presolvedP)
    GRBwrite(presolvedP[], presolvedfile)
    JuMP.write_to_file(JuMP.read_from_file(presolvedfile), presolvedfile)
end

# remove duplicate ref buses
function remove_duplicate_ref_buses(data)
    
    if haskey(data, "nw")
        for (nw_id, nw_data) in data["nw"]
            found_ref_bus = false

            for (b, bus) in nw_data["bus"]

                if bus["bus_type"] == 3
                    if !found_ref_bus
                        found_ref_bus = true
                    else
                        bus["bus_type"] = 1
                    end
                end

            end

        end
    else
        found_ref_bus = false
        
        for (b, bus) in data["bus"]

            if bus["bus_type"] == 3
                if !found_ref_bus
                    found_ref_bus = true
                else
                    bus["bus_type"] = 1
                end
            end

        end
    end

end

# Modified from PowerModels.jl (https://github.com/lanl-ansi/PowerModels.jl), see LICENSE
""
function solve_model(data::Dict{String,<:Any}, model_type::Type, optimizer, build_method;
        ref_extensions=[], solution_processors=[], relax_integrality=false,
        multinetwork=false, kwargs...)

    if multinetwork != _IM.ismultinetwork(data)
        model_requirement = multinetwork ? "multi-network" : "single-network"
        data_type = _IM.ismultinetwork(data) ? "multi-network" : "single-network"
        Memento.error(_LOGGER, "attempted to build a $(model_requirement) model with $(data_type) data")
    end

    start_time = time()
    pm = instantiate_model(data, model_type, build_method; ref_extensions=ref_extensions, kwargs...)
    Memento.debug(_LOGGER, "pm model build time: $(time() - start_time)")

    if PARAMS[:WRITING][:LP] || PARAMS[:WRITING][:MPS] || PARAMS[:WRITING][:PRESOLVED]
        rename_variables!(pm.model)
    end

    if PARAMS[:PRECISION][:SOLVER_PREC_DIGITS] >= 0 && PARAMS[:ROUNDING][:MODEL_COEFFS]
        pm.model = round_coeffs!(pm.model)
    end
    
    if PARAMS[:WRITING][:LP]
        if occursin("BENDERS", pm.data["application"])
            JuMP.MOI.set(pm.model, Gurobi.CallbackFunction(), nothing)
        end
        JuMP.write_to_file(pm.model, "$(pm.data["name"])_$(pm.data["application"]).lp")
        @info("LP file written to $(pm.data["name"])_$(pm.data["application"]).lp")
        if occursin("BENDERS", pm.data["application"])
            throw("Benders decomposition cannot be continued after writing LP file because the callback function is removed")
        end
    end
    
    if PARAMS[:WRITING][:MPS]
        JuMP.write_to_file(pm.model, "$(pm.data["name"])_$(pm.data["application"]).mps")
        @info("MPS file written to $(pm.data["name"])_$(pm.data["application"]).mps")
    end

    if PARAMS[:WRITING][:REW] 
        JuMP.write_to_file(pm.model, "$(pm.data["name"])_$(pm.data["application"]).rew")
        @info("MPS file written to $(pm.data["name"])_$(pm.data["application"]).rew")
    end

    if PARAMS[:WRITING][:DEC] 
        @assert(PARAMS[:WRITING][:LP] )
        # read cur_model.lp and create cur_model.dec from it
        write_decomposition("timesteps")
    end
    
    start_time = time()
    #JuMP.set_start_values(pm.model, variable_primal_start=(x -> JuMP.value(x))) # not required, since solver env stays the same
    result = optimize_model!(pm, relax_integrality=relax_integrality, optimizer=optimizer, solution_processors=solution_processors)
    Memento.debug(_LOGGER, "pm model solve and solution time: $(time() - start_time)")

    if PARAMS[:WRITING][:PRESOLVED] 
        write_presolved_to_file(pm.model, "$(pm.data["name"])_$(pm.data["application"]).presolved.lp")
        @info("Presolved file written to $(pm.data["name"])_$(pm.data["application"]).presolved.lp")
    end

    if PARAMS[:DEBUG][:RUN_IIS]
        #JuMP.write_to_file(pm.model, "model.lp")
        if result["termination_status"] == JuMP.MOI.INFEASIBLE || result["termination_status"] == JuMP.MOI.INFEASIBLE_OR_UNBOUNDED
            JuMP.compute_conflict!(pm.model)
            iis_model, _ = JuMP.copy_conflict(pm.model)
            print(iis_model)
            if PARAMS[:WRITING][:IIS]
                if occursin("BENDERS", pm.data["application"])
                    JuMP.MOI.set(iis_model, Gurobi.CallbackFunction(), nothing)
                end
                JuMP.write_to_file(iis_model, "conflict_model.lp")
            end
        end
    end

    return result
end