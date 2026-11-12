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

"data reading and preprocessing function"
function read_and_preprocess_data(casefile::String, 
                                  settings::String,
                                  busbars::Bool, 
                                  application::String, 
                                  slack_types::Vector{Symbol}, 
                                  timesteps::Int64, 
                                  start_time::String
                                )::Tuple{Dict{String, Any}, Bool}

    casename = get_casename(casefile)
    # this is a file containing multiple casefiles in the "files" field of the json
    if occursin(".json", casefile)
        files = JSON.parsefile(casefile)
        casefiles   = files["files"]
        linkingfile = files["ramps"]
        # define main file to use for the first data reading and busbar expansion
        casefile  = casefiles[1]
    else
        casefiles = [casefile]
        linkingfile = ""
    end
    # parse case file
    data = parse_file(casefile)

    # round input data to avoid numerical problems
    round_component_data!(data)

    data["start_time"] = start_time
    data["BUSBARS"] = busbars

    # add slack types
    if ! occursin("SLACK", application) 
        slack_types = Vector{Symbol}()
    end
    add_slack_types_to_data!(data, slack_types)


    # define contingencies, format ::Vector{Dict{Symbol, Vector{Int64}}}
    if occursin("TIME", application)
        data = add_time_steps(data, casefiles, linkingfile, timesteps, slack_types, busbars)
    end

    data["start_time"] = start_time
    data["name"]        = casename
    # get basename of settings
    data["settings"]    = split(basename(settings), ".")[1]
    data["cost_coeffs"] = Dict{String, Float64}()
    for (key, value) in PARAMS[:OBJECTIVE][:COST_COEFFS]
        data["cost_coeffs"][lowercase(string(key))] = value
    end
    data["BUSBARS"]     = busbars
    data["WARMSTART_PROVIDED"] = false
    data["application"] = application


    # adjust rate a
    if occursin(application, "DC_SCOTS_TIME_SLACK") || occursin("DC_SCOTS_TIME_SLACK_BENDERS", application)
        if haskey(data, "nw")
            for n in 1:timesteps
                adjust_ratea!(data["nw"]["$n"])
            end
        end
    end

    if application == "AC_OPF"
        convert_to_generic_cost_data!(data)
    end

    nested = application == "DC_SCOTS_TIME_SLACK_BENDERS_NESTED" || application == "DC_SCOTS_TIME_SLACK_BENDERS_BENDERS"

    # apply certain params directly to data (switchability...)
    apply_params!(data, busbars)

    return data, nested
end

"additional preprocessing function (doesn't require reading)"
function preprocess_data!(data::Dict{String, Any}, 
                          application::String, 
                          slack_types::Vector{Symbol},
                          timesteps::Int64,
                          nested::Bool, 
                          solver_convex::JuMP.MOI.OptimizerWithAttributes,
                          solver_nlp::JuMP.MOI.OptimizerWithAttributes
                         )

    # constant loss estimation
    PARAMS[:LOSS_ESTIMATION][:ESTIMATE_LOSSES] && estimate_losses!(data, solver_nlp)

    if ! PARAMS[:BENDERS][:BASE_PF_IN_MASTER] && occursin("BENDERS", application)
        push!(data["nw"]["1"]["contingencies"], Dict(:branch => Vector{Int64}())::Contingency)
    end
    
    PARAMS[:PRINTING][:INSTANCE_INFO] && print_pm_data_info(data)
    PARAMS[:WARMSTARTING][:ENABLED] && warm_start!(data, solver_convex, slack_types)

    # "wrong" params
    if occursin("BENDERS", application) && PARAMS[:CONTINGENCIES][:SLACK_EQUALS_SEC_SLACK]
        @info("PARAMS[:CONTINGENCIES][:SLACK_EQUALS_SEC_SLACK]=true is not supported in Benders mode and will be ignored.")
    end

    PARAMS[:HEURISTIC_SOLVING][:ENABLED] && enable_heuristic_solving!(data, application, timesteps)
    PARAMS[:PRINTING][:CONFIDENTIAL] && show_system_balance(data) 

    # fix
    remove_duplicate_ref_buses(data)
end

"postprocess results: additional calculations, printing, exporting"
function postprocess_results!(data::Dict{String, Any},
                              result::Dict{String, Any},
                              application::String,
                              timesteps::Int64,
                              nested::Bool,
                              solver_convex::JuMP.MOI.OptimizerWithAttributes,
                              solver_nlp::JuMP.MOI.OptimizerWithAttributes,
                              t_start::Float64
                            )
    casename = data["name"]
    if occursin("SCOTS", application) && PARAMS[:HEURISTIC_SOLVING][:ENABLED] && timesteps == 1
        reset_switchability!(data["nw"]["1"])
    end

    # postprocess, i.e., do all writing that should not be included in the measured runtime
    if result["termination_status"] == JuMP.MOI.OPTIMAL || result["termination_status"] == JuMP.MOI.LOCALLY_SOLVED || (haskey(result, "primal_status") && result["primal_status"] == JuMP.MOI.FEASIBLE_POINT)
        if (! PARAMS[:BENDERS][:BASE_PF_IN_MASTER] || occursin("NESTED", application) || occursin("BENDERS_BENDERS", application)) && ! PARAMS[:POSTPROCESSING][:CHECKS][:CHECK_WITH_COMPACT] && occursin("DC", application)
            # calculate power flow for base case
            calculate_power_flow!(data, result, solver_convex, base_only=true)
        elseif PARAMS[:POSTPROCESSING][:CHECKS][:CHECK_WITH_COMPACT] && occursin("DC", application)
            # calculate power flow for base case
            calculate_power_flow!(data, result, solver_convex, base_only=false)
        end
    end
    
    if result["termination_status"] == JuMP.MOI.OPTIMAL || result["termination_status"] == JuMP.MOI.LOCALLY_SOLVED || (haskey(result, "primal_status") && result["primal_status"] == JuMP.MOI.FEASIBLE_POINT)
        if occursin("OTS", application)
            #print_topology_changes(casename, data, result)
        end
        if occursin("SLACK", application)
            print_slack_values(casename, data, result, data["slack_types"])
            if ! occursin(CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:SLACKED_CASE_FILETYPE], casename) && PARAMS[:WRITING][:SLACKED_CASEFILE]
                write_slack_values(data["slack_types"], result, casename)
                write_matpower_with_slacks(data, result, casename)
            end
        end
        if occursin("TIME", application)
            update_data_safe!(data, result)
            print_generator_changes(casename, data)
            gen_changes = print_pg_delta(casename, data)
            line_switches = print_topology_changes(casename, data)
            result["n_switching_actions"] = sum([length(nw_ls) for (nw, nw_ls) in line_switches])
            pst_shifts = print_transformer_changes(casename, data)

            # print line overloads
            if ! occursin("BENDERS", application) || ! PARAMS[:BENDERS][:BASE_PF_IN_MASTER] || nested # if BASE_PF is not in master, we did calculate_power_flow, so we have at least a base flow here
                print_line_overloads(casename, data, result)
            end
            if occursin("BENDERS", application) && ! nested
                print_line_overloads_benders(casename, data, result)
            end

            # compare extracted cost components (collected in data["objective_costs"]) with the objective value
            #check_obj(data["objective_costs"], result["objective"], base_only=!nested)

            # check line overloads explicitly with compact model
            if PARAMS[:POSTPROCESSING][:OVERLOADS][:EXPLICIT] && ! PARAMS[:POSTPROCESSING][:OVERLOADS][:EXPLICIT_INDIVIDUAL] 
                @info("Checking overloads using compact model.")
            data, result, application = check_relevant_overloads(data, result, solver_convex, for_validation=false)
            elseif PARAMS[:POSTPROCESSING][:OVERLOADS][:EXPLICIT_INDIVIDUAL]
                @info("Checking individual contingencies (debugging mode).")

                # run compact and Benders formulation for this specific contingency
                check_relevant_overloads_individual(data, result, solver_convex, get_relevant_contingencies_from_benders(data, result))
            end

            # check AC-feasibility of DC solution
            PARAMS[:POSTPROCESSING][:AC_FEAS][:CHECK] && check_ac_feasibility!(result, data, solver_nlp)

            result["total_time"] = time() - t_start
            # write solving data
            if PARAMS[:WRITING][:SOLVING_DATA_CSV]
                write_solving_data(casename, data, line_switches, gen_changes)
            end
        end
    else
        result["total_time"] = time() - t_start
    end

    if occursin("BENDERS", application)
        add_solving_data_to_results!(data, result, casename, callback_stats=callback_stats, callback_times=callback_times)
    else
        add_solving_data_to_results!(data, result, casename)
    end

    # result summary
    print_solving_statistics(casename, result, application=application)
end


### Getters ###
"returns parallel branches of given orig_branch, excluding the orig_branch"
function get_parallel_branches(data_branch::Dict{String, Any}, orig_branch::Dict{String, Any})::Vector{Pair{String, Dict}}
    parallel_branches = Vector{Pair{String, Dict}}()
    o_f_bus = orig_branch["f_bus"]::Int64
    o_t_bus = orig_branch["t_bus"]::Int64
    o_index_str = string(orig_branch["index"])

    for (b_id, branch::Dict{String, Any}) in data_branch
        f_bus = branch["f_bus"]::Int64
        t_bus = branch["t_bus"]::Int64
        if (o_f_bus == f_bus && o_t_bus == t_bus) || 
           (o_f_bus == t_bus && o_t_bus == f_bus)
            # do not add the original branch
            if b_id ≠ o_index_str
               push!(parallel_branches, b_id => branch)
            end
        end
        @assert(t_bus ≠ f_bus)
    end
    return parallel_branches
end

### Printing ###
"print step (and message) with time difference to last print"
function printst(step::String, message::String)
    global CLOCK
    @printf("[ %5.2fs | %-21s ] %s\n", time() - CLOCK, step, message)
    CLOCK = time()    
end

"print basic information on the instance at hand"
function _print_pm_data_info(pm_data_hour::Dict{String, Any})
    name        = pm_data_hour["name"]
    hour        = get(pm_data_hour, "hour", "1")
    n_buses     = length(pm_data_hour["bus"])
    n_branches  = length(pm_data_hour["branch"])
    n_gens      = length(pm_data_hour["gen"])
    n_loads     = length(pm_data_hour["load"])
    n_contingencies = length(pm_data_hour["contingencies"])
    n_monitored = length(get(pm_data_hour, "monitoring", pm_data_hour["branch"]))
    n_switchable = count(branch -> branch["switchable"], values(pm_data_hour["branch"]))
    n_off = count(branch -> branch["status"] == 0.0, values(pm_data_hour["branch"]))
    n_substations = length(get(pm_data_hour, "substations", []))
    n_couplers = length(filter(branch -> branch["type"] == "K", collect(values(pm_data_hour["branch"]))))
    n_breakers = count(branch -> branch["type"] == "breaker", collect(values(pm_data_hour["branch"])))
    A = ([length(substation["busbar"]) for substation in values(get(pm_data_hour, "substations", []))])
    n_substation_sizes = [(i, count(==(i), A)) for i in unique(A)]  
    # sort n_substation_sizes by first entry of each tuple
    n_substation_sizes = sort(n_substation_sizes, by=x->x[1])
    # create string in the form 1: %3d / 2: %3d / 3: %3d / 4: %3d ... for each unique substation size
    n_substation_sizes_str = join(["$i: $(n)" for (i, n) in n_substation_sizes], " / ")

    network_data_str = @sprintf("Buses: %4d, Branches: %4d, Gens: %4d, Loads: %4d", n_buses, n_branches, n_gens, n_loads)
    substation_data_str = @sprintf("Substations: %d (%s), Couplers: %d, Breakers: %d", n_substations, n_substation_sizes_str, n_couplers, n_breakers)
    contingency_monitoring_str = @sprintf("Contingencies: %d, Monitored: %d", n_contingencies, n_monitored)
    switchable_str = @sprintf("Switchable: %d, Off: %d", n_switchable, n_off)
    @info("Data Info:\n Name: $name, Hour: $hour\n $network_data_str\n $substation_data_str\n $contingency_monitoring_str\n $switchable_str")
end
function print_pm_data_info(pm_data_hour::Dict{String, Any})
    if haskey(pm_data_hour, "nw")
        for (_, nw_data) in pm_data_hour["nw"]
            _print_pm_data_info(nw_data)
        end
    else
        _print_pm_data_info(pm_data_hour)
    end
end

### Data Consistency ###
"returns `gen`s minimum production when considering re-dispatch: `max((gen['pg_0'] - gen['RDP-']), gen['pmin'])`"
function get_gen_min_prod(gen::Dict)::Float64
    return max((gen["pg_0"] - gen["RDP-"]), gen["pmin"])
end

"returns `gen`s maximum production when considering re-dispatch: `min((gen['pg_0'] + gen['RDP+']), gen['pmax'])`"
function get_gen_max_prod(gen::Dict)::Float64
    return min((gen["pg_0"] + gen["RDP+"]), gen["pmax"])
end

function fix_gsk_switchability!(data::Dict)
    n_fixed = 0
    ## filter out all gens that are turned off and not switchable
    if haskey(data, "gsk_factors")
        ## filter gsks that can not be turned on & propagate to gsks
        # filter gsks that can not be turned on
        for (_, gsk_factors) in data["gsk_factors"]
            for (i, gsk_factor) in gsk_factors

                gen = data["gen"]["$i"]

                # check if factor can be reached
                min_prod = (1 - gsk_factor) * get_gen_min_prod(gen)
                max_prod = (1 - gsk_factor) * get_gen_max_prod(gen)
                min_total_prod = gsk_factor * sum(get_gen_min_prod(data["gen"]["$j"]) for (j, factor) in gsk_factors if j != i)
                max_total_prod = gsk_factor * sum(get_gen_max_prod(data["gen"]["$j"]) for (j, factor) in gsk_factors if j != i)
                if min_prod > max_total_prod || max_prod < min_total_prod # infeasible, this gen has to remain off
                    gen["status"] = 0
                    gen["switchable"] = false
                    @debug("Inconsistency: Gen $(gen["name"]) with gsk factor $gsk_factor cannot be turned on, setting status = 0 and switchable = false.")
                    n_fixed += 1

                    break # we'll propagate this
                end
            end

        end

        # propagate to gsks
        remove_gsks = Vector{String}()
        for (_, gen) ∈ data["gen"]

            # do we want to remove this?
            if !is_gen_switchable(gen) && gen["status"] == 0

                # get the gsk
                if haskey(gen, "gsk_name")
                    gsk_name = gen["gsk_name"]

                    # did we already check this?
                    if gsk_name ∉ remove_gsks
                        push!(remove_gsks, gsk_name)
                        gsk_factors = data["gsk_factors"][gsk_name]

                        # set all gens to non-switchable and off
                        for (other_gen_id, _) in gsk_factors
                            other_gen = data["gen"]["$other_gen_id"]

                            other_gen["switchable"] = false
                            other_gen["status"]     = 0
                        end

                    end

                end
            end

        end

        # remove all gsks that are switched off
        filter!(name_factors -> name_factors[1] ∉ remove_gsks, data["gsk_factors"])
    end
    return n_fixed
    printst("GSKs", "Added GSK names to gens and fixed switchability.")
end

"adds a 'gsk_name' field to each gen that is part of one"
function add_gsk_name_to_gens!(data::Dict)
    if haskey(data, "gsk_factors")

        for (gsk_name, gsk_factors) in data["gsk_factors"]
            for (gen_id, _) in gsk_factors
                @assert !haskey(data["gen"]["$gen_id"], "gsk_name")

                data["gen"]["$gen_id"]["gsk_name"] = gsk_name
            end
        end

    end
end

"check gen consistency"
function check_gen_consistency!(data)
    ics_logger = _LOGGER
    n_errors = 0
    for (nw, nw_data) in data["nw"]
        add_gsk_name_to_gens!(nw_data)
        
        for (g_id, gen) in nw_data["gen"]
            name = get(gen, "name", g_id)

            # negative generation is not permitted
            if (gen["pg_0"] - gen["RDP-"]) < 0
                @warn("Inconsistency: Gen $name. pg_0 - RDP- = $(gen["pg_0"]) - $(gen["RDP-"]) = $(gen["pg_0"] - gen["RDP-"]) < 0.")
            end
            if gen["pmin"] < 0
                @warn("Inconsistency: Gen $name. pmin = $(gen["pmin"]) < 0. Setting pmin = 0.")
                n_errors += 1
            end
            # check switchable flag
            if isapprox(gen["pg_0"] - gen["RDP-"], 0, atol=PARAMS[:PRECISION][:POWER_PREC]) && gen["pmin"] > 0 && ! gen["switchable"]
                @warn("Inconsistency: Gen $name. pg_0 - RDP- = $(gen["pg_0"]) - $(gen["RDP-"]) = 0, but gen not switchable. Setting switchable = true.")
                n_errors += 1
                gen["switchable"] = true
            end
            if ! isapprox(gen["pg_0"] - gen["RDP-"], 0, atol=PARAMS[:PRECISION][:POWER_PREC]) && gen["switchable"]
                @warn("Inconsistency: Gen $name. pg_0 - RDP- = $(gen["pg_0"]) - $(gen["RDP-"]) = $(gen["pg_0"] - gen["RDP-"]) > 0, but gen is switchable. Setting switchable = false.")
                n_errors += 1
                gen["switchable"] = false
            end

            # check status flag
            # if pmin = 0, the gen is always on and not switchable
            if gen["pmin"] == 0 && gen["switchable"] == true
                gen["switchable"] = false
                gen["status"] = 1
            end
            # correct pmin from data
            if ! isapprox(gen["pg_0"] - gen["RDP-"], 0, atol=PARAMS[:PRECISION][:POWER_PREC]) && (gen["pg_0"] - gen["RDP-"]) < gen["pmin"]
                if gen["pg_0"] - gen["RDP-"] < gen["pmin"] - PARAMS[:PRECISION][:POWER_PREC] 
                    @warn("Inconsistency: Gen $name. pg_0 - RDP- = $(gen["pg_0"]) - $(gen["RDP-"]) = $(gen["pg_0"] - gen["RDP-"]) < $(gen["pmin"]) = pmin. Reducing pmin = $(gen["pmin"]) to $(gen["pg_0"] - gen["RDP-"]).")
                    n_errors += 1
                end
                gen["pmin"] = gen["pg_0"] - gen["RDP-"]
            end
            # correct pmin if gen is off, but should be able to turn on
            if gen["status"] == 0 && gen["switchable"] == true && (gen["pg_0"] + gen["RDP+"]) < gen["pmin"] && (gen["pg_0"] + gen["RDP+"]) > 0
                if gen["pg_0"] + gen["RDP+"] < gen["pmin"] - PARAMS[:PRECISION][:POWER_PREC] 
                    @warn("Inconsistency: Gen $name. pg_0 + RDP+ = $(gen["pg_0"]) + $(gen["RDP+"]) = $(gen["pg_0"] + gen["RDP+"]) < $(gen["pmin"]) = pmin. Reducing pmin = $(gen["pmin"]) to $(gen["pg_0"] + gen["RDP+"]).")
                    n_errors += 1
                end
                gen["pmin"] = gen["pg_0"] + gen["RDP+"]
            end
            # correct pmax from data
            if (gen["pg_0"] + gen["RDP+"]) > gen["pmax"] 
                if gen["pg_0"] + gen["RDP+"] > gen["pmax"] + PARAMS[:PRECISION][:POWER_PREC] 
                    @warn("Inconsistency: Gen $name. pg_0 + RDP+ = $(gen["pg_0"]) + $(gen["RDP+"]) = $(gen["pg_0"] + gen["RDP+"]) > $(gen["pmax"]) = pmax. Increasing pmax = $(gen["pmax"]) to $(gen["pg_0"] + gen["RDP+"])")
                    n_errors += 1
                end
                gen["pmax"] = gen["pg_0"] + gen["RDP+"]
            end
            # correct pg_0
            if gen["pg_0"] > 0 && gen["pg_0"] < gen["pmin"]
                @warn("Inconsistency: Gen $name. pg_0 = $(gen["pg_0"]) < $(gen["pmin"]) = pmin. Setting pmin to pg_0.")
                n_errors += 1
                gen["pmin"] = gen["pg_0"]
            end
            if gen["pg_0"] == 0 && gen["switchable"] && gen["RDP+"] > 0 && gen["pmin"] > gen["ramp_60_pos"]
                @warn("Inconsistency: Gen $name. pg_0 = 0, switchable = true, RDP+ > 0, but pmin = $(gen["pmin"]) >  $(gen["ramp_60_pos"]) = ramp_60_pos. Setting pmin = $(gen["pmin"]) to $(gen["ramp_60_pos"]).")
                n_errors += 1
                gen["pmin"] = gen["ramp_60_pos"]
            end
            if ! isapprox(gen["pg"], gen["pg_0"], atol=PARAMS[:PRECISION][:POWER_PREC])
                @warn("Inconsistency: Gen $name. pg = $(gen["pg"]) ≠ pg_0 = $(gen["pg_0"]). Setting pg = pg_0.")
                n_errors += 1
                gen["pg"] = gen["pg_0"]
            end
            if gen["pg_0"] < 0
                @warn("Inconsistency: Gen $name. pg_0 = $(gen["pg_0"]) < 0. Setting pg_0 = 0.")
                n_errors += 1
                old_pg_0 = deepcopy(gen["pg_0"])
                bus_id = gen["gen_bus"]
                gen["pg_0"] = 0.0
                gen["pg"] = 0.0
                gen["pmax"] = gen["pmax"] - old_pg_0
                gen["pmin"] = gen["pmin"] - old_pg_0
                gen["RDP+"] = gen["RDP+"] - old_pg_0
                gen["RDP-"] = gen["RDP-"] - old_pg_0
                # add load to data with negative previous pg_0
                loads_max_idx = maximum([load["index"] for (_, load) in nw_data["load"]])
                idx = loads_max_idx + 1
                push!(nw_data["load"], string(idx) => create_bus_load(- old_pg_0, 0.0, bus_id, 1, idx))
                @info("Added slack load for gen $name with negative pg_0 = $old_pg_0. New load index: $idx.")
            end
        end
        printst("Gen Consistency", "Resolved $n_errors generator inconsistencies.")

        n_gsk_errors = fix_gsk_switchability!(nw_data)
        printst("GSK Consistency", "Resolved $n_gsk_errors GSK inconsistencies.")
    end

    # if on, apply gen override, and remove all gens that are not switchable and are off
    for (nw, nw_data) in data["nw"]
        n_fixed = 0
        if PARAMS[:LIMITS][:GENERATOR][:GEN_REDISPATCH_OVERRIDE]
            for (g_id, gen) in nw_data["gen"]
            # all gens are on, and no limit on pmin
                gen["status"] = 1
                gen["switchable"] = false
                gen["pmin"] = min(gen["pmin"], 0)

                n_fixed += 1
                continue
            end
            printst("Gen Override", "Set all gens to on due to activated gen redispatch override.")
        end

        # remove gens that are not switchable and are off
        n_before = length(nw_data["gen"])
        filter!(gid_gen -> is_gen_switchable(gid_gen[2]) || gid_gen[2]["status"] == 1, nw_data["gen"])
        printst("Gen Cleanup", "Removed $(n_before - length(nw_data["gen"])) non-switchable and off gens.")
    end
end

"check whether all non-electrical branches are classified correctly for LinEqSysBBSolver"
function check_substation_integrity!(data_branch::Dict{String, Any}, data_substations::Dict{String, Dict{String, Vector{Any}}})
    #if PARAMS[:BENDERS][:SUBPROBLEM_SOLVER_TYPE] == LinEqSysBBSolver
        # flatten all coupler branches in data_substations[...]["K"]
        coupler_branches = Vector{Int64}()
        for (_, substation) in data_substations
            for branch_id in substation["K"]
                push!(coupler_branches, parse(Int64, branch_id))
            end
        end
        @assert(length(coupler_branches) == length(unique(coupler_branches)), "Coupler branches are not unique: $(coupler_branches)")
        breaker_branches = Vector{Int64}()
        for (_, branch) in data_branch
            if branch["type"] == "breaker"
                push!(breaker_branches, branch["index"])
            end
        end
        @assert(length(breaker_branches) == length(unique(breaker_branches)), "Breaker branches are not unique: $(breaker_branches)")

        for (b_id, branch) ∈ data_branch
            if ! branch["electrical"]
                if branch["transformer"]
                    @error("Non-electrical branch $(branch["branch_name_intern"]) ($(branch["index"])) is a transformer.")
                    branch["transformer"] = false
                end
                if branch["type"] ∉ ["K", "breaker"]
                    @error("Non-electrical branch $(branch["branch_name_intern"]) ($(branch["index"])) is not a coupler or breaker.")
                else
                    if branch["type"] == "K" && branch["index"] ∉ coupler_branches
                        @error("Non-electrical branch $(branch["branch_name_intern"]) ($(branch["index"])) is not in any substation.")
                        f_bus = branch["f_bus"]
                        t_bus = branch["t_bus"]

                        found = false
                        for (_, substation) in data_substations

                            if "$f_bus" ∈ substation["busbar"] || "$t_bus" ∈ substation["busbar"]
                                push!(substation["K"], "$b_id")

                                if "$f_bus" ∉ substation["busbar"]
                                    push!(substation["busbar"], "$f_bus")
                                elseif "$t_bus" ∉ substation["busbar"]
                                    push!(substation["busbar"], "$t_bus")
                                end

                                found = true
                                break
                            end
                        end

                        if !found
                            data_substations["substation_$b_id"] = Dict("busbar" => ["$f_bus", "$t_bus"], "aux_bus" => Vector{Dict{String, Any}}(), "K" => ["$b_id"])
                        end
                    elseif branch["type"] == "breaker" && branch["index"] ∉ breaker_branches
                        @error("Non-electrical branch $(branch["branch_name_intern"]) ($(branch["index"])) is not a breaker")
                    end
                end
            end
        end
    #end
    printst("Substation Integrity", "Checked substation data integrity.")
end

### Slack ###
"add slack load to represent present system inbalance"
function add_single_slack_load!(data)
    # add slack node based on system inbalance
    total_pd = 0.0
    total_pg = 0.0
    for (load_id, load) in data["load"]
        if load["pd"] != 0 && load["status"] == 0
            @warn("Load $(load["name"]) is off but has a positive pd of $(load["pd"]).")
        end
        total_pd += load["pd"]
    end
    for (gen_id, gen) in data["gen"]
        if gen["pg"] > 0 && gen["status"] == 0
            @warn("Gen $(gen["name"]) is off but has a positive pg of $(gen["pg"]).")
        end
        total_pg += gen["pg"]
    end
    @debug("[System Balance] Total demand: $total_pd, total generation: $total_pg, imbalance: $(total_pd - total_pg).")
    data["slack"] = total_pd - total_pg
    bus = collect(values(filter(b -> b[2]["bus_type"] == 3, data["bus"])))[begin]

    # add load for the slack at the slack bus
    push!(data["load"], "-1" => create_bus_load(- data["slack"], 0, bus["index"], 1, -1))
    @info("[Slack] Added a slack load of $(data["slack"]) at bus $(bus["name"]).")
end

"for each bus, add a slack load. The relative distribution depends on the absolute bus load."
function add_slack_distributed!(data::Dict{String, Any}, data_bus::Dict{String, Any}, data_gen::Dict{String, Any}, data_load::Dict{String, Any})
    # calculate total slack
    total_pd = 0.0
    total_pg = 0.0
    for (_, load) in data_load
        total_pd += load["pd"]
    end
    for (_, gen) in data_gen
        if gen["pg"] > 0 && gen["status"] == 0
            @warn("Gen $(gen["name"]) is off but has a positive pg of $(gen["pg"]).")
        end
        total_pg += gen["pg"]
    end
    @debug("[System Balance] Total demand: $total_pd, total generation: $total_pg, imbalance: $(total_pd - total_pg).")
    data["slack"] = total_pd - total_pg

    # to each bus, add a slack load that is proportional to the absolute load at the bus (should be only one load, since executed after merge_loads!)
    total_abs_load = sum([abs(load["pd"]) for (_, load) in data_load])
    bus_loads = Dict{Int64, Float64}()
    # get all buses that have any load connected. at this point in the code, they must be connected to the main network, so we do not power yet unconnected parts here.
    for (_, bus::Dict{String, Any}) ∈ data_bus
        load_sum = 0.0
        bus_index = bus["index"]::Int64

        for (_, load::Dict{String, Any}) ∈ data_load
            if load["load_bus"]::Int64 == bus_index && load["pd"]::Float64 ≠ 0.0
                load_sum += abs(load["pd"])
            end
        end

        if load_sum > 0.0
            bus_loads[bus_index] = load_sum
        end
    end

    n_slack_loads = 0
    loads_max_idx = maximum([load["index"] for (_, load) in data_load])
    for (bus_id, bus_load) ∈ bus_loads
        # add load for the slack at the slack bus
        individual_distributed_slack = data["slack"] * (bus_load / total_abs_load)
        idx = loads_max_idx + 1 + n_slack_loads
        push!(data_load, string(idx) => create_bus_load(- individual_distributed_slack, 0.0, bus_id, 1, idx))
        n_slack_loads += 1
    end
    printst("Slack Load", "Added $n_slack_loads slack loads to all $(length(keys(bus_loads))) powered buses.")
end

### Merging ###
"merges buses if the branch between them is not switchable and has a very low resistance"
function merge_buses!(data_branch::Dict{String, Any}, data_bus::Dict{String, Any}, data_gen::Dict{String, Any}, data_load::Dict{String, Any}, data_substations::Dict{String, Dict{String, Vector{Any}}})
    # We can merge buses if _all_ of the following conditions are met for _all_ branches connecting the two buses:
    #   - branch is not switchable/shiftable
    #   - branch is not in contingencies
    #   - branch is not electrical (i.e., has a low resistance)
    # We merge by deleting the branch (A,B) (and all parallel branches), and moving all loads, gens and branches connected to bus B to bus A, and delete bus B.

    n_branches_before = length(data_branch)
    n_buses_before    = length(data_bus)
    n_consistency_errors = 0
    
    # checked branches (to keep track during parallel branches), buses to merge into first one, branches to delete
    checked = Vector{Int64}()
    ActionItem = Vector{Union{Int64, Vector{String}}}
    action_list = Vector{ActionItem}()

    for (b_id, branch) in data_branch
        # if the branch is switchable/shiftable, electrical, or in contingencies, or off (and not switchable), skip
        if is_branch_switchable(branch) || is_branch_shiftable(branch) || is_branch_electrical(branch) || branch["in_contingencies"] || branch["br_status"] == 0
            continue
        end
        @assert(branch["type"] ≠ "breaker", "Why was this breaker branch created if it is not switchable?")

        # if we already looked at this branch in the parallel branches of another branch, skip
        if branch["index"] ∈ checked
            continue
        else
            push!(checked, branch["index"])

            # check if there are parallel branches, delete them as well, or skip this merging if there are any that do not meet the conditions
            parallel_branches::Vector{Pair{String, Dict}} = get_parallel_branches(data_branch, branch)
            action_list_item::ActionItem = [branch["f_bus"], branch["t_bus"], [b_id]] # ActionItem = pairs of buses to merge together, and all buses that are between them

            fail = false # flag that becomes true if a parallel branch has different properties, and thus no merging may take place
            for (p_b_id, p_branch) in parallel_branches
                push!(checked, p_branch["index"])

                # check if the parallel branch meets the conditions, add to action item to later remove all parallel branches if no fails occured
                if is_branch_switchable(p_branch)
                    @debug("[Merging Elements] Parallel branch $(p_branch["branch_name_intern"]) ($(p_branch["index"])) is switchable.")
                    n_consistency_errors += 1
                    fail = true
                    push!(action_list_item[3], p_b_id)
                elseif is_branch_shiftable(p_branch)
                    @debug("[Merging Elements] Parallel branch $(p_branch["branch_name_intern"]) ($(p_branch["index"])) is shiftable.")
                    n_consistency_errors += 1
                    fail = true
                    push!(action_list_item[3], p_b_id)
                elseif is_branch_electrical(p_branch)
                    @debug("[Merging Elements] Parallel branch $(p_branch["branch_name_intern"]) ($(p_branch["index"])) is electrical.")
                    #@printf("  - [Parallel] Type = %s, R = %8.4f, br_r = %8.4f\n", p_branch["type"], p_branch["R"], p_branch["br_r"])
                    #printf("  - [Original] Type = %s, R = %8.4f, br_r = %8.4f\n",   branch["type"],   branch["R"],   branch["br_r"])
                    n_consistency_errors += 1
                    #fail = true
                    push!(action_list_item[3], p_b_id)
                elseif branch["in_contingencies"]
                    @debug("[Merging Elements] Parallel branch $(p_branch["branch_name_intern"]) ($(p_branch["index"])) is in contingencies.")
                    n_consistency_errors += 1
                    fail = true
                    push!(action_list_item[3], p_b_id)
                elseif branch["br_status"] == 0
                    @debug("[Merging Elements] Parallel branch $(p_branch["branch_name_intern"]) ($(p_branch["index"])) is off. Removing it as well.")
                    push!(action_list_item[3], p_b_id)
                else
                    push!(action_list_item[3], p_b_id)
                end
            end
        end

        # at least one parallel branch did not meet the requirements, so skip the merging here
        if ! fail
            push!(action_list, action_list_item)
        end
    end
    sort!(action_list, by=x->x[1])
    printst("Bus Merging", "Found $(length(action_list)) buses to merge (including circularity duplicates).")

    while ! isempty(action_list)
        (merge_bus_idx::Int64, delete_bus_idx::Int64, delete_branches_idxs::Vector{String}) = pop!(action_list)

        if merge_bus_idx == delete_bus_idx
            # delete all branches that were supposed to be deleted
            for b_id in delete_branches_idxs
                delete!(data_branch, b_id)
            end
            continue
        end

        merge_bus = data_bus[string(merge_bus_idx)]
        delete_bus = data_bus[string(delete_bus_idx)]

        str = "[Merging Elements] Merging buses $(merge_bus["index"]) ($(get(merge_bus, "bus_name_copy", ""))) and $(delete_bus["index"]) ($(get(delete_bus, "bus_name_copy", ""))) connected by: \n"
        for b_id in delete_branches_idxs
            branch = data_branch[b_id]
            str *= @sprintf("      - %-28s (%4d: %4d -> %4d, br_r = %6.2f)\n", get(branch, "branch_name_intern", ""), branch["index"], branch["f_bus"], branch["t_bus"], branch["br_r"])
            delete!(data_branch, b_id)
        end
        #@info(str)

        # check if reference bus
        if data_bus[string(delete_bus_idx)]["bus_type"] == 3
            @debug("Re-assigning reference bus from $(delete_bus["index"]) to $(merge_bus["index"])")
            data_bus[string(merge_bus_idx)]["bus_type"] = 3
        end 

        # delete bus, replace all references to this bus in the action_list by the merge bus
        delete!(data_bus, string(delete_bus_idx))
        for (j, (m_b, d_b, _)) in enumerate(action_list)
            if m_b == delete_bus_idx
                #println("Replace: $(m_b) <- $(d_b) by $(merge_bus_idx) <- $(d_b)")
                action_list[j][1] = merge_bus_idx
            end
            if d_b == delete_bus_idx
                #println("Replace: $(m_b) <- $(d_b) by $(m_b) <- $(merge_bus_idx)")
                action_list[j][2] = merge_bus_idx
            end
        end

        # delete bus from pm_data_hour["substations"]
        for (_, substation) in data_substations
            filter!(x -> x ≠ string(delete_bus_idx), substation["busbar"])
            # aux_buses do not have to be changed since we don't merge breakers
        end

        # move loads to the new bus
        for (_, load::Dict{String, Any}) in data_load
            if load["load_bus"]::Int64 == delete_bus_idx
                load["load_bus"] = merge_bus_idx
                load["source_id"][2] = merge_bus_idx # ["bus", load["load_bus"]]
            end
        end

        # move gens to the new bus
        for (_, gen::Dict{String, Any}) in data_gen
            if gen["gen_bus"]::Int64 == delete_bus_idx
                gen["gen_bus"] = merge_bus_idx
                gen["source_id"][2] = merge_bus_idx # ["bus", gen["gen_bus"]]
            end
        end

        # connect all branches that were previously connected to the old bus to the new bus
        for (_, branch::Dict{String, Any}) in data_branch
            if branch["f_bus"]::Int64 == delete_bus_idx || branch["t_bus"]::Int64 == delete_bus_idx
                if branch["f_bus"] == delete_bus_idx
                    branch["f_bus"] = merge_bus_idx
                end
                if branch["t_bus"] == delete_bus_idx
                    branch["t_bus"] = merge_bus_idx
                end
            end
        end
    end
    filter!(bid_branch -> bid_branch[2]["f_bus"] != bid_branch[2]["t_bus"], data_branch)

    # fix pm_data_hour["contingencies"]
    #n_contingencies_before = length(pm_data_hour["contingencies"])
    #filter!(x -> x[:branch][1] ∈ parse.(Int64, keys(data_branch)), pm_data_hour["contingencies"])

    printst("Bus Merging", "Merged buses with low resistance branches. Remaining branches: $(length(data_branch)) ($(n_branches_before - length(data_branch)) removed), remaining buses: $(length(data_bus)) ($(n_buses_before - length(data_bus)) removed), found $n_consistency_errors consistency errors.")

end

"for each bus, merge all its connected loads into one single load with index equal to the bus index"
function merge_loads!(data_bus::Dict{String, Any}, data_load::Dict{String, Any})
    n_loads     = length(data_load)
    new_loads   = Dict{String, Dict{String, Any}}()
    # TODO: build connected_loads as a dictionary
    connected_loads = Dict{Int64, Vector{Dict{String, Any}}}()
    for (_, load) in data_load
        push!(get!(connected_loads, load["load_bus"], Vector{Dict{String, Any}}()), load)
    end

    for (_, bus) in data_bus
        _connected_loads = get(connected_loads, bus["index"], [])
        if !isempty(_connected_loads)
            # merge loads
            new_load = deepcopy(_connected_loads[1])
            new_load["index"] = bus["index"]
            new_load["source_id"] = ["bus", bus["index"]]
            new_load["pd"] = sum(load["pd"] for load in _connected_loads)
            new_load["qd"] = sum(load["qd"] for load in _connected_loads)

            # add new load
            push!(new_loads, string(new_load["index"]) => new_load)

            # delete old loads
            for load in _connected_loads
                delete!(data_load, string(load["index"]))
            end
        end
    end

    # add new loads
    for (l_id, load) in new_loads
        if load["pd"] ≠ 0 #|| load["qd"] ≠ 0
            data_load[l_id] = load
        end
    end
    printst("Load Merging", "Aggregated $n_loads loads into $(length(data_load)) loads.")
end