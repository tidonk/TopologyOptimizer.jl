"""Process-local state for an optional helper-to-main incumbent mailbox.

The paired launcher owns a private directory and identical input/settings for
both workers. Only serializable names/values cross processes, never JuMP refs
or Gurobi environments. Publication is an atomic rename on the same filesystem.
"""
Base.@kwdef mutable struct IncumbentExchange
    directory::String
    role::Symbol
    variables::Vector{JuMP.VariableRef}
    names::Vector{String}
    contingency_keys::Vector{String}
    sequence::Int = 0
    published::Float64 = Inf
    last_poll::Float64 = -Inf
    acknowledged::Int = 0
    last_submission::Float64 = -Inf
    submitted_objective::Float64 = Inf
    last_log::Float64 = -Inf
    last_primal::Float64 = Inf
    last_dual::Float64 = -Inf
    packet_mtime::Float64 = -Inf
    packet::Any = nothing
    attempts::Int = 0
end

struct ParallelHelperFinished <: Exception
    diagnostic::Dict{String,Any}
end

function parallel_main_finished()
    return get(ENV, "INNOPTEM_EXCHANGE_STOP_WITH_MAIN", "0") == "1" &&
        get(ENV, "INNOPTEM_EXCHANGE_ROLE", "") == "donor" &&
        isfile(joinpath(get(ENV, "INNOPTEM_EXCHANGE_DIR", ""), "main.done"))
end

function stop_parallel_helper!(pm; solved=false)
    diagnostic = Dict{String,Any}("reason" => "main_finished",
        "tal" => PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP],
        "prior_warmstart_solver_seconds" => get(pm.data, "warmstart_time", 0.0),
        "rscots_stages" => get(pm.data, "rscots_stages", Any[]),
        "restricted_solver_controls" => get(pm.data, "restricted_solver_controls", nothing))
    if solved
        diagnostic["last_solve_seconds"] = solve_time(pm.model)
        diagnostic["last_status"] = string(termination_status(pm.model))
        diagnostic["last_objective"] = has_values(pm.model) ? objective_value(pm.model) : nothing
    end
    throw(ParallelHelperFinished(diagnostic))
end

function init_incumbent_exchange(pm; contingencies=Contingency[])
    ccall(:jl_generating_output, Cint, ()) == 1 && return nothing
    directory = get(ENV, "INNOPTEM_EXCHANGE_DIR", "")
    role = Symbol(get(ENV, "INNOPTEM_EXCHANGE_ROLE", "off"))
    tal = PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP]
    isempty(directory) && return nothing
    role in (:donor, :receiver, :control) || error("Invalid incumbent exchange role")
    helper = get(ENV, "INNOPTEM_EXCHANGE_HELPER", "tal1")
    main_tal = parse(Int, get(ENV, "INNOPTEM_EXCHANGE_MAIN_TAL", "2"))
    active = role == :donor ? (helper == "scopf" ? tal == 0 : helper == "rscots" ? tal in (0, 1) : helper == "rscots-sequential" ? 0 <= tal <= main_tal : tal == 1) : tal == main_tal
    active || return nothing
    parallel_main_finished() && stop_parallel_helper!(pm)
    if role == :donor && helper == "rscots" && tal == 1
        restrict_parallel_helper!(pm)
    elseif role == :donor && helper == "rscots-sequential" && tal == main_tal && !haskey(pm.data, "_rscots_stage")
        haskey(pm.data, "_rscots_final_stage") || error("Staged RSCOTS requires a feasible SCOPF initialization")
        stage = pm.data["_rscots_final_stage"]
        restrict_rscots_stage!(pm, stage)
        set_time_limit_sec(pm.model, min(time_limit_sec(pm.model), stage["time_limit"]))
    end
    if role == :donor && helper in ("rscots", "rscots-sequential") && tal > 0
        config = PARAMS[:WARMSTARTING][:RSCOTS]
        config[:HELPER_MIP_GAP] < 0 || set_optimizer_attribute(pm.model, "MIPGap", config[:HELPER_MIP_GAP])
        config[:HELPER_DISABLE_OBJECTIVE_STOP] && set_optimizer_attribute(pm.model, "BestObjStop", -1e100)
        pm.data["restricted_solver_controls"] = Dict(
            "mipgap" => get_optimizer_attribute(pm.model, "MIPGap"),
            "bestobjstop" => get_optimizer_attribute(pm.model, "BestObjStop"))
    end
    variables = JuMP.all_variables(pm.model)
    names = JuMP.name.(variables)
    all(!isempty, names) && allunique(names) || error("Exchange requires unique variable names")
    contingency_keys = [join([string(k) * ":" * join(sort(v), ";")
                        for (k, v) in sort(collect(c); by=p -> string(first(p)))], "|") for c in contingencies]
    state = IncumbentExchange(; directory, role, variables, names, contingency_keys)
    # RSCOTS continues the same publication sequence after its SCOPF stage.
    packet_path = joinpath(directory, "incumbent.bin")
    if role == :donor && isfile(packet_path)
        packet = deserialize(packet_path)
        state.sequence = packet.sequence
        state.published = packet.objective
    end
    return state
end

"Restrict the helper to its top-k branches while preserving the shared variable schema."
function restrict_parallel_helper!(pm)
    haskey(pm.data, "result") || error("RSCOTS requires a feasible SCOPF initialization")
    score = Symbol(get(ENV, "INNOPTEM_EXCHANGE_SCORE", "SENS"))
    methods = (:RND, :FULL, :SENS, :CONGPROX, :ETA, :FIXDUAL, :ADMITTANCE_SECANT)
    score in methods || error("Unknown RSCOTS score: $score")
    count = parse(Int, get(ENV, "INNOPTEM_EXCHANGE_SWITCHABLE", "16"))
    count > 0 || error("RSCOTS branch count must be positive")
    nw_data = pm.data["nw"]["1"]
    score_start = time()
    if score == :ADMITTANCE_SECANT
        config = PARAMS[:WARMSTARTING][:RSCOTS]
        scores = calc_br_sw_heur_scores(nw_data, config[:WEIGHTS])
        JuMP.set_time_limit_sec(pm.model, min(JuMP.time_limit_sec(pm.model), config[:TIME_LIMIT]))
    elseif score == :FIXDUAL
        remaining = parse(Float64, ENV["INNOPTEM_EXCHANGE_LIMIT"]) - (time() - exchange_epoch(ENV["INNOPTEM_EXCHANGE_DIR"]))
        solver = JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(CONSTS[:OTHER][:GUROBI_ENV]),
            "Threads" => PARAMS[:SOLVER][:MULTITHREADING][:N_THREADS_MIP_MP],
            "TimeLimit" => max(0.0, remaining),
            "LogFile" => joinpath(CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:RESULTS_DIR], "fixing-duals.log"))
        # This auxiliary LP computes scores, not transferable incumbents/bounds.
        scores, diagnostic = withenv("INNOPTEM_EXCHANGE_DIR" => nothing) do
            parallel_fixed_topology_duals(pm.data, solver;
                perturbation_check=get(ENV, "INNOPTEM_CHECK_DUAL_PERTURBATIONS", "0") == "1")
        end
        diagnostic["scoring_wall_seconds"] = time() - score_start
        open(joinpath(CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:RESULTS_DIR], "fixing-duals.json"), "w") do io
            JSON.print(io, diagnostic, 2)
        end
        isnothing(scores) && error("FIXDUAL requires an optimal fixing LP; see fixing-duals.json")
    else
        funcs = Dict(:RND => calc_br_sw_heur_scores_rnd, :FULL => calc_br_sw_heur_scores_full,
            :SENS => calc_br_sw_heur_scores_sens, :CONGPROX => calc_br_sw_heur_scores_congprox,
            :ETA => calc_br_sw_heur_scores_eta)
        # Single-method rankings need no normalization (including all-zero scores).
        scores = funcs[score](nw_data; nw_sol=pm.data["result"]["solution"]["nw"]["1"],
            PARAMS[:HEURISTIC_SOLVING][:RSCOTS][:HEUR]...)
    end
    z = var(pm, 1, :z_branch)
    ids = [id for (id, branch) in ref(pm, 1, :branch) if is_branch_switchable(branch)]
    all(id -> isfinite(scores[id]), ids) || error("Nonfinite RSCOTS scores")
    if score == :ADMITTANCE_SECANT
        ids = admittance_branch_order(nw_data, scores, ids;
            initially_open_first=get(PARAMS[:WARMSTARTING][:RSCOTS][:HEUR], :INITIALLY_OPEN_FIRST, false))
    else
        sort!(ids; by=id -> (-scores[id], id))
    end
    selected = Set(first(ids, min(count, length(ids))))
    for id in ids
        id in selected && continue
        branch = ref(pm, 1, :branch, id)
        JuMP.fix(z[id], get(branch, "status", branch["br_status"]); force=true)
    end
    # Keep all variables and contingency identities so the receiver can check
    # the original objective and every imported value without a schema conversion.
    path = joinpath(CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:RESULTS_DIR], "helper-selection.json")
    open(path, "w") do io
        JSON.print(io, Dict("score" => string(score), "requested" => count,
            "selected" => sort(collect(selected)), "eligible" => length(ids),
            "scores" => Dict(string(id) => scores[id] for id in ids),
            "solver_time_limit" => JuMP.time_limit_sec(pm.model),
            "initially_open_first" => score == :ADMITTANCE_SECANT && get(PARAMS[:WARMSTARTING][:RSCOTS][:HEUR], :INITIALLY_OPEN_FIRST, false),
            "tie_break" => "ascending branch id", "scoring_wall_seconds" => time() - score_start), 2)
    end
end

function exchange_epoch(directory)
    return parse(Float64, read(joinpath(directory, "epoch"), String))
end

function exchange_event(state, event; objective=NaN, dual=NaN, sequence=state.sequence, detail="")
    t = time() - exchange_epoch(state.directory)
    path = joinpath(state.directory, "$(state.role).events.csv")
    open(path, "a") do io
        println(io, join((t, event, objective, dual, sequence, detail), ','))
    end
end

"Publish only after the caller has checked every contingency without errors/cuts."
function publish_incumbent!(state, pm, cb_data)
    isnothing(state) && return
    state.role == :donor || return
    values = JuMP.callback_value.(cb_data, state.variables)
    objective = JuMP.value(v -> JuMP.callback_value(cb_data, v), JuMP.objective_function(pm.model))
    publish_incumbent_values!(state, values, objective)
end

"Publish a fully checked recourse repair before the donor solve terminates."
function publish_repaired_incumbent!(state, pm, cb_data, repair)
    isnothing(state) && return
    state.role == :donor && !isnothing(repair) || return
    variables, repaired_values = repair
    by_variable = Dict(zip(variables, repaired_values))
    values = [get(by_variable, v, JuMP.callback_value(cb_data, v)) for v in state.variables]
    by_variable = Dict(zip(state.variables, values))
    objective = JuMP.value(v -> by_variable[v], JuMP.objective_function(pm.model))
    publish_incumbent_values!(state, values, objective)
end

"Publish a validated callback solution or the solver's final accepted incumbent."
function publish_incumbent_values!(state, values, objective)
    state.role == :donor || return
    length(values) == length(state.names) || error("Shared incumbent value count differs")
    all(isfinite, values) && isfinite(objective) || return
    objective < state.published - 1e-8 * max(1.0, abs(objective)) || return
    state.sequence += 1
    packet = (version=1, sequence=state.sequence, objective=objective,
              names=state.names, values=values, contingency_keys=state.contingency_keys)
    temp = joinpath(state.directory, "incumbent.tmp")
    serialize(temp, packet)
    mv(temp, joinpath(state.directory, "incumbent.bin"); force=true)
    state.published = objective
    if get(ENV, "INNOPTEM_EXCHANGE_AUDIT", "0") == "1"
        open(joinpath(state.directory, "incumbent.sol"), "w") do io
            for (name, value) in zip(state.names, values)
                println(io, name, " ", value)
            end
        end
    end
    exchange_event(state, "published"; objective)
end

"Poll independently of the fractional-separation budget; never wait for the donor."
function exchange_callback!(state, pm, cb_data, cb_where)
    isnothing(state) && return
    if state.role == :donor && time() - state.last_poll >= 0.1
        state.last_poll = time()
        if parallel_main_finished()
            GRBterminate(cb_data.model)
            return :stop
        end
    end
    if cb_where == GRB_CB_MIPSOL && get(ENV, "INNOPTEM_EXCHANGE_AUDIT", "0") == "1"
        obj = Ref{Cdouble}()
        GRBcbget(cb_data, cb_where, GRB_CB_MIPSOL_OBJ, obj)
        exchange_event(state, "candidate"; objective=obj[])
    end
    codes = if cb_where == GRB_CB_MIP
        (GRB_CB_MIP_OBJBST, GRB_CB_MIP_OBJBND)
    elseif cb_where == GRB_CB_MIPNODE
        (GRB_CB_MIPNODE_OBJBST, GRB_CB_MIPNODE_OBJBND)
    elseif cb_where == GRB_CB_MIPSOL
        (GRB_CB_MIPSOL_OBJBST, GRB_CB_MIPSOL_OBJBND)
    else
        return
    end
    now = time()
    primal, dual = Ref{Cdouble}(), Ref{Cdouble}()
    GRBcbget(cb_data, cb_where, codes[1], primal)
    GRBcbget(cb_data, cb_where, codes[2], dual)
    # Bound logging must not consume the receiver's independent polling budget.
    if primal[] != state.last_primal || dual[] != state.last_dual || now - state.last_log >= 1.0
        exchange_event(state, "bounds"; objective=primal[], dual=dual[])
        state.last_primal, state.last_dual, state.last_log = primal[], dual[], now
    end
    state.role == :receiver || return
    # Queue external solutions outside the current candidate's separation.
    # Immediate heuristic return statuses do not establish a persistent update.
    cb_where == GRB_CB_MIP || return
    now - state.last_poll >= 0.05 || return
    state.last_poll = now
    path = joinpath(state.directory, "incumbent.bin")
    isfile(path) || return
    mtime = stat(path).mtime
    if isnothing(state.packet) || mtime != state.packet_mtime
        packet = deserialize(path)
        packet.version == 1 || error("Unsupported incumbent packet")
        length(packet.names) == length(packet.values) == length(state.names) &&
            allunique(packet.names) && Set(packet.names) == Set(state.names) ||
            error("Donor and receiver master variable schemas differ")
        packet.contingency_keys == state.contingency_keys || error("Donor and receiver contingency identities/order differ")
        isfinite(packet.objective) && all(isfinite, packet.values) || error("Nonfinite shared incumbent")
        state.packet, state.packet_mtime = packet, mtime
    end
    packet = state.packet
    packet.sequence > state.acknowledged || return
    is_new = packet.sequence > state.sequence
    state.sequence = packet.sequence
    is_new && (state.attempts = 0)
    target = is_new ? packet.objective : state.submitted_objective
    target < primal[] - 1e-8 * max(1.0, abs(target)) || begin
        exchange_event(state, is_new ? "not_improving" : "confirmed"; objective=primal[])
        state.acknowledged = packet.sequence
        return
    end
    # Back off persistent rejections instead of injecting the same point four
    # times a second for hours. A new donor incumbent is always tried immediately.
    !is_new && now - state.last_submission < min(30.0, 0.25 * 2.0^min(state.attempts - 1, 7)) && return
    by_name = Dict(zip(packet.names, packet.values))
    values = [by_name[name] for name in state.names]
    all(isfinite, values) || error("Nonfinite shared incumbent")
    receiver_objective = JuMP.value(v -> by_name[JuMP.name(v)], JuMP.objective_function(pm.model))
    isapprox(receiver_objective, packet.objective; rtol=1e-8, atol=1e-6) ||
        error("Shared incumbent objective differs in receiver model")
    state.last_submission = now
    state.attempts += 1
    variables = state.variables
    mode = "full"
    if !is_new && pm isa AbstractPowerModel
        # Independently computed dual cuts can differ at the feasibility
        # tolerance. A conservative increase of recourse estimates lets the
        # receiver validate the same decisions without repeatedly rejecting a
        # numerically marginal point. Never weaken or bypass a lazy cut.
        eta = Set(get_master_vars(pm; nw=1)[:η])
        for (i, v) in enumerate(variables)
            if v in eta
                margin = max(100 * SOLVER_PRECISION, 100 * eps(max(1.0, abs(values[i]))))
                values[i] = min(values[i] + margin, JuMP.has_upper_bound(v) ? JuMP.upper_bound(v) : Inf)
            end
        end
        mode = "recourse_margin"
    end
    submitted_values = Dict(zip(variables, values))
    state.submitted_objective = JuMP.value(v -> submitted_values[v], JuMP.objective_function(pm.model))
    status = JuMP.MOI.submit(pm.model, JuMP.MOI.HeuristicSolution(cb_data), variables, values)
    exchange_event(state, "submitted"; objective=state.submitted_objective, detail=string(status) * "|" * mode)
end

"Audit whether newly generated receiver cuts exclude the donor's point."
function audit_exchange_cuts!(state, cuts)
    isnothing(state) && return
    state.role == :receiver && get(ENV, "INNOPTEM_EXCHANGE_AUDIT", "0") == "1" || return
    path = joinpath(state.directory, "incumbent.bin")
    isfile(path) || return
    packet = deserialize(path)
    by_name = Dict(zip(packet.names, packet.values))
    for con in cuts
        violation = JuMP.value(v -> by_name[JuMP.name(v)], con.func) - con.set.upper
        if violation > SOLVER_PRECISION
            eta = join([JuMP.name(v) for v in keys(con.func.terms) if startswith(JuMP.name(v), "benders_eta")], "|")
            exchange_event(state, "cut_rejects_donor"; objective=packet.objective,
                           dual=violation, sequence=packet.sequence, detail=eta)
        end
    end
end

"Synchronize the first solve (SCOPF) of both processes on a common wall clock."
function parallel_solve_start!(pm)
    ccall(:jl_generating_output, Cint, ()) == 1 && return
    directory = get(ENV, "INNOPTEM_EXCHANGE_DIR", "")
    isempty(directory) && return
    role = get(ENV, "INNOPTEM_EXCHANGE_ROLE", "off")
    parallel_main_finished() && stop_parallel_helper!(pm)
    if PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP] > 0 && get(ENV, "INNOPTEM_EXCHANGE_AUDIT", "0") == "1"
        JuMP.MOI.Utilities.attach_optimizer(JuMP.backend(pm.model))
        grb = JuMP.unsafe_backend(pm.model)
        GRBupdatemodel(grb)
        GRBwrite(grb, joinpath(directory, "$(role).lp"))
    end
    ready = joinpath(directory, "$(role).ready")
    if !isfile(ready)
        write(ready, "ready\n")
        deadline = time() + parse(Float64, get(ENV, "INNOPTEM_STARTUP_TIMEOUT", "3600"))
        while !isfile(joinpath(directory, "epoch"))
            time() < deadline || error("Timed out waiting for paired solve")
            sleep(0.05)
        end
    end
    if haskey(ENV, "INNOPTEM_EXCHANGE_LIMIT")
        remaining = parse(Float64, ENV["INNOPTEM_EXCHANGE_LIMIT"]) - (time() - exchange_epoch(directory))
        current_limit = JuMP.time_limit_sec(pm.model)
        JuMP.set_time_limit_sec(pm.model, min(something(current_limit, Inf), max(0.0, remaining)))
    end
end

"Record solver-reported final bounds, including the feasible SCOPF anchor."
function parallel_solve_finish!(pm)
    ccall(:jl_generating_output, Cint, ()) == 1 && return
    directory = get(ENV, "INNOPTEM_EXCHANGE_DIR", "")
    isempty(directory) && return
    role = Symbol(get(ENV, "INNOPTEM_EXCHANGE_ROLE", "off"))
    state = pop!(pm.data, "_parallel_incumbent_exchange", nothing)
    if isnothing(state)
        state = IncumbentExchange(; directory, role, variables=JuMP.VariableRef[], names=String[], contingency_keys=String[])
    elseif JuMP.has_values(pm.model)
        # Gurobi can accept a repaired heuristic incumbent without another raw
        # candidate reaching our no-cut publication path. Never lose that final
        # (solver-accepted) solution when the helper exits.
        publish_incumbent_values!(state, JuMP.value.(state.variables), JuMP.objective_value(pm.model))
    end
    tal = PARAMS[:LIMITS][:BRANCH][:TAL_TIMESTEP]
    objective = JuMP.has_values(pm.model) ? JuMP.objective_value(pm.model) : Inf
    if role == :donor && haskey(pm.data, "_rscots_final_stage") && !haskey(pm.data, "_rscots_stage")
        stage = pop!(pm.data, "_rscots_final_stage")
        merge!(stage, Dict("objective" => objective, "solve_seconds" => solve_time(pm.model),
            "status" => string(termination_status(pm.model))))
        push!(pm.data["rscots_stages"], stage)
    end
    dual = tal == 0 ? -Inf : JuMP.objective_bound(pm.model)
    exchange_event(state, tal == 0 ? "scopf" : "final"; objective, dual,
                   detail=string(JuMP.termination_status(pm.model)))
    if get(ENV, "INNOPTEM_EXCHANGE_STOP_WITH_MAIN", "0") == "1"
        main_tal = parse(Int, get(ENV, "INNOPTEM_EXCHANGE_MAIN_TAL", "2"))
        if role in (:receiver, :control) && tal == main_tal
            write(joinpath(directory, "main.done"), string(termination_status(pm.model)))
        elseif parallel_main_finished()
            stop_parallel_helper!(pm; solved=true)
        end
    end
end
