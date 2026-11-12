#!/usr/bin/env julia
# Compare runs of TopologyOptimizer across a set of settings (against the default baseline)
# on all instances under check/instances/, with TAL=1.
#
# Usage:
#   julia --project=. check/scripts/compare_runs.jl [settings_folder]
#
# settings_folder defaults to settings/5-partial_decomp.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using TopologyOptimizer
using Statistics
using Printf

const REPO_ROOT      = normpath(joinpath(@__DIR__, "..", ".."))
const INSTANCE_ROOT  = joinpath(REPO_ROOT, "check", "instances")
const DEFAULT_SETTS  = joinpath(REPO_ROOT, "settings", "default.toml")
const TAL            = 1

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

"Collect all instance .m files under check/instances/."
function collect_instances()
    files = String[]
    for d in sort(readdir(INSTANCE_ROOT))
        sub = joinpath(INSTANCE_ROOT, d)
        isdir(sub) || continue
        for f in sort(readdir(sub))
            endswith(f, ".m") && push!(files, joinpath(sub, f))
        end
    end
    return files
end

"Wilcoxon signed-rank test (two-sided) on paired differences x .- y.
Returns p-value via normal approximation with tie + zero correction."
function wilcoxon_signed_rank_p(x::AbstractVector, y::AbstractVector)
    @assert length(x) == length(y)
    d = x .- y
    d = d[d .!= 0]                       # drop zeros
    n = length(d)
    n < 2 && return NaN
    a = abs.(d)
    p = sortperm(a)
    ranks = zeros(Float64, n)
    i = 1
    while i <= n
        j = i
        while j < n && a[p[j+1]] == a[p[i]]
            j += 1
        end
        avg = (i + j) / 2
        for k in i:j
            ranks[p[k]] = avg
        end
        i = j + 1
    end
    W_plus  = sum(ranks[d .> 0])
    W_minus = sum(ranks[d .< 0])
    W = min(W_plus, W_minus)
    μ = n * (n + 1) / 4
    σ = sqrt(n * (n + 1) * (2n + 1) / 24)
    σ == 0 && return NaN
    z = (W - μ) / σ
    # two-sided p-value via standard normal CDF (Abramowitz & Stegun 7.1.26)
    return 2 * (1 - 0.5 * (1 + erf_approx(abs(z) / sqrt(2))))
end

"Abramowitz–Stegun erf approximation (max abs error ≈ 1.5e-7)."
function erf_approx(x::Float64)
    sign_x = sign(x)
    x = abs(x)
    a1, a2, a3, a4, a5, p = 0.254829592, -0.284496736, 1.421413741,
                            -1.453152027, 1.061405429, 0.3275911
    t = 1.0 / (1.0 + p * x)
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x)
    return sign_x * y
end

stars(p) = isnan(p) ? "" : (p < 0.001 ? "***" : p < 0.01 ? "**" : p < 0.05 ? "*" : "")

geomean(v) = isempty(v) ? NaN : exp(mean(log.(max.(v, eps()))))

# ---------------------------------------------------------------------------
# Solve
# ---------------------------------------------------------------------------

"Run TopologyOptimizer on `casefile` with the given settings file. Returns (objective, time_total)."
function run_one(casefile::String, settings_file::String)
    add_par = Dict{Symbol,Any}(:LIMITS => Dict(:BRANCH => Dict(:TAL_TIMESTEP => TAL)))
    use_bb  = endswith(casefile, "_bb.m")
    @info "Running" instance=basename(casefile) settings=basename(settings_file) busbars=use_bb
    res = TopologyOptimizer.main(
        casefile;
        application = "DC_SCOTS_TIME_SLACK_BENDERS",
        busbars     = use_bb,
        params_file = settings_file,
        additional_params = add_par,
    )
    obj = get(res, "objective", NaN)
    t   = haskey(res, "time") && isa(res["time"], Dict) ? get(res["time"], "total", NaN) : NaN
    return obj, t
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main(settings_folder::String = joinpath(REPO_ROOT, "settings", "5-partial_decomp"))
    @assert isdir(settings_folder) "Settings folder not found: $settings_folder"

    instances = collect_instances()
    @assert !isempty(instances) "No instances found under $INSTANCE_ROOT"

    # baseline = default, then each .toml in folder
    settings_files = vcat(DEFAULT_SETTS,
                          sort([joinpath(settings_folder, f)
                                for f in readdir(settings_folder) if endswith(f, ".toml")]))
    setting_names  = vcat("baseline",
                          [splitext(basename(f))[1] for f in settings_files[2:end]])

    n_inst = length(instances)
    n_set  = length(settings_files)
    objs   = fill(NaN, n_inst, n_set)
    times  = fill(NaN, n_inst, n_set)

    for (j, sf) in enumerate(settings_files), (i, inst) in enumerate(instances)
        try
            objs[i, j], times[i, j] = run_one(inst, sf)
        catch e
            @warn "Run failed" instance=basename(inst) settings=basename(sf) exception=(e, catch_backtrace())
        end
    end

    inst_names = [basename(splitdir(f)[1]) * "/" * basename(f) for f in instances]

    print_table("Objective", inst_names, setting_names, objs; baseline_col = 1)
    print_table("Time [s]",  inst_names, setting_names, times; baseline_col = 1)

    # Objective equivalence check
    base_obj = objs[:, 1]
    println("\n=== Objective equivalence vs baseline (rel. diff > 1e-3) ===")
    any_diff = false
    for j in 2:n_set
        for i in 1:n_inst
            (isnan(objs[i, j]) || isnan(base_obj[i])) && continue
            rel = abs(objs[i, j] - base_obj[i]) / max(abs(base_obj[i]), 1.0)
            if rel > 1e-3
                @warn "Objective mismatch" instance=inst_names[i] setting=setting_names[j] baseline=base_obj[i] got=objs[i, j] rel
                any_diff = true
            end
        end
    end
    any_diff || println("  All objectives match baseline within tolerance.")
end

function print_table(title::String, rows::Vector{String}, cols::Vector{String},
                     M::Matrix{Float64}; baseline_col::Int = 1)
    println("\n=== ", title, " ===")
    name_w = max(maximum(length, rows), length("instance"))
    col_w  = max(12, maximum(length, cols))
    fmt(x) = isnan(x) ? "—" : @sprintf("%.3g", x)
    # header
    print(rpad("instance", name_w), " | ")
    println(join([lpad(c, col_w) for c in cols], " | "))
    println(repeat("-", name_w + 3 + length(cols) * (col_w + 3)))
    # body
    for i in eachindex(rows)
        print(rpad(rows[i], name_w), " | ")
        println(join([lpad(fmt(M[i, j]), col_w) for j in 1:length(cols)], " | "))
    end
    println(repeat("-", name_w + 3 + length(cols) * (col_w + 3)))
    # footer rows
    for (label, fn) in (("sum", v -> sum(skipmissing(filter(!isnan, v)))),
                        ("arith mean", v -> mean(filter(!isnan, v))),
                        ("geo mean",   v -> geomean(filter(!isnan, v))))
        print(rpad(label, name_w), " | ")
        println(join([lpad(fmt(fn(M[:, j])), col_w) for j in 1:length(cols)], " | "))
    end
    # Wilcoxon vs baseline
    print(rpad("Wilcoxon p", name_w), " | ")
    ps = [j == baseline_col ? NaN :
          wilcoxon_signed_rank_p(filter(!isnan, M[:, j]), filter(!isnan, M[:, baseline_col]))
          for j in 1:length(cols)]
    println(join([lpad(isnan(p) ? "—" : @sprintf("%.3g", p), col_w) for p in ps], " | "))
    print(rpad("signif.", name_w), " | ")
    println(join([lpad(stars(p), col_w) for p in ps], " | "))
end

if abspath(PROGRAM_FILE) == @__FILE__
    folder = length(ARGS) >= 1 ? ARGS[1] : joinpath(REPO_ROOT, "settings", "5-partial_decomp")
    main(folder)
end
