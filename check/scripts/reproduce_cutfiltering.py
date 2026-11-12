#!/usr/bin/env python3
"""Run the six paper configurations and write a performance profile.

Python 3.9+ and Julia 1.12.6 are required. No Python packages are needed.
The default three instances illustrate the main comparison, not a new benchmark.
Run with --help for paths, the full instance set, and plot-only operation.
"""

import argparse
import csv
import hashlib
import html
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]
DEMO_CASES = ("pglib_scopf_case116", "pglib_scopf_case218_1", "pglib_scopf_case248")
METHODS = {
    "unfiltered": ("No filtering", "000000", ""),
    "scip": ("SCIP cut selection", "000000", "dashed"),
    "random": ("Random", "000000", "dotted"),
    "violation": ("Violation", "FE6100", ""),
    "diversity": ("Diversity", "DC267F", ""),
    "hybrid": ("Hybrid", "785EF0", ""),
}

# Each solve uses a fresh Julia process, as in the reported experiments.
WORKER = r'''
using TopologyOptimizer, JSON, LinearAlgebra
VERSION == v"1.12.6" || error("Use Julia 1.12.6 to reproduce the experiments")
LinearAlgebra.BLAS.set_num_threads(1)
casefile, parameter_file, output = ARGS
params = TopologyOptimizer.convert_params(JSON.parsefile(parameter_file))
consts = Dict{Symbol,Any}(:USER_ENVIRONMENT => Dict(:DIRECTORIES => Dict(
    :RESULTS_DIR => output, :LOG_DIR => output,
    :SOLVER_LOG_DIR => joinpath(output, "logs"))))
safe(x) = x
safe(x::AbstractFloat) = isfinite(x) ? x : nothing
safe(x::Union{Enum,Symbol,Type,Function}) = string(x)
safe(x::AbstractDict) = Dict(string(k) => safe(v) for (k,v) in x)
safe(x::AbstractArray) = map(safe, x)
save(name, value) = open(io -> JSON.print(io, safe(value), 2), joinpath(output, name), "w")
refs = (Ref{Cint}(), Ref{Cint}(), Ref{Cint}())
TopologyOptimizer.Gurobi.GRBversion(refs...)
solver_version = join(getindex.(refs), ".")
solver_version == "13.0.2" || error("Use Gurobi 13.0.2, found $solver_version")
save("versions.json", Dict("julia" => string(VERSION), "gurobi" => solver_version))
try
    result = TopologyOptimizer.main(casefile; additional_params=params,
        additional_consts=consts, plot=false, busbars=false, timesteps=1,
        slack_types=[:ratea])
    keys = ("termination_status", "primal_status", "solve_time", "objective", "objective_lb")
    small = Dict(k => get(result, k, nothing) for k in keys)
    small["callback"] = get(result, "callback", Dict())
    save("result.json", small)
finally
    save("parameters.json", TopologyOptimizer.PARAMS)
end
'''


def parameters(method, limit, threads):
    """Changes to settings/default.toml used in the main comparison."""
    selection = dict(ENABLED=method != "unfiltered", METHOD=":VIOLATION",
                     FEASCUT_ONLY=False, FEASCUT_INF_VIOL=method in ("violation", "hybrid"),
                     MIN_OPT_CUTS=0, VIOLATION=dict(LIMIT_ABS_NUM=10))
    if method == "scip":
        selection.update(METHOD=":SCIP_HYBRID", SCIP_HYBRID=dict(
            LIMIT_MODE=":scip_default", FORCE_FEASIBILITY=False, MIN_ORTHO=0.9,
            MAXCUTSROOT=2000, MAXCUTS=100, POST_EFFICACY_TOL=1e-6, SEED=42))
    elif method == "random":
        selection.update(METHOD=":RANDOM", RANDOM=dict(
            NUM_CUTS_MODE=":rel", NUM_CUTS_REL_ALPHA=0.05))
    elif method == "violation":
        selection["VIOLATION"].update(LIMIT_MODE=":rel_num", LIMIT_REL_ALPHA=0.05)
    elif method in ("diversity", "hybrid"):
        selection.update(METHOD=":DIVERSITY", DIVERSITY=dict(
            USE_DBSCAN=False, DISTANCE_METRIC=":cosine", NUM_CLUSTERS_MODE=":rel",
            NUM_CLUSTERS_REL_ALPHA=0.05,
            REPRESENTATIVE_MODE=":centroid" if method == "diversity" else ":max_violation"))
    return dict(APPLICATION="DC_SCOTS_TIME_SLACK_BENDERS",
        LIMITS=dict(BRANCH=dict(TAL_TIMESTEP=1)), TECHNICAL=dict(SEED=42),
        SOLVER=dict(SOLVING_TIME_LIMIT=limit, MULTITHREADING=dict(
            N_THREADS_MIP_MP=threads, N_THREADS_BENDERS=1, N_THREADS_MIP_SP=1)),
        BENDERS=dict(PARTIAL=dict(ENABLED=False),
            CUTS=dict(CUT_STAB_ROUND_FREQ=-1, SELECTION=selection)))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def finite(x):
    return isinstance(x, (int, float)) and math.isfinite(x)


def read_result(method, case, directory, returncode, limit):
    row = dict(method=method, case=case, status="ERROR", seconds=None,
               objective=None, bound=None, solved=False)
    path = directory / "result.json"
    if returncode or not path.is_file():
        return row
    result = json.loads(path.read_text())
    row.update(status=result.get("termination_status", "ERROR"),
               seconds=result.get("solve_time"), objective=result.get("objective"),
               bound=result.get("objective_lb"))
    if row["status"] == "OPTIMAL":
        if (not all(finite(row[k]) for k in ("seconds", "objective", "bound"))
                or row["seconds"] <= 0 or result.get("primal_status") != "FEASIBLE_POINT"
                or abs(row["objective"] - row["bound"]) > 1e-4 * max(1., abs(row["objective"]))):
            row["status"] = "INCONSISTENT"
        else:
            row["solved"] = row["seconds"] <= limit
    return row


def validate_pairs(rows):
    for case in dict.fromkeys(row["case"] for row in rows):
        solved = [row for row in rows if row["case"] == case and row["solved"]]
        if len(solved) > 1:
            tolerance = 1e-4 * max(1., *(abs(row["objective"]) for row in solved))
            if max(row["bound"] for row in solved) > min(row["objective"] for row in solved) + tolerance:
                raise ValueError(f"Inconsistent objective intervals for {case}; inspect the saved results")


def profile_points(values, n):
    values = sorted(v for v in values if math.isfinite(v) and v <= 16)
    y = sum(v <= 1 for v in values) / n
    points = [(1., y)]
    for v in values:
        if v > 1:
            points.append((v, y))
            y += 1 / n
            points.append((v, y))
    points.append((16., y))
    return points


def plots(rows, output):
    validate_pairs(rows)
    cases = list(dict.fromkeys(row["case"] for row in rows))
    methods = [m for m in METHODS if any(row["method"] == m for row in rows)]
    lookup = {(row["case"], row["method"]): row for row in rows}
    if len(lookup) != len(rows) or len(rows) != len(cases) * len(methods):
        raise ValueError("A profile requires exactly one result for every instance and method")
    best = {c: min((lookup[c, m]["seconds"] for m in methods if lookup[c, m]["solved"]),
                   default=math.inf) for c in cases}
    curves = {m: profile_points([lookup[c, m]["seconds"] / best[c]
        if lookup[c, m]["solved"] else math.inf for c in cases], len(cases)) for m in methods}
    tikz = [r"\begin{tikzpicture}",
        r"\begin{axis}[width=.95\linewidth,height=.58\linewidth,",
        r"xmin=1,xmax=16,ymin=0,ymax=1.025,xmode=log,log basis x=2,",
        r"xtick={1,2,4,8,16},xticklabels={1,2,4,8,16},ytick={0,.2,.4,.6,.8,1},",
        r"xlabel={Runtime ratio to fastest method, $\tau$},ylabel={Fraction of instances},",
        r"grid=major,grid style={black!15},tick align=outside,tick pos=left,",
        r"legend style={at={(.98,.04)},anchor=south east,draw=none,font=\small},",
        r"legend cell align=left,every axis plot/.append style={no marks,line width=1.2pt}]"]
    svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="800" height="500" viewBox="0 0 800 500">',
           '<rect width="800" height="500" fill="white"/>',
           '<g font-family="sans-serif" font-size="15" fill="black">',
           f'<text x="400" y="28" text-anchor="middle">Performance profile ({len(cases)} instances)</text>']
    for i in range(5):
        x = 85 + i * 165
        svg += [f'<path d="M{x} 50V410" stroke="#ddd"/>',
                f'<text x="{x}" y="435" text-anchor="middle">{2**i}</text>']
    for i in range(6):
        y = 410 - i * 72
        svg += [f'<path d="M85 {y}H745" stroke="#ddd"/>',
                f'<text x="73" y="{y+5}" text-anchor="end">{i/5:.1f}</text>']
    for m in methods:
        label, color, style = METHODS[m]
        tikz += [rf"\definecolor{{profile{m}}}{{HTML}}{{{color}}}",
                 rf"\addplot[profile{m},{style}] coordinates {{" +
                 " ".join(f"({x:.14g},{y:.14g})" for x,y in curves[m]) + "};",
                 rf"\addlegendentry{{{label}}}"]
        dash = {'dashed': 'stroke-dasharray="8 5"', 'dotted': 'stroke-dasharray="2 5"'}.get(style, '')
        points = " ".join(f"{85+165*math.log2(x):.3f},{410-360*y:.3f}" for x,y in curves[m])
        svg.append(f'<polyline points="{points}" fill="none" stroke="#{color}" stroke-width="2.5" {dash}/>')
    svg.append('<rect x="500" y="244" width="235" height="157" fill="white" fill-opacity=".95"/>')
    for i, m in enumerate(methods):
        label, color, style = METHODS[m]
        y = 266 + i*24
        dash = {'dashed': 'stroke-dasharray="8 5"', 'dotted': 'stroke-dasharray="2 5"'}.get(style, '')
        svg += [f'<path d="M510 {y-5}H548" stroke="#{color}" stroke-width="2.5" {dash}/>',
                f'<text x="560" y="{y}">{html.escape(label)}</text>']
    svg += ['<text x="415" y="475" text-anchor="middle">Runtime ratio to fastest method</text>',
            '<text transform="translate(23,230) rotate(-90)" text-anchor="middle">Fraction of instances</text>',
            '</g></svg>']
    (output / "performance_profile.svg").write_text("\n".join(svg) + "\n")
    (output / "performance_profile.tikz").write_text("\n".join(tikz + [r"\end{axis}", r"\end{tikzpicture}"]) + "\n")
    (output / "performance_profile.tex").write_text(
        "\\documentclass[border=3pt]{standalone}\n\\usepackage{pgfplots}\n"
        "\\pgfplotsset{compat=1.18}\n\\begin{document}\n"
        "\\input{performance_profile.tikz}\n\\end{document}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--instances", type=Path, help="Directory containing the supplied .m files")
    parser.add_argument("--output", type=Path, default=ROOT / "check/results/cutfilter-demo")
    parser.add_argument("--julia", default=shutil.which("julia"), help="Julia 1.12.6 executable")
    choice = parser.add_mutually_exclusive_group()
    choice.add_argument("--cases", nargs="+", default=list(DEMO_CASES), help="Instance names without .m")
    choice.add_argument("--all", action="store_true", help="Run all pglib_scopf_case*.m files in --instances")
    parser.add_argument("--methods", nargs="+", choices=METHODS, default=list(METHODS))
    parser.add_argument("--time-limit", type=float, default=600., help="Seconds per solve (default: 600; paper: 10800)")
    parser.add_argument("--threads", type=int, default=4, help="Master threads (paper: 4)")
    parser.add_argument("--plot-only", action="store_true", help="Regenerate plots from the saved runs.json")
    args = parser.parse_args()
    output = args.output.resolve()
    if args.plot_only:
        plots(json.loads((output / "runs.json").read_text()), output)
        return 0
    if not args.julia or not args.instances or not args.instances.is_dir():
        parser.error("Supply an instance directory and a Julia 1.12.6 executable")
    if not math.isfinite(args.time_limit) or args.time_limit <= 0 or args.threads < 1:
        parser.error("Time limit and thread count must be positive")
    cases = sorted(p.stem for p in args.instances.glob("pglib_scopf_case*.m")) if args.all else args.cases
    if not cases or len(set(cases)) != len(cases) or len(set(args.methods)) != len(args.methods):
        parser.error("Provide a nonempty instance set and no duplicate instances or methods")
    inputs = {c: (args.instances / f"{c}.m").resolve() for c in cases}
    for c, path in inputs.items():
        if Path(c).name != c or not path.is_file():
            parser.error(f"Missing or invalid instance: {c}")
    if output.exists():
        parser.error("Output directory already exists; choose a new directory or use --plot-only")
    output.mkdir(parents=True)
    sources = sorted((ROOT / "src").rglob("*.jl")) + [ROOT / "src/constants.toml",
        ROOT / "settings/default.toml", ROOT / "Project.toml", ROOT / "Manifest.toml"]
    write_json(output / "setup.json", dict(cases=cases, methods=args.methods,
        time_limit=args.time_limit, master_threads=args.threads, seed=42,
        description="Illustrative three-instance demonstration" if cases == list(DEMO_CASES) else "User-selected instances",
        input_sha256={c: digest(p) for c,p in inputs.items()},
        source_sha256={str(p.relative_to(ROOT)): digest(p) for p in sources},
        runner_sha256=digest(Path(__file__))))
    rows = []
    env = dict(os.environ, JULIA_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               OMP_NUM_THREADS="1", JULIA_PKG_PRECOMPILE_AUTO="0")
    for case in cases:
        for method in args.methods:
            directory = output / case / method
            directory.mkdir(parents=True)
            parameter_file = directory / "requested-parameters.json"
            write_json(parameter_file, parameters(method, args.time_limit, args.threads))
            print(f"[{len(rows)+1}/{len(cases)*len(args.methods)}] {case}: {METHODS[method][0]}", flush=True)
            command = [args.julia, "--startup-file=no", "--history-file=no", f"--project={ROOT}",
                       "-e", WORKER, str(inputs[case]), str(parameter_file), str(directory)]
            with (directory / "run.log").open("w") as log:
                result = subprocess.run(command, cwd=directory, env=env, stdout=log, stderr=subprocess.STDOUT)
            row = read_result(method, case, directory, result.returncode, args.time_limit)
            rows.append(row)
            write_json(output / "runs.json", rows)
            with (output / "runtimes.csv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(row))
                writer.writeheader()
                writer.writerows(rows)
            print(f"  {row['status']}, solve time {row['seconds']} s", flush=True)
            if row["status"] in ("ERROR", "INCONSISTENT"):
                print(f"  Inspect {directory / 'run.log'}", file=sys.stderr)
    plots(rows, output)
    print(f"Results and performance profiles: {output}")
    print("Unsuccessful runs have infinite ratios. This small demonstration does not replace the full benchmark.")
    return int(any(row["status"] in ("ERROR", "INCONSISTENT") for row in rows))


if __name__ == "__main__":
    sys.exit(main())
