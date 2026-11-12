#!/usr/bin/env python3
"""Run two isolated TAL workers in one allocation, retaining reproducible evidence."""
from __future__ import annotations
import argparse
import hashlib
import fcntl
import getpass
import json
import math
import os
from pathlib import Path
import shutil
import signal
import re
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[3]


def settings_time_limit(path):
    # Only [SOLVER].SOLVING_TIME_LIMIT controls the master. Do not silently
    # substitute a short pilot limit for a settings file that inherits defaults.
    def read_limit(file):
        section = None
        found = []
        for line in Path(file).read_text().splitlines():
            line = line.split("#", 1)[0].strip()
            if line.startswith("["):
                section = line
            elif section == "[SOLVER]" and re.match(r"SOLVING_TIME_LIMIT\s*=", line):
                found.append(float(line.split("=", 1)[1].strip().replace("_", "")))
        if len(found) > 1:
            raise ValueError("Duplicate SOLVER.SOLVING_TIME_LIMIT")
        return found[0] if found else None
    value = read_limit(path)
    if value is None:
        value = read_limit(ROOT / "settings/default.toml")
    if value is None or not math.isfinite(value) or value <= 0:
        raise ValueError("A finite positive solver time limit is required")
    return value


def persist_events(mailbox, output, offsets):
    """Append complete new records once, rather than recopying growing files."""
    for path in mailbox.glob("*.events.csv"):
        offset = offsets.get(path.name, 0)
        with path.open("rb") as stream:
            stream.seek(offset)
            chunk = stream.read()
        end = chunk.rfind(b"\n") + 1
        if end:
            with (output / path.name).open("ab") as stream:
                stream.write(chunk[:end])
            offsets[path.name] = offset + end


def write_json_atomic(path, data):
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(data, indent=2) + "\n")
    temp.replace(path)


def prepare_receiver_depot(donor, receiver):
    """Seed the persistent node-local cache once; keep packages in donor fallback."""
    if donor == receiver:
        raise ValueError("Donor and receiver must use distinct primary depots")
    if not donor.is_dir():
        raise ValueError(f"Existing node-local Julia depot is missing: {donor}")
    receiver.mkdir(parents=True, exist_ok=True)
    with (receiver / ".seed.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        marker = receiver / ".seed-complete.json"
        if not marker.exists():
            compiled = donor / "compiled"
            if compiled.exists():
                shutil.copytree(compiled, receiver / "compiled", dirs_exist_ok=True,
                                ignore=shutil.ignore_patterns("*.pidfile", "jl_*"))
            write_json_atomic(marker, dict(donor=str(donor), seeded_epoch=time.time()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("settings", type=Path)
    parser.add_argument("experiment")
    parser.add_argument("--time-limit", type=float)
    parser.add_argument("--helper", choices=["tal1", "scopf", "rscots", "rscots-sequential"], default="tal1")
    parser.add_argument("--main-tal", type=int, default=2)
    parser.add_argument("--score", choices=["RND", "FULL", "SENS", "CONGPROX", "ETA", "FIXDUAL", "ADMITTANCE_SECANT"], default="SENS")
    parser.add_argument("--switchable", type=int, default=16)
    parser.add_argument("--threads", type=int, nargs=2, default=[4, 4])
    parser.add_argument("--no-sharing", action="store_true")
    parser.add_argument("--stop-helper-with-main", action="store_true")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--audit-models", action="store_true")
    parser.add_argument("--startup-timeout", type=float, default=3600)
    parser.add_argument("--donor-depot", type=Path)
    parser.add_argument("--receiver-depot", type=Path)
    args = parser.parse_args()
    if args.time_limit is None:
        args.time_limit = settings_time_limit(args.settings)
    if not math.isfinite(args.time_limit) or args.time_limit <= 0 or min(args.threads) < 1 or sum(args.threads) > 8:
        parser.error("Use a positive time limit and at most eight total threads")
    if not math.isfinite(args.startup_timeout) or args.startup_timeout <= 0:
        parser.error("Use a finite positive startup timeout")
    if args.main_tal < 1 or args.switchable < 1:
        parser.error("Main TAL and switchable branch count must be positive")
    if args.helper == "rscots-sequential" and args.main_tal == 1:
        # Preserve the previously tested single-stage helper at TAL=1.
        args.helper = "rscots"
    case, settings = args.case.resolve(), args.settings.resolve()
    output = args.output or ROOT / "check/results" / args.experiment / settings.stem / case.stem
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if (output / "manifest.json").exists():
        parser.error(f"Refusing to overwrite an existing run: {output}")
    try:
        with (output / "run.claim").open("x") as claim:
            claim.write(f"pid={os.getpid()} host={os.uname().nodename}\n")
    except FileExistsError:
        parser.error(f"Another launch already claimed this run: {output}")
    # Both processes are on one machine. The mailbox is local, private and fresh;
    # only logs/results are written to the shared experiment directory.
    local = Path(tempfile.mkdtemp(prefix="innoptem-pair-"))
    mailbox = local / "exchange"
    mailbox.mkdir()
    julia = os.environ.get("INNOPTEM_JULIA", shutil.which("julia") or "julia")
    depot = os.environ.get("INNOPTEM_BASE_DEPOT", os.environ.get("JULIA_DEPOT_PATH", str(Path.home() / ".julia")))
    donor_depot = (args.donor_depot or Path(depot.split(os.pathsep)[0])).expanduser().resolve()
    receiver_depot = (args.receiver_depot or donor_depot.with_name(donor_depot.name + "-tal2")).expanduser().resolve()
    if donor_depot == receiver_depot:
        parser.error("Donor and receiver must use distinct primary depots")
    prepare_receiver_depot(donor_depot, receiver_depot)
    fallback_depot = os.environ.get("INNOPTEM_FALLBACK_DEPOT", str(Path.home() / ".julia"))
    worker = ROOT / "check/scripts/cluster/parallel_worker.jl"
    roles = ["donor", "control" if args.no_sharing else "receiver"]
    revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True)
    manifest = dict(case=str(case), settings=str(settings), threads=args.threads,
                    time_limit=args.time_limit, roles=roles, julia=julia,
                    helper=args.helper, main_tal=args.main_tal, score=args.score if args.helper.startswith("rscots") else None,
                    switchable=args.switchable if args.helper.startswith("rscots") else None,
                    main_warmstart=False if args.helper != "tal1" else None,
                    stop_helper_with_main=args.stop_helper_with_main,
                    hostname=os.uname().nodename, start_epoch=time.time(),
                    startup_timeout=args.startup_timeout,
                    donor_depot=str(donor_depot), receiver_depot=str(receiver_depot),
                    git_head=revision.stdout.strip() or "snapshot",
                    source_root=str(ROOT),
                    case_sha256=hashlib.sha256(case.read_bytes()).hexdigest(),
                    settings_sha256=hashlib.sha256(settings.read_bytes()).hexdigest(),
                    clock="Shared elapsed time since both workers reach their first optimize call.")
    files = [str(p.relative_to(ROOT)) for folder in ("src", "ext")
             for p in (ROOT / folder).rglob("*") if p.is_file() and p.suffix == ".jl"]
    files += ["Project.toml", "Manifest.toml", "settings/default.toml",
              "check/scripts/cluster/run_parallel_tal.py", "check/scripts/cluster/parallel_worker.jl"]
    manifest["source_sha256"] = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()
                                  for p in files}
    write_json_atomic(output / "manifest.json", manifest)
    shutil.copy2(settings, output / "settings.toml")
    procs, streams = [], []
    offsets = {}
    last_persist = 0.0

    def terminate(signum=None, frame=None):
        for proc in procs:
            if proc.poll() is None:
                try:
                    os.killpg(proc.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
        if signum is not None:
            raise SystemExit(128 + signum)

    signal.signal(signal.SIGTERM, terminate)
    signal.signal(signal.SIGINT, terminate)
    try:
        for role, threads in zip(roles, args.threads):
            role_out = output / role
            role_out.mkdir()
            (mailbox / f"{role}.events.csv").write_text("time_s,event,primal,dual,sequence,detail\n")
            env = dict(os.environ, INNOPTEM_EXCHANGE_DIR=str(mailbox),
                       INNOPTEM_EXCHANGE_ROLE=role,
                       INNOPTEM_EXCHANGE_HELPER=args.helper,
                       INNOPTEM_EXCHANGE_MAIN_TAL=str(args.main_tal),
                       INNOPTEM_EXCHANGE_SCORE=args.score,
                       INNOPTEM_EXCHANGE_SWITCHABLE=str(args.switchable),
                       INNOPTEM_EXCHANGE_STOP_WITH_MAIN="1" if args.stop_helper_with_main else "0",
                       INNOPTEM_EXCHANGE_AUDIT="1" if args.audit_models else "0",
                       JULIA_DEPOT_PATH=(f"{donor_depot}:{fallback_depot}:" if role == "donor"
                                         else f"{receiver_depot}:{donor_depot}:{fallback_depot}:"),
                       INNOPTEM_STARTUP_TIMEOUT=str(args.startup_timeout),
                       JULIA_PKG_PRECOMPILE_AUTO="0", JULIA_NUM_PRECOMPILE_TASKS="1",
                       OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
            # Each worker has one Julia thread for separation and its own Gurobi
            # master thread budget. Gurobi pauses its search during callbacks.
            if args.helper != "tal1":
                env["INNOPTEM_EXCHANGE_LIMIT"] = str(args.time_limit)
            command = [julia, "--startup-file=no", "--history-file=no", "--compiled-modules=existing", "--threads=1",
                       f"--project={ROOT}", str(worker), str(case), str(settings),
                       str(role_out), role, str(args.time_limit), str(threads)]
            log = open(role_out / "worker.log", "w")
            streams.append(log)
            scratch = local / role
            (scratch / "check/results/logs").mkdir(parents=True)
            procs.append(subprocess.Popen(command, cwd=scratch, env=env, stdout=log,
                                          stderr=subprocess.STDOUT, start_new_session=True))
        deadline = time.monotonic() + args.startup_timeout
        # A helper is optional: even a failure before its barrier must not
        # prevent the main solver from starting. Main-worker failures stay fatal.
        while not (mailbox / f"{roles[1]}.ready").exists() or (
                procs[0].poll() is None and not (mailbox / "donor.ready").exists()):
            if procs[1].poll() is not None:
                raise RuntimeError("Main worker exited before the start barrier; inspect worker.log")
            if time.monotonic() > deadline:
                raise TimeoutError("Workers did not reach the start barrier")
            time.sleep(0.1)
        manifest["solve_epoch"] = time.time()
        (mailbox / "epoch.tmp").write_text(str(manifest["solve_epoch"]))
        (mailbox / "epoch.tmp").replace(mailbox / "epoch")
        deadline = time.time() + (2 if args.helper == "tal1" else 1) * args.time_limit + 900
        while any(p.poll() is None for p in procs):
            if procs[1].poll() not in (None, 0):
                raise RuntimeError("Main worker failed; inspect worker.log")
            if time.time() > deadline:
                raise TimeoutError("Paired solve exceeded its watchdog")
            # Persist intermediate evidence, including when Condor terminates.
            if time.time() - last_persist >= 2.0:
                persist_events(mailbox, output, offsets)
                last_persist = time.time()
            time.sleep(0.25)
        manifest["returncodes"] = [p.returncode for p in procs]
        if manifest["returncodes"][1] != 0:
            raise RuntimeError("Main worker failed; inspect worker.log")
        manifest["main_complete"] = True
        manifest["complete"] = manifest["returncodes"][0] == 0
        stopped = output / "donor" / "helper-stopped.json"
        if stopped.exists():
            manifest["helper_stopped_with_main"] = json.loads(stopped.read_text())
        if not manifest["complete"]:
            manifest["helper_error"] = "Helper failed; main solve completed independently. Inspect donor/worker.log."
    except BaseException as exc:
        manifest["error"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        terminate()
        for proc in procs:
            try:
                proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
        for stream in streams:
            stream.close()
        persist_events(mailbox, output, offsets)
        for pattern in ("*.bin", "*.lp", "*.sol"):
            for path in mailbox.glob(pattern):
                shutil.copyfile(path, output / path.name)
        manifest["returncodes"] = [p.returncode for p in procs]
        manifest["end_epoch"] = time.time()
        write_json_atomic(output / "manifest.json", manifest)
        shutil.rmtree(local)
    # The dispatcher marker is published only after final evidence/manifest.
    if not manifest["complete"]:
        # A usable main result is not a successful paired experiment. Retain
        # its diagnostics, but do not publish the dispatcher's PAIR_COMPLETE.
        print(output)
        return
    marker = output.parent / f"{case.stem}.{settings.stem}.csv"
    temp = marker.with_suffix(".csv.tmp")
    temp.write_text(f"instance_file,settings,status\n{case.stem},{settings.stem},PAIR_COMPLETE\n")
    temp.replace(marker)
    print(output)


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--print-watchdog":
        print(int(2 * settings_time_limit(sys.argv[2]) + 3600 + 900))
    else:
        main()
