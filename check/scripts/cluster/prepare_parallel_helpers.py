#!/usr/bin/env python3
"""Freeze the exp5 parallel-helper benchmark and prepare (but do not submit) Condor jobs."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[3]
METHODS = ['scopf', 'rscots-rnd', 'rscots-full', 'rscots-sens', 'rscots-congprox', 'rscots-eta']

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('campaign', type=Path)
    p.add_argument('--baseline', type=Path, required=True)
    p.add_argument('--priority', type=int, default=100)
    args = p.parse_args()
    dest = args.campaign.resolve()
    dest.mkdir(parents=True, exist_ok=False)
    source = dest / 'source'
    source.mkdir()
    for name in ['src', 'ext']:
        if (ROOT/name).exists():
            shutil.copytree(ROOT/name, source/name)
    for name in ['Project.toml', 'Manifest.toml', 'settings/default.toml',
                 'check/scripts/cluster/run_parallel_tal.py', 'check/scripts/cluster/parallel_worker.jl']:
        (source/name).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT/name, source/name)
    # Required by the package's optional precompile workload; no jobs write here.
    source.joinpath('check/instances').symlink_to(ROOT/'check/instances', target_is_directory=True)
    settings = dest/'exp5.toml'
    shutil.copy2(args.baseline/'settings/partialdec.toml', settings)
    cases = (args.baseline/'active-cohort.txt').read_text().splitlines()
    assert len(cases) == len(set(cases)) == 131 and 'pglib_scopf_case197_2' not in cases
    (dest/'active-cohort.txt').write_text('\n'.join(cases)+'\n')
    (dest/'instances').mkdir()
    for case in cases:
        shutil.copy2(args.baseline/'instances'/f'{case}.m', dest/'instances'/f'{case}.m')
    for name in ['results', 'cluster_logs']:
        (dest/name).mkdir()
    revision = subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert revision == subprocess.check_output(['git','rev-parse','origin/main'],cwd=ROOT,text=True).strip()
    patch = subprocess.check_output(['git','diff','HEAD','--','src','Project.toml'],cwd=ROOT)
    (dest/'source.patch').write_bytes(patch)
    files = {str(f.relative_to(source)):digest(f) for f in source.rglob('*') if f.is_file() and not f.is_symlink()}
    spec = dict(revision=revision, source_checkout=str(ROOT), source_files=files,
                patch_sha256=hashlib.sha256(patch).hexdigest(),
                settings_sha256=digest(settings), settings_reference=str(args.baseline/'settings/partialdec.toml'),
                input_sha256={case:digest(dest/'instances'/f'{case}.m') for case in cases},
                methods=METHODS, jobs=len(cases)*len(METHODS), main_tal=1, helper_tal={'scopf':0,'rscots':1},
                switchable_branches=16, main_warmstart=False, rscots_scopf_initialization=True,
                time_limit=10800, clock='Shared elapsed wall time from first optimize barrier; SCOPF and scoring count inside the helper budget',
                threads={'main':4,'helper':4,'julia_per_process':1}, memory_MB=[14000,20000], max_retries=1,
                priority=args.priority,
                selection='Top 16 by descending raw score, ascending branch id for ties; FULL is constant-score top 16',
                )
    (dest/'source-provenance.json').write_text(json.dumps(spec,indent=2)+'\n')
    shutil.copy2(Path(__file__), dest/'prepare.py')
    (dest/'launch.py').write_text('''#!/usr/bin/env python3
import argparse, json, os, shutil, subprocess, time, uuid
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('case');p.add_argument('method');p.add_argument('--limit',type=float,default=10800)
p.add_argument('--threads',type=int,default=4)
a=p.parse_args();root=Path(__file__).resolve().parent
out=root/'results'/a.method/a.case
if (out/'manifest.json').exists():
    old=json.loads((out/'manifest.json').read_text())
    if old.get('complete') and (out/'receiver/result.json').is_file():
        raise SystemExit(0)
if out.exists():
    archive=root/'retry-archive'/a.method/a.case
    archive.mkdir(parents=True,exist_ok=True)
    out.rename(archive/(time.strftime('%Y%m%dT%H%M%S',time.gmtime())+'-'+uuid.uuid4().hex[:8]))
cmd=['/usr/bin/python3',str(root/'source/check/scripts/cluster/run_parallel_tal.py'),
     str(root/'instances'/f'{a.case}.m'),str(root/'exp5.toml'),root.name,
     '--helper','scopf' if a.method=='scopf' else 'rscots','--main-tal','1',
     '--switchable','16','--time-limit',str(a.limit),'--threads',str(a.threads),str(a.threads),
     '--output',str(out)]
if a.method!='scopf': cmd+=['--score',a.method.split('-',1)[1].upper()]
raise SystemExit(subprocess.call(cmd))
''')
    lines = [
        'universe = vanilla', 'executable = /usr/bin/python3', 'transfer_executable = false',
        f'initialdir = {dest}', f'arguments = "{dest}/launch.py $(case) $(method)"',
        'getenv = HOME,GRB_LICENSE_FILE,INNOPTEM_JULIA,INNOPTEM_BASE_DEPOT,TMPDIR',
        'request_cpus = 8',
        'request_memory = ifThenElse(((NumJobStarts =!= UNDEFINED) && (NumJobStarts > 1)) || ((NumJobCompletions =!= UNDEFINED) && (NumJobCompletions > 0)), 20000, 14000)',
        'max_retries = 1', 'success_exit_code = 0', f'priority = {args.priority}', 'requirements = True',
        f'+Campaign = "{dest.name}"', '+ExperimentCase = "$(case)"', '+ExperimentMethod = "$(method)"',
        '+TargetTAL = 1', '+SolverTimeLimit = 10800', '+HelperSwitchableBranches = 16',
        '+ParallelHelper = True', '+ValidationPending = True', f'+SourceRevision = "{revision}"',
        f'batch_name = {dest.name}',
        *[f'{key} = {dest}/cluster_logs/$(case).$(method).{ext}' for key,ext in [('output','out'),('error','err'),('log','log')]],
        'queue case,method from (']
    # Interleave settings, and start with small and formerly disagreeing cases.
    first=['pglib_scopf_case112','pglib_scopf_case105','pglib_scopf_case147','pglib_scopf_case191_2','pglib_scopf_case115']
    ordered=first+[c for c in cases if c not in first]
    lines += [f'{case} {method}' for case in ordered for method in METHODS]+[')']
    (dest/'jobs.sub').write_text('\n'.join(lines)+'\n')
    print(json.dumps(dict(campaign=str(dest),jobs=spec['jobs'],revision=revision,submitted=False)))

if __name__=='__main__': main()
