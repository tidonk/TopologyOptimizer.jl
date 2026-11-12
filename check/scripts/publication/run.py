#!/usr/bin/env python3
"""Run a paper configuration in a fresh Julia process (Python 3.11+)."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import tomllib

ROOT = Path(__file__).resolve().parents[3]
SETTINGS = ROOT / 'settings/publication'
MAIN_SETTINGS = {'compact': '1-compact.toml', 'benders': '2-benders.toml', 'benders-sleu': '3-benders-sleu.toml', 'spsel': '4a-sp-selection.toml', 'cutfilt': '4b-cut-filtering.toml', 'spselcutfilt': '5-sp-cf-combined.toml', 'partialdec': '6-partial-decomposition.toml', 'branchfix': '7-branch-fixing.toml', 'sequential': '8a-heuristic-warmstart.toml', 'hybrid': '8b-parallel-heuristic.toml', 'benders-threads': '9-multithreading.toml'}


def merge(a, b):
    result = dict(a)
    for k, v in b.items():
        result[k] = merge(result.get(k, {}), v) if isinstance(v, dict) else v
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--suite', choices=['main', 'practical'], required=True)
    parser.add_argument('--method', required=True)
    parser.add_argument('--tal', type=int, choices=range(1, 6))
    parser.add_argument('--instances', type=Path, required=True)
    parser.add_argument('--cases', nargs='+', required=True, help='Case names without .m, or all')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--julia', default=os.environ.get('INNOPTEM_JULIA', shutil.which('julia') or 'julia'))
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    if (args.suite == 'practical') != (args.tal is not None):
        parser.error('--tal is required only for the practical suite')
    name = args.method + (f'-tal{args.tal}' if args.tal else '')
    if args.suite == 'main':
        if args.method not in MAIN_SETTINGS:
            parser.error(f'Unknown configuration: {name}')
        settings = ROOT / 'settings' / MAIN_SETTINGS[args.method]
    else:
        settings = SETTINGS / 'practical' / f'{name}.toml'
    if not settings.is_file():
        parser.error(f'Unknown configuration: {name}')
    hashes = {line.split()[1][:-2]: line.split()[0] for line in (SETTINGS/'instances.sha256').read_text().splitlines()}
    cases = sorted(hashes) if args.cases == ['all'] else args.cases
    if len(cases) != len(set(cases)):
        parser.error('Repeated cases')
    defaults = tomllib.loads((ROOT/'settings/default.toml').read_text())
    effective = merge(defaults, tomllib.loads(settings.read_text()))
    parallel = args.method in ('hybrid', 'parallel')
    env = dict(os.environ, INNOPTEM_JULIA=args.julia, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
    for case in cases:
        casefile = (args.instances / (case + '.m')).resolve()
        if case not in hashes or hashlib.sha256(casefile.read_bytes()).hexdigest() != hashes[case]:
            parser.error(f'Case does not match the published input: {case}')
        output = (args.output / name / case).resolve()
        if output.exists():
            parser.error(f'Refusing to overwrite {output}')
        if parallel:
            threads = str(effective['SOLVER']['MULTITHREADING']['N_THREADS_MIP_MP'])
            restricted = effective['WARMSTARTING']['RSCOTS']
            scores = [key for key, weight in restricted['WEIGHTS'].items() if weight != 0]
            if len(scores) != 1:
                parser.error('The parallel launcher requires one active scoring method')
            command = [sys.executable, str(ROOT/'check/scripts/cluster/run_parallel_tal.py'),
                       str(casefile), str(settings), name, '--helper', 'rscots-sequential',
                       '--main-tal', str(effective['LIMITS']['BRANCH']['TAL_TIMESTEP']),
                       '--score', scores[0], '--switchable', str(restricted['N_SW_BR_ABS']),
                       '--threads', threads, threads, '--stop-helper-with-main', '--output', str(output)]
        else:
            command = [args.julia, '--startup-file=no', '--history-file=no', '--compiled-modules=existing',
                       f"--threads={effective['SOLVER']['MULTITHREADING']['N_THREADS_BENDERS']}",
                       f'--project={ROOT}', str(Path(__file__).with_name('worker.jl')),
                       str(casefile), str(settings), str(output)]
        print(json.dumps(dict(case=case, command=command)), flush=True)
        if args.dry_run:
            continue
        if parallel:
            subprocess.run(command, cwd=ROOT, env=env, check=True)
            continue
        output.mkdir(parents=True)
        files = [p for folder in ('src', 'ext') for p in (ROOT/folder).rglob('*.jl')]
        files += [ROOT/p for p in ('Project.toml', 'Manifest.toml', 'settings/default.toml',
                                   'check/scripts/publication/worker.jl', 'check/scripts/publication/run.py')]
        files.append(settings)
        record = dict(command=command, case_sha256=hashes[case],
                      git_head=subprocess.check_output(['git','rev-parse','HEAD'], cwd=ROOT, text=True).strip(),
                      source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
        (output/'run.json').write_text(json.dumps(record, indent=2)+'\n')
        shutil.copy2(settings, output/'settings.toml')
        with tempfile.TemporaryDirectory(prefix='innoptem-paper-') as scratch, (output/'worker.log').open('w') as log:
            (Path(scratch)/'check/results/logs').mkdir(parents=True)
            subprocess.run(command, cwd=scratch, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        if not (output/'result.json').is_file():
            raise RuntimeError(f'Missing result: {output}')


if __name__ == '__main__':
    main()
