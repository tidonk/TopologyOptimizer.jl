#!/usr/bin/env python3
"""Validate completed SCIP-hybrid runs against the matched paper controls.

Only derived analysis files are written. Production results and preflight
records are left unchanged. This is a recorded-result consistency audit, not
an independent mathematical certificate or a new solver experiment.
"""
import argparse
import collections
import concurrent.futures
import copy
import csv
import datetime
import hashlib
import json
import math
import re
import statistics
from pathlib import Path
import os

import export_full_cutfilter as export

METHOD = 'scip_hybrid_all_005'
NATIVE_METHOD = 'scip_hybrid_native_full'
ROOT = (Path(os.environ.get("TOPOLOGYOPTIMIZER_RESULTS_ROOT", Path(__file__).resolve().parents[2] / "check/results")))
DEFAULT_CAMPAIGN = ROOT / '20261026-cutfilter-scip-hybrid-efficacy-b8a05381'
DEFAULT_REFERENCE = ROOT / '20261023-cutfilter-full-005-fast-b8a05381'


def read(path):
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matched_control_solved(row):
    return (row['status'] == 'OPTIMAL' and row['primal_status'] == 'FEASIBLE_POINT'
        and row['objective_agrees_within_mip_tolerance'] and row['paired_bound_intervals_overlap']
        and row['reported_gap_valid'] and not row['inverted_bound'] and not row['incumbent_warning']
        and not row['n_benders_cut_value_fail'] and not row['n_questionable_sol'])


def export_eligible(row, method):
    # A comparison disagreement is distinct from a failed comparison. Keep it
    # visible in the native audit rather than interpreting it as a failed solve.
    return matched_control_solved(row) and (
        method == NATIVE_METHOD or not row['n_benders_cut_value_disagree'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign', type=Path, default=DEFAULT_CAMPAIGN)
    parser.add_argument('--reference', type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--method', choices=(METHOD, NATIVE_METHOD), default=METHOD)
    args = parser.parse_args()
    campaign, reference = args.campaign, args.reference
    provenance = read(campaign / 'source-provenance.json')
    ref_provenance = read(reference / 'source-provenance.json')
    method = args.method
    expected = [e for e in read(campaign / 'expected.json') if e['method'] == method]
    assert len(expected) == len({e['case'] for e in expected}) == 131
    assert {e['method'] for e in expected} == {method}
    reference_analysis = reference / 'analysis-paper-20261026-main-six-common'
    with (reference_analysis / 'records.csv').open() as stream:
        controls = {(r['case'], r['method']): r for r in csv.DictReader(stream)}
    assert {e['case'] for e in expected} == {key[0] for key in controls}
    common = set(read(reference_analysis / 'summary.json')['main_common_cases'])
    scripts = {name: digest(campaign / name) for name in ('worker.jl', 'launch.py')}
    assert scripts['worker.jl'] == digest(reference / 'worker.jl')
    if method == NATIVE_METHOD:
        preparation = read(campaign / 'preparation.json')
        assert preparation['source_archive_sha256'] == provenance['source_archive_sha256']
        assert all(value == preparation['hashes'][name] for name, value in scripts.items())
    else:
        assert scripts['launch.py'] == digest(reference / 'launch.py')

    # Verify both selector-only source changes, not just their descriptions.
    source_checks = []
    current = campaign
    while current != reference:
        p = read(current / 'source-provenance.json')
        parent = Path(p['reference_campaign'])
        assert parent in {reference, ROOT / '20261026-cutfilter-scip-hybrid-default-b8a05381',
                          ROOT / '20261026-cutfilter-scip-hybrid-efficacy-b8a05381',
                          ROOT / '20261027-cutfilter-grid-random20-s42'}
        assert digest(current / 'source.tar') == p['source_archive_sha256']
        assert digest(parent / 'source.tar') == p['reference_source_archive_sha256']
        new, old = export.archive_files(current / 'source.tar'), export.archive_files(parent / 'source.tar')
        changed = sorted(n for n in new.keys() | old.keys() if new.get(n) != old.get(n))
        assert changed == sorted(p['source_changed_files']), changed
        assert all(n.startswith(('src/core/', 'test/')) for n in changed)
        if 'source_overrides' in p:
            assert all(new[n] == p['source_overrides'][n] for n in changed)
        else:
            assert all(new[n] == digest(current / 'source_changes' / n) for n in changed)
        source_checks.append(dict(campaign=str(current), changed_files=changed))
        current = parent
    preflight = read(campaign / 'preflight/selector-validation.json')
    assert preflight['passed'] and preflight['source_archive_sha256'] == provenance['source_archive_sha256']
    for command in preflight.get('commands', []):
        assert command['returncode'] == 0
        assert digest(campaign / 'preflight' / (command['name'] + '.log')) == command['log_sha256']
    if method == NATIVE_METHOD:
        assert digest(campaign / 'preflight/selector-tests.log') == preflight['log_sha256']
        for name in ('settings-validation.json', 'smoke-validation.json'):
            check = read(campaign / 'preflight' / name)
            assert check['passed'] and check['source_archive_sha256'] == provenance['source_archive_sha256']

    def load(entry):
        case = entry['case']
        output = Path(entry['output'])
        assert output == campaign / 'results' / method / case
        assert digest(campaign / entry['instance']) == entry['instance_sha256'] == ref_provenance['input_sha256'][case + '.m']
        assert digest(campaign / entry['settings']) == provenance['settings_sha256'][method]
        row = export.load_record(entry, provenance, scripts)
        assert row['complete'] and row['returncode'] == 0, case
        result = read(output / 'result.json')
        params = read(output / 'resolved-params.json')
        original = read(reference / 'results/scip_dynamic_all_005' / case / 'resolved-params.json')
        selection = params['BENDERS']['CUTS']['SELECTION']
        assert selection['METHOD'] == 'SCIP_HYBRID'
        assert selection['MIN_OPT_CUTS'] == 0 and not selection['FEASCUT_INF_VIOL']
        expected_selection = dict(FORCE_FEASIBILITY=True, LIMIT_MODE='rel_num',
            LIMIT_REL_ALPHA=.05, MIN_ORTHO=.9, POST_EFFICACY_TOL=1e-6, SEED=42)
        if method == NATIVE_METHOD:
            expected_selection.update(FORCE_FEASIBILITY=False, LIMIT_MODE='scip_default',
                                      MAXCUTSROOT=2000, MAXCUTS=100)
        assert selection['SCIP_HYBRID'] == expected_selection
        normalized = copy.deepcopy(params)
        normalized['BENDERS']['CUTS']['SELECTION'].pop('SCIP_HYBRID')
        normalized['BENDERS']['CUTS']['SELECTION']['METHOD'] = 'SCIP_DYNAMIC'
        assert normalized == original, (case, 'Non-selector settings changed')
        runtime = read(output / 'runtime.json')
        assert runtime['success'] and runtime['julia_version'] == '1.12.6'
        assert runtime['gurobi_version'] == '13.0.2' and runtime['master_threads'] == 4
        assert result['callback']['SP_optimizer_type'] in ('InnOpTEM.LinEqSysRfSolver', 'TopologyOptimizer.LinEqSysRfSolver')
        assert all(math.isfinite(row[k]) for k in ('objective', 'bound', 'seconds'))
        assert row['seconds'] > 0
        text = (output / 'worker.log').read_text(errors='replace')
        tolerances = re.findall(r'Optimal solution found \(tolerance ([0-9.eE+-]+)\)', text)
        assert tolerances, case
        mip_tolerance = float(tolerances[-1])
        other = [controls[case, m] for m in ('violation_all_005', 'hybrid_all_005')]
        differences = [abs(row['objective'] - float(r['objective'])) for r in other]
        scale = max(abs(row['objective']), *(abs(float(r['objective'])) for r in other))
        numerical_tolerance = 1e-5 + 1e-6 * scale
        lb = max(row['bound'], *(float(r['bound']) for r in other))
        ub = min(row['objective'], *(float(r['objective']) for r in other))
        row.update(reported_mip_tolerance=mip_tolerance,
            objective_absolute_difference=max(differences),
            objective_relative_difference=max(differences) / max(scale, 1e-12),
            strict_objective_agreement=max(differences) <= numerical_tolerance,
            objective_agrees_within_mip_tolerance=max(differences) <= 1e-5 + mip_tolerance * scale,
            paired_bound_intervals_overlap=lb <= ub + numerical_tolerance,
            inverted_bound=row['bound'] > row['objective'] + numerical_tolerance,
            reported_gap_valid=row['objective'] - row['bound'] <= 1e-5 + mip_tolerance * abs(row['objective']),
            resolved_params_sha256=digest(output / 'resolved-params.json'),
            worker_log_sha256=digest(output / 'worker.log'))
        for name in ('n_benders_cut_value_disagree', 'n_benders_cut_value_fail', 'n_questionable_sol'):
            row[name] = result['callback'][name]
        for name, value in result['scip_hybrid_acceptance'].items():
            row['acceptance_' + name] = value
        row['matched_control_solved'] = matched_control_solved(row)
        row['validated_solved'] = row['matched_control_solved'] and not row['n_benders_cut_value_disagree']
        row['export_eligible'] = export_eligible(row, method)
        return row

    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        rows = list(pool.map(load, expected))
    diagnostics = dict(completed=len(rows), validated_solved=sum(r['validated_solved'] for r in rows),
        matched_control_solved=sum(r['matched_control_solved'] for r in rows),
        export_eligible=sum(r['export_eligible'] for r in rows),
        status_counts=dict(collections.Counter(r['status'] for r in rows)),
        strict_objective_agreements=sum(r['strict_objective_agreement'] for r in rows),
        objective_agreements_within_reported_mip_tolerance=sum(r['objective_agrees_within_mip_tolerance'] for r in rows),
        overlapping_bound_intervals=sum(r['paired_bound_intervals_overlap'] for r in rows),
        max_absolute_objective_difference=max(r['objective_absolute_difference'] for r in rows),
        max_relative_objective_difference=max(r['objective_relative_difference'] for r in rows),
        incumbent_warnings=sum(r['incumbent_warning'] for r in rows),
        cut_value_disagreements=sum(r['n_benders_cut_value_disagree'] for r in rows),
        cut_value_failures=sum(r['n_benders_cut_value_fail'] for r in rows),
        questionable_solutions=sum(r['n_questionable_sol'] for r in rows),
        acceptance_totals={k: sum(r.get(k, 0) for r in rows)
                           for k in sorted({k for r in rows for k in r if k.startswith('acceptance_')})})
    summaries = []
    for label, cases in [('All 131', {r['case'] for r in rows}), ('Common 126', common)]:
        rr = [r for r in rows if r['case'] in cases]
        summaries.append(dict(group=label, n=len(rr), mean_seconds=statistics.mean(r['seconds'] for r in rr),
            sgm_seconds=export.sgm([r['seconds'] for r in rr]),
            mean_callbacks=statistics.mean(r['callback_rounds'] for r in rr),
            mean_cuts=statistics.mean(r['added_benders_cuts'] for r in rr)))
    audit = dict(passed=all(r['export_eligible'] for r in rows),
        checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), campaign=str(campaign),
        method=method, script_sha256=scripts,
        reference_campaign=str(reference), reference_analysis=str(reference_analysis),
        source_sha256=provenance['source_archive_sha256'], expected_sha256=digest(campaign / 'expected.json'),
        source_checks=source_checks, diagnostics=diagnostics, summaries=summaries,
        checks=['All 131 identities, successful exits, inputs, source and settings hashes',
                'Resolved parameters identical to matched controls outside SCIP selector settings',
                'Objectives within reported MIP tolerance and paired bound intervals overlap',
                'No incumbent warnings, cut-value failures or questionable solutions'],
        diagnostic_policy=('Native results require matched objectives and bounds and no recorded cut-value failures. '
            'Cut-value comparison disagreements are counted separately and do not themselves classify a solve as failed. '
            'They are not asserted to be harmless, and independent solution or cut validity is not established.'
            if method == NATIVE_METHOD else 'All cut-value comparison and failure counters must be zero.'),
        scope='Recorded results and matched-control consistency, not an independent proof of optimality')
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    for name, values in [('records.csv', rows), ('summary.csv', summaries)]:
        with (args.output / name).open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(values[0]))
            writer.writeheader()
            writer.writerows(values)
    print(json.dumps(dict(passed=audit['passed'], diagnostics=diagnostics, summaries=summaries), indent=2))
    if not audit['passed']:
        raise SystemExit('Review flagged runs before updating the paper.')


if __name__ == '__main__':
    main()
