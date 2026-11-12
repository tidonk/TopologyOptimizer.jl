#!/usr/bin/env python3
"""Export the frozen eight-method campaign; retain failures in the full test set.

Objective screening is deliberately the historical 1% consensus rule, not a
certificate of correctness. Unfinished runs remain explicitly provisional.
"""
import argparse
import collections
import concurrent.futures
import csv
import datetime
import hashlib
import json
import math
import statistics
import tarfile
from pathlib import Path

METHODS = ['unfiltered', 'random_all_005', 'violation_all_005',
           'diversity_all_005', 'hybrid_all_005', 'hybrid_plus_all_005',
           'scip_dynamic_all_005', 'violation_pb_110']
PB_REPLACEMENT = 'violation_pb_allfeas_110'
LABELS = ['No filtering', 'Random', 'Violation', 'Diversity', 'Hybrid',
          'Hybrid+', 'SCIP cut selection', 'PB']
# Comparison methods use black; the other methods use the IBM palette.
# Hybrid+ and PB share their base method's color.
COLORS = ['000000', '000000', 'FE6100', 'DC267F', '785EF0', '785EF0',
          '000000', 'FE6100']
PLOT_STYLES = {'scip_dynamic_all_005': 'dashed', 'random_all_005': 'dotted'}
LIMIT = 10800
REFERENCE = 'scip_dynamic_all_005'
FAST_METHODS = {'hybrid_all_005', 'hybrid_plus_all_005', 'diversity_all_005', 'scip_dynamic_all_005'}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def instance_label(case):
    return r'\texttt{' + case.split('case')[1].replace('_', r'\_') + '}'


def sgm(values):
    return math.exp(statistics.mean(math.log(x + 10) for x in values)) - 10


def load_record(e, provenance, script_hashes):
    p = Path(e['output'])
    raw = (p / 'manifest.json').read_bytes()
    m = json.loads(raw)
    expected = dict(method=e['method'], case=e['case'],
                    source_sha256=provenance['source_archive_sha256'],
                    input_sha256=e['instance_sha256'],
                    settings_sha256=provenance['settings_sha256'][e['method']],
                    worker_sha256=script_hashes['worker.jl'],
                    launcher_sha256=script_hashes['launch.py'],
                    master_threads=4, time_limit=LIMIT, seed=42, tal=1,
                    julia_threads=1)
    for key, value in expected.items():
        assert m[key] == value, (e['method'], e['case'], key)
    row = dict(method=e['method'], case=e['case'], complete=m['complete'],
               returncode=m.get('returncode'), status='RUNNING', seconds=None,
               objective=None, bound=None, primal_status=None,
               manifest_sha256=sha(raw), result_sha256=None,
               callback_rounds=None, added_benders_cuts=None,
               incumbent_warning=False, kmedoids_fallback=False)
    if m['complete']:
        raw = (p / 'result.json').read_bytes()
        r = json.loads(raw)
        row.update(status=r['termination_status'], seconds=r['solve_time'],
                   objective=r.get('objective'), bound=r.get('objective_lb'),
                   primal_status=r.get('primal_status'), result_sha256=sha(raw))
        row.update(callback_rounds=r.get('callback', {}).get('n_callbacks'),
                   added_benders_cuts=r.get('callback', {}).get('n_benders_cuts_added'))
        assert m['returncode'] == 0, (e['method'], e['case'], m['returncode'])
        assert json.loads((p / 'runtime.json').read_text())['success']
        log = (p / 'worker.log').read_text(errors='replace')
        row['incumbent_warning'] = 'New incumbent solution accepted even though' in log
        row['kmedoids_fallback'] = 'falling back to top-k by violation' in log
    return row


def archive_files(path):
    with tarfile.open(path) as archive:
        return {m.name.removeprefix('./'): sha(archive.extractfile(m).read())
                for m in archive.getmembers() if m.isfile()}


def paired_metrics(reference, candidate):
    assert reference.keys() == candidate.keys()
    result = {'instances': len(reference)}
    for metric in ('seconds', 'callback_rounds', 'added_benders_cuts'):
        pairs = [(reference[c][metric], candidate[c][metric]) for c in reference]
        assert all(a is not None and b is not None for a, b in pairs)
        positive = [(a, b) for a, b in pairs if a > 0]
        result[metric] = dict(reference_mean=statistics.mean(a for a, _ in pairs),
            candidate_mean=statistics.mean(b for _, b in pairs),
            candidate_lower=sum(b < a for a, b in pairs),
            equal=sum(b == a for a, b in pairs), candidate_higher=sum(b > a for a, b in pairs),
            relative_reduction_instances=len(positive),
            mean_instance_relative_reduction=statistics.mean((a-b)/a for a, b in positive) if positive else None,
            reduction_of_means=1-statistics.mean(b for _, b in pairs)/statistics.mean(a for a, _ in pairs)
                if statistics.mean(a for a, _ in pairs) > 0 else None)
    return result


def step_points(ratios, xmax, n):
    counts = collections.Counter(x for x in ratios if math.isfinite(x))
    y = counts.pop(1.0, 0) / n
    points = [(1., y)]
    for x, count in sorted(counts.items()):
        points.append((x, y))
        y += count / n
        points.append((x, y))
    points.append((xmax, y))
    assert abs(y - sum(math.isfinite(x) for x in ratios) / n) < 1e-12
    return points


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--campaign', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--paper-output', type=Path, required=True)
    ap.add_argument('--exclude-method', action='append', default=[], choices=METHODS)
    ap.add_argument('--main-exclude-method', action='append', default=[], choices=METHODS + [PB_REPLACEMENT])
    ap.add_argument('--all-only', action='store_true', help='show only all-instance runtime statistics in the main table')
    ap.add_argument('--common-only', action='store_true', help='show statistics on instances solved by every main method; retain full-set solved counts and profile')
    ap.add_argument('--with-counts', action='store_true', help='append mean callback/cut counts; requires --all-only or --common-only')
    ap.add_argument('--exclude-case', action='append', default=[])
    ap.add_argument('--method-order', nargs='+', choices=METHODS + [PB_REPLACEMENT])
    ap.add_argument('--pb-campaign', type=Path,
                    help='replace PB results with the completed all-feasibility-cuts campaign')
    ap.add_argument('--replacement-campaign', type=Path,
                    help='replace the four affected selectors with a completed, matched 528-run campaign')
    ap.add_argument('--scip-campaign', type=Path,
                    help='replace the SCIP comparison with the validated hybrid-selector campaign')
    ap.add_argument('--scip-audit', type=Path,
                    help='completed-result audit.json produced by audit_scip_hybrid_completed.py')
    args = ap.parse_args()
    if args.all_only and args.common_only:
        ap.error('--all-only and --common-only are mutually exclusive')
    if args.with_counts and not (args.all_only or args.common_only):
        ap.error('--with-counts requires --all-only or --common-only')
    if bool(args.scip_campaign) != bool(args.scip_audit):
        ap.error('--scip-campaign and --scip-audit must be supplied together')
    c, out, paper = args.campaign, args.output, args.paper_output
    out.mkdir(parents=True, exist_ok=True)
    paper.mkdir(parents=True, exist_ok=True)
    expected_raw = (c / 'expected.json').read_bytes()
    expected = json.loads(expected_raw)
    provenance = json.loads((c / 'source-provenance.json').read_text())
    scripts = json.loads((c / 'submission-remaining.json').read_text())['script_sha256']
    assert len(expected) == len({(e['method'], e['case']) for e in expected}) == 1056
    assert set(args.exclude_case) <= {e['case'] for e in expected}
    selected = [(m, label, color) for m, label, color in
                zip(METHODS, LABELS, COLORS) if m not in args.exclude_method]
    if args.pb_campaign:
        assert 'pglib_scopf_case245' in args.exclude_case
        selected = [(PB_REPLACEMENT if m == 'violation_pb_110' else m,
                     label, color) for m, label, color in selected]
    if args.method_order:
        assert len(args.method_order) == len(set(args.method_order)) == len(selected)
        assert set(args.method_order) == {s[0] for s in selected}
        selected.sort(key=lambda s: args.method_order.index(s[0]))
    methods, labels, colors = map(list, zip(*selected))
    main_methods = [m for m in methods if m not in args.main_exclude_method]
    assert main_methods, 'Main comparison must include at least one method'
    assert REFERENCE in methods
    base_expected = [e for e in expected if e['method'] in methods and e['case'] not in args.exclude_case]
    records_to_load = [(e, provenance, scripts) for e in base_expected]
    fast_rows = None
    fast_metadata = None
    if args.replacement_campaign:
        rc = args.replacement_campaign
        replacement_raw = (rc / 'expected.json').read_bytes()
        replacement_expected = json.loads(replacement_raw)
        rp = json.loads((rc / 'source-provenance.json').read_text())
        rs = {name: sha((rc / name).read_bytes()) for name in ('worker.jl', 'launch.py')}
        assert len(replacement_expected) == 528
        assert set(rp['methods']) == FAST_METHODS
        identities = {(e['method'], e['case']) for e in replacement_expected}
        assert len(identities) == 528
        assert identities == {(m, e['case']) for m in FAST_METHODS for e in expected}
        assert sha((rc / 'source.tar').read_bytes()) == rp['source_archive_sha256']
        assert sha((c / 'source.tar').read_bytes()) == provenance['source_archive_sha256']
        original_files, fast_files = archive_files(c / 'source.tar'), archive_files(rc / 'source.tar')
        changed = sorted(n for n in original_files.keys() | fast_files.keys()
                         if original_files.get(n) != fast_files.get(n))
        assert changed == ['src/core/cutsel.jl'], changed
        for name, digest in rs.items():
            assert digest == scripts[name] == sha((c / name).read_bytes()), name
        for method in FAST_METHODS:
            assert (rp['settings_sha256'][method] == provenance['settings_sha256'][method]
                    == sha((rc / 'settings' / (method + '.toml')).read_bytes())
                    == sha((c / 'settings' / (method + '.toml')).read_bytes()))
        for e in replacement_expected:
            assert e['instance_sha256'] == provenance['input_sha256'][e['case'] + '.m'] == rp['input_sha256'][e['case'] + '.m']
            assert sha((rc / e['instance']).read_bytes()) == e['instance_sha256']
        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
            fast_rows = list(pool.map(lambda e: load_record(e, rp, rs), replacement_expected))
        assert all(r['complete'] and r['returncode'] == 0 for r in fast_rows)
        records_to_load = [(e, pp, ss) for e, pp, ss in records_to_load if e['method'] not in FAST_METHODS]
        fast_metadata = dict(campaign=str(rc), expected_sha256=sha(replacement_raw),
            source_sha256=rp['source_archive_sha256'], validated_complete_runs=528,
            source_changed_files=changed, retained_controls=[m for m in methods if m not in FAST_METHODS])
    replacement_expected_sha256 = None
    if args.pb_campaign and PB_REPLACEMENT in methods:
        pc = args.pb_campaign
        replacement_raw = (pc / 'expected.json').read_bytes()
        replacement_expected_sha256 = sha(replacement_raw)
        replacement_expected = json.loads(replacement_raw)
        replacement_provenance = json.loads((pc / 'source-provenance.json').read_text())
        replacement_scripts = json.loads((pc / 'submission.json').read_text())['script_sha256']
        assert replacement_provenance['source_archive_sha256'] == provenance['source_archive_sha256']
        assert {e['case'] for e in replacement_expected} == {
            e['case'] for e in expected if e['method'] == 'violation_pb_110' and e['case'] not in args.exclude_case}
        assert len(replacement_expected) == 131
        assert all(e['method'] == PB_REPLACEMENT for e in replacement_expected)
        records_to_load += [(e, replacement_provenance, replacement_scripts)
                            for e in replacement_expected]
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        rows = list(pool.map(lambda item: load_record(*item), records_to_load))
    if fast_rows is not None:
        rows += [r for r in fast_rows if r['method'] in methods and r['case'] not in args.exclude_case]
    scip_metadata = None
    if args.scip_campaign:
        sc = args.scip_campaign
        audit = json.loads(args.scip_audit.read_text())
        sp = json.loads((sc / 'source-provenance.json').read_text())
        se_raw = (sc / 'expected.json').read_bytes()
        scip_method = audit.get('method', 'scip_hybrid_all_005')
        assert scip_method in {'scip_hybrid_all_005', 'scip_hybrid_native_full'}
        se = [e for e in json.loads(se_raw) if e['method'] == scip_method]
        ss = {name: sha((sc / name).read_bytes()) for name in ('worker.jl', 'launch.py')}
        assert audit['passed'] and Path(audit['campaign']) == sc
        assert audit['expected_sha256'] == sha(se_raw)
        assert audit['source_sha256'] == sp['source_archive_sha256'] == sha((sc / 'source.tar').read_bytes())
        assert ss['worker.jl'] == scripts['worker.jl']
        if scip_method == 'scip_hybrid_native_full':
            assert ss == audit['script_sha256']
        else:
            assert ss['launch.py'] == scripts['launch.py']
        old_cases = {r['case'] for r in rows if r['method'] == REFERENCE}
        assert len(se) == len({e['case'] for e in se}) == len(old_cases)
        assert {e['case'] for e in se} == old_cases
        assert {e['method'] for e in se} == {scip_method}
        with (args.scip_audit.parent / 'records.csv').open() as stream:
            checked = {r['case']: r for r in csv.DictReader(stream)}
        eligibility = 'export_eligible' if scip_method == 'scip_hybrid_native_full' else 'validated_solved'
        assert set(checked) == old_cases and all(r[eligibility] == 'True' for r in checked.values())
        def load_scip(entry):
            row = load_record(entry, sp, ss)
            previous = checked[entry['case']]
            assert row['complete'] and row['returncode'] == 0
            assert row['manifest_sha256'] == previous['manifest_sha256']
            assert row['result_sha256'] == previous['result_sha256']
            output = Path(entry['output'])
            assert sha((output / 'resolved-params.json').read_bytes()) == previous['resolved_params_sha256']
            assert sha((output / 'worker.log').read_bytes()) == previous['worker_log_sha256']
            # Preserve the existing table/profile slot, while retaining the
            # actual method identity in every exported record and summary.
            row['source_method'] = row['method']
            row['method'] = REFERENCE
            return row
        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
            scip_rows = list(pool.map(load_scip, se))
        rows = [r for r in rows if r['method'] != REFERENCE] + scip_rows
        scip_metadata = dict(campaign=str(sc), method=scip_method,
            selector=('SCIP 10.0.1 hybrid with default cut limits and post-selection efficacy acceptance'
                      if scip_method == 'scip_hybrid_native_full'
                      else 'SCIP 10.0.1 hybrid with post-selection efficacy acceptance'),
            source_sha256=sp['source_archive_sha256'], expected_sha256=sha(se_raw),
            audit=str(args.scip_audit), audit_sha256=sha(args.scip_audit.read_bytes()),
            validated_runs=len(scip_rows), legacy_comparison_slot=REFERENCE)
    for row in rows:
        row.setdefault('source_method', row['method'])
    cases = sorted({r['case'] for r in rows})
    assert len(cases) == 132 - len(set(args.exclude_case))
    assert len(rows) == len(cases) * len(methods)
    by_case = {case: [r for r in rows if r['case'] == case] for case in cases}
    for case, records in by_case.items():
        opts = [r['objective'] for r in records if r['status'] == 'OPTIMAL'
                and r['primal_status'] == 'FEASIBLE_POINT']
        assert len(opts) > len(records) / 2, (case, opts)
        center = statistics.median(opts)
        tol = max(1e-5, abs(center) * .01)
        assert sum(abs(x - center) <= tol for x in opts) > len(records) / 2
        for r in records:
            r['consensus_objective'] = center
            r['objective_outlier'] = (r['status'] == 'OPTIMAL'
                                      and abs(r['objective'] - center) > tol)
            r['inverted_bound'] = (r['objective'] is not None and r['bound'] is not None
                                  and r['bound'] > r['objective'] + 1e-5 + 1e-6 * abs(r['objective']))
            r['screened_solved'] = (r['status'] == 'OPTIMAL' and not r['objective_outlier']
                                     and not r['inverted_bound']
                                     and r['primal_status'] == 'FEASIBLE_POINT')
            if r['screened_solved']:
                assert math.isfinite(r['seconds']) and r['seconds'] > 0
            r['penalized_seconds'] = r['seconds'] if r['screened_solved'] else LIMIT
    common = [case for case, rr in by_case.items() if all(r['screened_solved'] for r in rr)]
    main_common = [case for case, rr in by_case.items()
                   if all(r['screened_solved'] for r in rr if r['method'] in main_methods)]
    main_common_summary = []
    if args.common_only:
        assert main_common, 'No instances are solved by every main method'
        for method, label in zip(methods, labels):
            if method not in main_methods:
                continue
            rr = [r for r in rows if r['method'] == method and r['case'] in main_common]
            assert len(rr) == len(main_common) and all(r['screened_solved'] for r in rr)
            assert all(r['callback_rounds'] is not None and r['added_benders_cuts'] is not None for r in rr)
            main_common_summary.append(dict(method=method, label=label, instances=len(rr),
                mean_seconds=statistics.mean(r['seconds'] for r in rr),
                sgm_seconds=sgm([r['seconds'] for r in rr]),
                mean_callback_rounds=statistics.mean(r['callback_rounds'] for r in rr),
                mean_added_benders_cuts=statistics.mean(r['added_benders_cuts'] for r in rr)))
    common_lookup = {s['method']: s for s in main_common_summary}
    summary = []
    ratios = {m: [] for m in main_methods}
    appendix_ratios = {m: [] for m in methods}
    for case, rr in by_case.items():
        best = min(r['seconds'] for r in rr if r['screened_solved'])
        main_best = min((r['seconds'] for r in rr if r['screened_solved']
                         and r['method'] in main_methods), default=math.inf)
        for r in rr:
            appendix_ratios[r['method']].append(r['seconds'] / best if r['screened_solved'] else math.inf)
            if r['method'] in main_methods:
                ratios[r['method']].append(r['seconds'] / main_best if r['screened_solved'] else math.inf)
    for method, label in zip(methods, labels):
        rr = [r for r in rows if r['method'] == method]
        summary.append(dict(method=method, label=label,
            solved=sum(r['screened_solved'] for r in rr),
            time_limit=sum(r['status'] == 'TIME_LIMIT' for r in rr),
            cut_error=sum(r['status'] == 'CUT_ERROR' for r in rr),
            objective_outlier=sum(r['objective_outlier'] for r in rr),
            unfinished=sum(not r['complete'] for r in rr),
            incumbent_warnings=sum(r['incumbent_warning'] for r in rr),
            kmedoids_fallbacks=sum(r['kmedoids_fallback'] for r in rr),
            sgm_all=sgm([r['penalized_seconds'] for r in rr]),
            sgm_common=sgm([r['seconds'] for r in rr if r['case'] in common]),
            mean_all=statistics.mean(r['penalized_seconds'] for r in rr),
            mean_callback_rounds=statistics.mean(r['callback_rounds'] for r in rr if r['callback_rounds'] is not None),
            mean_added_benders_cuts=statistics.mean(r['added_benders_cuts'] for r in rr if r['added_benders_cuts'] is not None),
            mean_common=statistics.mean(r['seconds'] for r in rr if r['case'] in common),
            profile_scope='main' if method in main_methods else 'appendix_all_methods',
            fastest=sum(x == 1 for x in ratios.get(method, appendix_ratios[method])),
            within_twice=sum(x <= 2 for x in ratios.get(method, appendix_ratios[method])),
            appendix_fastest=sum(x == 1 for x in appendix_ratios[method]),
            appendix_within_twice=sum(x <= 2 for x in appendix_ratios[method])))
    ref = next(s for s in summary if s['method'] == REFERENCE)
    for s in summary:
        s['ratio_all'] = s['sgm_all'] / ref['sgm_all']
        s['ratio_common'] = s['sgm_common'] / ref['sgm_common']
    metadata = dict(snapshot_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        campaign=str(c), expected_sha256=sha(expected_raw), source_sha256=provenance['source_archive_sha256'],
        pb_campaign=str(args.pb_campaign) if args.pb_campaign else None,
        pb_expected_sha256=replacement_expected_sha256,
        replacement_campaign=fast_metadata,
        scip_replacement=scip_metadata,
        reference=REFERENCE, total_instances=len(cases), common_solved=len(common),
        excluded_methods=args.exclude_method, excluded_cases=args.exclude_case,
        main_methods=main_methods, main_excluded_methods=args.main_exclude_method,
        all_only=args.all_only, common_only=args.common_only, with_counts=args.with_counts,
        main_common_cases=main_common, main_common_summary=main_common_summary,
        table_count_policy=('identical instances solved by all main methods; no runtime cutoff'
                            if args.common_only else 'counts shown only for methods solving every included instance; otherwise omitted'),
        profile_denominator='fastest screened-solved runtime among main_methods on each instance',
        appendix_profile_denominator='fastest screened-solved runtime among all included methods on each instance',
        common_runtime_instances=len(common), common_cases=common, shift=10, penalty=LIMIT,
        objective_screen='within 1% of median reported optimal objective; strict majority required; not a correctness certificate',
        unfinished_policy='10800 s in provisional all-instance SGM, infinity in profile',
        common_policy=f'screened solved by all {len(methods)} included methods; no runtime cutoff',
        methods=summary)
    (out / 'summary.json').write_text(json.dumps(metadata, indent=2) + '\n')
    for filename, data in [('records.csv', rows), ('summary.csv', summary)]:
        with (out / filename).open('w') as f:
            w = csv.DictWriter(f, fieldnames=list(data[0])); w.writeheader(); w.writerows(data)
    compact_methods = (['random_all_005'] if 'random_all_005' in methods else []) + [
        'violation_all_005', 'diversity_all_005', 'hybrid_all_005', 'hybrid_plus_all_005']
    compact_methods = [m for m in compact_methods if m not in args.main_exclude_method]
    if args.common_only:
        compact_methods = main_methods
    if args.replacement_campaign and all(m in methods for m in compact_methods):
        compact_rows = main_common_summary if args.common_only else [next(s for s in summary if s['method'] == m) for m in compact_methods]
        if not args.common_only:
            assert all(s['solved'] == len(cases) for s in compact_rows), 'Compact arithmetic-means table requires all included fixed-budget methods solved on the full cohort'
        compact = [r'% Generated from the same records and count scope as the main runtime table.',
            r'\begin{tabular}{lrr}', r'\toprule',
            r'Strategy & Mean callback rounds & Mean added Benders cuts \\', r'\midrule']
        for index, s in enumerate(compact_rows):
            compact.append(f"{s['label']} & {s['mean_callback_rounds']:.2f} & {s['mean_added_benders_cuts']:.2f}" + r' \\')
            if s['label'] == 'Random' and index < len(compact_rows) - 1:
                compact.append(r'\addlinespace[4pt]')
        compact += [r'\bottomrule', r'\end{tabular}']
        (paper / 'callback_cut_table.tex').write_text('\n'.join(compact) + '\n')
        count_cases = main_common if args.common_only else cases
        lookup = {m: {r['case']: r for r in rows if r['method'] == m and r['case'] in count_cases} for m in compact_methods}
        diagnostics = dict(instances=len(count_cases), compact_methods=compact_methods,
            scope=('All main methods on identical jointly solved instances' if args.common_only else
                   'All solved fixed-budget methods listed in compact_methods; main exclusions applied; SCIP and PB excluded'),
            arithmetic_means=compact_rows,
            versus_violation={m: paired_metrics(lookup['violation_all_005'], lookup[m])
                              for m in ('hybrid_all_005', 'hybrid_plus_all_005') if m in compact_methods})
        if fast_rows is not None:
            old_expected = [e for e in base_expected if e['method'] in FAST_METHODS]
            with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
                old_rows = list(pool.map(lambda e: load_record(e, provenance, scripts), old_expected))
            diagnostics['old_vs_new'] = {}
            for m in sorted(FAST_METHODS & set(methods)):
                old = {r['case']: r for r in old_rows if r['method'] == m}
                new = {r['case']: r for r in rows if r['method'] == m}
                diagnostics['old_vs_new'][m] = paired_metrics(old, new)
            diagnostics['old_vs_new_raw'] = [dict(method=r['method'], case=r['case'],
                old_seconds=r['seconds'], old_callback_rounds=r['callback_rounds'], old_added_benders_cuts=r['added_benders_cuts'],
                new_seconds=next(n['seconds'] for n in rows if n['method']==r['method'] and n['case']==r['case']),
                new_callback_rounds=next(n['callback_rounds'] for n in rows if n['method']==r['method'] and n['case']==r['case']),
                new_added_benders_cuts=next(n['added_benders_cuts'] for n in rows if n['method']==r['method'] and n['case']==r['case'])) for r in old_rows]
        (out / 'callback-cut-diagnostics.json').write_text(json.dumps(diagnostics, indent=2) + '\n')
    bold_best = args.with_counts and (args.all_only or args.common_only)
    if bold_best:
        displayed = [s for s in summary if s['method'] in main_methods]
        displayed_stats = lambda s: common_lookup[s['method']] if args.common_only else s
        count_rows = [s for s in displayed if args.common_only or s['solved'] == len(cases)]
        best_mean = min((displayed_stats(s)['mean_seconds'] if args.common_only else s['mean_all']) for s in displayed)
        best_sgm = min((displayed_stats(s)['sgm_seconds'] if args.common_only else s['sgm_all']) for s in displayed)
        best_callbacks = min(displayed_stats(s)['mean_callback_rounds'] for s in count_rows)
        best_cuts = min(displayed_stats(s)['mean_added_benders_cuts'] for s in count_rows)

    def table_number(value, best, places=2):
        number = f'{value:.{places}f}'
        return r'\bfseries ' + number if bold_best and number == f'{best:.{places}f}' else number

    lines = [r'% Generated by check/scripts/export_full_cutfilter.py.',
             r'\begin{tabular}{lrrrrrrr}', r'\toprule',
             r' & & \multicolumn{3}{c}{All %d instances} & \multicolumn{3}{c}{Common, $n=%d$} \\' % (len(cases), len(common)),
             r'\cmidrule(lr){3-5}\cmidrule(lr){6-8}',
             r'Strategy & Solved & Mean (\unit{\second}) & SGM (\unit{\second}) & Ratio & Mean (\unit{\second}) & SGM (\unit{\second}) & Ratio \\', r'\midrule']
    if args.all_only or args.common_only:
        lines = [r'% Generated by check/scripts/export_full_cutfilter.py.',
                 r'\begin{tabular}{lrrr}', r'\toprule',
                 r'Strategy & Solved & Mean (\unit{\second}) & SGM (\unit{\second}) \\', r'\midrule']
        if args.with_counts:
            lines = [r'% Generated by check/scripts/export_full_cutfilter.py.',
                     r'\begingroup',
                     r'\setlength{\tabcolsep}{4pt}',
                     r'\renewcommand{\arraystretch}{1.15}',
                     r'\sisetup{group-separator={\,},group-minimum-digits=4,detect-weight=true}',
                     r'\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}l S[table-format=3.0] *{3}{S[table-format=3.2]} S[table-format=5.2]@{}}', r'\toprule',
                     r' & & \multicolumn{2}{c}{Solve time (\unit{\second})} & \multicolumn{2}{c}{Mean counts} \\',
                     r'\cmidrule(lr){3-4}\cmidrule(lr){5-6}',
                     r'Strategy & {Solved} & {Mean} & {SGM} & {\shortstack{Callback\\rounds}} & {\shortstack{Benders\\cuts}} \\', r'\midrule']
    for index, s in enumerate(summary):
        if s['method'] not in main_methods:
            continue
        if args.all_only or args.common_only:
            stats = common_lookup[s['method']] if args.common_only else s
            mean_time = stats['mean_seconds'] if args.common_only else s['mean_all']
            sgm_time = stats['sgm_seconds'] if args.common_only else s['sgm_all']
            row = (f"{s['label']} & {s['solved']}"
                   f" & {table_number(mean_time, best_mean) if bold_best else f'{mean_time:.2f}'}"
                   f" & {table_number(sgm_time, best_sgm) if bold_best else f'{sgm_time:.2f}'}")
            if args.with_counts:
                if args.common_only or s['solved'] == len(cases):
                    count_cases = main_common if args.common_only else cases
                    rr = [r for r in rows if r['method'] == s['method'] and r['case'] in count_cases]
                    assert len(rr) == len(count_cases) and all(
                        r['callback_rounds'] is not None and r['added_benders_cuts'] is not None for r in rr
                    ), 'Displayed count means require complete counts on every instance'
                    row += (f" & {table_number(stats['mean_callback_rounds'], best_callbacks)}"
                            f" & {table_number(stats['mean_added_benders_cuts'], best_cuts)}")
                else:
                    row += ' & {--} & {--}'
            lines.append(row + r' \\')
        else:
            lines.append(f"{s['label']} & {s['solved']} & {s['mean_all']:.2f} & {s['sgm_all']:.2f} & {s['ratio_all']:.3f} & {s['mean_common']:.2f} & {s['sgm_common']:.2f} & {s['ratio_common']:.3f}" + r' \\')
        if s['label'] == 'Random' and any(t['method'] in main_methods for t in summary[index + 1:]):
            lines.append(r'\addlinespace[4pt]')
    if args.with_counts and (args.all_only or args.common_only):
        lines += [r'\bottomrule', r'\end{tabular*}', r'\endgroup']
    else:
        lines += [r'\bottomrule', r'\end{tabular}']
    (paper / 'runtime_table.tex').write_text('\n'.join(lines) + '\n')
    # The appendix uses the same records and classification, avoiding stale
    # per-instance data beside the refreshed summary.
    case_order = sorted(cases, key=lambda x: tuple(int(v) for v in x.split('case')[1].split('_')))
    # Align decimal markers directly when normal and bold numbers share a column.
    for part in range(3):
        lines = [r'\begin{tabular}{l*{' + str(len(methods)) + r'}{S}}', r'\toprule',
                 'Instance & ' + ' & '.join('{' + ('SCIP' if label == 'SCIP cut selection' else label) + '}' for label in labels) + r' \\',
                 r'\midrule']
        for case in case_order[part * 44:(part + 1) * 44]:
            rr = {r['method']: r for r in by_case[case]}
            best = min(r['seconds'] for r in rr.values() if r['screened_solved'])
            cells = [instance_label(case)]
            for method in methods:
                r = rr[method]
                if r['screened_solved']:
                    cell = f"{r['seconds']:.2f}"
                    if r['seconds'] == best:
                        cell = r'\bfseries ' + cell
                elif r['objective_outlier']:
                    cell = '{OBJ}'
                else:
                    cell = '{' + {'TIME_LIMIT': 'TL', 'CUT_ERROR': 'CE', 'RUNNING': 'RUN'}.get(r['status'], 'ERR') + '}'
                cells.append(cell)
            lines.append(' & '.join(cells) + r' \\')
        lines += [r'\bottomrule', r'\end{tabular}']
        (paper / f'per_instance_part{part + 1}.tex').write_text('\n'.join(lines) + '\n')
    xmax = max(4, 2 ** math.ceil(math.log2(max(x for rr in ratios.values() for x in rr if math.isfinite(x)))))
    ticks = [2 ** k for k in range(int(math.log2(xmax)) + 1)]
    lines = [r'% Generated by check/scripts/export_full_cutfilter.py.', r'\begin{tikzpicture}']
    for i, (method, col) in enumerate(zip(methods, colors)):
        if method in main_methods:
            lines.append(r'\definecolor{fullprofile%d}{HTML}{%s}' % (i, col))
    lines += [r'\begin{axis}[width=.90\linewidth,height=.49\linewidth,',
        f'xmin=1,xmax={xmax},ymin=0,ymax=1.025,',
        r'xmode=log,log basis x=2,xtick={' + ','.join(map(str, ticks)) + '},',
        r'xticklabels={' + ','.join(map(str, ticks)) + '},',
        r'ytick={0,.2,.4,.6,.8,1},yticklabels={' + ','.join(rf'\qty{{{value}}}{{\percent}}' for value in range(0, 101, 20)) + '},',
        r'xlabel={Runtime ratio to fastest method, $\tau$},ylabel={Fraction of instances},',
        r'grid=major,grid style={black!15},tick align=outside,tick pos=left,',
        r'tick label style={font=\small},label style={font=\small},',
        r'legend columns=1,legend cell align=left,',
        r'legend style={at={(.98,.04)},anchor=south east,draw=none,',
        r'fill=white,fill opacity=.95,text opacity=1,font=\small,inner sep=2pt},',
        r'every axis plot/.append style={no marks,line width=1.2pt}]']
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'DejaVu Serif', 'pdf.fonttype': 42})
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for i, (method, label) in enumerate(zip(methods, labels)):
        if method not in main_methods:
            continue
        points = step_points(ratios[method], xmax, len(cases))
        line_style = PLOT_STYLES.get(method, 'solid')
        plot_options = f'fullprofile{i}' + (f',{line_style}' if line_style != 'solid' else '')
        lines.append(r'\addplot[%s] coordinates {' % plot_options)
        lines.extend(f'({x:.15g},{y:.15g})' for x, y in points)
        lines += ['};', r'\addlegendentry{%s}' % label]
        ax.plot([p[0] for p in points], [p[1] for p in points], color='#' + colors[i],
                linestyle={'solid': '-', 'dashed': '--', 'dotted': ':'}[line_style],
                label=label, linewidth=1.8)
    lines += [r'\end{axis}', r'\end{tikzpicture}']
    (paper / 'performance_profile.tikz').write_text('\n'.join(lines) + '\n')
    ax.set(xscale='log', xlim=(1, xmax), ylim=(0, 1.025),
           xlabel='Runtime ratio to fastest method', ylabel='Fraction of instances')
    ax.set_xticks(ticks, [str(x) for x in ticks]); ax.grid(alpha=.2)
    ax.legend(loc='upper center', bbox_to_anchor=(.5, -.18), ncol=3, frameon=False)
    fig.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(out / ('performance_profile.' + ext), dpi=160, bbox_inches='tight')
    print(json.dumps(metadata, indent=2))


if __name__ == '__main__':
    main()
