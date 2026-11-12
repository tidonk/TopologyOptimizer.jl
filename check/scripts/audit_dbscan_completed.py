#!/usr/bin/env python3
"""Read-only campaign audit; write summaries separately from production runs."""
import collections
import concurrent.futures
import csv
import hashlib
import importlib.util
import json
import math
import re
import statistics
from pathlib import Path
import os

ROOT=(Path(os.environ.get("TOPOLOGYOPTIMIZER_RESULTS_ROOT", Path(__file__).resolve().parents[2] / "check/results")))
C=ROOT/'20261026-cutfilter-dbscan-eps050-min5-b8a05381'
B=ROOT/'20261023-cutfilter-full-005-fast-b8a05381'
OUT=C/'analysis-20261026'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')

def main():
    spec=importlib.util.spec_from_file_location('export',ROOT.parent/'scripts/export_full_cutfilter.py')
    export=importlib.util.module_from_spec(spec);spec.loader.exec_module(export)
    expected=read(C/'expected.json');p=read(C/'source-provenance.json');bp=read(B/'source-provenance.json')
    assert len(expected)==len({e['case'] for e in expected})==131
    assert sha(C/'source.tar')==p['source_archive_sha256']==bp['source_archive_sha256']
    scripts={n:sha(C/n) for n in ('launch.py','worker.jl')}
    assert all(scripts[n]==sha(B/n) for n in scripts)
    reference=list(csv.DictReader((B/'analysis-paper-20261026-main-six-common/records.csv').open()))
    lookup={(r['case'],r['method']):r for r in reference}
    size={r['case']:int(r['buses']) for r in csv.DictReader((B/'analysis-size-20261026/per_instance_size_effects.csv').open())}
    common=set(read(B/'analysis-paper-20261026-main-six-common/summary.json')['main_common_cases'])
    reference_selection=read(B/'results/hybrid_all_005/pglib_scopf_case100/resolved-params.json')['BENDERS']['CUTS']['SELECTION']
    def load(e):
        assert e['instance_sha256']==bp['input_sha256'][e['case']+'.m']==sha(C/e['instance'])
        assert sha(C/e['settings'])==p['settings_sha256'][e['method']]
        row=export.load_record(e,p,scripts);assert row['complete'],e['case']
        result=read(Path(e['output'])/'result.json');resolved=read(Path(e['output'])/'resolved-params.json')
        sel=resolved['BENDERS']['CUTS']['SELECTION']
        assert sel['DIVERSITY']['USE_DBSCAN'] is True and sel['MIN_OPT_CUTS']==0
        assert sel['METHOD']=='DIVERSITY' and sel['FEASCUT_INF_VIOL'] is True
        assert sel['DIVERSITY']['REPRESENTATIVE_MODE']=='max_violation'
        assert sel['DIVERSITY']['DBSCAN_EPS']==.5 and sel['DIVERSITY']['DBSCAN_MIN_PTS']==5
        sel['DIVERSITY']['USE_DBSCAN']=False
        assert sel==reference_selection,'Selection settings differ beyond USE_DBSCAN'
        ref=float(lookup[e['case'],'hybrid_all_005']['consensus_objective'])
        row['objective_absolute_difference']=abs(row['objective']-ref)
        row['objective_relative_difference']=row['objective_absolute_difference']/max(abs(ref),1e-12)
        row['strict_objective_agreement']=row['objective_absolute_difference']<=1e-5+1e-6*abs(ref)
        row['objective_outlier_1pct']=row['objective_absolute_difference']>max(1e-5,.01*abs(ref))
        row['inverted_bound']=row['bound']>row['objective']+1e-5+1e-6*abs(row['objective'])
        log=(Path(e['output'])/'worker.log').read_text()
        tolerances=re.findall(r'Optimal solution found \(tolerance ([0-9.eE+-]+)\)',log)
        assert tolerances,e['case']
        row['reported_mip_relative_tolerance']=float(tolerances[-1])
        row['objective_agrees_within_mip_tolerance']=row['objective_absolute_difference']<=1e-5+row['reported_mip_relative_tolerance']*max(abs(ref),abs(row['objective']))
        other=[lookup[e['case'],m] for m in ('violation_all_005','hybrid_all_005')]
        lb=max([row['bound']]+[float(r['bound']) for r in other])
        ub=min([row['objective']]+[float(r['objective']) for r in other])
        row['paired_bound_intervals_overlap']=lb<=ub+1e-5+1e-6*abs(ub)
        for k in ('n_benders_cut_value_disagree','n_benders_cut_value_fail'):
            row[k]=result['callback'][k]
        row['buses']=size[e['case']]
        row['screened_solved']=row['status']=='OPTIMAL' and row['primal_status']=='FEASIBLE_POINT' and not row['objective_outlier_1pct'] and not row['inverted_bound']
        return row
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:rows=list(pool.map(load,expected))
    assert all(r['screened_solved'] for r in rows),'Review incomplete/invalid results before comparing'
    summaries=[]
    groups=[('All 131',lambda r:True),('Common 126',lambda r:r['case'] in common),
            ('Small <200',lambda r:r['buses']<200),('Medium 200–299',lambda r:200<=r['buses']<300),('Large >=300',lambda r:r['buses']>=300)]
    for group,test in groups:
        rr=[r for r in rows if test(r)]
        for method,label in [('violation_all_005','Violation'),('hybrid_all_005','Hybrid'),('dbscan','DBSCAN')]:
            def values(key):return [r[key] if method=='dbscan' else float(lookup[r['case'],method][key]) for r in rr]
            t=values('seconds')
            summaries.append(dict(group=group,method=label,n=len(rr),mean_seconds=statistics.mean(t),sgm_seconds=export.sgm(t),
                mean_callbacks=statistics.mean(values('callback_rounds')),mean_cuts=statistics.mean(values('added_benders_cuts'))))
    diagnostics=dict(completed=len(rows),status_counts=dict(collections.Counter(r['status'] for r in rows)),
        strict_objective_agreements=sum(r['strict_objective_agreement'] for r in rows),
        objective_agreements_within_reported_mip_tolerance=sum(r['objective_agrees_within_mip_tolerance'] for r in rows),
        overlapping_bound_intervals_with_violation_and_hybrid=sum(r['paired_bound_intervals_overlap'] for r in rows),
        max_absolute_objective_difference=max(r['objective_absolute_difference'] for r in rows),
        max_relative_objective_difference=max(r['objective_relative_difference'] for r in rows),
        objective_outliers=sum(r['objective_outlier_1pct'] for r in rows),inverted_bounds=sum(r['inverted_bound'] for r in rows),
        cut_value_disagreements=sum(r['n_benders_cut_value_disagree'] for r in rows),cut_value_failures=sum(r['n_benders_cut_value_fail'] for r in rows),
        incumbent_warnings=sum(r['incumbent_warning'] for r in rows),clustering_fallbacks=sum(r['kmedoids_fallback'] for r in rows),
        faster_than_hybrid=sum(r['seconds']<float(lookup[r['case'],'hybrid_all_005']['seconds']) for r in rows),
        faster_than_violation=sum(r['seconds']<float(lookup[r['case'],'violation_all_005']['seconds']) for r in rows))
    OUT.mkdir(exist_ok=True)
    save(OUT/'audit.json',dict(campaign=str(C),source_sha256=p['source_archive_sha256'],settings=read(C/'settings-audit.json'),
        diagnostics=diagnostics,summaries=summaries,scope='Existing generated benchmark, case245 excluded; no independent proof of optimality beyond recorded solver/cross-method checks'))
    for name,rr in [('records.csv',rows),('summary.csv',summaries)]:
        with (OUT/name).open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
    print(json.dumps(dict(diagnostics=diagnostics,summaries=summaries),indent=2))

if __name__=='__main__':main()
