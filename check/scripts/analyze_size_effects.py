#!/usr/bin/env python3
"""Exploratory size/runtime analysis of the exact six-method paper records.

Does not edit manuscript assets or mix in the pending SCIP hybrid rerun.
"""
import argparse
import concurrent.futures
import csv
import hashlib
import json
import re
from pathlib import Path
import os
import numpy as np
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter, NullFormatter

BASE = (Path(os.environ.get("TOPOLOGYOPTIMIZER_RESULTS_ROOT", Path(__file__).resolve().parents[2] / "check/results")) / '20261023-cutfilter-full-005-fast-b8a05381')
V, H = 'violation_all_005', 'hybrid_all_005'
LABELS = {'unfiltered':'No filtering', 'scip_dynamic_all_005':'SCIP dynamic',
          'random_all_005':'Random', V:'Violation', 'diversity_all_005':'Diversity', H:'Hybrid'}
ORANGE, PURPLE = '#D55E00', '#785EF0'

def sgm(t, shift=10):
    return float(np.expm1(np.mean(np.log(np.asarray(t)+shift)))+1-shift)

def write_csv(path, rows):
    with path.open('w') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def size_metadata(case, provenance):
    raw=(BASE/'instances'/f'{case}.m').read_bytes()
    assert hashlib.sha256(raw).hexdigest()==provenance['input_sha256'][f'{case}.m']
    text=raw.decode()
    def count(name):
        block=re.search(r'mpc\.'+name+r'\s*=\s*\[(.*?)\];',text,re.S).group(1)
        return sum(bool(line.split('%')[0].strip()) for line in block.splitlines())
    # Use the actual solver-reported contingency count, not a branch-count proxy.
    log=(BASE/'results'/H/case/'worker.log').read_text()
    matches=re.findall(r'Contingencies:\s*(\d+)',log)
    assert matches and len(set(matches))==1,case
    return case,dict(buses=count('bus'),branches=count('branch'),scenarios=int(matches[0]))

def paired_stats(cases, lookup, sizes):
    v=np.array([lookup[c,V]['seconds'] for c in cases]); h=np.array([lookup[c,H]['seconds'] for c in cases])
    result=dict(n=len(cases),violation_mean=float(v.mean()),hybrid_mean=float(h.mean()),
        mean_reduction_pct=float(100*(1-h.mean()/v.mean())),violation_sgm=sgm(v),hybrid_sgm=sgm(h),
        sgm_reduction_pct=float(100*(1-sgm(h)/sgm(v))),hybrid_faster=int(sum(h<v)),
        median_shifted_speedup=float(np.median((v+10)/(h+10))),
        saved_seconds=float(sum(v-h)),bus_min=min(sizes[c]['buses'] for c in cases),
        bus_max=max(sizes[c]['buses'] for c in cases))
    for metric in ('callback_rounds','added_benders_cuts'):
        before=np.mean([lookup[c,V][metric] for c in cases]);after=np.mean([lookup[c,H][metric] for c in cases])
        result[metric+'_reduction_pct']=float(100*(1-after/before))
    return result

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args(); out=args.output;out.mkdir(parents=True,exist_ok=True)
    source=BASE/'analysis-paper-20261026-main-six-common'
    summary=json.loads((source/'summary.json').read_text())
    methods=summary['main_methods'];common=summary['main_common_cases']
    assert len(common)==126 and summary['shift']==10
    raw=(source/'records.csv').read_bytes()
    rows=list(csv.DictReader(raw.decode().splitlines()))
    lookup={}
    for r in rows:
        for k in ['seconds','penalized_seconds','callback_rounds','added_benders_cuts']:r[k]=float(r[k])
        r['screened_solved']=r['screened_solved']=='True';lookup[r['case'],r['method']]=r
    cases=sorted({r['case'] for r in rows});assert len(cases)==131
    for s in summary['main_common_summary']:
        rr=[lookup[c,s['method']]['seconds'] for c in common]
        assert abs(np.mean(rr)-s['mean_seconds'])<1e-9 and abs(sgm(rr)-s['sgm_seconds'])<1e-9
    provenance=json.loads((BASE/'source-provenance.json').read_text())
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        sizes=dict(pool.map(lambda c:size_metadata(c,provenance),cases))
    assert all(sizes[c]['buses']==int(c.split('case')[1].split('_')[0]) for c in cases)
    groups=[('Small (<200 buses)',lambda n:n<200),('Medium (200–299)',lambda n:200<=n<300),('Large (≥300)',lambda n:n>=300)]
    common_groups=[(label,[c for c in common if test(sizes[c]['buses'])]) for label,test in groups]
    all_groups=[(label,[c for c in cases if test(sizes[c]['buses'])]) for label,test in groups]
    # Alternate, size-only boundaries: empirical terciles of all 131 sizes.
    cuts=np.quantile([sizes[c]['buses'] for c in cases],[1/3,2/3])
    terciles=[(f'Tercile {i+1}',[c for c in common if np.searchsorted(cuts,sizes[c]['buses'],side='right')==i]) for i in range(3)]
    pairs=[dict(group=label,**paired_stats(cc,lookup,sizes)) for label,cc in common_groups]
    tables=[]
    for label,cc in [('All jointly solved',common)]+common_groups:
        for m in methods:
            t=[lookup[c,m]['seconds'] for c in cc]
            tables.append(dict(group=label,n=len(cc),method=LABELS[m],mean_seconds=float(np.mean(t)),sgm_seconds=sgm(t)))
    v=np.array([lookup[c,V]['seconds'] for c in common]);h=np.array([lookup[c,H]['seconds'] for c in common])
    delta=v-h; logr=np.log((v+10)/(h+10));rng=np.random.default_rng(20261026)
    corr={}
    for metric in ('buses','branches','scenarios'):
        x=np.array([sizes[c][metric] for c in common]);rho,p=spearmanr(x,logr)
        boots=[spearmanr(x[i],logr[i])[0] for i in rng.integers(0,len(common),(2000,len(common)))]
        corr[metric]=dict(spearman=float(rho),p_unadjusted=float(p),bootstrap_95=np.quantile(boots,[.025,.975]).tolist())
    order=np.argsort(-delta)
    contributions=[dict(case=common[i],**sizes[common[i]],violation_seconds=float(v[i]),hybrid_seconds=float(h[i]),
                        saved_seconds=float(delta[i]),share_net_savings_pct=float(100*delta[i]/delta.sum())) for i in order]
    robustness={}
    for k in (0,1,3,5,10):
        keep=order[k:]
        robustness[f'omit_{k}_largest_savings']=dict(n=len(keep),violation_mean=float(v[keep].mean()),
            hybrid_mean=float(h[keep].mean()),mean_reduction_pct=float(100*(1-h[keep].mean()/v[keep].mean())),
            violation_sgm=sgm(v[keep]),hybrid_sgm=sgm(h[keep]))
    effects=[]
    for c in cases:
        r=dict(case=c,**sizes[c],jointly_solved=c in common)
        for m in methods:
            r[LABELS[m]+'_seconds']=lookup[c,m]['seconds'];r[LABELS[m]+'_solved']=lookup[c,m]['screened_solved']
        r['hybrid_shifted_speedup']=(lookup[c,V]['seconds']+10)/(lookup[c,H]['seconds']+10)
        r['hybrid_saved_seconds']=lookup[c,V]['seconds']-lookup[c,H]['seconds']
        effects.append(r)
    # Sanity check all 131: V/H solve every case, so this comparison needs no penalty.
    assert all(lookup[c,m]['screened_solved'] for c in cases for m in (V,H))
    all_pairs=[dict(group=label,**paired_stats(cc,lookup,sizes)) for label,cc in all_groups]
    all_method_sensitivity=[]
    for label,cc in all_groups:
        for m in methods:
            t=[lookup[c,m]['penalized_seconds'] for c in cc]
            all_method_sensitivity.append(dict(group=label,method=LABELS[m],n=len(cc),
                solved=sum(lookup[c,m]['screened_solved'] for c in cc),mean_seconds=float(np.mean(t)),sgm_seconds=sgm(t)))
    paired_boot=rng.integers(0,len(common),(10000,len(common)))
    mean_boot=100*(1-h[paired_boot].mean(axis=1)/v[paired_boot].mean(axis=1))
    log_boot=logr[paired_boot].mean(axis=1)
    report=dict(source=str(source),records_sha256=hashlib.sha256(raw).hexdigest(),
        purpose='Exploratory size stratification and explanation of arithmetic vs shifted geometric mean; no paper edits',
        shift=10,solver_limit=10800,main_instances=126,full_instances=131,
        old_scip_policy='Existing dynamic selector only; pending default-hybrid rerun is excluded',
        size_definition='Bus count from hashed input; branch count from input; scenarios from solver logs',
        primary_boundaries=[200,300],boundary_policy='Round bus-count thresholds; no runtime-based boundary optimization',
        global_stats=paired_stats(common,lookup,sizes),size_groups=pairs,correlations=corr,
        alternate_size_tercile_boundaries=cuts.tolist(),
        alternate_size_terciles=[dict(group=l,**paired_stats(cc,lookup,sizes)) for l,cc in terciles],
        all_131=paired_stats(cases,lookup,sizes),all_131_groups=all_pairs,
        largest_savings=contributions[:10],largest_losses=contributions[-5:],robustness=robustness,
        bootstrap_95_mean_reduction_pct=np.quantile(mean_boot,[.025,.975]).tolist(),
        bootstrap_95_geometric_shifted_speedup=np.exp(np.quantile(log_boot,[.025,.975])).tolist(),
        caveat='Single runs per instance; bootstrap concerns sampled test instances, not measurement noise. Correlations and all splits are exploratory, not causal or confirmatory.')
    for name,data in [('size_groups.csv',pairs),('all_method_size_table.csv',tables),('per_instance_size_effects.csv',effects),
                      ('savings_contributions.csv',contributions),('all_131_size_groups.csv',all_pairs),
                      ('all_131_penalized_methods.csv',all_method_sensitivity)]:write_csv(out/name,data)
    (out/'analysis.json').write_text(json.dumps(report,indent=2)+'\n')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,
                         'axes.spines.top':False,'axes.spines.right':False})
    figures=[]
    buses=np.array([sizes[c]['buses'] for c in common]);ratios=np.exp(logr)
    fig,ax=plt.subplots(1,2,figsize=(11.2,4.5),gridspec_kw={'width_ratios':[1,1.12]})
    lim=[min(v.min(),h.min())*.75,max(v.max(),h.max())*1.3]
    pts=ax[0].scatter(v,h,c=buses,cmap='viridis',s=30,alpha=.8,edgecolor='white',linewidth=.3,zorder=3)
    ax[0].plot(lim,lim,color='.4',lw=1);ax[0].set(xscale='log',yscale='log',xlim=lim,ylim=lim,
        xlabel='Violation solve time (s)',ylabel='Hybrid solve time (s)',title='(a) Paired runtimes')
    ax[0].text(.04,.95,'Above line: Violation faster',transform=ax[0].transAxes,va='top',color=ORANGE,fontsize=9)
    ax[0].text(.96,.04,'Below line: Hybrid faster',transform=ax[0].transAxes,ha='right',color=PURPLE,fontsize=9)
    for i in order[:3]:ax[0].annotate(common[i].split('case')[1],(v[i],h[i]),xytext=(-5,-13),
        textcoords='offset points',ha='right',fontsize=8)
    cb=fig.colorbar(pts,ax=ax[0],pad=.02,fraction=.05);cb.set_label('Buses')
    ax[1].scatter(buses,ratios,c=np.where(ratios>1,PURPLE,ORANGE),s=27,alpha=.65,edgecolor='white',linewidth=.3)
    ax[1].axhline(1,color='.35',lw=1);ax[1].set(yscale='log',ylim=(.5,2.7),xlabel='Number of buses',
        ylabel='Shifted speedup of Hybrid over Violation',title='(b) Size and relative runtime')
    ax[1].set_yticks([.5,1,2]);ax[1].yaxis.set_major_formatter(FuncFormatter(lambda t,_:f'{t:g}'))
    ax[1].yaxis.set_minor_formatter(NullFormatter())
    for boundary in [200,300]:ax[1].axvline(boundary,color='.8',linestyle=':',lw=.8)
    rho=corr['buses']['spearman'];lo,hi=corr['buses']['bootstrap_95']
    ax[1].text(.02,.98,f'Spearman ρ = {rho:.2f}',transform=ax[1].transAxes,va='top',fontsize=9)
    ax[1].text(.98,.02,'Above 1: Hybrid faster',transform=ax[1].transAxes,ha='right',fontsize=9,color=PURPLE)
    for a in ax:a.grid(alpha=.15,which='major');a.set_axisbelow(True)
    fig.suptitle('Hybrid versus Violation — 126 jointly solved instances',fontsize=12)
    fig.text(.5,.01,'Shifted speedup = (Violation time + 10 s) / (Hybrid time + 10 s).',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.05,1,1]);figures.append(('runtime_and_size',fig))
    # Distribution of relative effects beside transparent absolute contributions.
    fig,ax=plt.subplots(1,2,figsize=(11.2,4.3),gridspec_kw={'width_ratios':[1,1.25]})
    sort=np.sort(ratios);ax[0].step(sort,np.arange(1,len(sort)+1)/len(sort),where='post',color=PURPLE,lw=2)
    ax[0].axvline(1,color='.4',lw=1);ax[0].axhline(.5,color='.6',ls=':',lw=1)
    ax[0].set(xscale='log',xlim=(.5,2.7),ylim=(0,1.02),xlabel='Shifted speedup of Hybrid over Violation',
              ylabel='Fraction of instances',title='(a) Distribution of relative runtimes')
    ax[0].set_xticks([.5,1,2]);ax[0].xaxis.set_major_formatter(FuncFormatter(lambda t,_:f'{t:g}'))
    ax[0].xaxis.set_minor_formatter(NullFormatter())
    ax[0].text(.04,.94,f'Hybrid faster: {sum(h<v)}/{len(v)}\nMedian ratio: {np.median(ratios):.3f}',transform=ax[0].transAxes,va='top')
    chosen=order[:5];bars=[delta[i] for i in chosen]+[float(sum(delta[order[5:]]))]
    labels=[common[i].split('case')[1] for i in chosen]+[f'Other {len(common)-5}\ninstances']
    ax[1].barh(np.arange(6),bars,color=[PURPLE if x>0 else ORANGE for x in bars]);ax[1].invert_yaxis()
    ax[1].set_yticks(np.arange(6));ax[1].set_yticklabels(labels);ax[1].axvline(0,color='.4',lw=.8)
    ax[1].set(xlabel='Total seconds saved by Hybrid (Violation − Hybrid)',title='(b) Contributions to the mean advantage')
    for i,val in enumerate(bars):ax[1].annotate(f'{val:+,.0f}',(val,i),xytext=(4 if val>=0 else -4,0),
        textcoords='offset points',va='center',ha='left' if val>=0 else 'right',fontsize=9)
    ax[1].margins(x=.20)
    for a in ax:a.grid(alpha=.15,axis='x');a.set_axisbelow(True)
    fig.suptitle('Why the arithmetic and shifted geometric means differ',fontsize=12)
    fig.text(.5,.01,'126 jointly solved instances. Shifted speedup = (Violation time + 10 s) / (Hybrid time + 10 s).',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.05,1,1]);figures.append(('relative_effects_and_savings',fig))
    # Compact method table; every cell uses identical cases within a size group.
    fig,ax=plt.subplots(figsize=(11.2,4.4));ax.axis('off')
    cols=['Method']+[f'{label}\nn={len(cc)}' for label,cc in common_groups]+['All\nn=126']
    cells=[]
    for m in methods:
        vals=[]
        for label,cc in common_groups+[('All jointly solved',common)]:
            t=[lookup[c,m]['seconds'] for c in cc];vals.append(f'{np.mean(t):.1f} / {sgm(t):.1f}')
        cells.append([LABELS[m]]+vals)
    tab=ax.table(cellText=cells,colLabels=cols,loc='center',cellLoc='center',colWidths=[.21,.19,.19,.19,.19])
    tab.auto_set_font_size(False);tab.set_fontsize(10);tab.scale(1,2.0)
    for (r,c),cell in tab.get_celld().items():
        cell.set_edgecolor('#dddddd');cell.set_linewidth(.5)
        if r==0:cell.set_facecolor('#eeeeee');cell.set_text_props(weight='bold');cell.set_height(.16)
        elif r==6:cell.set_facecolor('#f1edfe')
    ax.set_title('Solve times by network size: arithmetic mean / shifted geometric mean (s)',pad=18)
    fig.text(.5,.04,'Same 126 jointly solved instances as the current paper; shift 10 s; three-hour limit.\nSCIP denotes the existing dynamic selector. The replacement rerun is not included.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.10,1,1]);figures.append(('size_table',fig))
    with PdfPages(out/'size_analysis_preview.pdf') as pdf:
        for name,fig in figures:
            fig.savefig(out/f'{name}.png',dpi=170,bbox_inches='tight')
            fig.savefig(out/f'{name}.pdf',bbox_inches='tight');pdf.savefig(fig,bbox_inches='tight');plt.close(fig)
    # Optional LaTeX table stays with diagnostic outputs, not in the manuscript.
    tex=['% Exploratory preview only: mean / shifted geometric mean, in seconds.',r'\begin{tabular}{lrrrr}',r'\toprule',
         'Method & Small & Medium & Large & All '+r'\\',r'\midrule']
    tex += [' & '.join(c)+r' \\' for c in cells]
    tex += [r'\bottomrule',r'\end{tabular}']
    (out/'size_table.tex').write_text('\n'.join(tex)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
