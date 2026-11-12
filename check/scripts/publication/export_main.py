"""Synchronize the approved 130-instance ablation table and existing profile style."""
import json,re,sys,math,statistics
from pathlib import Path
import pandas as pd
import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--input', type=Path, required=True)
parser.add_argument('--threads-input', type=Path, required=True)
parser.add_argument('--paper', type=Path, required=True)
parser.add_argument('--diagnostics', type=Path, required=True)
args = parser.parse_args()
HERE = args.input
PAPER = args.paper
TEMPLATES = Path(__file__).with_name('templates')
(PAPER/'data/across_all/FINAL/full').mkdir(parents=True, exist_ok=True)
args.diagnostics.mkdir(parents=True, exist_ok=True)
from plot_style import _compute_slowdown, _ecdf
import numpy as np
rows=json.loads((HERE/'snapshot.json').read_text())
report=json.loads((HERE/'summary-excluding-case162.json').read_text())
methods=[r['method'] for r in report['summary']]
# Keep raw solver output intact; apply reviewed validity decisions to publication statistics.
overrides=json.loads((HERE/'publication-overrides.json').read_text())
for override in overrides:
 matched=[r for r in rows if (r['method'],r['case'])==(override['method'],override['case'])]
 assert len(matched)==1
 matched[0]['issues'].append(override['reason'])
summary=[]
for method in methods:
 rs=[r for r in rows if r['method']==method and r['case']!='pglib_scopf_case162']
 valid=[r['result']['termination_status']=='OPTIMAL' and not r['issues'] for r in rs]
 times=[r['clock'] if ok else 10800 for r,ok in zip(rs,valid)]
 errors=sum(bool(r['issues']) for r in rs)
 summary.append(dict(method=method,n=len(rs),solved=sum(valid),limits=len(rs)-sum(valid)-errors,errors=errors,
                     mean=statistics.mean(times),geomean=math.exp(statistics.mean(math.log(t+10) for t in times))-10))
report['summary']=summary
report['publication_overrides']=overrides
(args.diagnostics/'summary-excluding-case162.json').write_text(json.dumps(report,indent=2)+'\n')
# Preserve the curated table's configurations, order, and typography.
p=PAPER/'data/across_all/FINAL/full/comparison_short_all.tex'
s=(TEMPLATES/p.name).read_text();s=s.replace('Sol. & TL & Err.', 'Sol. & Lim. & Err.')
lines=s.splitlines();i=0
for k,line in enumerate(lines):
 if re.match(r'^(?:[1-8][ab]? &|\\textcolor\{gray\}\{8[ab]\})',line):
  r=summary[i];i+=1;parts=line.split(' & ');assert len(parts)==7
  nums=[str(r['solved']),f"{r['mean']:.2f}",f"{r['geomean']:.2f}"]
  if r['method'] in ['sequential','hybrid']:
   nums=[r'\textcolor{gray}{'+x+'}' for x in nums]
   parts[:4]=[x if x.startswith(r'\textcolor{gray}{') else r'\textcolor{gray}{'+x+'}' for x in parts[:4]]
  lines[k]=' & '.join(parts[:4]+nums)+r' \\'
 elif line.startswith(r'\textcolor{gray}{9}'):
  parts=line.split(' & ')
  thread_summary=args.threads_input
  if thread_summary.exists():
   r=json.loads(thread_summary.read_text());assert r['instances']==r['solved']==130 and r['publication_acceptance']['all_optimal']
   parts[1]=r'\textcolor{gray}{7}'
   nums=[str(r['solved']),f"{r['mean']:.2f}",f"{r['geomean']:.2f}"]
   nums=[r'\textcolor{gray}{'+x+'}' for x in nums]
   lines[k]=' & '.join(parts[:4]+nums)+r' \\'
  else:
   lines[k]=' & '.join(parts[:4])+r' & \multicolumn{3}{c}{\textcolor{gray}{Pending}} \\'
assert i==10
s='\n'.join(lines)+'\n';s=s.replace('% n = 131;', '% n = 130; excluded pglib_scopf_case162;')
s=re.sub(r'% Multithreading source: .*', '% Multithreading source: 20261022-final-3h-ablation-branchfix-benders-threads/status-audit/summary.json', s)
s = re.sub(r'% Source: 20261021-final-3h-ablation-heuristics-budget/status-audit/summary-excluding-case162.json\n', '', s)
s += '% Source: 20261021-final-3h-ablation-heuristics-budget/status-audit/summary-excluding-case162.json\n'
p.write_text(s)
# Use the existing slowdown implementation and retain all TikZ styling.
shown=[m for m in methods if m not in ['sequential','hybrid']]
df=pd.DataFrame([{'instance':r['case'],'settings':r['method'],'status':r['result']['termination_status'] if not r['issues'] else 'ERROR','time_solve':r['clock']} for r in rows if r['case']!='pglib_scopf_case162' and r['method'] in shown])
assert len(df)==130*len(shown)
slow=_compute_slowdown(df,'time_solve',shown)
p=PAPER/'data/across_all/FINAL/full/performance_profile_time_solve.tikz.tex'
s=(TEMPLATES/p.name).read_text()
s=re.sub(r'\\addplot\+\[[^\n]+\] coordinates \{\n.*?\n\};\n\\addlegendentry\{([^}]+)\}\n', lambda m: '' if m[1]=='Heuristic warmstart' else m[0], s, flags=re.S)
s=re.sub(r'^\\node[^\n]*\$n=\d+\$[^\n]*\n', '', s, flags=re.M)
pattern=r'(\\addplot\+\[[^\n]+\] coordinates \{\n).*?(\n\};)';matches=list(re.finditer(pattern,s,re.S));assert len(matches)==len(shown)
for m,method in reversed(list(zip(matches,shown))):
 vals=slow.loc[slow.settings==method,'slowdown'].to_numpy();assert len(vals)==130
 grid=np.unique(np.concatenate(([1.0],vals[(vals>=1)&(vals<=100)],[100.0])))
 ys=_ecdf(vals,grid,n_total=130)
 coords='\n'.join(f'  ({x:.12g},{y:.12g})' for x,y in zip(grid,ys))
 s=s[:m.start()]+m[1]+coords+m[2]+s[m.end():]
s=s.replace('% Generated by check/scripts/visu/perf_profile.py; do not edit by hand.','% Generated by status-audit/export_paper_130.py using check/scripts/visu/perf_profile.py.\n% n = 130; excluded pglib_scopf_case162; sequential times include warmstart.')
s=s.replace('$n=131$', '$n=130$')
p.write_text(s)
slow.to_csv(args.diagnostics/'profile-130.csv',index=False)
print('Updated table and performance profile for 130 instances.')
