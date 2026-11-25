#!/usr/bin/env python3
"""Summarize five-wheel exposed 3D probe without excluding failed inputs."""
import argparse
import hashlib
import json
import statistics
from collections import Counter
from pathlib import Path
from summarize_cod_candidate_v2 import recovered


def main():
    p=argparse.ArgumentParser(); p.add_argument('input',type=Path); p.add_argument('--output',type=Path,required=True); args=p.parse_args()
    records=[json.loads(line) for line in args.input.read_text().splitlines()]; head=records[0]; rows=records[1:-1]; marker=records[-1]
    if marker.get('type')!='complete' or marker['rows']!=len(rows) or head['total']!=len(rows): raise ValueError('incomplete run')
    if len({r['index'] for r in rows})!=len(rows): raise ValueError('duplicate inputs')
    out={'source_sha256':hashlib.sha256(args.input.read_bytes()).hexdigest(),'script_sha256':head['script_sha256'],'rows':len(rows),'engines':{}}
    for engine in head['engines']:
        vals=[row.get('engines',{}).get(engine,{}) for row in rows]
        info={'statuses':dict(Counter(v.get('status','missing') for v in vals)), 'input_errors':sum('input_error' in row for row in rows), 'scoring_errors':sum('scoring_error' in v for v in vals),'rmsd_missing':sum(len(v.get('rmsd_errors',[])) for v in vals), 'attempts':sum(v.get('attempted',0) for v in vals),'conformers':sum(v.get('accepted',0) for v in vals),'converged':sum(s==0 for v in vals for s in v.get('statuses',[]))}
        for label, cutoff, require_complete in [('fixed',None,True),('timed_prefix',head['timed_prefix_seconds'],False)]:
            count=0; converged=0
            for row in rows:
                raw=row.get('engines',{}).get(engine,{})
                if require_complete and raw.get('status')!='ok': continue
                item=dict(raw)
                if cutoff is not None:
                    keep=[i for i,t in enumerate(item.get('arrival_seconds',[])) if t<=cutoff]
                    for key in ('rmsd','stereo','geometry','statuses'):
                        arr=item.get(key,[]); item[key]=[arr[i] for i in keep if i<len(arr)]
                copy=dict(row); copy['engines']={engine:item}
                count+=recovered(copy,engine); converged+=recovered(copy,engine,True)
            info[label]={'recovered':count,'converged_recovered':converged}
        times=[v['elapsed_with_cleanup'] for v in vals if 'elapsed_with_cleanup' in v]; mem=[v['worker_peak_rss'] for v in vals if v.get('worker_peak_rss') is not None]
        info['elapsed_seconds']={'count':len(times),'median':statistics.median(times) if times else None,'max':max(times,default=None)}
        info['worker_peak_rss_bytes']={'count':len(mem),'median':statistics.median(mem) if mem else None,'max':max(mem,default=None)}
        stereo=[s for v in vals for s in v.get('stereo',[])]
        info['stereo']={'checked':len(stereo),'violated':sum(bool(s.get('violated_atoms') or s.get('violated_bonds')) for s in stereo),'unevaluable':sum(bool(s.get('unevaluable_atoms') or s.get('unevaluable_bonds')) for s in stereo),'coverage_missing':sum(v.get('accepted',0)-len(v.get('stereo',[])) for v in vals)}
        out['engines'][engine]=info
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n'); print(json.dumps(out,indent=2))

if __name__=='__main__': main()
