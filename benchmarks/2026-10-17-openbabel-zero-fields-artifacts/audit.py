from pathlib import Path
import argparse,collections,csv,gzip,hashlib,itertools,json,random,statistics
E=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--archive-only',action='store_true');args=parser.parse_args()
plan=json.loads((E/'measurement-plan.json').read_text());decision=json.loads((E/'decision.json').read_text());pre=json.loads((E/'preflight.json').read_text())
assert decision['experiment_complete'] and pre['passed'] and pre['records']==4384 and pre['formal_observations']==0
MODES=('parse','write','roundtrip')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def corpus_sha(p):
    return sha(p) if p.exists() else hashlib.sha256(gzip.decompress(p.with_name(p.name+'.gz').read_bytes())).hexdigest()
def stats(a,b,seed):
    assert len(a)==len(b)==24 and all(type(x) is int and x>0 for x in a+b)
    ratios=[x/y for x,y in zip(a,b,strict=True)];rng=random.Random(seed)
    boot=sorted(statistics.median([ratios[rng.randrange(24)] for _ in range(24)]) for _ in range(10000))
    return statistics.median(ratios),[boot[249],boot[9749]],sum(r>1 for r in ratios)
lanes=[];rows=0;repeats=0;processed=0;checks={};all_journals=set()
for profile in plan['profiles']:
    kind,name,count=profile['kind'],profile['profile'],profile['records'];key=kind+'-'+name
    report=json.loads((E/f'paired-{key}.json').read_text());journal=E/(key+'-journals');records=pre['profiles'][key]
    assert report['complete'] and len(report['blocks'])==24 and report['corpus']==profile
    assert report['preflight']==records and corpus_sha(E/profile['path'])==profile['sha256']
    names=tuple(profile['orders'][0]);assert set(names)==({'baseline','candidate','openbabel'} if kind=='sdf' else {'baseline','candidate'})
    wanted=collections.Counter(itertools.permutations(names))
    actual=collections.Counter(tuple(x) for x in profile['orders'])
    assert set(wanted)==set(actual) and all(v==24//len(wanted) for v in actual.values())
    expected_files=set();units={}
    for arm in names:
        for mode in MODES:
            trial=records['instrumentation'][arm][mode];assert trial['output_units']%trial['repeats']==0
            units[arm,mode]=trial['output_units']//trial['repeats']
    for mode in MODES:assert units['baseline',mode]==units['candidate',mode]
    for index,block in enumerate(report['blocks']):
        assert block['order']==profile['orders'][index]
        for arm in block['order']:
            for mode in MODES:
                path=journal/f'{index:02}-{arm}-{mode}.json';expected_files.add(path.name);all_journals.add(str(path.relative_to(E)))
                raw=json.loads(path.read_text());row=raw['timing'];assert raw=={'profile':key,'block':index,'arm':arm,'mode':mode,'timing':block[arm][mode]}
                integer_fields=('elapsed_ns_raw','elapsed_ns','elapsed_ns_per_batch','repeats','records','output_units','records_per_batch','minimum_repeats','warmup_repeats','warmup_elapsed_ns')
                assert all(type(row[field]) is int and row[field]>0 for field in integer_fields)
                assert row['mode']==mode and row['repeats']>=3 and row['minimum_repeats']==3
                assert row['records']==row['repeats']*count and row['records_per_batch']==count
                assert row['elapsed_ns_raw']>=100_000_000 and row['elapsed_ns']==row['elapsed_ns_per_batch']==row['elapsed_ns_raw']//row['repeats']
                assert row['minimum_warmup_ns']==500_000_000 and row['minimum_timed_ns']==100_000_000 and row['warmup_elapsed_ns']>=500_000_000
                assert row['output_units']==units[arm,mode]*row['repeats'],(key,arm,mode,'different work units')
                rows+=1;repeats+=row['repeats'];processed+=row['records']
    assert {p.name for p in journal.iterdir()}==expected_files
    for suffix in ('snapshot','dump'):
        a=(E/f'preflight-{key}-baseline-{suffix}.gz');b=(E/f'preflight-{key}-candidate-{suffix}.gz')
        assert gzip.decompress(a.read_bytes())==gzip.decompress(b.read_bytes())
    checks[key]={'complete_blocks':24,'all_permutation_counts':{'/'.join(k):v for k,v in actual.items()},'every_journal_matches_full_result':True,'all_budgets_counts_normalizations_and_units_valid':True,'source_snapshots_and_warmed_rewrites_equal':True}
    for mode in MODES:
        candidate=[block['candidate'][mode]['elapsed_ns'] for block in report['blocks']]
        for arm,seed in [('baseline',0)]+([('openbabel',1)] if kind=='sdf' else []):
            base=[block[arm][mode]['elapsed_ns'] for block in report['blocks']];median,ci,wins=stats(base,candidate,seed)
            saved=report['summary'][mode]['source'] if arm=='baseline' else report['summary'][mode]['openbabel']['speedup']
            allowed=arm=='baseline' or report['summary'][mode]['openbabel']['semantic_prerequisite']
            if allowed:assert saved=={'median':round(median,3),'ci95':[round(ci[0],3),round(ci[1],3)],'winning_blocks':wins,'blocks':24}
            else:assert saved is None
            if arm=='baseline':assert report['summary'][mode]['source_unrounded_median']==median
            lanes.append({'profile':key,'mode':mode,'comparator':arm,'ranked':allowed,'paired_median':median,'ci95_low':ci[0],'ci95_high':ci[1],'winning_blocks':wins,'blocks':24,'baseline_median_ns':statistics.median(base),'candidate_median_ns':statistics.median(candidate),'strict_external_win':arm=='openbabel' and allowed and ci[0]>1 and wins==24})
assert rows==len(all_journals)==1584 and len(lanes)==36 and len(checks)==10
source=[r for r in lanes if r['comparator']=='baseline'];primary=[r for r in source if r['profile']+'/'+r['mode'] in plan['primary']];assert len(source)==30 and len(primary)==2
floor=[r for r in source if r['paired_median']<.99];adopt=not floor and all(r['ci95_low']>1 for r in primary)
assert adopt==decision['local_adoption_gate_passed']
if not args.archive_only:
    seal=json.loads((E/'seal.json').read_text())
    for path,h in seal['files_sha256'].items():assert sha(Path(path))==h,path
    root=Path('/Users/k_tanabe/.codex/worktrees/rdkit-cosmolkit-speed-20261017/chematic')
    for path,h in seal['source_files_sha256'].items():assert sha(root/path)==h,path
else:
    inventory=json.loads((E/'manifest.json').read_text())
    for path,h in inventory['files_sha256'].items():assert sha(E/path)==h,path
result={'experiment_complete':True,'independent_gate_passed':adopt,'checks':checks,'timing_rows':rows,'timed_complete_batches':repeats,'timed_processed_records':processed,'source_lanes':30,'all_source_records_equal':4384,'source_floor_failures':floor,'primary':primary,'minimum_source_median':min(source,key=lambda x:x['paired_median']),'all_source_medians_below_one':[r for r in source if r['paired_median']<1],'all_source_ci_upper_below_one':[r for r in source if r['ci95_high']<1],'ranked_external_lanes':sum(r['comparator']=='openbabel' and r['ranked'] for r in lanes),'strict_external_wins':[r for r in lanes if r['strict_external_win']],'seal_still_valid':True if not args.archive_only else 'original seal was validated before packaging; archive SHA inventory rechecked here','lanes':lanes}
(E/('archive-audit.json' if args.archive_only else 'audit.json')).write_text(json.dumps(result,indent=2)+'\n')
with (E/('archive-lanes.csv' if args.archive_only else 'lanes.csv')).open('w') as output:
    writer=csv.DictWriter(output,fieldnames=list(lanes[0]));writer.writeheader();writer.writerows(lanes)
print(json.dumps({k:v for k,v in result.items() if k not in ('checks','lanes','all_source_medians_below_one','strict_external_wins')},indent=2))
