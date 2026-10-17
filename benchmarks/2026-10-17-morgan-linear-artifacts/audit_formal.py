from pathlib import Path
import collections, csv, gzip, hashlib, itertools, json, random, statistics
E=Path(__file__).resolve().parent
ARMS=('baseline','candidate','rdkit','cosmolkit')
PRIMARY=('morgan2_2048','morgan2_chiral','morgan3_2048')
CONTROLS=('tpsa','labute_asa','num_rings','maccs')
MODES=('first_use','repeated_use','parse_inclusive')
def sha(raw):return hashlib.sha256(raw).hexdigest()
def encoded(v):return json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
def statistics_of(a,b,seed):
    ratios=[x/y for x,y in zip(a,b,strict=True)]
    rng=random.Random(seed)
    boot=sorted(statistics.median([ratios[rng.randrange(24)] for _ in range(24)]) for _ in range(10000))
    return statistics.median(ratios),[boot[249],boot[9749]],sum(x>1 for x in ratios)
lanes=[];journals=0;accounting=0;total_batches=0;total_rows=0;checks={}
for offset in (0,5000):
    directory=E/f'profile-{offset}'
    result=json.loads((directory/'result.json').read_text())
    plan=json.loads((directory/'plan.json').read_text())
    assert result['complete'] and len(result['blocks'])==24 and plan['rows']==1000
    orders=list(itertools.permutations(ARMS));random.Random(20261017+offset).shuffle(orders)
    assert plan['orders']==[list(x) for x in orders]
    assert len(set(orders))==24
    for name in ARMS:
        assert [sum(order[position]==name for order in orders) for position in range(4)]==[6]*4
    expected_files={f'{i:02}-{name}.json.gz' for i in range(24) for name in ARMS}
    assert {p.name for p in directory.glob('*.json.gz')}==expected_files
    assert not list(directory.glob('*.partial'))
    reference={}
    raw={name:collections.defaultdict(list) for name in ARMS}
    for i,order in enumerate(orders):
        block=result['blocks'][i]
        assert block['order']==list(order)
        for name in order:
            sample=json.loads(gzip.decompress((directory/f'{i:02}-{name}.json.gz').read_bytes()));journals+=1
            assert sample['label']==name and sample['offset']==offset and sample['rows']==1000 and not sample['probe_only']
            assert sample['corpus_sha256']==plan['corpus_sha256']
            assert {k:sample[k] for k in ('version','module_path','native_sha256')}==plan['source_manifest']['expected_metadata'][name]
            compact=json.loads(json.dumps(sample))
            values={op:record['values'] for op,record in sample['operations'].items()}
            if name not in reference:reference[name]=values
            assert reference[name]==values
            expected_ops=PRIMARY+CONTROLS if name in ARMS[:2] else PRIMARY
            assert set(values)==set(expected_ops)
            for op,record in sample['operations'].items():
                assert len(values[op])==1000
                if op in PRIMARY or op=='maccs':
                    assert all(bits==sorted(set(bits)) and all(type(bit) is int and 0<=bit<2048 for bit in bits) for bits in values[op])
                value_hash=sha(encoded(values[op]))
                for mode in MODES:
                    detail=record['accounting'][mode];times=detail['batch_ns'];count=len(times)
                    budget=50_000_000 if op in PRIMARY else 20_000_000
                    assert detail['mode']==mode and detail['minimum_accumulated_ns']==budget
                    assert count>0 and all(type(t) is int and t>0 for t in times) and sum(times)>=budget
                    assert detail['elapsed_ns_raw']==sum(times) and detail['batches']==count and detail['processed_rows']==1000*count
                    assert detail['elapsed_ns_per_batch']==sum(times)//count==record[mode]
                    assert detail['output_sha256_per_batch']==[value_hash]*count
                    raw[name][op+'/'+mode].append(record[mode]);accounting+=1;total_batches+=count;total_rows+=1000*count
                compact['operations'][op].pop('values')
            assert compact==block[name]
        for op in PRIMARY+CONTROLS:assert reference['baseline'][op]==reference['candidate'][op]
        for name in ('baseline','candidate','cosmolkit'):
            for op in PRIMARY:assert reference[name][op]==reference['rdkit'][op]
    checks[str(offset)]={'all_24_permutations_once':True,'each_arm_position_count':6,'all_source_and_external_outputs_equal':True,'no_missing_rows_or_batches':True,'journal_count':96}
    for op in PRIMARY+CONTROLS:
        for mode in MODES:
            key=op+'/'+mode;b=raw['candidate'][key]
            for seed,name in enumerate(('baseline','rdkit','cosmolkit')):
                if name!='baseline' and op not in PRIMARY:continue
                a=raw[name][key];assert len(a)==len(b)==24
                median,ci,wins=statistics_of(a,b,seed)
                saved=result['summary'][key][name]
                assert saved['outputs_equal']
                assert saved['unrounded_paired_median']==median
                assert saved['speedup']=={'median':round(median,3),'ci95':[round(ci[0],3),round(ci[1],3)],'winning_blocks':wins,'blocks':24}
                lanes.append({'offset':offset,'operation':op,'mode':mode,'comparator':name,'paired_median':median,'ci95_low':ci[0],'ci95_high':ci[1],'winning_blocks':wins,'blocks':24,'comparator_median_ns':statistics.median(a),'candidate_median_ns':statistics.median(b),'strict_external_win':name!='baseline' and ci[0]>1 and wins==24})
source=[r for r in lanes if r['comparator']=='baseline'];assert len(source)==42 and len(lanes)==78
primaries=[r for r in source if r['operation']=='morgan2_chiral' and r['mode']=='repeated_use'];assert len(primaries)==2
floor=[r for r in source if r['paired_median']<.99]
broad=json.loads((E/'broad-audit.json').read_text());assert broad['source_rows_equal'] and all(record=={'equal':10000 if key.startswith('rdkit-js') else 4072} for key,record in broad['comparisons']['candidate'].items())
report={'complete':True,'checks':checks,'journals':journals,'accounting_cells':accounting,'total_batches':total_batches,'total_processed_rows_in_timed_batches':total_rows,'source_lane_count':len(source),'all_lane_count':len(lanes),'primary_source_lanes':primaries,'source_floor_violations':floor,'source_ci_upper_below_one':[r for r in source if r['ci95_high']<1],'external_strict_wins':[r for r in lanes if r['strict_external_win']],'broad_source_and_rdkit_equal_cells':56288,'source_adoption_passed':all(r['ci95_low']>1 for r in primaries) and not floor,'lanes':lanes}
(E/'formal-audit.json').write_text(json.dumps(report,indent=2)+'\n')
with (E/'formal-lanes.csv').open('w') as output:
    writer=csv.DictWriter(output,fieldnames=list(lanes[0]));writer.writeheader();writer.writerows(lanes)
print(json.dumps({k:v for k,v in report.items() if k not in ('lanes','external_strict_wins')},indent=2))
print('external strict wins',len(report['external_strict_wins']),'of36')
