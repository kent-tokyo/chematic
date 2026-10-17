from pathlib import Path
import collections, hashlib, json
E=Path(__file__).resolve().parent
engines=('baseline','candidate','rdkit','cosmolkit')
data={name:[json.loads(line) for line in (E/f'broad-{name}.jsonl').read_text().splitlines()] for name in engines}
headers={name:[r for r in rows if 'sha256' in r] for name,rows in data.items()}
rows={name:[r for r in values if 'index' in r] for name,values in data.items()}
assert all(len(value)==14072 for value in rows.values())
assert rows['baseline']==rows['candidate'], 'source output or refusal changed'
assert headers['baseline']==headers['candidate']
report={'source_exact_rows_including_refusals':14072,'source_rows_equal':True,'headers':headers,'comparisons':{}}
for name in ('candidate','cosmolkit'):
    counts={};details=[]
    for src,ref in zip(rows[name],rows['rdkit'],strict=True):
        assert (src['corpus'],src['index'])==(ref['corpus'],ref['index'])
        for op in ('morgan_r2_chiral0','morgan_r2_chiral1','morgan_r3_chiral0','morgan_r3_chiral1'):
            a=src['values'].get(op,{});b=ref['values'].get(op,{})
            state=('both_unavailable' if 'value' not in a and 'value' not in b else 'source_unavailable' if 'value' not in a else 'rdkit_unavailable' if 'value' not in b else 'equal' if a['value']==b['value'] else 'mismatch')
            counts.setdefault(src['corpus']+'/'+op,collections.Counter())[state]+=1
            if state!='equal':details.append({'corpus':src['corpus'],'index':src['index'],'op':op,'state':state,'source':src,'rdkit':ref})
    report['comparisons'][name]={k:dict(v) for k,v in counts.items()}
    (E/f'broad-{name}-differences.json').write_text(json.dumps(details,indent=2)+'\n')
report['raw_sha256']={name:hashlib.sha256((E/f'broad-{name}.jsonl').read_bytes()).hexdigest() for name in engines}
(E/'broad-audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report['comparisons'],indent=2))
