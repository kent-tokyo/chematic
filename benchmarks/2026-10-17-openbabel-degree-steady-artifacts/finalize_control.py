from pathlib import Path
import hashlib,json,shutil
O=Path(__file__).resolve().parent
R=Path('/private/tmp/chematic-degree-only-steady-20261017')
NAME='2026-10-17-openbabel-degree-steady'
D=R/'benchmarks'/(NAME+'-artifacts')
assert json.loads((O/'decision.json').read_text())['source_valid_epochs']==10
rows=[]
for design,count in [('cyclic',21),('balanced',24)]:
 data=json.loads((O/('same-binary-'+design+'.json')).read_text())
 assert data['complete'] and len(data['blocks'])==count and data['same_rust_executable_for_both_labels']
 for mode,ratio in data['source_label_ratios'].items():
  rows.append(f"| {design} | {mode} | {ratio['median']:.3f} | {ratio['ci95'][0]:.3f} .. {ratio['ci95'][1]:.3f} | {ratio['winning_blocks']}/{count} |")
 for block in data['blocks']:
  for arm in ['baseline','candidate','openbabel']:
   for mode in ['parse','write','roundtrip']:
    x=block[arm][mode]
    assert x['records']==1000*x['repeats'] and x['elapsed_ns']==x['elapsed_ns_raw']//x['repeats']
 for name,digest in data['binary_sha256'].items():
  assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==digest
section='''## Same-binary diagnostic, excluded from adoption

After all ten formal candidate epochs completed, both Rust labels execute the exact same `baseline-sdf` path and bytes. The same Open Babel third arm, unfiltered SDF-0 corpus, modes, warmup and time budgets are retained. The two designs were fixed in `control-plan.json` before any control samples: 21 cyclic rotations and 24 blocks spanning all six permutations, each used four times with one fixed seed. No control is repeated and no control observations enter the candidate ratios, retention gates or external rankings.

A nonunit ratio here indicates possible order/environment effects for identical executable bytes. A ratio consistent with one does not prove the candidate slowdown arose from code layout. The failed candidate remains unadopted regardless of these diagnostic observations.

| Design | Mode | First label / second label median | 95% CI | Second label wins |
|---|---|---:|---|---:|
'''+ '\n'.join(rows)+'\n\n'
p=R/'benchmarks'/(NAME+'.md');text=p.read_text();marker='## Evidence and reproduction';assert text.count(marker)==1 and '## Same-binary diagnostic' not in text;p.write_text(text.replace(marker,section+marker))
for name in ['control-plan.json','control.py','control.log','same-binary-cyclic.json','same-binary-balanced.json','finalize_control.py']:
 shutil.copy2(O/name,D/name)
packet={str(p.relative_to(D)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(D.rglob('*')) if p.is_file() and p.name!='packet.sha256.json'}
(D/'packet.sha256.json').write_text(json.dumps(packet,indent=2)+'\n')
for name,digest in packet.items():assert hashlib.sha256((D/name).read_bytes()).hexdigest()==digest
print('Final packet verified:',len(packet),'files')
print('\n'.join(rows))
