from pathlib import Path
import json,hashlib,shutil,time,subprocess
old=Path('/private/tmp/chematic-openbabel-count-rounded-20261015')
while True:
 status=(old/'measurement.log').read_text() if (old/'measurement.log').exists() else ''
 if 'Traceback' in status:raise RuntimeError('Rounded-capacity measurement failed')
 if 'source_diagnostics passed' in status:break
 time.sleep(2)
negative={p.stem:{k:v['speedup'] for k,v in json.loads(p.read_text())['summary'].items() if v['speedup']['ci95'][1]<1} for p in old.glob('focused21-*.json')};assert any('write' in lanes for lanes in negative.values()),'Review rounded candidate before starting another implementation'
(old/'decision.json').write_text(json.dumps({'adopt_source':False,'reason':'Reject rounded capacities alone: original profile parser/roundtrip gains are retained, but prepared-writer regressions remain in at least one prescribed separate follow-up. Rounded capacities did not eliminate the prepared-operation cost. All output and 848 native/10 runner gates pass. Archive all initial/continuation/follow-up results separately. The next implementation adds a direct writer optimization for the common zero atom-map field; it will be tested and measured as a separate full candidate, never pooled with this run.','reproduced_negative_lanes':negative,'merged':False,'published':False},indent=2)+'\n')
r=Path('/private/tmp/chematic-openbabel-count-zero-map-20261015');s=Path('/private/tmp/chematic-openbabel-count-zero-map-source-20261015');assert not r.exists() and not s.exists();r.mkdir();s.mkdir();base=json.loads((old/'baseline-source.json').read_text());active=Path('/Users/k_tanabe/.codex/worktrees/relevant-cycle-speed/chematic')
for n,h in base['sha256'].items():
 data=(active/n).read_bytes();assert hashlib.sha256(data).hexdigest()==h,n;p=s/n;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data)
(r/'baseline-source.json').write_text(json.dumps(base,indent=2)+'\n')
for n in ['baseline-v3000','baseline-build.log','baseline-mol3000.rs','baseline-CHANGELOG.md','baseline-artifact-origin.json','continuation-plan.json','continuation-corpus-manifest.json','corpus-continued-chains.molblocks']:shutil.copyfile(old/n,r/n)
(r/'baseline-v3000').chmod(0o755);(r/'baseline-writer.rs').write_bytes((s/'crates/chematic-mol/src/mol3000/writer.rs').read_bytes())
for n in ['optimize.py','round.py','resume-validation.py','freeze.py','run_formal.py','source_diagnostics.py','measure.py','after-validation.py','deliver.py','package.py']:
 t=(old/n).read_text().replace('count-rounded-source-20261015','count-zero-map-source-20261015').replace('count-rounded-target-20261015','count-zero-map-target-20261015').replace('count-rounded-20261015','count-zero-map-20261015').replace('openbabel-v3000-count-rounded','openbabel-v3000-count-zero-map')
 if n=='resume-validation.py':t=t.replace("('rounding.log',['python3',str(r/'round.py')]),", "('rounding.log',['python3',str(r/'round.py')]),('zero-map.log',['python3',str(r/'zero-map.py')]),")
 if n=='freeze.py':
  t=t.replace("assert changed==['CHANGELOG.md','crates/chematic-mol/src/mol3000.rs']","assert changed==['CHANGELOG.md','crates/chematic-mol/src/mol3000.rs','crates/chematic-mol/src/mol3000/writer.rs']")
  t=t.replace("'crates/chematic-mol/src/mol3000.rs':'baseline-mol3000.rs'","'crates/chematic-mol/src/mol3000.rs':'baseline-mol3000.rs','crates/chematic-mol/src/mol3000/writer.rs':'baseline-writer.rs'")
  t=t.replace("'preceding_rejected_trial':{", "'preceding_rounded_trial':{'manifest_sha256':sha(Path('/private/tmp/chematic-openbabel-count-rounded-20261015/source-manifest.json')),'decision_sha256':sha(Path('/private/tmp/chematic-openbabel-count-rounded-20261015/decision.json'))},'preceding_rejected_trial':{")
 if n=='package.py':
  t=t.replace("('crates/chematic-mol/src/mol3000.rs','baseline-mol3000.rs')]", "('crates/chematic-mol/src/mol3000.rs','baseline-mol3000.rs'),('crates/chematic-mol/src/mol3000/writer.rs','baseline-writer.rs')]")
  t=t.replace('V3000 bounded rounded capacity — 2026-10-15','V3000 bounded capacities and zero atom-map writing — 2026-10-15')
  t=t.replace('The writer is unchanged; writer timing movement is not attributed to a new writer algorithm.','The 2D writer directly appends the exact constant zero-map suffix for absent or zero atom maps; nonzero maps retain the original formatter. The parser and writer effects are combined in this measured candidate, without an isolated attribution of each effect.')
  t=t.replace('No public API, dependency, tolerance or unsafe code is added.','The 2D writer also bypasses integer formatting for the common zero atom-map field, appending the same exact bytes; nonzero maps keep the existing formatter. Existing frozen-writer comparisons exercise absent, zero, nonzero and maximum atom maps. No public API, dependency, tolerance or unsafe code is added.')
 (r/n).write_text(t)
p=json.loads((old/'measurement-plan.json').read_text());p['revision_reason']='Rounded parser capacities alone retain prepared-writer regressions. Add a direct constant suffix for absent/zero atom maps in the 2D writer; all preceding rejected trial outcomes retained separately.';p['intended_improvement']='V3000 parse/write/roundtrip via bounded reservation and zero atom-map output shortcut; full fixed matrix rerun';(r/'measurement-plan.json').write_text(json.dumps(p,indent=2)+'\n')
audit={};corpora=Path('/private/tmp/chematic-openbabel-v3000-20261014')
for profile in ['0','5000','chains']:
 path=corpora/f'corpus-{profile}.molblocks';inside=False;counts={}
 for raw in path.read_text().splitlines():
  if raw=='M  V30 BEGIN ATOM':inside=True;continue
  if raw=='M  V30 END ATOM':inside=False;continue
  if inside:
   tokens=raw.removeprefix('M  V30 ').split();counts[tokens[5]]=counts.get(tokens[5],0)+1
 audit[profile]={'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'atom_map_lexemes':counts,'atoms':sum(counts.values())}
(r/'atom-map-audit.json').write_text(json.dumps(audit,indent=2)+'\n')
(r/'zero-map.py').write_text('''from pathlib import Path
s=Path('/private/tmp/chematic-openbabel-count-zero-map-source-20261015')
p=s/'crates/chematic-mol/src/mol3000/writer.rs';t=p.read_text();old='''+repr('''        write!(out, " 0.0000 {}", atom.atom_map.unwrap_or(0))
            .expect("writing to String cannot fail");''')+''';assert old in t;t=t.replace(old,''' + repr('''        if atom.atom_map.unwrap_or(0) == 0 {
            out.push_str(" 0.0000 0");
        } else {
            write!(out, " 0.0000 {}", atom.atom_map.unwrap_or(0))
                .expect("writing to String cannot fail");
        }''')+''',1);p.write_text(t)
p=s/'CHANGELOG.md';p.write_text(p.read_text().replace('  state transitions, errors and diagnostics.','  state transitions, errors and diagnostics. Write absent or zero V3000 2D\n  atom maps directly, retaining the original nonzero-map formatting.'))
print('Zero atom-map direct suffix written')
''')
print('Combined writer candidate prepared after completed rounded-trial rejection',flush=True)
with (r/'validation-pipeline.log').open('w') as log:subprocess.run(['python3',str(r/'resume-validation.py')],stdout=log,stderr=subprocess.STDOUT,check=True)
print('Combined candidate all validation passed',flush=True)
