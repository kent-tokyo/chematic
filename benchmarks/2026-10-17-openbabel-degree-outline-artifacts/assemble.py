from pathlib import Path
import gzip,hashlib,json,statistics,shutil,subprocess,sys
O=Path(__file__).resolve().parent;C=Path('/private/tmp/chematic-degree-only-steady-evidence-20261017');R=Path('/private/tmp/chematic-degree-only-steady-20261017');NAME='2026-10-17-openbabel-degree-outline';D=R/'benchmarks'/(NAME+'-artifacts');D.mkdir(exist_ok=True)
def sha(raw):return hashlib.sha256(raw).hexdigest()
plan=json.loads((O/'measurement-plan.json').read_text());decision=json.loads((O/'decision.json').read_text());validation=json.loads((O/'validation.json').read_text());manifest=json.loads((O/'source-manifest.json').read_text());proof=json.loads((O/'binary-source-proof.json').read_text());base=manifest['base_commit'];candidate=manifest['candidate_commit']
sys.path.insert(0,str(R/'scripts'));from bench_openbabel_sdf_throughput import paired_median_ratio_ci,eligible
src=subprocess.check_output(['git','show',base+':crates/chematic-perception/src/stereo_validation.rs'],cwd=R,text=True);a=src.index('pub fn stereo_centers(');b=src.index('\n/// Summarise',a);old=src[a:b].strip()
test=(R/'crates/chematic-perception/src/stereo_degree_tests.rs').read_text();a=test.index('fn frozen_stereo_centers(');b=test.index('\nfn compare(',a);assert old==test[a:b].replace('fn frozen_stereo_centers(','pub fn stereo_centers(',1).strip()
(O/'frozen-classifier-proof.json').write_text(json.dumps({'base_commit':base,'complete_frozen_classifier_exact':True,'body_sha256':sha(old.encode()),'pinned_graphs':14072,'five_atom_graph_H_variants':6144,'necessary_condition':'Every emitted shape requires three or four heavy neighbors; their count cannot exceed graph degree. Ring-cache preparation and valid branches remain.'},indent=2)+'\n')
rows=[];audits=[];floors=[];negative=[];primary={};samples=0;external_count=0
for p in plan['profiles']:
 key=p['kind']+'-'+p['profile'];data=json.loads((O/('paired-'+key+'.json')).read_text());assert data['complete'] and len(data['blocks'])==24 and data['schema']=='balanced-steady-throughput/v1'
 assert data['preflight']['snapshot_sha256']['baseline']==data['preflight']['snapshot_sha256']['candidate'] and data['preflight']['rewrite_sha256']['baseline']==data['preflight']['rewrite_sha256']['candidate']
 assert [r['order'] for r in data['blocks']]==p['orders']
 arms=p['orders'][0]
 for arm in arms:
  for pos in range(len(arms)):assert sum(r['order'][pos]==arm for r in data['blocks'])==24//len(arms)
 for block in data['blocks']:
  for arm in arms:
   for mode in ['parse','write','roundtrip']:
    x=block[arm][mode];n=x['repeats'];raw=x['elapsed_ns_raw']
    assert type(n) is int and n>=3 and raw>=100000000 and x['elapsed_ns']==x['elapsed_ns_per_batch']==raw//n>0
    assert x['records']==n*p['records'] and x['records_per_batch']==p['records']
    assert x['warmup_repeats']>=1 and x['warmup_elapsed_ns']>=500000000
    assert x['minimum_warmup_ns']==500000000 and x['minimum_timed_ns']==100000000
    samples+=1
 for mode in ['parse','write','roundtrip']:
  values={arm:[b[arm][mode]['elapsed_ns'] for b in data['blocks']] for arm in arms}
  ratio=paired_median_ratio_ci(values['baseline'],values['candidate'],seed=0);median=statistics.median(a/b for a,b in zip(values['baseline'],values['candidate'],strict=True))
  assert ratio==decision['source_ratios'][key+'/'+mode]==data['summary'][mode]['source'] and median==data['summary'][mode]['source_unrounded_median']
  if median<.99:floors.append(key+'/'+mode)
  if ratio['ci95'][1]<1:negative.append(key+'/'+mode)
  if key+'/'+mode in plan['primary']:primary[key+'/'+mode]=ratio
  def fmt(x):return f"{x['median']:.3f} [{x['ci95'][0]:.3f}, {x['ci95'][1]:.3f}]; {x['winning_blocks']}/24"
  external='not measured in this candidate'
  if 'openbabel' in arms:
   entry=data['summary'][mode]['openbabel'];assert entry==decision['fresh_SDF_Open_Babel'][key+'/'+mode]
   allowed=eligible(data['preflight']['semantics'],True,mode,'openbabel');assert allowed==entry['semantic_prerequisite']
   speed=paired_median_ratio_ci(values['openbabel'],values['candidate'],seed=1) if allowed else None;assert speed==entry['speedup'];external_count+=1
   external=fmt(speed)+(' strict win' if entry['strict_speed_win'] else '') if speed else 'unranked: semantic prerequisite fails'
  rows.append('| '+key+' | '+mode+' | '+f"{statistics.median(values['baseline'])/1e6:.3f} | {statistics.median(values['candidate'])/1e6:.3f} | "+fmt(ratio)+' | '+external+' |')
  audits.append({'lane':key+'/'+mode,'unrounded_source_median':median,'ratio_recomputed':True})
assert samples==1584 and len(audits)==30 and external_count==6
assert floors==[x['lane'] for x in decision['source_floor_failures']] and negative==[x['lane'] for x in decision['negative_intervals']]
assert decision['local_adoption_gate_passed']==(len(primary)==2 and all(x['ci95'][0]>1 for x in primary.values()) and not floors)
(O/'audit.json').write_text(json.dumps({'passed':True,'all_1584_sample_counts_budgets_normalizations_checked':True,'all_30_source_ratios_recomputed':True,'all_SDF_external_eligibility_and_ranked_ratios_recomputed':True,'planned_orders_and_position_balance_verified':True,'no_filtered_rows_or_reruns':True,'comparisons':audits},indent=2)+'\n')
(O/'candidate.patch').write_bytes(subprocess.check_output(['git','diff',base,candidate,'--','crates','tools','scripts'],cwd=R))
for name in ['baseline-build-command.json']:
 shutil.copy2(C/name,O/('reused-'+name))
for p in sorted(O.iterdir()):
 if not p.is_file() or p.name=='assemble.log':continue
 name=p.name
 if p.suffix in ['.json','.py'] or name=='common-harness.patch':shutil.copy2(p,D/name)
 elif p.suffix in ['.log','.txt'] or name=='candidate.patch':(D/(name+'.gz')).write_bytes(gzip.compress(p.read_bytes(),mtime=0))
for archive in sorted(O.glob('unmeasured-*-preflight')):
 if archive.is_dir():shutil.copytree(archive,D/archive.name,dirs_exist_ok=True)
for p in plan['profiles']:
 key=p['kind']+'-'+p['profile'];(D/(p['path']+'.gz')).write_bytes(gzip.compress((O/p['path']).read_bytes(),mtime=0))
 for mode in ['snapshot','dump']:
  left=O/f'preflight-{key}-baseline-{mode}.gz';right=O/f'preflight-{key}-candidate-{mode}.gz';raw=gzip.decompress(left.read_bytes());assert raw==gzip.decompress(right.read_bytes())
  shutil.copy2(left,D/(key+'-source-common-'+mode+'.gz'))
 for mode in ['dump','smiles']:
  path=O/f'preflight-{key}-openbabel-{mode}.gz'
  if path.exists():shutil.copy2(path,D/path.name)
common=manifest['common_probe_files_sha256']
for name,h in common.items():
 assert sha((R/name).read_bytes())==sha((C/'baseline-source'/name).read_bytes())==h
 path=D/'common-instrumentation'/name;path.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(R/name,path)
(D/'README.md').write_text('Fresh balanced source experiment: ten profiles, 30 source lanes, 24 blocks; fresh Open Babel timings only on two SDF profiles. Exactly 1584 native samples retain raw elapsed, actual repeats, processed records, budgets and normalized ns/full-batch. Prior C outcomes are referenced only for semantic provenance, never for elapsed statistics. Common source outcomes are byte-identical for both source arms. Full binaries remain at '+str(O)+'. Scripts retain original paths and require adaptation.\n')
packet={str(p.relative_to(D)):sha(p.read_bytes()) for p in sorted(D.rglob('*')) if p.is_file() and p.name!='packet.sha256.json'};(D/'packet.sha256.json').write_text(json.dumps(packet,indent=2)+'\n')
status='PASS: local source adoption gate only; hosted CI and main integration are separate' if decision['local_adoption_gate_passed'] else 'FAIL: outlined candidate not adopted'
env=json.loads((O/'environment.json').read_text())
report=f'''# Balanced throughput: outlined stereocenter classification

Decision: **{status}**. Baseline `{base}`; measured source `{candidate}`. One macOS host, native source builds, steady-state throughput. This does not establish universal superiority, cold-process latency, physical-conformer quality or published-wheel performance.

## Source and correctness

The production change is the necessary degree>=3 precondition in `stereo_centers`, plus `inline(never)` on that classifier. Every shape emitted by the encoded baseline classifier needs three or four heavy neighbors, whose count cannot exceed graph degree. The guard preserves global ring-cache preparation, valid-center branches and final shape tests. Numeric formatting is unchanged. No API, dependency or production unsafe additions.

The full baseline classifier is frozen byte-exact after renaming: 14,072 pinned parsed graphs, 6,144 five-atom/H variants and directed boundaries match. All 4,384 fresh native parse snapshots and warmed rewrites are byte-identical across source arms and match the preceding C baseline streams. Fresh SDF observer assessment also reproduces all accepted rows and mismatch counts; semantic failures remain unranked. These establish preservation of existing behavior within the stated tests, not chemical correctness of all possible inputs.

Fresh native validation: {validation['native']['passed']} passed, {validation['native']['failed']} failed, {validation['native']['ignored']} existing ignored in {validation['native']['suites']} suites (perception, MOL, chemistry, fingerprints). All-target release Clippy, formatting, unchanged refactor ratchet and 24 runner/accounting tests pass. The nine-mode instrumentation check is excluded from formal timings. Three preparation failures stopped before the first formal sample (an import error, an inappropriate comparison to Open Babel's time-bearing program header, and a shadowed observer variable). All original drivers, logs, seals and unmeasured outputs are retained. A standalone ten-profile preflight then passed before any timing epoch. Open Babel SDF semantic prerequisites are recomputed independently from fresh output; its program-header timestamp is outside the unchanged declared molecular fields. Neither source byte-exact conditions nor performance thresholds/orders were changed.

Actual generated code calls a separate classifier from the 3D reader. The new reader owner has {proof['candidate_reader']['instructions']} instructions and the separate classifier {proof['candidate_classifier']['instructions']}; the prior inlined baseline owner had 3,881 instructions. Instruction counts explain source realization, not total performance or the cause of previous writer penalties. Full disassembly is retained.

## Prospective protocol and limits

The formatter-only and compound attempts in #984 and the degree-only C attempt in [its full report](2026-10-17-openbabel-degree-steady.md) failed their fixed gates. C retains its SDF write penalties and all negative observations. Its same-binary controls showed a nonunit cyclic write interval, while all balanced-control intervals included one. Those controls do not establish causation and do not erase C's failed result. This is a distinct outlined source candidate with a newly fixed protocol, not an unchanged-candidate rerun or pooled analysis.

Before any D timing: 24 blocks/profile; all six three-arm permutations used four times for the two SDF profiles, and twelve AB/twelve BA source orders for the remaining eight profiles, with one fixed shuffle per profile. Orders and position balance are audited. Source adoption requires both predeclared branched-3D parse CI lower bounds>1 and every one of 30 unrounded paired source medians>=0.99, plus exact outputs. No threshold exceptions, filtering, performance-informed reruns, retiming or reuse of previous elapsed samples.

The same ten complete corpus hashes total 4,384 records: two 1,000-row SDF slices, corresponding V3000 slices, 64 linear chains, continued variants and four arbitrary-coordinate 3D chain/F-branched variants spanning 16..1024 carbons. These synthetic coordinates are deterministic nonplanar inputs, not physical conformers or a representative holdout. Source-only comparisons cover all ten; fresh external Open Babel timings cover only the two SDF slices. There are no new external speed rankings for the other eight profiles, and old speed ratios are never multiplied.

Each mode receives >=500ms warmup, then >=100ms timed work and >=3 complete batches. Actual repeats can differ by engine; raw elapsed, repeats, processed records, output units and budgets remain in every row. Normalization is raw elapsed/actual repeats, rounded down by <1ns per batch. All 1,584 samples are independently audited. Parse/roundtrip create fresh objects and retain the full batch; prepared write reuses its batch. SDF Rust reuses a String; V3000 and Open Babel allocate. Input loading and startup/plugin discovery are excluded; SDF framing is inside timing and V3000 framing outside. Whole-process RSS includes preparation/warmup/batch retention, is not isolated-operation or streaming memory, and is not ranked. Host load is recorded but external activity and CPU frequency are uncontrolled.

The baseline reuses exact C sealed binaries, actual empty-target production compiler logs and source proof at current main `{base}`; no C timing samples are used. Baseline source and every copied binary still match that seal. The candidate uses a newly empty target with opt3/full LTO/one CGU/one build job/incremental=false/strip=none. Both production libraries are recompiled; seven common native probe sources are byte-identical in the two arms. Compiler versions are checked against C; package-local Open Babel files, sources, protocol, inputs and binaries are sealed before/after profiles. System/package-external runtime data are outside the seal.

Open Babel `{env['openbabel'].strip()}`, native C++ -O3; observer RDKit `{env['rdkit']}`. External parse ranking requires accepted equivalent canonical graphs. Write/roundtrip additionally require declared atom/bond fields, coordinates to 1e-4, names/properties and explicit unknown stereo. Strict external wins require CI lower>1 and 24/24 faster blocks. Paired bootstrap uses 10,000 resamples, seed0 source/seed1 Open Babel; intervals display three decimals, adoption uses unrounded source medians.

## Complete matrix

Milliseconds per whole unfiltered batch, not process duration. Paired medians need not equal ratios of marginal medians.

| Profile | Mode | Base ms | Candidate ms | Base/candidate [95% CI]; wins | Fresh Open Babel/candidate [95% CI]; wins |
|---|---|---:|---:|---|---|
'''+ '\n'.join(rows)+'\n\n## Decision and all negative observations\n\n```json\n'+json.dumps(decision,indent=2)+'\n```\n\n## Evidence and reproduction\n\n[Hashed evidence packet]('+NAME+'-artifacts/README.md) contains every raw native sample, counts/budgets, input and output streams, exact source patch/common instrumentation, reused baseline provenance, fresh candidate code, native test/compiler logs, full disassembly and the independently recomputed audit. Full local evidence: `'+str(O)+'`. Rebuild baseline main with the retained common probe patch and the measured candidate in independent targets with recorded flags; run `build_native.py`, `instrument.py`, `seal.py`, `formal.py`, `assemble.py` after adapting original paths. The fresh process after native validation must be free of own competing builds/profiling. Retain all failures and do not pool previous elapsed data. Hosted CI completion, integration revisions and main merge are separate from this measured source.\n'
(R/'benchmarks'/(NAME+'.md')).write_text(report)
index=R/'benchmarks/README.md';text=index.read_text();marker='## What counts as a speed win';entry='- [Outlined classifier with balanced native timing]('+NAME+'.md): 24-block source retention across 30 lanes; fresh SDF-only Open Babel comparison.\n\n';assert marker in text
if entry not in text:index.write_text(text.replace(marker,entry+marker))
for name,h in packet.items():assert sha((D/name).read_bytes())==h
print('Report and',len(packet),'hashed evidence files assembled and verified',flush=True)
