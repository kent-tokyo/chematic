from pathlib import Path
import gzip,hashlib,json,statistics,shutil,subprocess,sys
O=Path(__file__).resolve().parent;R=Path('/private/tmp/chematic-degree-only-steady-20261017');NAME='2026-10-17-openbabel-degree-steady';D=R/'benchmarks'/(NAME+'-artifacts');D.mkdir(exist_ok=True)
def sha(data):return hashlib.sha256(data).hexdigest()
plan=json.loads((O/'measurement-plan.json').read_text());decision=json.loads((O/'decision.json').read_text());validation=json.loads((O/'validation.json').read_text());manifest=json.loads((O/'source-manifest.json').read_text());base=manifest['base_commit'];candidate=manifest['candidate_commit']
src=subprocess.check_output(['git','show',base+':crates/chematic-perception/src/stereo_validation.rs'],cwd=R,text=True);a=src.index('pub fn stereo_centers(');b=src.index('\n/// Summarise',a);old=src[a:b].strip()
test=(R/'crates/chematic-perception/src/stereo_degree_tests.rs').read_text();a=test.index('fn frozen_stereo_centers(');b=test.index('\nfn compare(',a);assert old==test[a:b].replace('fn frozen_stereo_centers(','pub fn stereo_centers(',1).strip()
(O/'frozen-classifier-proof.json').write_text(json.dumps({'base_commit':base,'complete_frozen_classifier_exact':True,'body_sha256':sha(old.encode()),'pinned_graphs':14072,'five_atom_graph_hydrogen_variants':6144,'necessary_condition':'Every return shape requires heavy_nbs length 3 or 4; heavy_nbs.len <= graph degree. find_sssr preparation remains before the loop.'},indent=2)+'\n')
patch=subprocess.check_output(['git','diff',base,candidate,'--','crates','tools','scripts'],cwd=R);(O/'candidate.patch').write_bytes(patch)
common=json.loads((O/'source-proof.json').read_text())['common_instrumentation_sha256'];common_patch=subprocess.check_output(['git','diff',base,candidate,'--',*common],cwd=R);(O/'common-harness.patch').write_bytes(common_patch)
sys.path.insert(0,str(R/'validation/cosmolkit_comparison'));from bench_corpus import paired_median_ratio_ci
rows=[];audit=[];failures=[];negative=[];primary={};native_samples=0
for p in plan['profiles']:
 key=p['kind']+'-'+p['profile'];data=json.loads((O/('paired-'+key+'.json')).read_text());assert data['complete'] and data['schema']=='openbabel-steady-throughput/v1' and len(data['blocks'])==21 and data['source_parse_snapshots_equal'] and data['rewrite_sha256']['baseline']==data['rewrite_sha256']['candidate']
 for block in data['blocks']:
  for arm in ['baseline','candidate','openbabel']:
   for mode in ['parse','write','roundtrip']:
    row=block[arm][mode];n=row['repeats'];raw=row['elapsed_ns_raw']
    assert type(n) is int and n>=3 and type(raw) is int and raw>=100000000
    assert row['records']==n*p['records'] and row['records_per_batch']==p['records']
    assert row['warmup_repeats']>=1 and row['warmup_elapsed_ns']>=500000000
    assert row['minimum_warmup_ns']==500000000 and row['minimum_timed_ns']==100000000
    assert row['elapsed_ns']==row['elapsed_ns_per_batch']==raw//n>0
    native_samples+=1
 for mode in ['parse','write','roundtrip']:
  values={arm:[x[arm][mode]['elapsed_ns'] for x in data['blocks']] for arm in ['baseline','candidate','openbabel']}
  source=paired_median_ratio_ci(values['baseline'],values['candidate'],seed=0);assert source==decision['source_ratios'][key+'/'+mode]
  raw_median=statistics.median(x/y for x,y in zip(values['baseline'],values['candidate'],strict=True))
  if raw_median<.99:failures.append(key+'/'+mode)
  if source['ci95'][1]<1:negative.append(key+'/'+mode)
  if key+'/'+mode in plan['primary']:primary[key+'/'+mode]=source
  reported=data['summary'][mode]['baseline']['speedup']
  if reported is not None:assert reported==source
  ob=data['summary'][mode]['openbabel'];speed=ob['speedup']
  if speed is not None:assert speed==paired_median_ratio_ci(values['openbabel'],values['candidate'],seed=1)
  def fmt(x):return f"{x['median']:.3f} [{x['ci95'][0]:.3f}, {x['ci95'][1]:.3f}]; {x['winning_blocks']}/21"
  rows.append('| '+key+' | '+mode+' | '+f"{statistics.median(values['baseline'])/1e6:.3f} | {statistics.median(values['candidate'])/1e6:.3f} | "+fmt(source)+' | '+(fmt(speed)+(' strict win' if ob['strict_speed_win'] else '') if speed else 'unranked: semantic prerequisite fails')+' |')
  audit.append({'lane':key+'/'+mode,'source_unrounded_median':raw_median,'source_ratio_recomputed':True,'openbabel_ratio_recomputed':speed is not None})
assert native_samples==1890 and len(audit)==30
assert failures==[x['lane'] for x in decision['source_floor_failures']];assert negative==[x['lane'] for x in decision['negative_intervals']]
assert decision['local_adoption_gate_passed']==(len(primary)==2 and all(x['ci95'][0]>1 for x in primary.values()) and not failures)
(O/'audit.json').write_text(json.dumps({'passed':True,'all_1890_native_sample_counts_budgets_and_normalizations_checked':True,'all_30_source_ratios_recomputed':True,'all_ranked_comparator_ratios_recomputed':True,'adoption_recomputed_from_unrounded_medians':True,'no_filtered_rows_or_reruns':True,'comparisons':audit},indent=2)+'\n')
for name in ['measurement-plan.json','prospective-plan-original.json','source-proof.json','decision.json','frozen-classifier-proof.json','binary-source-proof.json','source-manifest.json','binary-seal.json','environment.json','validation.json','native-check-commands.json','preflight-source-equality.json','instrumentation-check.json','audit.json','diagnostic-profile-scope.json']:
 shutil.copy2(O/name,D/name)
for name in ['setup.py','patch_harness.py','build.py','seal.py','formal.py','native.py','profile.py','instrument.py','assemble.py','run.py']:
 shutil.copy2(O/name,D/name)
for name in ['baseline-build.log','candidate-build.log','openbabel-sdf-build.log','openbabel-v3000-build.log','source-commit.log','native-tests.log','clippy.log','format.log','refactor.log','runner-tests.log','instrument.log','profile.log','seal.log','driver.log','diagnostic-profile.txt','diagnostic-profile-worker.json','baseline-stereo-owner-disassembly.txt','candidate-stereo-owner-disassembly.txt','candidate.patch','common-harness.patch','formal.log']:
 (D/(name+'.gz')).write_bytes(gzip.compress((O/name).read_bytes(),mtime=0))
for name,h in common.items():
 assert sha((R/name).read_bytes())==h==sha((O/'baseline-source'/name).read_bytes())
 path=D/'common-instrumentation'/name;path.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(R/name,path)
for p in plan['profiles']:
 key=p['kind']+'-'+p['profile']
 for suffix in ['json','log']:
  name='paired-'+key+'.'+suffix;(D/(name+'.gz')).write_bytes(gzip.compress((O/name).read_bytes(),mtime=0))
 (D/(p['path']+'.gz')).write_bytes(gzip.compress((O/p['path']).read_bytes(),mtime=0))
outcomes=D/'outcomes';outcomes.mkdir(exist_ok=True);streams={}
for p in plan['profiles']:
 key=p['kind']+'-'+p['profile'];data=json.loads((O/('paired-'+key+'.json')).read_text())
 for mode in ['snapshot','dump']:
  left=O/f'preflight-{key}-baseline-{mode}.gz';right=O/f'preflight-{key}-candidate-{mode}.gz';raw=gzip.decompress(left.read_bytes());assert raw==gzip.decompress(right.read_bytes())
  hashes=data['snapshot_sha256'] if mode=='snapshot' else data['rewrite_sha256'];assert hashes['baseline']==hashes['candidate']==sha(raw)
  name=key+'-source-common-'+mode+'.gz';shutil.copy2(left,outcomes/name);streams[name]={'engines':['baseline','candidate'],'source_streams_identical':True,'records':p['records'],'raw_sha256':sha(raw)}
 for path in sorted((O/('paired-'+key+'.artifacts')).glob('openbabel.*')):
  raw=path.read_bytes()
  if path.suffix in ['.sdf','.molblocks']:assert sha(raw)==data['rewrite_sha256']['openbabel']
  name=key+'-'+path.name+'.gz';(outcomes/name).write_bytes(gzip.compress(raw,mtime=0));streams[name]={'engine':'openbabel','expected_records':p['records'],'raw_sha256':sha(raw),'artifact':path.name}
(D/'outcome-streams.json').write_text(json.dumps(streams,indent=2)+'\n')
(D/'README.md').write_text('Fixed steady-state source experiment on one macOS host. Every raw native sample retains actual repeats, processed records, raw elapsed and warmup/budget fields; elapsed_ns is raw elapsed divided by actual repeats, in ns per full unfiltered batch. All raw outcomes and unranked Open Babel artifacts are included. Common instrumentation is byte-identical in both arms, and the complete classifier is frozen. Prior fixed-three-repeat observations are not reused. Scripts retain original absolute paths and require adaptation; full binaries remain at '+str(O)+'.\n')
packet={str(p.relative_to(D)):sha(p.read_bytes()) for p in sorted(D.rglob('*')) if p.is_file() and p.name!='packet.sha256.json'};(D/'packet.sha256.json').write_text(json.dumps(packet,indent=2)+'\n')
status='PASS' if decision['local_adoption_gate_passed'] else 'FAIL: degree-only candidate not adopted'
report=f'''# Open Babel comparison: degree-only stereocenter filtering with steady native batches

Decision: **{status}**. Both branched-3D parse CI lower bounds must exceed 1; every one of 30 unrounded source median ratios must be >=0.99; source outcomes must be exact. This is source-build evidence on one macOS host, not universal or published-package superiority.

## Change and correctness

Baseline `{base}`; measured candidate `{candidate}`. The production change is only the degree>=3 precondition in `stereo_centers`. Every old return shape requires three or four heavy neighbors, including constrained ring nitrogen, and heavy-neighbor count cannot exceed graph degree. Rejecting lower-degree atoms before Vec allocation, H counting and rooted branch traversal preserves every valid-center path and the final shape check. Global `find_sssr` cache preparation remains before the loop. Numerical formatting is byte-identical to the baseline source. No public API, dependency or production unsafe additions.

The complete old classifier is byte-exact after only renaming its declaration; 14,072 parsed graphs, all 6,144 five-vertex/H-override variants and directed boundaries match it. All 4,384 native source snapshots and warmed rewrites match byte for byte. Native validation: {validation['native']['passed']} passed, {validation['native']['failed']} failed, {validation['native']['ignored']} existing ignored across {validation['native']['suites']} suites covering perception, MOL, chemistry and fingerprints. Release all-target Clippy, formatting and the unchanged refactor ratchet pass; all {validation['runner_log'].split('passed')[0].split()[-1]} runner/accounting tests pass. An additional nine-mode native instrumentation preflight verifies raw counts, elapsed time and budgets in all three engines; those observations and the baseline CPU diagnostic are excluded from every formal sample.

## New fixed throughput protocol

The rejected formatter and compound candidates remain in PR #984, with their complete observations. This candidate removes numeric changes and isolates the necessary classifier condition. Before its timing began, the source baseline was updated to latest main including #981/#986. The original prospective plan is retained. No earlier elapsed samples are reused, pooled or used as new observations. The two parser primary lanes and the 0.99 retention threshold remain fixed.

Ten unchanged corpus profiles/hashes: two unfiltered 1,000-row diverse SDF slices, corresponding V3000 slices, 64 linear chains, continued versions, two synthetic 3D chain profiles and two F-branched variants; total 4,384 records. Synthetic chains span 16..1024 carbons with arbitrary deterministic finite nonplanar coordinates. They are not physical conformers, embedding-quality evidence or a representative holdout.

Each profile has 21 rotating three-arm fresh-process blocks, with parse, prepared write and roundtrip. All engines first run the measured mode for at least 500ms, then time at least 100ms and at least three whole batches. Faster engines may process more complete batches. Every sample records raw elapsed, actual iterations, actual records and units; ns per batch is raw elapsed divided by actual iterations, rounded down by less than 1ns. Paired ratios use those per-batch values. All 1,890 sample counts, budgets and normalizations are independently rechecked. This is steady-state throughput, not cold-process or startup latency. Absolute times and whole-process RSS are not directly comparable to prior fixed-three-repeat experiments.

Parse and roundtrip create fresh result objects on every iteration; prepared write reuses its batch. Rust SDF reuses a String; V3000 and Open Babel allocate. Rewrites are checked after the same write warmup/time budget, including all unranked rows. Input loading and native startup/plugin discovery are excluded; V3000 framing is outside timing and SDF framing inside. Whole-process RSS includes preparation, warmup and batch retention; it is not isolated-operation/streaming memory and is not ranked. External host load remains uncontrolled and is logged at each block start.

Both initially empty independent targets use opt3/full LTO/one CGU/one build job/incremental=false/strip=none. All baseline files match Git except seven explicitly listed probe-only instrumentation files, whose bytes match candidate files exactly. Both real production compilations are verified in verbose logs; the actual 3D owner containing the classifier has different machine code. The instrumentation patch and exact sources are retained. Comparator {json.loads((O/'environment.json').read_text())['openbabel'].strip()}, native C++ probe -O3; observer RDKit {json.loads((O/'environment.json').read_text())['rdkit']}. Sources, baseline files, inputs, binaries, protocol and package-local Open Babel files are sealed before/after profiles; system/package-external runtime data are outside the seal.

Paired bootstrap: 10,000 resamples, seed0 for source, seed1 for Open Babel, intervals reported to three decimals. Raw source medians determine retention. External parse ranking requires accepted equivalent canonical graphs; write/roundtrip also require declared atom/bond fields, coordinates to 1e-4, names/properties and explicit unknown stereo. Rejections/mismatches are retained and unranked. Strict external wins additionally require 21/21 faster blocks. Earlier source ratios are never multiplied into these measurements.

## Complete matrix

Milliseconds per complete batch from actual repetition counts, not total native process duration. Paired medians need not equal ratios of displayed marginal medians.

| Profile | Mode | Base ms | Candidate ms | Base/candidate [95% CI]; wins | Open Babel/candidate [95% CI]; wins |
|---|---|---:|---:|---|---|
'''+ '\n'.join(rows)+'\n\n## Decision and every negative result\n\n```json\n'+json.dumps(decision,indent=2)+'\n```\n\n## Evidence and reproduction\n\n[Complete packet]('+NAME+'-artifacts/README.md) includes raw samples/counts/budgets, complete input and outcome streams, exact common probe patch, source/binary/native seals, compiler logs, disassembly, frozen-classifier proof, native test logs and recomputed statistical audit. Rebuild the two named source revisions in independent roots, apply the same retained common probe patch to the baseline, run `build.py`, `native.py`, `instrument.py`, `profile.py`, `seal.py`, then `formal.py`; adapt the original absolute paths. `assemble.py` verifies the raw accounting and all statistics. Full local evidence: `'+str(O)+'`. A later integration revision requires correctness validation and is not the timed candidate.\n'
(R/'benchmarks'/(NAME+'.md')).write_text(report)
index=R/'benchmarks/README.md';text=index.read_text();marker='## What counts as a speed win';assert marker in text
entry='- [Degree-only steady native throughput]('+NAME+'.md): parser-only production change, symmetric warmup/time budgets, raw counts and all 30 retention lanes.\n\n'
if entry not in text:index.write_text(text.replace(marker,entry+marker))
print('Report and',len(packet),'hashed evidence files assembled',flush=True)
