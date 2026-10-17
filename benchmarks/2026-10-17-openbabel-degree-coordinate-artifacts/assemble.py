from pathlib import Path
import gzip,hashlib,json,statistics,shutil,subprocess,sys
O=Path(__file__).resolve().parent;R=Path('/private/tmp/chematic-openbabel-writer-20261017')
NAME='2026-10-17-openbabel-degree-coordinate';D=R/'benchmarks'/(NAME+'-artifacts');D.mkdir(exist_ok=True)
plan=json.loads((O/'measurement-plan.json').read_text());decision=json.loads((O/'decision.json').read_text());validation=json.loads((O/'validation.json').read_text());manifest=json.loads((O/'source-manifest.json').read_text())
base=manifest['base_commit'];candidate=manifest['candidate_commit']
src=subprocess.check_output(['git','show',base+':crates/chematic-mol/src/mol2000/writer.rs'],cwd=R,text=True)
old=src[src.index('pub(crate) fn push_fixed4'):src.index('\npub(crate) fn push_right_aligned_u32')].strip()
test=(R/'crates/chematic-mol/src/mol2000/fixed4_tests.rs').read_text();frozen=test[test.index('fn frozen_push_fixed4'):test.index('\nfn check(')].replace('fn frozen_push_fixed4','pub(crate) fn push_fixed4',1).strip();assert old==frozen
(O/'frozen-function-proof.json').write_text(json.dumps({'base_commit':base,'complete_frozen_function_exact':True,'trimmed_function_sha256':hashlib.sha256(old.encode()).hexdigest(),'tests':'Every fractional pair 0000..9999, 100000 fixed-seed random f64 bit patterns, signs, padding widths 0/1/6/10/24/30, decimal carries, neighboring f64 values around rounding ties, zero/subnormal/large/nonfinite boundaries; both frozen and standard-format oracles; existing 400000 generated coordinates plus directed formatter cases retained'},indent=2)+'\n')
classifier=subprocess.check_output(['git','show',base+':crates/chematic-perception/src/stereo_validation.rs'],cwd=R,text=True)
a=classifier.index('pub fn stereo_centers(');b=classifier.index('\n/// Summarise',a);old_classifier=classifier[a:b].strip()
t=(R/'crates/chematic-perception/src/stereo_degree_tests.rs').read_text();a=t.index('fn frozen_stereo_centers(');b=t.index('\nfn compare(',a)
assert old_classifier==t[a:b].replace('fn frozen_stereo_centers(','pub fn stereo_centers(',1).strip()
(O/'frozen-classifier-proof.json').write_text(json.dumps({'base_commit':base,'complete_frozen_classifier_exact':True,'body_sha256':hashlib.sha256(old_classifier.encode()).hexdigest(),'pinned_graphs':14072,'five_atom_graph_hydrogen_variants':6144,'proof':'Every emitted branch needs heavy_nbs length 3 or 4, including constrained ring nitrogen; heavy_nbs.len <= degree. find_sssr remains before the loop.'},indent=2)+'\n')
patch=subprocess.check_output(['git','diff',base,candidate,'--','crates'],cwd=R);(O/'candidate.patch').write_bytes(patch)
sys.path.insert(0,str(R/'validation/cosmolkit_comparison'))
from bench_corpus import paired_median_ratio_ci
rows=[];audit=[];failures=[];negative=[];primary={}
for p in plan['profiles']:
 key=p['kind']+'-'+p['profile'];d=json.loads((O/('paired-'+key+'.json')).read_text())
 assert d['complete'] and len(d['blocks'])==21 and d['source_parse_snapshots_equal'] and d['rewrite_sha256']['baseline']==d['rewrite_sha256']['candidate']
 for mode in ['parse','write','roundtrip']:
  values={arm:[x[arm][mode]['elapsed_ns'] for x in d['blocks']] for arm in ['baseline','candidate','openbabel']}
  assert all(len(v)==21 and all(x>0 for x in v) for v in values.values())
  source=paired_median_ratio_ci(values['baseline'],values['candidate'],seed=0);assert source==decision['source_ratios'][key+'/'+mode]
  raw=statistics.median(x/y for x,y in zip(values['baseline'],values['candidate'],strict=True))
  if raw<.99:failures.append(key+'/'+mode)
  if source['ci95'][1]<1:negative.append(key+'/'+mode)
  if key+'/'+mode in plan['primary']:primary[key+'/'+mode]=source
  if d['summary'][mode]['baseline']['speedup'] is not None:assert source==d['summary'][mode]['baseline']['speedup']
  ob=d['summary'][mode]['openbabel'];speed=ob['speedup']
  if speed is not None:assert speed==paired_median_ratio_ci(values['openbabel'],values['candidate'],seed=1)
  def fmt(x):return f"{x['median']:.3f} [{x['ci95'][0]:.3f}, {x['ci95'][1]:.3f}]; {x['winning_blocks']}/21"
  rows.append('| '+key+' | '+mode+' | '+f"{statistics.median(values['baseline'])/3e6:.3f} | {statistics.median(values['candidate'])/3e6:.3f} | "+fmt(source)+' | '+(fmt(speed)+(' strict win' if ob['strict_speed_win'] else '') if speed else 'unranked: semantic prerequisite fails')+' |')
  audit.append({'lane':key+'/'+mode,'source_ratio_recomputed':True,'source_unrounded_median':raw,'openbabel_ratio_recomputed':speed is not None})
assert len(audit)==30
assert failures==[x['lane'] for x in decision['source_floor_failures']]
assert negative==[x['lane'] for x in decision['negative_intervals']]
assert decision['local_adoption_gate_passed']==(len(primary)==4 and all(x['ci95'][0]>1 for x in primary.values()) and not failures)
(O/'audit.json').write_text(json.dumps({'passed':True,'all_30_source_ratios_recomputed':True,'all_ranked_comparator_ratios_recomputed':True,'adoption_recomputed_from_unrounded_source_medians':True,'no_filtered_rows_or_timing_reruns':True,'complete_frozen_function_exact':True,'comparisons':audit},indent=2)+'\n')
for name in ['measurement-plan.json','decision.json','source-proof.json','baseline-reuse.json','frozen-classifier-proof.json','frozen-function-proof.json','binary-source-proof.json','source-manifest.json','binary-seal.json','environment.json','validation.json','native-check-commands.json','preflight-source-equality.json','audit.json','diagnostic-profile-scope.json','baseline-prerequisite-command.json']:
 shutil.copy2(O/name,D/name)
for name in ['setup.py','build.py','seal.py','formal.py','native.py','profile.py','assemble.py','run.py']:shutil.copy2(O/name,D/name)
for name in ['baseline-build.log','candidate-build.log','source-commit.log','native-tests.log','clippy.log','format.log','refactor.log','runner-tests.log','diagnostic-profile.txt','diagnostic-profile-worker.json','profile-sampler-output.txt','baseline-fixed4-disassembly.txt','candidate-fixed4-disassembly.txt','baseline-stereo-owner-disassembly.txt','candidate-stereo-owner-disassembly.txt','candidate.patch','formal.log','baseline-prerequisite-tests.log','driver.log','build-driver.log','native-driver.log','profile.log','seal.log']:
 (D/(name+'.gz')).write_bytes(gzip.compress((O/name).read_bytes(),mtime=0))
for p in plan['profiles']:
 key=p['kind']+'-'+p['profile']
 for suffix in ['json','log']:
  name='paired-'+key+'.'+suffix;(D/(name+'.gz')).write_bytes(gzip.compress((O/name).read_bytes(),mtime=0))
 (D/(p['path']+'.gz')).write_bytes(gzip.compress((O/p['path']).read_bytes(),mtime=0))
outcomes=D/'outcomes';outcomes.mkdir(exist_ok=True);outcome_manifest={}
for profile in plan['profiles']:
 key=profile['kind']+'-'+profile['profile'];matrix=json.loads((O/('paired-'+key+'.json')).read_text())
 for mode in ['snapshot','dump']:
  left=O/f'preflight-{key}-baseline-{mode}.gz';right=O/f'preflight-{key}-candidate-{mode}.gz'
  raw=gzip.decompress(left.read_bytes());assert raw==gzip.decompress(right.read_bytes())
  digest=hashlib.sha256(raw).hexdigest()
  check=matrix['snapshot_sha256'] if mode=='snapshot' else matrix['rewrite_sha256']
  assert check['baseline']==check['candidate']==digest
  name=key+'-source-common-'+mode+'.gz';shutil.copy2(left,outcomes/name)
  outcome_manifest[name]={'engines':['baseline','candidate'],'source_streams_identical':True,'records':profile['records'],'raw_sha256':digest}
 for path in sorted((O/('paired-'+key+'.artifacts')).glob('openbabel.*')):
  raw=path.read_bytes();digest=hashlib.sha256(raw).hexdigest()
  if path.suffix in ['.sdf','.molblocks']:assert digest==matrix['rewrite_sha256']['openbabel']
  name=key+'-'+path.name+'.gz';(outcomes/name).write_bytes(gzip.compress(raw,mtime=0))
  outcome_manifest[name]={'engine':'openbabel','records':profile['records'],'raw_sha256':digest,'artifact':path.name}
(D/'outcome-streams.json').write_text(json.dumps(outcome_manifest,indent=2)+'\n')
(D/'README.md').write_text('One-host source experiment. The fixed plan, every raw timing block, exact inputs, source/native/executable seals, production disassembly, verbose compiler logs, full frozen-function proof and recomputed statistical audit are retained. The outcomes directory contains shared complete source snapshots/rewrites and every retained Open Babel rewrite/graph artifact, including unranked rows. Full executables and duplicate per-arm streams remain at '+str(O)+'. Recipes retain original local absolute paths and require adaptation. No universal or published-package performance claim.\n')
packet={str(p.relative_to(D)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(D.rglob('*')) if p.is_file() and p.name!='packet.sha256.json'}
(D/'packet.sha256.json').write_text(json.dumps(packet,indent=2)+'\n')
status='PASS' if decision['local_adoption_gate_passed'] else 'FAIL: compound candidate not adopted'
report=f"""# Open Babel comparison: coordinate pairs and impossible-center degree filtering

Decision: **{status}**. Both diverse SDF write lanes and both branched-3D parse lanes must improve with paired CI lower bounds above 1. Every one of the 30 unrounded source median ratios must be >=0.99 and all source outputs must match. This is one macOS source experiment, not universal or published-package superiority.

## Change and correctness

Baseline `{base}`; measured compound candidate `{candidate}`. Both contain the identical independent ring-partner repair #986 over main `f20304f4`. The formatter converts decimal digits in pairs with a 200-byte ASCII table and combines common padding into the same stack buffer. Rounding, negative zero and general-format fallbacks are retained.

The classifier now rejects graph degree <3 before allocating heavy-neighbor vectors, calculating H counts or traversing branch signatures. Every former return shape requires three or four heavy neighbors, including constrained ring nitrogen, and heavy-neighbor count cannot exceed graph degree. The original global `find_sssr` preparation still runs first; all valid-center paths and the final output-shape check are unchanged. No public API, dependency or unsafe additions. The complete frozen old classifier matches the base source exactly after declaration renaming and is compared on 14,072 pinned graphs, all 6,144 five-vertex/H-override graphs and directed boundaries. The frozen numeric formatter and Rust standard formatting remain independent byte oracles for every fractional pair, 100,000 deterministic random f64 bits, padding/sign/rounding boundaries and the pre-existing 400,000 generated coordinate cases.

All 4,384 source native snapshots and rewritten records are byte-exact. Native validation: {validation['native']['passed']} passed, {validation['native']['failed']} failed, {validation['native']['ignored']} existing ignored across {validation['native']['suites']} suites covering perception, MOL, chemistry and fingerprints. Release all-target Clippy, formatting, unchanged refactor budget and ten runner tests pass. A separate baseline branched-3D parse CPU trace retains classifier/branch-traversal ancestry; its worker elapsed value is excluded from all formal samples.

## Separate candidate and fixed protocol

The [formatter-only experiment](2026-10-17-openbabel-coordinate-pairs.md) improved SDF writes but failed its retention gate at branched-spiral 3D parse, median 0.983815. This compound candidate adds a necessary precondition to remove that workload's known impossible-center traversal. It has a new fixed plan and a separate timing epoch. No old elapsed observations are reused, pooled or filtered, and no retention threshold is relaxed. The previous candidate and its initial pre-measurement main-regression failure remain fully recorded in that earlier packet.

Ten unchanged input profiles/hashes: two unfiltered 1,000-row diverse SDF slices; corresponding V3000 slices; 64 linear chains; continued lines; two synthetic 3D chain profiles and two fluorine-branched versions. Synthetic chains span 16..1024 carbons and use arbitrary deterministic finite nonplanar coordinates. They are not physical conformers, embedding-quality evidence or a representative holdout. Four predeclared primary lanes and all 30 retention lanes are recorded below, including every negative interval and unranked arm.

Twenty-one rotating three-arm fresh-process blocks per profile, three native batch repeats, parse/prepared write/roundtrip. Native startup, input loading and plugin discovery are excluded. V3000 framing is outside the timer and SDF framing inside. Diagnostic readers retain full batches. Prepared writing has an untimed warm write; Rust SDF reuses a String while V3000 and Open Babel allocate. Whole-process high-water RSS includes preparation, retained batches and warmup, not isolated-operation or streaming memory. External host load is uncontrolled and logged at each block start.

The exact independent baseline binaries, source tree and verbose build logs are reused from the sealed formatter-only experiment; every source file still matches Git `{base}` and every binary hash matches that prior seal. The candidate uses a new initially empty independent target. Rust compiler identity is verified unchanged; opt3/full LTO/one CGU/one build job/incremental=false/strip=none. Both production libraries are verified in the verbose logs. Actual coordinate-formatter bodies and the 3D perception owner containing the classifier have different machine code in the arms. `baseline-reuse.json` records the original baseline build command and source root; all timing blocks are new. Native comparator: {json.loads((O/'environment.json').read_text())['openbabel'].strip()}, C++ probe -O3; independent observer RDKit {json.loads((O/'environment.json').read_text())['rdkit']}. Seals cover tracked source, baseline files, corpora, binaries and package-local Open Babel files before/after profiles. System/package-external runtime data remain outside the seal.

Ratios are reference/candidate time, deterministic paired 10,000-resample bootstrap, seed 0 for source and 1 for Open Babel, reported to three decimals. Retention uses unrounded medians. Parse ranking requires accepted equivalent canonical graphs; write/roundtrip ranking also requires declared atom/bond fields, coordinates to 1e-4, names/properties and explicit unknown stereo. Mismatches/rejections remain unranked. Strict external wins additionally require 21/21 faster blocks. Ratios from the predecessor or older experiments are never multiplied into these measurements.

## Complete matrix

Milliseconds per full batch: native elapsed divided by three. Paired medians need not equal ratios of displayed marginal medians.

| Profile | Mode | Base ms | Candidate ms | Base/candidate [95% CI]; wins | Open Babel/candidate [95% CI]; wins |
|---|---|---:|---:|---|---|
"""+'\n'.join(rows)+'\n\n## Decision and every negative result\n\n```json\n'+json.dumps(decision,indent=2)+'\n```\n\n## Evidence and reproduction\n\n[Complete packet]('+NAME+'-artifacts/README.md) contains every raw matrix, exact input, complete common source output streams, retained Open Babel graph/rewrite outcomes, source/binary seals, compiler logs, disassembly, frozen-function proofs, full test logs and recomputed statistical audit. Original absolute paths in recipes require adaptation. Rebuild the named baseline and candidate in independent roots/targets using the original command in `baseline-reuse.json` and `build.py`; run `native.py`, `profile.py`, `seal.py`, then `formal.py`. `assemble.py` independently checks all source/ranked comparator ratios, the raw median floors, complete frozen functions and outcome hashes. Full local evidence: `'+str(O)+'`. A later integration revision requires correctness validation and is not the timed candidate.\n'
(R/'benchmarks'/(NAME+'.md')).write_text(report)
index=R/'benchmarks/README.md';s=index.read_text();marker='## What counts as a speed win';assert marker in s
entry='- [Coordinate pairs and degree filtering]('+NAME+'.md): separate compound-source experiment with four primary lanes, all 30 retention results and complete outcome streams.\n\n'
if entry not in s:s=s.replace(marker,entry+marker);index.write_text(s)
print('Report and',len(packet),'hashed files assembled',flush=True)
