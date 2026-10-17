from pathlib import Path
import gzip,hashlib,json,statistics,shutil,subprocess,sys
O=Path(__file__).resolve().parent;R=Path('/private/tmp/chematic-openbabel-writer-20261017')
NAME='2026-10-17-openbabel-coordinate-pairs';D=R/'benchmarks'/(NAME+'-artifacts');D.mkdir(exist_ok=True)
plan=json.loads((O/'measurement-plan.json').read_text());decision=json.loads((O/'decision.json').read_text());validation=json.loads((O/'validation.json').read_text());manifest=json.loads((O/'source-manifest.json').read_text())
base=manifest['base_commit'];candidate=manifest['candidate_commit']
src=subprocess.check_output(['git','show',base+':crates/chematic-mol/src/mol2000/writer.rs'],cwd=R,text=True)
old=src[src.index('pub(crate) fn push_fixed4'):src.index('\npub(crate) fn push_right_aligned_u32')].strip()
test=(R/'crates/chematic-mol/src/mol2000/fixed4_tests.rs').read_text();frozen=test[test.index('fn frozen_push_fixed4'):test.index('\nfn check(')].replace('fn frozen_push_fixed4','pub(crate) fn push_fixed4',1).strip();assert old==frozen
(O/'frozen-function-proof.json').write_text(json.dumps({'base_commit':base,'complete_frozen_function_exact':True,'trimmed_function_sha256':hashlib.sha256(old.encode()).hexdigest(),'tests':'Every fractional pair 0000..9999, 100000 fixed-seed random f64 bit patterns, signs, padding widths 0/1/6/10/24/30, decimal carries, neighboring f64 values around rounding ties, zero/subnormal/large/nonfinite boundaries; both frozen and standard-format oracles; existing 400000 generated coordinates plus directed formatter cases retained'},indent=2)+'\n')
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
  if key in plan['primary'] and mode=='write':primary[key]=source
  if d['summary'][mode]['baseline']['speedup'] is not None:assert source==d['summary'][mode]['baseline']['speedup']
  ob=d['summary'][mode]['openbabel'];speed=ob['speedup']
  if speed is not None:assert speed==paired_median_ratio_ci(values['openbabel'],values['candidate'],seed=1)
  def fmt(x):return f"{x['median']:.3f} [{x['ci95'][0]:.3f}, {x['ci95'][1]:.3f}]; {x['winning_blocks']}/21"
  rows.append('| '+key+' | '+mode+' | '+f"{statistics.median(values['baseline'])/3e6:.3f} | {statistics.median(values['candidate'])/3e6:.3f} | "+fmt(source)+' | '+(fmt(speed)+(' strict win' if ob['strict_speed_win'] else '') if speed else 'unranked: semantic prerequisite fails')+' |')
  audit.append({'lane':key+'/'+mode,'source_ratio_recomputed':True,'source_unrounded_median':raw,'openbabel_ratio_recomputed':speed is not None})
assert len(audit)==30
assert failures==[x['lane'] for x in decision['source_floor_failures']]
assert negative==[x['lane'] for x in decision['negative_intervals']]
assert decision['local_adoption_gate_passed']==(len(primary)==2 and all(x['ci95'][0]>1 for x in primary.values()) and not failures)
(O/'audit.json').write_text(json.dumps({'passed':True,'all_30_source_ratios_recomputed':True,'all_ranked_comparator_ratios_recomputed':True,'adoption_recomputed_from_unrounded_source_medians':True,'no_filtered_rows_or_timing_reruns':True,'complete_frozen_function_exact':True,'comparisons':audit},indent=2)+'\n')
for name in ['measurement-plan.json','decision.json','source-proof.json','source-copy-correction.json','frozen-function-proof.json','binary-source-proof.json','source-manifest.json','binary-seal.json','environment.json','validation.json','native-check-commands.json','preflight-source-equality.json','audit.json','diagnostic-profile-scope.json','baseline-prerequisite-command.json']:
 shutil.copy2(O/name,D/name)
for name in ['setup.py','prepare.py','build.py','seal.py','formal.py','native.py','profile.py','assemble.py','run_corrected.py']:shutil.copy2(O/name,D/name)
for name in ['baseline-build.log','candidate-build.log','source-commit.log','native-tests.log','clippy.log','format.log','refactor.log','runner-tests.log','diagnostic-profile.txt','diagnostic-profile-worker.json','profile-sampler-output.txt','baseline-fixed4-disassembly.txt','candidate-fixed4-disassembly.txt','candidate.patch','formal.log','prerequisite-targeted-tests.log','prerequisite-commit.log','baseline-prerequisite-tests.log','corrected-driver.log','build-driver.log','native-driver.log','profile.log','seal.log']:
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
A=D/'unmeasured-main-failure';A.mkdir(exist_ok=True)
for path in (O/'unmeasured-main-failure').iterdir():
 if path.name.startswith(('baseline-sdf','baseline-v3000','candidate-sdf','candidate-v3000')):continue
 if path.suffix in ['.log','.txt']:(A/(path.name+'.gz')).write_bytes(gzip.compress(path.read_bytes(),mtime=0))
 else:shutil.copy2(path,A/path.name)
(D/'README.md').write_text('One-host source experiment. The fixed plan, every raw timing block, exact inputs, source/native/executable seals, production disassembly, verbose compiler logs, full frozen-function proof and recomputed statistical audit are retained. The outcomes directory contains shared complete source snapshots/rewrites and every retained Open Babel rewrite/graph artifact, including unranked rows. Full executables and duplicate per-arm streams remain at '+str(O)+'. Recipes retain original local absolute paths and require adaptation. No universal or published-package performance claim.\n')
packet={str(p.relative_to(D)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(D.rglob('*')) if p.is_file() and p.name!='packet.sha256.json'}
(D/'packet.sha256.json').write_text(json.dumps(packet,indent=2)+'\n')
status='PASS' if decision['local_adoption_gate_passed'] else 'FAIL: formatter-only candidate was not adopted'
report=f'''# Open Babel comparison: decimal-pair coordinate formatting

Decision: **{status}**. The two diverse SDF write lanes are the predeclared primary. All source outputs must match and all 30 source median ratios must be at least 0.99. This is one macOS source experiment, not universal superiority or published-package evidence.

## Change and correctness

Measured baseline `{base}` (the independent ring-partner repair #986 over original main `f20304f4`); measured production candidate `{candidate}`. Both arms contain the identical correctness repair, so their production-source difference is the coordinate formatter. The shared fixed-four-decimal coordinate formatter converts integer and fractional digits two at a time with a 200-byte ASCII pair table. Common field padding is included in the same stack buffer and appended once. The rounding arithmetic, negative sign including negative zero, and standard-format fallback for ties, large and nonfinite values are unchanged. The helper moves into its own module to keep the existing writer within its refactor budget. No public API, dependency or unsafe additions. A separate baseline own-CPU profile observes prepared SDF write mode; its elapsed value is excluded from all formal results and preparation ancestry must not be attributed to the coordinate formatter.

The complete frozen old function is exact after renaming its declaration. Both it and Rust standard formatting serve as independent byte oracles for every fractional pair, 100,000 deterministic random f64 bit patterns and directed rounding/sign/padding boundaries. The pre-existing oracle for 400,000 generated coordinates plus directed cases remains. Native MOL library/integration validation: {validation['native']['passed']} passed, {validation['native']['failed']} failed, {validation['native']['ignored']} existing ignored across {validation['native']['suites']} suites. Release all-target Clippy, formatting, unchanged refactor ratchet and ten benchmark-runner tests pass. All 4,384 native source parse snapshots and rewritten records match byte for byte.

## Fixed protocol and provenance

Ten fixed profiles: two unfiltered 1,000-row diverse SDF slices; corresponding V3000 slices; 64 linear chains; continued versions; two synthetic 3D chain profiles and two fluorine-branched versions. Synthetic chains span 16..1024 carbon atoms. Their coordinates are arbitrary deterministic finite nonplanar coordinates, not physical conformers, embedding quality or representative holdout molecules. Both SDF write source CI lower bounds must exceed 1; every one of the 30 unrounded paired source median ratios must be >=0.99. Negative intervals, semantic rejections and floor failures remain visible. One valid epoch per profile, no adaptive timing reruns or pooling.

Twenty-one rotating three-arm fresh-process blocks per profile; three native batch repetitions, parse/prepared write/roundtrip. Input loading, process startup and plugin discovery are excluded. V3000 framing is outside the timer, SDF framing inside. Diagnostic Rust readers retain full batches. Prepared writing receives one untimed warm write; SDF reuses a Rust String while V3000 and Open Babel allocate. Whole-process high-water RSS includes preparation, retained inputs and warmup; it is not operation-only or streaming memory. External host load is uncontrolled and logged per block. Bootstrap intervals describe these paired measurements, not every host or molecule.

Separate source roots and initially empty targets, Rust opt3/full LTO/one codegen unit/one job/incremental=false/strip=none. Baseline files are verified against the named Git base. Both verbose build logs compile the production MOL library; actual `push_fixed4` machine-code bodies differ. Before any formal samples, a mandatory native test reproduced an existing main regression in both baseline and candidate: an explicit-H ring CH counted itself as its own partner. The separately committed #986 repair excludes self-partners and retains distinct partners. Baseline MOL2 tests and the complete candidate MOL suite pass after applying it identically to both arms. Original build/seal records and both failures remain in `unmeasured-main-failure/`. Inputs/hashes, primary modes, block counts and all adoption thresholds are unchanged. Two unused candidate files were removed from the initial baseline copy before source sealing; baseline production never referenced them, no source code changed and no formal samples existed. {json.loads((O/'environment.json').read_text())['openbabel'].strip()} native C++ probes use -O3; RDKit {json.loads((O/'environment.json').read_text())['rdkit']} is the independent output observer. Seals cover tracked sources, baseline files, exact corpora, executables and package-local Open Babel files before/after each profile. System and package-external runtime data are outside that seal.

Ratios are reference/candidate elapsed times, deterministic paired 10,000-resample bootstrap, seeds 0 (source) and 1 (Open Babel), reported to three decimals. Source retention is assessed using unrounded medians. Open Babel parse ranking requires equivalent accepted canonical graphs. Write/roundtrip ranking additionally requires declared atom/bond fields, coordinates to 1e-4, names/properties and explicit unknown stereo. Rejected or mismatched arms remain unranked; a strict external win requires a CI lower bound above 1 and 21/21 faster blocks. Earlier source optimization ratios are not multiplied into these measurements.

## Complete matrix

Milliseconds per full batch, native elapsed divided by three repetitions. Paired ratios need not equal ratios of the displayed marginal medians.

| Profile | Mode | Base ms | Candidate ms | Base/candidate [95% CI]; wins | Open Babel/candidate [95% CI]; wins |
|---|---|---:|---:|---|---|
'''+ '\n'.join(rows)+'\n\n## Decision and all negative results\n\n```json\n'+json.dumps(decision,indent=2)+'\n```\n\n## Evidence and reproduction\n\n[Complete frozen packet]('+NAME+'-artifacts/README.md), including all ten raw matrices, exact inputs, complete shared source outcome streams and retained Open Babel outcomes, source/binary/native seals, compiler logs, disassembly, test logs and independently recomputed statistics. Full local evidence: `'+str(O)+'`. Rebuild the named revisions in independent roots/targets with `build.py`; run `native.py`, `seal.py` and `formal.py` with retained original absolute paths adapted. `assemble.py` independently verifies source ratios, eligible competitor ratios and the adoption decision. A later merge revision requires correctness validation and is not the timed candidate.\n'
(R/'benchmarks'/(NAME+'.md')).write_text(report)
index=R/'benchmarks/README.md';s=index.read_text();marker='## What counts as a speed win';assert marker in s
entry='- [Decimal-pair coordinate formatting]('+NAME+'.md): ten fixed Open Babel/source profiles, full byte oracles and all source retention results.\n\n'
if entry not in s:s=s.replace(marker,entry+marker);index.write_text(s)
print('Report and',len(packet),'hashed files assembled',flush=True)
