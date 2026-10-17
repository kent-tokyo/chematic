from pathlib import Path
import gzip,hashlib,json,statistics,shutil,subprocess
O=Path(__file__).resolve().parent;R=Path('/private/tmp/chematic-stereo-shape-filter-20261017')
D=R/'benchmarks/2026-10-17-openbabel-stereo-shape-filter-artifacts';D.mkdir(exist_ok=True)
plan=json.loads((O/'measurement-plan.json').read_text()); decision=json.loads((O/'decision.json').read_text());validation=json.loads((O/'validation.json').read_text());manifest=json.loads((O/'source-manifest.json').read_text())
base=manifest['base_commit'];candidate=manifest['candidate_commit']
src=subprocess.check_output(['git','show',base+':crates/chematic-perception/src/stereo_validation.rs'],cwd=R,text=True)
old=src[src.index('pub fn stereo_centers('):src.index('\n/// Summarise',src.index('pub fn stereo_centers('))].strip()
test=(R/'crates/chematic-perception/src/stereo_shape_tests.rs').read_text();frozen=test[test.index('fn frozen_stereo_centers('):test.index('\nfn compare(')].replace('fn frozen_stereo_centers(','pub fn stereo_centers(',1).strip();assert old==frozen
(O/'frozen-function-proof.json').write_text(json.dumps({'base_commit':base,'complete_frozen_function_exact':True,'trimmed_function_sha256':hashlib.sha256(old.encode()).hexdigest(),'pinned_graphs':14072,'five_atom_graph_hydrogen_variants':6144},indent=2)+'\n')
patch=subprocess.check_output(['git','diff',base,candidate,'--','crates'],cwd=R);(O/'candidate.patch').write_bytes(patch)
rows=[]; audit=[]
import sys
sys.path.insert(0,str(R/'validation/cosmolkit_comparison'))
from bench_corpus import paired_median_ratio_ci
for p in plan['profiles']:
 key=p['kind']+'-'+p['profile'];d=json.loads((O/('paired-'+key+'.json')).read_text())
 assert d['complete'] and len(d['blocks'])==21 and d['source_parse_snapshots_equal'] and d['rewrite_sha256']['baseline']==d['rewrite_sha256']['candidate']
 for mode in ['parse','write','roundtrip']:
  values={arm:[x[arm][mode]['elapsed_ns'] for x in d['blocks']] for arm in ['baseline','candidate','openbabel']}
  source=paired_median_ratio_ci(values['baseline'],values['candidate'],seed=0)
  assert source==decision['source_ratios'][key+'/'+mode]
  if d['summary'][mode]['baseline']['speedup'] is not None:
   assert source==d['summary'][mode]['baseline']['speedup']
  ob=d['summary'][mode]['openbabel'];speed=ob['speedup']
  if speed is not None: assert speed==paired_median_ratio_ci(values['openbabel'],values['candidate'],seed=1)
  def fmt(x): return f"{x['median']:.3f} [{x['ci95'][0]:.3f}, {x['ci95'][1]:.3f}]; {x['winning_blocks']}/21"
  rows.append('| '+key+' | '+mode+' | '+f"{statistics.median(values['baseline'])/3e6:.3f} | {statistics.median(values['candidate'])/3e6:.3f} | "+fmt(source)+' | '+(fmt(speed)+(' strict win' if ob['strict_speed_win'] else '') if speed else 'unranked: semantic prerequisite fails')+' |')
  audit.append({'lane':key+'/'+mode,'source_ratio_recomputed':True,'openbabel_ratio_recomputed':speed is not None})
assert len(audit)==30
(O/'audit.json').write_text(json.dumps({'passed':True,'all_30_source_ratios_recomputed':True,'all_ranked_comparator_ratios_recomputed':True,'no_filtered_rows_or_timing_reruns':True,'complete_frozen_function_exact':True,'comparisons':audit},indent=2)+'\n')
for name in ['measurement-plan.json','decision.json','source-proof.json','frozen-function-proof.json','binary-source-proof.json','source-manifest.json','binary-seal.json','environment.json','validation.json','native-check-commands.json','preflight-source-equality.json','audit.json','diagnostic-profile-scope.json']:
 shutil.copy2(O/name,D/name)
for name in ['setup.py','prepare.py','build.py','seal.py','formal.py','native.py','profile.py','assemble.py']:
 shutil.copy2(O/name,D/name)
for name in ['baseline-build.log','candidate-build.log','source-commit.log','native-tests.log','clippy.log','format.log','refactor.log','runner-tests.log','diagnostic-profile.txt','diagnostic-profile-worker.json','baseline-stereo-owner-disassembly.txt','candidate-stereo-owner-disassembly.txt','candidate.patch','formal.log']:
 (D/(name+'.gz')).write_bytes(gzip.compress((O/name).read_bytes(),mtime=0))
for p in plan['profiles']:
 key=p['kind']+'-'+p['profile']
 for suffix in ['json','log']:
  name='paired-'+key+'.'+suffix;(D/(name+'.gz')).write_bytes(gzip.compress((O/name).read_bytes(),mtime=0))
 (D/(p['path']+'.gz')).write_bytes(gzip.compress((O/p['path']).read_bytes(),mtime=0))
(D/'README.md').write_text('Frozen one-host source experiment. Plan, every raw timing block, exact corpora, source/native/executable seals, actual production disassembly, verbose compiler logs, frozen-function proof and independent statistical audit are retained. Full executable files and preflight snapshots/rewrites remain at '+str(O)+'. Build the named revisions in isolated per-arm targets; scripts retain original local paths and require path adaptation. No package or universal performance claim.\n')
failure=D/'setup-failure';failure.mkdir(exist_ok=True)
for n in ['failure.json','formal.log','measurement-plan.json','binary-seal.json','prepare.py','v3000-continued.molblocks.gz']: shutil.copy2(O/'setup-failure'/n,failure/n)
packet={str(p.relative_to(D)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(D.rglob('*')) if p.is_file() and p.name!='packet.sha256.json'}
(D/'packet.sha256.json').write_text(json.dumps(packet,indent=2)+'\n')
status='PASS' if decision['local_adoption_gate_passed'] else 'FAIL: retain as an unpromoted draft'
report=f'''# Open Babel comparison: reject impossible local stereocenter shapes before branch traversal

Decision: **{status}**. Exact source outputs are retained; primary and retention results below determine adoption. This is a source experiment on one macOS host, not universal superiority or a published-package result.

## Change and correctness

Measured base `{base}`; measured candidate `{candidate}`. In `stereo_centers`, CH₂ atoms can satisfy `groups == 4` but fail the later local output-shape condition. The unchanged condition now rejects those atoms before `substituent_signatures_distinct` builds rooted branch signatures and clones visited HashSets. The branch algorithm and valid-candidate order remain unchanged. No public API/dependency/unsafe additions. A baseline own-CPU sample confirms branch traversal/allocation is a major load on the branched 3D primary; that diagnostic worker's timing is excluded from every formal sample.

The complete frozen classifier is byte-exact to the base after only renaming its declaration. All 14,072 pinned graphs, all 6,144 five-atom/hydrogen-override graph combinations, and directed ring-nitrogen/explicit-H/charged/isotope/annotated boundaries match. Native validation: {validation['native']['passed']} passed, {validation['native']['failed']} failed, {validation['native']['ignored']} ignored across {validation['native']['suites']} suites, covering perception, MOL, chemistry and fingerprints. Release all-target Clippy, formatting, refactor ratchet and ten benchmark-runner tests pass. All 4,384 profile source snapshots and rewrites are byte-exact. The frozen classifier supplies direct public-result retention; this experiment does not add a Python/RDKit/COS performance claim.

## Fixed protocol and build provenance

Ten profiles: two 1,000-row diverse SDF slices, two corresponding V3000 slices, 64 linear chains, continued versions, two synthetic nonplanar linear-chain profiles, and two branched versions with F attached to carbon 2. Each synthetic set spans 16..1024 carbons. Coordinates are deterministic finite arbitrary coordinates, not physical conformers, embedding quality, or a fresh representative holdout. Both branched 3D parse lanes are the predeclared primary. Adoption requires exact source results, both primary paired 95% CI lower bounds above 1, and all 30 source medians at least 0.99. Every negative interval is disclosed. No sample filtering, epoch pooling or adaptive reruns.

Twenty-one rotating three-arm fresh-process blocks per profile, three native batch repetitions, parse/prepared write/roundtrip. Startup, file loading and plugin discovery are excluded; V3000 delimiter splitting is outside the timer, SDF framing inside. Parse retains the full diagnostic batch. Prepared writing has one untimed warm write. Rust SDF writing reuses a String; V3000 and Open Babel allocate strings. RSS is whole-process high-water usage including setup, batch retention and warmup, not isolated-operation or streaming memory. External host load is uncontrolled and recorded per block.

Rust release opt3/full LTO/one codegen unit/one job, incremental=false, strip=none. Separate initially empty targets and separate source roots for both arms. Every baseline source file is checked against the base Git revision. Both verbose logs explicitly compile the production perception library. Full LTO inlines `stereo_centers` into `perceive_tetrahedral_from_3d_with_mode`; that actual production owner has different machine-code bodies in the two arms, independently of example executable hashes. The initial preflight rejected a generated continuation-line prefix with one space instead of two. That failed setup supplied zero formal timing samples; the defective corpus, plan, seal and log remain in `setup-failure/`. Only that framing prefix was corrected, with unchanged graphs, rows, source, hypothesis and gates. No reused timing samples. Open Babel 3.2.1 native C++ probes are built with -O3; RDKit {json.loads((O/'environment.json').read_text())['rdkit']} observes graph and rewrite contracts. Source, baseline, corpus, binary and package-local Open Babel files are sealed before/after profiles; system and package-external runtime data are outside that seal.

Ratios are reference/candidate times with the unchanged deterministic 10,000-resample paired bootstrap. Open Babel parse ranking requires equivalent canonical graphs; write/roundtrip ranking requires accepted equivalent rewrites with declared atom/bond fields, coordinates to 1e-4, names/properties and unknown stereo. Mismatches/rejections remain present and unranked. A strict comparator win also requires 21/21 faster blocks. The earlier PR #966 whole-molecule guard and PR #960 writer optimization are already in this base; these ratios must not be multiplied by their historical ratios.

## Complete matrix

Milliseconds per full batch (native elapsed divided by three repetitions); ratios are paired medians, not ratios of the displayed marginal medians.

| Profile | Mode | Base ms | Candidate ms | Base/candidate [95% CI]; wins | Open Babel/candidate [95% CI]; wins |
|---|---|---:|---:|---|---|
'''+ '\n'.join(rows)+ '\n\n## Decision, all negative intervals and floor failures\n\n```json\n'+json.dumps(decision,indent=2)+'\n```\n\n## Evidence and reproduction\n\n[Complete frozen packet](2026-10-17-openbabel-stereo-shape-filter-artifacts/README.md), including all ten raw timing matrices, inputs, source and executable seals, compiler logs, disassembly, tests and audit. Full local evidence: `'+str(O)+'`. Build the two named source revisions using `build.py` in independent roots/targets, run `native.py`, `seal.py`, then `formal.py`; adapt retained original absolute paths. `assemble.py` independently recomputes every source/ranked comparator ratio and records every lane before assembling this report. A later merge revision, if any, requires fresh correctness validation and is not the timed candidate.\n'
(R/'benchmarks/2026-10-17-openbabel-stereo-shape-filter.md').write_text(report)
index=R/'benchmarks/README.md';s=index.read_text();marker='## What counts as a speed win';assert marker in s
s=s.replace(marker,'- [Local stereocenter shape filter](2026-10-17-openbabel-stereo-shape-filter.md): ten fixed Open Babel/source profiles with exact outputs and all retention results.\n\n'+marker);index.write_text(s)
print('report and',len(packet),'hashed evidence files assembled',flush=True)
