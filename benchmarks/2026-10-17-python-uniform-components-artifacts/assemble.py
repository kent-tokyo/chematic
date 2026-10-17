from pathlib import Path
import json,gzip,hashlib,subprocess,statistics,shutil,sys,difflib
O=Path(__file__).resolve().parent;R=Path('/private/tmp/chematic-python-asa-preparation-20261017');D=R/'benchmarks/2026-10-17-python-uniform-components-artifacts';D.mkdir(exist_ok=True)
plan=json.loads((O/'measurement-plan.json').read_text());decision=json.loads((O/'decision.json').read_text());native=json.loads((O/'native-validation.json').read_text());parity=json.loads((O/'parity.json').read_text());proof=json.loads((O/'binary-source-proof.json').read_text())
sys.path.insert(0,str(R/'validation/cosmolkit_comparison'));from bench_corpus import paired_median_ratio_ci
rows=[];audits=[];comparators={}
def fmt(r): return f"{r['median']:.3f} [{r['ci95'][0]:.3f}, {r['ci95'][1]:.3f}]; {r['winning_blocks']}/20"
for offset in plan['offsets']:
 data=json.loads((O/f'frontier-{offset}.json').read_text());assert data['complete'] and len(data['blocks'])==20
 for op in plan['operations']:
  assert data['agreement'][op]['candidate_vs_baseline']['passed']
  for mode in plan['modes']:
   key=f'{offset}/{op}/{mode}';summary=data['summary'][op+'/'+mode]
   values={arm:[b[arm]['operations'][op][mode] for b in data['blocks']] for arm in ['baseline','candidate','rdkit','cosmolkit']}
   ratio=paired_median_ratio_ci(values['baseline'],values['candidate'],seed=0);assert ratio==decision['source_ratios'][key]==summary['baseline']['speedup']
   cells=[]
   for seed,arm in [(1,'rdkit'),(2,'cosmolkit')]:
    result=summary[arm];r=result['speedup']
    if r is not None: assert r==paired_median_ratio_ci(values[arm],values['candidate'],seed=seed)
    cells.append(fmt(r)+(' strict win' if result['demonstrated_win'] else '') if r is not None else 'unranked: output mismatch')
   copies=plan['fixed_batch_copies'][op]
   rows.append(f"| {offset} | {op}/{mode} | {statistics.median(values['baseline'])/(copies*1e6):.3f} | {statistics.median(values['candidate'])/(copies*1e6):.3f} | {fmt(ratio)} | {cells[0]} | {cells[1]} |")
   audits.append({'lane':key,'source_ratio_recomputed':True,'ranked_comparator_ratios_recomputed':True})
   if op=='labute_asa' and mode=='first_use':comparators[str(offset)]=summary
assert len(audits)==42
(O/'audit.json').write_text(json.dumps({'passed':True,'all_42_source_ratios_recomputed':True,'all_ranked_comparator_ratios_recomputed':True,'full_frozen_function_and_binary_source_proof':proof['complete_frozen_function_exact'] and proof['source_contrast_valid'],'comparisons':audits},indent=2)+'\n')
(O/'primary-comparators.json').write_text(json.dumps(comparators,indent=2)+'\n')
patch=subprocess.check_output(['git','diff',proof['base_commit'],proof['candidate_commit'],'--','crates'],cwd=R);(O/'candidate.patch').write_bytes(patch)
old=(R/'scripts/bench_python_speed_frontier.py').read_text().splitlines(keepends=True);new=(O/'bench-fixed-batches.py').read_text().splitlines(keepends=True)
(O/'common-instrumentation.patch').write_text(''.join(difflib.unified_diff(old,new,fromfile='scripts/bench_python_speed_frontier.py',tofile='bench-fixed-batches.py')))
json_files=['measurement-plan.json','decision.json','source-manifest.json','comparator-manifest.json','binary-source-proof.json','protocol-seal.json','environment.json','frozen-reference.json','native-validation.json','native-check-commands.json','parity.json','instrumentation-check.json','audit.json','primary-comparators.json','baseline-build-command.json','candidate-build-command.json','baseline-python-test-command.json','candidate-python-test-command.json','profile-5000-retry1-status.json']
for name in json_files:shutil.copy2(O/name,D/name)
for name in ['build.py','implement.py','prepare-protocol.py','check-instrumentation.py','bench-fixed-batches.py','native.py','verify.py','verify_maccs.py','seal.py','formal.py','assemble.py','profile.py','profile-worker.py','profile-5000-retry.py']:
 shutil.copy2(O/name,D/name)
logs=['baseline-build.log','candidate-build.log','source-commit.log','shortcut-tests-attempt1.log','shortcut-tests.log','native-tests.log','clippy.log','format.log','refactor.log','runner-tests.log','baseline-python-tests.log','candidate-python-tests.log','profile.log','profile-0.txt','profile-5000-retry1.txt','profile-worker-0.json','profile-worker-5000-retry1.json','profile-worker-0.stderr','profile-worker-5000.stderr','profile-5000-retry1.log','candidate.patch','common-instrumentation.patch','formal.log','verify.log']
for arm in ['baseline','candidate']:
 for k in ['descriptor','morgan','maccs']:logs.append(f'{arm}-{k}.log')
 for idx in [0,1]:logs.append(f'{arm}-shortcut-{idx}-disassembly.txt')
for name in logs:(D/(name+'.gz')).write_bytes(gzip.compress((O/name).read_bytes(),mtime=0))
for offset in plan['offsets']:
 for suffix in ['json','log']:
  name=f'frontier-{offset}.{suffix}';(D/(name+'.gz')).write_bytes(gzip.compress((O/name).read_bytes(),mtime=0))
# One retained copy per byte-identical two-arm complete output stream.
for kind in ['descriptor','morgan','maccs']:
 assert (O/f'baseline-{kind}.jsonl').read_bytes()==(O/f'candidate-{kind}.jsonl').read_bytes()
 (D/f'equal-source-{kind}.jsonl.gz').write_bytes(gzip.compress((O/f'baseline-{kind}.jsonl').read_bytes(),mtime=0))
(D/'README.md').write_text('One-host source experiment. Includes every raw block, plan, all public-output streams (one copy of each byte-identical two-arm stream), source/native/wheel seals, actual production disassembly, tests, original setup errors, common-instrumentation patch and independent ratio recomputation. Exact wheels, installed extensions and separate compiler targets remain locally at '+str(O)+'. Retained scripts use original absolute paths; adapt paths before reproduction. No package or universal speed claim.\n')
packet={str(p.relative_to(D)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(D.rglob('*')) if p.is_file() and p.name!='packet.sha256.json'}
(D/'packet.sha256.json').write_text(json.dumps(packet,indent=2)+'\n')
status='PASS' if decision['local_adoption_gate_passed'] else 'FAIL: keep an unpromoted draft'
compact={k:v for k,v in decision.items() if k!='source_ratios'}
report=f'''# Python uniform-component preparation — 2026-10-17

Decision: **{status}**. This source experiment measures first-use Labute ASA, the complete 42-lane retention frontier, and output-qualified RDKit/COSMolKit comparisons. It does not establish toolkit-wide or published-package superiority.

## Change

Measured base `{proof['base_commit']}`; measured candidate `{proof['candidate_commit']}`. Baseline own-CPU samples attribute significant initial ASA work through the actual Python getter to aromatic descriptor preparation. The profiling process also parses/frees batches; only getter ancestry supports this finding. The first attempt to attach the second diagnostic sampler failed during process startup; its log remains retained, and a separate retry after startup succeeded. No diagnostic worker supplies formal timings.

For a molecule already accepted by the original RDKit-model safety gate, the new structural predicate proves that each cyclic component is fully explicit aromatic or saturated without an incident multiple bond. Aromatic edges must be cyclic with flagged endpoints; non-aromatic cyclic edges must have unflagged endpoints; cyclic multiple bonds and exocyclic multiple bonds touching a non-aromatic ring atom keep the original path. Neutral hypercoordinate group-16 states and the existing macrocycle donor rule retain their checks. The shared cyclic flags/components are prepared as before, but temporary component fact tables can be skipped when the original complete-preperceived result is provably Identity. The full fallback is unchanged. No ring algorithm, public API, dependency, unsafe code, or eager descriptor parsing is added.

## Correctness and provenance

The complete frozen base shortcut is verified byte-exact after renaming its declaration. It matches cold/repeated decisions for all 14,072 pinned graphs and their attempted Kekulé forms, all 4,096 five-atom/flag/H combinations, and directed fused/mixed/exocyclic/macrocycle/unusual cases. The first new-test compilation used a nonexistent convenience Kekulé API; the retained failed log supplied no formal timing, and the test was corrected to the existing core API before committing/building the candidate. Production source did not change during that correction.

Independent wheel processes retain exact complete descriptor, Morgan and MACCS streams: {parity['counts']['descriptor_outcomes']:,} descriptor outcomes (float bits), {parity['counts']['morgan_outputs']:,} Morgan outputs and {parity['counts']['maccs_outcomes']:,} MACCS outcomes over 14,072 inputs, zero parse/evaluation errors. Six existing nonfinite descriptor outcomes remain represented by their exact bits; they are not counted as valid external chemistry agreement. The ordered four-property dump and the separate fresh-per-operation timing workers check different preparation contexts.

Rust release validation: {native['native']['passed']:,} passed, {native['native']['failed']} failed, {native['native']['ignored']} existing ignored across {native['native']['suites']} suites (perception, SMILES, chemistry, fingerprints and MOL). Release all-target Clippy, formatting, refactor ratchet and benchmark-runner tests pass. Full Python suites: `{parity['python_tests']['baseline']}` / `{parity['python_tests']['candidate']}`. Thirty skips reflect APFS EILSEQ for non-UTF-8 filenames; the other two reflect optional elixcee/OGB dependencies unavailable in this environment.

Both wheels use separate initially empty targets and separate source roots, Rust opt3/full LTO/one CGU/one job, incremental=false, strip=none, native-inchi. Every baseline source file matches its Git revision. Both verbose logs compile the production perception library, and the actual shortcut/callback instruction hashes differ. Wheel/extension and source hashes are sealed, as are RDKit, COSMolKit and NumPy package-local native files. System libraries and package-external runtime data are outside that seal. All builds, checks and diagnostics finish before timing; no old wheel or timing samples are reused.

## Fixed protocol

CPython 3.13.6, RDKit 2026.09.1, COSMolKit 0.5.0, NumPy 2.5.3 on one shared macOS arm64 host. Two unfiltered 1,000-record slices (offsets 0 and 5,000) of the pinned browser 10k corpus. Twenty rotating four-arm fresh-process blocks per slice, sequential, one complete source-valid epoch each. Primary: both first-use Labute ASA paired CI lower bounds >1. Adoption also requires exact source outcomes and all 42 source median ratios >=0.99. Every negative interval is disclosed; no adaptive timing reruns, sample pooling or substitution.

Scalar operations (TPSA, ASA, ring count) use **20 copies of each full slice, parsed into 20,000 separate objects**. Fingerprints use the original single 1,000-record batch. This fixed choice was made before formal timing to reduce noise in historically roughly 40-microsecond warm scalar batches. A state-recording factory spy independently confirms 61,001 object creations, 40,001 objects evaluated once and 20,000 evaluated twice for one scalar worker: the initial parse probe/disjoint warm object, fresh first-use objects, once-warmed repeated-use objects, and genuine parse-inclusive objects. All preflight elapsed fields are discarded. These are repeated input records, not additional independent chemical cases or statistical blocks; comparisons with historical 1,000-object absolute timings are not claimed.

First-use parsing occurs outside the operation timer and may compare different eager/lazy preparation. Repeated-use warms the operation once per object. Parse-inclusive performs actual parsing and evaluation inside the same timer. Native result objects are retained; normalization/output comparisons run outside timing. Whole-process peak RSS includes imports, preparation and all operations; scalar batch sizes are larger than historical matrices. No isolated-operation memory claim. External host load is uncontrolled and recorded.

Ratios are reference/candidate elapsed times. Intervals use the unchanged deterministic paired bootstrap (10,000 resamples, three-decimal reporting). A strict external win requires equivalent outputs, CI lower>1 and all 20 blocks faster; mismatches remain unranked. Scalar mismatch counts include the 20 input copies. Source equivalence is exact; external scalar tolerance is 1e-9, fingerprints/rings require exact values. No Linux, Windows, WASM, independent holdout or published-artifact inference.

## Complete matrix

Displayed milliseconds normalize marginal batch medians to 1,000 calls by dividing scalar timings by 20. They are averages within a 20,000-object batch, not separate 1,000-object experiments. Ratios/intervals use the original raw paired batches.

| Offset | Operation / mode | Base ms per 1k | Candidate ms per 1k | Base/candidate [95% CI]; wins | RDKit/candidate [95% CI]; wins | COSMolKit/candidate [95% CI]; wins |
|---|---|---:|---:|---|---|---|
'''+ '\n'.join(rows)+'\n\n## Decision and every negative interval/floor failure\n\n```json\n'+json.dumps(compact,indent=2)+'\n```\n\n## Evidence and reproduction\n\n[Complete frozen packet](2026-10-17-python-uniform-components-artifacts/README.md): plan, both raw matrices, every retained public output, original validation/profile errors, compiler logs, production disassembly, common instrumentation patch, manifests and independent recomputation. Full local evidence: `'+str(O)+'`. Build the named revisions with `build.py` in separate targets; run `native.py`, `verify.py`, `check-instrumentation.py`, `seal.py`, then `formal.py`. Adapt the retained original absolute paths. `assemble.py` independently recomputes every source and ranked external ratio before writing this report. A later integration head, if any, requires fresh correctness checks and is not the timed candidate.\n'
(R/'benchmarks/2026-10-17-python-uniform-components.md').write_text(report)
p=R/'benchmarks/README.md';s=p.read_text();marker='## What counts as a speed win';assert marker in s
link='- [Uniform cyclic-component descriptor preparation](2026-10-17-python-uniform-components.md): exact full-source outcomes and all 42 fixed RDKit/COSMolKit retention lanes.\n\n'
assert '2026-10-17-python-uniform-components.md' not in s;s=s.replace(marker,link+marker);p.write_text(s)
print('Report and',len(packet),'recursively hashed evidence files assembled',flush=True)
