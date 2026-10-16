from pathlib import Path
import gzip, hashlib, json, shutil, statistics
ROOT=Path('/private/tmp/chematic-sdf-buffer-scan-20261016')
OUT=Path(__file__).resolve().parent
STEM='2026-10-16-python-scalar-model-reuse'
PACKET=ROOT/'benchmarks'/(STEM+'-artifacts')
plan=json.loads((OUT/'measurement-plan.json').read_text())
parity=json.loads((OUT/'parity.json').read_text())
native=json.loads((OUT/'native-validation.json').read_text())
matrix={f'{lane}-{offset}':json.loads((OUT/f'{lane}-{offset}.json').read_text()) for lane in ['frontier','bundle'] for offset in plan['offsets']}
assert all(d['complete'] and len(d['blocks'])==plan['blocks'] for d in matrix.values())
assert parity['outputs']['baseline']==parity['outputs']['candidate']
primary={k:d['summary']['tpsa_then_asa/first_use']['baseline']['speedup'] for k,d in matrix.items() if k.startswith('bundle')}
guard={k:{op:v['baseline']['speedup'] for op,v in d['summary'].items() if v['baseline']['speedup'] is None or v['baseline']['speedup']['median']<0.99} for k,d in matrix.items()}
negative={k:{op:v['baseline']['speedup'] for op,v in d['summary'].items() if v['baseline']['speedup'] and v['baseline']['speedup']['ci95'][1]<1} for k,d in matrix.items()}
passed=all(v and v['ci95'][0]>1 for v in primary.values()) and not any(guard.values())
counts={k:{arm:{'strict_wins':sum(v[arm]['demonstrated_win'] for v in d['summary'].values()),'eligible':sum(v[arm]['speedup'] is not None for v in d['summary'].values())} for arm in ['rdkit','cosmolkit']} for k,d in matrix.items()}
decision={'adoption_passed':passed,'primary':primary,'guard_failures':guard,'all_negative_source_intervals':negative,'comparator_strict_wins':counts,'source_lanes':sum(len(d['summary']) for d in matrix.values()),'universal_superiority':False,'epochs':'Four complete single 20-block matrices; no timing reruns, pooling or substitution'}
assert decision['source_lanes']==48
(OUT/'decision.json').write_text(json.dumps(decision,indent=2)+'\n')
def fmt(v):
 if v is None:return 'unranked: unequal output'
 return f"{v['median']:.3f} [{v['ci95'][0]:.3f}, {v['ci95'][1]:.3f}]; {v['winning_blocks']}/{v['blocks']}"
text=['# Scalar descriptor model reuse — 2026-10-16\n',f'Local adoption gate: **{"PASS" if passed else "FAIL"}**. RDKit/COSMolKit-wide superiority remains unproved. All measured lanes are retained.\n']
if not passed:text+=['This candidate must remain unpromoted: the predeclared primary or retention gate did not pass. Correctness CI does not override this performance decision.\n']
text+=['## Change and correctness\n',f"Base `{plan['base']}`, source version 1.1.5. Only scalar RDKit-model callbacks memoize `Option<Molecule>` after sanitization. Different scalar descriptors on the same original immutable graph can reuse that model. Generic callbacks keep their prior uncached model and aromatic-view behavior. The conversion algorithm, chemical gate, failure fallback and nested TLS bypass are unchanged. Clone and graph edits clear the derived model slot. No unsafe code, new dependency or altered descriptor formula.\n",
'This adds one `OnceLock<Arc<…>>` slot to every molecule and retains the corrected model (including any derived views) until that original molecule is dropped or edited. It trades memory for potential reuse. No universal or per-operation memory reduction is claimed.\n',
'The regression oracle freezes the complete original model callback function; `frozen-reference.json` verifies tokens after only renaming and removing whitespace. All 14,072 pinned input graphs preserve atom/bond data, adjacency, stereo neighbor order, bond directions, stereo groups/group bonds, tags and R-group labels, with both callback modes and repeated reads. Directed tests cover hydrogen/isotope handling, specified stereo, radicals/charges, dative cleanup, failed sanitization fallback, clone/edit invalidation, shared scalar model identity, ordinary graphs that need no model, generic callback isolation, nested calls, panic unwinding and concurrent cold readers.\n',
'The two source wheels also preserve all public descriptor float bits/integers, Morgan bytes and MACCS bytes across those 14,072 inputs: 112,576 descriptor outcomes, 56,288 Morgan outputs and 28,144 MACCS outcomes. Six pre-existing nonfinite descriptor outcomes are retained and excluded from accuracy claims. Development inputs are exposed, not a sealed holdout.\n',
'## Frozen protocol\n',
'Two unfiltered 1,000-row slices at offsets 0 and 5,000 of the pinned 10,000-input corpus. Four-arm rotating fresh processes: baseline, candidate, RDKit and COSMolKit. One 20-block standard matrix and one 20-block TPSA-then-ASA bundle matrix per slice. Standard operations: TPSA, Labute ASA, ring count, Morgan2 2048, chiral Morgan2, Morgan3 2048 and MACCS. Each has first-use, repeated-use and actual parse-inclusive modes: 42 standard plus six bundle source lanes. No adaptive timing rerun, sample pooling or substitution.\n',
'The bundle calls RDKit-default N/O-only TPSA followed by Labute ASA on the same molecule and retains both native return values inside timing. First use parses fresh objects outside timing; repeated use warms the entire operation/bundle once per object; parse inclusive parses and computes inside timing. All normalization and validation happen afterwards. The ordinary seven-operation runner is unchanged. Eager perception during comparator parsing differs from lazy first-use work. Imports/startup are excluded.\n',
'Rust 1.97.0; fresh CPython 3.13.6 wheels; opt3/full LTO/one codegen unit/one job, native-inchi, strip=none. RDKit 2026.09.1 and COSMolKit 0.5.0 native package files are rehashed before/after profiles. System/transitive libraries and package-external data are excluded, so this is not a full runtime closure. One shared macOS arm64 host, uncontrolled external activity; per-block load averages and whole-process RSS are retained. These are source results, not released-package or whole-toolkit wins.\n',
'Ratios are reference time / candidate time; brackets are deterministic paired-bootstrap 95% intervals (10,000 resamples, three-decimal rounding). Adoption requires exact public source outputs, CI lower >1 for both bundle first-use lanes, and all 48 source medians >=0.99. Every negative source interval is disclosed. A strict competitor win also requires matching outputs and 20/20 faster blocks. Different outputs remain unranked.\n',
'## Primary result\n','| Profile | Source ratio |\n|---|---:|']
for k,v in primary.items():text.append(f'| {k} / tpsa_then_asa/first_use | {fmt(v)} |')
for k,d in matrix.items():
 text += [f'\n## Complete matrix: {k}\n','| Operation / mode | Baseline | RDKit | COSMolKit |\n|---|---:|---:|---:|']
 for op,values in d['summary'].items():
  cells=[fmt(values[a]['speedup'])+(' **strict win**' if a!='baseline' and values[a]['demonstrated_win'] else '') for a in ['baseline','rdkit','cosmolkit']]
  text.append('| '+op+' | '+' | '.join(cells)+' |')
text+=['\n## Every negative source interval\n','```json\n'+json.dumps(negative,indent=2)+'\n```','\n## Every source median below 0.99\n','```json\n'+json.dumps(guard,indent=2)+'\n```','\n## Absolute times\n','Medians in milliseconds per complete 1,000-row batch.\n','| Profile / operation / mode | Baseline | Candidate | RDKit | COSMolKit |\n|---|---:|---:|---:|---:|']
for k,d in matrix.items():
 for key in d['summary']:
  op,mode=key.split('/')
  values=[statistics.median(b[a]['operations'][op][mode] for b in d['blocks'])/1e6 for a in ['baseline','candidate','rdkit','cosmolkit']]
  text.append('| '+k+' / '+key+' | '+' | '.join(f'{v:.3f}' for v in values)+' |')
text+=['\n## Whole-process memory\n','RSS includes imports, preparation, warmup and all worker operations; it does not isolate model retention or operation memory.\n','| Profile | Baseline MiB | Candidate MiB | RDKit MiB | COSMolKit MiB |\n|---|---:|---:|---:|---:|']
for k,d in matrix.items():
 values=[statistics.median(b[a]['peak_process_rss_bytes'] for b in d['blocks'])/2**20 for a in ['baseline','candidate','rdkit','cosmolkit']]
 text.append('| '+k+' | '+' | '.join(f'{v:.1f}' for v in values)+' |')
text+=['\n## Verification and retained evidence\n',f"Native core/smiles/chemistry tests: {native['passed']:,} passed, {native['failed']} failed, {native['ignored']} ignored across {native['suites']} suites (`--lib --tests`, debug profile with debug info disabled). Core/smiles all-target Clippy, format and refactor-budget checks passed. Python: 2,302 passed, two skipped. Three bundle semantic claim-gate tests passed. These checks are distinct from the optimized timed wheels.\n",
'Pre-timing failures were corrected before any formal worker: a new test initially referenced the wrong path for DerivedScalar; disk exhaustion interrupted candidate compilation and native test linking; two incomplete, nonexecutable linker outputs were removed and relinked. The stopped build/test logs and targeted cleanup records are retained. No failed validation supplied timing observations. The baseline snapshot includes only named Git-tracked build inputs, excluding ignored Python caches.\n',
f'[{STEM}-artifacts]({STEM}-artifacts/README.md) contains all raw matrices, fresh post-timing observer outputs and all COSMolKit mismatch rows, decisions, frozen plan, source/native/wheel manifests, candidate patch, drivers and logs. Full wheels, installed extensions, source snapshot and 14,072-input compressed output streams remain at `{OUT}` with hashes. Source reconstruction is verified; bit-reproducible rebuilds are not claimed.\n',
'More hosts, sealed holdout data, other operations, 3D/reactions, CLI startup and published Python/WASM artifacts remain open. A negative interval on an untouched API does not alone establish a causal code regression, but still counts under the fixed retention rule.\n']
(ROOT/'benchmarks'/(STEM+'.md')).write_text('\n'.join(text)+'\n')
PACKET.mkdir(exist_ok=True)
for p in OUT.iterdir():
 if p.is_file() and (p.suffix in ('.json','.py','.log','.patch') or p.name.startswith('observer-')):
  if p.name.startswith(('frontier-','bundle-')) and p.suffix=='.json':
   (PACKET/(p.name+'.gz')).write_bytes(gzip.compress(p.read_bytes(),mtime=0))
  else:shutil.copy2(p,PACKET/p.name)
(PACKET/'README.md').write_text(f'# Scalar model reuse evidence\n\nRead [the report](../{STEM}.md), frozen `measurement-plan.json`, and `decision.json`. Every matrix retains all 20 raw blocks and output hashes. Fresh post-timing observer values revalidate the recorded outputs; observer elapsed fields were discarded. `parity.json` hashes the complete 14,072-input source dumps, retained locally with the wheels at `{OUT}`.\n\n`candidate.patch` reconstructs all five changed/added candidate files from base `{plan["base"]}`. `packet-audit.json` verifies patch reconstruction, all source files, every paired statistic, output hashes and adoption. Compiler/build/test logs and provenance manifests are retained. This packet does not establish published-package or universal superiority.\n')
print(json.dumps(decision,indent=2))
