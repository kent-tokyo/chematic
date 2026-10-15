from pathlib import Path
import json,hashlib,gzip,tarfile,shutil,collections,importlib.util,re,statistics
r=Path('/private/tmp/chematic-descriptor-profile-20261015');s=Path('/private/tmp/chematic-descriptor-profile-source-20261015');m=json.loads((r/'source-manifest.json').read_text());sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert all(sha(s/n)==h for n,h in m['candidate_sha256'].items())
assert 'ALL VALIDATION PASSED' in (r/'validation.log').read_text()
assert '1389 passed, 1 skipped' in (r/'python-tests.log').read_text()
decision=json.loads((r/'decision.json').read_text());adopt=decision['adopt_source']
d=r/'delivery';a=d/'benchmarks/2026-10-15-python-normalization-storage-artifacts';a.mkdir(parents=True,exist_ok=True)
changed=m['changed_documentation_files']+m['changed_production_files'];new=m['added_test_files']
if adopt:
 for n in changed+new:
  p=d/n;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(s/n,p)
for n in ['source-manifest.json','baseline-source.json','normalization-storage.patch','adaptive-selection.json','measurement-plan.json','decision.json','comparator-environment.txt','optimize.py','validate.py','freeze.py','run-frontiers.py','run-followups.py','package.py','profile.py','profile-worker.py','baseline-cache-cleanup.json','native-cache-cleanup.json','completed-build-cache-cleanup.json','measure.py','deliver.py','prepare.py','prepare-package.py','resume-validation.py','profile-summary.json','after-validation.py']:
 shutil.copyfile(r/n,a/n)
for p in sorted(r.glob('*.log')):shutil.copyfile(p,a/p.name)
for glob in ['frontier-*.json','focused-*.json','sample-*.txt','sample-*-worker.jsonl']:
 for p in sorted(r.glob(glob)):shutil.copyfile(p,a/p.name)
for arm in ['baseline','candidate']:
 p=next((r/(arm+'-wheel')).glob('*.whl'));q=a/(arm+'-wheel')/p.name;q.parent.mkdir();shutil.copyfile(p,q)
oldnames={'CHANGELOG.md':'baseline-CHANGELOG.md','crates/chematic-perception/src/aromaticity/model.rs':'baseline-model.rs'}
for n,f in oldnames.items():
 p=a/'baseline-overrides'/n;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(r/f,p)
for arm in ['baseline','candidate']:
 for kind in ['descriptors','morgan']:
  src=r/(arm+'-'+kind+'.jsonl');dst=a/(src.name+'.gz')
  with dst.open('wb') as out,gzip.GzipFile(filename='',mode='wb',fileobj=out,mtime=0) as z:z.write(src.read_bytes())
  assert hashlib.sha256(gzip.decompress(dst.read_bytes())).hexdigest()==sha(src)
with tarfile.open(a/'candidate-build-source.tar.gz','w:gz') as tar:
 for n in m['candidate_sha256']:tar.add(s/n,arcname='source/'+n,recursive=False)
with tarfile.open(a/'candidate-build-source.tar.gz') as tar:
 observed={member.name.removeprefix('source/'):hashlib.sha256(tar.extractfile(member).read()).hexdigest() for member in tar.getmembers()}
assert observed==m['candidate_sha256']
reconstructed=observed.copy()
for n,f in oldnames.items():reconstructed[n]=sha(a/'baseline-overrides'/n)
for n in new:del reconstructed[n]
assert reconstructed==m['baseline_source']['sha256']
(a/'source-reconstruction.json').write_text(json.dumps({'patch_exact':m['patch_reconstruction_exact'],'archive_exact':True,'baseline_overrides_exact':True,'candidate_files':len(observed),'baseline_files':len(reconstructed)},indent=2)+'\n')
raw={offset:json.loads((r/f'frontier-{offset}.json').read_text()) for offset in [0,5000]};focused={offset:json.loads((r/f'focused-{offset}.json').read_text()) for offset in [0,5000] if (r/f'focused-{offset}.json').exists()};selection=json.loads((r/'adaptive-selection.json').read_text())
for j in list(raw.values())+list(focused.values()):
 assert j['complete'] and len(j['blocks'])==20
 assert all(v['candidate_vs_baseline']['passed'] for v in j['agreement'].values())
def fmt(v):
 x=v.get('speedup')
 if not x:return 'unranked (output mismatch)'
 return f"{x['median']:.3f} [{x['ci95'][0]:.3f}, {x['ci95'][1]:.3f}]; {x['winning_blocks']}/{x['blocks']}"+(' **strict win**' if v['demonstrated_win'] else '')
def table(j):
 lines=['| Operation / preparation | Immediate baseline | RDKit | COSMolKit |','|---|---:|---:|---:|']
 for k,v in j['summary'].items():lines.append('| '+k+' | '+' | '.join(fmt(v[name]) for name in ['baseline','rdkit','cosmolkit'])+' |')
 return '\n'.join(lines)
counts={}
for mode in ['first_use','repeated_use','parse_inclusive']:
 vals=[v[name] for j in raw.values() for k,v in j['summary'].items() if k.endswith('/'+mode) for name in ['rdkit','cosmolkit']]
 counts[mode]={'strict_wins':sum(v['demonstrated_win'] for v in vals),'eligible':sum(bool(v.get('speedup')) for v in vals),'suppressed':sum(not v.get('speedup') for v in vals)}
(a/'strict-win-counts.json').write_text(json.dumps(counts,indent=2)+'\n')
parity={}
for kind in ['descriptors','morgan']:
 counter=collections.Counter()
 for line in (r/('candidate-'+kind+'.jsonl')).open():
  row=json.loads(line)
  if 'index' not in row:continue
  counter['inputs']+=1;counter['parse_errors']+='parse_error' in row
  for group in ['first','repeated','values']:
   for out in row.get(group,{}).values():counter['outcomes']+=1;counter['evaluation_errors']+='error' in out
 parity[kind]=dict(counter)
assert parity['descriptors']=={'inputs':14072,'parse_errors':0,'outcomes':112576,'evaluation_errors':0},parity
assert parity['morgan']=={'inputs':14072,'parse_errors':0,'outcomes':56288,'evaluation_errors':0},parity
(a/'parity-counts.json').write_text(json.dumps(parity,indent=2)+'\n')
results=re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored;', (r/'native-tests.log').read_text());native={k:sum(int(x[i]) for x in results) for i,k in enumerate(['passed','failed','ignored'])};assert native['failed']==0
(a/'native-test-counts.json').write_text(json.dumps(native,indent=2)+'\n')
lines=['# Python aromatic normalization storage — 2026-10-15','',
('The production change is adopted in the attached source candidate.' if adopt else 'The production experiment is rejected for adoption; only its report and evidence packet are delivered.')+' It is unmerged and does not establish complete superiority over RDKit or COSMolKit, or verify current main.','',
'CPU sampling identified aromatic normalization as a material part of cold descriptors. The candidate reserves known atom/bond capacities in both normalization builders and reuses the original per-atom implicit-H vector as an indexed correction map. This removes a second correction-vector allocation and per-atom linear searches. The implicit-H calls, aromaticity model, bond-order decisions, atom/bond enumeration and metadata-copy sequence remain unchanged, including historical behavior in the patched rebuild. Two pre-existing formatting fragments in the same included source file are normalized. No public API, dependency or unsafe code is added.','',
'The immediate baseline includes all eight previously accepted source candidates: cycle-root pruning, V3000 writer buffers, shared Kekulé descriptor views, V3000 parser token reuse, ring-component counting, the no-double-bond E/Z skip, aromaticity component reuse and V3000 positive-zero parsing. It excludes all rejected experiments. Both Python wheels were built freshly before/after this single production change; this avoids attributing inherited wins to the new optimization.','',
'Decision: '+decision['reason'],'',
'All 14,072 fixed inputs have byte-identical descriptor and Morgan dumps: 112,576 descriptor outcomes (first/repeated, exact float bits), 56,288 Morgan outputs (R2/R3 with both chirality settings), zero parse/evaluation errors. The added differential test freezes the original normalization finish function and compares every one of the 10,000 existing development inputs under both Hückel and RDKit-like models, on parsed graphs and their successful Kekulé forms. It compares atoms, bonds, atom tags, R-group labels, adjacency order, implicit H, stereo groups with bond memberships, neighbor order, bond-direction markers/anchors and SMILES closure metadata. It asserts that hydrogen corrections are exercised. Existing native/Python suites separately cover tags and other metadata. No new chemical input records were added.','',
f"Scoped release Rust tests: {native['passed']:,} passed, {native['failed']} failed, {native['ignored']} existing ignored. Full Python suite: 1,389 passed, one skipped. Release all-target Clippy for perception and chemistry, formatting and the refactor budget passed. The complete logs are retained. The initial native-test build stopped on an incorrect include-module path; that failure log is retained with the corrected successful rerun. Only this run owned regenerable targets were removed between phases. Timing starts after all builds and validation finish.",'',
'## Diagnostic profile','', 'Four separate 8-second workers repeatedly parse all 1,000 rows and evaluate TPSA or Labute ASA, with fresh descriptor state. macOS sample records five seconds at a one-millisecond requested interval. These are diagnostic CPU samples, not benchmark timings or attribution of wall time. NCI Labute ASA has 2,063 of 3,601 main-thread samples under descriptor aromaticity and 585 under its normalized molecule build, summed across nonrecursive call sites; full trees and worker checksums are retained. The no-N/O input audit uses pinned RDKit over all unfiltered rows and is retained in profile-summary.json. No chemical rows were added or filtered.', '', '## Protocol','',
'Pinned measured environment: RDKit 2026.09.1, COSMolKit 0.5.0-rc.22, CPython 3.13.6, Rust 1.97.0, maturin 1.13.3, macOS arm64. The fixed, unfiltered slices contain 1,000 inputs each, offsets 0 and 5,000. Each uses 20 rotating fresh-process blocks across baseline, candidate and the two comparators. This is an exposed development corpus on one host, without a sealed holdout or toolkit-wide claim.','',
'Ratios are reference elapsed time / candidate elapsed time; greater than 1 means faster. Brackets are paired median-ratio 95% bootstrap intervals. Strict wins require semantic agreement, CI lower > 1 and all 20 blocks faster. Every operation, mode, block, error and mismatch is retained. Unchanged-operation movement is not evidence of a new algorithm.','',
'First-use parses fresh objects outside the timer and may compare different eager/lazy work; repeated-use warms the operation once per object; parse-inclusive times actual fresh parsing and evaluation. Native return objects remain retained during timing; normalization and semantic checks occur outside timing. Whole-process peak RSS includes imports and preparation and is not per-operation allocation. Source, corpus, wheels and installed extension hashes are verified.','',
'Initial comparator strict wins: '+', '.join(f"{mode}: {v['strict_wins']}/{v['eligible']}" for mode,v in counts.items())+'. COSMolKit produces six unequal ring counts at offset 0, so all three corresponding modes remain unranked. The full remaining descriptor deficits are retained below.','']
for offset,j in raw.items():lines += [f'## Initial matrix: offset {offset}','',table(j),'']
lines += ['## Adaptive negative observations','',
'The plan was saved before timing: any initial baseline-relative lane with CI upper < 1 selects its entire operation for one separate 20-block follow-up in all modes. Initial observations remain above; follow-ups are not pooled or substituted. New follow-up negatives remain recorded without recursive retries.','', '```json',json.dumps(selection,indent=2),'```','']
for offset,j in focused.items():lines += [f'### Separate follow-up: offset {offset}','',table(j),'']
lines += ['## Reproduction and status','',
'[The evidence packet](2026-10-15-python-normalization-storage-artifacts/README.md) contains both measured wheels, every raw matrix/follow-up, exact compressed output dumps, build/test logs, the original implementation, the patch and all frozen source/input files. The archive is read back and all '+str(len(observed))+' candidate file hashes verified; baseline overrides and removal of the added test reconstruct all '+str(len(reconstructed))+' baseline hashes. The patch also reconstructs the changed files in an independent tiny Git repository. Hashes verify source reconstruction, not bit-reproducibility of rebuilt wheels.','',
'Last verified Git base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`. Shared metadata was not refreshed; this is source-version 1.1.4 evidence with no hosted CI, merge or publication. Previous unresolved regressions against older candidates, including chiral Morgan parse-inclusive, ring/Morgan repeated-use and MACCS observations, are not resolved by an immediate-baseline comparison. Separate run win counts must not be treated as attributable code improvements.','']
report=d/'benchmarks/2026-10-15-python-normalization-storage.md';report.write_text('\n'.join(lines))
(a/'README.md').write_text('''# Aromatic normalization storage evidence

Read `decision.json`, `measurement-plan.json`, `source-manifest.json`, both `frontier-*.json` reports and `adaptive-selection.json`. All blocks, modes, output gates, load observations and installed extension hashes remain. Any `focused-*.json` reports are separate adaptive diagnostics.

Extract `candidate-build-source.tar.gz` to obtain `source/`. To reconstruct the immediate baseline, replace paths using `baseline-overrides/` and remove the added test listed in the source manifest. `source-reconstruction.json` checks every source/input hash. The patch reproduces the candidate changes. Source version is 1.1.4; this is not main or a release.

The two wheel directories contain the exact measured binaries. Install each separately using Python 3.13 and `pip install --no-deps --target <directory> <wheel>`. Pin comparator dependencies from `comparator-environment.txt`. Run the extracted `scripts/bench_python_speed_frontier.py` with `--corpus source/validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi --offset 5000 --limit 1000 --blocks 20` and supply `--baseline-wheel`, `--baseline-path`, `--candidate-wheel`, `--candidate-path`, `--source-manifest`, `--base-revision` and `--output`. Repeat at offset 0. Use `run-followups.py` as the selection-policy reference; its run paths are host-specific.

Rebuild with the recorded Rust/maturin toolchains, release LTO and `extension-module,native-inchi` features. Gzipped exact dumps decompress to hashes in the source manifest. `artifact-sha256.json` indexes all other packet files. Read the decision before adopting an archived experiment.
''')
spec=importlib.util.spec_from_file_location('names',s/'scripts/check_product_names.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
paths=[report,a/'README.md']+([d/'CHANGELOG.md'] if adopt else [])
assert not module.scan_paths(d,paths)
index={p.relative_to(a).as_posix():sha(p) for p in sorted(a.rglob('*')) if p.is_file()};(a/'artifact-sha256.json').write_text(json.dumps(index,indent=2)+'\n')
(r/'delivery-sha256.json').write_text(json.dumps({p.relative_to(d).as_posix():sha(p) for p in sorted(d.rglob('*')) if p.is_file()},indent=2)+'\n')
print('Packaged',len(index),'artifacts; adopted=',adopt,'validation=',native, 'strict_wins=',counts)
