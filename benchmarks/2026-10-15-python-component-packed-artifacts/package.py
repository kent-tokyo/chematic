from pathlib import Path
import json,hashlib,gzip,tarfile,shutil,collections,importlib.util,re,statistics
r=Path('/private/tmp/chematic-component-packed-20261015');s=Path('/private/tmp/chematic-component-packed-source-20261015');m=json.loads((r/'source-manifest.json').read_text());sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert all(sha(s/n)==h for n,h in m['candidate_sha256'].items())
assert 'ALL VALIDATION PASSED' in (r/'validation.log').read_text()
assert '1389 passed, 1 skipped' in (r/'python-tests.log').read_text()
decision=json.loads((r/'decision.json').read_text());adopt=decision['adopt_source']
d=r/'delivery';a=d/'benchmarks/2026-10-15-python-component-packed-artifacts';a.mkdir(parents=True,exist_ok=True)
changed=m['changed_documentation_files']+m['changed_production_files'];new=m['added_test_files']
if adopt:
 for n in changed+new:
  p=d/n;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(s/n,p)
for p in sorted(r.iterdir()):
 if p.is_file() and p.suffix in ['.py','.json','.patch','.log','.txt'] and p.name not in ['delivery-sha256.json','delivery-verification.json']:shutil.copyfile(p,a/p.name)
for p in sorted(r.glob('*.log')):shutil.copyfile(p,a/p.name)
for glob in ['frontier-*.json','focused-*.json','sample-*.txt','sample-*-worker.jsonl']:
 for p in sorted(r.glob(glob)):shutil.copyfile(p,a/p.name)
for arm in ['baseline','candidate']:
 p=next((r/(arm+'-wheel')).glob('*.whl'));q=a/(arm+'-wheel')/p.name;q.parent.mkdir(exist_ok=True);shutil.copyfile(p,q)
oldnames={'CHANGELOG.md':'baseline-CHANGELOG.md','crates/chematic-perception/src/rdkit_parity/view.rs':'baseline-view.rs'}
for n,f in oldnames.items():
 p=a/'baseline-overrides'/n;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(r/f,p)
for arm in ['baseline','candidate']:
 for kind in ['descriptors','morgan']:
  src=r/(arm+'-'+kind+'.jsonl');dst=a/(src.name+'.gz')
  with dst.open('wb') as out,gzip.GzipFile(filename='',mode='wb',fileobj=out,mtime=0) as z:z.write(src.read_bytes())
  assert hashlib.sha256(gzip.decompress(dst.read_bytes())).hexdigest()==sha(src)
for src in r.glob('*-partial.jsonl'):
 dst=a/(src.name+'.gz')
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
lines=['# Python contiguous aromatic-component scratch — 2026-10-15','',
('The production change is adopted in the attached source candidate.' if adopt else 'The production experiment is rejected for adoption; only its report and evidence packet are delivered.')+' It is unmerged and does not establish complete superiority over RDKit or COSMolKit, or verify current main.','',
'The candidate combines three private u32 component-fact vectors into one allocation with adjacent root, cyclic-edge-count and atom-count slices. Their total element count stays n + 2c; the kind vector and public cached flags/components are unchanged. Traversal, component numbering, donor/candidate/Hückel predicates and arithmetic are unchanged. No public API, dependency or unsafe code is added. Both preceding FIFO and inline SmallVec trials are rejected, archived separately and excluded from this candidate.','',
'The immediate baseline includes all ten previously accepted source candidates, including aromatic normalization storage and the V3000 bounded-capacity/zero-map writer changes. All rejected experiments remain excluded. The exact baseline Python wheel is reused after matching all 10,057 baseline source/input hashes and its wheel/installed-extension hashes. It was built before either trial; the candidate is built freshly. New timing observations are separate and never pooled. The origin is recorded in baseline-artifact-origin.json.','',
'Decision: '+decision['reason'],'',
'All 14,072 fixed inputs have byte-identical descriptor and Morgan dumps: 112,576 descriptor outcomes (first/repeated, exact float bits), 56,288 Morgan outputs (R2/R3 with both chirality settings), zero parse/evaluation errors. A differential test freezes the original ComponentInfo constructor. All 14,072 inputs compare exact root labels, component kinds, edge/atom counts, cyclic flags, booleans and root-iteration order on repeated constructor calls. The frozen reference constructor populates shared ring caches before the candidate constructor; fresh/cold versus repeated descriptor output parity is separately checked by the two independent Python dump processes. No new chemical input records were added.','',
f"Scoped release Rust tests: {native['passed']:,} passed, {native['failed']} failed, {native['ignored']} existing ignored. Full Python suite: 1,389 passed, one skipped. Release all-target Clippy for perception and chemistry, formatting and the refactor budget passed. The complete logs are retained. The reused baseline build originally hit ENOSPC; its failure and successful retry logs are retained with the baseline origin record. Successful current candidate validation logs are retained. Only this run owned regenerable targets were removed between phases. Timing starts after all builds and validation finish.",'',
'## Measurement interruption','',
'ChEMBL initial collection stopped after 16 saved blocks when opening the block-17 checkpoint failed with transient ENOSPC. The just-computed four-arm block 17 was not persisted and its timings could not be recovered. Both the 16-block incomplete matrix and error logs are retained. The existing runner --resume revalidated the exact source, wheels, corpus, protocol and saved output hashes in untimed fresh workers, then recollected block 17 and continued through block 20. The original 16 saved blocks are unchanged. NCI completed without interruption and was not repeated. This collection boundary limits interpretation of small differences; no performance-based sample was removed and no measurements from another candidate were pooled. See measurement-recovery.json and frontier-0-resume.log.','',
'## Protocol','',
'Pinned measured environment: RDKit 2026.09.1, COSMolKit 0.5.0-rc.22, CPython 3.13.6, Rust 1.97.0, maturin 1.13.3, macOS arm64. The fixed, unfiltered slices contain 1,000 inputs each, offsets 0 and 5,000. Each uses 20 rotating fresh-process blocks across baseline, candidate and the two comparators. This is an exposed development corpus on one host, without a sealed holdout or toolkit-wide claim.','',
'Ratios are reference elapsed time / candidate elapsed time; greater than 1 means faster. Brackets are paired median-ratio 95% bootstrap intervals. Strict wins require semantic agreement, CI lower > 1 and all 20 blocks faster. Every successfully checkpointed operation, mode, block, error and mismatch is retained; the unavailable observation is disclosed below. Unchanged-operation movement is not evidence of a new algorithm.','',
'First-use parses fresh objects outside the timer and may compare different eager/lazy work; repeated-use warms the operation once per object; parse-inclusive times actual fresh parsing and evaluation. Native return objects remain retained during timing; normalization and semantic checks occur outside timing. Whole-process peak RSS includes imports and preparation and is not per-operation allocation. Source, corpus, wheels and installed extension hashes are verified.','',
'Initial comparator strict wins: '+', '.join(f"{mode}: {v['strict_wins']}/{v['eligible']}" for mode,v in counts.items())+'. COSMolKit produces six unequal ring counts at offset 0, so all three corresponding modes remain unranked. The full remaining descriptor deficits are retained below.','']
for offset,j in raw.items():lines += [f'## Initial matrix: offset {offset}','',table(j),'']
lines += ['## Adaptive negative observations','',
'The plan was saved before timing: any initial baseline-relative lane with CI upper < 1 selects its entire operation for one separate 20-block follow-up in all modes. Initial observations remain above; follow-ups are not pooled or substituted. New follow-up negatives remain recorded without recursive retries.','', '```json',json.dumps(selection,indent=2),'```','']
for offset,j in focused.items():lines += [f'### Separate follow-up: offset {offset}','',table(j),'']
lines += ['## Reproduction and status','',
'[The evidence packet](2026-10-15-python-component-packed-artifacts/README.md) contains both measured wheels, every raw matrix/follow-up, exact compressed output dumps, build/test logs, the original implementation, the patch and all frozen source/input files. The archive is read back and all '+str(len(observed))+' candidate file hashes verified; baseline overrides and removal of the added test reconstruct all '+str(len(reconstructed))+' baseline hashes. The patch also reconstructs the changed files in an independent tiny Git repository. Hashes verify source reconstruction, not bit-reproducibility of rebuilt wheels.','',
'Last verified Git base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`. Shared metadata was not refreshed; this is source-version 1.1.4 evidence with no hosted CI, merge or publication. Previous unresolved regressions against older candidates, including chiral Morgan parse-inclusive, ring/Morgan repeated-use and MACCS observations, are not resolved by an immediate-baseline comparison. Separate run win counts must not be treated as attributable code improvements.','']
report=d/'benchmarks/2026-10-15-python-component-packed.md';report.write_text('\n'.join(lines))
(a/'README.md').write_text('''# Contiguous component scratch evidence

Read `decision.json`, `measurement-plan.json`, `source-manifest.json`, both `frontier-*.json` reports and `adaptive-selection.json`. All blocks, modes, output gates, load observations and installed extension hashes remain. Any `focused-*.json` reports are separate adaptive diagnostics.

Extract `candidate-build-source.tar.gz` to obtain `source/`. To reconstruct the immediate baseline, replace paths using `baseline-overrides/` and remove the added test listed in the source manifest. `source-reconstruction.json` checks every source/input hash. The patch reproduces the candidate changes. Source version is 1.1.4; this is not main or a release.

The two wheel directories contain the exact measured binaries. Install each separately using Python 3.13 and `pip install --no-deps --target <directory> <wheel>`. Pin comparator dependencies from `comparator-environment.txt`. Run the extracted `scripts/bench_python_speed_frontier.py` with `--corpus source/validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi --offset 5000 --limit 1000 --blocks 20` and supply `--baseline-wheel`, `--baseline-path`, `--candidate-wheel`, `--candidate-path`, `--source-manifest`, `--base-revision` and `--output`. Repeat at offset 0. Use `run-followups.py` as the selection-policy reference; its run paths are host-specific.

Rebuild with the recorded Rust/maturin toolchains, release LTO and `extension-module,native-inchi` features. Gzipped exact dumps decompress to hashes in the source manifest. `artifact-sha256.json` indexes all other packet files. Read the decision before adopting an archived experiment.
''')
spec=importlib.util.spec_from_file_location('names',s/'scripts/check_product_names.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
paths=[report,a/'README.md']+([d/'CHANGELOG.md'] if adopt else [])
assert not module.scan_paths(d,paths)
index={p.relative_to(a).as_posix():sha(p) for p in sorted(a.rglob('*')) if p.is_file() and p.name!='artifact-sha256.json'};(a/'artifact-sha256.json').write_text(json.dumps(index,indent=2)+'\n')
(r/'delivery-sha256.json').write_text(json.dumps({p.relative_to(d).as_posix():sha(p) for p in sorted(d.rglob('*')) if p.is_file()},indent=2)+'\n')
print('Packaged',len(index),'artifacts; adopted=',adopt,'validation=',native, 'strict_wins=',counts)
