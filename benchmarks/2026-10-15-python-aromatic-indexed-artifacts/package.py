from pathlib import Path
import json,hashlib,gzip,tarfile,shutil,collections,importlib.util,re,statistics
r=Path('/private/tmp/chematic-parity-lazy-20261015');s=Path('/private/tmp/chematic-parity-lazy-source-20261015');m=json.loads((r/'source-manifest.json').read_text());sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert all(sha(s/n)==h for n,h in m['candidate_sha256'].items())
assert 'ALL VALIDATION PASSED' in (r/'validation.log').read_text()
assert '1389 passed, 1 skipped' in (r/'python-tests.log').read_text()
decision=json.loads((r/'decision.json').read_text());adopt=decision['adopt_source']
d=r/'delivery';a=d/'benchmarks/2026-10-15-python-aromatic-indexed-artifacts';a.mkdir(parents=True,exist_ok=True)
changed=m['changed_documentation_files']+m['changed_production_files']+m['changed_test_files'];new=m['added_test_files']
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
oldnames={'CHANGELOG.md':'baseline-CHANGELOG.md','crates/chematic-perception/src/aromaticity/model.rs':'baseline-model.rs','crates/chematic-perception/src/aromaticity/normalization_storage_tests.rs':'baseline-normalization-storage-tests.rs'}
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

lines=['# Python indexed aromatic normalization — 2026-10-15','',
('The production change is adopted in the attached source candidate.' if adopt else 'The production experiment is rejected for adoption; only its report and evidence are delivered.')+' This is unmerged source-version 1.1.4 evidence, without current-main verification or hosted CI. Complete superiority over RDKit or COSMolKit is not established.','',
'Decision: '+decision['reason'],'',
'## Change and correctness','',
'The candidate stores temporary atom flags and bond orders in dense index-aligned arrays instead of a second atom hash set and bond hash map. Molecule iterators enumerate dense indices in order, including after removal. Model membership predicates, promote-only versus authoritative atom semantics, aromatic-bond endpoint fallback, traversal order, implicit-H preservation, direction stashing, metadata copying and the existing patched-rebuild behavior are unchanged. The model itself still uses its original sets; no ring or aromaticity algorithm is changed. No public API, dependency or unsafe code is added.','',
'The immediate baseline includes all fourteen previously accepted source candidates, including the inactive E/Z guard and inline small-ring queue. Rejected prototypes are excluded. Both current wheels were built freshly; the baseline was built before edits and all 10,061 baseline source/input hashes were checked. No previous wheel or timing was reused.','',
'The frozen-oracle proof confirms that the complete immediate-baseline builders, bond-order map and finish function differ only by rustfmt whitespace. This independent copy checks 14,072 pinned inputs in parsed and attempted Kekule forms (original clone on Kekule error), with both aromaticity algorithms and both promote-only/authoritative builders: 112,576 full graph/metadata comparisons. Additional empty, reordered-adjacency, stale-flag, arbitrary-model and metadata cases cover authoritative aromatic-bond endpoint fallback and ignored out-of-range model entries. The older normalization storage oracle remains intact; its call adapter uses the new production builder. Independent baseline/candidate processes additionally produce byte-identical full descriptor and Morgan dumps: 112,576 descriptor outcomes including exact floating bits and 56,288 Morgan outputs, with zero parse/evaluation errors. Existing licensed input corpora are reused; no new molecular records were added.','',
f"Scoped release Rust tests: {native['passed']:,} passed, {native['failed']} failed, {native['ignored']} existing ignored. Python: 1,389 passed, one skipped. Release all-target Clippy for chemistry/perception, formatting and the unchanged refactor budget (487 production files, 8,063 functions, 18 file allowances, 88 function allowances) passed. Full logs are retained. Builds and CPU sampling completed before formal timing; only owned completed build targets were removed.",'',
'Fresh baseline CPU samples are retained in sample-*.txt and getter-profile-review.json. At offset 5,000, the normalization builder appears beneath 468/2,084 TPSA getter samples and 462/2,194 Labute ASA getter samples. At offset 0 it has no separately named samples. Counts are filtered to the descriptor getter ancestry; inclusive function counts overlap and are not additive, and LTO inlining can hide names. Sampling identifies a plausible target, not a predicted speedup. Initial sandbox sampling, the corrected shared-import compilation failure and a separate SIGTERM-interrupted validation compilation are retained. The interruption had no completed test results; the same source/tests resumed from Cargo fingerprints before any formal timing.','',
'## Protocol','',
'Pinned environment: CPython 3.13.6, RDKit 2026.09.1, COSMolKit 0.5.0-rc.22, Rust 1.97.0, maturin 1.13.3, macOS arm64. Two fixed unfiltered 1,000-input slices, offsets 5,000 then 0, use 20 rotating fresh-process four-arm blocks, seven operations and all three preparation modes. This is one exposed development corpus on one host, without a sealed holdout or toolkit-wide claim.','',
'Ratios are reference elapsed / candidate elapsed; greater than 1 is faster. Brackets are paired median-ratio 95% bootstrap intervals. A strict win requires semantic agreement, CI lower > 1 and all 20 blocks faster. Output mismatches suppress ranking; all observations and failures are retained. Whole-process RSS includes imports and preparation.','',
'First-use prepares fresh parsed objects outside timing and can compare different eager/lazy preparation. Repeated-use warms the operation once on each object; parse-inclusive times fresh parsing and evaluation. Return objects are retained during timing; normalization and semantic checks occur outside timing. Every block validates output determinism and wheel/extension identity.','',
'An external checkpoint wrapper leaves frozen benchmark workers, timers, output checks and selection policy unchanged. It reserves 16 MiB of owned storage and on ENOSPC releases that reserve and retries writing the same computed checkpoint for at most 60 seconds. It never removes or remeasures a performance observation. Wrapper logs retain any such event. Matrix completion and any checkpoint events are recorded in the packet; all initial and selected adaptive matrices are retained.','',
'Initial comparator strict wins: '+', '.join(f"{mode}: {v['strict_wins']}/{v['eligible']}" for mode,v in counts.items())+'. COSMolKit has six unequal ring counts at offset 0; the corresponding modes remain unranked. Every remaining deficit is shown below.','']
for offset,j in raw.items():lines += [f'## Initial matrix: offset {offset}','',table(j),'']
lines += ['## Adaptive observations','',
'The plan was saved before timing: an initial source-relative CI upper < 1 selects its entire operation for one separate 20-block follow-up in all modes. Initial and follow-up observations remain separate, without pooling or recursive retries. All new follow-up negatives remain recorded.','', '```json',json.dumps(selection,indent=2),'```','']
for offset,j in focused.items():lines += [f'### Separate follow-up: offset {offset}','',table(j),'']
lines += ['## Reproduction and status','',
'[The evidence packet](2026-10-15-python-aromatic-indexed-artifacts/README.md) contains both measured wheels, every initial/adaptive matrix, compressed exact output dumps, build/test and profile logs, original implementation, exact patch, and all frozen source/input files. Archive readback checks all '+str(len(observed))+' candidate hashes; baseline overrides and removal of the added test reconstruct all '+str(len(reconstructed))+' baseline hashes. The patch independently reconstructs changed files in a tiny Git repository. Hashes demonstrate source reconstruction, not bit-reproducibility of rebuilt wheels.','',
'Last verified Git base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`. Shared Git metadata was not refreshed; this is not merged, published or hosted-CI-validated. Prior unresolved regressions against older source candidates, including V3000 reading and Morgan/MACCS observations, are not resolved by an immediate-baseline comparison. OpenBabel format lanes are outside this descriptor experiment. Separate-run win counts do not establish attributable algorithm gains.','']
report=d/'benchmarks/2026-10-15-python-aromatic-indexed.md';report.write_text('\n'.join(lines))
(a/'README.md').write_text("""# Indexed aromatic normalization evidence

Read decision.json, measurement-plan.json, source-manifest.json, both frontier matrices and adaptive-selection.json. Any focused matrices are separate follow-ups. All blocks, modes, gates and hashes are retained.

Extract candidate-build-source.tar.gz to obtain source/. Replace paths with baseline-overrides/ and remove the added test listed in the manifest to reconstruct the immediate baseline. Source-reconstruction.json checks every hash. Read the adoption decision before using an archived experiment.

The wheel directories contain exact measured binaries. Install separately with Python 3.13 and pip install --no-deps --target <directory> <wheel>. Comparator dependencies are pinned in comparator-environment.txt. Build with recorded release LTO, extension-module/native-inchi features, Rust and maturin versions.

Run the extracted scripts/bench_python_speed_frontier.py using --corpus source/validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi --offset 5000 --limit 1000 --blocks 20, and supply baseline/candidate wheel and installation paths, source manifest, base revision and output path. Repeat at offset 0. run-frontiers.py and run-followups.py record exact host commands and selection policy; paths need adaptation on another host. Gzipped exact dumps match manifest hashes. artifact-sha256.json indexes all other packet files.
""")
spec=importlib.util.spec_from_file_location('names',s/'scripts/check_product_names.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
assert not module.scan_paths(d,[report,a/'README.md']+([d/'CHANGELOG.md'] if adopt else []))
index={p.relative_to(a).as_posix():sha(p) for p in sorted(a.rglob('*')) if p.is_file() and p.name!='artifact-sha256.json'};(a/'artifact-sha256.json').write_text(json.dumps(index,indent=2)+'\n')
(r/'delivery-sha256.json').write_text(json.dumps({p.relative_to(d).as_posix():sha(p) for p in sorted(d.rglob('*')) if p.is_file()},indent=2)+'\n')
print('Packaged',len(index),'artifacts; adopted=',adopt,'validation=',native,'strict_wins=',counts)
