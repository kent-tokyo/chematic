from pathlib import Path
import json,hashlib,gzip,tarfile,shutil,statistics,re,importlib.util
r=Path('/private/tmp/chematic-openbabel-ring-queue-20261015');s=Path('/private/tmp/chematic-openbabel-ring-queue-source-20261015');d=r/'delivery';pre='2026-10-15-openbabel-ring-queue';a=d/'benchmarks'/(pre+'-artifacts');a.mkdir(parents=True,exist_ok=True);sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();m=json.loads((r/'source-manifest.json').read_text());assert all(sha(s/n)==h for n,h in m['frozen_build_source_sha256'].items());assert 'ALL VALIDATION PASSED' in (r/'validation-pipeline.log').read_text()
decision=json.loads((r/'decision.json').read_text());adopt=decision['adopt_source'];changed=m['only_changed_existing_source'];added=m['added_source_files']+m['added_test_files']
if adopt:
 for n in changed+added:
  p=d/n;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(s/n,p)
reports={profile:json.loads((r/f'paired21-{profile}.json').read_text()) for profile in ['0','5000','chains']};diagnostic=json.loads((r/'continued21.json').read_text());selection=json.loads((r/'adaptive-selection.json').read_text());focused={p.stem:json.loads(p.read_text()) for p in r.glob('focused21-*.json')}
for profile,j in reports.items():
 assert j['complete'] and len(j['blocks'])==21 and j['source_parse_snapshots_equal'] and j['source_manifest']==m
 assert (r/f'paired21-{profile}.artifacts/baseline.molblocks').read_bytes()==(r/f'paired21-{profile}.artifacts/candidate.molblocks').read_bytes()
for j in [diagnostic,*focused.values()]:assert j['complete'] and len(j['blocks'])==21
for p in sorted(r.iterdir()):
 if p.is_file() and p.suffix in ['.json','.patch','.log','.py','.rs','.partial'] and p.name not in ['delivery-sha256.json','delivery-verification.json']:shutil.copyfile(p,a/p.name)
for arm in ['baseline','candidate']:
 shutil.copyfile(r/(arm+'-v3000'),a/(arm+'-v3000'));(a/(arm+'-v3000')).chmod(0o755)
for folder in list(r.glob('paired21-*.artifacts'))+[r/'continuation.artifacts',r/'large-neighborhood.artifacts']:
 for p in sorted(folder.iterdir()):
  if p.is_file():
   q=a/folder.name/(p.name+'.gz');q.parent.mkdir(exist_ok=True)
   with q.open('wb') as out,gzip.GzipFile(filename='',mode='wb',fileobj=out,mtime=0) as z:z.write(p.read_bytes())
   assert hashlib.sha256(gzip.decompress(q.read_bytes())).hexdigest()==sha(p)
old=Path('/private/tmp/chematic-openbabel-v3000-20261014');oldfiles={'CHANGELOG.md':'baseline-CHANGELOG.md','crates/chematic-perception/src/stereo2d_ez_direction.rs':'baseline-ez.rs'}
with tarfile.open(a/'candidate-build-source.tar.gz','w:gz') as tar:
 for n in m['frozen_build_source_sha256']:tar.add(s/n,arcname='candidate/'+n,recursive=False)
 for n,f in oldfiles.items():tar.add(r/f,arcname='baseline-overrides/'+n,recursive=False)
 for profile in reports:tar.add(old/f'corpus-{profile}.molblocks',arcname='corpora/'+profile+'.molblocks',recursive=False)
 tar.add(r/'corpus-continued-chains.molblocks',arcname='corpora/continued-chains.molblocks',recursive=False)
 tar.add(r/'corpus-large-neighborhood.molblocks',arcname='corpora/large-neighborhood.molblocks',recursive=False)
with tarfile.open(a/'candidate-build-source.tar.gz') as tar:
 hashes={member.name.removeprefix('candidate/'):hashlib.sha256(tar.extractfile(member).read()).hexdigest() for member in tar.getmembers() if member.name.startswith('candidate/')};assert hashes==m['frozen_build_source_sha256'];baseline=hashes.copy()
 for n in changed:baseline[n]=hashlib.sha256(tar.extractfile('baseline-overrides/'+n).read()).hexdigest()
 for n in added:del baseline[n]
 assert baseline==m['baseline_source_files_sha256']
 assert hashlib.sha256(tar.extractfile('corpora/large-neighborhood.molblocks').read()).hexdigest()==sha(r/'corpus-large-neighborhood.molblocks')
counts={'passed':0,'failed':0,'ignored':0}
for name in ['native-library-tests.log','native-integration-tests.log']:
 for values in re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored;', (r/name).read_text()):
  for k,v in zip(counts,values):counts[k]+=int(v)
assert counts['passed']>=1138 and counts['failed']==0,counts
(a/'verification.json').write_text(json.dumps({'source_archive_exact':True,'baseline_reconstruction_exact':True,'candidate_file_count':len(hashes),'baseline_file_count':len(baseline),'native_tests':counts,'source_parse_and_rewrite_equality':{profile:{'source_parse_snapshots_equal':j['source_parse_snapshots_equal'],'snapshot_sha256':j['snapshot_sha256'],'rewrite_sha256':j['rewrite_sha256']} for profile,j in reports.items()},'continued_input_equivalence':json.loads((r/'continuation-equivalence.json').read_text()),'adaptive_selection':selection},indent=2)+'\n')
def fmt(v):
 z=v['speedup'];return 'unranked (semantic mismatch)' if z is None else f"{z['median']:.3f} [{z['ci95'][0]:.3f}, {z['ci95'][1]:.3f}]; {z['winning_blocks']}/21"+(' **strict win**' if v.get('strict_speed_win') else '')
def ms(j,arm,mode):return statistics.median(b[arm][mode]['elapsed_ns'] for b in j['blocks'])/1e6

oracle_match=re.search(r'ring_queue_oracle (\{[^\n]+\})',(r/'native-library-tests.log').read_text());assert oracle_match
oracle=json.loads(oracle_match[1]);assert oracle['rows']==14072 and oracle['queries']>1_000_000 and oracle['hashed_cases']>0 and oracle['inline_cases']>0
(a/'ring-search-oracle-counts.json').write_text(json.dumps(oracle,indent=2)+'\n')
lines=['# Inline E/Z ring-search queue — 2026-10-15','',
('The production optimization is adopted in the attached source candidate.' if adopt else 'The production experiment is rejected for adoption; only its report and evidence are delivered.')+' This does not establish toolkit-wide superiority, current main or a release.','', 'Decision: '+decision['reason'],'',
'## Change and correctness','',
'The small-ring breadth-first search uses one inline queue of up to 32 atom indices as both the visited sequence and the level frontiers. Captured begin/end positions preserve the exact depth boundary and adjacency traversal order. Before adding a 33rd visited atom, the candidate seeds a HashSet from the queue and uses hashed membership thereafter. Large neighborhoods therefore retain hashed lookup complexity. An exhausted frontier returns false without unused empty iterations. No chemical predicate, coordinate arithmetic, tolerance, input filtering, public API, dependency or unsafe code is changed. The preceding inactive E/Z guard and every fallback classification/diagnostic remain unchanged.','',
'The immediate baseline contains all thirteen accepted source candidates, including the inactive E/Z topology guard. Rejected prototypes are excluded. Its exact native executable is reused only after all 10,060 current baseline source/input hashes and its binary hash match the accepted source13 build. Baseline build/test logs and origin are retained. The candidate is built freshly in a dedicated target; no timing samples are reused or pooled. Every build and CPU sample completes before formal timing.','',
f"The independent BFS oracle freezes the original traversal; only its return values add a visited-count diagnostic. Every bond of all 14,072 existing pinned inputs is checked at 15 thresholds (0 through 12, 33 and 65), totaling {oracle['queries']:,} queries. Results match exactly: {oracle['inline_cases']:,} queries stay within 32 visited atoms and {oracle['hashed_cases']:,} exceed that boundary, exercising both storage paths. True and false outcomes both occur. Existing 506,592 frozen public E/Z comparisons also pass again with coordinate failures, unspecified-bond masks, annotations, parsed/Kekule representations and repeated calls. Independent native snapshots and rewritten bytes match on all 2,064 fixed MOL inputs and all 64 continued chains. These oracle comparisons reuse existing licensed input records. A separate project-authored synthetic diagnostic is described below.",'',
f"Scoped release Rust tests: {counts['passed']:,} passed, {counts['failed']} failed, {counts['ignored']} existing ignored. Full MOL and perception libraries plus six integration suites cover E/Z/wedges, 3D, square planar, semantic preservation and source-ID lookup. Ten runner gates, release all-target Clippy for MOL/perception, formatting and unchanged refactor budget (487 production files, 8,063 functions, 18 file allowances, 88 function allowances) pass. Complete logs remain in the packet.",'',
'Fresh source13 CPU profiles over the fixed 1,000-input batches are retained. The harness also performs untimed preparation; samples diagnose the executed read path, not pure timed latency. Inclusive counts overlap and are not additive; LTO can hide symbols. Sampling confirms ring-search work as a target, not a predicted speedup.','',
'## Protocol','',
'Fixed unfiltered profiles: 1,000 records at offsets 0 and 5,000, plus 64 controlled linear chains with lengths 16 through 1,024. One exposed development corpus on one macOS arm64 host, without a sealed holdout or toolkit-wide coverage. Open Babel 3.2.1 libraries/plugins and RDKit 2026.09.1 observer are pinned by hashes. Each profile uses 21 rotating three-arm fresh-process blocks, five batch repeats and all parse/prepared-write/roundtrip modes.','',
'Ratios are medians of per-block paired reference/candidate elapsed ratios; above 1 is faster. Milliseconds are separate marginal medians, so their quotient need not equal the paired ratio. Brackets are paired median-ratio 95% intervals using 10,000 deterministic bootstrap resamples. A strict win requires semantic agreement, CI lower > 1 and all 21 blocks faster. Original and follow-up observations remain separate, without pooling or substitution. Movement in unchanged paths may include host noise and binary layout; it is not proof of a new algorithm.','',
'Timers exclude startup, loading/framing and plugin discovery. Preparation retains a full parsed batch and performs one warm write. Parse evaluates new batches with default diagnostics, prepared-write creates strings from prepared records, roundtrip parses and writes fresh records. RSS includes preparation/plugins and is whole-process high-water usage. Output, source, input, executable and native library hashes are checked.','',
'An external checkpoint/artifact wrapper reserves 16 MiB of owned storage, releases it on ENOSPC and retries writing the same computed report for at most 60 seconds. The frozen native timers, observer/bootstrap functions, selection policy and measured observations are unchanged. Logs retain any such event; collection-completion.json records interruptions if any.','',
'## Original full matrix','', '| Profile | Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks | Open Babel / candidate [95% CI]; faster blocks |','|---|---|---:|---:|---|---|']
for profile,j in reports.items():
 for mode,v in j['summary'].items():lines.append(f"| {profile} | {mode} | {ms(j,'baseline',mode):.3f} | {ms(j,'candidate',mode):.3f} | {fmt(v['baseline'])} | {fmt(v['openbabel'])} |")
lines+=['','## Semantic exclusions','',
'Diverse Open Babel lanes remain unranked whenever indexed graphs, stereo or rewritten fields differ. No records or stereo markers are removed to obtain a ranking. Controlled-chain agreement establishes only that checked workload. Source-only snapshots include indexed atoms/bonds, stereo order, exact coordinate bits, names/metadata and geometry/stereo diagnostics; rewrite bytes are also required to match.','', '| Profile | Open Babel parse graph mismatches | Parse rejected | Rewrite rejected | Rewrite graph mismatches | Unknown-stereo mismatches |','|---|---:|---:|---:|---:|---:|']
for profile,j in reports.items():
 q=j['semantics']['openbabel_parse'];w=j['semantics']['openbabel'];c=w.get('mismatch_counts',{});lines.append(f"| {profile} | {q['mismatches']} | {len(q['rejected_rows'])} | {len(w['rejected_rows'])} | {c.get('graph',0)} | {c.get('explicit_unknown_stereo',0)} |")
lines+=['','## Separate continuation diagnostic','',
'All numeric atom/bond rows of the 64 chains are physically continued after the id. Original and continued source snapshots/rewrites must be byte-identical in both arms. This source-only diagnostic has 21 alternating fresh-process blocks, five repeats and all modes. It does not rank Open Babel or replace an original profile.','', '| Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks |','|---|---:|---:|---|']
for mode,v in diagnostic['summary'].items():lines.append(f"| {mode} | {ms(diagnostic,'baseline',mode):.3f} | {ms(diagnostic,'candidate',mode):.3f} | {fmt(v)} |")
lines+=['','## Adaptive regression checks','',
'The rule saved before timing selects an entire profile for one separate 21-block source-only all-mode follow-up when any initial source-relative CI is wholly below 1. Initial outcomes remain; newly appearing follow-up negatives do not cause recursive retries.','', '```json',json.dumps(selection,indent=2),'```','']
for label,j in focused.items():
 lines += [f'### {label}','','| Lane | Baseline / candidate [95% CI]; faster blocks |','|---|---:|']
 for mode,v in j['summary'].items():lines.append('| '+mode+' | '+fmt(v)+' |')
 lines.append('')
large=json.loads((r/'large-neighborhood21.json').read_text());assert large['complete'] and len(large['blocks'])==21
large_records=json.loads((r/'large-neighborhood-records.json').read_text());assert len(large_records)==64 and sum(v['uses_hash_overflow'] for v in large_records)==32
lines += ['## Additional large-neighborhood diagnostic','',
'This exploratory source-only diagnostic was declared after gains were observed in the two initial fixed slices, before generating or timing its inputs. It neither replaces nor pools with any original matrix and does not compare Open Babel. All outcomes remain, without recursive follow-up.','',
'The project-authored recipe creates 64 valid acyclic carbon graphs: a central C=C with equal balanced binary subtrees at both ends, eight depths and eight fixed atom permutations, spanning 6 through 1,022 atoms. All graphs sanitize in pinned RDKit; no record is removed. Production ring search visits more than 32 atoms on 32 records and stays inline on 32. Full provenance, license, minimization, exact recipe and pinned oracle are in large-neighborhood-corpus-manifest.json.','',
'Before timing, native snapshots and rewrite bytes match between both source arms for every record. RDKit verifies the generated versus native canonical/rewritten graph, name and exact coordinates. The diagnostic uses 21 alternating two-arm fresh-process blocks, five repeats and all modes through the unchanged source-only driver. It is exposed synthetic evidence, not a sealed holdout or a guarantee for every large graph.','',
'| Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks |','|---|---:|---:|---|']
for mode,v in large['summary'].items():lines.append(f"| {mode} | {ms(large,'baseline',mode):.3f} | {ms(large,'candidate',mode):.3f} | {fmt(v)} |")
lines.append('')
lines+=['## Reproduction and status','',
'[The evidence packet]('+pre+'-artifacts/README.md) retains exact measured binaries, every matrix/snapshot/rewrite, failures if any, build/test/profile logs, full frozen source/input archive and patch. Archive readback verifies all '+str(len(hashes))+' candidate hashes and reconstructs all '+str(len(baseline))+' baseline hashes. An independent tiny Git repository checks patch reconstruction. Source reconstruction is not bit-reproducibility of rebuilt binaries.','',
'Last verified shared Git base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`; refs were not refreshed. Source version is 1.1.4, without merge, hosted CI or publication. Immediate-baseline improvements do not establish recovery of every earlier regression; older writer-induced read regressions, Python Morgan/MACCS observations and other unresolved lanes remain separate evidence.','']
report=d/'benchmarks'/(pre+'.md');report.write_text('\n'.join(lines))
(a/'README.md').write_text("""# Inline ring-search evidence

Read decision.json, measurement-plan.json, source-manifest.json, all paired21 matrices, continued21.json and adaptive-selection.json. Focused reports are separate diagnostics. Original inputs and semantic exclusions remain.

Extract candidate-build-source.tar.gz: candidate/ is the complete build/input source, baseline-overrides/ contains pre-change files, corpora/ contains exact fixed streams. Reconstruct baseline by applying overrides and removing the added test listed in the manifest. verification.json checks every source/input hash. candidate.patch reconstructs changes independently. Existing input provenance is in chemistry-speed-followups-v1.manifest.json; Main fixed corpora are unchanged; the separate synthetic workload has its own provenance manifest and plan.

The two native macOS arm64 executables are exactly measured. Rebuild with the recorded Rust toolchain and cargo build --release --locked -p chematic-mol --example openbabel_v3000_throughput, opt3/fat-LTO/codegen-units1. Bit identity of a rebuild is not promised. Open Babel harness/library/plugin hashes are pinned and C++ harness source is archived.

Run scripts/bench_openbabel_v3000_throughput.py with --baseline, --candidate, --openbabel, --corpus, --records 1000, --blocks 21, --repeats 5, --source-manifest and --output in the pinned observer environment. Controlled chains use --records 64. Darwin time -l requires OS-statistics access. Driver paths are host-specific. Exact gzip streams are verified; artifact-sha256.json indexes every other packet file. Read the adoption decision before using an archived prototype.
""")
spec=importlib.util.spec_from_file_location('names',s/'scripts/check_product_names.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);assert not mod.scan_paths(d,[report,a/'README.md']+([d/'CHANGELOG.md'] if adopt else []))
index={p.relative_to(a).as_posix():sha(p) for p in sorted(a.rglob('*')) if p.is_file() and p.name!='artifact-sha256.json'};(a/'artifact-sha256.json').write_text(json.dumps(index,indent=2)+'\n');assert all(sha(a/n)==h for n,h in index.items());(r/'delivery-sha256.json').write_text(json.dumps({p.relative_to(d).as_posix():sha(p) for p in sorted(d.rglob('*')) if p.is_file()},indent=2)+'\n');print('Packaged',len(index),'artifacts; native tests',counts,'adopted',adopt)
