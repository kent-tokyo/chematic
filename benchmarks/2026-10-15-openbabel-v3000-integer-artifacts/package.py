from pathlib import Path
import json,hashlib,gzip,tarfile,shutil,statistics,re,importlib.util
r=Path('/private/tmp/chematic-openbabel-integer-20261015');s=Path('/private/tmp/chematic-openbabel-integer-source-20261015');d=r/'delivery';pre='2026-10-15-openbabel-v3000-integer';a=d/'benchmarks'/(pre+'-artifacts');a.mkdir(parents=True,exist_ok=True);sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();m=json.loads((r/'source-manifest.json').read_text());assert all(sha(s/n)==h for n,h in m['frozen_build_source_sha256'].items());assert 'ALL VALIDATION PASSED' in (r/'validation-pipeline.log').read_text()
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
for folder in list(r.glob('paired21-*.artifacts'))+[r/'continuation.artifacts']:
 for p in sorted(folder.iterdir()):
  if p.is_file():
   q=a/folder.name/(p.name+'.gz');q.parent.mkdir(exist_ok=True)
   with q.open('wb') as out,gzip.GzipFile(filename='',mode='wb',fileobj=out,mtime=0) as z:z.write(p.read_bytes())
   assert hashlib.sha256(gzip.decompress(q.read_bytes())).hexdigest()==sha(p)
old=Path('/private/tmp/chematic-openbabel-v3000-20261014');oldfiles={'CHANGELOG.md':'baseline-CHANGELOG.md','crates/chematic-mol/src/mol2000.rs':'baseline-mol2000.rs','crates/chematic-mol/src/mol2000/writer.rs':'baseline-v2000-writer.rs','crates/chematic-mol/src/mol3000/writer.rs':'baseline-v3000-writer.rs'}
with tarfile.open(a/'candidate-build-source.tar.gz','w:gz') as tar:
 for n in m['frozen_build_source_sha256']:tar.add(s/n,arcname='candidate/'+n,recursive=False)
 for n,f in oldfiles.items():tar.add(r/f,arcname='baseline-overrides/'+n,recursive=False)
 for profile in reports:tar.add(old/f'corpus-{profile}.molblocks',arcname='corpora/'+profile+'.molblocks',recursive=False)
 tar.add(r/'corpus-continued-chains.molblocks',arcname='corpora/continued-chains.molblocks',recursive=False)
with tarfile.open(a/'candidate-build-source.tar.gz') as tar:
 hashes={member.name.removeprefix('candidate/'):hashlib.sha256(tar.extractfile(member).read()).hexdigest() for member in tar.getmembers() if member.name.startswith('candidate/')};assert hashes==m['frozen_build_source_sha256'];baseline=hashes.copy()
 for n in changed:baseline[n]=hashlib.sha256(tar.extractfile('baseline-overrides/'+n).read()).hexdigest()
 for n in added:del baseline[n]
 assert baseline==m['baseline_source_files_sha256']
counts={'passed':0,'failed':0,'ignored':0}
for name in ['native-library-tests.log','native-integration-tests.log']:
 for values in re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored;', (r/name).read_text()):
  for k,v in zip(counts,values):counts[k]+=int(v)
assert counts['passed']>=849 and counts['failed']==counts['ignored']==0,counts
(a/'verification.json').write_text(json.dumps({'source_archive_exact':True,'baseline_reconstruction_exact':True,'candidate_file_count':len(hashes),'baseline_file_count':len(baseline),'native_tests':counts,'source_parse_and_rewrite_equality':{profile:{'source_parse_snapshots_equal':j['source_parse_snapshots_equal'],'snapshot_sha256':j['snapshot_sha256'],'rewrite_sha256':j['rewrite_sha256']} for profile,j in reports.items()},'continued_input_equivalence':json.loads((r/'continuation-equivalence.json').read_text()),'adaptive_selection':selection},indent=2)+'\n')
def fmt(v):
 z=v['speedup'];return 'unranked (semantic mismatch)' if z is None else f"{z['median']:.3f} [{z['ci95'][0]:.3f}, {z['ci95'][1]:.3f}]; {z['winning_blocks']}/21"+(' **strict win**' if v.get('strict_speed_win') else '')
def ms(j,arm,mode):return statistics.median(b[arm][mode]['elapsed_ns'] for b in j['blocks'])/1e6
lines=['# V3000 direct integer writing — 2026-10-15','',('The production optimization is adopted in the attached source candidate.' if adopt else 'The production optimization is rejected; only its report and evidence are delivered.')+' No toolkit-wide superiority, current-main or release claim follows.','', 'Decision: '+decision['reason'],'',
'The V3000 2D writer appends atom/bond ids using the existing private decimal formatter with zero padding width; the existing three-column field table now also supplies its unpadded suffix. Bond-order codes, all currently 0 through 9, are appended as one decimal digit. Existing whitespace, element symbols, atom maps, properties, wedge orientation, CFG fields and floating-point formatting remain unchanged. Values above 999 retain the original general integer path. Three-column V2000 output remains unchanged. The formatter is exposed only within the crate. No external public API, dependency, unsafe code or tolerance is added. The rejected ASCII lexer is excluded and the reader is unchanged; parse timing movement is not attributed to a new reader algorithm.','',
'The immediate baseline contains all eleven previously accepted source candidates, including packed private aromatic-component facts. All rejected experiments remain excluded. The exact baseline11 native release/LTO executable is reused from the preceding rejected ASCII trial, where it was built before either edit. All 10058 baseline source hashes and its binary hash match; baseline-artifact-origin.json records the origin. Candidate is built freshly in a separate target. No timing observation from the ASCII trial is pooled. All builds, validation and cleanup of owned completed targets preceded timing. All 2,064 original records retain byte-identical snapshots and rewritten bytes between source arms; all 64 continued chains match their original records in both arms.','',
f"Scoped Rust tests: {counts['passed']:,} passed, zero failed/ignored. The full MOL library and six integration suites cover E/Z and wedge stereo, 3D, square planar, V3000 semantic preservation and source-ID lookup. The decimal formatter is checked against the standard integer formatter for 143,702 fields: widths 0 through 12, all values 0 through 1,024, decimal-power boundaries, maximum u32 values, and 10,000 deterministic u32 samples. Existing frozen V3000 writer tests and full raw corpus snapshots/rewrites check unchanged output, including maps, properties, stereo and coordinates. No new molecular input records are added. Ten benchmark runner gates, release all-target Clippy, formatting and the unchanged refactor budget pass.",'',
'## Protocol','',
'Fixed unfiltered development profiles: 1,000 records at offsets 0 and 5,000, plus 64 controlled linear chains of lengths 16 through 1,024. One macOS arm64 host; Open Babel 3.2.1 native libraries/plugins and RDKit 2026.09.1 observer pinned by hashes. Each profile has 21 rotating three-arm fresh-process blocks, five batch repeats, and all parse/prepared-write/roundtrip modes. No sealed holdout or toolkit-wide coverage.','',
'Ratios are comparator elapsed time / candidate elapsed time; above 1 is faster. Brackets are paired median-ratio 95% intervals from 10,000 deterministic bootstrap resamples. Strict wins require semantic agreement, CI lower > 1 and all 21 blocks faster. Every original and follow-up observation is retained without pooling or substitution. An external checkpoint/artifact-write wrapper retries ENOSPC using an owned 16MiB reserve, preserving computed reports in memory; source benchmark timing/observer/bootstrap functions and all native timers are unchanged. Driver source and any retry messages are retained.','',
'Timers exclude startup, input loading/framing and plugin discovery. Preparation retains a full parsed batch and includes one untimed warm write; parse evaluates new batches with default diagnostics, prepared-write creates output strings, roundtrip parses and writes fresh records. RSS is whole-process high-water usage including preparation/plugins; it is not per-operation allocation evidence. Source, input, native executable and library hashes are checked.','',
'## Original full matrix','', '| Profile | Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks | Open Babel / candidate [95% CI]; faster blocks |','|---|---|---:|---:|---|---|']
for profile,j in reports.items():
 for mode,v in j['summary'].items():lines.append(f"| {profile} | {mode} | {ms(j,'baseline',mode):.3f} | {ms(j,'candidate',mode):.3f} | {fmt(v['baseline'])} | {fmt(v['openbabel'])} |")
lines+=['','## Semantic exclusions','', 'Diverse Open Babel lanes remain unranked when internal graphs, stereo or rewrite fields differ. No row or stereo marker is removed. Controlled-chain agreement establishes only that checked workload.','', '| Profile | Open Babel parse graph mismatches | Parse rejected | Rewrite rejected | Rewrite graph mismatches | Unknown-stereo mismatches |','|---|---:|---:|---:|---:|---:|']
for profile,j in reports.items():
 p=j['semantics']['openbabel_parse'];w=j['semantics']['openbabel'];c=w.get('mismatch_counts',{});lines.append(f"| {profile} | {p['mismatches']} | {len(p['rejected_rows'])} | {len(w['rejected_rows'])} | {c.get('graph',0)} | {c.get('explicit_unknown_stereo',0)} |")
lines+=['','## Separate continuation diagnostic','', 'Every numeric atom/bond row of all 64 controlled chains is split after its id into continued physical lines. All snapshots and rewrites match the original chains in both arms. The source-only diagnostic uses 21 alternating two-arm fresh-process blocks, five repeats and all modes, without Open Babel ranking or RSS measurement. It never replaces an original profile.','', '| Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks |','|---|---:|---:|---|']
for mode,v in diagnostic['summary'].items():lines.append(f"| {mode} | {ms(diagnostic,'baseline',mode):.3f} | {ms(diagnostic,'candidate',mode):.3f} | {fmt(v)} |")
lines+=['','## Adaptive regression checks','', 'The predeclared rule selects an entire profile for one separate 21-block source-only follow-up in all modes if any initial source-relative CI is wholly below 1. Initial and follow-up outcomes remain separate. Newly appearing negatives do not trigger recursive retries.','', '```json',json.dumps(selection,indent=2),'```','']
for label,j in focused.items():
 lines += [f'### {label}','','| Lane | Baseline / candidate [95% CI]; faster blocks |','|---|---:|']
 for mode,v in j['summary'].items():lines.append('| '+mode+' | '+fmt(v)+' |')
 lines.append('')
lines+=['## Reproduction and status','', '[The evidence packet]('+pre+'-artifacts/README.md) retains exact measured binaries, all matrices/snapshots/rewrites, failures if any, test/build logs, full frozen source/input archive and exact patch. Archive readback verifies all '+str(len(hashes))+' candidate hashes and reconstructs all '+str(len(baseline))+' baseline hashes. A separate tiny Git repository verifies patch reconstruction. Source reconstruction is not bit-reproducibility of rebuilt binaries.','', 'Last verified shared Git base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`; refs were not refreshed. These are version 1.1.4 unmerged source candidates, with no hosted CI or publication. Prior unrelated source-candidate regressions against older baselines remain unresolved by this incremental comparison.','']
report=d/'benchmarks'/(pre+'.md');report.write_text('\n'.join(lines))
(a/'README.md').write_text('''# V3000 integer writer evidence

Read decision.json, measurement-plan.json, source-manifest.json, all paired21 reports, continued21.json and adaptive-selection.json. Focused reports remain separate diagnostics. All original inputs and semantic exclusions are retained.

Extract candidate-build-source.tar.gz: candidate/ contains the full build/input source, baseline-overrides/ contains original changed files, corpora/ contains fixed exact benchmark streams. Reconstruct baseline by applying overrides and deleting added files listed in the manifest. verification.json checks all source hashes. candidate.patch reconstructs the candidate in an independent tiny Git repository. No new molecular corpus rows are added; original corpus provenance remains in chemistry-speed-followups-v1.manifest.json.

The two native macOS arm64 executables are exact measured binaries. Rebuild with the recorded Rust toolchain and cargo build --release --locked -p chematic-mol --example openbabel_v3000_throughput, release opt3/LTO/codegen-units1. Bit identity of a rebuild is not promised. Open Babel native harness/library/plugin paths and hashes are pinned; its C++ harness source is archived.

Use extracted scripts/bench_openbabel_v3000_throughput.py with --baseline, --candidate, --openbabel, --corpus, --records 1000, --blocks 21, --repeats 5, --source-manifest and --output in the pinned observer environment. Use --records 64 for controlled chains. Darwin time -l requires OS-statistics access. Driver paths are host-specific. Gzip streams are verified against every raw file. artifact-sha256.json indexes other packet files. These are unmerged version 1.1.4 source candidates.
''')
spec=importlib.util.spec_from_file_location('names',s/'scripts/check_product_names.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);assert not mod.scan_paths(d,[report,a/'README.md']+([d/'CHANGELOG.md'] if adopt else []))
index={p.relative_to(a).as_posix():sha(p) for p in sorted(a.rglob('*')) if p.is_file() and p.name!='artifact-sha256.json'};(a/'artifact-sha256.json').write_text(json.dumps(index,indent=2)+'\n');assert all(sha(a/n)==h for n,h in index.items());(r/'delivery-sha256.json').write_text(json.dumps({p.relative_to(d).as_posix():sha(p) for p in sorted(d.rglob('*')) if p.is_file()},indent=2)+'\n');print('Packaged',len(index),'artifacts; native tests',counts,'adopted',adopt)
