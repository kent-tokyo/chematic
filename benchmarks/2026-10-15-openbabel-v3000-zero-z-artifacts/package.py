from pathlib import Path
import json,hashlib,gzip,tarfile,shutil,statistics,re,importlib.util
r=Path('/private/tmp/chematic-openbabel-zero-z-20261015');s=Path('/private/tmp/chematic-openbabel-zero-z-source-20261015');d=r/'delivery';pre='2026-10-15-openbabel-v3000-zero-z';a=d/'benchmarks'/(pre+'-artifacts');a.mkdir(parents=True);sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();m=json.loads((r/'source-manifest.json').read_text());assert all(sha(s/n)==h for n,h in m['frozen_build_source_sha256'].items());assert 'ALL VALIDATION PASSED' in (r/'validation-pipeline.log').read_text()
decision=json.loads((r/'decision.json').read_text());adopt=decision['adopt_source'];changed=m['only_changed_existing_source'];added=m['added_source_files']+m['added_test_files']
if adopt:
 for n in changed+added:
  p=d/n;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(s/n,p)
reports={profile:json.loads((r/f'paired21-{profile}.json').read_text()) for profile in ['0','5000','chains']};diagnostic=json.loads((r/'continued21.json').read_text());selection=json.loads((r/'adaptive-selection.json').read_text());focused={p.stem:json.loads(p.read_text()) for p in r.glob('focused21-*.json')}
for profile,x in reports.items():
 assert x['complete'] and len(x['blocks'])==21 and x['source_parse_snapshots_equal'] and x['source_manifest']==m
 assert (r/f'paired21-{profile}.artifacts/baseline.molblocks').read_bytes()==(r/f'paired21-{profile}.artifacts/candidate.molblocks').read_bytes()
for j in [diagnostic,*focused.values()]:assert j['complete'] and len(j['blocks'])==21
for p in sorted(r.iterdir()):
 if p.is_file() and p.suffix in ['.json','.patch','.log','.py']:shutil.copyfile(p,a/p.name)
for arm in ['baseline','candidate']:
 shutil.copyfile(r/(arm+'-v3000'),a/(arm+'-v3000'));(a/(arm+'-v3000')).chmod(0o755)
for folder in list(r.glob('paired21-*.artifacts'))+[r/'continuation.artifacts']:
 for p in sorted(folder.iterdir()):
  if p.is_file():
   q=a/folder.name/(p.name+'.gz');q.parent.mkdir(exist_ok=True)
   with q.open('wb') as out,gzip.GzipFile(filename='',mode='wb',fileobj=out,mtime=0) as z:z.write(p.read_bytes())
   assert hashlib.sha256(gzip.decompress(q.read_bytes())).hexdigest()==sha(p)
old=Path('/private/tmp/chematic-openbabel-v3000-20261014')
with tarfile.open(a/'candidate-build-source.tar.gz','w:gz') as tar:
 for n in m['frozen_build_source_sha256']:tar.add(s/n,arcname='candidate/'+n,recursive=False)
 for n,f in [('CHANGELOG.md','baseline-CHANGELOG.md'),('crates/chematic-mol/src/mol3000.rs','baseline-mol3000.rs')]:tar.add(r/f,arcname='baseline-overrides/'+n,recursive=False)
 for profile in reports:tar.add(old/f'corpus-{profile}.molblocks',arcname='corpora/'+profile+'.molblocks',recursive=False)
 tar.add(r/'corpus-continued-chains.molblocks',arcname='corpora/continued-chains.molblocks',recursive=False)
with tarfile.open(a/'candidate-build-source.tar.gz') as tar:
 hashes={member.name.removeprefix('candidate/'):hashlib.sha256(tar.extractfile(member).read()).hexdigest() for member in tar.getmembers() if member.name.startswith('candidate/')}
 assert hashes==m['frozen_build_source_sha256']
 baseline=hashes.copy()
 for n in changed:baseline[n]=hashlib.sha256(tar.extractfile('baseline-overrides/'+n).read()).hexdigest()
 for n in added:del baseline[n]
 assert baseline==m['baseline_source_files_sha256']
counts={'passed':0,'failed':0,'ignored':0}
for n in ['native-library-tests.log','native-integration-tests.log']:
 for values in re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored;', (r/n).read_text()):
  for k,v in zip(counts,values):counts[k]+=int(v)
assert counts=={'passed':844,'failed':0,'ignored':0},counts
(a/'verification.json').write_text(json.dumps({'source_archive_exact':True,'baseline_reconstruction_exact':True,'candidate_file_count':len(hashes),'baseline_file_count':len(baseline),'native_tests':counts,'source_parse_and_rewrite_equality':{profile:{'source_parse_snapshots_equal':x['source_parse_snapshots_equal'],'snapshot_sha256':x['snapshot_sha256'],'rewrite_sha256':x['rewrite_sha256']} for profile,x in reports.items()},'continued_input_equivalence':json.loads((r/'continuation-equivalence.json').read_text()),'adaptive_selection':selection},indent=2)+'\n')
def fmt(v):
 z=v['speedup'];return 'unranked (semantic mismatch)' if z is None else f"{z['median']:.3f} [{z['ci95'][0]:.3f}, {z['ci95'][1]:.3f}]; {z['winning_blocks']}/21"+(' **strict win**' if v.get('strict_speed_win') else '')
def ms(j,arm,mode):return statistics.median(b[arm][mode]['elapsed_ns'] for b in j['blocks'])/1e6
lines=['# V3000 positive-zero z parsing — 2026-10-15','',
('The production optimization is adopted in the attached source candidate.' if adopt else 'The production optimization is rejected for adoption; only its experiment report and evidence are delivered.')+' No toolkit-wide superiority, current-main or release claim follows from this run.','',
'Decision: '+decision['reason'],'',
'The V3000 coordinate reader now handles the exact finite positive-zero z lexemes `0` and `0.0000` directly, avoiding the general float parser for those common 2D values. Every other spelling, including signed zero, still uses the original parser and finite-value check. x and y are evaluated first in the original order, preserving the first typed error, its axis, input text and source line. No public API, dependency, tolerance or unsafe code is added.','',
'The immediate baseline contains all seven preceding accepted candidates, including the latest aromaticity component reuse. It excludes the rejected direct-line experiment. Both native release/LTO executables were built freshly before and after this change; source, binaries, corpus and Open Babel libraries are hashed. All original 2,064 records retain byte-identical parser snapshots and rewrites between the two source arms, and all 64 physical-continuation variants match the original chains in both arms.','',
f"Scoped Rust validation: {counts['passed']:,} passed, no failures or ignored tests. This includes the full MOL library and four integration suites for 2D E/Z, wedge stereo, 3D conformers and square-planar I/O. Runner semantic gates: 10 passed. Release all-target Clippy, formatting and refactor checks passed. The frozen original coordinate parser is compared on 55,296 Cartesian lexeme/source-line cases and 10,000 deterministic f64 fallback cases. Exact tuple bits and complete typed-error representations match, including negative-zero bits and x/y/z failure precedence. The original numeric recipe, license, minimization, oracle source hash and normalized payload hash are recorded in the new zero-z fixture manifest. No new molecular records were added.",'',
'## Protocol','',
'All fixed profiles are retained without filtering: 1,000 records each at offsets 0 and 5,000, and 64 linear chains of lengths 16 through 1,024 in steps of 16. This is an exposed development corpus on one macOS arm64 host, without a sealed holdout. Open Babel 3.2.1 native library/plugins and RDKit 2026.09.1 observer are pinned. Each original profile has 21 rotating three-arm fresh-process blocks and five full-batch repeats per mode. All builds and validation precede timing; profiles and diagnostics run sequentially. Only this run owned, completed build target was removed after retaining both measured binaries.','',
'Ratios are comparator time / candidate time; above 1 means faster. Brackets are paired median-ratio 95% intervals from 10,000 deterministic bootstrap resamples. Strict wins require semantic agreement, CI lower > 1 and all 21 blocks faster. Every record, block, mismatch, failure and initial/adaptive negative is retained. The writer is unchanged; writer timing movement is not attributed to a new writer algorithm.','',
'Timers exclude startup, input loading/framing and plugin discovery. The native harness prepares a parsed batch and performs an untimed warm write before timing. Parse retains the complete diagnostic batch, prepared writing creates output strings, and roundtrip parses and writes fresh records. RSS is whole-process high-water usage including preparation and plugins; it is not per-operation allocation evidence. The general-parser bypass follows from the implementation; no per-operation allocation counter is measured.','',
'## Original full matrix','', '| Profile | Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks | Open Babel / candidate [95% CI]; faster blocks |','|---|---|---:|---:|---|---|']
for profile,j in reports.items():
 for mode,v in j['summary'].items():lines.append(f"| {profile} | {mode} | {ms(j,'baseline',mode):.3f} | {ms(j,'candidate',mode):.3f} | {fmt(v['baseline'])} | {fmt(v['openbabel'])} |")
lines+=['','## Semantic exclusions','',
'Open Babel disagreements remain unranked. No stereo marker is removed and no row is dropped to obtain a speed result. Diverse inputs include graph mismatches, explicit unknown-double-stereo loss and rejected rewrites; those lanes cannot support equivalent-output speed wins. Controlled-chain agreement does not establish broad chemical correctness or superiority.','', '| Profile | Open Babel parse graph mismatches | Parse rejected | Rewrite rejected | Rewrite graph mismatches | Unknown-stereo mismatches |','|---|---:|---:|---:|---:|---:|']
for profile,j in reports.items():
 p=j['semantics']['openbabel_parse'];w=j['semantics']['openbabel'];c=w.get('mismatch_counts',{});lines.append(f"| {profile} | {p['mismatches']} | {len(p['rejected_rows'])} | {len(w['rejected_rows'])} | {c.get('graph',0)} | {c.get('explicit_unknown_stereo',0)} |")
lines+=['','## Separate continuation diagnostic','',
'The predeclared source-only diagnostic splits every numeric V30 atom and bond row of all controlled chains after its id into continued physical lines. All original and continued snapshots and rewrites agree between both source arms. It uses 21 alternating two-arm fresh-process blocks, all modes and five repeats; it has no Open Babel ranking or RSS measurement and never replaces an original observation.','', '| Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks |','|---|---:|---:|---|']
for mode,v in diagnostic['summary'].items():lines.append(f"| {mode} | {ms(diagnostic,'baseline',mode):.3f} | {ms(diagnostic,'candidate',mode):.3f} | {fmt(v)} |")
lines+=['','## Adaptive regression checks','',
'Any initial baseline-relative interval wholly below 1 selects that profile for one separate 21-block source-only follow-up in all three modes. All initial results remain above; follow-ups are not pooled or substituted, and new follow-up negatives do not cause recursive retries.','', '```json',json.dumps(selection,indent=2),'```','']
for label,j in focused.items():
 lines += [f'### {label}','','| Lane | Baseline / candidate [95% CI]; faster blocks |','|---|---:|']
 for mode,v in j['summary'].items():lines.append('| '+mode+' | '+fmt(v)+' |')
 lines.append('')
lines+=['## Reproduction and status','',
'[The evidence packet]('+pre+'-artifacts/README.md) contains both measured source executables, every initial/follow-up result and semantic output, full build/test logs, source manifests, fixture recipe/provenance, the exact patch and a complete build/input archive. Every '+str(len(hashes))+' candidate archive hash and '+str(len(baseline))+' reconstructed baseline hash was verified. The patch also reconstructs the changed files in an independent tiny Git repository. The preceding [flat-geometry allocation experiment](2026-10-15-openbabel-v3000-flat-geometry.md) was rejected because none of its original-profile parse intervals established a gain. Its source, outcomes and logs remain archived separately. Both its point-allocation changes and shared geometry helper are excluded from this trial; measurements are not pooled.','',
'Last verified shared Git base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`. Shared refs were not refreshed; hashes identify version 1.1.4 source candidates, not current main. This run has no hosted CI, merge or publication. Previous source-candidate negatives in unrelated Python or native lanes are not resolved by an incremental comparison. The earlier rejected direct-line candidate remains archived and excluded.','']
report=d/'benchmarks'/(pre+'.md');report.write_text('\n'.join(lines))
(a/'README.md').write_text('''# Flat-coordinate geometry speed evidence

Start with `decision.json`, `measurement-plan.json`, `source-manifest.json`, all three `paired21-*.json` results, `continued21.json` and `adaptive-selection.json`. Any focused follow-ups are separate diagnostics and do not replace initial results.

The archive contains the complete frozen candidate build/input tree under `candidate/`, two baseline source overrides under `baseline-overrides/` and the unchanged exact benchmark inputs under `corpora/`. Reconstruct the baseline by applying the overrides to a copy and deleting all added files listed in the source manifest. `verification.json` validates every archived candidate hash and reconstructed baseline hash. `candidate.patch` reconstructs changed files exactly. The numeric recipe and fixture manifest record original provenance, license, minimization, pinned source oracle and payload hash.

`baseline-v3000` and `candidate-v3000` are the exact measured macOS arm64 native executables. Open Babel harness/libraries/plugins remain at the pinned manifest paths; the archived C++ harness source reproduces that harness against the recorded installation. Build Rust executables with the recorded toolchain and `cargo build --release --locked -p chematic-mol --example openbabel_v3000_throughput`. Rebuilt binary bit identity is not promised.

Run archived `scripts/bench_openbabel_v3000_throughput.py` with `--baseline`, `--candidate`, `--openbabel`, `--corpus`, `--records 1000`, `--blocks 21`, `--repeats 5`, `--source-manifest` and `--output` in the pinned RDKit observer environment. Use `--records 64` for chains. Darwin `/usr/bin/time -l` needs OS-statistics access. Driver script working paths are host-specific and must be adapted on another host. Original snapshots/rewrite streams are gzip-compressed and read-back verified. All initial failures and successful validation logs are retained. `artifact-sha256.json` indexes all other packet files. These are unmerged version 1.1.4 source candidates, not released or toolkit-wide proof.
''')
spec=importlib.util.spec_from_file_location('names',s/'scripts/check_product_names.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);assert not mod.scan_paths(d,[report,a/'README.md']+([d/'CHANGELOG.md'] if adopt else []))
index={p.relative_to(a).as_posix():sha(p) for p in sorted(a.rglob('*')) if p.is_file()};(a/'artifact-sha256.json').write_text(json.dumps(index,indent=2)+'\n');(r/'delivery-sha256.json').write_text(json.dumps({p.relative_to(d).as_posix():sha(p) for p in sorted(d.rglob('*')) if p.is_file()},indent=2)+'\n')
print('Packaged',len(index),'artifacts; native tests',counts,'adopted',adopt)
