from pathlib import Path
import json,hashlib,gzip,tarfile,shutil,statistics,re,importlib.util
r=Path('/private/tmp/chematic-openbabel-lines-20261015');s=Path('/private/tmp/chematic-openbabel-lines-source-20261015');d=r/'delivery';pre='2026-10-15-openbabel-v3000-direct-lines';a=d/'benchmarks'/(pre+'-artifacts');a.mkdir(parents=True);sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();m=json.loads((r/'source-manifest.json').read_text());assert all(sha(s/n)==h for n,h in m['frozen_build_source_sha256'].items());assert 'ALL VALIDATION PASSED' in (r/'validation-final.log').read_text()
changed=m['only_changed_existing_source'];added=m['added_source_files']+m['added_test_files']
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
results=[tuple(map(int,x)) for x in re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored', (r/'native-tests-retry.log').read_text())];assert len(results)>10 and sum(x[1] for x in results)==0
counts={'passed':sum(x[0] for x in results),'failed':sum(x[1] for x in results),'ignored':sum(x[2] for x in results),'runner_gate_tests':10}
(a/'verification.json').write_text(json.dumps({'native_tests':counts,'patch_reconstruction_exact':True,'source_archive_sha256':sha(a/'candidate-build-source.tar.gz'),'archive_candidate_hashes_exact':True,'baseline_reconstruction_hashes_exact':True,'candidate_file_count':len(hashes),'baseline_file_count':len(baseline),'source_parse_and_rewrite_equality':{profile:{'source_parse_snapshots_equal':x['source_parse_snapshots_equal'],'snapshot_sha256':x['snapshot_sha256'],'rewrite_sha256':x['rewrite_sha256']} for profile,x in reports.items()},'continued_input_equivalence':json.loads((r/'continuation-equivalence.json').read_text()),'adaptive_selection':selection},indent=2)+'\n')
def fmt(v):
 z=v['speedup'];return 'unranked (semantic mismatch)' if z is None else f"{z['median']:.3f} [{z['ci95'][0]:.3f}, {z['ci95'][1]:.3f}]; {z['winning_blocks']}/21"+(' **strict win**' if v.get('strict_speed_win') else '')
def ms(j,arm,mode):return statistics.median(b[arm][mode]['elapsed_ns'] for b in j['blocks'])/1e6
lines=['# V3000 direct payload consumption — 2026-10-15','', 'Rejected for adoption: the direct-line candidate did not improve original-profile parsing and repeated source-only follow-ups reproduced writer regressions of roughly 3–4%. All results and the source patch are retained. No production source or changelog from this candidate is applied to the attached worktree.','',
'Ordinary V3000 payloads now go directly from validated physical lines into the parser and reused token buffer. The logical-line vector is deferred until the first continuation; that entire suffix uses the unchanged continuation folder and a reused suffix token buffer. This removes the intermediate logical-line array on inputs without continuations. Physical-line validation still inspects the whole input before parsing, preserving error precedence. Unicode token splitting, source line numbers, metadata, diagnostics and chemical interpretation are unchanged. This is not a constant-memory file-streaming parser.','',
'Unmerged source version 1.1.4 candidate. The immediate baseline includes all five previous candidates: cycle-root pruning, V3000 writer buffers, shared Kekulé descriptor views, V3000 parser token reuse and ring-component reuse. Both native release/LTO binaries were freshly built; the baseline source was verified before editing. Last verified Git base is `cac480fb23cd8adeec6e5e9bed1f18677f87b392`; shared Git refs were not refreshed after earlier stalls. Source and binary hashes identify this experiment, not current main, a release or a registry artifact.','',
'All 2,064 original fixed records retain byte-identical baseline/candidate internal snapshots and rewrites. All 64 continued-chain variants also retain the original-chain snapshots and rewrite bytes in both source arms. The continuation diagnostic transforms physical framing only; it is not an additional diverse-chemistry or Open Babel ranking claim.','',
f"Full MOL tests: {counts['passed']:,} passed, {counts['failed']} failed, {counts['ignored']} ignored. Runner semantic-gate tests: 10 passed. Release all-target MOL Clippy, workspace formatting and the refactor budget passed. After Clippy identified one redundant borrow in an added test, that style-only fix passed all five added tests again; the full 1,018-test run preceded the fix and its log is retained. Frozen-reader differential tests compare complete graphs, coordinate bits, metadata, stereo/group fields, 2D/3D diagnostics, rewrites and displayed errors. They cover every continuation split position, malformed/dangling continuations, first-physical-line error precedence, CRLF/no final newline, Unicode whitespace, malformed numbers and duplicate/opaque attributes.",'',
'## Protocol','',
'Unchanged and unfiltered fixed profiles: 1,000 records each at offsets 0 and 5,000, and all 64 linear chains of lengths 16, 32, ..., 1,024. Exposed development corpora on one macOS arm64 host; no sealed holdout or toolkit-wide claim. Open Babel 3.2.1 native C++ library/plugins and RDKit 2026.09.1 observer are pinned. Each profile has 21 rotating three-arm fresh-process blocks, five full-batch repeats per mode. Builds/tests precede all timing; profiles, diagnostics and any adaptive follow-ups run sequentially. Source, binaries, libraries and inputs are verified before/after measurement.','',
'Native API timers exclude startup, input loading, delimiter splitting and plugin discovery. The unchanged harness prepares a parsed batch and performs an untimed warm write even before parse timing. Parse retains the full diagnostic batch; prepared writing allocates strings after one warm write; roundtrip actually parses and writes a fresh batch. RSS is the whole-process high-water mark including preparation/plugins/allocator state, not per-operation allocation.','',
'Ratios are comparator elapsed time / candidate elapsed time; greater than 1 means faster. Brackets show paired median-ratio 95% intervals from 10,000 deterministic bootstrap resamples. A strict speed win requires semantic agreement, CI lower > 1 and all 21 blocks faster. All input records, blocks, errors, mismatches and output markers are retained. Intervals crossing 1 do not establish an improvement. The writer is unchanged; writer movement is recorded without attributing it to a new writer algorithm.','',
'## Original full matrix','', '| Profile | Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks | Open Babel / candidate [95% CI]; faster blocks |','|---|---|---:|---:|---|---|']
for profile,j in reports.items():
 for mode,v in j['summary'].items():lines.append(f"| {profile} | {mode} | {ms(j,'baseline',mode):.3f} | {ms(j,'candidate',mode):.3f} | {fmt(v['baseline'])} | {fmt(v['openbabel'])} |")
lines+=['','## Semantic exclusions','',
'Open Babel disagreements remain unranked. No marker was removed or input filtered to make its results pass. The diverse-profile original observations include unknown-double-bond stereo loss and rejected rewrites; they are neither equivalent-output speed wins nor failures of the new source candidate.','', '| Profile | Open Babel parse graph mismatches | Parse rejected | Rewrite rejected | Rewrite graph mismatches | Unknown-stereo mismatches |','|---|---:|---:|---:|---:|---:|']
for profile,j in reports.items():
 p=j['semantics']['openbabel_parse'];w=j['semantics']['openbabel'];c=w.get('mismatch_counts',{});lines.append(f"| {profile} | {p['mismatches']} | {len(p['rejected_rows'])} | {len(w['rejected_rows'])} | {c.get('graph',0)} | {c.get('explicit_unknown_stereo',0)} |")
lines+=['','## Continuation fallback diagnostic','',
'The separately declared source-only diagnostic splits every numeric V30 atom and bond record of all 64 controlled chains after its id into two physical lines. Original and continued snapshots and rewrites match byte for byte in both source arms. It runs 21 alternating two-source-arm fresh-process blocks with the same five-repeat native timers and all three modes, without RSS or an Open Babel arm. It diagnoses the legacy suffix path rather than replacing any original result.','', '| Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks |','|---|---:|---:|---|']
for mode,v in diagnostic['summary'].items():lines.append(f"| {mode} | {ms(diagnostic,'baseline',mode):.3f} | {ms(diagnostic,'candidate',mode):.3f} | {fmt(v)} |")
lines+=['','## Adaptive negative observations','',
'Every original baseline-relative interval wholly below 1 selects that profile for one separate source-only 21-block follow-up, retaining all three modes. The same rule applies to the initial continuation diagnostic. No pooling, substitution of original results or recursive retries. Newly appearing follow-up negatives remain unresolved observations.','', '```json',json.dumps(selection,indent=2),'```','']
for label,j in focused.items():
 lines += [f'### {label}','', '| Lane | Baseline / candidate [95% CI]; faster blocks |','|---|---|']
 for mode,v in j['summary'].items():lines.append(f'| {mode} | {fmt(v)} |')
 lines+=['']
lines+=['## Limits and evidence','',
'Complete Open Babel-wide speed superiority is not established. Controlled-chain wins apply to that profile; diverse semantic exclusions, other formats, rich query/SGROUP chemistry, cold startup, filesystem/CLI/Python overhead, general 3D, independent corpora, other hosts and published artifacts remain separate evaluations. Earlier candidates’ regressions against older checkpoints are not resolved by this incremental V3000 comparison. No main merge, hosted CI, release or publication was performed.','',
f'[Evidence packet]({pre}-artifacts/README.md): every original/adaptive block, load observation, mismatch, source and binary/library/input hash, exact output streams, logs, measured source executables, reconstructable patch and complete frozen build/input archive. Archive bytes and both reconstructed source trees are verified. The [existing corpus provenance](../validation/benchmark_corpora/chemistry-speed-followups-v1.manifest.json) and [third-party notices](../THIRD_PARTY_NOTICES.md) apply.','']
report=d/'benchmarks'/(pre+'.md');report.write_text('\n'.join(lines))
(a/'README.md').write_text('''# Direct V3000 payload evidence

Read `source-manifest.json`, the three `paired21-*.json` original reports, `continued21.json`, `adaptive-selection.json` and any `focused21-*.json` follow-ups. All observations and semantic mismatches remain in their original lane. Strict wins require agreement, CI lower > 1 and 21/21 faster blocks. Source-only continuation/follow-up diagnostics never rank Open Babel.

`candidate-build-source.tar.gz` contains the complete frozen build/input snapshot under `candidate/`, exact fixed inputs under `corpora/` and the two previous source files under `baseline-overrides/`. Reconstruct the baseline by copying the candidate, applying these overrides and removing both added files listed in `source-manifest.json`. `verification.json` verifies every candidate archive hash and every reconstructed baseline hash. The supplied patch was separately applied and compared to the candidate.

`baseline-v3000` and `candidate-v3000` are the exact measured native source executables (macOS arm64). Open Babel's separately pinned native harness, libraries and plugins remain at manifest paths. The archived `tools/openbabel_v3000_throughput.cpp` reproduces the harness against that Open Babel installation. Rebuild source executables with `cargo build --release --locked -p chematic-mol --example openbabel_v3000_throughput`, using the recorded Rust toolchain and release/LTO settings.

Run the archived `scripts/bench_openbabel_v3000_throughput.py` in the pinned RDKit observer environment with `--baseline`, `--candidate`, `--openbabel`, `--corpus`, `--records 1000`, `--blocks 21`, `--repeats 5`, `--source-manifest` and `--output`; use `--records 64` for controlled chains. Darwin `/usr/bin/time -l` requires OS-statistics access. The continuation/source-only driver records its exact original private working paths; adapt those paths if reproducing elsewhere. All three original corpora are unchanged. Exact output streams are gzip-compressed, and `artifact-sha256.json` indexes the packet.

The initial full-test build failed with ENOSPC; its log and the complete successful retry are retained. Builds and validation finished before timing. These are version 1.1.4 source candidates. Hashes prove measured artifact identity and source reconstruction; they do not promise bit-identical rebuilt executable files or toolkit-wide speed superiority.
''')
spec=importlib.util.spec_from_file_location('names',s/'scripts/check_product_names.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);assert not mod.scan_paths(d,[report,d/'CHANGELOG.md',a/'README.md'])
index={p.relative_to(a).as_posix():sha(p) for p in sorted(a.rglob('*')) if p.is_file()};(a/'artifact-sha256.json').write_text(json.dumps(index,indent=2)+'\n');(r/'delivery-sha256.json').write_text(json.dumps({p.relative_to(d).as_posix():sha(p) for p in sorted(d.rglob('*')) if p.is_file()},indent=2)+'\n')
print('Packaged',len(index),'artifacts; native tests',counts)
