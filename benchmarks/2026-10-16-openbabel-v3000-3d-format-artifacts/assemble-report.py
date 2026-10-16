from pathlib import Path
import gzip,hashlib,json,shutil,statistics
ROOT=Path('/private/tmp/chematic-sdf-buffer-scan-20261016');OUT=Path(__file__).resolve().parent
STEM='2026-10-16-openbabel-v3000-3d-format';PACKET=ROOT/'benchmarks'/(STEM+'-artifacts')
plan=json.loads((OUT/'measurement-plan.json').read_text());validation=json.loads((OUT/'validation.json').read_text())
matrix={v['kind']+'-'+v['profile']:json.loads((OUT/f"paired-{v['kind']}-{v['profile']}.json").read_text()) for v in plan['profiles']}
assert all(d['complete'] and len(d['blocks'])==21 for d in matrix.values())
negative={};below={};unranked=[]
for p,d in matrix.items():
 negative[p]={};below[p]={}
 for mode,row in d['summary'].items():
  v=row['baseline']['speedup']
  if v is None:unranked.append([p,mode]);continue
  if v['ci95'][1]<1:negative[p][mode]=v
  if v['median']<.99:below[p][mode]=v
primary={p:matrix[p]['summary']['write']['baseline'] for p in ['v3000-3d-spiral','v3000-3d-offset']}
primary_pass=all(v['speedup'] and v['speedup']['ci95'][0]>1 for v in primary.values())
equal=all(d['source_parse_snapshots_equal'] and d['rewrite_sha256']['baseline']==d['rewrite_sha256']['candidate'] for d in matrix.values())
passed=validation['passed'] and equal and primary_pass and not any(below.values()) and not unranked
decision={'local_adoption_gate_passed':passed,'all_source_snapshots_and_rewrites_exact':equal,'primary':primary,'primary_ci_pass':primary_pass,'all_24_source_median_floor':.99,'below_floor':below,'unranked_source_lanes':unranked,'all_negative_source_intervals':negative,'openbabel_strict_wins':{p:{mode:v['openbabel']['strict_speed_win'] for mode,v in d['summary'].items()} for p,d in matrix.items()},'eight_complete_single_epochs_no_sample_pooling':True,'initial_failed_rss_launch_completed_blocks':0,'universal_superiority':False}
(OUT/'decision.json').write_text(json.dumps(decision,indent=2)+'\n')
native=validation['native']
text=['# V3000 3D coordinate formatting — 2026-10-16/17 (JST)\n',f'Local adoption gate: **{"PASS" if passed else "FAIL"}**. Open Babel-wide superiority remains unproved. These are native source results on one shared macOS arm64 host.\n']
if not passed:text+=['Keep this candidate unpromoted: the frozen primary or retention rule below failed. Correctness CI alone does not change that performance decision.\n']
text+=['## Implementation and preservation\n',f"Base `{plan['base_revision']}`, source version 1.1.6. Only `append_atoms_3d` changes in production: it appends atom/map integers and x/y/z using the existing integer and `push_fixed4` helpers already used by V2000 and 2D V3000. Ties, large magnitudes and nonfinite values still use the exact Rust formatter fallback. Atom properties, missing-coordinate zeros, checked-writer validation, bond writing, stereo, collections, headers, allocation policy, parsing and public APIs are unchanged. No unsafe code or new dependency. Earlier unpromoted candidates are excluded.\n",
'The complete original 3D atom function is frozen in a private regression test and independently checked against base tokens after its name is restored. Tests compare whole atom sections, retained prefixes, signed zeros, numeric/map/id boundaries, missing/extra points, pseudoatoms and Unicode/duplicate opaque properties. An additional 20,000 deterministic f64 bit patterns and 60,000 rounding-tie/adjacent-value cases match exactly. Existing full-writer 2D/3D/checked-writer differential, rich metadata and stereo tests remain active. All 4,256 profile-input snapshots and rewrite byte streams match the baseline; chemical structures overlap across profiles.\n',
f"MOL validation: {native['passed']:,} passed, {native['failed']} failed, {native['ignored']} ignored across {native['suites']} suites (`--lib --tests`). All-target release Clippy, formatting, refactor budget and ten runner semantic gates passed. Correctness tests use opt3/LTO=false/16 codegen units; timing binaries use full LTO/one codegen unit. Source/build snapshots contain Git-tracked inputs plus the named 3D probe, excluding ignored Python caches.\n",
'## Frozen workload and limits\n',
'The six existing retention profiles are unchanged: two unfiltered 1,000-record SDF slices, two unfiltered 1,000-record V3000 slices (offsets 0 and 5,000), 64 linear alkanes with 16–1,024 atoms, and continuation-only copies of those chains. The two added primary profiles each reuse the 64 achiral chains with explicit deterministic nonplanar finite coordinates. `controlled-3d-recipe.json` records source/output hashes and formulas. These are synthetic coordinate-serialization controls, not physically optimized conformers, diverse 3D chemistry, embedding quality or an independent holdout. No rows were filtered after measuring outputs.\n',
'The new 3D probe calls the real public conformer writer and refuses an absent parsed conformer before timing during warmup. The same probe is built against both source trees. Raw source differs only through Rustfmt normalization (closure blocks/trailing commas); `frozen-writer-proof.json` verifies equality after applying that formatter and retains both raw hashes. Two-dimensional probes and the Open Babel harnesses are unchanged.\n',
'Preparation and the plan were frozen on 2026-10-16 JST; formal collection crossed midnight into 2026-10-17 without changing the workload or binaries.\n',
'## Protocol and adoption rule\n',
'Twenty-one rotating three-arm fresh-process blocks per profile; three native batch repetitions for parse, prepared write and parse-plus-write. Startup/input loading/plugin discovery are excluded; V3000 framing is outside the timer, SDF framing inside. Parse retains the whole diagnostic batch. Prepared write follows one warm write. SDF reuses a Rust String; V3000 and Open Babel return allocating strings. The 3D probe retains all xyz fields. RSS is whole-process high-water memory including preparation and warmup, not isolated operation or streaming memory. External host activity is uncontrolled; load averages are retained.\n',
'Rust 1.97.0, opt3/full LTO/one codegen unit/one job, incremental=false, strip=none. Open Babel 3.2.1 C++ -O3 executables are reused by exact hashes; historical timing samples are not reused. RDKit 2026.09.1 independently observes graph, indexed atoms/bonds, xyz within 1e-4, names, properties and explicit-unknown stereo. Internal canonical graphs gate parse rankings; accepted equivalent rewrites gate write/roundtrip. Rejected/mismatching rows stay in place and unranked. The 281 package-local native-file entries are rehashed before/after profiles; aliases may repeat contents. System/transitive libraries and package-external data are outside this seal.\n',
'Ratios are reference time / candidate time; intervals use the unchanged deterministic 10,000-resample paired bootstrap, rounded to three decimals. Adoption was fixed before timing: exact source bytes/full errors, both controlled 3D write CI lower>1, every one of 24 source-relative medians >=0.99, and disclosure of all negative intervals. A strict Open Babel win additionally needs matching declared outputs and 21/21 faster blocks. Single-host source results do not establish published-package or toolkit-wide superiority.\n',
'An initial startup attempt reached one baseline SDF parse child but `/usr/bin/time -l` failed to obtain `sysctl kern.clockrate` inside the sandbox; no block or valid timing sample was accepted. Its incomplete report and failure logs are retained separately. After obtaining process-statistics permission, exactly one complete 21-block epoch per profile ran with the unchanged frozen sources, binaries, corpus and plan. No performance-informed rerun, adaptive filtering, sample pooling or substitution occurred. A prior integrity-only launch also stopped before workers because Rustfmt completed after a test-file hash was sealed; that test-only manifest was re-sealed before measurement.\n',
'## Complete matrix\n','Median milliseconds per full batch (native elapsed divided by three repetitions). Ratios are paired, not ratios of displayed medians.\n','| Profile | Mode | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; wins | Open Babel / candidate [95% CI]; wins |\n|---|---|---:|---:|---|---|']
def cell(row):
 v=row['speedup']
 return 'unranked: semantic prerequisite fails' if v is None else f"{v['median']:.3f} [{v['ci95'][0]:.3f}, {v['ci95'][1]:.3f}]; {v['winning_blocks']}/21"+(' **strict win**' if row['strict_speed_win'] else '')
for p,d in matrix.items():
 for mode,row in d['summary'].items():
  ms=[statistics.median(b[a][mode]['elapsed_ns']/b[a][mode]['repeats']/1e6 for b in d['blocks']) for a in ['baseline','candidate']]
  text.append(f'| {p} | {mode} | {ms[0]:.3f} | {ms[1]:.3f} | {cell(row["baseline"])} | {cell(row["openbabel"])} |')
text+=['\n## Decision, every negative interval and every failed floor\n','```json\n'+json.dumps(decision,indent=2)+'\n```','\n## Whole-process RSS\n','| Profile | Mode | Baseline MiB | Candidate MiB | Open Babel MiB |\n|---|---|---:|---:|---:|']
for p,d in matrix.items():
 for mode in d['summary']:
  values=[statistics.median(b[a][mode]['peak_process_rss_bytes'] for b in d['blocks'])/2**20 for a in ['baseline','candidate','openbabel']]
  text.append('| '+p+' | '+mode+' | '+' | '.join(f'{v:.2f}' for v in values)+' |')
text+=['\n## Retained proof\n',f'The [{STEM}-artifacts]({STEM}-artifacts/README.md) packet preserves every raw matrix, decision, mismatch and negative interval, frozen plan/recipe, source/library/corpus/binary hashes, reconstructable patch and complete logs/drivers. Full binaries, generated profile inputs, parse snapshots and rewrites remain locally at `{OUT}` with hashes. The audit reopens all saved outputs, checks fresh post-timing native dump/snapshot/SMILES observers against the same declared gates, and recomputes all 48 paired comparisons; source reconstruction does not imply bit-reproducible rebuilds.\n',
'Pre-timing test setup corrections are retained: the empty graph constructor was corrected and an unnecessary trait import removed. Tiny fixtures absent from sparse checkout were restored unchanged from Git. Only completed dedicated build targets were removed after process/open-file checks; the eight measured executables, all source and evidence were preserved and rehashed. The checkpoint wrapper retries only the same in-memory write after ENOSPC, never a worker or timer.\n',
'Open items include diverse/realistic 3D conformers, other formats, larger/streaming memory, additional hosts, sealed holdout chemistry, cold CLI startup, Python/WASM and published artifacts. Unchanged lanes moving in the matrix do not isolate the cause, but still count under the frozen retention rule.\n']
(ROOT/'benchmarks'/(STEM+'.md')).write_text('\n'.join(text)+'\n');PACKET.mkdir(exist_ok=True)
for p in OUT.iterdir():
 if not p.is_file():continue
 if p.suffix in ('.log','.patch') or (p.name.startswith('paired-') and p.suffix=='.json'):
  (PACKET/(p.name+'.gz')).write_bytes(gzip.compress(p.read_bytes(),mtime=0))
 elif p.suffix in ('.json','.py'):shutil.copy2(p,PACKET/p.name)
(PACKET/'README.md').write_text(f'# 3D V3000 formatting evidence\n\nRead [the report](../{STEM}.md), frozen plan and `decision.json`. All eight complete 21-block matrices and every mismatch/negative interval are retained. Failed startup reports contain zero accepted blocks and are separate from the complete epochs. Source/binary/corpus/library hashes, full validation logs, reconstructable patch and audit are included. Full native executables and saved output streams remain locally at `{OUT}`. The 3D workload is synthetic achiral chains, not physical conformer or toolkit-wide proof.\n')
print(json.dumps(decision,indent=2))
