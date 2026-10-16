from pathlib import Path
import gzip, hashlib, json, shutil, statistics
ROOT=Path('/private/tmp/chematic-sdf-buffer-scan-20261016')
OUT=Path(__file__).resolve().parent
STEM='2026-10-16-openbabel-v3000-borrowed-rows'
PACKET=ROOT/'benchmarks'/(STEM+'-artifacts')
plan=json.loads((OUT/'measurement-plan.json').read_text())
validation=json.loads((OUT/'validation.json').read_text())
matrices={v['kind']+'-'+v['profile']:json.loads((OUT/f"paired-{v['kind']}-{v['profile']}.json").read_text()) for v in plan['profiles']}
assert all(d['complete'] and len(d['blocks'])==21 for d in matrices.values())
negative={};below={};unranked=[]
for profile,data in matrices.items():
 negative[profile]={};below[profile]={}
 for mode,row in data['summary'].items():
  value=row['baseline']['speedup']
  if value is None:unranked.append([profile,mode]);continue
  if value['ci95'][1]<1:negative[profile][mode]=value
  if value['median']<.99:below[profile][mode]=value
primary={p:matrices[p]['summary']['parse']['baseline'] for p in ['v3000-0','v3000-5000']}
primary_pass=all(v['speedup'] and v['speedup']['ci95'][0]>1 for v in primary.values())
equal=all(d['source_parse_snapshots_equal'] and d['rewrite_sha256']['baseline']==d['rewrite_sha256']['candidate'] for d in matrices.values())
passed=validation['passed'] and equal and primary_pass and not any(below.values()) and not unranked
decision={'local_adoption_gate_passed':passed,'all_source_snapshots_and_rewrites_exact':equal,'primary':primary,'primary_ci_pass':primary_pass,'all_18_source_median_floor':.99,'below_floor':below,'unranked_source_lanes':unranked,'all_negative_source_intervals':negative,'openbabel_strict_wins':{p:{mode:r['openbabel']['strict_speed_win'] for mode,r in d['summary'].items()} for p,d in matrices.items()},'no_timing_samples_reused':True,'universal_superiority':False}
(OUT/'decision.json').write_text(json.dumps(decision,indent=2)+'\n')
native=validation['native']
text=[f'# Borrowed V3000 logical rows — 2026-10-16\n',f'Local adoption gate: **{"PASS" if passed else "FAIL"}**. Complete Open Babel-wide superiority remains unproved. Native source measurements on one shared macOS arm64 host are separate from published artifacts.\n']
if not passed:text+=['Both primary parse intervals improve (1.020 [1.018, 1.023] and 1.041 [1.032, 1.049]), but the fixed retention gate fails: SDF writing has medians 0.986 and 0.983, and continued V3000 parsing has median 0.989, each below the 0.99 floor. The SDF writer algorithm is unchanged; these measurements do not isolate the cause. This candidate must remain unpromoted even if correctness CI succeeds.\n']
text += ['## Implementation and preservation\n',f"Baseline: clean source main `{plan['base_revision']}`, source version 1.1.5. Both Rust arms are freshly built from this base and the one-file production delta; earlier unpromoted sparse-ID and descriptor trials are excluded. The baseline includes requested-only writer losses and the dative-carbonyl E/Z correction.\n",
'V3000 reading previously built a complete vector of logical-line wrappers even when all payloads were already borrowed from the raw input. The candidate detects that no exact V30 payload ends in a continuation hyphen, counts the same logical rows, and feeds their borrowed payloads directly to the unchanged parser. Any continuation takes the original collection/joining path. Physical-line validation, error precedence, Unicode tokenization, exact allocation hints, chemical interpretation, full diagnostics, public signatures and dependencies are unchanged. No unsafe code. Physical lines and the full input remain resident; this is not a streaming-file API or memory bound.\n',
'A frozen complete original reader from the named base compares full returned graphs, atom/bond/tag/R-group/stereo-neighbor and group-membership data, coordinate and conformer float bits, metadata, every diagnostic family, and rewrites. It checks all 2,064 pinned V3000 records twice, long/short/empty/Unicode/continued rows, exact malformed-input error variants and fields, unusual headers, huge counts, allocation capacities, 3D stereo and enhanced collections. The original reader function tokens are independently checked after restoring only its name/visibility. Existing rich/query/opaque/3D tests remain active.\n',
f"Full MOL checks: {native['passed']:,} passed, {native['failed']} failed, {native['ignored']} ignored across {native['suites']} suites. Ten runner semantic gates, release all-target MOL Clippy, formatting and refactor budget pass. Native correctness uses opt3/LTO=false/16 codegen units, separate from the full-LTO/one-codegen-unit timing binaries. Complete preflight snapshots and rewrite bytes match across all 4,128 profile inputs. These inputs overlap in chemical structures.\n",
'## Frozen protocol\n',
'Two unfiltered 1,000-record SDF slices and two unfiltered 1,000-record V3000 slices at offsets 0 and 5,000, plus 64 controlled chains of 16–1,024 atoms and their continuation-only variants. Every atom row in the latter is split after its x coordinate; all seven declared observer fields match the original 64 chains. The two chain profiles represent the same chemistry with different physical formatting. All inputs are exposed development data, not a sealed holdout. Existing corpus provenance/licenses and third-party notices apply.\n',
'Twenty-one rotating three-arm fresh-process blocks per profile and three native batch repetitions per lane. Startup, input loading and initial plugin discovery are outside timers. SDF framing is inside parsing; V3000 delimiter splitting is outside. Parse retains the full diagnostic batch. Prepared write follows an untimed warm write; roundtrip parses fresh inputs and writes inside the timer. Rust SDF reuses its String; V3000 and C++ use allocating native return strings. Whole-process RSS includes loading, preparation, warmup, allocator context and plugins; it is not operation-only or streaming memory. External host activity is uncontrolled and load averages are retained.\n',
'Rust 1.97.0, release opt3/full LTO/one codegen unit/one build job, incremental=false, strip=none. Open Babel 3.2.1 C++ -O3 harnesses are reused by exact binary/library hashes; no historical timings are reused. RDKit 2026.09.1 is an independent observer, not a timed competitor. The 281 package-local native-file entries are pinned before/after each profile; aliases can repeat contents. System/transitive libraries outside the named package directories and native data files are outside this seal.\n',
'Ratios are reference time / candidate time. Paired median intervals use the unchanged deterministic 10,000-resample bootstrap and three-decimal reporting. Strict Open Babel wins require semantic prerequisites, interval lower>1 and 21/21 faster blocks. Internal canonical graph equality gates parsing; all seven declared rewrite fields gate write/roundtrip. Rejected/mismatching records remain in their original positions and unranked without filtering or repair.\n',
'The local plan was fixed before timing: exact source snapshots, rewrites and full reports/errors; V3000 parse interval lower>1 on both diverse slices; all 18 source-relative medians >=0.99; disclose every negative interval. No adaptive timing reruns, pooled epochs, substituted samples or relaxed gates. An unchanged lane moving does not establish its cause.\n',
'## Complete matrix\n','Times are median milliseconds per complete batch, native elapsed divided by three repetitions. Ratios use paired observations.\n',
'| Profile | Mode | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; wins | Open Babel / candidate [95% CI]; wins |\n|---|---|---:|---:|---|---|']
def cell(value):
 r=value['speedup']
 return 'unranked: semantic prerequisite fails' if r is None else f"{r['median']:.3f} [{r['ci95'][0]:.3f}, {r['ci95'][1]:.3f}]; {r['winning_blocks']}/21"+(' **strict win**' if value['strict_speed_win'] else '')
for profile,data in matrices.items():
 for mode,row in data['summary'].items():
  times={arm:statistics.median(b[arm][mode]['elapsed_ns']/b[arm][mode]['repeats']/1e6 for b in data['blocks']) for arm in ['baseline','candidate']}
  text.append(f"| {profile} | {mode} | {times['baseline']:.3f} | {times['candidate']:.3f} | {cell(row['baseline'])} | {cell(row['openbabel'])} |")
text += ['\n## Full decision and every negative interval\n','```json\n'+json.dumps(decision,indent=2)+'\n```\n','## Whole-process peak RSS\n','| Profile | Mode | Baseline MiB | Candidate MiB | Open Babel MiB |\n|---|---|---:|---:|---:|']
for profile,data in matrices.items():
 for mode in data['summary']:
  memory=[statistics.median(b[arm][mode]['peak_process_rss_bytes'] for b in data['blocks'])/2**20 for arm in ['baseline','candidate','openbabel']]
  text.append(f'| {profile} | {mode} | '+' | '.join(f'{v:.2f}' for v in memory)+' |')
text += ['\n## Retention and limits\n',
'The immutable file manifest also includes pre-existing Python bytecode caches copied with the baseline tree. Those files are unused by both native Rust executables and the external RDKit observer environment. `manifest-entry-classification.json` separates them from Git-tracked code/build inputs. Cache hashes are incidental files, not reconstructable source. All native source/build inputs match committed files; no timing observation or adoption decision changed.\n',
'Existing tiny test fixtures excluded by sparse checkout were restored from unchanged Git contents before testing. The timing wrapper reserves 64MiB and can retry only the same in-memory owned output write after ENOSPC; it never repeats a worker or timer.\n',
f'The [{STEM}-artifacts]({STEM}-artifacts/README.md) packet retains complete raw matrices and every mismatch/negative interval, decisions, source/compiler/corpus/library/binary hashes, frozen plan, reconstructable patch, full-report oracle proof, build/test logs and launch/verification drivers. Full binaries, generated profile inputs and complete source snapshot/rewrite streams remain locally under `{OUT}` with versioned hashes. Source reconstruction does not imply bit-reproducible compiler output.\n',
'Additional hosts, sealed holdout chemistry, richer query/SGROUP throughput, cold CLI startup, Python/WASM throughput, published packages, other formats and FF/3D performance remain open. Correctness checks do not prove toolkit-wide superiority, and source measurements are separate from main integration, hosted CI and publication.\n']
(ROOT/'benchmarks'/(STEM+'.md')).write_text('\n'.join(text)+'\n')
PACKET.mkdir(exist_ok=True)
for p in sorted(OUT.iterdir()):
 if not p.is_file():continue
 if p.suffix in ('.log','.patch') or p.name.startswith('paired-') and p.suffix=='.json':
  (PACKET/(p.name+'.gz')).write_bytes(gzip.compress(p.read_bytes(),mtime=0))
 elif p.suffix in ('.json','.py'):shutil.copy2(p,PACKET/p.name)
(PACKET/'README.md').write_text('Read [the report](../'+STEM+'.md) with the fixed plan and decision. Complete raw matrices, all negative observations and mismatch rows are retained. Binaries/full generated streams stay local with versioned hashes; no historical timing is reused. The candidate patch reconstructs measured source files from the pinned base.\n')
with (PACKET/'README.md').open('a') as output:
 output.write('\nThe immutable measurement manifest includes incidental pre-existing Python bytecode caches. `manifest-entry-classification.json` separates those unused entries from Git-tracked native source/build inputs. Cache bytes are not a source-reconstruction claim. Timing samples and the failed adoption decision are unchanged.\n')
seal={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(PACKET.iterdir()) if p.is_file() and p.name!='packet.sha256.json'}
(PACKET/'packet.sha256.json').write_text(json.dumps(seal,indent=2)+'\n')
print(json.dumps(decision,indent=2))
