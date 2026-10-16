from pathlib import Path
import gzip, hashlib, json, shutil, statistics
ROOT=Path('/private/tmp/chematic-sdf-buffer-scan-20261016')
OUT=Path(__file__).resolve().parent
STEM='2026-10-16-python-model-hook-inline'
PACKET=ROOT/'benchmarks'/(STEM+'-artifacts')
plan=json.loads((OUT/'measurement-plan.json').read_text())
parity=json.loads((OUT/'parity.json').read_text())
native=json.loads((OUT/'native-validation.json').read_text())
matrix={str(offset):json.loads((OUT/f'frontier-{offset}.json').read_text()) for offset in plan['offsets']}
assert all(d['complete'] and len(d['blocks'])==plan['blocks'] for d in matrix.values())
assert parity['outputs']['baseline']==parity['outputs']['candidate']
primary={offset:{key:data['summary'][key]['baseline']['speedup'] for key in ['tpsa/first_use','labute_asa/first_use']} for offset,data in matrix.items()}
guard_failures={offset:{key:values['baseline']['speedup'] for key,values in data['summary'].items() if values['baseline']['speedup'] is None or values['baseline']['speedup']['median']<0.99} for offset,data in matrix.items()}
negative={offset:{key:values['baseline']['speedup'] for key,values in data['summary'].items() if values['baseline']['speedup'] and values['baseline']['speedup']['ci95'][1]<1} for offset,data in matrix.items()}
passed=all(v['ci95'][0]>1 for rows in primary.values() for v in rows.values()) and not any(guard_failures.values())
counts={offset:{arm:{'strict_wins':sum(row[arm]['demonstrated_win'] for row in data['summary'].values()),'eligible':sum(row[arm]['speedup'] is not None for row in data['summary'].values())} for arm in ['rdkit','cosmolkit']} for offset,data in matrix.items()}
decision={'adoption_passed':passed,'primary':primary,'guard_failures':guard_failures,'all_negative_source_intervals':negative,'comparator_strict_wins':counts,'universal_superiority':False,'epochs':'Two complete single 20-block epochs; no timing reruns, pooling or substitutions'}
(OUT/'decision.json').write_text(json.dumps(decision,indent=2)+'\n')
def fmt(value):
 if value is None:return 'unranked: output difference'
 return f"{value['median']:.3f} [{value['ci95'][0]:.3f}, {value['ci95'][1]:.3f}]; {value['winning_blocks']}/{value['blocks']}"
text=[f'# Inline aromatic model-gate scratch — 2026-10-16\n',f'Local adoption gate: **{"PASS" if passed else "FAIL"}**. RDKit/COSMolKit-wide superiority remains unproved. Every measured lane is retained.\n']
if not passed:text+=['This candidate must remain unpromoted: the fixed primary or retention gate below did not pass. CI correctness alone cannot override this performance decision.\n']
text += ['## Source change and preservation\n',f"Baseline: `{plan['base']}` (source version 1.1.5). In `rdkit_model_hook.rs`, common aromaticity-model validation uses inline scratch storage for the visited state, electron counts, component traversal, local adjacency and BFS. The FIFO queue uses a contiguous buffer and a read cursor. Atom/edge visitation order, shortest-path parents, electron counting, rejection rules, sanitization, caches and public APIs are unchanged. Larger buffers spill to the heap. No dependency, unsafe code or new chemical shortcut.\n",
'The inline buffers trade bounded stack storage for fewer heap allocations. Whole-process RSS is recorded below, rather than presented as per-operation or large-input memory proof.\n',
'The test oracle freezes three complete original gate functions from the named base and checks their saved tokens after restoring renamed identifiers. All 14,072 pinned inputs preserve cold/repeated gate decisions. Directed empty, dummy, radical, charged, dative, disconnected, aromatic and macrocycle cases are covered, along with fused graphs around the inline/spill boundaries, alternative electron assignments and non-component spectator atoms.\n',
'The two measured source wheels preserve all descriptor float bits/integers, Morgan bytes and MACCS bytes on the same 14,072 inputs: 112,576 descriptor outcomes, 56,288 Morgan outputs and 28,144 MACCS outcomes. Six pre-existing nonfinite descriptor outcomes are preserved, not scored as accuracy passes. These are exposed development inputs, not an independent holdout.\n',
'## Fixed protocol and gates\n',
'Two unfiltered 1,000-row slices (offsets 0 and 5,000) of the pinned 10,000-input corpus. Twenty rotating four-arm fresh-process blocks per slice, seven operations and three preparation modes: 42 source-relative lanes. No adaptive timing reruns or pooled epochs. Every worker retains native return objects during timing and normalizes values afterwards. First use parses fresh objects outside the timer and includes lazy work; repeated use warms each object once; parse inclusive runs actual parsing and the operation together. Competitor parsing may eagerly perceive chemistry. The timed TPSA lane uses RDKit-default N/O contributions.\n',
'Both source wheels were freshly built with Rust 1.97.0, opt3/full LTO/one codegen unit/one job, native-inchi and strip=none. Comparators: RDKit 2026.09.1, COSMolKit 0.5.0 and CPython 3.13.6 on one shared macOS arm64 host. Package-local native files are hashed before and after each profile; system/transitive libraries and native data files outside those package directories are excluded, so this is not a full runtime closure. Startup/import time is excluded. External host activity is uncontrolled; load averages and whole-process peak RSS are retained.\n',
'Ratios are reference time / candidate time. Brackets are paired-bootstrap 95% intervals using the unchanged deterministic 10,000-resample procedure, rounded to three decimals. Local adoption requires exact source outputs, positive interval lower bounds (>1) for all four first-use TPSA/ASA primary lanes, and all 42 source-relative medians >=0.99. Strict competitor wins additionally require output equivalence and all 20 blocks faster. Mismatches remain unranked. Negative intervals and budget failures are disclosed below.\n',
'## Primary decision\n','| Offset | Operation | Source ratio |\n|---|---|---:|']
for offset,rows in primary.items():
 for key,value in rows.items():text.append(f'| {offset} | {key} | {fmt(value)} |')
for offset,data in matrix.items():
 text += [f'\n## Complete matrix: offset {offset}\n','| Operation / mode | Baseline | RDKit | COSMolKit |\n|---|---:|---:|---:|']
 for key,values in data['summary'].items():
  cells=[fmt(values[arm]['speedup'])+(' **strict win**' if arm!='baseline' and values[arm]['demonstrated_win'] else '') for arm in ['baseline','rdkit','cosmolkit']]
  text.append('| '+key+' | '+' | '.join(cells)+' |')
text += ['\n## Every negative source interval\n','```json\n'+json.dumps(negative,indent=2)+'\n```\n','## Every median budget failure\n','```json\n'+json.dumps(guard_failures,indent=2)+'\n```\n','## Absolute batch times and whole-process memory\n','Milliseconds below are medians for a complete 1,000-input batch.\n','| Offset / operation / mode | Baseline ms | Candidate ms | RDKit ms | COSMolKit ms |\n|---|---:|---:|---:|---:|']
for offset,data in matrix.items():
 for key in data['summary']:
  op,mode=key.split('/')
  times=[statistics.median(b[arm]['operations'][op][mode] for b in data['blocks'])/1e6 for arm in ['baseline','candidate','rdkit','cosmolkit']]
  text.append(f'| {offset} / {key} | '+' | '.join(f'{v:.3f}' for v in times)+' |')
text += ['\nPeak RSS includes imports, preparation, warmup and every operation in the worker. It is not isolated operation memory.\n','| Offset | Baseline MiB | Candidate MiB | RDKit MiB | COSMolKit MiB |\n|---|---:|---:|---:|---:|']
for offset,data in matrix.items():
 values=[statistics.median(b[arm]['peak_process_rss_bytes'] for b in data['blocks'])/2**20 for arm in ['baseline','candidate','rdkit','cosmolkit']]
 text.append(f'| {offset} | '+' | '.join(f'{v:.1f}' for v in values)+' |')
text += ['\n## Validation, retention and limits\n',f"Rust perception/chemistry: {native['passed']:,} passed, {native['failed']} failed, {native['ignored']} ignored across {native['suites']} suites. Release all-target perception Clippy, formatting and the refactor budget passed. Python: 2,302 passed, 2 skipped. Native correctness uses opt3/LTO=false/16 codegen units in a separate target; these are not the timed full-LTO wheels.\n",
'Pre-timing failures are retained: the initial baseline build snapshot lacked a workspace tools directory; the initial SmallVec array sizes were unsupported by the existing dependency features; and the new empty-input test incorrectly asked the parser to accept an empty SMILES. These were corrected before formal timing. No failed build/test supplied observations.\n',
f'[{STEM}-artifacts]({STEM}-artifacts/README.md) retains complete raw matrices, all slice observer outputs and mismatch rows, decisions, the frozen plan, source/library/wheel manifests, exact reconstructable candidate patch, launch/validation drivers and build/test logs. Full wheels, installed extensions, original source snapshots and all 14,072-input output streams remain under `{OUT}` with versioned hashes. Reconstructable source does not imply bit-reproducible rebuilt binaries.\n',
'Additional hosts, sealed holdout inputs, other operations, 3D/reaction workflows, cold CLI startup, Python/WASM publication and toolkit-wide superiority remain open. Source-relative interval movement in an untouched API does not establish a causal code regression. A source result also does not establish a published-package win.\n']
(ROOT/'benchmarks'/(STEM+'.md')).write_text('\n'.join(text)+'\n')
PACKET.mkdir(exist_ok=True)
for p in OUT.iterdir():
 if p.is_file() and (p.suffix in ('.json','.py','.log','.patch') or p.name.startswith('observer-')):
  if p.name.startswith('frontier-') and p.suffix=='.json':
   (PACKET/(p.name+'.gz')).write_bytes(gzip.compress(p.read_bytes(),mtime=0))
  else:shutil.copy2(p,PACKET/p.name)
(PACKET/'README.md').write_text(f'# Inline model-gate evidence packet\n\nRead [the report](../{STEM}.md) with the predeclared measurement plan and `decision.json`. Complete single-epoch matrices are compressed without dropping any raw lane. `observer-*` files contain fresh post-timing values; their elapsed fields were discarded. Source output dumps and wheel hashes are in `parity.json`; full binaries and 14,072-input streams remain local at `{OUT}`.\n\n`candidate.patch` reconstructs the two changed/added Rust files from base `{plan["base"]}`. Compiler flags and build logs, the complete source/input manifest, package-local native library hashes, frozen-function extraction proof and validation commands are included. The cleanup record lists only unused regenerable compiler metadata/static libraries removed after process/open-file checks. No timing worker was restarted.\n')
print(json.dumps(decision,indent=2))
