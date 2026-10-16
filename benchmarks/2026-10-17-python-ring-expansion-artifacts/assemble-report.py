from pathlib import Path
import gzip,hashlib,json,shutil,statistics
ROOT=Path('/private/tmp/chematic-python-ring-first-use-20261017')
OUT=Path(__file__).resolve().parent
STEM='2026-10-17-python-ring-expansion'
PACKET=ROOT/'benchmarks'/(STEM+'-artifacts')
plan=json.loads((OUT/'measurement-plan.json').read_text())
manifest=json.loads((OUT/'source-manifest.json').read_text())
parity=json.loads((OUT/'parity.json').read_text())
native=json.loads((OUT/'native-validation.json').read_text())
matrix={str(o):json.loads((OUT/f'frontier-{o}.json').read_text()) for o in plan['offsets']}
assert all(d['complete'] and len(d['blocks'])==20 for d in matrix.values())
assert parity['outputs']['baseline']==parity['outputs']['candidate']
primary=matrix['0']['summary']['num_rings/first_use']['baseline']['speedup']
guards={o:{k:v['baseline']['speedup'] for k,v in d['summary'].items() if not v['baseline']['speedup'] or v['baseline']['speedup']['median']<.99} for o,d in matrix.items()}
negative={o:{k:v['baseline']['speedup'] for k,v in d['summary'].items() if v['baseline']['speedup'] and v['baseline']['speedup']['ci95'][1]<1} for o,d in matrix.items()}
passed=primary is not None and primary['ci95'][0]>1 and not any(guards.values())
counts={o:{arm:{'strict_wins':sum(v[arm]['demonstrated_win'] for v in d['summary'].values()),'eligible':sum(v[arm]['speedup'] is not None for v in d['summary'].values())} for arm in ['rdkit','cosmolkit']} for o,d in matrix.items()}
decision={'adoption_passed':passed,'primary':{'offset':0,'lane':'num_rings/first_use','speedup':primary},'guard_failures':guards,'all_negative_source_intervals':negative,'comparator_strict_wins':counts,'universal_superiority':False,'epochs':'Exactly one complete 20-block epoch per profile; no timing reruns or pooling'}
(OUT/'decision.json').write_text(json.dumps(decision,indent=2)+'\n')
def fmt(x):
 return 'unranked: unequal outputs' if x is None else f"{x['median']:.3f} [{x['ci95'][0]:.3f}, {x['ci95'][1]:.3f}]; {x['winning_blocks']}/{x['blocks']}"
text=[f'# Relevant-cycle expansion — 2026-10-17\n',f'Local performance adoption gate: **{"PASS" if passed else "FAIL"}**. RDKit/COSMolKit-wide superiority is not established.\n']
if not passed:text+=['This candidate remains unpromoted. Correctness checks do not override the failed performance gate.\n']
text += ['## Change and correctness boundary\n',f"Measured baseline `{plan['base']}`; measured candidate `{manifest['candidate_commit']}`; source version 1.1.6. The candidate compares two shortest paths at equal distances from their shared root. A shared vertex has one root distance, so this changes the intersection check from a nested membership scan to a single aligned scan without changing the predicate.\n",
'`relevant_cycles_with_limits` separates the output-ring limit from the per-endpoint shortest-path limit. Existing relevant-cycle APIs continue to use the same limit for both. The RDKit ring descriptor bounds output expansion at its existing three-times-legacy rejection threshold, retaining the original one-million path limit and legacy fallback. Above that output threshold the final descriptor was already the legacy count. No ring rule, atom order, sorting, graph eligibility, sanitization or cache rule changes. The broader relevant-cycle algorithm and its legacy compatibility limitations remain.\n',
'The test oracle freezes the complete original expansion and intersection functions from the baseline, verified after restoring renamed identifiers. All 1,024 simple five-vertex graphs compare exact ring sequence, every output cap through count+1, equal-limit errors, and every shortest-path pair. Cages, disconnected graphs, long multiword macrocycles, and a path/output-cap separation regression are also checked. Existing fixed-10k ring-count, edit-invalidation and arbitrary-graph differential tests remain active.\n',
'All 14,072 exposed inputs preserve exact source descriptor float bits/integers, Morgan bytes, and MACCS bytes: 112,576 descriptor outcomes, 56,288 Morgan outputs, and 28,144 MACCS outcomes. Six pre-existing nonfinite descriptor outcomes are retained and excluded from accuracy claims. This is source preservation, not a claim that every outcome matches RDKit or COSMolKit.\n',
'## Frozen protocol\n',
'Two unfiltered 1,000-row slices of the pinned 10k corpus, at offsets 0 and 5,000. Each has one fixed epoch of 20 rotating four-arm fresh-process blocks: baseline wheel, candidate wheel, RDKit and COSMolKit. Seven operations and three modes give 42 source lanes. Native return values stay alive during timing; normalization is outside. First-use molecules are freshly parsed outside timing; repeated use warms each molecule once; parse-inclusive executes actual parsing and the operation together. RDKit parsing eagerly prepares ring information, while chematic can defer work; first-use timings alone do not represent equal end-to-end work. TPSA is RDKit-default N/O.\n',
'Fresh source wheels use Rust 1.97.0, opt3/full LTO/one codegen unit/one job, native-inchi, strip=none. Comparators are RDKit 2026.09.1 and COSMolKit 0.5.0 on CPython 3.13.6, one shared macOS arm64 host. Complete source/input hashes, wheel/extension hashes and package-local comparator native hashes are checked before/after each profile. These hashes exclude system libraries and package-external runtime data. Load averages and whole-process peak RSS are retained; external host activity is uncontrolled.\n',
'Ratios are reference time / candidate time. Intervals use the unchanged deterministic paired-bootstrap procedure with 10,000 resamples. The primary is offset-0 first-use ring count, with 95% interval lower bound >1. All 42 source-relative medians must also be >=0.99 and all 14,072 source outcomes exact. Every negative interval and failed floor is disclosed. Strict competitor wins additionally require equivalent outputs and 20/20 faster blocks; unequal outputs are unranked. No adaptive reruns, pooled epochs or sample substitutions.\n',
'The initial pre-timing plan covered output limiting alone. Code inspection then identified the quadratic path intersection; the candidate description was amended before the candidate wheel and before any formal timing. The primary, corpora, block count and acceptance gates did not change. Both plans are retained. An older source wheel was used only for preliminary per-row diagnosis.\n',
'## Primary decision\n',f"Offset 0, `num_rings/first_use`, source ratio: **{fmt(primary)}**.\n"]
for o,d in matrix.items():
 text += [f'## Complete matrix: offset {o}\n','| Operation / mode | Baseline | RDKit | COSMolKit |\n|---|---:|---:|---:|']
 for k,values in d['summary'].items():
  cells=[fmt(values[a]['speedup'])+(' **strict win**' if a!='baseline' and values[a]['demonstrated_win'] else '') for a in ['baseline','rdkit','cosmolkit']]
  text.append('| '+k+' | '+' | '.join(cells)+' |')
text += ['\n## Every negative source interval\n','```json\n'+json.dumps(negative,indent=2)+'\n```\n','## Every retention-floor failure\n','```json\n'+json.dumps(guards,indent=2)+'\n```\n','## Absolute batch times and process memory\n','Median milliseconds per complete 1,000-input batch.\n','| Offset / lane | Baseline | Candidate | RDKit | COSMolKit |\n|---|---:|---:|---:|---:|']
for o,d in matrix.items():
 for k in d['summary']:
  op,mode=k.split('/');ns=[statistics.median(b[a]['operations'][op][mode] for b in d['blocks'])/1e6 for a in ['baseline','candidate','rdkit','cosmolkit']]
  text.append(f'| {o} / {k} | '+' | '.join(f'{x:.3f}' for x in ns)+' |')
text += ['\nRSS includes imports, preparation, warmup and all operations; it is not isolated operation memory.\n','| Offset | Baseline MiB | Candidate MiB | RDKit MiB | COSMolKit MiB |\n|---|---:|---:|---:|---:|']
for o,d in matrix.items():
 rss=[statistics.median(b[a]['peak_process_rss_bytes'] for b in d['blocks'])/2**20 for a in ['baseline','candidate','rdkit','cosmolkit']]
 text.append(f'| {o} | '+' | '.join(f'{x:.1f}' for x in rss)+' |')
text += ['\n## Validation and retained evidence\n',f"Final native perception/chemistry tests: {native['passed']:,} passed, {native['failed']} failed, {native['ignored']} ignored in {native['suites']} suites. Both affected crates passed release all-target Clippy. Format and refactor budget passed. Python: {native['python_result']}. Native tests use opt3/LTO=false/16 codegen units; these binaries supply no timing samples.\n",
'Pre-timing setup failures are retained: the initial sparse build-input snapshot lacked data files and was completed before a successful wheel build; native/Python validation needed restoration of tracked nested fixtures. Both wheels first hit 30 sandbox denials while creating non-UTF8 filename fixtures; rerunning with appropriate permissions passed all supported cases and explicitly skipped those 30 filesystem-unsupported cases. Initial denied logs are retained. The initial focused test pass covered output limiting; final native validation was rerun after the linear intersection and expanded differential tests. No failed or superseded check supplied observations.\n',
f'[{STEM}-artifacts]({STEM}-artifacts/README.md) retains the frozen plans, complete raw timing matrices, every post-timing observer output and mismatch row, exact source dumps, source/wheel/native manifests, reconstructable candidate patch, build/validation/launch drivers and logs. Local source snapshots and the exact wheel binaries remain at `{OUT}`. Source reconstruction does not guarantee bit-reproducible binaries.\n',
'Additional hosts, sealed holdouts, other operations, 3D/reaction workloads and published Python/WASM packages remain unmeasured. Source intervals in untouched APIs do not by themselves establish a causal regression. A finite source comparison cannot prove universal or package-wide superiority.\n']
(ROOT/'benchmarks').mkdir(exist_ok=True)
(ROOT/'benchmarks'/(STEM+'.md')).write_text('\n'.join(text)+'\n')
PACKET.mkdir(exist_ok=True)
for p in OUT.iterdir():
 if p.is_file() and (p.suffix in ['.json','.py','.log','.patch','.gz']):
  if p.name.startswith('frontier-') and p.suffix=='.json':(PACKET/(p.name+'.gz')).write_bytes(gzip.compress(p.read_bytes(),mtime=0))
  else:shutil.copy2(p,PACKET/p.name)
(PACKET/'README.md').write_text(f'# Relevant-cycle expansion evidence\n\nRead [the complete report](../{STEM}.md) and `decision.json`. Each `frontier-*.json.gz` contains every raw observation from one 20-block epoch. `observer-*` files are fresh post-timing outputs; their elapsed fields were discarded. Full 14,072-input source streams are compressed without filtering. `candidate.patch` reconstructs the three modified Rust files and one added oracle test file from `{plan["base"]}`. Source snapshots, wheels and installed extensions remain at `{OUT}`. No timing reruns or pooling.\n')
print(json.dumps(decision,indent=2))
