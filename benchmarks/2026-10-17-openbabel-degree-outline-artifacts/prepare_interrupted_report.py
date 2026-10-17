from pathlib import Path
import hashlib,json,statistics,sys
O=Path(__file__).resolve().parent;R=Path('/private/tmp/chematic-degree-only-steady-20261017')
sys.path.insert(0,str(R/'scripts'));from bench_openbabel_sdf_throughput import paired_median_ratio_ci
plan=json.loads((O/'measurement-plan.json').read_text());continuation=json.loads((O/'continuation-seal.json').read_text());finished=json.loads((O/'continuation-finished.json').read_text());assert finished['all_preexisting_raw_files_unchanged'] and finished['interrupted_profile_not_rerun']
for name,h in continuation['existing_raw_files_sha256'].items():assert hashlib.sha256((O/name).read_bytes()).hexdigest()==h
source={};primary={};external={};floors=[];negative=[];missing=[];complete=[];partial=[];persisted=0
for p in plan['profiles']:
 key=p['kind']+'-'+p['profile'];data=json.loads((O/('paired-'+key+'.json')).read_text());count=len(data['blocks']);arms=p['orders'][0];persisted+=count*len(arms)*3
 if not data['complete']:
  assert key=='v3000-continued' and count==5
  partial.append({'profile':key,'persisted_blocks':5,'planned_blocks':24,'reason':'ENOSPC during persistence of sixth block; original five blocks retained byte-exact','unpersisted_native_rows':6})
  missing += [key+'/'+mode for mode in ['parse','write','roundtrip']];continue
 assert count==24;complete.append(key)
 for mode in ['parse','write','roundtrip']:
  a=[row['baseline'][mode]['elapsed_ns'] for row in data['blocks']];b=[row['candidate'][mode]['elapsed_ns'] for row in data['blocks']];ratio=paired_median_ratio_ci(a,b,seed=0);median=statistics.median(x/y for x,y in zip(a,b,strict=True))
  assert ratio==data['summary'][mode]['source'] and median==data['summary'][mode]['source_unrounded_median'];source[key+'/'+mode]=ratio
  if median<.99:floors.append({'lane':key+'/'+mode,'speedup':ratio,'unrounded_median':median})
  if ratio['ci95'][1]<1:negative.append({'lane':key+'/'+mode,'speedup':ratio})
  if key+'/'+mode in plan['primary']:primary[key+'/'+mode]=ratio
  if p['kind']=='sdf':external[key+'/'+mode]=data['summary'][mode]['openbabel']
assert len(complete)==9 and len(partial)==1 and len(source)==27 and len(missing)==3 and persisted==1470
assert floors and not (O/'decision.json').exists()
decision={'experiment_complete':False,'local_adoption_gate_passed':False,'reasons':['Two completed SDF source write medians below fixed0.99 floor','v3000-continued stopped during save; three source lanes lack a complete24-block epoch'],'primary':primary,'source_ratios':source,'negative_intervals':negative,'source_floor_failures':floors,'missing_lanes':missing,'completed_profiles':complete,'incomplete_profiles':partial,'source_valid_epochs':9,'planned_native_rows':1584,'persisted_native_rows':1470,'unpersisted_native_rows':6,'never_measured_native_rows':108,'fresh_SDF_Open_Babel':external,'all_4384_preflight_source_snapshots_and_rewrites_exact':True,'interrupted_profile_not_rerun':True,'four_remaining_profiles_measured_for_first_time_with_identical_original_timing_loop':True,'universal_superiority':False}
(O/'decision.json').write_text(json.dumps(decision,indent=2)+'\n')
s=(O/'assemble.py').read_text()
def replace(a,b):
 global s
 assert s.count(a)==1,(a,s.count(a));s=s.replace(a,b)
replace("assert data['complete'] and len(data['blocks'])==24 and data['schema']=='balanced-steady-throughput/v1'","assert data['schema']=='balanced-steady-throughput/v1'\n assert (data['complete'] and len(data['blocks'])==24) if key!='v3000-continued' else (not data['complete'] and len(data['blocks'])==5)")
replace("assert [r['order'] for r in data['blocks']]==p['orders']","assert [r['order'] for r in data['blocks']]==p['orders'][:len(data['blocks'])]")
replace("for pos in range(len(arms)):assert sum(r['order'][pos]==arm for r in data['blocks'])==24//len(arms)","for pos in range(len(arms)):\n   if data['complete']:assert sum(r['order'][pos]==arm for r in data['blocks'])==24//len(arms)")
replace("    samples+=1\n for mode", "    samples+=1\n if not data['complete']:\n  for mode in ['parse','write','roundtrip']:rows.append('| '+key+' | '+mode+' | — | — | incomplete 5/24; no ratio or CI | not measured |')\n  continue\n for mode")
replace('assert samples==1584 and len(audits)==30 and external_count==6','assert samples==1470 and len(audits)==27 and external_count==6')
replace("and not floors)\n(O/'audit.json')","and not floors and not decision['missing_lanes'])\n(O/'audit.json')")
replace("'all_1584_sample_counts_budgets_normalizations_checked':True,'all_30_source_ratios_recomputed':True", "'all_1470_persisted_sample_counts_budgets_normalizations_checked':True,'all_27_completed_source_ratios_recomputed':True,'planned_rows':1584,'unpersisted_rows':6,'never_measured_rows':108,'missing_lanes_unranked':decision['missing_lanes']")
replace("'no_filtered_rows_or_reruns':True", "'no_performance_informed_filtering_or_reruns':True")
replace("if not p.is_file() or p.name=='assemble.log':continue", "if not p.is_file() or p.name in ['assemble.log','assemble-interrupted.log']:continue")
replace('Fresh balanced source experiment: ten profiles, 30 source lanes, 24 blocks;', 'Interrupted experiment: nine complete profiles and one five-block partial profile; 27 completed source lanes, three unranked incomplete lanes;')
replace('Exactly 1584 native samples retain', 'Exactly 1470 persisted native samples retain')
replace('FAIL: outlined candidate not adopted','INCOMPLETE AND FAILED: outlined candidate not adopted')
replace('# Balanced throughput: outlined stereocenter classification','# Interrupted balanced throughput: outlined stereocenter classification')
replace('## Source and correctness', '''## Interruption and missing data

The original run stopped with transient ENOSPC while opening the sixth block result for `v3000-continued`. The original JSON remains readable with five blocks. Six completed native calls in the failed block have no persisted timing values; they are not reconstructed. Of 1,584 planned native rows, 1,470 are retained, six failed persistence and 108 were never measured. Nine profiles have complete 24-block epochs; the three continued-profile lanes remain unranked. This experiment is ineligible for adoption even apart from the failed SDF floors.

After checking process/open-file use, unused owned build caches were removed with all copied formal binaries, source, inputs, logs and data preserved. The exact allocator/process causing the transient allocation shortage was not identified. Only the four profiles never reached were then measured for the first time with the original timing-loop bytes and fixed orders, while all six preexisting raw files were SHA-checked unchanged. The completed profiles and interrupted profile were never rerun. No observations within a profile were pooled across this interruption. The original full-plan assembler is retained; a separate interrupted-report assembler checks the available coverage and forces the overall decision false.

## Source and correctness''')
replace('Orders and position balance are audited.', 'Planned-order prefixes are audited for every saved row; complete-profile position balance is audited. The interrupted five-block prefix is not ranked.')
replace('Source-only comparisons cover all ten;', 'Source-only timing was planned for all ten, but only nine epochs completed;')
replace('All 1,584 samples are independently audited.', 'All 1,470 persisted rows are independently audited; six unavailable values are not reconstructed and 108 planned rows were not measured.')
replace('source proof at current main', 'source proof at the then-current main snapshot')
replace('## Complete matrix','## Complete coverage matrix, including the incomplete profile')
replace('contains every raw native sample', 'contains every persisted raw native sample')
replace('`formal.py`, `assemble.py` after adapting original paths.', '`formal.py` after adapting original paths. The recorded interrupted run additionally uses `resume_remaining.py` for only the four previously unmeasured profiles and `prepare_interrupted_report.py` / `assemble_interrupted.py` to report the incomplete result; the original full-plan `assemble.py` cannot validate this incomplete run.')
replace('24-block source retention across 30 lanes; fresh SDF-only Open Babel comparison.', 'interrupted: nine complete profiles, one five-block partial; failed SDF retention gates, fresh SDF-only Open Babel comparison.')
compile(s,str(O/'assemble_interrupted.py'),'exec');(O/'assemble_interrupted.py').write_text(s)
print(json.dumps({'complete_profiles':len(complete),'incomplete_profiles':partial,'persisted_rows':persisted,'primary':primary,'floor_failures':floors,'negative_intervals':negative,'adoption':False}),flush=True)
