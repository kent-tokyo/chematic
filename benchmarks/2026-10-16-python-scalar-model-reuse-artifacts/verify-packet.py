from pathlib import Path
import gzip, hashlib, importlib, json, re, subprocess, sys, tempfile
ROOT=Path('/private/tmp/chematic-sdf-buffer-scan-20261016')
OUT=Path(__file__).resolve().parent
PACKET=ROOT/'benchmarks/2026-10-16-python-scalar-model-reuse-artifacts'
sys.path.insert(0,str(ROOT/'scripts'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
manifest=json.loads((OUT/'source-manifest.json').read_text())
assert all(sha(ROOT/p)==v for p,v in manifest['candidate_sha256'].items())
seal=json.loads((OUT/'protocol-seal.json').read_text())
assert all(sha(OUT/p)==v for p,v in seal.items())
changed=manifest['changed_build_files'];added=manifest['added_files']
baseline=json.loads((OUT/'baseline-source-manifest.json').read_text())
assert all(sha(OUT/'baseline-source'/p)==v for p,v in baseline['sha256'].items())
# Every frozen candidate source/input must be a base blob or a named change/addition.
paths=list(manifest['candidate_sha256'])
raw=subprocess.check_output(['git','cat-file','--batch'],input=''.join(manifest['base']+':'+p+'\n' for p in paths).encode(),cwd=ROOT)
at=0
for path in paths:
 end=raw.index(b'\n',at);header=raw[at:end].split();at=end+1
 if header[-1]==b'missing':
  assert path in added,path
  continue
 size=int(header[2]);blob=raw[at:at+size];at+=size+1
 if path not in changed:
  assert hashlib.sha256(blob).hexdigest()==manifest['candidate_sha256'][path],path
assert at==len(raw)
base_model=subprocess.check_output(['git','show',manifest['base']+':crates/chematic-smiles/src/rdkit/model_view.rs'],cwd=ROOT,text=True)
original=base_model[base_model.index('fn with_rdkit_model_molecule_impl<R>('):base_model.index('/// Descriptor/fingerprint callbacks')].strip()
frozen=(ROOT/'crates/chematic-smiles/src/rdkit/model_view_cache_tests.rs').read_text()
frozen=frozen[frozen.index('fn frozen_model_callback<R>('):frozen.index('fn snapshot(')].strip().replace('frozen_model_callback','with_rdkit_model_molecule_impl')
assert re.sub(r'\s+','',original)==re.sub(r'\s+','',frozen)
with tempfile.TemporaryDirectory(prefix='chematic-scalar-model-patch-',dir='/private/tmp') as directory:
 tmp=Path(directory)
 subprocess.run(['git','init','-q',str(tmp)],check=True)
 for path in changed:
  original=subprocess.check_output(['git','show',manifest['base']+':'+path],cwd=ROOT)
  p=tmp/path;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(original)
 subprocess.run(['git','apply',str(OUT/'candidate.patch')],cwd=tmp,check=True)
 assert all(sha(tmp/p)==manifest['candidate_sha256'][p] for p in changed+added)
checks=0;strict=0;unranked=[];primary=[];medians=[];negative={}
for lane in ['frontier','bundle']:
 bench=importlib.import_module('bench_python_scalar_bundle' if lane=='bundle' else 'bench_python_speed_frontier')
 for offset in [0,5000]:
  profile=f'{lane}-{offset}'
  data=json.loads((OUT/f'{profile}.json').read_text())
  assert data['complete'] and len(data['blocks'])==20
  assert data['base_revision']==manifest['base'] and data['source_manifest']==manifest
  assert data['corpus']['sha256']==sha(ROOT/'validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi')
  observers={}
  for arm in ['baseline','candidate','rdkit','cosmolkit']:
   observer=json.loads(gzip.decompress((OUT/f'observer-{profile}-{arm}.json.gz').read_bytes()))
   values=observer['values'];observers[arm]=values
   assert hashlib.sha256(json.dumps(values,sort_keys=True).encode()).hexdigest()==observer['output_sha256']
   assert all(b[arm]['output_sha256']==observer['output_sha256'] and b[arm]['extension_sha256']==observer['extension_sha256'] and b[arm]['version']==observer['version'] for b in data['blocks'])
  negative[profile]={}
  for op in bench.OPS:
   agreement={'candidate_vs_baseline':bench.agreement(observers['baseline'][op],observers['candidate'][op],exact=True)}
   agreement.update({name+'_vs_rdkit':bench.agreement(observers['rdkit'][op],observers[name][op],exact=False) for name in ['baseline','candidate','cosmolkit']})
   assert agreement==data['agreement'][op]
   assert agreement['candidate_vs_baseline']['passed']
   for mode in bench.MODES:
    key=op+'/'+mode
    for seed,arm in enumerate(['baseline','rdkit','cosmolkit']):
     value=data['summary'][key][arm]
     allowed=bench.comparable(agreement,arm)
     expected=bench.paired_median_ratio_ci([b[arm]['operations'][op][mode] for b in data['blocks']],[b['candidate']['operations'][op][mode] for b in data['blocks']],seed=seed) if allowed else None
     assert value['speedup']==expected,(profile,key,arm)
     win=bool(expected and expected['ci95'][0]>1 and expected['winning_blocks']==20)
     assert value['demonstrated_win']==win
     checks+=1;strict+=int(win)
     if not allowed:unranked.append([profile,key,arm])
     if arm=='baseline':
      medians.append(expected['median'])
      if lane=='bundle' and mode=='first_use':primary.append(expected['ci95'][0])
      if expected['ci95'][1]<1:negative[profile][key]=expected
parity=json.loads((OUT/'parity.json').read_text())
assert parity['outputs']['baseline']==parity['outputs']['candidate']
for arm,records in parity['outputs'].items():
 for kind,record in records.items():assert hashlib.sha256(gzip.decompress((OUT/f'{arm}-{kind}.jsonl.gz').read_bytes())).hexdigest()==record['sha256']
decision=json.loads((OUT/'decision.json').read_text())
assert len(medians)==48 and len(primary)==2 and checks==144
assert decision['adoption_passed']==(all(x>1 for x in primary) and all(x>=.99 for x in medians))
assert decision['all_negative_source_intervals']==negative
proof={'source_files_verified':len(manifest['candidate_sha256']),'only_tracked_build_snapshot_files':True,'all_unnamed_source_files_match_base_blobs':True,'frozen_original_callback_reverified_from_git':True,'protocol_seal_verified':True,'patch_reconstructed_from_base':manifest['base'],'changed_and_added_files_verified':changed+added,'all_matrix_comparisons_recomputed':checks,'source_median_guards_verified':len(medians),'strict_wins_including_source':strict,'all_unranked_lanes':unranked,'all_20_hashes_match_fresh_observers':True,'full_source_dump_hashes_verified':True,'adoption_passed':decision['adoption_passed'],'single_complete_epochs_no_reruns':True,'universal_superiority':False}
(OUT/'packet-audit.json').write_text(json.dumps(proof,indent=2)+'\n')
(PACKET/'packet-audit.json').write_bytes((OUT/'packet-audit.json').read_bytes())
# Final packet contents may include this checker and the final report assembly driver.
for p in OUT.glob('*.py'):(PACKET/p.name).write_bytes(p.read_bytes())
files={p.name:sha(p) for p in sorted(PACKET.iterdir()) if p.is_file() and p.name!='packet.sha256.json'}
(PACKET/'packet.sha256.json').write_text(json.dumps(files,indent=2)+'\n')
assert all(sha(PACKET/p)==v for p,v in files.items())
print(json.dumps(proof,indent=2))
