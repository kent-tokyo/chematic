from pathlib import Path
import hashlib,json,shutil,subprocess
ROOT=Path('/private/tmp/chematic-sdf-buffer-scan-20261016');OUT=Path(__file__).resolve().parent
BASE=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip();assert BASE=='faa8092e6050ad6250c68498970245cfbc6579e9'
# Same instrumented 3D API probe is compiled against both source trees.
probe=(ROOT/'crates/chematic-mol/examples/openbabel_v3000_throughput.rs').read_text()
probe=probe.replace('//! Native V3000 MOL API throughput over a batch split before timing.','//! Native 3D V3000 MOL APIs; fail closed if a parsed conformer is absent.')
probe=probe.replace('write_mol_v3000};','write_mol_v3000_with_conformer};')
probe=probe.replace('write_mol_v3000(&r.mol, &r.metadata, &r.coords)','write_mol_v3000_with_conformer(&r.mol, &r.metadata, r.conformer.as_ref().expect("3D corpus requires conformer"))')
name='crates/chematic-mol/examples/openbabel_v3000_3d_throughput.rs';(ROOT/name).write_text(probe)
paths=subprocess.check_output(['git','ls-files'],cwd=ROOT,text=True).splitlines()
paths=[p for p in paths if (p.startswith(('crates/','data/','tools/')) or '/' not in p) and (ROOT/p).is_file() and p!='.git']+[name]
baseline=OUT/'baseline-source';baseline.mkdir(exist_ok=True)
for name in paths:
 dest=baseline/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,dest)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
(OUT/'baseline-source-manifest.json').write_text(json.dumps({'base_revision':BASE,'source_files_sha256':{p:sha(baseline/p) for p in sorted(paths)},'only_git_tracked_plus_identical_3d_probe':True,'instrumentation':[name]},indent=2)+'\n')
old=Path('/private/tmp/chematic-v3000-record-stream-evidence-20261016')
for name in ['sdf-0.sdf','sdf-5000.sdf','v3000-0.molblocks','v3000-5000.molblocks','v3000-chains.molblocks','v3000-continued.molblocks','openbabel-sdf','openbabel-v3000']:
 shutil.copy2(old/name,OUT/name)
# Explicit controlled 3D files: no embedding/energy or realistic conformer claim.
recipe={}
for profile in ['spiral','offset']:
 source=(OUT/'v3000-chains.molblocks').read_text();lines=[];record=0;in_atoms=False;atoms=0
 for line in source.splitlines(keepends=True):
  if line=='  chematic\n' or line.lstrip().startswith('RDKit'):
   line='  chematic          3D\n'
  if line.rstrip('\n')=='M  V30 BEGIN ATOM':in_atoms=True
  elif line.rstrip('\n')=='M  V30 END ATOM':in_atoms=False
  elif line.rstrip('\n')=='$$$$':record+=1
  if in_atoms and line.startswith('M  V30 ') and line[7:8].isdigit():
   words=line[7:].split();i=int(words[0]);
   if profile=='spiral':x=i*.1251;y=((i*37)%101-50)*.1234;z=((i*73)%97-48)*.2345
   else:x=(i-512)*1.2345;y=((i*19)%53-26)*.3456;z=-((i*43)%79-39)*.4567
   words[2:5]=[f'{x:.4f}',f'{y:.4f}',f'{z:.4f}'];line='M  V30 '+' '.join(words)+'\n';atoms+=1
  lines.append(line)
 assert record==64
 payload=''.join(lines)
 # Confirm the source dimensional marker instead of silently projecting 2D.
 assert payload.count('          3D\n')==64
 target=OUT/f'v3000-3d-{profile}.molblocks';target.write_text(payload)
 recipe[profile]={'records':record,'atoms':atoms,'sha256':sha(target),'source_sha256':sha(OUT/'v3000-chains.molblocks'),'coordinate_formula':'x=i*.1251;y=((i*37)%101-50)*.1234;z=((i*73)%97-48)*.2345' if profile=='spiral' else 'x=(i-512)*1.2345;y=((i*19)%53-26)*.3456;z=-((i*43)%79-39)*.4567'}
(OUT/'controlled-3d-recipe.json').write_text(json.dumps({'scope':'Synthetic finite nonplanar coordinates on pinned achiral linear chains, not physical conformer/embedding validation','profiles':recipe},indent=2)+'\n')
profiles=[{'kind':k,'profile':p,'records':n} for k,p,n in [('sdf','0',1000),('sdf','5000',1000),('v3000','0',1000),('v3000','5000',1000),('v3000','chains',64),('v3000','continued',64),('v3000-3d','spiral',64),('v3000-3d','offset',64)]]
plan={'base_revision':BASE,'candidate':'Use existing byte-compatible fixed4 and integer append helpers in 3D V3000 atom rows only; unchanged fallback, other writer and parser paths','primary':'prepared 3D V3000 write CI lower>1 on both controlled 64-record coordinate profiles','adoption':'All source snapshots/rewrites and frozen formatter output bytes exact; both primary write CI lower>1; all 24 source-relative medians>=0.99; disclose every negative interval','profiles':profiles,'blocks':21,'native_repeats':3,'comparator_win':'semantic prerequisite, paired CI lower>1 and 21/21 faster blocks','build':{'opt_level':3,'lto':True,'codegen_units':1,'jobs':1,'incremental':False,'strip':'none'},'test_build':{'release':True,'lto':False,'codegen_units':16,'jobs':1},'no_adaptive_timing_reruns':True,'scope':'Shared macOS arm64 host, exposed development corpora; primary profiles are synthetic achiral chains with finite nonplanar coordinates; source 3D writer only, not full Open Babel or realistic conformer superiority','previous_unpromoted_candidates':'excluded','checkpoint_policy':'64MiB owned reserve; retry same checkpoint, no worker/timer restarts'}
(OUT/'measurement-plan.json').write_text(json.dumps(plan,indent=2)+'\n');print('Prepared',len(paths),'tracked/instrumentation inputs; frozen 8-profile plan')
