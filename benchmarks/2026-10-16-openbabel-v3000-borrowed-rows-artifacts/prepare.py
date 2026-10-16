from pathlib import Path
import gzip, hashlib, json, shutil, subprocess
ROOT=Path('/private/tmp/chematic-sdf-buffer-scan-20261016')
OUT=Path(__file__).resolve().parent
OLD=Path('/private/tmp/chematic-v3000-loss-work-evidence-20261016')
BASE=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
assert BASE=='2251a9f5d2e706087efaf55cb0f8ab1ff20b8016'
baseline=OUT/'baseline-source';baseline.mkdir()
for name in ['crates','data','tools']:shutil.copytree(ROOT/name,baseline/name)
for p in ROOT.iterdir():
 if p.is_file() and p.name!='.git':shutil.copy2(p,baseline/p.name)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
source={p.relative_to(baseline).as_posix():sha(p) for p in sorted(baseline.rglob('*')) if p.is_file()}
(OUT/'baseline-source-manifest.json').write_text(json.dumps({'base_revision':BASE,'source_files_sha256':source},indent=2)+'\n')
for name in ['sdf-0.sdf','sdf-5000.sdf','v3000-0.molblocks','v3000-5000.molblocks','v3000-chains.molblocks','openbabel-sdf','openbabel-v3000']:
 shutil.copy2(OLD/name,OUT/name)
original=(OUT/'v3000-chains.molblocks').read_text()
continued=[];in_atoms=False;changes=0
for line in original.splitlines(keepends=True):
 if line.rstrip('\n')=='M  V30 BEGIN ATOM':in_atoms=True
 elif line.rstrip('\n')=='M  V30 END ATOM':in_atoms=False
 if in_atoms and line.startswith('M  V30 ') and line[7:8].isdigit():
  words=line[7:].split();assert len(words)>=6
  line='M  V30 '+' '.join(words[:3])+' -\nM  V30 '+' '.join(words[3:])+'\n';changes+=1
 continued.append(line)
(OUT/'v3000-continued.molblocks').write_text(''.join(continued))
(OUT/'continuation-corpus.json').write_text(json.dumps({'recipe':'Split every atom row after the x coordinate with exact V30 continuation; no field changes','source_sha256':sha(OUT/'v3000-chains.molblocks'),'continued_sha256':sha(OUT/'v3000-continued.molblocks'),'changed_atom_rows':changes,'records':64},indent=2)+'\n')
profiles=[{'kind':k,'profile':p,'records':n} for k,p,n in [('sdf','0',1000),('sdf','5000',1000),('v3000','0',1000),('v3000','5000',1000),('v3000','chains',64),('v3000','continued',64)]]
plan={'base_revision':BASE,'candidate':'Borrow plain V30 rows directly after unchanged physical validation; retain original collection for continuations and exact allocation hints','primary':'parse CI lower>1 on both unfiltered V3000 slices at offsets 0 and 5000','adoption':'Exact source snapshots/rewrites and frozen full reports/errors; both primary parse paired CI lower>1; all 18 source-relative medians >=0.99; disclose all negative intervals','profiles':profiles,'blocks':21,'native_repeats':3,'comparator_win':'semantic prerequisite, paired CI lower>1 and 21/21 faster blocks','build':{'release':True,'opt_level':3,'lto':True,'codegen_units':1,'jobs':1,'incremental':False,'strip':'none'},'test_build':{'release':True,'lto':False,'codegen_units':16,'jobs':1,'strip':'none'},'no_adaptive_timing_reruns':True,'scope':'Shared macOS arm64 host, exposed development corpora, current native source; no toolkit-wide or publication claim','previous_unpromoted_candidates':'excluded','checkpoint_policy':'64MiB owned reserve; retry identical in-memory checkpoint only, never worker or timer'}
(OUT/'measurement-plan.json').write_text(json.dumps(plan,indent=2)+'\n')
print('Prepared fresh baseline',len(source),'files; continued atom rows',changes)
