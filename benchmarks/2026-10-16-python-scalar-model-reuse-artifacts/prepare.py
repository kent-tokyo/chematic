from pathlib import Path
import hashlib, json, shutil, subprocess
ROOT=Path('/private/tmp/chematic-sdf-buffer-scan-20261016')
OUT=Path(__file__).resolve().parent
base=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
assert base=='2251a9f5d2e706087efaf55cb0f8ab1ff20b8016'
tracked=subprocess.check_output(['git','ls-files','-z'],cwd=ROOT).decode().split('\0')
paths=sorted(p for p in tracked if p and (p.startswith(('crates/','data/','tools/')) or '/' not in p) and (ROOT/p).is_file())
assert all('/__pycache__/' not in p for p in paths)
snapshot=OUT/'baseline-source';snapshot.mkdir()
manifest={}
for path in paths:
 target=snapshot/path;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/path,target)
 manifest[path]=hashlib.sha256(target.read_bytes()).hexdigest()
raw=subprocess.check_output(['git','cat-file','--batch'],cwd=ROOT,input=''.join(base+':'+p+'\n' for p in paths).encode());at=0
for p in paths:
 end=raw.index(b'\n',at);header=raw[at:end].decode().split();assert len(header)==3,p;size=int(header[2]);at=end+1
 assert hashlib.sha256(raw[at:at+size]).hexdigest()==manifest[p],p;at+=size+1
assert at==len(raw)
(OUT/'baseline-source-manifest.json').write_text(json.dumps({'base':base,'sha256':manifest,'only_git_tracked_files':True},indent=2)+'\n')
plan={'base':base,'candidate':'Memoize the sanitized core model only for scalar descriptor callbacks; original generic model callbacks remain uncached; graph edits/clones clear the new derived slot','primary':['tpsa_then_asa/first_use at offsets 0 and 5000'],'bundle_operations':['tpsa','labute_asa'],'bundle_modes':['first_use','repeated_use','parse_inclusive'],'corpus':'validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi','offsets':[0,5000],'rows_each':1000,'blocks':20,'operations':['tpsa','labute_asa','num_rings','morgan2_2048','morgan2_chiral','morgan3_2048','maccs'],'modes':['first_use','repeated_use','parse_inclusive'],'adoption':'Both bundle first-use paired CI lower>1; all 14072 source outputs exact; all 48 source-relative medians >=0.99; disclose all negative intervals; no adaptive timing reruns','comparator_win':'Equivalent output, paired CI lower>1 and 20/20 faster blocks; mismatches unranked','build':{'opt_level':3,'lto':True,'codegen_units':1,'jobs':1,'incremental':False,'features':['native-inchi'],'strip':'none'},'scope':'One shared macOS arm64 host, exposed development corpus, source wheels; not toolkit-wide or published-package superiority','memory':'Model copies are retained until their original molecule is dropped; whole-process RSS recorded, no operation-memory or universal-memory win required or claimed','epochs':'One complete 20-block matrix per slice for ordinary operations and the descriptor bundle, with no timing pooling or substitution'}
(OUT/'measurement-plan.json').write_text(json.dumps(plan,indent=2)+'\n')
print('Prepared',len(paths),'Git-verified baseline files; no ignored caches copied')
