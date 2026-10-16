from pathlib import Path
import gzip,hashlib,json,os,shutil,subprocess
ROOT=Path('/private/tmp/chematic-3d-stereo-guard-20261017');OUT=Path(__file__).resolve().parent
MEASURED=Path('/private/tmp/chematic-3d-stereo-guard-valid-evidence-20261017')
env=os.environ.copy();env.update(json.loads((OUT/'native-check-commands.json').read_text())['env'])
cmd=['cargo','build','--offline','--locked','--release','-p','chematic-mol','--example','openbabel_sdf_throughput','--example','openbabel_v3000_throughput','--example','openbabel_v3000_3d_throughput']
(OUT/'probe-build-command.json').write_text(json.dumps({'command':cmd,'env':{k:v for k,v in env.items() if k.startswith('CARGO_')}},indent=2)+'\n')
with (OUT/'probe-build.log').open('w') as log:subprocess.run(cmd,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
proof={'head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'scope':'Untimed combined-head preservation only; no new timing observations','profiles':{}}
for kind in ['sdf','v3000','v3000-3d']:
 shutil.copy2(Path(env['CARGO_TARGET_DIR'])/'release/examples'/('openbabel_'+kind.replace('-','_')+'_throughput'),OUT/('integrated-'+kind))
for item in json.loads((MEASURED/'measurement-plan.json').read_text())['profiles']:
 kind=item['kind'];key=kind+'-'+item['profile'];corpus=MEASURED/(key+('.sdf' if kind=='sdf' else '.molblocks'));values={}
 for mode in ['snapshot','dump']:
  raw=subprocess.check_output([str(OUT/('integrated-'+kind)),str(corpus),mode,'1',str(item['records'])])
  expected=gzip.decompress((MEASURED/f'preflight-{key}-candidate-{mode}.gz').read_bytes());assert raw==expected,(key,mode)
  (OUT/f'{key}-{mode}.gz').write_bytes(gzip.compress(raw,mtime=0));values[mode+'_sha256']=hashlib.sha256(raw).hexdigest()
 proof['profiles'][key]={'records':item['records'],'both_outputs_exact':True,**values}
proof['all_4256_profile_inputs_preserved']=sum(v['records'] for v in proof['profiles'].values())==4256
assert proof['all_4256_profile_inputs_preserved']
(OUT/'output-equality.json').write_text(json.dumps(proof,indent=2)+'\n');print('All 4256 combined-head snapshots and rewrites exact')
