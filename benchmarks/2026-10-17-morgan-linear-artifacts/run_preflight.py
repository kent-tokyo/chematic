from pathlib import Path
import hashlib, importlib.metadata, json, os, subprocess, sys, zipfile
E=Path(__file__).resolve().parent
R=Path('/Users/k_tanabe/.codex/worktrees/rdkit-cosmolkit-speed-20261017/chematic')
meta={};wheels={}
for name in ('baseline','candidate'):
    wheel,=list((E/(name+'-wheel')).glob('*.whl'))
    assert (E/(name+'-site')).is_dir()
    wheels[name]={'path':str(wheel),'sha256':hashlib.sha256(wheel.read_bytes()).hexdigest()}
    with zipfile.ZipFile(wheel) as archive:
        wheels[name]['native_sha256']={Path(p).name:hashlib.sha256(archive.read(p)).hexdigest() for p in archive.namelist() if p.endswith(('.so','.dylib'))}
assert wheels['baseline']['native_sha256']!=wheels['candidate']['native_sha256']
for offset in (0,5000):
    outputs={}
    for name in ('baseline','candidate','rdkit','cosmolkit'):
        env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1'
        if name in wheels:env['PYTHONPATH']=str(E/(name+'-site'))
        else:env.pop('PYTHONPATH',None)
        cmd=[sys.executable,str(R/'scripts/bench_python_morgan_linear.py'),'--root',str(R),'--worker',name,'--probe','--corpus',str(R/'validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi'),'--offset',str(offset),'--limit','1000']
        run=subprocess.run(cmd,env=env,capture_output=True,check=True)
        (E/f'probe-{name}-{offset}.json').write_bytes(run.stdout)
        sample=json.loads(run.stdout);assert sample['probe_only'] and sample['rows']==1000
        metadata={k:sample[k] for k in ['version','module_path','native_sha256']}
        if name in meta:assert meta[name]==metadata
        meta[name]=metadata;outputs[name]=sample
        if name in wheels:assert sample['native_sha256']==wheels[name]['native_sha256']
    for op in ('morgan2_2048','morgan2_chiral','morgan3_2048'):
        assert all(outputs[name]['operations'][op]['values']==outputs['rdkit']['operations'][op]['values'] for name in outputs),op
    for op in ('tpsa','labute_asa','num_rings','maccs'):
        assert outputs['baseline']['operations'][op]['values']==outputs['candidate']['operations'][op]['values']
assert meta['rdkit']['version']=='2026.09.1'
assert meta['cosmolkit']['version']=='0.5.0'
runtime={};distributions={}
for name in ('rdkit','cosmolkit','numpy'):
    dist=importlib.metadata.distribution(name);count=0
    for file in dist.files or []:
        if str(file).endswith(('.pyc','.pyo')) or '__pycache__' in str(file):continue
        path=Path(dist.locate_file(file)).resolve()
        if path.is_file():runtime[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest();count+=1
    distributions[name]={'version':dist.version,'checked_files':count}
files=subprocess.check_output(['git','ls-files'],cwd=R,text=True).splitlines()
source={p:hashlib.sha256((R/p).read_bytes()).hexdigest() for p in files if (R/p).is_file() and (p.startswith('crates/') or p in ['Cargo.toml','Cargo.lock','scripts/bench_python_morgan_linear.py','scripts/bench_python_speed_frontier.py','scripts/verify_python_morgan_frontier.py','validation/cosmolkit_comparison/bench_corpus.py','validation/cosmolkit_comparison/corpus_engines.py'])}
build=json.loads((E/'candidate-build-source.json').read_text())
assert all(source[p]==h for p,h in build['source_files_sha256'].items())
manifest={'baseline_commit':'c84c68f297f95562d1d29c75234f2190cd5f9d75','candidate_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),'candidate_source_sha256':source,'wheels':wheels,'driver_sha256':hashlib.sha256((R/'scripts/bench_python_morgan_linear.py').read_bytes()).hexdigest(),'runtime_files_sha256':runtime,'distributions':distributions,'expected_metadata':meta,'preflight':{'both_profiles_all_four_arms_all_modes_passed':True,'elapsed_samples_excluded_from_formal':True},'measurement_plan_sha256':hashlib.sha256((E/'measurement-plan.json').read_bytes()).hexdigest()}
(E/'source-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Preflight passed; sealed',len(runtime),'runtime files',flush=True)
