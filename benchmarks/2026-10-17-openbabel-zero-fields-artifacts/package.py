from pathlib import Path
import gzip,hashlib,json,shutil,subprocess
E=Path(__file__).resolve().parent;R=Path('/Users/k_tanabe/.codex/worktrees/rdkit-cosmolkit-speed-20261017/chematic');A=R/'benchmarks/2026-10-17-openbabel-zero-fields-artifacts';A.mkdir(exist_ok=False)
seal=json.loads((E/'seal.json').read_text())
for path,h in seal['files_sha256'].items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==h,path
for path,h in seal['source_files_sha256'].items():assert hashlib.sha256((R/path).read_bytes()).hexdigest()==h,path
for p in E.iterdir():
    if p.is_file() and p.suffix in ('.json','.jsonl','.py','.log','.txt','.rs','.csv','.gz') and p.name!='baseline-source-tree.tar.gz':shutil.copyfile(p,A/p.name)
    elif p.is_file() and p.suffix in ('.sdf','.molblocks'):(A/(p.name+'.gz')).write_bytes(gzip.compress(p.read_bytes(),compresslevel=6,mtime=0))
    elif p.is_dir() and p.name.endswith('-journals'):shutil.copytree(p,A/p.name)
patch=subprocess.check_output(['git','diff','c84c68f297f95562d1d29c75234f2190cd5f9d75','a0b47272beb93af19bed1f2fe303237affbc0970','--','crates','scripts','tools'],cwd=R)
(A/'candidate.patch').write_bytes(patch)
files={str(p.relative_to(A)):hashlib.sha256(p.read_bytes()).hexdigest() for p in A.rglob('*') if p.is_file()}
local={p.name:{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in E.iterdir() if p.is_file() and (p.name in ['baseline-source-tree.tar.gz','baseline-sdf','baseline-v3000','baseline-v3000_3d','candidate-sdf','candidate-v3000','candidate-v3000_3d','openbabel-sdf','openbabel-v3000'])}
(A/'manifest.json').write_text(json.dumps({'schema':'openbabel-zero-fields-artifacts/v1','base_commit':'c84c68f297f95562d1d29c75234f2190cd5f9d75','candidate_commit':'a0b47272beb93af19bed1f2fe303237affbc0970','independent_local_gate_passed':True,'file_count_excluding_manifest':len(files),'files_sha256':files,'local_native_artifacts_retained_not_vendored':local,'scope':'one macOSarm64 host; native source builds; no whole-toolkit, holdout or released-package claim'},indent=2)+'\n')
print('Packaged',len(files),'files;',sum(p.stat().st_size for p in A.rglob('*') if p.is_file()),'bytes')
