from pathlib import Path
import hashlib, json, os, shutil, subprocess

ROOT = Path('/private/tmp/chematic-sdf-buffer-scan-20261016')
OUT = Path(__file__).resolve().parent
BASELINE = OUT / 'baseline-source'
BASE = subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
assert BASE == '3e7f04b20593af4b07569d886b56afc1184468fe'
BASELINE.mkdir()
for name in ['crates','data']:
    shutil.copytree(ROOT / name, BASELINE / name)
for p in ROOT.iterdir():
    if p.is_file() and p.name != '.git': shutil.copy2(p, BASELINE / p.name)
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
paths = sorted(p.relative_to(BASELINE).as_posix() for p in BASELINE.rglob('*') if p.is_file())
manifest = {p:sha(BASELINE / p) for p in paths}
(OUT / 'baseline-source-manifest.json').write_text(json.dumps({'base':BASE,'sha256':manifest},indent=2)+'\n')
plan = {
 'base':BASE,
 'candidate':'Inline small aromatic-model gate scratch buffers; preserve graph traversal, pi-electron rules, BFS queue order and large-input heap fallback',
 'primary':['tpsa/first_use and labute_asa/first_use at offsets 0 and 5000'],
 'corpus':'validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi','offsets':[0,5000],'rows_each':1000,'blocks':20,
 'operations':['tpsa','labute_asa','num_rings','morgan2_2048','morgan2_chiral','morgan3_2048','maccs'],
 'modes':['first_use','repeated_use','parse_inclusive'],
 'adoption':'All four primary paired 95% CI lower >1; exact source outputs on all 14072 pinned inputs; all 42 source-relative medians >=0.99; retain all negative intervals; no adaptive timing reruns',
 'comparator_win':'Equivalent output, paired CI lower >1 and 20/20 wins; mismatches unranked',
 'build':{'opt_level':3,'lto':True,'codegen_units':1,'jobs':1,'incremental':False,'features':['native-inchi'],'strip':'none'},
 'scope':'One shared macOS arm64 host, exposed development corpus, fresh source wheels; not published or universal superiority',
}
(OUT / 'measurement-plan.json').write_text(json.dumps(plan,indent=2)+'\n')
print('Baseline frozen',len(paths),'files',flush=True)
