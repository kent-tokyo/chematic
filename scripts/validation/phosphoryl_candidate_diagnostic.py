#!/usr/bin/env python3
"""Independent stereo check of two exposed phosphorus development inputs."""
import gzip
import json
from pathlib import Path
import sys
from rdkit import Chem, rdBase
import chematic
sys.path.insert(0, str(Path.cwd() / 'scripts/benchmarks'))
from platinum_development_probe import stereo_violations
with gzip.open('validation/results/cod-linux-development-v4-2026-11-22.json.gz','rt') as source:
    rows=json.load(source)['results']
names={'VJUNGNXBONODDR-UHFFFAOYSA-N','HRPUQROAFUIJQN-UHFFFAOYSA-N'}
out=[]
for row in rows:
    if row['name'] not in names: continue
    ref=Chem.AddHs(Chem.MolFromSmiles(row['smiles']))
    m=chematic.from_smiles(row['smiles']).add_hydrogens()
    candidates=m.etkdg_embed_multiple_confs_mmff_corrected(20,random_seed=42,max_mmff_iters=200,repair_mmff_iters=400)
    bad=stereo_violations(ref,[coords for _,_,coords in candidates])
    out.append({'name':row['name'],'baseline_count':row['chematic_count'],'candidate_count':len(candidates),
                'unfinished':sum(status != 0 for status,_,_ in candidates),
                'violations':sum(bool(x['atoms'] or x['bonds']) for x in bad),'details':bad})
assert len(out)==2
result={'rdkit':rdBase.rdkitVersion,'results':out}
Path('phosphoryl-candidate-diagnostic.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))

assert all(item["candidate_count"] > 0 and item["violations"] == 0 for item in out)
