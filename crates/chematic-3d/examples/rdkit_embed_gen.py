"""RDKit reference: EmbedMolecule(AddHs(m), randomSeed=42) coordinates.
usage: gen_embed.py SMI N MODE   (MODE: dg = no ET/K, etkdg = defaults)"""
import sys, json
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")
path, n, mode = sys.argv[1], int(sys.argv[2]), sys.argv[3]
for smi in [l.split()[0] for l in open(path) if l.strip()][:n]:
    m = Chem.MolFromSmiles(smi)
    if m is None:
        continue
    mh = Chem.AddHs(m)
    kw = {"randomSeed": 42}
    if mode == "dg":
        kw.update(useExpTorsionAnglePrefs=False, useBasicKnowledge=False)
    try:
        cid = AllChem.EmbedMolecule(mh, **kw)
        coords = [repr(float(x)) for x in mh.GetConformer().GetPositions().flatten()] if cid == 0 else None
        err = None
    except Exception as e:
        coords, err = None, repr(e)[:200]
    print(json.dumps({"smiles": smi, "mode": mode, "coords": coords, "error": err}))
