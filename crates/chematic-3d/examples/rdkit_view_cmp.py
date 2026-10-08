"""Compare chematic's RDKit view dump with RDKit's AddHs(MolFromSmiles(s))."""
import sys, json
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
bad = {}
n = ok = 0
for line in open(sys.argv[1]):
    d = json.loads(line)
    m = Chem.MolFromSmiles(d["smiles"])
    if m is None:
        continue
    n += 1
    if "error" in d:
        bad.setdefault("error", []).append((d["smiles"], d["error"]))
        continue
    mh = Chem.AddHs(m)
    ri = mh.GetRingInfo()
    ref = {
        "atoms": [[a.GetAtomicNum(), int(a.GetHybridization()), int(a.GetChiralTag()), a.GetTotalNumHs(True)] for a in mh.GetAtoms()],
        "bonds": [[b.GetBeginAtomIdx(), b.GetEndAtomIdx(), int(b.GetBondType()), b.GetIsConjugated(), int(b.GetStereo()), list(b.GetStereoAtoms())] for b in mh.GetBonds()],
        "atom_bonds": [[b.GetIdx() for b in a.GetBonds()] for a in mh.GetAtoms()],
        "atom_rings": [list(r) for r in ri.AtomRings()],
        "bond_rings": [list(r) for r in ri.BondRings()],
    }
    good = True
    for k, v in ref.items():
        if d[k] != v:
            good = False
            bad.setdefault(k, []).append(d["smiles"])
    ok += good
print(f"{ok}/{n}")
for k, v in bad.items():
    print(k, len(v), v[:3])
