"""Dump plain-DG embedding inputs/outputs from RDKit for the Rust checker."""
import sys, json
from rdkit import Chem, RDLogger
from rdkit.Chem import rdDistGeom
RDLogger.DisableLog("rdApp.*")
path, n = sys.argv[1], int(sys.argv[2])
rows = [l.split()[0] for l in open(path) if l.strip()][:n]
for smi in rows:
    m0 = Chem.MolFromSmiles(smi)
    if m0 is None:
        continue
    mh = Chem.AddHs(m0)
    if len(Chem.GetMolFrags(mh)) != 1:
        continue
    bm = rdDistGeom.GetMoleculeBoundsMatrix(mh)
    ps = rdDistGeom.EmbedParameters()
    ps.useExpTorsionAnglePrefs = False
    ps.useBasicKnowledge = False
    ps.enforceChirality = True
    ps.randomSeed = 42
    ps.SetBoundsMat(bm)
    piece = Chem.Mol(mh)
    Chem.AssignStereochemistry(piece)
    ri = piece.GetRingInfo()
    chiral, tet = [], []
    for at in piece.GetAtoms():
        if at.GetAtomicNum() == 1:
            continue
        ct = at.GetChiralTag()
        isch = ct in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        if not (isch or (at.GetAtomicNum() in (6, 7) and at.GetDegree() == 4)):
            continue
        nb = [b.GetOtherAtomIdx(at.GetIdx()) for b in at.GetBonds()]
        lo, up = 5.0, 100.0
        if len(nb) < 4:
            lo = 2.0
            nb.append(at.GetIdx())
        fused = sum(1 for s in ri.AtomRingSizes(at.GetIdx()) if s < 5) > 1
        idx = [at.GetIdx()] + nb
        if ct == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            chiral.append([idx, lo, up, fused])
        elif ct == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            chiral.append([idx, -up, -lo, fused])
        else:
            if ri.NumAtomRings(at.GetIdx()) < 2 or ri.IsAtomInRingOfSize(at.GetIdx(), 3):
                continue
            tet.append([idx, 0.0, 0.0, fused])
    dbe, sdb = [], []
    for b in piece.GetBonds():
        if b.GetBondType() != Chem.BondType.DOUBLE:
            continue
        for atm in (b.GetBeginAtom(), b.GetEndAtom()):
            if atm.GetDegree() < 2:
                continue
            o = b.GetOtherAtom(atm)
            for nbr in atm.GetNeighbors():
                if nbr.GetIdx() == o.GetIdx():
                    continue
                ob = piece.GetBondBetweenAtoms(atm.GetIdx(), nbr.GetIdx())
                if ob.GetBondType() != Chem.BondType.SINGLE and atm.GetDegree() == 2:
                    continue
                dbe.append([nbr.GetIdx(), atm.GetIdx(), o.GetIdx()])
        if b.GetStereo() > Chem.BondStereo.STEREOANY:
            sign = -1 if b.GetStereo() in (Chem.BondStereo.STEREOCIS, Chem.BondStereo.STEREOZ) else 1
            sa = list(b.GetStereoAtoms())
            sdb.append([[sa[0], b.GetBeginAtomIdx(), b.GetEndAtomIdx(), sa[1]], sign])
    cid = rdDistGeom.EmbedMolecule(mh, ps)
    coords = mh.GetConformer().GetPositions().flatten().tolist() if cid == 0 else None
    N = mh.GetNumAtoms()
    print(json.dumps({"smiles": smi, "n": N, "bounds": [repr(x) for x in bm.flatten().tolist()], "coords": None if coords is None else [repr(x) for x in coords],
                      "chiral": chiral, "tet": tet, "dbe": dbe, "sdb": sdb}))
