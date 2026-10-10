"""Freeze raw cis/trans and fragment ranking inputs with RDKit 2026.03.1."""
import itertools
import json
from pathlib import Path
from rdkit import Chem, rdBase

assert rdBase.rdkitVersion == "2026.03.1", rdBase.rdkitVersion

def describe(mol):
    rings = mol.GetRingInfo()
    atoms = [dict(atomic_num=a.GetAtomicNum(), isotope=a.GetIsotope(),
                  formal_charge=a.GetFormalCharge(), atom_map=a.GetAtomMapNum(),
                  chiral_tag=min(int(a.GetChiralTag()), 3), total_num_hs=a.GetTotalNumHs(),
                  num_rings=rings.NumAtomRings(a.GetIdx()), ring_stereo=False)
             for a in mol.GetAtoms()]
    bonds = [dict(begin=b.GetBeginAtomIdx(), end=b.GetEndAtomIdx(),
                  bond_type=int(b.GetBondType()), stereo=int(b.GetStereo()),
                  stereo_atoms=list(b.GetStereoAtoms()) or None)
             for b in mol.GetBonds()]
    return dict(atoms=atoms, bonds=bonds)

rows = []
for source in ["FC(Cl)=C(Br)I.FC(Cl)=C(Br)I", "C/C=C/C=C/C",
               "F/C=C(/Cl)C=C(/Br)I", "F/C=C1/CCCCCCC1"]:
    base = Chem.MolFromSmiles(source)
    doubles = [b.GetIdx() for b in base.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE]
    assert len(doubles) <= 2
    for states in itertools.product(range(6), repeat=len(doubles)):
        for swap in [False, True]:
            mol = Chem.Mol(base)
            for index, state in zip(doubles, states):
                bond = mol.GetBondWithIdx(index)
                a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                left = [n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() != b]
                right = [n.GetIdx() for n in mol.GetAtomWithIdx(b).GetNeighbors() if n.GetIdx() != a]
                bond.SetStereoAtoms(left[-1 if swap else 0], right[0])
                bond.SetStereo(Chem.BondStereo.values[state])
            rows.append(dict(label=f"{source}/{states}/{swap}", kind="molecule",
                             **describe(mol), ranks=list(Chem.CanonicalRankAtoms(mol))))

for source in ["CC(=O)Oc1ccccc1C(=O)O", "[13CH3][C@@H](F)[NH2+:7]C",
               "C1CC2CCC1C2", "C/C=C(/F)C=C/C"]:
    for chirality in [False, True]:
        mol = Chem.MolFromSmiles(source)
        if not chirality:
            Chem.RemoveStereochemistry(mol)
        n, m = mol.GetNumAtoms(), mol.GetNumBonds()
        for selected in [list(range(n)), list(range(1, n)), list(range(0, n, 2))]:
            bonds = [b.GetIdx() for b in mol.GetBonds()
                     if b.GetBeginAtomIdx() in selected and b.GetEndAtomIdx() in selected]
            # In this pinned Python wrapper, includeAtomMaps is not forwarded:
            # Wrap/rdmolfiles.cpp passes includeChiralPresence in its position.
            # True here selects the C++ includeAtomMaps=true contract that the
            # Rust raw ranking API ports. C++ includeChiralPresence stays false.
            ranks = list(Chem.CanonicalRankAtomsInFragment(mol, selected, bonds,
                         includeChirality=chirality, includeChiralPresence=True))
            rows.append(dict(label=f"{source}/fragment/{selected}/{chirality}", kind="fragment",
                             **describe(mol), ranks=[r & 0xffffffff for r in ranks],
                             atoms_in_play=[i in selected for i in range(n)],
                             bonds_in_play=[i in bonds for i in range(m)], include_chirality=chirality))
    mol = Chem.MolFromSmiles(source)
    symbols = ["same" if i % 2 == 0 else "other" for i in range(mol.GetNumAtoms())]
    ranks = list(Chem.CanonicalRankAtomsInFragment(mol, list(range(mol.GetNumAtoms())),
                 list(range(mol.GetNumBonds())), atomSymbols=symbols, breakTies=False,
                 includeChirality=False, includeIsotopes=False, includeAtomMaps=False))
    rows.append(dict(label=f"{source}/symbols", kind="symbols", **describe(mol),
                     ranks=ranks, atom_symbols=symbols))

root = Path(__file__).resolve().parents[1]
(root / "validation/rdkit-2026.03.1-ranking-boundary.json").write_text(
    json.dumps(dict(rdkit_version=rdBase.rdkitVersion, rows=rows), separators=(",", ":")) + "\n")
print("wrote", len(rows), "cases")
