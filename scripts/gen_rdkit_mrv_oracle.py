#!/usr/bin/env python3
"""Real RDKit oracle + fixture generator for the IO-3 (MRV) acceptance gate.

Unlike IO-1/IO-2 (chematic-authored multi-record files with RDKit run
against them), MRV fixtures are individual single-molecule files that only
RDKit's own `Chem.MolToMrvBlock` can author at scale -- so this script both
generates the `.mrv` fixture files (via RDKit) AND records RDKit's own
ground truth for each (atom/bond counts, elements, charges, isotopes,
atom-maps, RDKit's own canonical SMILES, 2D/3D coordinates), rather than
splitting fixture-authoring and oracle-running into two scripts.

Draws real, chemically diverse SMILES from the Morgan M4-A0 corpus for the
bulk "general" category, plus small hand-picked scenarios for isotope/
charge/radical/stereo/atom-map/3D/disconnected-fragment coverage that a
random corpus draw would under-represent. Only molecules from chematic's
OWN supported-feature set are included -- S-groups/polymers/reactions/
R-groups/enhanced-stereo/query atoms are never in this fixture set, since
chematic's port deliberately errors on them by design (see
`chematic_mol::mrv`'s module docs) and RDKit-succeeds-where-chematic-errors
would be a documented divergence, not something to gate on here.

Usage:
    python scripts/gen_rdkit_mrv_oracle.py --corpus <SMILES.csv> \\
        --out-dir <fixtures_dir> --manifest-out <manifest.json>
"""

from __future__ import annotations

import argparse
import json
import sys

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDepictor
except ImportError:
    print("rdkit is required to author fixtures", file=sys.stderr)
    raise


def load_corpus(path, limit=None):
    smis = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s:
                smis.append(s)
    if limit:
        smis = smis[:limit]
    return smis


def rdkit_mol(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        return mol
    except Exception:
        return None


def build_entries(corpus):
    valid = [s for s in corpus if rdkit_mol(s) is not None]
    entries = []

    # -- bulk "general" draw: real chemical diversity ------------------------
    general = valid[:160]
    for i, smi in enumerate(general):
        entries.append({"id": f"general_{i}", "category": "general", "smiles": smi, "with_3d": False})

    # -- hand-picked category coverage ---------------------------------------
    acyclic = ["CCCCCC", "CC(C)CC(=O)O", "CCOCC", "CC#CC", "CCN(CC)CC"]
    for i, smi in enumerate(acyclic):
        entries.append({"id": f"acyclic_{i}", "category": "acyclic", "smiles": smi, "with_3d": False})

    aromatic = ["c1ccccc1", "c1ccncc1", "c1ccoc1", "c1cc[nH]c1", "c1ccc2ccccc2c1"]
    for i, smi in enumerate(aromatic):
        entries.append({"id": f"aromatic_{i}", "category": "aromatic", "smiles": smi, "with_3d": False})

    fused_ring = ["c1ccc2ccccc2c1", "c1ccc2[nH]ccc2c1", "C1CCC2CCCCC2C1", "c1ccc2c(c1)ccc1ccccc12"]
    for i, smi in enumerate(fused_ring):
        entries.append({"id": f"fused_ring_{i}", "category": "fused_ring", "smiles": smi, "with_3d": False})

    isotope = ["[13CH4]", "[2H]C([2H])([2H])O", "c1ccc(cc1)[15NH2]", "[13cH]1ccccc1"]
    for i, smi in enumerate(isotope):
        entries.append({"id": f"isotope_{i}", "category": "isotope", "smiles": smi, "with_3d": False})

    charge = ["[NH4+]", "CC(=O)[O-]", "[Na+].[Cl-]", "c1ccc(cc1)[NH3+]", "C(=O)([O-])[O-]"]
    for i, smi in enumerate(charge):
        entries.append({"id": f"charge_{i}", "category": "charge", "smiles": smi, "with_3d": False})

    radical = ["[CH3]", "[O][O]", "c1cc[c]cc1"]
    for i, smi in enumerate(radical):
        entries.append({"id": f"radical_{i}", "category": "radical", "smiles": smi, "with_3d": False})

    tetrahedral_stereo = ["[C@H](N)(C)C(=O)O", "[C@@H](N)(C)C(=O)O", "C[C@H](O)[C@@H](N)C", "C1C[C@H]2CC[C@@H]1C2"]
    for i, smi in enumerate(tetrahedral_stereo):
        entries.append({"id": f"tetrahedral_stereo_{i}", "category": "tetrahedral_stereo", "smiles": smi, "with_3d": False})

    ez_stereo = ["C/C=C/C", "C/C=C\\C", "CC(=C(\\C)Cl)/F", "C(/C=C/Cl)Br"]
    for i, smi in enumerate(ez_stereo):
        entries.append({"id": f"ez_stereo_{i}", "category": "ez_stereo", "smiles": smi, "with_3d": False})

    atom_maps = ["[CH3:1][OH:2]", "[CH:1](=[O:2])[NH2:3]", "[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1"]
    for i, smi in enumerate(atom_maps):
        entries.append({"id": f"atom_map_{i}", "category": "atom_map", "smiles": smi, "with_3d": False})

    coords_3d = ["CCO", "c1ccccc1", "CC(=O)O", "C1CCCCC1", "CC(N)C(=O)O"]
    for i, smi in enumerate(coords_3d):
        entries.append({"id": f"coords_3d_{i}", "category": "coords_3d", "smiles": smi, "with_3d": True})

    disconnected = ["CC(=O)O.[Na+]", "[Na+].[Cl-]", "c1ccccc1.c1ccncc1", "CC(=O)[O-].[NH4+]"]
    for i, smi in enumerate(disconnected):
        entries.append({"id": f"disconnected_{i}", "category": "disconnected", "smiles": smi, "with_3d": False})

    return [e for e in entries if rdkit_mol(e["smiles"]) is not None]


def build_fixture(entry, out_dir):
    mol = Chem.MolFromSmiles(entry["smiles"])
    rdDepictor.Compute2DCoords(mol)

    coords_3d = None
    if entry["with_3d"]:
        mol3d = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol3d, randomSeed=42) == 0:
            AllChem.MMFFOptimizeMolecule(mol3d)
            mol3d = Chem.RemoveHs(mol3d)
            conf = mol3d.GetConformer()
            coords_3d = [list(conf.GetAtomPosition(i)) for i in range(mol3d.GetNumAtoms())]
            mol = mol3d  # use the (possibly atom-reordered-by-Hs) 3D-embedded mol consistently

    block = Chem.MolToMrvBlock(mol)
    fixture_path = f"{out_dir}/{entry['id']}.mrv"
    with open(fixture_path, "w") as f:
        f.write(block)

    conf = mol.GetConformer()
    coords_2d = [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y] for i in range(mol.GetNumAtoms())]

    atoms = []
    for atom in mol.GetAtoms():
        atoms.append({
            "symbol": atom.GetSymbol(),
            "charge": atom.GetFormalCharge(),
            "isotope": atom.GetIsotope() or None,
            "atom_map": atom.GetAtomMapNum() or None,
            "radical_electrons": atom.GetNumRadicalElectrons(),
        })

    bonds = []
    for bond in mol.GetBonds():
        bonds.append({
            "begin": bond.GetBeginAtomIdx(),
            "end": bond.GetEndAtomIdx(),
            "order": str(bond.GetBondType()),
        })

    return {
        "id": entry["id"],
        "category": entry["category"],
        "file": f"{entry['id']}.mrv",
        "known_smiles": entry["smiles"],
        "rdkit_canonical_smiles": Chem.MolToSmiles(mol),
        "atom_count": mol.GetNumAtoms(),
        "bond_count": mol.GetNumBonds(),
        "atoms": atoms,
        "bonds": bonds,
        "coords_2d": coords_2d,
        "coords_3d": coords_3d,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--corpus", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--manifest-out", required=True)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    import os
    os.makedirs(args.out_dir, exist_ok=True)

    corpus = load_corpus(args.corpus, args.limit)
    entries = build_entries(corpus)

    fixtures = []
    failed = []
    for entry in entries:
        try:
            fixtures.append(build_fixture(entry, args.out_dir))
        except Exception as e:
            failed.append({"id": entry["id"], "smiles": entry["smiles"], "error": str(e)})

    manifest = {
        "rdkit_version": __import__("rdkit").__version__,
        "total_fixtures": len(fixtures),
        "total_failed_to_generate": len(failed),
        "failed": failed,
        "fixtures": fixtures,
    }
    with open(args.manifest_out, "w") as f:
        json.dump(manifest, f, indent=1)

    print(f"total_fixtures={len(fixtures)} failed_to_generate={len(failed)}", file=sys.stderr)


if __name__ == "__main__":
    main()
