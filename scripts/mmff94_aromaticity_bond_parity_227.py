#!/usr/bin/env python3
"""Issue #227 Phase 0.2: RDKit oracle for MMFF94-specific aromaticity, atom
AND bond level, on a 12-molecule fixture matrix chosen to stress fused-ring
aromatic-bond promotion (the exact class of bug Phase 0.2 fixes in
`compute_mmff94_aromatic_view`: a ring's bonds must only be promoted to
`BondOrder::Aromatic` when that specific ring's own Huckel check passed --
never reconstructed by checking whether every atom in the ring happens to be
aromatic via *other*, unrelated accepted rings).

Ground truth source: `Chem.SetAromaticity(mol, Chem.AROMATICITY_MMFF94)` on a
freshly-Kekulized copy -- this calls RDKit's real `setMMFFAromaticity`
(`Code/GraphMol/Aromaticity.cpp`) directly, independent of
`MMFFGetMoleculeProperties`'s own internal usage of it. Verified on caffeine
against the already-established, independently-confirmed pattern (carbonyl
ring atoms non-aromatic, imidazole ring atoms aromatic, fusion bond
aromatic-via-imidazole-only) before trusting it for the other 11 fixtures.

Run: .venv/bin/python scripts/mmff94_aromaticity_bond_parity_227.py \
  > validation/results/mmff94_aromaticity_bond_parity_227_oracle.json
"""

import json
import sys

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

FIXTURES = [
    ("benzene", "c1ccccc1"),
    ("naphthalene", "c1ccc2ccccc2c1"),
    ("anthracene", "c1ccc2cc3ccccc3cc2c1"),
    ("indole", "c1ccc2[nH]ccc2c1"),
    ("quinoline", "c1ccc2ncccc2c1"),
    ("purine", "c1nc2[nH]cnc2cn1"),
    ("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
    ("azulene", "c1ccc2cccccc12"),
    # Fused aromatic + non-aromatic ring: benzo ring stays aromatic,
    # fused saturated cyclohexane ring must NOT be promoted even though
    # both of its fusion atoms are independently aromatic via the benzo
    # ring -- the exact fused-ring hazard Phase 0.2 targets.
    ("tetralin_fused_nonaromatic", "c1ccc2c(c1)CCCC2"),
    # Exocyclic-carbonyl-bearing fused system, distinct skeleton from
    # caffeine: the pyridone-like ring's own exocyclic C=O blocks only
    # that ring's aromaticity; the fused benzo ring is unaffected.
    ("carbostyril_exocyclic_carbonyl", "O=c1ccc2ccccc2[nH]1"),
    # Spiro system: aromatic benzo ring, non-aromatic indane 5-ring, and a
    # non-aromatic cyclohexane ring joined at a single spiro atom.
    ("spiro_indane_cyclohexane", "c1ccc2c(c1)CC3(CCCCC3)C2"),
    # Bridged polycycle, no aromaticity at all -- a negative control for
    # SSSR/ring-perception edge cases bridged systems are known to stress.
    ("norbornane_bridged", "C1CC2CCC1C2"),
]


def process(name, smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"name": name, "smiles": smiles, "status": "parse_failure"}

    n_heavy = mol.GetNumAtoms()

    # MMFF-specific atom+bond aromaticity ground truth.
    amol = Chem.Mol(mol)
    Chem.Kekulize(amol, clearAromaticFlags=True)
    Chem.SetAromaticity(amol, Chem.AROMATICITY_MMFF94)
    atom_aromatic = {a.GetIdx(): a.GetIsAromatic() for a in amol.GetAtoms()}
    bond_aromatic = {}
    for b in amol.GetBonds():
        key = f"{min(b.GetBeginAtomIdx(), b.GetEndAtomIdx())}-{max(b.GetBeginAtomIdx(), b.GetEndAtomIdx())}"
        bond_aromatic[key] = b.GetIsAromatic()

    # Ring-level accept/reject, derived from the bond flags above: a ring
    # is MMFF-accepted iff every one of its own bonds is aromatic.
    ri = mol.GetRingInfo()
    ring_accept = []
    for ring_bonds in ri.BondRings():
        accepted = all(
            bond_aromatic.get(
                f"{min(b.GetBeginAtomIdx(), b.GetEndAtomIdx())}-{max(b.GetBeginAtomIdx(), b.GetEndAtomIdx())}",
                False,
            )
            for b in (amol.GetBondWithIdx(i) for i in ring_bonds)
        )
        ring_accept.append(
            {
                "atoms": sorted(
                    {a for i in ring_bonds for a in (amol.GetBondWithIdx(i).GetBeginAtomIdx(), amol.GetBondWithIdx(i).GetEndAtomIdx())}
                ),
                "accepted": accepted,
            }
        )

    # Numeric atom types + MMFF bond type via the standard embed+properties
    # path (same pattern as mmff94_rdkit_type_oracle.py).
    molH = Chem.AddHs(mol)
    seed = 20260801
    embed_ret = AllChem.EmbedMolecule(molH, randomSeed=seed, useRandomCoords=True)
    if embed_ret != 0:
        embed_ret = AllChem.EmbedMolecule(molH, randomSeed=seed + 1, useRandomCoords=True)
    props = AllChem.MMFFGetMoleculeProperties(molH) if embed_ret == 0 else None

    atom_types = {}
    if props is not None:
        for atom in molH.GetAtoms():
            if atom.GetIdx() < n_heavy:
                atom_types[atom.GetIdx()] = props.GetMMFFAtomType(atom.GetIdx())

    bond_types = {}
    if props is not None:
        for b in mol.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            key = f"{min(i, j)}-{max(i, j)}"
            try:
                bond_types[key] = props.GetMMFFBondType(mol.GetBondBetweenAtoms(i, j).GetIdx())
            except Exception:
                bond_types[key] = None

    return {
        "name": name,
        "smiles": smiles,
        "status": "ok",
        "n_heavy_atoms": n_heavy,
        "element_sequence": [a.GetSymbol() for a in mol.GetAtoms()],
        "atom_aromatic": atom_aromatic,
        "bond_aromatic": bond_aromatic,
        "ring_accept": ring_accept,
        "atom_types": atom_types,
        "bond_types": bond_types,
    }


def main():
    results = [process(name, smi) for name, smi in FIXTURES]
    json.dump(results, sys.stdout, indent=2, sort_keys=True)
    print()


if __name__ == "__main__":
    main()
