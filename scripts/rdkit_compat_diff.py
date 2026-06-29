#!/usr/bin/env python3
"""SMARTS + aromaticity differential validation: chematic vs RDKit.

Two case families over a molecule corpus:

1. **SMARTS substructure** — for each (molecule, query) compare chematic's
   ``GetSubstructMatches`` to RDKit's as *sets of frozensets* (order-invariant;
   both index heavy atoms in input-SMILES order so atom indices correspond).
   A query chematic cannot handle is recorded as ``unsupported`` (it raises;
   never silently ignored).

2. **Aromaticity** — order-invariant counts of aromatic atoms and aromatic
   bonds (atom ordering differs post-perception, so compare counts not indices).

Usage:
    python scripts/rdkit_compat_diff.py [SMILES.csv] [--limit N]

Writes validation/results/rdkit_compat_diff.jsonl (one row per divergence).
"""
import json
import os
import sys

from chematic import rdkit_compat as Chem

try:
    from rdkit import Chem as RDChem
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
except ImportError:
    sys.exit("RDKit is required for rdkit_compat_diff.py (pip install rdkit)")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "validation", "results", "rdkit_compat_diff.jsonl")

# Curated query set (molecules come from the corpus, queries do not).
SMARTS = [
    "[OH]", "c", "[#7]", "C=O", "[NX3;H2,H1;!$(NC=O)]", "[r5]", "[r6]",
    "c1ccccc1", "[CX4]", "[!#6;!#1]", "C(=O)[OH]", "[nH]", "[#6]=[#6]",
    "[F,Cl,Br,I]", "[OX2H]", "[#16]",
]


def match_set(matches):
    """A set of frozensets of atom indices (order-invariant)."""
    return frozenset(frozenset(m) for m in matches)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-") \
        else os.path.expanduser("~/Downloads/SMILES.csv")
    limit = None
    if "--limit" in sys.argv:
        limit = int(sys.argv[sys.argv.index("--limit") + 1])

    smis = [l.strip() for l in open(path) if l.strip()]
    if limit:
        smis = smis[:limit]

    # Pre-parse RDKit queries once.
    rd_queries = {s: RDChem.MolFromSmarts(s) for s in SMARTS}

    n_mol = 0
    smarts_ok = smarts_mismatch = smarts_unsupported = 0
    arom_atom_ok = arom_atom_bad = 0
    arom_bond_ok = arom_bond_bad = 0
    rows = []

    for smi in smis:
        rm = RDChem.MolFromSmiles(smi)
        if rm is None:
            continue
        cm = Chem.MolFromSmiles(smi)
        if cm is None:
            continue
        n_mol += 1

        # --- SMARTS ---
        for s in SMARTS:
            rq = rd_queries[s]
            if rq is None:
                continue
            rd_set = match_set(rm.GetSubstructMatches(rq, uniquify=True))
            try:
                cm_set = match_set(cm.GetSubstructMatches(s, uniquify=True))
            except Exception as e:
                smarts_unsupported += 1
                rows.append({"case": "smarts", "smiles": smi, "query": s,
                             "status": "unsupported", "note": str(e)})
                continue
            if cm_set == rd_set:
                smarts_ok += 1
            else:
                smarts_mismatch += 1
                rows.append({"case": "smarts", "smiles": smi, "query": s,
                             "status": "count_differs",
                             "chematic_n": len(cm_set), "rdkit_n": len(rd_set)})

        # --- aromaticity (counts) ---
        cm_arom_atoms = sum(1 for a in cm._mol.atom_table if a[3])
        rd_arom_atoms = sum(1 for a in rm.GetAtoms() if a.GetIsAromatic())
        if cm_arom_atoms == rd_arom_atoms:
            arom_atom_ok += 1
        else:
            arom_atom_bad += 1
            rows.append({"case": "aromatic_atoms", "smiles": smi,
                         "status": "count_differs",
                         "chematic": cm_arom_atoms, "rdkit": rd_arom_atoms})

        cm_arom_bonds = sum(1 for b in cm._mol.bond_table if b[3])
        rd_arom_bonds = sum(1 for b in rm.GetBonds() if b.GetIsAromatic())
        if cm_arom_bonds == rd_arom_bonds:
            arom_bond_ok += 1
        else:
            arom_bond_bad += 1
            rows.append({"case": "aromatic_bonds", "smiles": smi,
                         "status": "count_differs",
                         "chematic": cm_arom_bonds, "rdkit": rd_arom_bonds})

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")

    smarts_total = smarts_ok + smarts_mismatch + smarts_unsupported
    print(f"corpus (parsed by both):  {n_mol}")
    print(f"SMARTS comparisons:       {smarts_total}")
    print(f"  match-set agreement:    {smarts_ok}/{smarts_total} "
          f"({100 * smarts_ok / max(smarts_total, 1):.2f}%)")
    print(f"  count differs:          {smarts_mismatch}")
    print(f"  unsupported (raised):   {smarts_unsupported}")
    print(f"aromatic-atom-count agreement: {arom_atom_ok}/{n_mol} "
          f"({100 * arom_atom_ok / max(n_mol, 1):.2f}%)  bad={arom_atom_bad}")
    print(f"aromatic-bond-count agreement: {arom_bond_ok}/{n_mol} "
          f"({100 * arom_bond_ok / max(n_mol, 1):.2f}%)  bad={arom_bond_bad}")
    print(f"\nwrote {len(rows)} divergence rows to {os.path.relpath(OUT, ROOT)}")


if __name__ == "__main__":
    main()
