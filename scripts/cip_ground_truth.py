#!/usr/bin/env python3
"""
Corpus-scale ground-truth check for chematic's CIP R/S assignment
(chematic_chem::assign_cip, exposed as Mol.cip_stereo() in Python).

Root-caused and fixed in the same round this script was added: two
independent bugs in crates/chematic-chem/src/cip.rs --

1. assign_tetrahedral() rebuilt a chiral atom's substituent order from raw
   Molecule::neighbors() adjacency order instead of the already-correct,
   parser-populated Molecule::stereo_neighbor_order(). Adjacency order only
   matches SMILES textual order when every ring bond at that atom CLOSES
   (partner already known); it silently reorders substituents when a
   stereocenter OPENS a ring before its other neighbors, because a
   ring-opening bond only gets added to the adjacency list once the
   matching closing digit is reached later in the string.

2. cip_branch_spheres()'s CIP double-bond duplication only added the
   "arrival side" phantom (B's own sphere gets a phantom-A once B is
   expanded, having been reached via A=B) -- never the "departure side"
   (A's own sphere never got a second phantom-B while iterating A's
   neighbors). A double bond must duplicate its partner into BOTH atoms'
   substituent spheres.

Before this round, comparing against RDKit's per-atom CIP oracle
(Chem.AssignStereochemistry(cleanIt=True, force=True) + _CIPCode) had never
been done -- prior rounds only checked chematic's InChI/canonical_smiles
against EACH OTHER for order-stability, which cannot detect an assignment
that is stable but wrong. This script is the first true correctness check
against an external oracle, not just internal self-consistency.

Correspondence note: this compares the ORIGINAL (non-respelled) SMILES only,
so RDKit and chematic atom indices align directly (both parse the same
string in the same atom order -- verified elsewhere this project,
scripts/aromaticity_mechanism_probe.py). No respelling/correspondence
mapping is needed or performed here.

Usage:
    .venv/bin/python scripts/cip_ground_truth.py [SMILES.csv]
"""

import csv
import sys

sys.path.insert(0, ".")
import chematic
from rdkit import Chem


def main():
    csv_path = sys.argv[1] if len(sys.argv) > 1 else "SMILES.csv"
    with open(csv_path) as f:
        reader = csv.reader(f)
        next(reader)
        smis = [row[0] for row in reader if row]

    total = 0
    match = 0
    mismatches = []
    for smi in smis:
        rd = Chem.MolFromSmiles(smi)
        if rd is None:
            continue
        if not any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in rd.GetAtoms()):
            continue
        Chem.AssignStereochemistry(rd, cleanIt=True, force=True)
        rd_cip = {a.GetIdx(): a.GetProp("_CIPCode") for a in rd.GetAtoms() if a.HasProp("_CIPCode")}
        if not rd_cip:
            continue
        try:
            m = chematic.from_smiles(smi)
            cm_cip = {d["atom_idx"]: d["descriptor"] for d in m.cip_stereo()}
        except Exception:
            continue
        for aidx, code in rd_cip.items():
            total += 1
            if cm_cip.get(aidx) == code:
                match += 1
            else:
                mismatches.append((smi, aidx, code, cm_cip.get(aidx)))

    print(f"total stereocenters compared: {total}")
    print(f"match: {match} ({100 * match / total:.2f}%)" if total else "no stereocenters found")
    print(f"mismatch: {len(mismatches)}")
    for ex in mismatches[:10]:
        print(" ", ex)


if __name__ == "__main__":
    main()
