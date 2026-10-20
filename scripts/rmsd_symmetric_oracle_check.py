#!/usr/bin/env python3
"""Independent RDKit oracle check for `chematic_3d::conformer::rmsd_symmetric`
(crates/chematic-3d/src/conformer.rs).

Reads the JSONL dump produced by
`crates/chematic-3d/examples/rmsd_symmetric_oracle_dump.rs` (one row per
molecule: SMILES + two heavy-atom-only conformers + chematic's own
`rmsd_symmetric` result) and independently recomputes each pair's
symmetry-aware RMSD via RDKit's `rdMolAlign.GetBestRMS`, the function
chematic's `rmsd_symmetric` was ported from.

The comparison uses RDKit's default terminal-group symmetrization. Acetate
and nitro exercise interchange of formally different terminal oxygens;
isotopically labelled isobutane checks that distinct labels stay distinct.

Run:
    cargo run --release -p chematic-3d --example rmsd_symmetric_oracle_dump \
      > /tmp/rmsd_symmetric_oracle_dump.jsonl
    .venv/bin/python scripts/rmsd_symmetric_oracle_check.py \
      /tmp/rmsd_symmetric_oracle_dump.jsonl
"""

import json
import sys

from rdkit import Chem
from rdkit.Chem import rdMolAlign
from rdkit.Geometry import Point3D

AGREE_TOL = 1e-3  # Å; two independent Kabsch/automorphism-search implementations


def mol_with_conformer(smiles, coords):
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        raise ValueError(f"RDKit failed to parse {smiles!r}")
    if mol.GetNumAtoms() != len(coords):
        raise ValueError(
            f"{smiles!r}: RDKit atom count {mol.GetNumAtoms()} != "
            f"dump coord count {len(coords)} (SMILES atom-order mismatch?)"
        )
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(i, Point3D(x, y, z))
    mol.AddConformer(conf, assignId=True)
    return mol


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/rmsd_symmetric_oracle_dump.jsonl"
    rows = [json.loads(line) for line in open(path)]

    print(f"{'name':<14} {'chematic':>12} {'rdkit':>12} {'diff':>10}  status")
    mismatches = []
    for row in rows:
        name = row["name"]
        probe = mol_with_conformer(row["smiles"], row["conformer_a"])
        ref = mol_with_conformer(row["smiles"], row["conformer_b"])
        rdkit_rmsd = rdMolAlign.GetBestRMS(probe, ref)
        chematic_rmsd = row["chematic_rmsd_symmetric"]
        diff = abs(chematic_rmsd - rdkit_rmsd)

        if diff <= AGREE_TOL:
            status = "OK"
        else:
            status = "MISMATCH"
            mismatches.append(name)

        print(f"{name:<14} {chematic_rmsd:>12.6f} {rdkit_rmsd:>12.6f} {diff:>10.6f}  {status}")

    if mismatches:
        print(f"\n{len(mismatches)} unexplained mismatch(es): {mismatches}", file=sys.stderr)
        sys.exit(1)
    print(f"\nAll {len(rows)} cases agree with RDKit within {AGREE_TOL} Å.")


if __name__ == "__main__":
    main()
