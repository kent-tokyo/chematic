#!/usr/bin/env python3
"""Issue #227 Phase 2: paired Torsion Fingerprint Deviation (TFD) between one
chematic arm and one RDKit arm, over molecules that succeeded on both sides.

chematic has no TFD implementation of its own (grepped, confirmed absent) --
per the Phase 2 directive's explicit instruction not to hand-roll one inside
a measurement PR, this uses RDKit's own `TorsionFingerprints.
GetTFDBetweenMolecules` directly, the same "reuse the real library, don't
reimplement" principle `rmsd_symmetric_oracle_check.py` already established
for RMSD (there via `rdMolAlign.GetBestRMS`).

Reads the SAME two dump JSONL files the Rust paired-RMSD addition to
`pipeline_v2_vs_rdkit_common_scorer.rs` reads (one row per (molecule, arm),
`coords` = heavy-atom-only, same manifest SMILES/atom-order as chematic's
own parse -- atom mapping already verified 265/265 in the Wave 1 aggregate).
Builds two RDKit `Mol` objects per joined molecule (one per engine), each
with a single conformer set from that engine's saved coordinates (heavy
atoms only, no explicit Hs added -- matches how both engines' coordinates
were saved), then calls `GetTFDBetweenMolecules(mol_chematic, mol_rdkit)`.

Usage:
    .venv/bin/python scripts/pipeline_v2_vs_rdkit_tfd_227.py \\
        <chematic_rows.jsonl> <rdkit_rows.jsonl> <chematic_arm> <rdkit_arm> \\
        > <output>.jsonl
"""

import json
import sys

from rdkit import Chem
from rdkit.Chem import TorsionFingerprints
from rdkit.Geometry import Point3D
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

MANIFESTS = [
    ("A", "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_a.json"),
    ("B", "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json"),
]


def load_manifest_smiles():
    out = {}
    for tier, path in MANIFESTS:
        with open(path) as f:
            data = json.load(f)
        out[tier] = {m["name"]: m["smiles"] for m in data["molecules"]}
    return out


def load_jsonl(path):
    with open(path) as f:
        return [json.loads(line) for line in f if line.strip()]


def coords_by_key(rows, arm):
    out = {}
    for r in rows:
        if r.get("arm") == arm and r.get("status") == "success" and r.get("coords") is not None:
            out[(r["tier"], r["name"])] = r["coords"]
    return out


def mol_with_conformer(smiles, coords):
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        return None, "parse_failure"
    if mol.GetNumAtoms() != len(coords):
        return None, "coords_count_mismatch"
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(i, Point3D(x, y, z))
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return mol, None


def main():
    chematic_path, rdkit_path, chematic_arm, rdkit_arm = sys.argv[1:5]
    smiles_by_name = load_manifest_smiles()
    chematic_rows = load_jsonl(chematic_path)
    rdkit_rows = load_jsonl(rdkit_path)

    ch_coords = coords_by_key(chematic_rows, chematic_arm)
    rd_coords = coords_by_key(rdkit_rows, rdkit_arm)

    n_ok = 0
    n_fail = 0
    for tier, name in sorted(ch_coords):
        if (tier, name) not in rd_coords:
            continue
        smiles = smiles_by_name.get(tier, {}).get(name)
        if smiles is None:
            print(json.dumps({"tier": tier, "name": name, "status": "integrity_error",
                               "reason": "smiles_not_found_in_manifest"}))
            n_fail += 1
            continue

        mol_ch, err_ch = mol_with_conformer(smiles, ch_coords[(tier, name)])
        mol_rd, err_rd = mol_with_conformer(smiles, rd_coords[(tier, name)])
        if err_ch or err_rd:
            print(json.dumps({"tier": tier, "name": name, "status": "integrity_error",
                               "reason": f"chematic:{err_ch} rdkit:{err_rd}"}))
            n_fail += 1
            continue

        try:
            tfd = TorsionFingerprints.GetTFDBetweenMolecules(mol_ch, mol_rd)
        except Exception as e:  # noqa: BLE001 -- typed, recorded, not swallowed
            print(json.dumps({"tier": tier, "name": name, "status": "tfd_exception",
                               "reason": str(e)}))
            n_fail += 1
            continue

        print(json.dumps({
            "tier": tier,
            "name": name,
            "chematic_arm": chematic_arm,
            "rdkit_arm": rdkit_arm,
            "status": "paired_tfd",
            "tfd": tfd,
        }))
        n_ok += 1

    print(f"tfd_ok={n_ok} tfd_fail={n_fail}", file=sys.stderr)


if __name__ == "__main__":
    main()
