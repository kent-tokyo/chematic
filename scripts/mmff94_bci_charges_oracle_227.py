#!/usr/bin/env python3
"""Live-RDKit-oracle dump of MMFF94 partial charges (issue #227 Phase 2, BCI
investigation). Measurement-only -- no chematic-ff code is touched by this
script.

Question this answers: does RDKit's real `computeMMFFCharges` use the same
per-bond `getMMFFBondType(bond)` RDKit's own bond-stretch/angle/torsion code
uses (the Kekulized/reperceived view chematic-ff's Phase 1 fix already
threads through those four term kinds), or something else? Source-level
answer already confirmed by direct read of the pinned RDKit commit
(`e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f`, `Code/GraphMol/
ForceFieldHelpers/MMFF/AtomTyper.cpp:3071-3488`): `computeMMFFCharges` calls
`this->getMMFFBondType(bond)` at line 3472 -- textually the SAME method
`getMMFFBondStretchParams` calls at line 3500, both keyed off the identical
`bond` object's `getBondType()` (i.e. the sanitized/Kekulized RDKit `mol`
built once per `MMFFMolProperties` construction). This script supplies the
empirical half of the falsification: does chematic's current (pre-Phase-2)
`mmff94_charges_numeric` output already match this oracle, or not?

No embedding/geometry needed -- MMFF charge/type assignment is purely
topological, same precedent as `scripts/mmff94_torsion_oracle_validate_227.py`
and `scripts/mmff94_stbn_oracle_validate_227.py` (both call
`MMFFGetMoleculeProperties` directly on the implicit-H `Chem.MolFromSmiles`
result, no `AddHs`/no conformer).

Run:
    .venv/bin/python scripts/mmff94_bci_charges_oracle_227.py \\
        > validation/results/mmff94_bci_charges_oracle_227.jsonl
"""

import json
import sys

from rdkit import Chem
from rdkit.Chem.rdForceFieldHelpers import MMFFGetMoleculeProperties
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

MANIFESTS = [
    ("A", "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_a.json"),
    ("B", "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json"),
]


def load_manifest(path):
    with open(path) as f:
        data = json.load(f)
    return data["molecules"]


def process(tier, name, smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"tier": tier, "name": name, "smiles": smiles, "status": "parse_failure"}

    try:
        props = MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
    except Exception as e:  # noqa: BLE001 -- typed failure, recorded not swallowed
        return {
            "tier": tier,
            "name": name,
            "smiles": smiles,
            "status": "mmff_properties_exception",
            "reason": str(e),
        }

    if props is None:
        return {"tier": tier, "name": name, "smiles": smiles, "status": "mmff_properties_none"}

    n_heavy = mol.GetNumAtoms()
    charges = []
    for idx in range(n_heavy):
        atom = mol.GetAtomWithIdx(idx)
        charges.append(
            {
                "index": idx,
                "element": atom.GetSymbol(),
                "rdkit_mmff_type": props.GetMMFFAtomType(idx),
                "rdkit_partial_charge": props.GetMMFFPartialCharge(idx),
            }
        )

    return {
        "tier": tier,
        "name": name,
        "smiles": smiles,
        "status": "ok",
        "n_heavy": n_heavy,
        "charges": charges,
    }


def main():
    rows = []
    for tier, path in MANIFESTS:
        for m in load_manifest(path):
            rows.append((tier, m["name"], m["smiles"]))

    n_ok = 0
    n_fail = 0
    for tier, name, smiles in rows:
        row = process(tier, name, smiles)
        print(json.dumps(row))
        if row["status"] == "ok":
            n_ok += 1
        else:
            n_fail += 1

    print(f"total={len(rows)} ok={n_ok} fail={n_fail}", file=sys.stderr)


if __name__ == "__main__":
    main()
