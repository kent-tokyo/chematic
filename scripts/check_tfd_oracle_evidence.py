#!/usr/bin/env python3
"""Validate the checked-in RDKit TFD comparison boundary.

The current result explicitly records rigid, zero-rotatable-torsion probes as
``not_applicable``. All other rows must contain finite TFD values and unique
corpus identities.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import TorsionFingerprints


ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "validation" / "results" / "pipeline_v2_vs_rdkit_tfd_current-v1.0.10.jsonl"
MANIFESTS = {
    "A": ROOT / "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_a.json",
    "B": ROOT / "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json",
}


def main() -> int:
    smiles: dict[tuple[str, str], str] = {}
    for tier, path in MANIFESTS.items():
        document = json.loads(path.read_text(encoding="utf-8"))
        for molecule in document["molecules"]:
            smiles[(tier, molecule["name"])] = molecule["smiles"]

    rows = [json.loads(line) for line in RESULT.read_text(encoding="utf-8").splitlines() if line.strip()]
    keys = [(row.get("tier"), row.get("name")) for row in rows]
    errors: list[str] = []
    if len(rows) != 250:
        errors.append(f"expected 250 corpus rows, got {len(rows)}")
    if len(set(keys)) != len(keys):
        errors.append("duplicate tier/name rows")

    paired = [row for row in rows if row.get("status") == "paired_tfd"]
    not_applicable = [row for row in rows if row.get("status") == "not_applicable"]
    other = [row for row in rows if row.get("status") not in {"paired_tfd", "not_applicable"}]
    if len(paired) != 201:
        errors.append(f"expected 201 finite TFD rows, got {len(paired)}")
    if len(not_applicable) != 49:
        errors.append(f"expected 49 rigid not-applicable rows, got {len(not_applicable)}")
    if other:
        errors.append(f"unexpected statuses: {sorted({row.get('status') for row in other})}")

    for row in paired:
        value = row.get("tfd")
        if not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0.0:
            errors.append(f"non-finite or negative TFD: {row.get('tier')}/{row.get('name')}")

    for row in not_applicable:
        key = (row.get("tier"), row.get("name"))
        smi = smiles.get(key)
        if smi is None:
            errors.append(f"not-applicable row is not in the pinned manifest: {key}")
            continue
        mol = Chem.MolFromSmiles(smi, sanitize=True)
        if mol is None:
            errors.append(f"not-applicable SMILES does not parse: {key}")
            continue
        torsions, _ = TorsionFingerprints.CalculateTorsionLists(mol)
        if torsions:
            errors.append(f"not-applicable row has rotatable torsions: {key}")

    if errors:
        print("TFD oracle evidence failures:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(
        "TFD oracle evidence OK: 250 unique corpus rows, "
        "201 finite paired values, 49 explicitly rigid/not-applicable rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
