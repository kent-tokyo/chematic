#!/usr/bin/env python3
"""Compare v1.0.8 Rust ``rdkit_torsion`` with RDKit's hashed torsion bits."""
from __future__ import annotations
import argparse, hashlib, json, subprocess
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""): h.update(block)
    return h.hexdigest()
def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("corpus", type=Path); p.add_argument("--json", type=Path, required=True)
    a = p.parse_args()
    result = subprocess.run(["cargo", "run", "-p", "chematic-cli", "--release", "--", "batch-fingerprints", "--input", str(a.corpus), "--algorithm", "rdkit_torsion"], cwd=ROOT, check=True, capture_output=True, text=True)
    batch = json.loads(result.stdout); matches = mismatches = 0; cf = rf = 0; examples = []
    from rdkit import Chem
    from rdkit.Chem import AllChem
    for row in batch["records"]:
        if row.get("error") is not None: cf += 1; continue
        mol = Chem.MolFromSmiles(row["input_smiles"])
        if mol is None: rf += 1; continue
        expected = sorted(AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect(mol, nBits=2048).GetOnBits())
        actual = row["fingerprint"]["set_bits"]
        if actual == expected: matches += 1
        else:
            mismatches += 1
            if len(examples) < 25: examples.append({"smiles": row["input_smiles"], "chematic_set_bits": actual, "rdkit_set_bits": expected})
    report = {"schema_version": 1, "target_version": "1.0.8", "corpus": str(a.corpus), "corpus_sha256": digest(a.corpus), "operation": "rdkit_torsion", "comparison_boundary": "same SMILES inputs; v1.0.8 Rust CLI versus RDKit AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect", "configuration": {"nBits": 2048, "torsionAtomCount": 4}, "rows": len(batch["records"]), "chematic_valid": batch["valid_count"], "chematic_failures": cf, "rdkit_failures": rf, "exact_matches": matches, "exact_mismatches": mismatches, "exact_match_pct_of_common_success": 100.0 * matches / (matches + mismatches) if matches + mismatches else None, "rdkit_version": __import__("rdkit").rdBase.rdkitVersion, "mismatch_examples": examples}
    a.json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8"); print(json.dumps(report, indent=2)); return 0
if __name__ == "__main__": raise SystemExit(main())
