#!/usr/bin/env python3
"""Measure the explicit Phase 1 fragment-selection profile on held-out fixtures."""
from __future__ import annotations
import argparse, json, subprocess
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fixtures", type=Path); p.add_argument("--json", type=Path, required=True)
    a = p.parse_args()
    rows = [json.loads(line) for line in a.fixtures.read_text().splitlines() if line.strip()]
    text = "\n".join(row["input"] for row in rows) + "\n"
    result = subprocess.run(["cargo", "run", "-p", "chematic-cli", "--release", "--", "batch-standardize", "--largest-fragment-only"], cwd=ROOT, input=text, text=True, capture_output=True, check=True)
    batch = json.loads(result.stdout)
    from rdkit import Chem
    outcomes = []
    for fixture, record in zip(rows, batch["records"]):
        standardization = record.get("standardization")
        actual = standardization.get("output_smiles") if standardization else None
        expected = fixture.get("expected_kept_fragment")
        def canonical(value):
            mol = Chem.MolFromSmiles(value) if value else None
            return Chem.MolToSmiles(mol, canonical=True) if mol else None
        outcomes.append({"id": fixture["id"], "input": fixture["input"], "expected_kept_fragment": expected, "output_smiles": actual, "expected_canonical_rdkit": canonical(expected), "output_canonical_rdkit": canonical(actual), "match": canonical(expected) == canonical(actual), "status": standardization.get("status") if standardization else None, "error": record.get("error")})
    matched = sum(row["match"] for row in outcomes)
    report = {"schema_version": 1, "target_version": "1.0.8", "fixture": str(a.fixtures), "operation": "standardize_largest_fragment_only", "comparison_boundary": "10 committed Phase 1 standardization holdouts; v1.0.8 Rust CLI explicit --largest-fragment-only profile, with RDKit canonical identity used only to compare expected and returned fragment", "rows": len(outcomes), "valid": batch["valid_count"], "errors": batch["error_count"], "matches": matched, "mismatches": len(outcomes) - matched, "match_pct": 100.0 * matched / len(outcomes) if outcomes else None, "outcomes": outcomes}
    a.json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8"); print(json.dumps(report, indent=2)); return 0
if __name__ == "__main__": raise SystemExit(main())
