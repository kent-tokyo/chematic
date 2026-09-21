#!/usr/bin/env python3
"""Measure TPSA parity on exposed public corpora without sealed-row reuse."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CORPORA = (
    ROOT / "scripts/chembl_accuracy_corpus_4999.smi",
    ROOT / "scripts/descriptor_census_corpus.smi",
)
DEFAULT_BINARY = ROOT / "target/debug/examples/descriptor_binding_dump"


def sha256(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def read_smiles(path: Path) -> list[str]:
    return [
        line.split()[0]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    ]


def evaluate(path: Path, binary: Path, rdkit_modules: dict[str, object]) -> dict[str, object]:
    smiles = read_smiles(path)
    run = subprocess.run(
        [str(binary)],
        input="\n".join(smiles) + "\n",
        text=True,
        capture_output=True,
        check=True,
    )
    records = [json.loads(line) for line in run.stdout.splitlines() if line.strip()]
    if len(records) != len(smiles):
        raise RuntimeError(f"{path}: implementation emitted {len(records)} rows for {len(smiles)} inputs")

    parse_failures = 0
    mismatches: list[dict[str, object]] = []
    Chem = rdkit_modules["Chem"]
    rdMolDescriptors = rdkit_modules["rdMolDescriptors"]
    for source, record in zip(smiles, records, strict=True):
        rd_mol = Chem.MolFromSmiles(source)
        if rd_mol is None or record.get("status") != "ok":
            parse_failures += 1
            continue
        expected = rdMolDescriptors.CalcTPSA(rd_mol, includeSandP=True)
        actual = record["descriptors"]["tpsa"]
        delta = abs(float(actual) - float(expected))
        if delta > 1e-6:
            mismatches.append(
                {"smiles": source, "chematic": actual, "rdkit": expected, "absolute_error": delta}
            )
    return {
        "path": str(path.relative_to(ROOT)),
        "sha256": sha256(path),
        "rows": len(smiles),
        "parse_failures": parse_failures,
        "strict_matches": len(smiles) - parse_failures - len(mismatches),
        "mismatches": len(mismatches),
        "mismatch_examples": mismatches[:20],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--corpus", type=Path, action="append")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    binary = args.binary.resolve()
    corpora = [path.resolve() for path in (args.corpus or DEFAULT_CORPORA)]
    if not binary.is_file():
        parser.error(f"descriptor binding dump not found: {binary}")

    import rdkit
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    results = [
        evaluate(path, binary, {"Chem": Chem, "rdMolDescriptors": rdMolDescriptors})
        for path in corpora
    ]
    total_rows = sum(int(result["rows"]) for result in results)
    strict_matches = sum(int(result["strict_matches"]) for result in results)
    report = {
        "schema_version": 1,
        "profile": "tpsa_exposed_public_corpus_parity_v1",
        "status": "development_regression",
        "scope": "Exposed public corpora; not sealed, unused-data evidence, or an adoption decision.",
        "oracle": {
            "engine": "RDKit",
            "version": rdkit.__version__,
            "operation": "rdMolDescriptors.CalcTPSA(includeSandP=True)",
        },
        "implementation": {
            "engine": "CheMatic Rust source",
            "operation": "descriptor_binding_dump tpsa",
            "binary_sha256": sha256(binary),
        },
        "strict_tolerance": 1e-6,
        "rows": total_rows,
        "strict_matches": strict_matches,
        "mismatches": total_rows - strict_matches,
        "corpora": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: report[key] for key in ("rows", "strict_matches", "mismatches")}))
    return 0 if report["mismatches"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
