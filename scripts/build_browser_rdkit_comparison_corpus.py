#!/usr/bin/env python3
"""Build the fixed, exposed candidate corpus for the browser RDKit.js comparison.

This is deliberately a performance corpus, not a sealed accuracy holdout. It
combines already-versioned inputs in a documented order, removes exact SMILES
duplicates, and writes both the selected rows and an input/output provenance
manifest. The hash assertions make a changed source fail closed instead of
silently changing the comparison population.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "validation" / "benchmark_corpora" / "rdkit-js-browser-candidates-v1.smi"
MANIFEST = OUTPUT.with_suffix(".manifest.json")
TARGET_ROWS = 12_000

SOURCES = (
    {
        "path": "scripts/descriptor_census_corpus.smi",
        "sha256": "d6f2ba3f128296f935007f0b0813aa97b6ebcc2457e014ddca2213ddd655276c",
        "description": "ChEMBL-derived descriptor census; CC BY-SA 3.0 provenance is documented in docs/rfcs/descriptor_census_rfc.md.",
    },
    {
        "path": "scripts/nci_first_5k_smiles_only.smi",
        "sha256": "1d9fa413f33c4d2e2f2da9742cdbd7ce8f39d33b0b9fd78ce012a0b7f24e031d",
        "description": "RDKit bundled NCI Diversity Set first_5K regression input; it is existing, exposed project data.",
    },
    {
        "path": "scripts/chembl_accuracy_corpus_4999.smi",
        "sha256": "1c47371dcbe37f4e0a141bf545b72bf238de2761fa3894fa251a552d84728d3e",
        "description": "Existing ChEMBL-derived accuracy corpus; CC BY-SA 3.0 provenance is recorded in validation/manifests/dataset_provenance.json.",
    },
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def main() -> int:
    selected: list[str] = []
    seen: set[str] = set()
    source_records: list[dict[str, object]] = []
    for source in SOURCES:
        path = ROOT / source["path"]
        actual = digest(path)
        if actual != source["sha256"]:
            raise SystemExit(
                f"source hash changed for {source['path']}: expected {source['sha256']}, got {actual}"
            )
        input_rows = rows(path)
        before = len(selected)
        for smiles in input_rows:
            if smiles in seen:
                continue
            seen.add(smiles)
            selected.append(smiles)
            if len(selected) == TARGET_ROWS:
                break
        source_records.append({
            **source,
            "input_rows": len(input_rows),
            "selected_unique_rows": len(selected) - before,
        })
        if len(selected) == TARGET_ROWS:
            break

    if len(selected) != TARGET_ROWS:
        raise SystemExit(f"only {len(selected)} unique rows available; expected {TARGET_ROWS}")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text("\n".join(selected) + "\n", encoding="utf-8")
    manifest = {
        "schema_version": 1,
        "id": "rdkit-js-browser-10k-v1",
        "purpose": "Fixed candidate input for browser compatibility qualification; not the final common-parse comparison population.",
        "status": "exposed_not_sealed",
        "not_for": "accuracy tuning, model selection, or a claim of independent generalization",
        "deduplication": "exact input SMILES, first occurrence wins in declared source order",
        "candidate_rows": TARGET_ROWS,
        "sources": source_records,
        "output": {
            "path": str(OUTPUT.relative_to(ROOT)),
            "rows": len(selected),
            "sha256": digest(OUTPUT),
        },
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest["output"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
