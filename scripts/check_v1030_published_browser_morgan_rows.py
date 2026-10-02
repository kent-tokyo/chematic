#!/usr/bin/env python3
"""Check direct browser bit parity on pinned published v1.0.30 artifacts.

This checks the exposed rows and typed refusal, not a general RDKit parity
claim. The browser runner compares every supported fingerprint bit in-page.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

if __package__:
    from .check_published_browser_paired import EXPECTED
else:
    from check_published_browser_paired import EXPECTED

ROOT = Path(__file__).resolve().parents[1]
CORPORA = {
    250: ROOT / "scripts/descriptor_census_corpus.smi",
    10_000: ROOT / "validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi",
}
PACKAGE_JSON_SHA256 = "6a18eb0fa882b4a2772751c9484db8200deca21764f0a9418d8a94196d6b5b90"
NPM_TARBALL_SHA256 = "fd427a161c11b14e6cca51e78e0c35ff05bbc52faf0539a2ddc2b7090e534201"
RDKIT_PACKAGE_JSON_SHA256 = "bb938eb6b28bf39baa9778f56aefcdea746d28060fe2c542209f877a71468fed"
UNSUPPORTED_INDEX = 8341


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(path: Path, operation: str, rows: int, engine: str = "chromium") -> dict:
    if (operation not in {"direct", "prepared"} or rows not in CORPORA
            or engine not in {"chromium", "firefox", "webkit"}):
        raise ValueError("unrecognized browser Morgan lane")
    document = json.loads(path.read_bytes())
    config = document["configuration"]
    artifacts = document["artifacts"]
    result = document["result"]
    if (document["schema_version"] != 1
            or document["gate"] != "browser-rdkit-ecfp4-bit-parity"
            or config["engine"] != engine
            or config["schematic_operation"] != operation
            or config["rows"] != rows
            or config["allow_unsupported"] != (rows == 10_000)
            or document["corpus"]["rows"] != rows
            or document["corpus"]["sha256"] != sha256(CORPORA[rows])):
        raise ValueError("browser corpus or protocol mismatch")
    for key, expected in EXPECTED.items():
        if artifacts[key]["sha256"] != expected:
            raise ValueError(f"browser artifact hash mismatch: {key}")
    package = artifacts["schematic_package"]
    if (artifacts["schematic_wasm"]["kind"] != "published"
            or package["name"] != "@kent-tokyo/chematic"
            or package["version"] != "1.0.30"
            or package["package_json_sha256"] != PACKAGE_JSON_SHA256
            or package["tarball_sha256"] != NPM_TARBALL_SHA256
            or artifacts["rdkit_package"]["version"] != "2026.03.6"
            or artifacts["rdkit_package"]["package_json_sha256"] != RDKIT_PACKAGE_JSON_SHA256
            or result["rdkit_version"] != "2026.03.6"):
        raise ValueError("browser package identity mismatch")
    unsupported = 1 if rows == 10_000 else 0
    if (result["compared_rows"] != rows
            or result["supported_rows"] != rows - unsupported
            or result["exact_matches"] != rows - unsupported
            or result["unsupported_rows"] != unsupported
            or result["first_mismatch"] is not None
            or result["mismatch_sample"]):
        raise ValueError("browser row-level Morgan parity failed")
    sample = result["unsupported_sample"]
    if unsupported:
        corpus_rows = [line.strip() for line in CORPORA[rows].read_text().splitlines() if line.strip()]
        if (len(sample) != 1 or sample[0]["index"] != UNSUPPORTED_INDEX
                or sample[0]["smiles"] != corpus_rows[UNSUPPORTED_INDEX]
                or "unsupported RDKit coordination sanitization" not in sample[0]["error"]):
            raise ValueError("browser typed-refusal identity changed")
    elif sample:
        raise ValueError("unexpected browser typed refusal")
    return {"rows": rows, "exact": result["exact_matches"],
            "typed_unsupported": unsupported, "operation": operation}


def main() -> int:
    results = []
    for rows, prefix in ((250, "row"), (10_000, "10k")):
        for operation in ("direct", "prepared"):
            path = ROOT / "benchmarks" / f"2026-10-03-v1030-published-chromium-morgan-{prefix}-{operation}.json"
            results.append(check(path, operation, rows))
    print(json.dumps({"browser_published_morgan_rows": results}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
