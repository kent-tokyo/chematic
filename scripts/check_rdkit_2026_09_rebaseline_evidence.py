#!/usr/bin/env python3
"""Check the published-npm 2026.03.6 -> 2026.09.1 rebaseline packet.

This does not turn the missing 2026.09.1 Python/native lanes into successes.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "validation" / "results"
DAY = "2026-10-04"
CORPUS = ROOT / "validation" / "benchmark_corpora" / "rdkit-js-browser-10k-v1.smi"
VERSIONS = ("2026.03.6", "2026.09.1")


def load(name: str) -> dict:
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def check() -> None:
    corpus_hash = sha256(CORPUS.read_bytes())
    availability = load(f"rdkit-2026-09-1-artifact-availability-{DAY}.json")
    require(availability["github_release"]["http_status"] == 200, "release missing")
    require(availability["pypi"]["http_status"] == 404, "Python availability snapshot changed")
    require(availability["npm"]["http_status"] == 200, "npm availability snapshot changed")

    artifacts: dict[str, dict] = {}
    for version in VERSIONS:
        provenance = load(f"rdkit-rebaseline-npm-provenance-v1.0.33-vs-{version}-{DAY}.json")
        require(provenance["inputs"]["sealed_accuracy_cohort_reused"] is False, "sealed corpus used")
        require(provenance["inputs"]["corpora"][0]["sha256"] == corpus_hash, "provenance corpus changed")
        require(len(provenance["lanes"]) == 1, "expected one provenance lane")
        lane = provenance["lanes"][0]
        require(lane["runtime"]["version"] == version, f"{version}: runtime mismatch")
        require(lane["missing_dimensions"] == [], f"{version}: missing provenance")
        require(lane["package"]["registry"]["status"] == "recorded", f"{version}: registry lock missing")
        require(lane["package"]["registry"]["integrity"].startswith("sha512-"), f"{version}: no SRI")
        artifact = lane["package"]["artifact"]
        require(artifact["sha256"] == lane["package"]["artifact_archive_sha256"], f"{version}: tarball hash mismatch")
        artifacts[version] = lane

        parity = load(f"rdkit-rebaseline-npm-morgan-parity-v1.0.33-vs-{version}-{DAY}.json")
        result = parity["result"]
        require(result["rdkit_version"] == version, f"{version}: parity runtime mismatch")
        require(result["compared_rows"] == 10000, f"{version}: partial parity corpus")
        require((result["exact_matches"], result["supported_rows"], result["unsupported_rows"]) == (9999, 9999, 1), f"{version}: Morgan counts changed")
        require(result["unsupported_sample"][0]["index"] == 8341, f"{version}: refusal row changed")
        require(result["first_mismatch"] is None, f"{version}: Morgan mismatch")

        browser = load(f"rdkit-rebaseline-npm-browser-runtime-v1.0.33-vs-{version}-{DAY}.json")
        require(browser["corpus"]["sha256"] == corpus_hash, f"{version}: browser corpus changed")
        require(browser["configuration"]["rows"] == 10000, f"{version}: partial browser corpus")
        require(browser["configuration"]["repetitions"] == 20, f"{version}: browser repetitions changed")
        require(browser["configuration"]["prepared_mode"] == "split", f"{version}: prepared operation changed")
        require(len(browser["raw_runs"]["rdkit"]) == 20, f"{version}: RDKit raw runs missing")
        require(len(browser["raw_runs"]["chematic"]) == 20, f"{version}: CheMatic raw runs missing")
        wasm = next(asset for asset in lane["package"]["assets"] if asset["path"] == "dist/RDKit_minimal.wasm")
        require(wasm["sha256"] == browser["artifacts"]["rdkit_wasm"]["sha256"], f"{version}: WASM artifact changed")
        require(browser["artifacts"]["schematic_package"]["version"] == "1.0.33", f"{version}: CheMatic package changed")
    old_browser = load(f"rdkit-rebaseline-npm-browser-runtime-v1.0.33-vs-{VERSIONS[0]}-{DAY}.json")
    new_browser = load(f"rdkit-rebaseline-npm-browser-runtime-v1.0.33-vs-{VERSIONS[1]}-{DAY}.json")
    require(old_browser["artifacts"]["schematic_wasm"]["sha256"] == new_browser["artifacts"]["schematic_wasm"]["sha256"], "CheMatic WASM differs between arms")

    comparison = load(f"rdkit-rebaseline-npm-oracle-delta-{VERSIONS[0]}-to-{VERSIONS[1]}-{DAY}.json")
    counts = comparison["counts"]
    require(comparison["old"]["runtime_version"] == VERSIONS[0], "old oracle runtime mismatch")
    require(comparison["new"]["runtime_version"] == VERSIONS[1], "new oracle runtime mismatch")
    require(comparison["corpus"]["sha256"] == corpus_hash, "oracle corpus changed")
    require((counts["rows"], counts["query_count"], counts["smarts_cells_compared"]) == (10000, 31, 310000), "oracle comparison incomplete")
    require((counts["parse_status_differences"], counts["canonical_spelling_differences"], counts["morgan_bit_differences"], counts["smarts_cell_differences"]) == (0, 0, 0, 12), "oracle delta counts changed")
    require(comparison["old"]["query_parse_failures"] == comparison["new"]["query_parse_failures"] == [], "query parser failures")
    for version, name in zip(VERSIONS, ("old", "new")):
        package = artifacts[version]["package"]
        require(comparison[name]["package_json"]["sha256"] == package["package_json"]["sha256"], f"{version}: oracle package changed")
        wasm = next(asset for asset in package["assets"] if asset["path"] == "dist/RDKit_minimal.wasm")
        require(comparison[name]["wasm"]["sha256"] == wasm["sha256"], f"{version}: oracle WASM changed")

    compressed = (ROOT / comparison["rows_output"]["path"]).read_bytes()
    require(sha256(compressed) == comparison["rows_output"]["compressed_sha256"], "compressed rows hash changed")
    body = gzip.decompress(compressed)
    require(sha256(body) == comparison["rows_output"]["uncompressed_sha256"], "rows content hash changed")
    status: Counter[tuple[str, str]] = Counter()
    changed: Counter[str] = Counter()
    for index, line in enumerate(body.splitlines()):
        row = json.loads(line)
        require(row["input_index"] == index, f"row {index}: index mismatch")
        status[(row["old_status"], row["new_status"])] += 1
        for difference in row["differences"]:
            require(difference["operation"] == "smarts", f"row {index}: unexpected difference")
            changed[difference["query"]] += 1
    require(index + 1 == 10000, "row file incomplete")
    require(status == {("ok", "ok"): 10000}, "row parse statuses changed")
    require(changed == {"[R2]": 6, "[R3]": 6}, "SMARTS difference families changed")
    print("RDKit 2026.09.1 npm rebaseline evidence OK: 10k rows, 310k SMARTS cells, 20 browser runs per arm")


if __name__ == "__main__":
    try:
        check()
    except (OSError, KeyError, IndexError, TypeError, ValueError) as exc:
        raise SystemExit(f"RDKit 2026.09.1 rebaseline evidence invalid: {exc}") from exc
