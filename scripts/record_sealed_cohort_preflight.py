#!/usr/bin/env python3
"""Create a check-in-safe evidence summary from a local sealed-cohort preflight."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.manifest.read_text(encoding="utf-8"))
    if source.get("protocol") != "rdkit-accuracy-sealed-cohort-v1":
        parser.error("not a sealed cohort preflight manifest")
    raw_source = source["source"]
    report = {
        "schema": "chematic.sealed-cohort-preflight-evidence.v1",
        "status": source["status"],
        "sealed": source["sealed"],
        "reason": source["reason"],
        "source": {
            "url": raw_source["source_url"],
            "release": raw_source["source_release"],
            "license": raw_source["license"],
            "license_url": raw_source.get("license_url"),
            "acquired_at": raw_source["acquired_at"],
            "source_sha256": raw_source["source_sha256"],
            "input_rows": raw_source["input_rows"],
            "response_hashes": [response for chunk in raw_source.get("chunks", []) for response in chunk["responses"]],
        },
        "identity_audit": {
            "oracle": {"engine": "RDKit", "version": "2025.09.3"},
            "sha256": source["identity_audit"]["sha256"],
            "reference_audit_sha256": [row["sha256"] for row in source["identity_audit"]["reference_audits"]],
            "keys": source["identity_audit"]["required_keys"],
        },
        "deduplication": source["deduplication"],
        "splits": source["splits"],
        "candidate_freeze": source["candidate_freeze"],
        "unused_data_attestation": source.get("unused_data_attestation"),
        "raw_input": "local-only; not committed",
    }
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "sealed_rows": report["splits"]["sealed_holdout"]["rows"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
