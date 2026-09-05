#!/usr/bin/env python3
"""Generate the small, versioned machine-readable release metadata document."""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", required=True, help="release version without the v prefix")
    parser.add_argument("--commit", required=True, help="40-character release commit")
    parser.add_argument("--released-at", required=True, help="RFC3339 release timestamp")
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not re.fullmatch(r"\d+\.\d+\.\d+", args.version):
        raise SystemExit("--version must be a semantic version such as 1.0.5")
    if not re.fullmatch(r"[0-9a-f]{40}", args.commit):
        raise SystemExit("--commit must be a 40-character lowercase git commit")
    if not args.released_at.endswith("Z") and "+" not in args.released_at:
        raise SystemExit("--released-at must be RFC3339 with a timezone")

    document = {
        "schema_version": 1,
        "product": "chematic",
        "release": {
            "version": args.version,
            "tag": f"v{args.version}",
            "commit": args.commit,
            "released_at": args.released_at,
        },
        "packages": {
            "rust": {
                "name": "chematic",
                "version": args.version,
                "registry_url": f"https://crates.io/crates/chematic/{args.version}",
            },
            "python": {
                "name": "chematic",
                "version": args.version,
                "registry_url": f"https://pypi.org/project/chematic/{args.version}/",
            },
            "npm": {
                "name": "@kent-tokyo/chematic",
                "version": args.version,
                "registry_url": f"https://www.npmjs.com/package/@kent-tokyo/chematic/v/{args.version}",
            },
        },
        "wasm": {
            "artifact": "crates/chematic-wasm/pkg/chematic_wasm_bg.wasm",
            "measurements": [
                {
                    "status": "historical",
                    "release_or_snapshot": "v1.0.2 release candidate",
                    "raw_bytes": 3300452,
                    "compressed_bytes": 1213008,
                    "command": "wasm-pack 0.13.1; wasm-opt 130 -O3; gzip -9",
                    "runtime": "macOS arm64",
                    "source": "benchmarks/2026-09-04-wasm-size.md",
                }
            ],
        },
        "mcp": {
            "tool_count": 20,
            "transport": ["stdio"],
            "network_enabled_tools": ["pubchem_lookup"],
        },
        "benchmarks": {
            "historical": [
                {
                    "path": "benchmarks/2026-09-04-wasm-size.md",
                    "status": "historical",
                    "versions": ["v1.0.2 release candidate"],
                    "commit": "not recorded in the source benchmark",
                    "corpus": "optimized chematic WASM artifact",
                    "hardware": "macOS arm64",
                },
                {
                    "path": "benchmarks/2026-09-04-streaming-formats.md",
                    "status": "historical",
                    "versions": ["chematic 1.0.3", "RDKit 2025.09.3"],
                    "commit": "not recorded in the source benchmark",
                    "corpus": "checked-in two-record SDF and two-frame XYZ fixtures",
                    "hardware": "macOS arm64",
                },
            ]
        },
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
