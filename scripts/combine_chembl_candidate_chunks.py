#!/usr/bin/env python3
"""Combine independently hashed local ChEMBL acquisition chunks without hiding provenance."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("chunks", type=Path, nargs="+")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    rows: list[str] = []
    chunks: list[dict[str, object]] = []
    for chunk in args.chunks:
        directory = chunk.resolve()
        source = directory / "source.smi"
        metadata_path = directory / "source-metadata.json"
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if metadata.get("source_sha256") != sha256(source):
            parser.error(f"source hash mismatch in {directory}")
        rows.extend(line for line in source.read_text(encoding="utf-8").splitlines() if line)
        chunks.append({"path": str(directory), "source_sha256": metadata["source_sha256"], "responses": metadata["responses"], "accepted_rows": metadata["accepted_single_fragment_rows"]})
    source = output / "source.smi"
    source.write_text("".join(f"{row}\n" for row in rows), encoding="utf-8")
    metadata = {
        "source_url": "https://www.ebi.ac.uk/chembl/api/data/molecule.json",
        "source_release": "ChEMBL REST API multi-chunk snapshot; every response URL/hash retained",
        "license": "ChEMBL data: CC BY-SA 3.0",
        "license_url": "https://creativecommons.org/licenses/by-sa/3.0/",
        "acquired_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "source_sha256": sha256(source),
        "accepted_single_fragment_rows": len(rows),
        "chunks": chunks,
        "redistribution": "local evaluation input only; raw source is not committed to this repository",
        "status": "acquired_not_sealed",
    }
    (output / "source-metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": metadata["status"], "accepted_rows": len(rows), "source_sha256": metadata["source_sha256"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
