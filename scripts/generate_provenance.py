#!/usr/bin/env python3
"""Create and optionally sign a deterministic artifact provenance manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def revision() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return "UNAVAILABLE"


def created_at() -> str:
    epoch = os.environ.get("SOURCE_DATE_EPOCH")
    if epoch is not None:
        return datetime.fromtimestamp(int(epoch), timezone.utc).replace(microsecond=0).isoformat()
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument("--artifact", action="append", type=Path, default=[])
    args = parser.parse_args()

    artifacts = []
    for path in sorted(args.artifact, key=lambda item: str(item)):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        artifacts.append({"path": str(path), "sha256": digest, "size": path.stat().st_size})
    manifest = {
        "schema": "chematic.provenance.v1",
        "created": created_at(),
        "source_revision": revision(),
        "artifacts": artifacts,
        "notes": [
            "This manifest records hashes; cryptographic signing is performed by the release CI key.",
            "Do not treat a locally generated key as release provenance.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"Provenance written: {args.output} ({len(artifacts)} artifacts)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
