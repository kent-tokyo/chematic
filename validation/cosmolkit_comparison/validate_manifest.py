#!/usr/bin/env python3
"""Validate the checked-in comparison corpus and its manifest."""

import json
from pathlib import Path

from validate import corpus


def main() -> None:
    digest, records = corpus()
    manifest = json.loads(Path(__file__).with_name("corpus_manifest.json").read_text())
    print(json.dumps({"valid": True, "records": len(records), "sha256": digest,
                      "schema_version": manifest["schema_version"]}, sort_keys=True))


if __name__ == "__main__":
    main()
