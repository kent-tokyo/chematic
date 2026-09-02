#!/usr/bin/env python3
"""Generate a deterministic SPDX 2.3 SBOM from Cargo's locked metadata.

This deliberately uses only Cargo and the Python standard library so the
release evidence can be produced offline after dependencies are cached.
Registry checksums are retained when Cargo reports them; workspace packages
are identified by their local manifest path.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def cargo_metadata() -> dict:
    return json.loads(
        subprocess.check_output(
            ["cargo", "metadata", "--format-version", "1", "--locked"],
            text=True,
        )
    )


def created_at() -> str:
    epoch = os.environ.get("SOURCE_DATE_EPOCH")
    if epoch is not None:
        return datetime.fromtimestamp(int(epoch), timezone.utc).replace(microsecond=0).isoformat()
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args()

    metadata = cargo_metadata()
    root = Path(metadata["workspace_root"]).resolve()
    packages = sorted(metadata["packages"], key=lambda package: package["id"])
    documents = []
    relationships = []
    refs_by_name = {}
    for index, package in enumerate(packages):
        refs_by_name.setdefault(package["name"], "SPDXRef-Package-" + str(index))
    for index, package in enumerate(packages):
        ref = "SPDXRef-Package-" + str(index)
        source = package.get("source") or "local:" + str(
            (Path(package["manifest_path"]).parent.resolve()).relative_to(root)
        )
        checksums = package.get("checksum")
        document = {
            "SPDXID": ref,
            "name": package["name"],
            "versionInfo": package["version"],
            "downloadLocation": source,
            "filesAnalyzed": False,
            "licenseConcluded": "NOASSERTION",
            "licenseDeclared": "NOASSERTION",
            "copyrightText": "NOASSERTION",
        }
        if checksums:
            document["checksums"] = [{"algorithm": "SHA256", "checksumValue": checksums}]
        documents.append(document)
        for dependency in package["dependencies"]:
            target = refs_by_name.get(dependency["name"])
            if target is None:
                # Cargo metadata normally includes every resolved package; keep
                # the document valid if a future Cargo mode omits one.
                continue
            relationships.append(
                {
                    "spdxElementId": ref,
                    "relationshipType": "DEPENDS_ON",
                    "relatedSpdxElement": target,
                }
            )

    output = {
        "spdxVersion": "SPDX-2.3",
        "dataLicense": "CC0-1.0",
        "SPDXID": "SPDXRef-DOCUMENT",
        "name": "chematic-cargo-lock-sbom",
        "documentNamespace": "https://github.com/kent-tokyo/chematic/sbom/cargo-lock",
        "creationInfo": {
            "created": created_at(),
            "creators": ["Tool: chematic/scripts/generate_sbom.py"],
        },
        "packages": documents,
        "relationships": sorted(
            relationships,
            key=lambda relationship: (
                relationship["spdxElementId"], relationship["relatedSpdxElement"]
            ),
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(f"SBOM written: {args.output} ({len(documents)} packages)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
