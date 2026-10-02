#!/usr/bin/env python3
"""Check that active roadmap packages and dependency classes stay explicit."""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ROADMAP = ROOT / "ROADMAP.md"
DISPOSITION = ROOT / "docs" / "roadmap-open-work.md"


def roadmap_release_line(roadmap: str) -> str:
    match = re.search(r"Release line: \*\*v([^*]+)\*\*", roadmap)
    if match is None:
        raise ValueError("roadmap release line is missing")
    return match.group(1)


ACCURACY_PACKAGES = {
    "A0": "Evaluation contract",
    "A1": "Perception and descriptors",
    "A2": "Stereo and identity",
    "A3": "Fingerprints and retrieval",
    "A4": "Workflows and interchange",
    "A5": "Independent adjudication",
    "A6": "3D and force fields",
}

DEPENDENCY_CLASSES = (
    "`local-open`",
    "`toolchain-open`",
    "`data-sealed`",
    "`external-open`",
    "`historical`",
)


def main() -> int:
    roadmap = ROADMAP.read_text(encoding="utf-8")
    disposition = DISPOSITION.read_text(encoding="utf-8")
    release_line = roadmap_release_line(roadmap)
    cargo = (ROOT / "Cargo.toml").read_text(encoding="utf-8")
    workspace_version = re.search(r'^version\s*=\s*"([^"]+)"\s*$', cargo, re.MULTILINE)
    errors: list[str] = []

    if workspace_version is None or release_line != workspace_version.group(1):
        errors.append("roadmap release line does not match Cargo workspace version")

    for document_name, document in (("ROADMAP", roadmap), ("open-work ledger", disposition)):
        package_rows = [
            match.groups()
            for line in document.splitlines()
            if (match := re.match(r"^\|\s*(A[0-6])\s+([^|]+?)\s*\|\s*([^|]+?)\s*\|", line))
        ]
        for package, title in ACCURACY_PACKAGES.items():
            matching = [row for row in package_rows if row[0] == package]
            if len(matching) != 1 or matching[0][1].strip() != title:
                errors.append(f"{document_name}: {package} package/title must appear exactly once")
                continue
            state = matching[0][2].strip().lower()
            if (package == "A0" and state != "complete") or (
                package != "A0" and state.startswith("complete")
            ):
                errors.append(f"{document_name}: {package} completion state is incorrect")

    for dependency_class in DEPENDENCY_CLASSES:
        if dependency_class not in disposition:
            errors.append(f"dependency class {dependency_class} is missing")

    for required_heading in (
        "## Current position",
        "## Priority order",
        "## Accuracy packages",
        "## Product phases",
    ):
        if required_heading not in roadmap:
            errors.append(f"roadmap heading is missing: {required_heading}")

    if errors:
        print("roadmap disposition failures:")
        print("\n".join(f"- {error}" for error in errors))
        return 1

    print(
        f"Roadmap disposition OK: {len(ACCURACY_PACKAGES)} accuracy packages (A0 complete) and "
        f"{len(DEPENDENCY_CLASSES)} dependency classes are explicit"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
