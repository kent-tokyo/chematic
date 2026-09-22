#!/usr/bin/env python3
"""Check that active roadmap packages and dependency classes stay explicit."""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ROADMAP = ROOT / "ROADMAP.md"
DISPOSITION = ROOT / "docs" / "roadmap-open-work.md"

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
    package_lines = [line for line in roadmap.splitlines() if re.match(r"^- \[[ x]\]", line)]
    errors: list[str] = []

    for package, title in ACCURACY_PACKAGES.items():
        marker = f"**{package} — {title}:**"
        matching = [line for line in package_lines if marker in line]
        if len(matching) != 1:
            errors.append(f"{package}: roadmap package must appear exactly once")
        elif package == "A0" and not matching[0].startswith("- [x]"):
            errors.append("A0: completed roadmap package must stay checked")
        elif package != "A0" and not matching[0].startswith("- [ ]"):
            errors.append(f"{package}: open roadmap package must stay unchecked")

    package_ids = {
        match.group(1)
        for line in package_lines
        if (match := re.search(r"\*\*(A[0-6]) —", line))
    }
    if package_ids != set(ACCURACY_PACKAGES):
        errors.append(f"unexpected accuracy package set: {sorted(package_ids)}")

    for dependency_class in DEPENDENCY_CLASSES:
        if dependency_class not in disposition:
            errors.append(f"dependency class {dependency_class} is missing")

    for required_heading in ("## Priority order", "## Product phases", "## v1.0.16 candidate boundary"):
        if required_heading not in roadmap:
            errors.append(f"roadmap heading is missing: {required_heading}")

    if errors:
        print("roadmap disposition failures:")
        print("\n".join(f"- {error}" for error in errors))
        return 1

    print(
        f"Roadmap disposition OK: {len(package_ids)} accuracy packages (A0 complete) and "
        f"{len(DEPENDENCY_CLASSES)} dependency classes are explicit"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
