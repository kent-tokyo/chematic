#!/usr/bin/env python3
"""Check that every unchecked roadmap area has an explicit disposition."""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ROADMAP = ROOT / "ROADMAP.md"
DISPOSITION = ROOT / "docs" / "roadmap-open-work.md"

# These are deliberately stable area markers rather than exact prose. Roadmap
# wording may grow, but an unchecked area must remain visible in the
# disposition table and in one dependency class.
AREA_MARKERS = {
    "equivalent-operation": ("equivalent operations", "P1"),
    "canonical-identity": ("canonical atom-order", "P2"),
    "performance-stretch": ("1.10x", "Performance stretch target"),
    "force-fields": ("MD/UFF/MMFF94", "P2"),
    "periodic-neighbors": ("periodic neighbor", "P2"),
    "symmetry-search": ("symmetry-heavy", "P2"),
    "binding-manifest": ("one fixture schema", "P3"),
    "browser-agent": ("browser and agent", "P3"),
    "reaction-breadth": ("reaction/SMARTS/medicinal", "P4"),
    "reaction-quality": ("curated reaction/query", "P4"),
    "force-field-quality": ("MMFF94/UFF typing", "P5"),
    "independent-review": ("S5 independent", "S5"),
    "maintenance": ("S6 continuous", "S6"),
}


def main() -> int:
    roadmap = ROADMAP.read_text(encoding="utf-8")
    disposition = DISPOSITION.read_text(encoding="utf-8")
    unchecked = [line for line in roadmap.splitlines() if re.match(r"^- \[ \]", line)]
    errors: list[str] = []
    for name, (roadmap_marker, disposition_marker) in AREA_MARKERS.items():
        if roadmap_marker not in roadmap:
            errors.append(f"{name}: roadmap marker is missing")
            continue
        if not any(roadmap_marker in line for line in unchecked):
            errors.append(f"{name}: marker is not represented by an unchecked item")
        if disposition_marker not in disposition:
            errors.append(f"{name}: disposition marker '{disposition_marker}' is missing")

    for required_class in ("`local-open`", "`local-toolchain`", "`external`", "`historical`"):
        if required_class not in disposition:
            errors.append(f"dependency class {required_class} is missing")

    if len(unchecked) != len(AREA_MARKERS):
        errors.append(
            f"unchecked roadmap count changed: expected {len(AREA_MARKERS)}, got {len(unchecked)}"
        )
    if errors:
        print("roadmap disposition failures:")
        print("\n".join(f"- {error}" for error in errors))
        return 1
    print(f"Roadmap disposition OK: {len(unchecked)} unchecked areas classified")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
