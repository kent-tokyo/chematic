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
    "canonical-identity": ("Exact canonical-SMILES", "P2"),
    "v3000-parity": ("Exact canonical-SMILES", "P2"),
    "force-fields": ("MD/UFF/MMFF94", "P2"),
    "periodic-neighbors": ("periodic neighbor", "P2"),
    "symmetry-search": ("symmetry-heavy", "P2"),
    "browser-agent": ("browser and agent", "P3"),
    "reaction-breadth": ("reaction/SMARTS/medicinal", "P4"),
    "reaction-quality": ("curated reaction/query", "P4"),
    "force-field-quality": ("MMFF94/UFF typing", "P5"),
    "independent-review": ("S5 independent", "S5"),
    "maintenance": ("S6 continuous", "S6"),
    "accuracy-contract": ("A0 — Accuracy evidence contract", "A0 / P0"),
    "accuracy-descriptors": ("A1 — Perception and descriptor", "A1 / P1/P2"),
    "accuracy-stereo": ("A2 — Stereo and identity", "A2 / P2"),
    "accuracy-retrieval": ("A3 — Fingerprint and retrieval", "A3 / P2/P3"),
    "accuracy-workflows": ("A4 — Workflow accuracy", "A4 / P1/P4"),
    "accuracy-gold": ("A5 — Independent accuracy", "A5 / P0/P6"),
    "accuracy-3d": ("A6 — 3D accuracy", "A6 / P5"),
}

# Abandoned work stays discoverable without inflating the active backlog.
HISTORICAL_MARKERS = {
    "performance-stretch": ("1.10x", "Performance stretch target"),
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

    for name, (roadmap_marker, disposition_marker) in HISTORICAL_MARKERS.items():
        if roadmap_marker not in roadmap:
            errors.append(f"{name}: historical roadmap marker is missing")
        if any(roadmap_marker in line for line in unchecked):
            errors.append(f"{name}: historical area must not be an unchecked item")
        if disposition_marker not in disposition:
            errors.append(f"{name}: disposition marker '{disposition_marker}' is missing")

    for required_class in ("`local-open`", "`local-toolchain`", "`external`", "`historical`"):
        if required_class not in disposition:
            errors.append(f"dependency class {required_class} is missing")

    # Multiple stable area markers may intentionally live in one compound
    # roadmap item (for example canonical output and cross-engine V3000
    # fixtures).  The important invariant is that every unchecked item is
    # represented in the disposition table, not a one-marker/one-line ratio.
    for line in unchecked:
        if not any(marker in line for marker, _ in AREA_MARKERS.values()):
            errors.append(f"unchecked roadmap item has no disposition marker: {line}")
    if errors:
        print("roadmap disposition failures:")
        print("\n".join(f"- {error}" for error in errors))
        return 1
    print(
        f"Roadmap disposition OK: {len(unchecked)} unchecked areas classified; "
        f"{len(HISTORICAL_MARKERS)} historical area excluded"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
