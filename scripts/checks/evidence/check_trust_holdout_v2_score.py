#!/usr/bin/env python3
"""Verify the immutable macOS score for the frozen v2 trust holdout."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCORE = ROOT / "validation/results/trust_holdout_v2_macos_arm64.json"
EXPECTED_SHA256 = "95635e8338b6fed3e3fc579364d23ba2635a9768dc1c93ada9fffcab9cc12c1e"


def main() -> None:
    assert hashlib.sha256(SCORE.read_bytes()).hexdigest() == EXPECTED_SHA256
    score = json.loads(SCORE.read_text())
    assert score["candidate_commit"] == "29b3fdc9ddf17c26cdb65b73be101e8e82046841"
    assert score["wheel"]["sha256"] == "fcba25a934c021ccc8ec5cdc5ded01191ef7484f369a99e624c39529333481f6"
    assert score["versions"] == {
        "chematic": "1.2.1",
        "python": "3.13.6",
        "rdkit": "2026.09.1",
    }
    lanes = score["lanes"]
    assert len(lanes["mcs"]["rows"]) == 64
    assert len(lanes["sdf"]["rows"]) == 64
    assert len(lanes["stereo_3d"]["rows"]) == 96
    assert lanes["sdf"]["columns"]["graph"]["success"] == 64
    assert lanes["sdf"]["columns"]["properties"]["success"] == 64
    assert lanes["sdf"]["columns"]["stereo"]["success"] == 64
    assert lanes["stereo_3d"]["columns"]["joint_embed_status"]["success"] == 96
    assert lanes["stereo_3d"]["columns"]["stereo_graph"]["success"] == 96
    print("Frozen trust holdout v2 macOS score passed")


if __name__ == "__main__":
    main()
