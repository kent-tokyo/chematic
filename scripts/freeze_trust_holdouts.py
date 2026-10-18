#!/usr/bin/env python3
"""Freeze outcome-blind trust holdouts from the exposed RDKit split.

The selection depends only on source text, fixed seeds, and input-side syntax.
No chematic or oracle operation is executed while choosing rows.  Generated
files intentionally contain inputs and protocols, never expected outputs.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "validation/rdkit_accuracy_split_v2/exposed_holdout.smi"
OUT = ROOT / "validation/trust_holdouts/v1"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def rank(seed: str, index: int, smiles: str) -> str:
    return sha(f"{seed}\0{index}\0{smiles}".encode())


def write_jsonl(path: Path, rows: list[dict[str, object]]) -> None:
    path.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in rows))


def main() -> None:
    source_bytes = SOURCE.read_bytes()
    molecules = [
        line.split()[0]
        for line in source_bytes.decode().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    OUT.mkdir(parents=True, exist_ok=True)

    ranked_mcs = sorted(
        ((rank("chematic-trust-mcs-v1", i, text), i, text)
         for i, text in enumerate(molecules))
    )[:128]
    mcs_rows = []
    modes = [("elements", "order"), ("any", "order"), ("elements", "any"), ("any", "any")]
    for pair_index in range(64):
        left = ranked_mcs[2 * pair_index]
        right = ranked_mcs[2 * pair_index + 1]
        atom_compare, bond_compare = modes[pair_index % len(modes)]
        mcs_rows.append({
            "id": f"mcs-{pair_index:03d}",
            "left_source_index": left[1],
            "left": left[2],
            "right_source_index": right[1],
            "right": right[2],
            "atom_compare": atom_compare,
            "bond_compare": bond_compare,
        })

    stereo_candidates = [
        (rank("chematic-trust-stereo3d-v1", i, text), i, text)
        for i, text in enumerate(molecules)
        if any(token in text for token in ("@", "/", "\\"))
    ]
    stereo_rows = [
        {"id": f"stereo3d-{n:03d}", "source_index": i, "smiles": text,
         "random_seed": 20261018 + n, "conformers": 8}
        for n, (_, i, text) in enumerate(sorted(stereo_candidates)[:96])
    ]

    ranked_sdf = sorted(
        ((rank("chematic-trust-sdf-attrs-v1", i, text), i, text)
         for i, text in enumerate(molecules))
    )[:64]
    sdf_rows = [
        {
            "id": f"sdf-{n:03d}",
            "source_index": i,
            "smiles": text,
            "title": f"holdout record {n}",
            "fields": {
                "row": str(i),
                "multiline": f"first line {n}\nsecond line {i}",
                "unicode": f"mu={n}\u03bc",
                "empty": "",
            },
            "enhanced_stereo_probe": any(token in text for token in ("@", "/", "\\")),
        }
        for n, (_, i, text) in enumerate(ranked_sdf)
    ]

    outputs = {
        "mcs-comparison.jsonl": mcs_rows,
        "sdf-attributes.jsonl": sdf_rows,
        "stereo-3d.jsonl": stereo_rows,
    }
    for name, rows in outputs.items():
        write_jsonl(OUT / name, rows)

    manifest = {
        "schema": "chematic.trust-holdout-manifest.v1",
        "frozen_on": "2026-10-18",
        "selection_policy": (
            "SHA-256 rank of source index and input text with lane-specific fixed seed; "
            "selection does not execute chematic, RDKit, COSMolKit, or inspect outcomes"
        ),
        "source": {
            "path": str(SOURCE.relative_to(ROOT)),
            "sha256": sha(source_bytes),
            "rows": len(molecules),
        },
        "lanes": {
            name: {
                "rows": len(rows),
                "sha256": sha((OUT / name).read_bytes()),
                "contains_expected_outputs": False,
            }
            for name, rows in outputs.items()
        },
        "use_policy": (
            "Freeze before measuring current or future candidates. Record every row; "
            "do not tune on these inputs, replace failures, or pool them with development data."
        ),
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
