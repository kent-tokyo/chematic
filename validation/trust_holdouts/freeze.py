#!/usr/bin/env python3
"""Freeze post-candidate, outcome-blind trust holdouts.

Inputs must be acquired after the candidate commit is frozen. Selection uses
only source identifiers, source text, syntax, and fixed SHA-256 ranks. No
chemistry engine is imported or executed by this script.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "validation/trust_holdouts/v1"
NCI_SHA256 = "031736b438a9596d89cf11802fa97934c3878ebbf021924c417a57e02c7930cd"
NCI_ARCHIVE_SHA256 = "3da792f624dc94caa9a3e9967fc7b1eed0ca79e0bd2d044f808cf916680a17ad"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def rank(seed: str, source_id: str, smiles: str) -> str:
    return sha(f"{seed}\0{source_id}\0{smiles}".encode())


def write_jsonl(path: Path, rows: list[dict[str, object]]) -> None:
    path.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in rows))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--nci-source", type=Path, required=True)
    parser.add_argument("--pubchem-source", type=Path, required=True)
    parser.add_argument("--candidate-commit", required=True)
    parser.add_argument("--candidate-frozen-at", required=True)
    return parser.parse_args()


def parse_nci(path: Path) -> tuple[bytes, list[tuple[str, str]]]:
    data = path.read_bytes()
    if sha(data) != NCI_SHA256:
        raise SystemExit("NCI source hash does not match pinned NCISMA99")
    rows: list[tuple[str, str]] = []
    for line in data.decode().splitlines():
        parts = line.split(maxsplit=1)
        if len(parts) == 2:
            rows.append((f"NCI-{parts[0]}", parts[1].strip()))
    return data, rows


def parse_pubchem(path: Path) -> tuple[bytes, dict[str, object], list[tuple[str, str]]]:
    data = path.read_bytes()
    payload = json.loads(data)
    if payload.get("schema") != "chematic.pubchem-holdout-source.v1":
        raise SystemExit("unexpected PubChem source schema")
    rows = [(f"CID-{row['cid']}", row["smiles"]) for row in payload["records"]]
    return data, payload, rows


def source_row(source_id: str, smiles: str) -> dict[str, object]:
    return {
        "source_id": source_id,
        "source_row_sha256": sha(f"{source_id}\t{smiles}\n".encode()),
        "smiles": smiles,
    }


def main() -> None:
    args = parse_args()
    frozen_at = datetime.fromisoformat(args.candidate_frozen_at.replace("Z", "+00:00"))
    if frozen_at.tzinfo is None:
        raise SystemExit("candidate-frozen-at must include a timezone")
    nci_bytes, nci = parse_nci(args.nci_source)
    pubchem_bytes, pubchem_payload, pubchem = parse_pubchem(args.pubchem_source)
    acquired_at = datetime.fromisoformat(pubchem_payload["acquired_at"].replace("Z", "+00:00"))
    if acquired_at <= frozen_at.astimezone(timezone.utc):
        raise SystemExit("PubChem source must be acquired after candidate freeze")
    OUT.mkdir(parents=True, exist_ok=True)

    ranked_mcs = sorted(
        (rank("chematic-trust-mcs-v1", source_id, text), source_id, text)
        for source_id, text in nci
    )[:128]
    modes = [("elements", "order"), ("any", "order"), ("elements", "any"), ("any", "any")]
    mcs_rows = []
    for pair_index in range(64):
        left = ranked_mcs[2 * pair_index]
        right = ranked_mcs[2 * pair_index + 1]
        atom_compare, bond_compare = modes[pair_index % len(modes)]
        mcs_rows.append({
            "id": f"mcs-{pair_index:03d}",
            "left": source_row(left[1], left[2]),
            "right": source_row(right[1], right[2]),
            "atom_compare": atom_compare,
            "bond_compare": bond_compare,
        })

    stereo_candidates = [
        (rank("chematic-trust-stereo3d-v1", source_id, text), source_id, text)
        for source_id, text in pubchem
        if any(token in text for token in ("@", "/", "\\"))
    ]
    if len(stereo_candidates) < 96:
        raise SystemExit(f"PubChem source has only {len(stereo_candidates)} stereo candidates")
    stereo_rows = [
        {
            "id": f"stereo3d-{n:03d}",
            **source_row(source_id, text),
            "random_seed": 20261018 + n,
            "conformers": 8,
        }
        for n, (_, source_id, text) in enumerate(sorted(stereo_candidates)[:96])
    ]

    ranked_sdf = sorted(
        (rank("chematic-trust-sdf-attrs-v1", source_id, text), source_id, text)
        for source_id, text in nci
    )[:64]
    sdf_rows = [
        {
            "id": f"sdf-{n:03d}",
            **source_row(source_id, text),
            "title": f"holdout record {n}",
            "fields": {
                "source_id": source_id,
                "multiline": f"first line {n}\nsecond line {source_id}",
                "unicode": f"mu={n}\u03bc",
                "empty": "",
            },
            "enhanced_stereo_probe": any(token in text for token in ("@", "/", "\\")),
        }
        for n, (_, source_id, text) in enumerate(ranked_sdf)
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
        "candidate_frozen_commit": args.candidate_commit,
        "candidate_frozen_at": args.candidate_frozen_at,
        "frozen_on": "2026-10-18",
        "selection_policy": (
            "Post-candidate external sources; SHA-256 rank of source identity and input text "
            "with lane-specific fixed seed. Selection does not execute chematic, RDKit, "
            "COSMolKit, or inspect outcomes."
        ),
        "sources": {
            "nci": {
                "url": "https://cactus.nci.nih.gov/download/nci/NCISMA99.sdz",
                "archive_sha256": NCI_ARCHIVE_SHA256,
                "decompressed_sha256": sha(nci_bytes),
                "rows": len(nci),
                "acquired_at": "2026-10-18T13:48:19+09:00",
                "acquired_after_candidate_freeze": True,
            },
            "pubchem": {
                "url_template": pubchem_payload["url_template"],
                "bundle_sha256": sha(pubchem_bytes),
                "records": len(pubchem),
                "acquired_at": pubchem_payload["acquired_at"],
                "requests": pubchem_payload["requests"],
                "acquired_after_candidate_freeze": True,
            },
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
            "Freeze before measuring this or future candidates. Record every row; do not tune "
            "on these inputs, replace failures, or pool them with development data."
        ),
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
