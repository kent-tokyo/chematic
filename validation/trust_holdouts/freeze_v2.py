#!/usr/bin/env python3
"""Freeze the post-candidate v2 trust holdout without inspecting outcomes.

The candidate commit must predate both external source acquisitions. Selection
uses only source identity, source text, syntax, and fixed SHA-256 ranks. Rows
already used by v1 are excluded so v2 is an independent evaluation set.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "validation/trust_holdouts/v2"
V1 = ROOT / "validation/trust_holdouts/v1"
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
    parser.add_argument("--nci-acquired-at", required=True)
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
    if payload.get("schema") != "chematic.pubchem-holdout-source.v2":
        raise SystemExit("unexpected PubChem source schema")
    rows = [(f"CID-{row['cid']}", row["smiles"]) for row in payload["records"]]
    return data, payload, rows


def source_row(source_id: str, smiles: str) -> dict[str, object]:
    return {
        "source_id": source_id,
        "source_row_sha256": sha(f"{source_id}\t{smiles}\n".encode()),
        "smiles": smiles,
    }


def v1_source_ids() -> set[str]:
    used: set[str] = set()
    for path in V1.glob("*.jsonl"):
        for line in path.read_text().splitlines():
            row = json.loads(line)
            for value in row.values():
                if isinstance(value, dict) and "source_id" in value:
                    used.add(value["source_id"])
            if "source_id" in row:
                used.add(row["source_id"])
    return used


def parse_time(value: str) -> datetime:
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        raise SystemExit("timestamps must include a timezone")
    return parsed.astimezone(timezone.utc)


def main() -> None:
    args = parse_args()
    frozen_at = parse_time(args.candidate_frozen_at)
    nci_acquired_at = parse_time(args.nci_acquired_at)
    nci_bytes, nci = parse_nci(args.nci_source)
    pubchem_bytes, pubchem_payload, pubchem = parse_pubchem(args.pubchem_source)
    pubchem_acquired_at = parse_time(pubchem_payload["acquired_at"])
    if nci_acquired_at <= frozen_at or pubchem_acquired_at <= frozen_at:
        raise SystemExit("all external sources must be acquired after candidate freeze")

    excluded = v1_source_ids()
    nci = [row for row in nci if row[0] not in excluded]
    pubchem = [row for row in pubchem if row[0] not in excluded]
    OUT.mkdir(parents=True, exist_ok=True)

    ranked_mcs = sorted(
        (rank("chematic-trust-mcs-v2", source_id, text), source_id, text)
        for source_id, text in nci
    )[:144]
    modes = [("elements", "order"), ("any", "order"), ("elements", "any"), ("any", "any")]
    mcs_rows: list[dict[str, object]] = []
    cursor = 0
    for row_index in range(64):
        width = 3 if row_index >= 48 else 2
        selected = ranked_mcs[cursor:cursor + width]
        cursor += width
        atom_compare, bond_compare = modes[row_index % len(modes)]
        mcs_rows.append({
            "id": f"mcs-{row_index:03d}",
            "inputs": [source_row(item[1], item[2]) for item in selected],
            "atom_compare": atom_compare,
            "bond_compare": bond_compare,
        })

    stereo_candidates = sorted(
        (rank("chematic-trust-stereo3d-v2", source_id, text), source_id, text)
        for source_id, text in pubchem
        if any(token in text for token in ("@", "/", "\\"))
    )
    if len(stereo_candidates) < 96:
        raise SystemExit(f"PubChem source has only {len(stereo_candidates)} new stereo candidates")
    stereo_rows = [
        {
            "id": f"stereo3d-{n:03d}",
            **source_row(source_id, text),
            "random_seed": 20261024 + n,
            "conformers": 8,
        }
        for n, (_, source_id, text) in enumerate(stereo_candidates[:96])
    ]

    ranked_sdf = sorted(
        (rank("chematic-trust-sdf-attrs-v2", source_id, text), source_id, text)
        for source_id, text in nci
    )[:64]
    sdf_rows = [
        {
            "id": f"sdf-{n:03d}",
            **source_row(source_id, text),
            "title": f"holdout v2 record {n}",
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
    for name, selected_rows in outputs.items():
        write_jsonl(OUT / name, selected_rows)

    manifest = {
        "schema": "chematic.trust-holdout-manifest.v2",
        "candidate_frozen_commit": args.candidate_commit,
        "candidate_frozen_at": args.candidate_frozen_at,
        "frozen_on": "2026-10-24",
        "selection_policy": (
            "Post-candidate external sources; SHA-256 rank of source identity and input text "
            "with lane-specific fixed v2 seed. v1 source IDs are excluded. Selection does not "
            "execute a chemistry engine or inspect outcomes."
        ),
        "sources": {
            "nci": {
                "url": "https://cactus.nci.nih.gov/download/nci/NCISMA99.sdz",
                "archive_sha256": NCI_ARCHIVE_SHA256,
                "decompressed_sha256": sha(nci_bytes),
                "rows": len(nci),
                "acquired_at": args.nci_acquired_at,
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
        "independence": {
            "excluded_v1_source_ids": len(excluded),
            "overlap_with_v1": 0,
        },
        "lanes": {
            name: {
                "rows": len(selected_rows),
                "sha256": sha((OUT / name).read_bytes()),
                "contains_expected_outputs": False,
            }
            for name, selected_rows in outputs.items()
        },
        "use_policy": (
            "One-shot scoring after this manifest is committed. Record every row; do not tune "
            "on these inputs, replace failures, or pool them with development data."
        ),
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
