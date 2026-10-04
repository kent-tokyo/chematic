#!/usr/bin/env python3
"""Gate checked Python reaction graphs, atom origins and template maps on 83 fixtures.

This is a source-candidate gate. It does not promote a local extension or wheel
to published-package evidence, and raw reaction embedding counts are diagnostic.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path
from zipfile import ZipFile

from rdkit import Chem, RDLogger, rdBase

if __package__:
    from .reaction_atom_provenance_gate import rdkit_sets, rust_sets
    from .reaction_template_map_gate import rdkit_map_sets, rust_map_sets
    from .run_reaction_compatibility_v2 import load_cases
    from .reaction_python_checked_candidates import python_row as candidate_row
else:
    from reaction_atom_provenance_gate import rdkit_sets, rust_sets
    from reaction_template_map_gate import rdkit_map_sets, rust_map_sets
    from run_reaction_compatibility_v2 import load_cases
    from reaction_python_checked_candidates import python_row as candidate_row

ROOT = Path(__file__).resolve().parents[1]
RDLogger.DisableLog("rdApp.*")
EXPECTED = {
    # v2_isotope_methanol_split and the three alanine stereo rows match
    # RDKit since #734 (were a refusal and typed unsupported).
    "graph_origin_map_match": 80,
    "joint_invalid_input": 3,
}


def source_artifact(chematic, wheel: Path | None) -> dict:
    if wheel is None:
        return {"kind": "editable_source_extension", "published": False}
    if wheel.suffix != ".whl" or not wheel.is_file():
        raise ValueError("--wheel must identify one existing wheel file")
    package_dir = Path(chematic.__file__).resolve().parent
    extensions = list(package_dir.glob("chematic*.so")) + list(package_dir.glob("chematic*.pyd"))
    if len(extensions) != 1:
        raise ValueError("installed chematic extension is missing or ambiguous")
    installed_digest = hashlib.sha256(extensions[0].read_bytes()).hexdigest()
    with ZipFile(wheel) as archive:
        entries = [name for name in archive.namelist()
                   if name.endswith("/" + extensions[0].name)]
        if len(entries) != 1:
            raise ValueError("wheel extension is missing or ambiguous")
        wheel_extension_digest = hashlib.sha256(archive.read(entries[0])).hexdigest()
    if installed_digest != wheel_extension_digest:
        raise ValueError("installed extension does not match the supplied wheel")
    return {"kind": "release_profile_source_wheel", "published": False,
            "wheel": wheel.name, "wheel_sha256": hashlib.sha256(wheel.read_bytes()).hexdigest(),
            "extension_sha256": installed_digest}


def rdkit_readable(case_id: str, smiles: str, atom_count: int) -> None:
    parsed = Chem.MolFromSmiles(smiles)
    if parsed is None or parsed.GetNumAtoms() != atom_count:
        raise ValueError(f"{case_id}: canonical product is not RDKit-readable")


def python_row(chematic, case: dict) -> dict:
    return candidate_row(chematic, case, rdkit_readable)


def candidates_rows(path: Path, cases: list[dict], hashes: dict) -> tuple[dict, list[dict]]:
    """Rows written by ``reaction_python_checked_candidates.py`` in another interpreter."""
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("schema") != "python-checked-reaction-candidates/v1":
        raise ValueError("unknown candidates schema")
    if data["fixtures"] != hashes:
        raise ValueError("candidates were produced from different fixtures")
    if sorted(data["rows"]) != sorted(case["id"] for case in cases):
        raise ValueError("candidate rows do not cover the fixture IDs")
    rows = [data["rows"][case["id"]] for case in cases]
    for case, row in zip(cases, rows):
        for product_set in row.get("sets", []):
            for item in product_set:
                rdkit_readable(case["id"], item["smiles"], len(item["atom_sources"]))
    return {**data["artifact"], "candidates_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}, rows


def classify(case: dict, candidate: dict) -> dict:
    status = candidate["status"]
    try:
        oracle_graph, oracle_origins, oracle_count = rdkit_sets(case)
        oracle_maps, map_count = rdkit_map_sets(case)
        if oracle_count != map_count:
            raise ValueError("RDKit reaction count changed between oracle passes")
    except ValueError as exc:
        outcome = "joint_invalid_input" if status == "typed_refusal" else "oracle_invalid"
        return {"id": case["id"], "outcome": outcome, "status": status, "detail": str(exc)}
    if status == "typed_unsupported":
        outcome = "typed_unsupported"
    elif status in {"typed_refusal", "partial_products"}:
        outcome = "typed_or_diagnosed_refusal" if candidate.get("reason") == "product_valence" else "unexpected_refusal"
    elif status in {"products", "no_match"}:
        graph, origins = rust_sets(case, {"sets_with_atom_sources": candidate["sets"]})
        maps = rust_map_sets({"sets": candidate["sets"]})
        outcome = ("graph_mismatch" if graph != oracle_graph else
                   "provenance_mismatch" if origins != oracle_origins else
                   "map_label_mismatch" if maps != oracle_maps else
                   "graph_origin_map_match")
    else:
        outcome = "unexpected_status"
    return {"id": case["id"], "outcome": outcome, "status": status,
            "reason": candidate.get("reason"), "diagnostics": candidate.get("diagnostics"),
            "rdkit_raw_product_sets": oracle_count}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=ROOT / "validation/reaction_product_parity_cases.json")
    parser.add_argument("--strata", type=Path, default=ROOT / "validation/reaction_product_parity_strata_v2.json")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--wheel", type=Path,
                        help="release-profile source wheel installed in this interpreter")
    parser.add_argument("--candidates", type=Path,
                        help="rows from reaction_python_checked_candidates.py (another interpreter)")
    parser.add_argument("--expected-rdkit", default="2026.03.6")
    args = parser.parse_args()
    if rdBase.rdkitVersion != args.expected_rdkit:
        parser.error(f"RDKit {rdBase.rdkitVersion} != {args.expected_rdkit}")
    if args.candidates and args.wheel:
        parser.error("--candidates and --wheel are exclusive")

    cases, hashes = load_cases(args.base, args.strata)
    if len(cases) != 83:
        raise ValueError(f"reaction fixture count changed: {len(cases)}")
    if args.candidates:
        artifact, candidates = candidates_rows(args.candidates, cases, hashes)
    else:
        import chematic

        artifact = source_artifact(chematic, args.wheel)
        candidates = [python_row(chematic, case) for case in cases]
    rows = [classify(case, candidate) for case, candidate in zip(cases, candidates)]
    counts = dict(sorted(Counter(row["outcome"] for row in rows).items()))
    report = {"schema": "source-python-checked-reaction-provenance/v1",
              "rdkit_version": rdBase.rdkitVersion, "fixtures": hashes,
              "artifact": artifact,
              "accounting": {"input": len(cases), "outcomes": counts}, "rows": rows,
              "limits": ["83 pinned cases only", "Raw embeddings are not required to match",
                         "No yield, selectivity or broad SMIRKS parity claim"]}
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(f"{args.output}: sha256={hashlib.sha256(rendered.encode()).hexdigest()}")
    else:
        print(rendered)
    if counts != EXPECTED:
        raise ValueError(f"checked Python outcome accounting differs: {counts} != {EXPECTED}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
