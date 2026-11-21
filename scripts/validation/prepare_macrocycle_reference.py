#!/usr/bin/env python3
"""Prepare a graph-audited reference SDF from a pinned macrocycle archive.

This checks input identity and coordinate finiteness only. It does not generate
or score conformers. Original SDF content is preserved except for the title,
which is replaced with the stable case ID required by the benchmark scorer.
"""
import argparse
from collections import Counter
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import re
import tarfile

from rdkit import Chem, RDLogger, rdBase


def sha(data):
    return hashlib.sha256(data).hexdigest()


def archive_sha(path):
    h = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_cases(archive, source_commit, case_ids):
    expected_root = f"MacroCycle-Benchmark-{source_commit}"
    files = {}
    with tarfile.open(archive, mode="r|gz") as source:
        for member in source:
            parts = Path(member.name).parts
            if len(parts) != 4 or parts[0] != expected_root or parts[1] != "cases" or parts[3] != "ligand_unprotonated.sdf" or parts[2] not in case_ids:
                continue
            case = parts[2]
            if case in files or not member.isfile() or member.size > 5_000_000:
                raise ValueError(f"duplicate, non-file or oversized reference: {case}")
            stream = source.extractfile(member)
            if stream is None:
                raise ValueError(f"reference could not be read: {case}")
            files[case] = stream.read()
    if set(files) != case_ids:
        raise ValueError(f"archive missing {len(case_ids - set(files))} indexed references")
    return files


def title_replaced(raw, name):
    lines = raw.splitlines(keepends=True)
    if not lines or not lines[0].endswith((b"\n", b"\r")) or raw.count(b"$$$$") != 1 or not raw.rstrip().endswith(b"$$$$"):
        raise ValueError("reference is not a single terminated SDF record")
    newline = b"\r\n" if lines[0].endswith(b"\r\n") else b"\n"
    result = name.encode("ascii") + newline + b"".join(lines[1:])
    if not result.endswith(b"\n"):
        result += b"\n"
    return result


def validate_reference(raw, indexed_smiles):
    molecules = list(Chem.ForwardSDMolSupplier(io.BytesIO(raw), removeHs=False))
    if len(molecules) != 1 or molecules[0] is None:
        raise ValueError("SDF parse failed or contains multiple records")
    mol = molecules[0]
    reference = Chem.MolFromSmiles(indexed_smiles)
    if reference is None:
        raise ValueError("indexed SMILES no longer parses")
    actual = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)
    expected = Chem.MolToSmiles(Chem.RemoveHs(reference), isomericSmiles=True)
    if actual != expected:
        raise ValueError("SDF and indexed isomeric graphs differ")
    if mol.GetNumConformers() != 1:
        raise ValueError("reference must have exactly one coordinate set")
    conf = mol.GetConformer()
    coords = [tuple(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    if not coords or not all(math.isfinite(x) for point in coords for x in point) or len(set(coords)) < 2:
        raise ValueError("reference coordinates are nonfinite or degenerate")
    return {"atoms": mol.GetNumAtoms(), "heavy_atoms": mol.GetNumHeavyAtoms(),
            "canonical_graph_sha256": sha(actual.encode())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("index", type=Path)
    parser.add_argument("split", type=Path)
    parser.add_argument("audit", type=Path)
    parser.add_argument("--sdf", type=Path, required=True)
    parser.add_argument("--prepared-split", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    index_bytes, split_bytes, audit_bytes = args.index.read_bytes(), args.split.read_bytes(), args.audit.read_bytes()
    split, audit = json.loads(split_bytes), json.loads(audit_bytes)
    if rdBase.rdkitVersion != split["rdkit_version"] or sha(index_bytes) != split["index_sha256"]:
        raise ValueError("index hash or RDKit version differs from frozen split")
    if (audit["index_sha256"] != sha(index_bytes) or audit["split_sha256"] != sha(split_bytes)
            or audit["invalid_exclusion_count"] != len(split["invalid_exclusions"]) or audit["overlaps"]):
        raise ValueError("invalid-exclusion audit is not clean or does not match")
    with args.index.open(encoding="utf-8", newline="") as source:
        index = list(csv.DictReader(source))
    if len(index) != len(split["rows"]):
        raise ValueError("index and frozen split have different record counts")
    case_ids = {r["case_id"] for r in index}
    if len(case_ids) != len(index) or any(not re.fullmatch(r"[A-Za-z0-9_.-]+", x) for x in case_ids):
        raise ValueError("duplicate or unsafe case IDs")
    RDLogger.DisableLog("rdApp.*")
    files = read_cases(args.archive, split["source_commit"], case_ids)
    records, pieces, failures = [], [], []
    for row in split["rows"]:
        source_index = row["index"]
        reference = index[source_index]
        case = reference["case_id"]
        if row["name"] != case or row["input_smiles_sha256"] != sha(reference["unprotonated_isomeric_smiles"].encode()):
            raise ValueError("index identity differs from frozen split")
        raw = files[case]
        try:
            details = validate_reference(raw, reference["unprotonated_isomeric_smiles"])
            normalized = title_replaced(raw, case)
        except Exception as exc:
            failures.append({"source_index": source_index, "name": case,
                             "split": row.get("split", row.get("status")), "error": repr(exc),
                             "source_sdf_sha256": sha(raw)})
            continue
        prepared = dict(row)
        prepared.update(index=len(records), source_index=source_index,
                        source_sdf_sha256=sha(raw), **details)
        records.append(prepared)
        pieces.append(normalized)
    if not records:
        raise ValueError("no valid reference records")
    sdf_bytes = b"".join(pieces)
    args.sdf.write_bytes(sdf_bytes)
    prepared_split = dict(split)
    prepared_split.update(source_sha256=sha(sdf_bytes), original_index_sha256=sha(index_bytes),
                          source_split_sha256=sha(split_bytes), source_audit_sha256=sha(audit_bytes),
                          archive_sha256=archive_sha(args.archive), rows=records,
                          counts=dict(Counter(r.get("split", r.get("status")) for r in records)))
    split_output = (json.dumps(prepared_split, indent=2) + "\n").encode()
    args.prepared_split.write_bytes(split_output)
    report = {"source_commit": split["source_commit"], "archive_sha256": prepared_split["archive_sha256"],
              "index_sha256": sha(index_bytes), "original_split_sha256": sha(split_bytes),
              "sdf_sha256": sha(sdf_bytes), "prepared_split_sha256": sha(split_output),
              "source_cases": len(index), "prepared_records": len(records),
              "prepared_counts": prepared_split["counts"], "input_failures": failures,
              "failure_counts": dict(Counter(f["split"] for f in failures))}
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: report[k] for k in ("source_cases", "prepared_records", "prepared_counts", "failure_counts")}, sort_keys=True))


if __name__ == "__main__":
    main()
