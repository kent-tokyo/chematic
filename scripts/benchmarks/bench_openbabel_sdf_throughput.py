#!/usr/bin/env python3
"""Three-arm SDF batch throughput with semantic prerequisites and raw blocks.

Input buffers are preloaded; executable startup and file loading are outside
native timers. Parsers retain the whole batch. RSS is whole-process peak,
including loading, preparation and warmup, not a streaming-memory guarantee.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import io
import json
import math
import os
from pathlib import Path
import platform
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "validation/cosmolkit_comparison"))
from bench_corpus import paired_median_ratio_ci  # noqa: E402

MODES = ("parse", "write", "roundtrip")
FIELDS = ("graph", "atoms", "bonds", "coordinates", "name", "properties", "explicit_unknown_stereo")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def observe(payload: bytes) -> list[dict]:
    return observe_rows(payload)


def observe_rows(payload: bytes, *, allow_rejected: bool = False) -> list:
    """Retain rejected row positions only when a caller explicitly asks."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    rows = []
    for molecule in Chem.ForwardSDMolSupplier(io.BytesIO(payload), sanitize=True, removeHs=False, strictParsing=True):
        if molecule is None:
            if allow_rejected:
                rows.append(None)
                continue
            raise ValueError("SDF observer rejected a row; no filtering allowed")
        unknown = sorted(tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx())))
                         for b in molecule.GetBonds() if b.GetStereo() == Chem.BondStereo.STEREOANY)
        rows.append({
            "graph": Chem.MolToSmiles(molecule),
            "atoms": [(a.GetAtomicNum(), a.GetIsotope(), a.GetFormalCharge(), a.GetNumRadicalElectrons(),
                       a.GetTotalNumHs()) for a in molecule.GetAtoms()],
            "bonds": sorted((min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                             b.GetBondTypeAsDouble()) for b in molecule.GetBonds()),
            "coordinates": [list(p) for p in molecule.GetConformer().GetPositions()],
            "name": molecule.GetProp("_Name") if molecule.HasProp("_Name") else "",
            "properties": molecule.GetPropsAsDict(includePrivate=False, includeComputed=False, autoConvertStrings=False),
            "explicit_unknown_stereo": unknown,
        })
    return rows


def semantic_check(reference: list[dict], output: list[dict]) -> dict:
    mismatches = {key: [] for key in FIELDS}
    if len(reference) != len(output):
        return {"passed": False, "expected_rows": len(reference), "observed_rows": len(output),
                "record_count_mismatch": True}
    for index, (a, b) in enumerate(zip(reference, output, strict=True)):
        for field in FIELDS:
            if field == "coordinates":
                passed = len(a[field]) == len(b[field]) and all(
                    len(pa) == len(pb) == 3 and all(
                        math.isfinite(x) and math.isfinite(y) and abs(x - y) <= 1e-4
                        for x, y in zip(pa, pb))
                    for pa, pb in zip(a[field], b[field]))
            elif field == "explicit_unknown_stereo":
                # Marking additional unspecified bonds "either" is a safe
                # normalization; every explicit input marker must survive.
                passed = set(a[field]).issubset(set(b[field]))
            else:
                passed = a[field] == b[field]
            if not passed:
                mismatches[field].append(index)
    return {"passed": not any(mismatches.values()), "expected_rows": len(reference),
            "observed_rows": len(output), "record_count_mismatch": False,
            "mismatch_counts": {key: len(values) for key, values in mismatches.items()},
            "mismatch_examples": {key: values[:8] for key, values in mismatches.items()}}


def timing_row(binary: Path, corpus: Path, mode: str, repeats: int, records: int) -> dict:
    process = subprocess.run(["/usr/bin/time", "-l", str(binary), str(corpus), mode, str(repeats), str(records)],
                             check=True, capture_output=True, timeout=180)
    row = json.loads(process.stdout)
    value = row.get("elapsed_ns")
    if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
        raise ValueError("invalid native timing")
    if row.get("records") != repeats * records or row.get("repeats") != repeats or row.get("mode") != mode:
        raise ValueError("native workload accounting differs")
    match = re.search(rb"(\d+)\s+maximum resident set size", process.stderr)
    if match is None:
        raise ValueError("Darwin peak RSS unavailable")
    row["peak_process_rss_bytes"] = int(match.group(1))
    return row


def eligible(semantics: dict, snapshots_equal: bool, mode: str, competitor: str) -> bool:
    if mode == "parse":
        return ((competitor != "baseline" or snapshots_equal)
                and semantics["candidate_parse"]["passed"] and semantics[competitor + "_parse"]["passed"])
    return semantics["candidate"]["passed"] and semantics[competitor]["passed"]


def main():
    from rdkit import Chem

    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("baseline", "candidate", "openbabel", "corpus", "output", "source_manifest"):
        parser.add_argument("--" + name.replace("_", "-"), type=Path, required=True)
    parser.add_argument("--records", type=int, default=1000)
    parser.add_argument("--blocks", type=int, default=21)
    parser.add_argument("--repeats", type=int, default=5)
    args = parser.parse_args()
    if platform.system() != "Darwin" or args.blocks < 21 or args.repeats < 1 or args.records < 1:
        parser.error("Darwin, at least 21 blocks and positive workload sizes required")
    binaries = {name: getattr(args, name) for name in ("baseline", "candidate", "openbabel")}
    source_payload = args.source_manifest.read_bytes()
    source_manifest = json.loads(source_payload)

    def check_source_files():
        for path, expected in source_manifest["source_files_sha256"].items():
            if sha((ROOT / path).read_bytes()) != expected:
                raise ValueError(f"source differs from manifest: {path}")

    check_source_files()
    payload = args.corpus.read_bytes()
    reference = observe(payload)
    if len(reference) != args.records:
        raise ValueError("source record count differs")
    artifact_dir = args.output.with_suffix(".artifacts")
    artifact_dir.mkdir(parents=True, exist_ok=True)
    semantics = {}
    snapshots = {}
    rewrite_hashes = {}
    for name, binary in binaries.items():
        data = subprocess.run([str(binary), str(args.corpus), "dump", "1", str(args.records)],
                              check=True, capture_output=True, timeout=180).stdout
        (artifact_dir / (name + ".sdf")).write_bytes(data)
        rewrite_hashes[name] = sha(data)
        semantics[name] = semantic_check(reference, observe(data))
        if name != "openbabel":
            snapshot = subprocess.run([str(binary), str(args.corpus), "snapshot", "1", str(args.records)],
                                      check=True, capture_output=True, timeout=180).stdout
            snapshots[name] = snapshot
            rows = json.loads(snapshot)
            canonical = [Chem.MolToSmiles(Chem.MolFromSmiles(row["canonical_smiles"])) for row in rows]
        else:
            smiles = subprocess.run([str(binary), str(args.corpus), "smiles", "1", str(args.records)],
                                    check=True, capture_output=True, timeout=180).stdout.decode().splitlines()
            canonical = [Chem.MolToSmiles(Chem.MolFromSmiles(row.split()[0])) for row in smiles]
        bad = [i for i, (a, b) in enumerate(zip(reference, canonical, strict=True)) if a["graph"] != b]
        semantics[name + "_parse"] = {"passed": not bad, "rows": len(canonical), "mismatches": len(bad), "examples": bad[:8]}
    snapshots_equal = snapshots["baseline"] == snapshots["candidate"]
    report = {
        "schema": "openbabel-sdf-throughput/v1", "created_utc": datetime.now(timezone.utc).isoformat(),
        "corpus": {"path": str(args.corpus), "sha256": sha(payload), "bytes": len(payload), "records": args.records},
        "binaries": {name: {"path": str(path), "sha256": sha(path.read_bytes())} for name, path in binaries.items()},
        "source_manifest": source_manifest, "source_manifest_sha256": sha(source_payload),
        "host": {"platform": platform.platform(), "python": sys.version, "cpus": os.cpu_count()},
        "openbabel_version": subprocess.run(["obabel", "-V"], check=True, capture_output=True, text=True).stdout.strip(),
        "rdkit_observer_version": Chem.rdBase.rdkitVersion,
        "protocol": {"blocks": args.blocks, "repeats": args.repeats, "modes": MODES,
                     "order": "rotating three-arm fresh processes, sequential",
                     "timing": "native in-process loop; input loading and warmup excluded",
                     "parse": "borrowed in-memory input, full batch retained; default diagnostic Rust reader",
                     "write": "prepared batch, one untimed write, reusable Rust String / native OB WriteString",
                     "roundtrip": "fresh batch parse plus write inside one timer",
                     "rss": "Darwin whole-process high-water RSS, including input, prepared batch and warmup"},
        "semantics": semantics, "source_parse_snapshots_equal": snapshots_equal,
        "snapshot_sha256": {name: sha(data) for name, data in snapshots.items()}, "rewrite_sha256": rewrite_hashes,
        "blocks": [], "summary": {}, "complete": False,
    }

    def save():
        temporary = args.output.with_suffix(".partial")
        temporary.write_text(json.dumps(report, indent=2) + "\n")
        temporary.replace(args.output)

    save()
    names = list(binaries)
    for block in range(args.blocks):
        order = names[block % 3:] + names[:block % 3]
        values = {"order": order, "load_average_start": list(os.getloadavg())}
        for name in order:
            values[name] = {mode: timing_row(binaries[name], args.corpus, mode, args.repeats, args.records) for mode in MODES}
        report["blocks"].append(values)
        save()
        print(f"block {block + 1}/{args.blocks}", flush=True)
    for mode in MODES:
        report["summary"][mode] = {}
        for index, name in enumerate(("baseline", "openbabel")):
            allowed = eligible(semantics, snapshots_equal, mode, name)
            speed = paired_median_ratio_ci(
                [row[name][mode]["elapsed_ns"] for row in report["blocks"]],
                [row["candidate"][mode]["elapsed_ns"] for row in report["blocks"]], seed=index) if allowed else None
            report["summary"][mode][name] = {"speedup": speed, "semantic_prerequisite": allowed,
                "strict_speed_win": bool(speed and speed["ci95"][0] > 1 and speed["winning_blocks"] == args.blocks)}
    if sha(args.corpus.read_bytes()) != report["corpus"]["sha256"]:
        raise ValueError("corpus changed during measurement")
    for name, binary in binaries.items():
        if sha(binary.read_bytes()) != report["binaries"][name]["sha256"]:
            raise ValueError("binary changed during measurement")
    if sha(args.source_manifest.read_bytes()) != sha(source_payload):
        raise ValueError("source manifest changed during measurement")
    check_source_files()
    report["complete"] = True
    save()


if __name__ == "__main__":
    main()
