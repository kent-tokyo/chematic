#!/usr/bin/env python3
"""V3000 native MOL APIs; retain every rejected rewrite and suppress its rank."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

from bench_openbabel_sdf_throughput import (
    ROOT, FIELDS, eligible, observe, observe_rows, semantic_check, sha, timing_row,
    paired_median_ratio_ci,
)

MODES = ("parse", "write", "roundtrip")


def rewrite_check(reference: list[dict], payload: bytes) -> dict:
    rows = observe_rows(payload, allow_rejected=True)
    result = {"passed": False, "expected_rows": len(reference), "observed_rows": len(rows),
              "record_count_mismatch": len(rows) != len(reference),
              "rejected_rows": [i for i, row in enumerate(rows) if row is None]}
    if result["record_count_mismatch"]:
        return result
    counts = dict.fromkeys(FIELDS, 0)
    examples = {field: [] for field in FIELDS}
    for index, (expected, row) in enumerate(zip(reference, rows, strict=True)):
        if row is None:
            continue  # Its position remains in rejected_rows; it never counts as a pass.
        check = semantic_check([expected], [row])
        for field, count in check["mismatch_counts"].items():
            counts[field] = counts.get(field, 0) + count
            examples.setdefault(field, [])
            if count and len(examples[field]) < 8:
                examples[field].append(index)
    result.update(mismatch_counts=counts, mismatch_examples=examples,
                  accepted_rows=len(rows) - len(result["rejected_rows"]))
    result["passed"] = not result["rejected_rows"] and not any(counts.values())
    return result


def inspect(binaries: dict[str, Path], corpus: Path, expected: int, artifacts: Path) -> dict:
    from rdkit import Chem

    reference = observe(corpus.read_bytes())
    if len(reference) != expected:
        raise ValueError("reference row count differs; no filtering allowed")
    artifacts.mkdir(parents=True, exist_ok=True)
    semantics, snapshots, rewritten = {}, {}, {}
    for name, binary in binaries.items():
        raw = subprocess.run([str(binary), str(corpus), "dump", "1", str(expected)],
                             capture_output=True, check=True, timeout=180).stdout
        (artifacts / (name + ".molblocks")).write_bytes(raw)
        rewritten[name] = sha(raw)
        semantics[name] = rewrite_check(reference, raw)
        mode = "smiles" if name == "openbabel" else "snapshot"
        raw = subprocess.run([str(binary), str(corpus), mode, "1", str(expected)],
                             capture_output=True, check=True, timeout=180).stdout
        (artifacts / (name + "." + mode)).write_bytes(raw)
        if name == "openbabel":
            strings = [line.split()[0] for line in raw.decode().splitlines()]
        else:
            snapshots[name] = raw
            strings = [row["canonical_smiles"] for row in json.loads(raw)]
        if len(strings) != expected:
            raise ValueError(f"{name}: internal snapshot row count differs")
        mismatches, rejected = [], []
        for index, (row, string) in enumerate(zip(reference, strings, strict=True)):
            mol = Chem.MolFromSmiles(string)
            if mol is None:
                rejected.append(index)
            elif Chem.MolToSmiles(mol) != row["graph"]:
                mismatches.append(index)
        semantics[name + "_parse"] = {"passed": not mismatches and not rejected,
            "rows": len(strings), "mismatches": len(mismatches), "mismatch_rows": mismatches,
            "rejected_rows": rejected}
    return {"semantics": semantics,
            "source_parse_snapshots_equal": snapshots["baseline"] == snapshots["candidate"],
            "snapshot_sha256": {name: sha(raw) for name, raw in snapshots.items()},
            "rewrite_sha256": rewritten}


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
        parser.error("Darwin, 21+ blocks and positive workload sizes required")
    binaries = {name: getattr(args, name) for name in ("baseline", "candidate", "openbabel")}
    source_bytes = args.source_manifest.read_bytes()
    source = json.loads(source_bytes)
    binary_hashes = {name: sha(path.read_bytes()) for name, path in binaries.items()}
    corpus_bytes = args.corpus.read_bytes()

    def integrity():
        if args.source_manifest.read_bytes() != source_bytes or args.corpus.read_bytes() != corpus_bytes:
            raise ValueError("manifest/corpus changed")
        for path, expected in source["source_files_sha256"].items():
            if sha((ROOT / path).read_bytes()) != expected:
                raise ValueError(f"source changed: {path}")
        for path, expected in source["native_libraries_sha256"].items():
            if sha(Path(path).read_bytes()) != expected:
                raise ValueError(f"Open Babel library changed: {path}")
        for name, path in binaries.items():
            if sha(path.read_bytes()) != binary_hashes[name]:
                raise ValueError(f"binary changed: {name}")

    integrity()
    checks = inspect(binaries, args.corpus, args.records, args.output.with_suffix(".artifacts"))
    report = {"schema": "openbabel-v3000-block-throughput/v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_manifest": source, "source_manifest_sha256": sha(source_bytes),
        "binaries": {name: {"path": str(path), "sha256": binary_hashes[name]} for name, path in binaries.items()},
        "corpus": {"path": str(args.corpus), "sha256": sha(corpus_bytes),
                   "bytes": len(corpus_bytes), "records": args.records},
        "host": {"platform": platform.platform(), "python": sys.version, "cpus": os.cpu_count()},
        "rdkit_observer_version": Chem.rdBase.rdkitVersion,
        "openbabel_version": subprocess.run(["obabel", "-V"], check=True, capture_output=True, text=True).stdout.strip(),
        "protocol": {"blocks": args.blocks, "repeats": args.repeats, "modes": MODES,
            "order": "rotating three-arm fresh processes, sequential",
            "timing": "native MOL API batch, delimiter splitting/loading/plugin discovery excluded for both engines",
            "parse": "borrowed pre-split MOL blocks; full parsed batch retained; default diagnostic Rust reader",
            "write": "prepared batch, one untimed warm write; allocating native MOL strings in both engines",
            "roundtrip": "fresh full-batch parse plus MOL writing inside one timer",
            "rss": "Darwin whole-process high-water RSS including loading, framing, prepared batch, warmup and plugins",
            "ranking": "equal internal graphs for parse; equal accepted rewrites for write/roundtrip; CI lower >1 and 21/21 faster"},
        **checks, "blocks": [], "summary": {}, "complete": False}

    def save():
        temporary = args.output.with_suffix(".partial")
        temporary.write_text(json.dumps(report, indent=2) + "\n")
        temporary.replace(args.output)

    save()
    names = list(binaries)
    for block in range(args.blocks):
        order = names[block % 3:] + names[:block % 3]
        row = {"order": order, "load_average_start": list(os.getloadavg())}
        for name in order:
            row[name] = {mode: timing_row(binaries[name], args.corpus, mode, args.repeats, args.records)
                         for mode in MODES}
        report["blocks"].append(row)
        save()
        print(f"block {block + 1}/{args.blocks}", flush=True)
    for mode in MODES:
        report["summary"][mode] = {}
        for index, name in enumerate(("baseline", "openbabel")):
            allowed = eligible(report["semantics"], report["source_parse_snapshots_equal"], mode, name)
            interval = paired_median_ratio_ci(
                [row[name][mode]["elapsed_ns"] for row in report["blocks"]],
                [row["candidate"][mode]["elapsed_ns"] for row in report["blocks"]], seed=index) if allowed else None
            report["summary"][mode][name] = {"speedup": interval,
                "semantic_prerequisite": allowed,
                "strict_speed_win": bool(interval and interval["ci95"][0] > 1
                                         and interval["winning_blocks"] == args.blocks)}
    integrity()
    report["complete"] = True
    save()


if __name__ == "__main__":
    main()
