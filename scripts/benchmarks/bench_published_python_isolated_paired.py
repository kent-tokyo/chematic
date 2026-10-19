#!/usr/bin/env python3
"""Paired fresh-process Python benchmark with separate time and peak-RSS axes.

Run two exact published wheels (or a wheel and RDKit) on one exposed corpus.
Only HBA and RDKit-compatible Morgan are supported. Each arm is checked
against its archived per-row output; timing and memory remain distinct claims.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import platform
import random
import resource
import statistics
import subprocess
import sys
import time
from pathlib import Path

OP_NAMES = {"hba": "hba", "morgan": "morgan_r2_2048(rdkit-compatible)"}
MODES = ("parse_inclusive", "prepared_first_use", "precomputed")
ARTIFACT_SHA256 = {
    ("chematic", "1.0.29"): "a9ceb3c685e6cf4abbbdb1ae1cdbaa76b9ba13c8015d53750ba617cba3f7c1d8",
    ("chematic", "1.0.30"): "f6b1b22dd898cf3f100001156d9504b3753cb6b392143cbc0e6f3be78ec723d5",
    ("rdkit", "2026.03.6"): "e16c467cb254a223e59a0cf81358c6b39da15a99d2909170d693e95778fddb41",
}
REFERENCE_SHA256 = {
    ("chematic", "1.0.29"): "13a122220f7add4c22f347684db953b53cd47a56c222c020d7dcf61c102c5be7",
    ("chematic", "1.0.30"): "a63bc20a8edce597904d9fa1d54cc21822af906b55282bd964b6747ce13feb9f",
    ("rdkit", "2026.03.6"): "a63bc20a8edce597904d9fa1d54cc21822af906b55282bd964b6747ce13feb9f",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(values: list[str]) -> str:
    return hashlib.sha256(("\n".join(values) + "\n").encode("ascii")).hexdigest()


def reference_digest(path: Path, operation: str, count: int) -> str:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            if row["op"] != OP_NAMES[operation]:
                continue
            values = []
            for index, output in enumerate(row["rows"][:count]):
                if output.get("input_index") != index or "error" in output:
                    raise ValueError(f"reference {path}: invalid row {index}")
                value = output["value"]
                values.append(str(value) if operation == "hba" else value["bytes_hex"])
            if len(values) != count:
                raise ValueError(f"reference {path}: expected {count} rows")
            return digest(values)
    raise ValueError(f"reference {path}: missing {OP_NAMES[operation]}")


def ci(ratios: list[float]) -> dict:
    if len(ratios) < 20 or any(r <= 0 for r in ratios):
        raise ValueError("20 positive paired observations required")
    rng = random.Random(0)
    n = len(ratios)
    bootstrap = sorted(statistics.median(ratios[rng.randrange(n)] for _ in range(n)) for _ in range(4000))
    return {"median_a_over_b": statistics.median(ratios),
            "bootstrap_95pct": [bootstrap[99], bootstrap[3899]],
            "resamples": 4000, "seed": 0}


def worker(args: argparse.Namespace) -> None:
    rows = [line.split()[0] for line in args.corpus.read_text(encoding="utf-8").splitlines() if line.strip()][:args.limit]
    if len(rows) != args.limit:
        raise ValueError("corpus shorter than requested limit")

    if args.library == "chematic":
        import chematic
        if chematic.__version__ != args.expected_version:
            raise ValueError(f"loaded chematic {chematic.__version__}, expected {args.expected_version}")
        parse = chematic.from_smiles
        if args.operation == "hba":
            evaluate = lambda mol: mol.rdkit_hba
            encode = str
        else:
            evaluate = lambda mol: mol.rdkit_ecfp4()
            encode = lambda value: bytes(value).hex()
        runtime_version = chematic.__version__
    else:
        from rdkit import Chem, DataStructs, RDLogger
        from rdkit.Chem import rdMolDescriptors, rdFingerprintGenerator
        RDLogger.DisableLog("rdApp.*")
        runtime_version = Chem.rdBase.rdkitVersion
        if runtime_version != args.expected_version:
            raise ValueError(f"loaded RDKit {runtime_version}, expected {args.expected_version}")
        parse = Chem.MolFromSmiles
        if args.operation == "hba":
            evaluate = rdMolDescriptors.CalcNumHBA
            encode = str
        else:
            generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
            evaluate = generator.GetFingerprint
            encode = lambda value: DataStructs.BitVectToBinaryText(value).hex()

    # Warm a disjoint molecule, never one of the timed input objects.
    warm = parse(rows[0])
    if warm is None:
        raise ValueError("warmup parse failed")
    evaluate(warm)
    setup_start = time.perf_counter_ns()
    molecules = None if args.mode == "parse_inclusive" else [parse(smiles) for smiles in rows]
    if molecules is not None and any(mol is None for mol in molecules):
        raise ValueError("prepared input parse failed")
    if args.mode == "precomputed":
        [evaluate(mol) for mol in molecules]
    setup_ns = time.perf_counter_ns() - setup_start
    rss_after_setup = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if args.mode == "parse_inclusive":
        def call_all():
            return [evaluate(parse(smiles)) for smiles in rows]
    else:
        def call_all():
            return [evaluate(mol) for mol in molecules]

    started = time.perf_counter_ns()
    values = call_all()
    operation_ns = time.perf_counter_ns() - started
    encoded = [encode(value) for value in values]
    highwater = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    scale = 1 if sys.platform == "darwin" else 1024
    print(json.dumps({"library": args.library, "version": runtime_version,
                      "operation": args.operation, "mode": args.mode,
                      "input_count": len(rows), "setup_ns": setup_ns,
                      "operation_ns": operation_ns, "output_sha256": digest(encoded),
                      "peak_process_rss_bytes": highwater * scale,
                      "rss_after_setup_highwater_bytes": rss_after_setup * scale,
                      "platform": platform.platform()}, sort_keys=True))


def parent(args: argparse.Namespace) -> None:
    if args.blocks < 20 or args.limit < 1:
        raise ValueError("at least 20 paired blocks and a positive limit required")
    arms = {}
    for label in ("a", "b"):
        python = getattr(args, f"python_{label}")
        artifact = getattr(args, f"artifact_{label}")
        reference = getattr(args, f"reference_{label}")
        library = getattr(args, f"library_{label}")
        version = getattr(args, f"version_{label}")
        identity = (library, version)
        if sha256(artifact) != ARTIFACT_SHA256.get(identity) or \
                sha256(reference) != REFERENCE_SHA256.get(identity):
            raise ValueError(f"{label}: unpinned artifact or output reference for {identity}")
        expected_digest = reference_digest(reference, args.operation, args.limit)
        arms[label] = {"python": str(python), "library": library,
                       "version": version, "artifact_sha256": sha256(artifact),
                       "reference_archive_sha256": sha256(reference),
                       "expected_output_sha256": expected_digest}
    pairs = []
    for block in range(args.blocks):
        order = ("a", "b") if block % 2 == 0 else ("b", "a")
        result = {"block": block, "order": list(order)}
        for label in order:
            arm = arms[label]
            command = [arm["python"], __file__, "--worker", "--corpus", str(args.corpus),
                       "--limit", str(args.limit), "--operation", args.operation,
                       "--mode", args.mode, "--library", arm["library"],
                       "--expected-version", arm["version"]]
            completed = subprocess.run(command, check=True, capture_output=True, text=True)
            value = json.loads(completed.stdout)
            if value["output_sha256"] != arm["expected_output_sha256"]:
                raise ValueError(f"block {block} {label}: published-output mismatch")
            result[label] = value
        pairs.append(result)
        print(f"block {block + 1}/{args.blocks}: {result['a']['operation_ns']} / {result['b']['operation_ns']} ns", flush=True)
    speed = ci([pair["a"]["operation_ns"] / pair["b"]["operation_ns"] for pair in pairs])
    memory = ci([pair["a"]["peak_process_rss_bytes"] / pair["b"]["peak_process_rss_bytes"] for pair in pairs])
    outputs_equal = arms["a"]["expected_output_sha256"] == arms["b"]["expected_output_sha256"]
    perception_mismatch = (args.mode == "prepared_first_use" and
                           arms["a"]["library"] != arms["b"]["library"])
    if not outputs_equal:
        speed_outcome = "not_counted_output_difference"
    elif perception_mismatch:
        speed_outcome = "not_counted_eager_lazy_perception"
    elif speed["bootstrap_95pct"][0] > 1:
        speed_outcome = "arm_b_faster_on_declared_lane"
    elif speed["bootstrap_95pct"][1] < 1:
        speed_outcome = "arm_a_faster_on_declared_lane"
    else:
        speed_outcome = "inconclusive_interval_crosses_parity"
    report = {"schema": "published-python-isolated-paired-v1",
              "artifacts": arms, "corpus": {"sha256": sha256(args.corpus), "limit": args.limit},
              "operation": args.operation, "mode": args.mode,
              "protocol": {"blocks": args.blocks, "order": "AB/BA alternating fresh interpreters",
                           "timing": "only operation calls; parse included only in parse_inclusive",
                           "memory": "whole-process high-water RSS, including imports, preparation and result objects; not operation allocation",
                           "warmup": "one separate molecule per arm"},
              "outputs_equal_across_arms": outputs_equal,
              "speed_a_over_b": speed, "speed_outcome": speed_outcome,
              "memory_a_over_b": memory,
              "memory_interpretation": "whole-process peak only; no per-operation allocation or cross-binding memory win",
              "pairs": pairs}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"speed_a_over_b": speed, "memory_a_over_b": memory,
                      "outputs_equal": report["outputs_equal_across_arms"]}))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=5000)
    parser.add_argument("--operation", choices=OP_NAMES, required=True)
    parser.add_argument("--mode", choices=MODES, required=True)
    parser.add_argument("--library", choices=("chematic", "rdkit"))
    parser.add_argument("--expected-version")
    for label in ("a", "b"):
        parser.add_argument(f"--python-{label}", type=Path)
        parser.add_argument(f"--artifact-{label}", type=Path)
        parser.add_argument(f"--reference-{label}", type=Path)
        parser.add_argument(f"--library-{label}", choices=("chematic", "rdkit"))
        parser.add_argument(f"--version-{label}")
    parser.add_argument("--blocks", type=int, default=20)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.worker:
        if args.library is None or args.expected_version is None:
            parser.error("worker requires --library and --expected-version")
        worker(args)
    else:
        if any(getattr(args, key) is None for label in ("a", "b") for key in
               (f"python_{label}", f"artifact_{label}", f"reference_{label}",
                f"library_{label}", f"version_{label}")) or args.output is None:
            parser.error("both pinned arms and --output are required")
        parent(args)


if __name__ == "__main__":
    main()
