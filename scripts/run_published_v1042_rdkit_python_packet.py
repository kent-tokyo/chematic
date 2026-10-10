#!/usr/bin/env python3
"""Reproduce the v1.0.42 Python/RDKit packet from published wheels.

The parent process performs the output-equivalence gate and launches one fresh
child per engine and block.  Child timings therefore do not load the competing
binding in the same interpreter.  A speed interval is emitted only when the
corresponding 10k output gate passes.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import random
import statistics
import subprocess
import sys
import time
from pathlib import Path

EXPECTED_CHEMATIC = "1.0.42"
EXPECTED_RDKIT = "2026.9.1"
OPS = ("tpsa", "labute_asa", "num_rings", "morgan2_chiral")


def packed_on_bits(raw: bytes) -> tuple[int, ...]:
    return tuple(
        index
        for index in range(len(raw) * 8)
        if raw[index // 8] >> (index % 8) & 1
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rows_from(path: Path, limit: int | None = None) -> list[str]:
    rows = [line.split()[0] for line in path.read_text().splitlines() if line.strip()]
    return rows if limit is None else rows[:limit]


def engine(engine_name: str):
    if engine_name == "chematic":
        import chematic

        parse = chematic.from_smiles
        operations = {
            "tpsa": lambda mol: mol.rdkit_tpsa,
            "labute_asa": lambda mol: mol.labute_asa,
            "num_rings": lambda mol: mol.num_rings,
            "morgan2_chiral": lambda mol: mol.rdkit_ecfp_config(
                2, 2048, include_chirality=True
            ),
        }
        outputs = operations | {
            "morgan2_chiral": lambda mol: packed_on_bits(operations["morgan2_chiral"](mol))
        }
        return parse, operations, outputs, chematic.__version__

    from rdkit import Chem, rdBase
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.Chem import rdMolDescriptors as descriptors

    morgan = rdFingerprintGenerator.GetMorganGenerator(
        radius=2, fpSize=2048, includeChirality=True
    )
    operations = {
        "tpsa": descriptors.CalcTPSA,
        "labute_asa": descriptors.CalcLabuteASA,
        "num_rings": descriptors.CalcNumRings,
        "morgan2_chiral": morgan.GetFingerprint,
    }
    outputs = operations | {
        "morgan2_chiral": lambda mol: tuple(operations["morgan2_chiral"](mol).GetOnBits())
    }
    return Chem.MolFromSmiles, operations, outputs, rdBase.rdkitVersion


def consume(value) -> int:
    if isinstance(value, (int, float)):
        return hash(value)
    if hasattr(value, "GetNumBits"):
        return int(value.GetNumBits())
    return len(value)


def measure_child(engine_name: str, corpus: Path, limit: int, block: int) -> dict:
    parse, operations, _, version = engine(engine_name)
    rows = rows_from(corpus, limit)
    ordered_ops = list(OPS[block % len(OPS) :] + OPS[: block % len(OPS)])
    result: dict[str, object] = {"engine": engine_name, "version": version}

    started = time.perf_counter()
    parsed = [parse(smiles) for smiles in rows]
    if any(mol is None for mol in parsed):
        raise RuntimeError(f"{engine_name} rejected a timing row")
    result["parse_seconds"] = time.perf_counter() - started

    for op in ordered_ops:
        fn = operations[op]
        prepared = [parse(smiles) for smiles in rows]
        if any(mol is None for mol in prepared):
            raise RuntimeError(f"{engine_name} rejected a prepared {op} row")
        checksum = 0
        started = time.perf_counter()
        for mol in prepared:
            value = fn(mol)
            checksum ^= consume(value)
        result[f"{op}_prepared_seconds"] = time.perf_counter() - started
        result[f"{op}_prepared_checksum"] = checksum

        checksum = 0
        started = time.perf_counter()
        for smiles in rows:
            mol = parse(smiles)
            if mol is None:
                raise RuntimeError(f"{engine_name} rejected a pipeline {op} row")
            value = fn(mol)
            checksum ^= consume(value)
        result[f"{op}_pipeline_seconds"] = time.perf_counter() - started
        result[f"{op}_pipeline_checksum"] = checksum
    return result


def paired_ratio_ci(numerator: list[float], denominator: list[float], seed: int) -> dict:
    ratios = [a / b for a, b in zip(numerator, denominator)]
    rng = random.Random(seed)
    estimates = [
        statistics.median(ratios[rng.randrange(len(ratios))] for _ in ratios)
        for _ in range(10_000)
    ]
    estimates.sort()
    return {
        "median": statistics.median(ratios),
        "ci95": [estimates[249], estimates[9749]],
        "winning_blocks": sum(ratio > 1.0 for ratio in ratios),
        "blocks": len(ratios),
    }


def compare_outputs(corpus: Path) -> dict:
    rows = rows_from(corpus)
    ch_parse, _, ch_ops, _ = engine("chematic")
    rd_parse, _, rd_ops, _ = engine("rdkit")
    report = {
        op: {"compared": 0, "matches": 0, "mismatches": [], "errors": []}
        for op in OPS
    }
    for index, smiles in enumerate(rows):
        try:
            ch_mol = ch_parse(smiles)
            rd_mol = rd_parse(smiles)
            if ch_mol is None or rd_mol is None:
                raise ValueError("one parser returned no molecule")
        except Exception as exc:  # pragma: no cover - retained in raw report
            for op in OPS:
                report[op]["errors"].append({"row": index, "error": str(exc)})
            continue
        for op in OPS:
            try:
                actual = ch_ops[op](ch_mol)
                expected = rd_ops[op](rd_mol)
                if op == "labute_asa":
                    matches = abs(actual - expected) <= 1e-9
                else:
                    matches = actual == expected
                report[op]["compared"] += 1
                report[op]["matches"] += int(matches)
                if not matches and len(report[op]["mismatches"]) < 50:
                    report[op]["mismatches"].append(
                        {"row": index, "smiles": smiles, "chematic": actual, "rdkit": expected}
                    )
            except Exception as exc:  # pragma: no cover - retained in raw report
                report[op]["errors"].append({"row": index, "error": str(exc)})
    for op, values in report.items():
        values["equivalent"] = (
            values["compared"] == len(rows)
            and values["matches"] == len(rows)
            and not values["errors"]
        )
        values["mismatch_count"] = values["compared"] - values["matches"]
        values["error_count"] = len(values["errors"])
    return report


def artifact_records(directory: Path) -> list[dict]:
    records = []
    for path in sorted(item for item in directory.iterdir() if item.is_file()):
        records.append(
            {"file": path.name, "bytes": path.stat().st_size, "sha256": sha256_file(path)}
        )
    return records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--artifact-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--blocks", type=int, default=21)
    parser.add_argument("--limit", type=int, default=1000)
    parser.add_argument("--child", choices=("chematic", "rdkit"))
    parser.add_argument("--block", type=int, default=0)
    args = parser.parse_args()
    if args.child:
        print(json.dumps(measure_child(args.child, args.corpus, args.limit, args.block)))
        return 0

    ch_version = importlib.metadata.version("chematic")
    rd_version = importlib.metadata.version("rdkit")
    if ch_version != EXPECTED_CHEMATIC or rd_version != EXPECTED_RDKIT:
        raise SystemExit(
            f"expected chematic {EXPECTED_CHEMATIC} / rdkit {EXPECTED_RDKIT}, "
            f"got {ch_version} / {rd_version}"
        )
    samples = {"chematic": [], "rdkit": []}
    for block in range(args.blocks):
        order = ("chematic", "rdkit") if block % 2 == 0 else ("rdkit", "chematic")
        for engine_name in order:
            command = [
                sys.executable,
                __file__,
                "--corpus",
                str(args.corpus),
                "--artifact-dir",
                str(args.artifact_dir),
                "--output",
                str(args.output),
                "--limit",
                str(args.limit),
                "--child",
                engine_name,
                "--block",
                str(block),
            ]
            completed = subprocess.run(command, check=True, capture_output=True, text=True)
            samples[engine_name].append(json.loads(completed.stdout))
        print(f"published Python block {block + 1}/{args.blocks}", file=sys.stderr)

    accuracy = compare_outputs(args.corpus)
    speed = {}
    for op_index, op in enumerate(OPS):
        speed[op] = {"eligible": accuracy[op]["equivalent"]}
        for lane in ("prepared", "pipeline"):
            key = f"{op}_{lane}_seconds"
            speed[op][lane] = (
                paired_ratio_ci(
                    [sample[key] for sample in samples["rdkit"]],
                    [sample[key] for sample in samples["chematic"]],
                    seed=0x1042 + op_index * 2 + (lane == "pipeline"),
                )
                if accuracy[op]["equivalent"]
                else None
            )

    corpus_bytes = args.corpus.read_bytes()
    report = {
        "schema": "published-v1.0.42-rdkit-packet/v1",
        "binding": "python",
        "versions": {"chematic": ch_version, "rdkit": rd_version},
        "artifacts": artifact_records(args.artifact_dir),
        "corpus": {
            "path": args.corpus.name,
            "rows": len(rows_from(args.corpus)),
            "sha256": hashlib.sha256(corpus_bytes).hexdigest(),
            "timing_rows": args.limit,
        },
        "environment": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "processor": platform.processor(),
            "python": platform.python_version(),
            "cpus": os.cpu_count(),
            "runner_os": os.environ.get("RUNNER_OS"),
            "runner_arch": os.environ.get("RUNNER_ARCH"),
            "runner_name": os.environ.get("RUNNER_NAME"),
        },
        "method": {
            "blocks": args.blocks,
            "order": "alternating fresh child process per engine and block",
            "prepared": "parse outside timing, then first operation call",
            "pipeline": "parse plus operation inside timing",
            "confidence_interval": "paired bootstrap of per-block RDKit/chematic ratios, 10000 resamples",
            "speed_claim_policy": "withhold both intervals unless the 10000-row output gate passes",
        },
        "accuracy": accuracy,
        "raw_samples_seconds": samples,
        "speedup_rdkit_over_chematic": speed,
    }
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
