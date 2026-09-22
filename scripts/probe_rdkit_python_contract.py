#!/usr/bin/env python3
"""Probe the runtime contract and call boundary of one distributed RDKit wheel.

This is deliberately not a chemistry-kernel benchmark.  It records the
runtime behavior of the installed Python binding and keeps scalar and batch
API timings separate so a future Boost.Python/nanobind change cannot be
silently folded into chemistry performance.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import statistics
import sys
import tempfile
import time
from pathlib import Path
from typing import Callable


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": sha256(path)}


def load_smiles(path: Path, limit: int) -> list[str]:
    rows: list[str] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.startswith("{"):
            record = json.loads(stripped)
            if record.get("_manifest") is True:
                continue
            value = record.get("smiles")
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{path}:{line_number}: missing smiles")
            rows.append(value.strip())
        else:
            rows.append(stripped.split()[0])
        if len(rows) == limit:
            break
    if not rows:
        raise ValueError(f"no SMILES rows in {path}")
    return rows


def summarize_value(value: object) -> object:
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, (list, tuple)):
        return {"type": type(value).__name__, "length": len(value)}
    return {"type": f"{type(value).__module__}.{type(value).__qualname__}"}


def observe(call: Callable[[], object]) -> dict[str, object]:
    try:
        return {"outcome": "returned", "value": summarize_value(call())}
    except Exception as exc:  # the exception type is the evidence
        return {
            "outcome": "raised",
            "exception_type": f"{type(exc).__module__}.{type(exc).__qualname__}",
            "message": str(exc),
        }


def timed(call: Callable[[], int], repetitions: int) -> dict[str, object]:
    call()
    samples: list[int] = []
    result = 0
    for _ in range(repetitions):
        started = time.perf_counter_ns()
        result = call()
        samples.append(time.perf_counter_ns() - started)
    ordered = sorted(samples)
    return {
        "repetitions": repetitions,
        "result_count": result,
        "samples_ns": samples,
        "median_ns": statistics.median(samples),
        # Nearest-rank percentile: ceil(0.95 * N), converted to zero-based indexing.
        "p95_ns": ordered[max(0, min(len(ordered) - 1, (95 * len(ordered) + 99) // 100 - 1))],
    }


def backend_record(callable_object: object) -> dict[str, object]:
    callable_type = type(callable_object)
    type_name = f"{callable_type.__module__}.{callable_type.__qualname__}"
    lowered = type_name.lower()
    if "boost.python" in lowered:
        backend = "boost_python"
    elif "nanobind" in lowered or "nb_func" in lowered:
        backend = "nanobind"
    else:
        backend = None
    return {"value": backend, "evidence": f"runtime callable type={type_name}"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-rdkit", required=True)
    parser.add_argument("--artifact", type=Path)
    parser.add_argument("--rows", type=int, default=1000)
    parser.add_argument("--repetitions", type=int, default=7)
    args = parser.parse_args()
    if args.rows < 1:
        parser.error("--rows must be positive")
    if args.repetitions < 3:
        parser.error("--repetitions must be at least 3")
    if not args.corpus.is_file():
        parser.error(f"corpus not found: {args.corpus}")
    if args.artifact is not None and not args.artifact.is_file():
        parser.error(f"artifact not found: {args.artifact}")
    return args


def main() -> int:
    args = parse_args()
    from rdkit import Chem, RDLogger, rdBase
    from rdkit.Chem import rdFingerprintGenerator

    RDLogger.DisableLog("rdApp.*")
    smiles = load_smiles(args.corpus, args.rows)
    parsed: list[object] = []
    failures: list[dict[str, object]] = []
    for index, value in enumerate(smiles):
        try:
            molecule = Chem.MolFromSmiles(value)
        except Exception as exc:
            failures.append(
                {"input_index": index, "smiles": value, "stage": "parse", "error": repr(exc)}
            )
            continue
        if molecule is None:
            failures.append(
                {"input_index": index, "smiles": value, "stage": "parse", "error": "returned_none"}
            )
        else:
            parsed.append(molecule)

    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    with tempfile.TemporaryDirectory(prefix="rdkit-contract-") as temporary:
        sdf_path = Path(temporary) / "probe.sdf"
        sdf_path.write_text(
            Chem.MolToMolBlock(Chem.MolFromSmiles("CCO")) + "\n$$$$\n",
            encoding="utf-8",
        )
        contract = {
            "str": observe(lambda: Chem.MolFromSmiles("CCO") is not None),
            "bytes": observe(lambda: Chem.MolFromSmiles(b"CCO") is not None),
            "path_like": observe(lambda: len(Chem.SDMolSupplier(sdf_path))),
            "keyword": observe(lambda: Chem.MolFromSmiles(SMILES="CCO") is not None),
            "invalid_type": observe(lambda: Chem.MolFromSmiles(123)),
            "invalid_smiles": observe(lambda: Chem.MolFromSmiles("C1(")),
            "morgan_iterable": {
                "list": observe(lambda: generator.GetFingerprints(parsed, numThreads=1)),
                "tuple": observe(lambda: generator.GetFingerprints(tuple(parsed), numThreads=1)),
                "generator": observe(
                    lambda: generator.GetFingerprints((mol for mol in parsed), numThreads=1)
                ),
            },
        }

    def parse_scalar() -> int:
        return sum(Chem.MolFromSmiles(value) is not None for value in smiles)

    def morgan_scalar() -> int:
        return sum(generator.GetFingerprint(mol) is not None for mol in parsed)

    def morgan_batch() -> int:
        return len(generator.GetFingerprints(parsed, numThreads=1))

    def empty_batch_boundary() -> int:
        return len(generator.GetFingerprints([], numThreads=1))

    distribution = importlib.metadata.distribution("rdkit")
    report = {
        "schema_version": 1,
        "profile": "rdkit_python_binding_contract_v1",
        "package": {
            "name": distribution.metadata.get("Name", "rdkit"),
            "version": distribution.version,
            "artifact": file_record(args.artifact) if args.artifact is not None else None,
        },
        "runtime": {
            "rdkit_version": rdBase.rdkitVersion,
            "backend": backend_record(Chem.MolFromSmiles),
            "python_version": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "executable": str(Path(sys.executable).resolve()),
            "platform": platform.platform(),
        },
        "corpus": file_record(args.corpus),
        "row_accounting": {
            "input_count": len(smiles),
            "success_count": len(parsed),
            "failed_count": len(failures),
            "failures": failures,
        },
        "contract": contract,
        "timing": {
            "classification": "python_binding_end_to_end; not chemistry_kernel_timing",
            "scalar_parse": timed(parse_scalar, args.repetitions),
            "scalar_morgan": timed(morgan_scalar, args.repetitions),
            "batch_morgan": timed(morgan_batch, args.repetitions),
            "empty_batch_boundary": timed(
                empty_batch_boundary, max(100, args.repetitions * 20)
            ),
            "parse_batch_api": {"availability": "unavailable"},
        },
    }
    report["expected_rdkit"] = args.expected_rdkit
    report["gate"] = {
        "runtime_matches": rdBase.rdkitVersion == args.expected_rdkit,
        "row_accounting_complete": len(smiles) == len(parsed) + len(failures),
    }
    report["gate_passed"] = all(report["gate"].values())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["gate_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
