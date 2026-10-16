#!/usr/bin/env python3
"""Measure RDKit-default TPSA followed by Labute ASA on each same molecule.

Uses the fixed four-arm, fresh-process protocol of the single-operation
frontier. Both native scalar results are retained inside timing; validation,
normalization and bit-exact source comparisons happen afterwards.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import resource
import subprocess
import struct
import sys
import time
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "validation/cosmolkit_comparison"))
from bench_corpus import paired_median_ratio_ci, table  # noqa: E402

OPS = ("tpsa_then_asa",)
MODES = ("first_use", "repeated_use", "parse_inclusive")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def corpus_rows(data: bytes, offset: int, limit: int) -> list[str]:
    if offset < 0 or limit < 1:
        raise ValueError("nonnegative offset and positive limit required")
    rows = [line.split()[0] for line in data.decode().splitlines() if line.strip()][offset:offset + limit]
    if len(rows) != limit:
        raise ValueError("corpus shorter than requested slice")
    return rows


def normal(value, engine: str, op: str):
    if len(value) != 2 or any(not math.isfinite(v) for v in value):
        raise ValueError("bundle requires two finite scalars")
    return list(value)


def worker(args) -> dict:
    engine = "chematic" if args.worker in ("baseline", "candidate") else args.worker
    parse, scalars = table(engine)

    def bundle(mol):
        return scalars["tpsa"](mol), scalars["labute_asa"](mol)

    functions = {"tpsa_then_asa": bundle}
    data = args.corpus.read_bytes()
    rows = corpus_rows(data, args.offset, args.limit)

    def parsed():
        molecules = [parse(s) for s in rows]
        if any(m is None for m in molecules):
            raise ValueError("input parse failed; no timing claim")
        return molecules

    # Warm disjoint objects, preserving first-use input state in the timed lane.
    for op in args.ops.split(","):
        functions[op](parse(rows[0]))
    start = time.perf_counter_ns()
    molecules = parsed()
    elapsed = time.perf_counter_ns() - start
    del molecules
    result = {"parse_ns": elapsed, "corpus_sha256": sha(data), "operations": {}}
    for op in args.ops.split(","):
        fn = functions[op]
        record = {}
        reference_values = None
        for mode in MODES:
            molecules = None if mode == "parse_inclusive" else parsed()
            if mode == "repeated_use":
                for m in molecules:
                    fn(m)
            start = time.perf_counter_ns()
            values = ([fn(parse(s)) for s in rows] if mode == "parse_inclusive"
                      else [fn(m) for m in molecules])
            ns = time.perf_counter_ns() - start
            if ns <= 0:
                raise ValueError("nonpositive timing")
            normalized = [normal(v, engine, op) for v in values]
            if reference_values is not None and reference_values != normalized:
                raise ValueError(f"{op}: output changes across preparation modes")
            reference_values = normalized
            record[mode] = ns
            del values, molecules
        record["values"] = reference_values
        result["operations"][op] = record
    scale = 1 if sys.platform == "darwin" else 1024
    result["peak_process_rss_bytes"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale
    module = __import__(engine)
    result["version"] = module.__version__
    result["module_path"] = module.__file__
    result["extension_sha256"] = {
        p.name: sha(p.read_bytes()) for p in Path(module.__file__).parent.glob("*.so")
    }
    return result


def agreement(reference: list, values: list, *, exact: bool) -> dict:
    if len(values) != len(reference):
        return {"passed": False, "rows": len(reference), "mismatches": len(reference)}
    def equal(a, b):
        if not isinstance(a, (list, tuple)) or not isinstance(b, (list, tuple)):
            return False
        if len(a) != 2 or len(b) != 2:
            return False
        for x, y in zip(a, b, strict=True):
            if (isinstance(x, bool) or isinstance(y, bool)
                    or not isinstance(x, (int, float)) or not isinstance(y, (int, float))
                    or not math.isfinite(x) or not math.isfinite(y)):
                return False
            if exact:
                if struct.pack(">d", x) != struct.pack(">d", y):
                    return False
            elif not math.isclose(x, y, rel_tol=0, abs_tol=1e-9):
                return False
        return True

    mismatches = sum(not equal(a, b) for a, b in zip(reference, values, strict=True))
    return {"passed": mismatches == 0, "rows": len(reference), "mismatches": mismatches}


def comparable(checks: dict, competitor: str) -> bool:
    if competitor == "baseline":
        return checks["candidate_vs_baseline"]["passed"]
    if competitor == "rdkit":
        return checks["candidate_vs_rdkit"]["passed"]
    return checks["candidate_vs_rdkit"]["passed"] and checks[competitor + "_vs_rdkit"]["passed"]


def parent(args) -> None:
    if args.blocks < 20 or args.limit < 1:
        raise ValueError("at least 20 blocks and a positive limit required")
    names = ("baseline", "candidate", "rdkit", "cosmolkit")
    wheels = {name: getattr(args, name + "_wheel") for name in names[:2]}
    installed_extensions = {}
    for name, path in wheels.items():
        with zipfile.ZipFile(path) as archive:
            installed_extensions[name] = {Path(p).name: sha(archive.read(p))
                                          for p in archive.namelist() if p.endswith(".so")}
    report = {
        "schema": "python-scalar-bundle/v1", "base_revision": args.base_revision,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "load_average_start": list(os.getloadavg()),
        "host": {"platform": platform.platform(), "python": sys.version, "cpus": os.cpu_count()},
        "wheels": {name: {"path": str(path), "sha256": sha(path.read_bytes())}
                   for name, path in wheels.items()},
        "corpus": {"path": str(args.corpus), "sha256": sha(args.corpus.read_bytes()),
                   "offset": args.offset, "limit": args.limit},
        "protocol": {"blocks": args.blocks, "order": "four rotating fresh-process arms, sequential",
                     "timing": "native return objects retained; normalization outside timing",
                     "first_use": "fresh parsed molecules; eager/lazy perception differs",
                     "repeated_use": "one untimed operation on every fresh molecule",
                     "parse_inclusive": "actual parse and operation inside the same timed loop",
                     "memory": "whole-process peak, including imports and all operations"},
        "blocks": [], "agreement": {}, "summary": {}, "complete": False,
    }
    report["source_manifest"] = json.loads(args.source_manifest.read_text()) if args.source_manifest else None
    if args.resume:
        previous = json.loads(args.output.read_text())
        for key in ("schema", "base_revision", "source_manifest", "wheels", "corpus", "protocol", "host"):
            if previous[key] != report[key]:
                raise ValueError(f"resume {key} differs")
        if len(previous["blocks"]) > args.blocks:
            raise ValueError("too many existing blocks")
        report = previous
    expected = {}
    args.output.parent.mkdir(parents=True, exist_ok=True)

    def save():
        tmp = args.output.with_suffix(".partial")
        tmp.write_text(json.dumps(report, indent=2) + "\n")
        tmp.replace(args.output)

    # Re-establish the full output values in fresh workers before reusing any
    # saved samples. Validation-only workers do not add timing observations.
    blocks = ([-1] if args.resume else []) + list(range(len(report["blocks"]), args.blocks))
    for block in blocks:
        order = names[block % 4:] + names[:block % 4]
        record = {"order": order, "load_average_start": list(os.getloadavg())}
        for name in order:
            env = os.environ.copy()
            if name in wheels:
                env["PYTHONPATH"] = str(getattr(args, name + "_path"))
            else:
                env.pop("PYTHONPATH", None)
            cmd = [sys.executable, __file__, "--worker", name, "--corpus", str(args.corpus),
                   "--limit", str(args.limit), "--offset", str(args.offset), "--ops", args.ops]
            process = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=300)
            if process.returncode:
                save()
                raise RuntimeError(f"{name}: worker failed, no completed speed claim\n{process.stderr[-2000:]}")
            sample = json.loads(process.stdout)
            if name in wheels and sample["extension_sha256"] != installed_extensions[name]:
                raise ValueError(f"{name}: installed extension does not match the supplied wheel")
            if sample["corpus_sha256"] != report["corpus"]["sha256"]:
                raise ValueError("corpus changed during measurement")
            outputs = {op: value.pop("values") for op, value in sample["operations"].items()}
            if name not in expected:
                expected[name] = outputs
            if outputs != expected[name]:
                raise ValueError(f"{name}: nondeterministic outputs across blocks")
            sample["output_sha256"] = sha(json.dumps(outputs, sort_keys=True).encode())
            if block == -1:
                for previous in report["blocks"]:
                    for key in ("output_sha256", "version", "extension_sha256"):
                        if previous[name][key] != sample[key]:
                            raise ValueError(f"{name}: saved {key} differs from fresh validation")
            record[name] = sample
        if block == -1:
            print("saved sample outputs and extensions revalidated", flush=True)
            continue
        report["blocks"].append(record)
        save()
        print(f"block {block + 1}/{args.blocks}", flush=True)

    for op in args.ops.split(","):
        checks = {"candidate_vs_baseline": agreement(expected["baseline"][op], expected["candidate"][op], exact=True)}
        checks.update({name + "_vs_rdkit": agreement(expected["rdkit"][op], expected[name][op], exact=False)
                       for name in ("baseline", "candidate", "cosmolkit")})
        report["agreement"][op] = checks
        for mode in MODES:
            comparisons = {}
            for index, name in enumerate(("baseline", "rdkit", "cosmolkit")):
                allowed = comparable(checks, name)
                interval = paired_median_ratio_ci(
                    [b[name]["operations"][op][mode] for b in report["blocks"]],
                    [b["candidate"]["operations"][op][mode] for b in report["blocks"]], seed=index,
                ) if allowed else None
                comparisons[name] = {"speedup": interval, "demonstrated_win": bool(
                    interval and interval["ci95"][0] > 1 and interval["winning_blocks"] == args.blocks)}
            report["summary"][op + "/" + mode] = comparisons
    report["complete"] = True
    report["load_average_end"] = list(os.getloadavg())
    save()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", choices=("baseline", "candidate", "rdkit", "cosmolkit"))
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=1000)
    parser.add_argument("--offset", type=int, default=0, help="zero-based nonempty SMILES row offset")
    parser.add_argument("--blocks", type=int, default=20)
    parser.add_argument("--ops", default=",".join(OPS))
    parser.add_argument("--base-revision")
    parser.add_argument("--source-manifest", type=Path)
    parser.add_argument("--resume", action="store_true", help="verify and reuse saved raw blocks")
    parser.add_argument("--output", type=Path)
    for name in ("baseline", "candidate"):
        parser.add_argument("--" + name + "-wheel", type=Path)
        parser.add_argument("--" + name + "-path", type=Path)
    args = parser.parse_args()
    if set(args.ops.split(",")) - set(OPS):
        parser.error("unsupported operation")
    if args.worker:
        print(json.dumps(worker(args)))
    else:
        if any(getattr(args, name) is None for name in (
            "baseline_wheel", "candidate_wheel", "baseline_path", "candidate_path", "output", "base_revision",
        )):
            parser.error("both source wheels, installed paths, base revision and output are required")
        parent(args)


if __name__ == "__main__":
    main()
