#!/usr/bin/env python3
"""Balanced source Morgan batches; controls source-only, every batch checked."""
from __future__ import annotations
import argparse, gzip, hashlib, itertools, json, os, platform, random, resource, shutil, statistics, subprocess, sys, time, zipfile
from pathlib import Path

PRIMARY = ("morgan2_2048", "morgan2_chiral", "morgan3_2048")
CONTROLS = ("tpsa", "labute_asa", "num_rings", "maccs")
MODES = ("first_use", "repeated_use", "parse_inclusive")
ARMS = ("baseline", "candidate", "rdkit", "cosmolkit")


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def encoded(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()


def atomic_json(path, value, compressed=False):
    raw = encoded(value)
    raw = gzip.compress(raw, compresslevel=1, mtime=0) if compressed else raw
    temporary = path.with_name(path.name + ".partial")
    temporary.write_bytes(raw)
    temporary.replace(path)


def rows_of(path, offset, limit):
    if offset < 0 or limit < 1:
        raise ValueError("positive corpus window required")
    rows = [line.split()[0] for line in path.read_text().splitlines() if line.strip()][
        offset : offset + limit
    ]
    if len(rows) != limit:
        raise ValueError("input shorter than requested; no filtering")
    return rows


def check_timing(record, mode, rows, budget):
    times = record["batch_ns"]
    count = len(times)
    if (
        mode not in MODES
        or record["mode"] != mode
        or not count
        or any(type(x) is not int or x <= 0 for x in times)
    ):
        raise ValueError("invalid raw batch timings")
    if record["minimum_accumulated_ns"] != budget or sum(times) < budget:
        raise ValueError("time budget differs or unmet")
    if (
        record["elapsed_ns_raw"] != sum(times)
        or record["batches"] != count
        or record["processed_rows"] != rows * count
    ):
        raise ValueError("native count accounting differs")
    if (
        record["elapsed_ns_per_batch"] != sum(times) // count
        or record["elapsed_ns_per_batch"] <= 0
    ):
        raise ValueError("batch normalization differs")
    hashes = record["output_sha256_per_batch"]
    if len(hashes) != count or len(set(hashes)) != 1:
        raise ValueError("per-batch output check differs")


def measure_mode(
    parse, fn, normalize, rows, mode, budget, expected=None, clock=time.perf_counter_ns
):
    if mode not in MODES or budget < 1 or not rows:
        raise ValueError("invalid workload")

    def parsed():
        result = [parse(s) for s in rows]
        if any(x is None for x in result):
            raise ValueError("parse refusal; no row filtering")
        return result

    prepared = parsed() if mode == "repeated_use" else None
    if prepared is not None:
        for mol in prepared:
            fn(mol)
    durations = []
    hashes = []
    total = 0
    reference = expected
    while total < budget:
        molecules = parsed() if mode == "first_use" else prepared
        started = clock()
        values = (
            [fn(parse(s)) for s in rows]
            if mode == "parse_inclusive"
            else [fn(m) for m in molecules]
        )
        elapsed = clock() - started
        if elapsed <= 0:
            raise ValueError("nonpositive native elapsed")
        normalized = [normalize(value) for value in values]
        if reference is None:
            reference = normalized
        if normalized != reference:
            raise ValueError("output changes across batches or modes")
        hashes.append(digest(encoded(normalized)))
        durations.append(elapsed)
        total += elapsed
        del values, normalized, molecules
    record = {
        "mode": mode,
        "minimum_accumulated_ns": budget,
        "batches": len(durations),
        "processed_rows": len(rows) * len(durations),
        "elapsed_ns_raw": total,
        "elapsed_ns_per_batch": total // len(durations),
        "batch_ns": durations,
        "output_sha256_per_batch": hashes,
    }
    check_timing(record, mode, len(rows), budget)
    return record, reference


def imports(root):
    sys.path.insert(0, str(root / "scripts"))
    sys.path.insert(0, str(root / "validation/cosmolkit_comparison"))
    from bench_corpus import table, paired_median_ratio_ci
    from bench_python_speed_frontier import normal, agreement

    return table, paired_median_ratio_ci, normal, agreement


def worker(args):
    table, _, normal, _ = imports(args.root)
    engine = "chematic" if args.worker in ARMS[:2] else args.worker
    parse, functions = table(engine)
    rows = rows_of(args.corpus, args.offset, args.limit)
    ops = PRIMARY + CONTROLS if engine == "chematic" else PRIMARY
    for op in ops:
        functions[op](parse(rows[0]))
    result = {
        "engine": engine,
        "label": args.worker,
        "offset": args.offset,
        "rows": len(rows),
        "corpus_sha256": digest(args.corpus.read_bytes()),
        "row_window_sha256": digest(encoded(rows)),
        "probe_only": args.probe,
        "operations": {},
    }
    for op in ops:
        budget = 1 if args.probe else (50_000_000 if op in PRIMARY else 20_000_000)
        reference = None
        record = {}
        for mode in MODES:
            detail, reference = measure_mode(
                parse,
                functions[op],
                lambda value: normal(value, engine, op),
                rows,
                mode,
                budget,
                reference,
            )
            record[mode] = detail["elapsed_ns_per_batch"]
            record.setdefault("accounting", {})[mode] = detail
        record["values"] = reference
        result["operations"][op] = record
    module = __import__(engine)
    root = Path(module.__file__).parent
    native = [
        p for p in root.rglob("*") if p.is_file() and p.suffix in (".so", ".dylib")
    ]
    assert len({p.name for p in native}) == len(native), "native basename collision"
    result.update(
        version=module.__version__,
        module_path=module.__file__,
        native_sha256={p.name: digest(p.read_bytes()) for p in native},
        peak_process_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        * (1 if sys.platform == "darwin" else 1024),
    )
    return result


def parent(args):
    _, ratio_fn, _, agreement = imports(args.root)
    out = args.output
    out.mkdir(exist_ok=True)
    assert not (out / "plan.json").exists(), "one formal epoch only"
    wheels = {name: getattr(args, name + "_wheel") for name in ARMS[:2]}
    native = {}
    for name, wheel in wheels.items():
        with zipfile.ZipFile(wheel) as archive:
            native[name] = {
                Path(p).name: digest(archive.read(p))
                for p in archive.namelist()
                if p.endswith((".so", ".dylib"))
            }
    source = json.loads(args.source_manifest.read_text())
    source_sha = digest(args.source_manifest.read_bytes())
    window_sha = digest(encoded(rows_of(args.corpus, args.offset, args.limit)))

    def integrity():
        assert digest(args.source_manifest.read_bytes()) == source_sha
        for name, h in source["candidate_source_sha256"].items():
            assert digest((args.root / name).read_bytes()) == h, name
        for name, wheel in wheels.items():
            assert digest(wheel.read_bytes()) == source["wheels"][name]["sha256"]
        assert digest(Path(__file__).read_bytes()) == source["driver_sha256"]
        for name, h in source["runtime_files_sha256"].items():
            assert digest(Path(name).read_bytes()) == h, name

    order = list(itertools.permutations(ARMS))
    random.Random(20261017 + args.offset).shuffle(order)
    plan = {
        "schema": "python-morgan-linear/v1",
        "source_manifest_sha256": source_sha,
        "source_manifest": source,
        "corpus_sha256": digest(args.corpus.read_bytes()),
        "offset": args.offset,
        "rows": args.limit,
        "orders": [list(x) for x in order],
        "host": {
            "platform": platform.platform(),
            "python": sys.version,
            "cpus": os.cpu_count(),
        },
        "timing": "accumulated per-full-batch operation time; actual counts and every batch duration retained; preparation/validation between batches excluded",
        "state": "24 fixed-order blocks planned; incomplete data cannot be adopted",
    }
    atomic_json(out / "plan.json", plan)
    blocks = []
    expected = {}
    checks = {}
    summary = {}
    expected_metadata = {}
    integrity()
    for index, names in enumerate(order):
        if shutil.disk_usage(out).free < 1024**3:
            raise RuntimeError("stop before block: less than1GiB free")
        record = {"order": list(names), "load_average_start": list(os.getloadavg())}
        for name in names:
            env = os.environ.copy()
            env["PYTHONDONTWRITEBYTECODE"] = "1"
            if name in wheels:
                env["PYTHONPATH"] = str(getattr(args, name + "_path"))
            else:
                env.pop("PYTHONPATH", None)
            command = [
                sys.executable,
                str(Path(__file__).resolve()),
                "--root",
                str(args.root),
                "--worker",
                name,
                "--corpus",
                str(args.corpus),
                "--offset",
                str(args.offset),
                "--limit",
                str(args.limit),
            ]
            process = subprocess.run(command, env=env, capture_output=True, timeout=600)
            if process.returncode:
                (out / f"{index:02}-{name}.stderr").write_bytes(process.stderr)
                atomic_json(
                    out / "failure.json",
                    {
                        "block": index,
                        "arm": name,
                        "returncode": process.returncode,
                        "expected_input_rows": args.limit,
                        "complete": False,
                    },
                )
                raise RuntimeError("worker failure; no completed claim")
            sample = json.loads(process.stdout)
            assert (
                sample["rows"] == args.limit
                and not sample["probe_only"]
                and sample["corpus_sha256"] == plan["corpus_sha256"]
            )
            assert sample["label"] == name and sample["row_window_sha256"] == window_sha
            assert set(sample["operations"]) == set(
                PRIMARY + CONTROLS if name in wheels else PRIMARY
            )
            if name in native:
                assert sample["native_sha256"] == native[name], name
            metadata = {
                k: sample[k] for k in ["version", "module_path", "native_sha256"]
            }
            assert metadata == source["expected_metadata"][name], name
            if name not in expected_metadata:
                expected_metadata[name] = metadata
            assert metadata == expected_metadata[name]
            values = {op: r["values"] for op, r in sample["operations"].items()}
            if name not in expected:
                expected[name] = values
            assert values == expected[name], "non-deterministic outputs"
            for op, r in sample["operations"].items():
                for mode in MODES:
                    check_timing(
                        r["accounting"][mode],
                        mode,
                        args.limit,
                        50_000_000 if op in PRIMARY else 20_000_000,
                    )
            atomic_json(out / f"{index:02}-{name}.json.gz", sample, True)
            for op in sample["operations"]:
                sample["operations"][op].pop("values")
            record[name] = sample
        blocks.append(record)
        atomic_json(out / "progress.json", {"complete": False, "blocks": blocks})
        print("block", index + 1, 24, flush=True)
    integrity()
    for op in PRIMARY + CONTROLS:
        test = {
            "candidate_vs_baseline": agreement(
                expected["baseline"][op], expected["candidate"][op], exact=True
            )
        }
        if op in PRIMARY:
            for name in ("baseline", "candidate", "cosmolkit"):
                test[name + "_vs_rdkit"] = agreement(
                    expected["rdkit"][op], expected[name][op], exact=False
                )
        checks[op] = test
        for mode in MODES:
            entry = {}
            for seed, name in enumerate(("baseline", "rdkit", "cosmolkit")):
                if name != "baseline" and op not in PRIMARY:
                    continue
                allowed = (
                    test["candidate_vs_baseline"]["passed"]
                    if name == "baseline"
                    else test["candidate_vs_rdkit"]["passed"]
                    and (name == "rdkit" or test["cosmolkit_vs_rdkit"]["passed"])
                )
                a = [r[name]["operations"][op][mode] for r in blocks]
                b = [r["candidate"]["operations"][op][mode] for r in blocks]
                assert len(a) == len(b) == 24 and all(x > 0 for x in a + b)
                speed = ratio_fn(a, b, seed=seed) if allowed else None
                entry[name] = {
                    "outputs_equal": allowed,
                    "speedup": speed,
                    "unrounded_paired_median": (
                        statistics.median(x / y for x, y in zip(a, b, strict=True))
                        if allowed
                        else None
                    ),
                    "strict_speed_win": bool(
                        speed and speed["ci95"][0] > 1 and speed["winning_blocks"] == 24
                    ),
                }
            summary[op + "/" + mode] = entry
    report = {
        **plan,
        "complete": True,
        "blocks": blocks,
        "agreement": checks,
        "summary": summary,
        "expected_metadata": expected_metadata,
        "journals": "Per-arm gzip journals retain full canonical values, every batch duration/count/digest; no filtered rows or elapsed reruns",
    }
    atomic_json(out / "result.json", report)
    print("profile complete", args.offset, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--worker", choices=ARMS)
    parser.add_argument("--probe", action="store_true")
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--offset", type=int, default=0)
    parser.add_argument("--limit", type=int, default=1000)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--source-manifest", type=Path)
    for name in ARMS[:2]:
        parser.add_argument("--" + name + "-wheel", type=Path)
        parser.add_argument("--" + name + "-path", type=Path)
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker(args), separators=(",", ":"), allow_nan=False))
    else:
        parent(args)


if __name__ == "__main__":
    main()
