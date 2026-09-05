#!/usr/bin/env python3
"""Alternate frozen source binaries; save samples and exact-output checks.

Each executable is hotpath_throughput built from the corresponding source.
Timing processes run sequentially. Interrupted reports retain finished pairs.
"""
import argparse
import hashlib
import json
import math
import platform
import statistics
import subprocess
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("baseline", "candidate", "smiles", "sdf", "output"):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument("--pairs", type=int, default=7)
    p.add_argument("--repeats", type=int, default=5)
    args = p.parse_args()
    if min(args.pairs, args.repeats) < 1:
        p.error("pairs and repeats must be positive")
    paths = {name: getattr(args, name).resolve() for name in ("baseline", "candidate", "smiles", "sdf")}
    provenance = {name: {"path": str(path), "sha256": sha(path)} for name, path in paths.items()}
    config = {"pairs": args.pairs, "repeats": args.repeats, "platform": platform.platform()}
    report = {"provenance": provenance, "configuration": config, "pairs": [], "complete": False}
    if args.output.exists():
        report = json.loads(args.output.read_text())
        if report["provenance"] != provenance or report["configuration"] != config:
            raise SystemExit("refusing to resume with changed binaries, inputs or configuration")
    def run(name, dump=False):
        cmd = [str(paths[name]), str(paths["smiles"]), str(paths["sdf"]), str(args.repeats)]
        if dump:
            cmd.append("--dump")
        return subprocess.run(cmd, check=True, capture_output=True, timeout=180).stdout
    a, b = run("baseline", True), run("candidate", True)
    if a != b:
        raise SystemExit("baseline/candidate output bytes differ")
    report["output_sha256"] = hashlib.sha256(a).hexdigest()
    def save():
        temporary = args.output.with_name(args.output.name + ".tmp")
        temporary.write_text(json.dumps(report, indent=2) + "\n")
        temporary.replace(args.output)
    save()
    for index in range(len(report["pairs"]), args.pairs):
        order = ["baseline", "candidate"] if index % 2 == 0 else ["candidate", "baseline"]
        pair = {name: json.loads(run(name)) for name in order}
        for result in pair.values():
            for lane in ("parse_us", "canonical_us", "sdf_read_us", "sdf_write_us"):
                if not math.isfinite(result[lane]) or result[lane] <= 0:
                    raise SystemExit(f"pair {index}: invalid timing {lane}")
        for field in ("output_fnv1a", "smiles", "records", "repeats"):
            if pair["baseline"][field] != pair["candidate"][field]:
                raise SystemExit(f"pair {index}: {field} mismatch")
        report["pairs"].append(pair)
        save()
        print(f"pair {index + 1}/{args.pairs} saved", flush=True)
    summary = {}
    for lane in ("parse_us", "canonical_us", "sdf_read_us", "sdf_write_us"):
        base = [p["baseline"][lane] for p in report["pairs"]]
        candidate = [p["candidate"][lane] for p in report["pairs"]]
        ratios = [a / b for a, b in zip(base, candidate)]
        summary[lane] = {"baseline_median": statistics.median(base),
                         "candidate_median": statistics.median(candidate),
                         "paired_speedup_median": statistics.median(ratios),
                         "paired_speedup_min": min(ratios), "paired_speedup_max": max(ratios)}
    report["summary"] = summary
    report["complete"] = True
    save()
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
