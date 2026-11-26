#!/usr/bin/env python3
"""Measure checkpoint overhead with deterministic synthetic conformers."""
import argparse
import hashlib
import importlib.util
import json
import platform
import tempfile
import time
from pathlib import Path

import chematic._timed_conformers as timed


class SyntheticMolecule:
    def __init__(self, atoms):
        self.coords = [[i * 0.1, i * 0.2, i * -0.1] for i in range(atoms)]

    def etkdg_embed_multiple_confs_mixed_checked_diagnostics(self, count, **kwargs):
        return [(1, 0.0, self.coords, None)], 1, None, None, None, None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--wrapper", type=Path, help="load a frozen wrapper source for the comparison")
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("repeats must be positive")
    global timed
    if args.wrapper:
        spec = importlib.util.spec_from_file_location("checkpoint_benchmark_wrapper", args.wrapper)
        timed = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(timed)
    wrapper = Path(timed.__file__)
    initial_hash = hashlib.sha256(wrapper.read_bytes()).hexdigest()
    rows = []
    for count in (100, 400, 800):
        for repeat in range(args.repeats):
            with tempfile.TemporaryDirectory(prefix="chematic-checkpoint-profile-") as directory:
                path = Path(directory) / "checkpoint.json"
                request = {"num_confs": count, "random_seed": 42, "max_iterations": 0,
                           "max_mmff_iters": 200, "repair_mmff_iters": 0, "checkpoint": str(path)}
                start = time.perf_counter()
                result = timed._run_partial(SyntheticMolecule(60), request)
                elapsed = time.perf_counter() - start
                start = time.perf_counter()
                recovered = timed._read_checkpoint(str(path))
                read_elapsed = time.perf_counter() - start
                assert recovered == json.loads(json.dumps(result))
                rows.append({"attempts": count, "repeat": repeat, "write_seconds": elapsed,
                             "read_seconds": read_elapsed, "bytes": path.stat().st_size})
    if hashlib.sha256(wrapper.read_bytes()).hexdigest() != initial_hash:
        raise RuntimeError("wrapper source changed during measurement")
    args.output.write_text(json.dumps({"platform": platform.platform(),
        "python": platform.python_version(), "wrapper_sha256": initial_hash,
        "atoms_per_conformer": 60, "synthetic": True, "rows": rows}, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
