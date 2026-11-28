#!/usr/bin/env python3
"""Compare the two checkpoint implementations on exposed development molecules."""
import argparse
import hashlib
import importlib.util
import json
import platform
import tempfile
from pathlib import Path

import chematic
from chematic import chematic as extension


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--before", type=Path, required=True)
    parser.add_argument("--after", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    modules = []
    for index, path in enumerate((args.before, args.after)):
        spec = importlib.util.spec_from_file_location(f"checkpoint_profile_{index}", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        modules.append(module)
    entries = []
    for row in json.loads(args.input.read_text(encoding="utf-8"))["results"]:
        molecule = chematic.from_smiles(row["smiles"]).add_hydrogens()
        results = []
        for module in modules:
            with tempfile.TemporaryDirectory(prefix="chematic-equivalence-") as directory:
                path = Path(directory) / "checkpoint"
                request = {"num_confs": 3, "random_seed": 42, "max_iterations": 0,
                           "max_mmff_iters": 200, "repair_mmff_iters": 0, "checkpoint": str(path)}
                result = module._run_partial(molecule, request)
                recovered = module._read_checkpoint(str(path))
                assert recovered == json.loads(json.dumps(result))
                results.append(result)
        encoded = [json.dumps(result, sort_keys=True, separators=(",", ":")) for result in results]
        entries.append({"index": row["index"], "name": row["name"],
                        "matches": encoded[0] == encoded[1],
                        "sha256": [hashlib.sha256(s.encode()).hexdigest() for s in encoded],
                        "conformers": [len(r["conformers"]) for r in results]})
        print(f"{len(entries)}", flush=True)
    output = {"input_sha256": digest(args.input), "before_wrapper_sha256": digest(args.before),
              "after_wrapper_sha256": digest(args.after), "extension_sha256": digest(Path(extension.__file__)),
              "python": platform.python_version(), "platform": platform.platform(), "results": entries}
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
