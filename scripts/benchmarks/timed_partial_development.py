#!/usr/bin/env python3
"""Compare batch and opt-in checkpointed timeouts on a development input."""
import argparse
import hashlib
import json
import platform
from pathlib import Path

import chematic
from chematic import chematic as extension


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seconds", type=float, default=0.5)
    args = parser.parse_args()
    records = json.loads(args.input.read_text())["results"]
    results = []
    for index, record in enumerate(records):
        result = {"name": record["name"], "index": record["index"]}
        # Alternate order to reduce systematic order effects.
        order = (False, True) if index % 2 == 0 else (True, False)
        result["order"] = ["partial" if flag else "batch" for flag in order]
        for partial in order:
            response = chematic.generate_conformers_timed(
                record["smiles"], 20, timeout_seconds=args.seconds,
                return_partial=partial,
            )
            result["partial" if partial else "batch"] = {
                "status": response["status"],
                "attempts_completed": response.get("attempts_completed"),
                "conformers": len(response.get("conformers", [])),
                "elapsed_seconds": response["elapsed_seconds"],
            }
        results.append(result)
        print(f"{index + 1}/{len(records)}", flush=True)
    module = Path(chematic.__file__).parent / "_timed_conformers.py"
    output = {"input_sha256": digest(args.input), "wrapper_sha256": digest(module),
              "extension_sha256": digest(Path(extension.__file__)),
              "platform": platform.platform(), "python": platform.python_version(),
              "package_path": str(Path(chematic.__file__).resolve()),
              "seconds": args.seconds, "num_confs": 20, "results": results}
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
