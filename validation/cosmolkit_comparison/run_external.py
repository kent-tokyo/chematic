#!/usr/bin/env python3
"""Run an installed competitor adapter against the shared smoke corpus.

The adapter is an executable (or command prefix) that receives ``--corpus`` and
must emit one common-schema JSON object per corpus row on stdout. Keeping the
adapter outside this repository makes the harness usable with installations
that expose different Python, CLI, or container entry points.
"""

import argparse
import shlex
import subprocess
import sys
from pathlib import Path

from validate import validate

CORPUS = Path(__file__).with_name("smoke_corpus.jsonl")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--engine", required=True, help="name recorded by the adapter")
    parser.add_argument("--adapter", required=True, help="adapter command prefix")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    command = shlex.split(args.adapter) + ["--corpus", str(CORPUS), "--engine", args.engine]
    try:
        result = subprocess.run(command, check=True, text=True, capture_output=True)
    except FileNotFoundError as exc:
        raise SystemExit(f"adapter executable not found: {exc.filename}") from exc
    except subprocess.CalledProcessError as exc:
        sys.stderr.write(exc.stderr)
        raise SystemExit(exc.returncode) from exc
    args.output.write_text(result.stdout)
    errors = validate(args.output)
    if errors:
        args.output.unlink(missing_ok=True)
        raise SystemExit("adapter emitted invalid common-schema output: " + "; ".join(errors))


if __name__ == "__main__":
    main()
