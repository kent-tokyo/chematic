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
    parser.add_argument("--timeout-seconds", type=float, default=120.0,
                        help="maximum adapter runtime (default: 120 seconds)")
    args = parser.parse_args()
    if args.timeout_seconds <= 0:
        parser.error("--timeout-seconds must be positive")
    command = shlex.split(args.adapter) + ["--corpus", str(CORPUS), "--engine", args.engine]
    try:
        result = subprocess.run(command, check=True, text=True, capture_output=True,
                                timeout=args.timeout_seconds)
    except FileNotFoundError as exc:
        raise SystemExit(f"adapter executable not found: {exc.filename}") from exc
    except subprocess.TimeoutExpired as exc:
        raise SystemExit(
            f"adapter exceeded timeout of {args.timeout_seconds:g} seconds"
        ) from exc
    except subprocess.CalledProcessError as exc:
        sys.stderr.write(exc.stderr)
        raise SystemExit(exc.returncode) from exc
    # Validate before publishing the output path. A failed adapter run must not
    # leave a misleading partial result that a later command could consume.
    candidate = args.output.with_name(args.output.name + ".candidate")
    candidate.write_text(result.stdout)
    errors = validate(candidate, expected_engine=args.engine)
    if errors:
        candidate.unlink(missing_ok=True)
        raise SystemExit("adapter emitted invalid common-schema output: " + "; ".join(errors))
    candidate.replace(args.output)


if __name__ == "__main__":
    main()
