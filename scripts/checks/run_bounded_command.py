#!/usr/bin/env python3
"""Run a CI command with a parent-process wall-clock deadline."""

from __future__ import annotations

import argparse
import subprocess
import sys


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout-seconds", type=int, required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if args.timeout_seconds <= 0 or not command:
        parser.error("provide a positive timeout and a command after --")
    try:
        return subprocess.run(command, timeout=args.timeout_seconds, check=False).returncode
    except subprocess.TimeoutExpired:
        print(
            f"CI command exceeded {args.timeout_seconds}s and was terminated: {command[0]}",
            file=sys.stderr,
            flush=True,
        )
        return 124


if __name__ == "__main__":
    raise SystemExit(main())
