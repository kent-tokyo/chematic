#!/usr/bin/env python3
"""Run a CI command with a parent-process wall-clock deadline."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys


def terminate_tree(process: subprocess.Popen[bytes]) -> None:
    if os.name == "nt":
        try:
            subprocess.run(
                ["taskkill", "/T", "/F", "/PID", str(process.pid)],
                timeout=30,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired):
            process.kill()
    else:
        process.kill()
    try:
        process.wait(timeout=10)
    except subprocess.TimeoutExpired:
        process.kill()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout-seconds", type=int, required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if args.timeout_seconds <= 0 or not command:
        parser.error("provide a positive timeout and a command after --")
    process = subprocess.Popen(command)
    try:
        return process.wait(timeout=args.timeout_seconds)
    except subprocess.TimeoutExpired:
        print(
            f"CI command exceeded {args.timeout_seconds}s; terminating process tree: {command[0]}",
            file=sys.stderr,
            flush=True,
        )
        terminate_tree(process)
        return 124


if __name__ == "__main__":
    raise SystemExit(main())
