#!/usr/bin/env python3
"""Run a CI command with a parent-process wall-clock deadline."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time


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
    parser.add_argument("--heartbeat-seconds", type=int, default=30)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if args.timeout_seconds <= 0 or args.heartbeat_seconds <= 0 or not command:
        parser.error("provide positive timeouts and a command after --")
    print(f"Starting bounded CI command: {command[0]}", flush=True)
    process = subprocess.Popen(command)
    start = time.monotonic()
    while True:
        remaining = args.timeout_seconds - (time.monotonic() - start)
        if remaining <= 0:
            print(
                f"CI command exceeded {args.timeout_seconds}s; terminating process tree: {command[0]}",
                file=sys.stderr,
                flush=True,
            )
            terminate_tree(process)
            return 124
        try:
            return process.wait(timeout=min(args.heartbeat_seconds, remaining))
        except subprocess.TimeoutExpired:
            print(
                f"CI command still running after {int(time.monotonic() - start)}s "
                f"(pid {process.pid})",
                flush=True,
            )


if __name__ == "__main__":
    raise SystemExit(main())
