#!/usr/bin/env python3
"""Run every Windows Python binding test file in bounded, logged batches."""

from __future__ import annotations

import argparse
from collections import deque
import os
from pathlib import Path
import subprocess
import sys
import time


def terminate_tree(process: subprocess.Popen[bytes]) -> None:
    if process.poll() is not None:
        return
    if os.name == "nt":
        try:
            subprocess.run(
                ["taskkill", "/T", "/F", "/PID", str(process.pid)],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
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


def run_with_deadline(
    command: list[str], log_path: Path, timeout_seconds: float, heartbeat_seconds: float
) -> int:
    """Capture child output to a file so the parent deadline cannot block on it."""
    start = time.monotonic()
    with log_path.open("wb") as log:
        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
        while True:
            elapsed = time.monotonic() - start
            remaining = timeout_seconds - elapsed
            if remaining <= 0:
                terminate_tree(process)
                return 124
            try:
                return process.wait(timeout=min(heartbeat_seconds, remaining))
            except subprocess.TimeoutExpired:
                elapsed = time.monotonic() - start
                if elapsed >= timeout_seconds:
                    terminate_tree(process)
                    return 124
                print(f"CI command still running after {int(elapsed)}s", flush=True)


def print_log_tail(log_path: Path, max_lines: int = 50) -> None:
    with log_path.open(encoding="utf-8", errors="replace") as log:
        lines = deque((line.rstrip()[:500] for line in log), maxlen=max_lines)
    for line in lines:
        print(line, file=sys.stderr)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tests-dir", type=Path, required=True)
    parser.add_argument("--report-dir", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=20)
    parser.add_argument("--batch-timeout-seconds", type=int, default=300)
    parser.add_argument("--pytest-timeout-seconds", type=int, default=90)
    parser.add_argument("--session-timeout-seconds", type=int, default=240)
    args = parser.parse_args()
    if min(
        args.batch_size,
        args.batch_timeout_seconds,
        args.pytest_timeout_seconds,
        args.session_timeout_seconds,
    ) <= 0:
        parser.error("batch size and timeouts must be positive")
    if args.session_timeout_seconds >= args.batch_timeout_seconds:
        parser.error("session timeout must be shorter than the parent deadline")
    files = sorted(args.tests_dir.rglob("test_*.py"))
    if not files:
        parser.error(f"no test_*.py files found in {args.tests_dir}")
    args.report_dir.mkdir(parents=True, exist_ok=True)
    batches = [files[i : i + args.batch_size] for i in range(0, len(files), args.batch_size)]
    for index, batch in enumerate(batches, 1):
        print(
            f"pytest batch {index}/{len(batches)}: {len(batch)} files "
            f"({batch[0].name} through {batch[-1].name})",
            flush=True,
        )
        report = args.report_dir / f"batch-{index:02d}.xml"
        log = args.report_dir / f"batch-{index:02d}.log"
        command = [
            sys.executable,
            "-m",
            "pytest",
            *(str(path) for path in batch),
            "-q",
            "-x",
            "--tb=short",
            f"--timeout={args.pytest_timeout_seconds}",
            "--timeout-method=thread",
            f"--session-timeout={args.session_timeout_seconds}",
            f"--junitxml={report}",
        ]
        result = run_with_deadline(command, log, args.batch_timeout_seconds, 30)
        if result:
            print_log_tail(log)
            print(f"pytest batch {index} failed with exit code {result}", file=sys.stderr)
            return result if result > 0 else 1
        print(f"pytest batch {index}/{len(batches)} passed", flush=True)
    print(f"All {len(files)} Python binding test files passed", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
