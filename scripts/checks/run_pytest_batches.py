#!/usr/bin/env python3
"""Run every Python binding test file in logged, time-bounded batches."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tests-dir", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=20)
    parser.add_argument("--batch-timeout-seconds", type=int, default=240)
    parser.add_argument("--pytest-timeout-seconds", type=int, default=90)
    parser.add_argument("--session-timeout-seconds", type=int, default=180)
    parser.add_argument("--report-dir", type=Path, required=True)
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
    batches = [files[i:i + args.batch_size] for i in range(0, len(files), args.batch_size)]
    for index, batch in enumerate(batches, 1):
        print(
            f"pytest batch {index}/{len(batches)}: {len(batch)} files, "
            f"{batch[0]} through {batch[-1]}",
            flush=True,
        )
        command = [
            sys.executable,
            "scripts/checks/run_bounded_command.py",
            "--timeout-seconds", str(args.batch_timeout_seconds),
            "--heartbeat-seconds", "30",
            "--",
            sys.executable, "-m", "pytest",
            *(str(path) for path in batch),
            "-v", "-x", "--capture=no",
            f"--timeout={args.pytest_timeout_seconds}",
            "--timeout-method=thread",
            f"--session-timeout={args.session_timeout_seconds}",
            f"--junitxml={args.report_dir / f'batch-{index:02d}.xml'}",
        ]
        result = subprocess.run(command, check=False)
        if result.returncode:
            print(f"pytest batch {index} failed with exit code {result.returncode}", file=sys.stderr)
            return result.returncode
    print(f"All {len(files)} Python binding test files passed", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
