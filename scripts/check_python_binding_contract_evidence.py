#!/usr/bin/env python3
"""Validate the checked-in Python binding contract evidence.

The test run itself belongs to the Python/toolchain lane. This dependency-free
check prevents a passing but stale report from being mistaken for current
workspace evidence.
"""

from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

from benchmark_version import workspace_version


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "validation" / "results" / f"python-binding-contract-v{workspace_version(ROOT)}.json"
SHA256 = re.compile(r"^[0-9a-f]{64}$")


def fail(message: str) -> int:
    print(f"Python binding contract evidence invalid: {message}", file=sys.stderr)
    return 1


def main() -> int:
    try:
        report = json.loads(REPORT.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        return fail(f"cannot read {REPORT}: {error}")

    version = workspace_version(ROOT)
    if report.get("schema_version") != 1:
        return fail("schema_version must be 1")
    if report.get("target_version") != version:
        return fail(f"target_version must be {version}")
    if report.get("status") != "local-verified":
        return fail("status must be local-verified")
    if report.get("gate") != "python_cross_binding_contract":
        return fail("gate name is incorrect")

    collected = report.get("tests_collected")
    passed = report.get("tests_passed")
    failed = report.get("tests_failed")
    skipped = report.get("tests_skipped")
    if not isinstance(collected, int) or collected <= 0:
        return fail("tests_collected must be a positive integer")
    if passed != collected or failed != 0 or skipped != 0:
        return fail("test counts do not describe a complete pass")

    python_version = report.get("python")
    pytest_version = report.get("pytest")
    if not isinstance(python_version, str) or not python_version:
        return fail("Python version is missing")
    if not isinstance(pytest_version, str) or not pytest_version:
        return fail("pytest version is missing")
    command = report.get("command")
    if not isinstance(command, str) or "pytest" not in command or "test" not in command:
        return fail("reproduction command is missing or does not run pytest")

    wheel = report.get("wheel")
    if not isinstance(wheel, dict):
        return fail("wheel metadata is missing")
    digest = wheel.get("sha256")
    if not isinstance(digest, str) or not SHA256.fullmatch(digest):
        return fail("wheel sha256 must be a lowercase 64-character digest")
    wheel_path = wheel.get("path")
    if not isinstance(wheel_path, str) or not wheel_path:
        return fail("wheel path is missing")
    candidate = Path(wheel_path)
    if candidate.is_file():
        actual = hashlib.sha256(candidate.read_bytes()).hexdigest()
        if actual != digest:
            return fail(f"wheel SHA-256 mismatch: expected {digest}, got {actual}")

    print(
        "Python binding contract evidence OK: "
        f"{collected} tests, 0 failures, 0 skips; target={version}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
