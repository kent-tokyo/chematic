#!/usr/bin/env python3
"""Require explicit text encodings in cross-platform Python binding tests."""

from __future__ import annotations

import ast
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PYTHON_TESTS = Path("crates/chematic-py/tests")


def unencoded_read_text_calls(root: Path) -> list[str]:
    """Return Path.read_text() calls that depend on the host locale."""
    failures: list[str] = []
    tests = root / PYTHON_TESTS
    for path in sorted(tests.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            if not isinstance(node.func, ast.Attribute) or node.func.attr != "read_text":
                continue
            has_encoding = bool(node.args) or any(
                keyword.arg == "encoding" for keyword in node.keywords
            )
            if not has_encoding:
                failures.append(f"{path.relative_to(root)}:{node.lineno}")
    return failures


def main() -> int:
    failures = unencoded_read_text_calls(ROOT)
    if failures:
        print("Python binding tests must read text fixtures with an explicit encoding:")
        for failure in failures:
            print(f"  {failure}")
        return 1
    print("Python binding fixture encodings are explicit")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
