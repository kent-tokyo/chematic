#!/usr/bin/env python3
"""Reject mutable GitHub Actions references in checked-in workflows."""

from pathlib import Path
import re
import sys


USES_RE = re.compile(r"^\s*(?:-\s*)?uses:\s*([^\s#]+)")
SHA_RE = re.compile(r"^[0-9a-f]{40}$")


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    workflows = sorted((root / ".github" / "workflows").glob("*.y*ml"))
    violations: list[str] = []
    checked = 0
    for workflow in workflows:
        for line_number, line in enumerate(
            workflow.read_text(encoding="utf-8").splitlines(), start=1
        ):
            match = USES_RE.match(line)
            if not match:
                continue
            checked += 1
            action_ref = match.group(1)
            if "@" not in action_ref:
                violations.append(f"{workflow}:{line_number}: missing @ref")
                continue
            action, ref = action_ref.rsplit("@", 1)
            if not action or not SHA_RE.fullmatch(ref):
                violations.append(
                    f"{workflow}:{line_number}: {action_ref!r} is not an immutable SHA"
                )
    if violations:
        print("Mutable GitHub Actions references found:", file=sys.stderr)
        print("\n".join(violations), file=sys.stderr)
        return 1
    print(f"Workflow action pins verified: {checked} immutable references")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
