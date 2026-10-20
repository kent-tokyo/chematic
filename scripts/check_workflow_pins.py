#!/usr/bin/env python3
"""Reject mutable GitHub Actions references in checked-in workflows."""

from pathlib import Path
import re
import sys


USES_RE = re.compile(r"^\s*(?:-\s*)?uses:\s*([^\s#]+)")
SHA_RE = re.compile(r"^[0-9a-f]{40}$")


def reference_error(reference: str, root: Path) -> str | None:
    if reference.startswith("./"):
        root = root.resolve()
        path = (root / reference).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            return f"{reference!r} is not a checked-in local workflow"
        if not path.is_relative_to(root / ".github" / "workflows"):
            return f"{reference!r} is outside .github/workflows"
        return None
    if "@" not in reference:
        return "missing @ref"
    action, ref = reference.rsplit("@", 1)
    if not action or not SHA_RE.fullmatch(ref):
        return f"{reference!r} is not an immutable SHA"
    return None


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
            error = reference_error(action_ref, root)
            if error:
                violations.append(f"{workflow}:{line_number}: {error}")
    if violations:
        print("Mutable GitHub Actions references found:", file=sys.stderr)
        print("\n".join(violations), file=sys.stderr)
        return 1
    print(f"Workflow action pins verified: {checked} immutable or same-commit references")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
