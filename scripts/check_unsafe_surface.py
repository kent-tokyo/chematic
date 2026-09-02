#!/usr/bin/env python3
"""Reject executable unsafe code outside the reviewed native-InChI FFI island."""

from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
ALLOWED = ROOT / "crates" / "chematic-inchi" / "src" / "native"
UNSAFE = re.compile(r"\bunsafe\s*(?:\{|fn\b|trait\b|impl\b|extern\b)")


def main() -> int:
    violations: list[str] = []
    for path in sorted((ROOT / "crates").rglob("*.rs")):
        allowed = ALLOWED in path.parents
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            code = line.split("//", 1)[0]
            if UNSAFE.search(code) and not allowed:
                violations.append(f"{path.relative_to(ROOT)}:{number}: {line.strip()}")
    if violations:
        print("unsafe code outside the reviewed native-InChI FFI boundary:", file=sys.stderr)
        print("\n".join(violations), file=sys.stderr)
        return 1
    print("unsafe surface OK: only the reviewed native-InChI FFI boundary is allowed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
