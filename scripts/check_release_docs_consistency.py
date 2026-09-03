#!/usr/bin/env python3
"""Check version, product-name, and release-key documentation invariants."""

from __future__ import annotations

import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DOCS = (
    ROOT / "README.md",
    ROOT / "README_ja.md",
    ROOT / "CHANGELOG.md",
    ROOT / "SECURITY.md",
    ROOT / "docs" / "v1.0-local-release-gate.md",
    ROOT / "docs" / "compatibility-scope.md",
    ROOT / "docs" / "use-cases" / "rust-server.md",
    ROOT / "docs" / "release-key-custody.md",
    ROOT / "crates" / "chematic-mcp" / "README.md",
    ROOT / "crates" / "chematic-inchi" / "README.md",
)


def main() -> int:
    errors: list[str] = []
    cargo = (ROOT / "Cargo.toml").read_text(encoding="utf-8")
    if not re.search(r'^version\s*=\s*"1\.0\.1"\s*$', cargo, re.MULTILINE):
        errors.append("Cargo.toml does not declare workspace version 1.0.1")

    for path in DOCS:
        try:
            text = path.read_text(encoding="utf-8")
        except OSError as exc:
            errors.append(f"{path.relative_to(ROOT)}: {exc}")
            continue
        relative = path.relative_to(ROOT)
        if "SCHEMATIC_RELEASE_PRIVATE_KEY" in text:
            errors.append(f"{relative}: stale SCHEMATIC release-key secret name")
        if relative != Path("README_ja.md") and "chematic" not in text:
            errors.append(f"{relative}: missing chematic product name")

    custody = (ROOT / "docs" / "release-key-custody.md").read_text(encoding="utf-8")
    if "CHEMATIC_RELEASE_PRIVATE_KEY" not in custody:
        errors.append("release-key custody document omits CHEMATIC_RELEASE_PRIVATE_KEY")
    gate = (ROOT / "docs" / "v1.0-local-release-gate.md").read_text(encoding="utf-8")
    for required in (
        "cargo test --workspace --all-targets --locked",
        "cargo test -p chematic-3d --lib -- --ignored",
        "validation/manifests/v1.0.0-long-run-evidence.json",
    ):
        if required not in gate:
            errors.append(f"release gate document omits {required}")
    for path in (
        ROOT / "docs" / "use-cases" / "rust-server.md",
        ROOT / "crates" / "chematic-mcp" / "README.md",
        ROOT / "crates" / "chematic-inchi" / "README.md",
    ):
        if 'version = "1.0.1"' not in path.read_text(encoding="utf-8"):
            errors.append(f"{path.relative_to(ROOT)}: current dependency example is not v1.0.1")

    if errors:
        print("Release documentation consistency failures:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        return 1
    print("Release documentation consistency OK: version, product name, key name, and gate references")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
