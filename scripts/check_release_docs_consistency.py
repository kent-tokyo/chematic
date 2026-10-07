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
    ROOT / "README_zh.md",
    ROOT / "CHANGELOG.md",
    ROOT / "SECURITY.md",
    ROOT / "docs" / "v1.0-local-release-gate.md",
    ROOT / "docs" / "compatibility-scope.md",
    ROOT / "docs" / "use-cases" / "rust-server.md",
    ROOT / "docs" / "release-key-custody.md",
    ROOT / "crates" / "chematic-mcp" / "README.md",
    ROOT / "crates" / "chematic-inchi" / "README.md",
)


# A relative Markdown link (not an image, URL or anchor) in the MkDocs tree.
LINK = re.compile(r"(?<!!)\[[^\]]*\]\(([^)\s]+)\)")


def site_link_errors() -> list[str]:
    """Links in docs/ (as MkDocs builds it, archive/ excluded) that leave
    docs/: `mkdocs build --strict` fails on them (the v1.0.36 Pages run)."""
    docs = ROOT / "docs"
    errors = []
    for path in sorted(docs.rglob("*.md")):
        if path.relative_to(docs).parts[0] == "archive":
            continue
        for target in LINK.findall(path.read_text(encoding="utf-8")):
            if "://" in target or target.startswith(("#", "mailto:")):
                continue
            resolved = (path.parent / target.split("#", 1)[0]).resolve()
            if not resolved.is_relative_to(docs.resolve()):
                errors.append(
                    f"{path.relative_to(ROOT)}: link {target} leaves docs/ "
                    "(MkDocs cannot resolve it; use the GitHub URL)"
                )
    return errors


def main() -> int:
    errors: list[str] = site_link_errors()
    cargo = (ROOT / "Cargo.toml").read_text(encoding="utf-8")
    version_match = re.search(r'^version\s*=\s*"([^"]+)"\s*$', cargo, re.MULTILINE)
    if version_match is None:
        print("Release documentation consistency failure: workspace version is missing", file=sys.stderr)
        return 1
    release_version = version_match.group(1)

    release_headings = {
        Path("README.md"): f"### v{release_version} release boundary",
        Path("README_ja.md"): f"### v{release_version} の対応範囲",
        Path("README_zh.md"): f"### v{release_version} 范围",
    }

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
        expected_heading = release_headings.get(relative)
        if expected_heading is not None and expected_heading not in text:
            errors.append(f"{relative}: missing current release boundary heading {expected_heading}")

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
        expected = f'version = "{release_version}"'
        if expected not in path.read_text(encoding="utf-8"):
            errors.append(
                f"{path.relative_to(ROOT)}: current dependency example is not v{release_version}"
            )

    if errors:
        print("Release documentation consistency failures:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(
        "Release documentation consistency OK: version, product name, key name, gate "
        "references and site-internal links"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
