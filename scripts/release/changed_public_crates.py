#!/usr/bin/env python3
"""List publishable workspace library crates whose public API may have changed."""

from __future__ import annotations

import argparse
import subprocess
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def public_crates() -> list[tuple[str, Path]]:
    crates = []
    for manifest in sorted((ROOT / "crates").glob("*/Cargo.toml")):
        data = tomllib.loads(manifest.read_text())
        package = data.get("package", {})
        # cargo-semver-checks only accepts crates with a library target.  A
        # publishable binary such as chematic-cli still belongs in the release
        # graph, but has no Rust API surface for this particular check.
        has_library_target = "lib" in data or (manifest.parent / "src/lib.rs").is_file()
        if package and package.get("publish", True) is not False and has_library_target:
            crates.append((package["name"], manifest.parent.relative_to(ROOT)))
    return crates


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", required=True)
    parser.add_argument("--head", default="HEAD")
    args = parser.parse_args()
    changed = set(subprocess.check_output(
        ["git", "diff", "--name-only", f"{args.base}...{args.head}"],
        cwd=ROOT, text=True,
    ).splitlines())
    workspace_wide = bool(changed & {"Cargo.toml", "Cargo.lock"})
    selected = [
        name for name, path in public_crates()
        if workspace_wide or any(item == str(path / "Cargo.toml") or item.startswith(f"{path}/src/") for item in changed)
    ]
    print(" ".join(selected))


if __name__ == "__main__":
    main()
