#!/usr/bin/env python3
"""Verify the crates.io publish graph is consistent before a release.

For every workspace crate that is publishable (no `publish = false`), every
crate it references via a normal or build dependency (NOT dev-dependency --
Cargo excludes those from the published manifest and crates.io's registry
constraints) must itself be:

  (a) publishable (not `publish = false`)
  (b) pinned to the current workspace version
  (c) reachable in a well-defined topological publish order (no cycles)

This exists because `chematic-cip` was `publish = false` while `chematic-chem`
gained a normal (non-dev) dependency on it in Milestone 5A -- nobody checked
publish-graph validity at merge time, so it wasn't caught until the 0.4.30
release was blocked at crates.io. See CHANGELOG.md's "Release" entry.

Usage:
    python scripts/check_publish_graph.py
"""

import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).parent.parent
CRATES_DIR = ROOT / "crates"


def workspace_version() -> str:
    text = (ROOT / "Cargo.toml").read_text()
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("version"):
            return line.split("=", 1)[1].strip().strip('"')
    sys.exit("Could not find workspace version in Cargo.toml")


def load_crate(path: Path) -> dict:
    with open(path / "Cargo.toml", "rb") as f:
        return tomllib.load(f)


def is_publishable(manifest: dict) -> bool:
    publish = manifest.get("package", {}).get("publish", True)
    # Cargo accepts `publish = false` or `publish = []` as "never publish";
    # any other value (true, or a list of registries) means publishable.
    if publish is False:
        return False
    if isinstance(publish, list) and len(publish) == 0:
        return False
    return True


def chematic_deps(manifest: dict, *sections: str) -> dict:
    """{name: declared version string} for chematic-* path deps in the given sections."""
    deps = {}
    for section in sections:
        for name, spec in manifest.get(section, {}).items():
            if name.startswith("chematic") and isinstance(spec, dict) and "path" in spec:
                deps[name] = spec.get("version")
    return deps


def main() -> int:
    ver = workspace_version()
    crate_dirs = sorted(p for p in CRATES_DIR.iterdir() if (p / "Cargo.toml").exists())

    manifests = {}
    for d in crate_dirs:
        m = load_crate(d)
        name = m["package"]["name"]
        manifests[name] = m

    errors = []
    # Normal-dependency-only graph, for cycle detection -- dev-dependencies are
    # excluded because Cargo strips them from the published manifest and they
    # never need to resolve against the registry.
    normal_graph = {
        name: chematic_deps(m, "dependencies", "build-dependencies")
        for name, m in manifests.items()
    }

    for name, manifest in manifests.items():
        if not is_publishable(manifest):
            continue
        for dep_name, declared_ver in normal_graph[name].items():
            dep_manifest = manifests.get(dep_name)
            if dep_manifest is None:
                errors.append(f"{name}: depends on unknown crate {dep_name}")
                continue
            if not is_publishable(dep_manifest):
                errors.append(
                    f"{name}: normal dependency on {dep_name}, which has "
                    f"`publish = false` -- crates.io cannot resolve this"
                )
            if declared_ver != ver:
                errors.append(
                    f"{name}: declares {dep_name} = version \"{declared_ver}\", "
                    f"workspace is {ver}"
                )

    # Cycle detection (normal-dependency graph only) via DFS.
    WHITE, GRAY, BLACK = 0, 1, 2
    color = {name: WHITE for name in manifests}
    order = []

    def visit(name, stack):
        if color[name] == GRAY:
            cycle = " -> ".join(stack[stack.index(name):] + [name])
            errors.append(f"cycle in normal-dependency graph: {cycle}")
            return
        if color[name] == BLACK:
            return
        color[name] = GRAY
        stack.append(name)
        for dep in sorted(normal_graph.get(name, ())):
            if dep in manifests:
                visit(dep, stack)
        stack.pop()
        color[name] = BLACK
        order.append(name)

    for name in sorted(manifests):
        visit(name, [])

    if errors:
        print("Publish graph check FAILED:")
        for e in errors:
            print(f"  - {e}")
        return 1

    publishable_order = [n for n in order if is_publishable(manifests[n])]
    print(f"Publish graph OK ({len(publishable_order)} publishable crates).")
    print("Safe publish order:")
    for n in publishable_order:
        print(f"  {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
