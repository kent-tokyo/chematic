#!/usr/bin/env python3
"""Check the frozen Python and WASM public surfaces for removals.

Rust compatibility is handled by cargo-semver-checks.  This checker covers
the two generated-language surfaces that cargo cannot see.  It intentionally
checks removals and newly-required Python arguments, not type widening or
runtime semantics.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import re
import subprocess
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PYI = Path("crates/chematic-py/python/chematic/__init__.pyi")
PYPROJECT = Path("crates/chematic-py/pyproject.toml")
WASM_BASELINE = ROOT / "demo/pkg/chematic_wasm.d.ts"
WASM_BASELINE_SHA256 = "bca23e85ad3ab1604160ba5cf6e84e39877b2121d0e5c10f4c3b2eacb264928a"


def required_arguments(node: ast.FunctionDef | ast.AsyncFunctionDef) -> int:
    positional = len(node.args.posonlyargs) + len(node.args.args)
    required_positional = positional - len(node.args.defaults)
    required_keyword_only = sum(default is None for default in node.args.kw_defaults)
    return required_positional + required_keyword_only


def public_python_api(text: str) -> dict[str, int | None]:
    tree = ast.parse(text)
    out: dict[str, int | None] = {}

    def insert(name: str, required: int | None) -> None:
        # An overload remains callable if any declared signature accepts the
        # old arity, so retain the least restrictive requirement.
        previous = out.get(name)
        if previous is None or required is None:
            out[name] = required
        else:
            out[name] = min(previous, required)

    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and not node.name.startswith("_"):
            insert(node.name, required_arguments(node))
        elif isinstance(node, ast.ClassDef) and not node.name.startswith("_"):
            out[node.name] = None
            for child in node.body:
                if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)) and (
                    not child.name.startswith("_") or child.name == "__init__"
                ):
                    insert(f"{node.name}.{child.name}", required_arguments(child))
    return out


def check_python(baseline: str) -> None:
    old = public_python_api(subprocess.check_output(["git", "show", f"{baseline}:{PYI}"], cwd=ROOT, text=True))
    new = public_python_api((ROOT / PYI).read_text())
    errors = []
    for name, required in old.items():
        if name not in new:
            errors.append(f"Python removed: {name}")
            continue
        if required is not None and new[name] is not None and new[name] > required:
            errors.append(f"Python added required arguments: {name} ({required} -> {new[name]})")
    if errors:
        raise SystemExit("\n".join(errors))


def dependency_contract(document: str) -> tuple[str | None, set[str], dict[str, set[str]]]:
    project = tomllib.loads(document)["project"]
    optional = {
        name: set(values)
        for name, values in project.get("optional-dependencies", {}).items()
    }
    return project.get("requires-python"), set(project.get("dependencies", [])), optional


def check_python_dependencies(baseline: str) -> None:
    old = dependency_contract(subprocess.check_output(
        ["git", "show", f"{baseline}:{PYPROJECT}"], cwd=ROOT, text=True
    ))
    new = dependency_contract((ROOT / PYPROJECT).read_text())
    errors = []
    if old[0] != new[0]:
        errors.append(f"Python requires-python changed: {old[0]!r} -> {new[0]!r}")
    removed = sorted(old[1] - new[1])
    if removed:
        errors.append("Python runtime dependency contract changed or removed: " + ", ".join(removed))
    for group, requirements in old[2].items():
        if group not in new[2]:
            errors.append(f"Python optional dependency group removed: {group}")
            continue
        changed = sorted(requirements - new[2][group])
        if changed:
            errors.append(f"Python optional dependency contract changed in {group}: " + ", ".join(changed))
    if errors:
        raise SystemExit("\n".join(errors))


def dts_surface(text: str) -> set[str]:
    surface = set(re.findall(r"^export (?:function|class|interface|type|enum) ([A-Za-z_$][\w$]*)", text, re.M))
    current = None
    depth = 0
    for line in text.splitlines():
        match = re.match(r"^export class ([A-Za-z_$][\w$]*)", line)
        if match:
            current = match.group(1)
            depth = line.count("{") - line.count("}")
            continue
        if current is not None:
            member = re.match(r"^\s+(?:readonly\s+)?([A-Za-z_$][\w$]*)\s*(?:\(|:)", line)
            if member and member.group(1) != "constructor":
                surface.add(f"{current}.{member.group(1)}")
            depth += line.count("{") - line.count("}")
            if depth <= 0:
                current = None
    return surface


def check_wasm(candidate: Path) -> None:
    digest = hashlib.sha256(WASM_BASELINE.read_bytes()).hexdigest()
    if digest != WASM_BASELINE_SHA256:
        raise SystemExit(f"WASM baseline drifted: {digest}")
    removed = sorted(dts_surface(WASM_BASELINE.read_text()) - dts_surface(candidate.read_text()))
    if removed:
        raise SystemExit("WASM declarations removed:\n" + "\n".join(removed))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", default="v1.1.9")
    parser.add_argument("--wasm-dts", type=Path)
    args = parser.parse_args()
    check_python(args.baseline)
    check_python_dependencies(args.baseline)
    if args.wasm_dts:
        check_wasm(args.wasm_dts)
    print("Python/WASM public API compatibility passed")


if __name__ == "__main__":
    main()
