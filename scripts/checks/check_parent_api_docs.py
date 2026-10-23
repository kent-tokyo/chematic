#!/usr/bin/env python3
"""Check the native Parent surface, typed contracts and rendered API anchors."""
from __future__ import annotations

import argparse
import ast
from html.parser import HTMLParser
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
METHODS = (
    "fragment_parent",
    "charge_parent",
    "isotope_parent",
    "stereo_parent",
    "tautomer_parent",
    "super_parent",
    "super_parent_report",
)
BUDGETED = METHODS[4:]
REPORT_TYPES = ("ParentStage", "SuperParentReport")


class Anchors(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids = set()

    def handle_starttag(self, tag, attrs):
        for name, value in attrs:
            if name == "id" and value is not None:
                self.ids.add(value)


def imported_names(tree):
    return {
        alias.asname or alias.name
        for node in tree.body
        if isinstance(node, ast.ImportFrom) and node.module == "_parent_types"
        for alias in node.names
    }


def audit(root, site_dir=None):
    errors = []
    package = root / "crates/chematic-py/python/chematic"
    native = (
        root / "crates/chematic-py/src/mol_methods/transformations.rs"
    ).read_text()
    stub = ast.parse((package / "__init__.pyi").read_text())
    mol = next(n for n in stub.body if isinstance(n, ast.ClassDef) and n.name == "Mol")
    methods = {
        name: [n for n in mol.body if isinstance(n, ast.FunctionDef) and n.name == name]
        for name in METHODS
    }
    for name in METHODS:
        if not re.search(r"\bfn\s+" + name + r"\b", native):
            errors.append(f"native method missing: {name}")
        if len(methods[name]) != 1:
            errors.append(f"expected one typed Mol.{name}; found {len(methods[name])}")
            continue
        method = methods[name][0]
        expected = (
            "SuperParentReport"
            if name == "super_parent_report"
            else ("tuple[Mol, str]" if name in BUDGETED else "Mol")
        )
        if method.returns is None or ast.unparse(method.returns) != expected:
            errors.append(f"wrong return annotation: Mol.{name}")
        doc = ast.get_docstring(method) or ""
        if "Returns:" not in doc:
            errors.append(f"return documentation missing: Mol.{name}")
        if name in BUDGETED:
            signature = re.search(
                r"#\[pyo3\(signature\s*=\s*\(([^)]*)\)\)\]\s*fn\s+" + name + r"\b",
                native,
            )
            if signature is None:
                errors.append(f"native Python signature missing: Mol.{name}")
            else:
                native_args = (
                    ast.parse("def f(" + signature.group(1) + "):\n    pass")
                    .body[0]
                    .args
                )
                args = method.args
                if [a.arg for a in args.args] != ["self"] + [
                    a.arg for a in native_args.args
                ] or [ast.literal_eval(d) for d in args.defaults] != [
                    ast.literal_eval(d) for d in native_args.defaults
                ]:
                    errors.append(f"wrong native defaults/signature: Mol.{name}")
            if "Args:" not in doc or "Raises:" not in doc:
                errors.append(f"argument/exception documentation missing: Mol.{name}")
    runtime = ast.parse((package / "__init__.py").read_text())
    for label, tree in [("runtime", runtime), ("stub", stub)]:
        for name in REPORT_TYPES:
            if name not in imported_names(tree):
                errors.append(f"{label} report export missing: {name}")
    guide = root / "docs/reference/standardization-parents.md"
    if not guide.is_file():
        errors.append("Parent guide missing")
    else:
        guide_text = guide.read_text()
        for name in METHODS:
            if f"#chematic.Mol.{name}" not in guide_text:
                errors.append(f"guide API link missing: {name}")
    api = (root / "docs/api/chematic.md").read_text()
    for name in ("Mol",) + REPORT_TYPES:
        if not re.search(r"^\s+- " + name + r"\s*$", api, re.M):
            errors.append(f"API reference member missing: {name}")
    checked_anchors = []
    if site_dir is not None:
        html = site_dir / "api/chematic/index.html"
        if not html.is_file():
            errors.append("rendered API HTML missing")
        else:
            parser = Anchors()
            parser.feed(html.read_text())
            checked_anchors = [f"chematic.Mol.{name}" for name in METHODS] + [
                f"chematic.{name}" for name in REPORT_TYPES
            ]
            for anchor in checked_anchors:
                if anchor not in parser.ids:
                    errors.append(f"rendered API anchor missing: {anchor}")
        if not (site_dir / "standardization-parents/index.html").is_file():
            errors.append("rendered Parent guide missing")
    return {
        "status": "failed" if errors else "passed",
        "methods": list(METHODS),
        "report_types": list(REPORT_TYPES),
        "checked_anchors": checked_anchors,
        "errors": errors,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--site-dir", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        report = audit(args.root, args.site_dir)
    except (OSError, SyntaxError, ValueError, StopIteration) as error:
        report = {"status": "failed", "errors": [str(error)]}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    for error in report["errors"]:
        print(error)
    if report["status"] == "passed":
        print(
            f"Parent API docs OK: {len(METHODS)} methods, {len(REPORT_TYPES)} report types; "
            + ("rendered HTML verified" if args.site_dir else "source verified")
        )
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
