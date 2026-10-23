#!/usr/bin/env python3
"""Fail on invalid or omitted Python type stubs, without importing the package.

Checks Python 3.9 grammar and compiler-level syntax (e.g. duplicate arguments).
Wheel checks inspect archive bytes, require py.typed, and compare the complete
chematic stub inventory with source. This is not a type/API or ABI parity check.
"""
from __future__ import annotations

import argparse
import ast
import glob
import hashlib
import json
import sys
import zipfile
import zlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PACKAGE = ROOT / "crates/chematic-py/python/chematic"
MAX_STUB_BYTES = 4 * 1024 * 1024
MAX_TOTAL_BYTES = 32 * 1024 * 1024


def diagnostic(path, kind, message, line=None, column=None):
    return dict(path=str(path), kind=kind, message=message, line=line, column=column)


def check_syntax(data, filename):
    """Compile to a code object only: never execute imports or annotations."""
    phase = "parse"
    try:
        tree = ast.parse(data, filename=filename, feature_version=(3, 9))
        phase = "compile"
        compile(tree, filename, "exec", dont_inherit=True)
    except (SyntaxError, ValueError, UnicodeError) as error:
        result = diagnostic(filename, "syntax", str(error),
                            getattr(error, "lineno", None), getattr(error, "offset", None))
        result["phase"] = phase
        return result
    return None


def inspect_stubs(blobs, label, errors):
    records = []
    total = 0
    for name, data in sorted(blobs.items()):
        total += len(data)
        if len(data) > MAX_STUB_BYTES or total > MAX_TOTAL_BYTES:
            errors.append(diagnostic(f"{label}:{name}", "limit", "stub byte limit exceeded"))
            continue
        issue = check_syntax(data, f"{label}:{name}")
        if issue:
            errors.append(issue)
        records.append(dict(path=name, bytes=len(data), sha256=hashlib.sha256(data).hexdigest(),
                            status="failed" if issue else "passed"))
    return records


def source_package(package):
    errors, blobs = [], {}
    if not package.is_dir():
        errors.append(diagnostic(package, "source", "source package directory is missing"))
    else:
        total = 0
        for path in sorted(package.rglob("*.pyi")):
            try:
                size = path.stat().st_size
                total += size
                if size > MAX_STUB_BYTES or total > MAX_TOTAL_BYTES:
                    raise ValueError("stub byte limit exceeded")
                blobs[path.relative_to(package).as_posix()] = path.read_bytes()
            except (OSError, ValueError) as error:
                errors.append(diagnostic(path, "source", str(error)))
    if "__init__.pyi" not in blobs:
        errors.append(diagnostic(package / "__init__.pyi", "missing", "required root stub is missing"))
    marker = package / "py.typed"
    try:
        if marker.stat().st_size > MAX_STUB_BYTES:
            raise ValueError("py.typed byte limit exceeded")
        marker_bytes = marker.read_bytes()
    except (OSError, ValueError) as error:
        marker_bytes = None
        errors.append(diagnostic(marker, "missing", f"required py.typed marker unavailable: {error}"))
    records = inspect_stubs(blobs, str(package), errors)
    return dict(path=str(package), stubs=records, errors=errors), blobs, marker_bytes


def check_wheel(path, source_blobs, source_marker):
    errors, blobs = [], {}
    result = dict(path=str(path), stubs=[], errors=errors)
    try:
        # Hashing is streamed; binary extension members are never loaded.
        with path.open("rb") as stream:
            digest = hashlib.sha256()
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
        result["sha256"] = digest.hexdigest()
        with zipfile.ZipFile(path) as archive:
            infos = archive.infolist()
            names, duplicates = set(), set()
            for info in infos:
                if info.filename in names:
                    duplicates.add(info.filename)
                names.add(info.filename)
            if duplicates:
                raise ValueError("duplicate archive members: " + ", ".join(sorted(duplicates)))
            total = 0
            for info in infos:
                if not info.filename.endswith(".pyi"):
                    continue
                total += info.file_size
                if info.file_size > MAX_STUB_BYTES or total > MAX_TOTAL_BYTES:
                    raise ValueError("stub byte limit exceeded")
                blobs[info.filename] = archive.read(info)
            marker_name = "chematic/py.typed"
            if marker_name not in names:
                errors.append(diagnostic(path, "missing", "wheel has no chematic/py.typed"))
            else:
                info = archive.getinfo(marker_name)
                if info.file_size > MAX_STUB_BYTES:
                    raise ValueError("py.typed byte limit exceeded")
                if archive.read(info) != source_marker:
                    errors.append(diagnostic(f"{path}:{marker_name}", "mismatch", "py.typed differs from source"))
    except (OSError, ValueError, zipfile.BadZipFile, RuntimeError, NotImplementedError, EOFError, zlib.error) as error:
        errors.append(diagnostic(path, "wheel", str(error)))
        return result
    result["stubs"] = inspect_stubs(blobs, str(path), errors)
    expected = {"chematic/" + name: data for name, data in source_blobs.items()}
    actual = {name: data for name, data in blobs.items() if name.startswith("chematic/")}
    for name in sorted(set(expected) - set(actual)):
        errors.append(diagnostic(f"{path}:{name}", "missing", "source stub is absent from wheel"))
    for name in sorted(set(actual) - set(expected)):
        errors.append(diagnostic(f"{path}:{name}", "unexpected", "wheel stub is absent from source"))
    for name in sorted(set(actual) & set(expected)):
        if actual[name] != expected[name]:
            errors.append(diagnostic(f"{path}:{name}", "mismatch", "wheel stub bytes differ from source"))
    return result


def audit(package, wheel_patterns):
    source, blobs, marker = source_package(package)
    errors, wheels, paths = list(source["errors"]), [], set()
    for pattern in wheel_patterns:
        matches = glob.glob(str(pattern))
        if not matches:
            errors.append(diagnostic(pattern, "missing", "wheel pattern matched no files"))
        paths.update(Path(p) for p in matches)
    for path in sorted(paths):
        wheel = check_wheel(path, blobs, marker)
        wheels.append(wheel)
        errors.extend(wheel["errors"])
    return dict(schema="chematic.python-stub-gate.v1", status="failed" if errors else "passed",
                grammar="Python 3.9 (ast feature_version; actual host compiler also checked)",
                host_python=sys.version.split()[0], source=source, wheels=wheels, errors=errors,
                scope="syntax, typed marker and exact source/wheel stub bytes; no type/API/ABI parity or execution")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-package", type=Path, default=DEFAULT_PACKAGE)
    parser.add_argument("--wheel", nargs="+", action="append", default=[], help="wheel paths or quoted glob patterns; no matches fails")
    parser.add_argument("--report", type=Path, help="write JSON report (including failures)")
    args = parser.parse_args(argv)
    result = audit(args.source_package, [p for group in args.wheel for p in group])
    payload = json.dumps(result, indent=2) + "\n"
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(payload, encoding="utf-8")
    for error in result["errors"]:
        where = error["path"]
        if error["line"] is not None:
            where += f':{error["line"]}:{error["column"] or 0}'
        print(f'{where}: {error["kind"]}: {error["message"]}', file=sys.stderr)
    print(f'Python stub gate {result["status"]}: {len(result["source"]["stubs"])} source stubs, {len(result["wheels"])} wheels')
    return 1 if result["errors"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
