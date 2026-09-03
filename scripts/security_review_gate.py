#!/usr/bin/env python3
"""Run the repository-local S4/S5 security review gate.

This gate is deliberately dependency-free and offline-friendly. It checks that
the reviewed security controls and their executable evidence are present; it
does not claim to replace an independent maintainer or external audit.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import tomllib
from pathlib import Path


def run(root: Path, command: list[str]) -> tuple[bool, str]:
    result = subprocess.run(command, cwd=root, capture_output=True, text=True)
    output = (result.stdout + result.stderr).strip()
    return result.returncode == 0, output


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    workspace = tomllib.loads((root / "Cargo.toml").read_text(encoding="utf-8"))
    workspace_version = workspace["workspace"]["package"]["version"]
    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "status": "pass" if passed else "fail", "detail": detail})

    fixture_path = root / "validation" / "cross_binding_contract.json"
    try:
        fixture = json.loads(fixture_path.read_text(encoding="utf-8"))
        adversarial = fixture["adversarial"]
        formats = {case["format"] for case in adversarial}
        check(
            "shared-adversarial-fixture",
            fixture["schema_version"] == 1 and len(adversarial) == 8 and len(formats) == 8,
            f"{len(adversarial)} cases across {len(formats)} formats",
        )
    except (OSError, KeyError, TypeError, ValueError) as error:
        check("shared-adversarial-fixture", False, str(error))

    required = [
        "crates/chematic-mol/tests/cross_binding_adversarial.rs",
        "crates/chematic-smiles/tests/cross_binding_contract.rs",
        "crates/chematic-py/tests/test_cross_binding_contract.py",
        "crates/chematic-wasm/tests/cross_binding_contract.test.mjs",
        ".github/workflows/miri.yml",
        ".github/workflows/sanitizers.yml",
        "scripts/generate_sbom.py",
        "scripts/generate_provenance.py",
        "scripts/check_workflow_pins.py",
    ]
    missing = [path for path in required if not (root / path).is_file()]
    check("reviewed-control-files", not missing, "all present" if not missing else ", ".join(missing))

    formats_target = (root / "fuzz" / "fuzz_targets" / "formats.rs").read_text(encoding="utf-8")
    selector = re.search(r"match selector % (\d+)", formats_target)
    check("parser-dispatch-fuzz", selector is not None and selector.group(1) == "22", "22 parser selectors")

    ok, output = run(root, [sys.executable, "scripts/check_workflow_pins.py"])
    check("immutable-workflow-pins", ok, output.splitlines()[-1] if output else "no output")
    ok, output = run(root, [sys.executable, "scripts/check_unsafe_surface.py"])
    check("unsafe-surface-allowlist", ok, output.splitlines()[-1] if output else "no output")

    failed = [item for item in checks if item["status"] != "pass"]
    report = {
        "schema": "chematic.security-review-gate.v1",
        "version": workspace_version,
        "checks": checks,
        "status": "fail" if failed else "pass",
        "external_review": "required-before-claiming-s5-exit",
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json_out:
        output_path = args.json_out if args.json_out.is_absolute() else root / args.json_out
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
