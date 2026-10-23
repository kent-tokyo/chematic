#!/usr/bin/env python3
"""Verify the installed public Python package's minimal runtime contract.

This script intentionally imports the package already installed into the active
environment. It neither builds nor imports the workspace binding, so a CI job
can use it after ``pip install chematic==<release>`` to attest to a published
wheel on its runner platform.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", required=True, help="expected published chematic version")
    parser.add_argument("--output", type=Path, required=True, help="JSON result path")
    args = parser.parse_args()

    import chematic
    from chematic import from_smiles

    actual_version = getattr(chematic, "__version__", None)
    if actual_version != args.version:
        raise SystemExit(
            f"installed chematic version mismatch: expected {args.version}, got {actual_version!r}"
        )
    benzene = from_smiles("c1ccccc1")
    formula = getattr(benzene, "formula", None)
    if formula != "C6H6":
        raise SystemExit(f"benzene formula mismatch: expected C6H6, got {formula!r}")

    result = {
        "schema_version": 1,
        "gate": "published_python_package_runtime_smoke",
        "expected_version": args.version,
        "installed_version": actual_version,
        "checks": {
            "import_schematic": True,
            "benzene_formula": formula,
        },
        "runtime": {
            "executable": sys.executable,
            "python_version": platform.python_version(),
            "implementation": platform.python_implementation(),
            "platform": platform.platform(),
        },
        "boundary": (
            "Verifies the already-installed public wheel's import, version, and one basic "
            "SMILES runtime operation only. It does not prove typing, all APIs, performance, "
            "or source/workspace behavior."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
