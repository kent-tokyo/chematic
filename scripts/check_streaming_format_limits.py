#!/usr/bin/env python3
"""Run the common streaming benchmark's negative-input format gate.

This is intentionally a small negative-input contract, not a throughput
benchmark. Every format accepted by ``streaming_benchmark`` gets one malformed
or resource-limit rejection and one input-size rejection. The gzip case
additionally proves that the limit is applied after decompression.
"""

from __future__ import annotations

import argparse
import gzip
import json
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "benchmarks" / "fixtures"


def run_runner(binary: list[str], fmt: str, path: Path, *extra: str) -> dict[str, object]:
    command = [*binary, "--format", fmt, "--path", str(path), "--repeats", "1", *extra]
    completed = subprocess.run(command, cwd=ROOT, check=True, text=True, capture_output=True)
    return json.loads(completed.stdout)


def require(condition: bool, message: str, errors: list[str]) -> None:
    if not condition:
        errors.append(message)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--binary",
        nargs="+",
        default=[
            "cargo",
            "run",
            "-p",
            "chematic-mol",
            "--example",
            "streaming_benchmark",
            "--offline",
            "--",
        ],
        help="runner command before --format (default: cargo run --example ...)",
    )
    args = parser.parse_args()

    malformed = {
        "sdf": "broken\n  chematic\n\n  NOTNUM  0  0 V2000\nM  END\n$$$$\n",
        "mol": "broken\n  chematic\n\n  NOTNUM  0  0 V2000\nM  END\n$$$$\n",
        "xyz": "2\nmissing second atom\nC 0 0 0\n",
        "v3000": "not a V3000 mol block\n",
        "mol2": "@<TRIPOS>MOLECULE\nmissing atom and bond sections\n",
        "cml": '<cml><molecule><atomArray><atom id="a1" elementType="C"/><atom id="a2" elementType="C"/></atomArray><bondArray><bond atomRefs2="a1 a2" order="bogus"/></bondArray></molecule></cml>',
        "cdxml": '<CDXML><page><fragment><n id="1" Element="6"/><b id="1"/></fragment></page></CDXML>',
        "mmcif": "data_empty\n",
        "pdb": "ATOM\n",
    }
    # CML/CDXML are deliberately lenient about unknown/empty structure and
    # PDB ignores non-record lines. Exercise their typed safety boundary with
    # an explicit line limit instead of claiming a malformed-record contract
    # they do not expose.
    malformed_options = {
        "cml": ("--max-line-bytes", "1"),
        "cdxml": ("--max-line-bytes", "1"),
        "pdb": ("--max-line-bytes", "1"),
    }
    valid = {
        "sdf": FIXTURES / "streaming.sdf",
        "mol": FIXTURES / "streaming.sdf",
        "xyz": FIXTURES / "streaming.xyz",
        "v3000": FIXTURES / "ethanol.v3000",
        "mol2": FIXTURES / "ethanol.mol2",
        "cml": FIXTURES / "ethanol.cml",
        "cdxml": FIXTURES / "ethanol.cdxml",
        "mmcif": FIXTURES / "minimal.mmcif",
        "pdb": FIXTURES / "minimal.pdb",
    }

    errors: list[str] = []
    with tempfile.TemporaryDirectory(prefix="chematic-streaming-limits-") as directory:
        temp = Path(directory)
        for fmt, content in malformed.items():
            path = temp / f"malformed.{fmt}"
            path.write_text(content, encoding="utf-8")
            result = run_runner(args.binary, fmt, path, *malformed_options.get(fmt, ()))
            require(
                result["records"] == 0 and result["failures"] == 1,
                f"{fmt} malformed case did not fail exactly once: {result}",
                errors,
            )

        for fmt, path in valid.items():
            result = run_runner(args.binary, fmt, path, "--max-input-bytes", "1")
            require(
                result["records"] == 0 and result["failures"] == 1,
                f"{fmt} oversized input was not rejected exactly once: {result}",
                errors,
            )

        gzip_path = temp / "streaming.sdf.gz"
        with gzip.open(gzip_path, "wb") as output:
            output.write(valid["sdf"].read_bytes())
        result = run_runner(args.binary, "sdf", gzip_path, "--gzip")
        require(
            result["records"] == 2 and result["failures"] == 0,
            f"gzip control case did not parse two records: {result}",
            errors,
        )
        result = run_runner(args.binary, "sdf", gzip_path, "--gzip", "--max-input-bytes", "1")
        require(
            result["records"] == 0 and result["failures"] == 1,
            f"gzip decompressed input limit was not enforced: {result}",
            errors,
        )

    if errors:
        print("streaming format limit failures:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        return 1
    print("streaming format limits OK: 9 negative, 9 oversized, 2 gzip cases")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
