#!/usr/bin/env python3
"""Run the common streaming benchmark's negative-input format gate.

This is intentionally a bounded negative-input contract, not a throughput
benchmark. Every format accepted by ``streaming_benchmark`` gets twelve cases
from the checked-in corpus and one input-size rejection. The gzip cases
additionally prove that the limit is applied after decompression for every
runner format.
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
SAFETY_CORPUS = ROOT / "validation" / "streaming_format_safety_cases.json"


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
        "sdf": [
            "broken\n  chematic\n\n  NOTNUM  0  0 V2000\nM  END\n$$$$\n",
            "broken\n  chematic\n\n  1  0  0  0  0  0            999 V2000\nthis is not an atom line\nM  END\n$$$$\n",
            "broken\n  chematic\n\n  1  1  0  0  0  0            999 V2000\n  0.0  0.0  0.0  C\nM  END\n$$$$\n",
        ],
        "mol": [
            "broken\n  chematic\n\n  NOTNUM  0  0 V2000\nM  END\n$$$$\n",
            "broken\n  chematic\n\n  1  0  0  0  0  0            999 V2000\nthis is not an atom line\nM  END\n$$$$\n",
            "broken\n  chematic\n\n  1  1  0  0  0  0            999 V2000\n  0.0  0.0  0.0  C\nM  END\n$$$$\n",
        ],
        "xyz": [
            "2\nmissing second atom\nC 0 0 0\n",
            "not-a-count\ncomment\nC 0 0 0\n",
            "-1\nnegative atom count\n",
        ],
        "extxyz": [
            "2\nProperties=species:S:1:pos:R:3\nC 0 0 0\n",
            "not-a-count\ncomment\nC 0 0 0\n",
            "1\nProperties=species:S:1:pos:R:3\nXx 0 0 0\n",
            "1\nProperties=species:S:1:pos:R:3\nC nan 0 0\n",
            "1\nProperties=species:S:1:pos:R:3:charge:R:1\nC 0 0 0 nope\n",
            "1\nProperties=species:S:1:pos:R:3\nC 0 0\n",
            "1\nProperties=species:S:1:pos:R:3\nC 0e 0 0\n",
            "1\nLattice=\"1 2 3\"\nC 0 0 0\n",
            "1\nProperties=species:S:1:pos:R\nC 0 0 0\n",
            "1\nProperties=species:S:1:pos:R:3\nC 0 0\nC 1 1 1\n",
        ],
        "v3000": [
            "not a V3000 mol block\n",
            "M  V30 COUNTS not-numbers\nM  END\n",
            "M  V30 BEGIN CTAB\nM  V30 COUNTS 1 1 0 0 0\nM  V30 END CTAB\n",
        ],
        "mol2": [
            "@<TRIPOS>MOLECULE\nmissing atom and bond sections\n",
            "@<TRIPOS>MOLECULE\nname\n1 1 0 0 0\n@<TRIPOS>ATOM\nnot-an-atom\n",
            "@<TRIPOS>MOLECULE\nname\n-1 0 0 0 0\n",
        ],
        "cml": [
            '<cml><molecule><atomArray><atom id="a1" elementType="C"/><atom id="a2" elementType="C"/></atomArray><bondArray><bond atomRefs2="a1 a2" order="bogus"/></bondArray></molecule></cml>',
            "<cml><molecule><atomArray>",
            "<cml><molecule><atomArray><atom id=\"a1\" elementType=\"Xx\"/>",
        ],
        "cdxml": [
            '<CDXML><page><fragment><n id="1" Element="6"/><b id="1"/></fragment></page></CDXML>',
            "<CDXML><page><fragment>",
            '<CDXML><fragment><n id="1" Element="999" p="0 0"/></fragment></CDXML>',
        ],
        "mmcif": [
            "data_empty\n",
            "data_empty\nloop_\n_atom_site.id\n",
            "data_empty\nloop_\n_atom_site.id\n1\n_atom_site.type_symbol\n",
        ],
        "pdb": [
            "ATOM\n",
            "ATOM      1  BAD\n",
            "HETATM not-a-pdb-record\n",
        ],
    }
    try:
        corpus = json.loads(SAFETY_CORPUS.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"streaming safety corpus read failure: {exc}", file=sys.stderr)
        return 1
    expected_formats = {"sdf", "mol", "xyz", "extxyz", "v3000", "mol2", "cml", "cdxml", "mmcif", "pdb"}
    if corpus.get("schema_version") != 1 or set(corpus.get("cases", {})) != expected_formats:
        print("streaming safety corpus has an invalid schema or format set", file=sys.stderr)
        return 1
    malformed = corpus["cases"]
    if any(
        not isinstance(cases, list) or len(cases) != 12 or any(not isinstance(case, str) for case in cases)
        for cases in malformed.values()
    ):
        print("streaming safety corpus must contain exactly twelve string cases for every format", file=sys.stderr)
        return 1

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
        "extxyz": FIXTURES / "streaming.extxyz",
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
        malformed_count = 0
        for fmt, cases in malformed.items():
            for case_index, content in enumerate(cases):
                path = temp / f"malformed-{case_index}.{fmt}"
                path.write_text(content, encoding="utf-8")
                result = run_runner(args.binary, fmt, path, *malformed_options.get(fmt, ()))
                malformed_count += 1
                require(
                    result["records"] == 0 and result["failures"] == 1,
                    f"{fmt} malformed case {case_index} did not fail exactly once: {result}",
                    errors,
                )

        for fmt, path in valid.items():
            result = run_runner(args.binary, fmt, path, "--max-input-bytes", "1")
            require(
                result["records"] == 0 and result["failures"] == 1,
                f"{fmt} oversized input was not rejected exactly once: {result}",
                errors,
            )

        expected_records = {
            "sdf": 2,
            "mol": 2,
            "xyz": 2,
            "extxyz": 2,
            "v3000": 1,
            "mol2": 1,
            "cml": 1,
            "cdxml": 1,
            "mmcif": 1,
            "pdb": 1,
        }
        gzip_cases = 0
        for fmt, source in valid.items():
            gzip_path = temp / f"{fmt}.gz"
            with gzip.open(gzip_path, "wb") as output:
                output.write(source.read_bytes())
            result = run_runner(args.binary, fmt, gzip_path, "--gzip")
            gzip_cases += 1
            require(
                result["records"] == expected_records[fmt] and result["failures"] == 0,
                f"gzip {fmt} control case did not parse expected records: {result}",
                errors,
            )
            result = run_runner(args.binary, fmt, gzip_path, "--gzip", "--max-input-bytes", "1")
            gzip_cases += 1
            require(
                result["records"] == 0 and result["failures"] == 1,
                f"gzip {fmt} decompressed input limit was not enforced: {result}",
                errors,
            )

    if errors:
        print("streaming format limit failures:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(f"streaming format limits OK: {malformed_count} negative, 10 oversized, {gzip_cases} gzip cases")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
