#!/usr/bin/env python3
"""Compare same-format rewrites through one semantic observer.

Both chematic and Open Babel outputs are parsed by the same Rust probe. This
avoids declaring one engine correct merely because it reproduces its own text.
The checked-in fixtures are deliberately small; this gate proves those fixture
lanes only, not broad format parity.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path
from typing import Any

try:
    from benchmark_version import workspace_version
except ModuleNotFoundError:  # imported as scripts.* by pytest
    from scripts.benchmark_version import workspace_version


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = {
    "v3000": ROOT / "benchmarks/fixtures/ethanol.v3000",
    "mol2": ROOT / "benchmarks/fixtures/ethanol.mol2",
    "cml": ROOT / "benchmarks/fixtures/ethanol.cml",
    "cdxml": ROOT / "benchmarks/fixtures/ethanol.cdxml",
}


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, cwd=ROOT, check=False, text=True, capture_output=True)


def snapshot(probe: Path, fmt: str, path: Path) -> dict[str, Any]:
    completed = run([str(probe), "--format", fmt, "--path", str(path)])
    if completed.returncode != 0:
        raise RuntimeError(f"semantic probe failed for {path}: {completed.stderr.strip()}")
    result = json.loads(completed.stdout)
    if not isinstance(result, dict):
        raise RuntimeError("semantic probe returned a non-object")
    return result


def rewrite_with_chematic(probe: Path, fmt: str, source: Path, output: Path) -> None:
    completed = run(
        [
            str(probe),
            "--format",
            fmt,
            "--path",
            str(source),
            "--rewrite-output",
            str(output),
        ]
    )
    if completed.returncode != 0:
        raise RuntimeError(f"chematic rewrite failed for {fmt}: {completed.stderr.strip()}")


def rewrite_with_openbabel(executable: str, fmt: str, source: Path, output: Path) -> None:
    if fmt == "v3000":
        command = [executable, "-imol", str(source), "-omol", "-x3", "-O", str(output)]
    else:
        command = [executable, f"-i{fmt}", str(source), f"-o{fmt}", "-O", str(output)]
    completed = run(command)
    if completed.returncode != 0:
        raise RuntimeError(f"Open Babel rewrite failed for {fmt}: {completed.stderr.strip()}")


def coordinate_comparison(
    reference: list[Any], candidate: list[Any], *, allow_uniform_scale: bool
) -> tuple[float | None, float | None]:
    if len(reference) != len(candidate):
        return None, None
    if not reference:
        return 0.0, 1.0
    scale = 1.0
    if allow_uniform_scale:
        ratios = [
            float(got) / float(want)
            for want, got in zip(reference, candidate)
            if abs(float(want)) > 1e-12
        ]
        if ratios:
            ratios.sort()
            scale = ratios[len(ratios) // 2]
    if abs(scale) <= 1e-12:
        return None, scale
    return (
        max(
            abs(float(want) - float(got) / scale)
            for want, got in zip(reference, candidate)
        ),
        scale,
    )


def selected_metadata(fmt: str, metadata: dict[str, Any]) -> Any:
    if fmt == "v3000":
        return {
            "sgroups": metadata.get("sgroups"),
            "bond_properties": metadata.get("bond_properties"),
            "atom_properties": metadata.get("atom_properties"),
        }
    if fmt == "mol2":
        atoms = metadata.get("atoms", [])
        bonds = metadata.get("bonds", [])
        return {
            "molecule_type": metadata.get("molecule_type"),
            "charge_type": metadata.get("charge_type"),
            "atoms": [
                {
                    "type": atom.get("type"),
                    "subst_id": atom.get("subst_id"),
                    "subst_name": atom.get("subst_name"),
                    "partial_charge": atom.get("partial_charge"),
                    "status_bits": atom.get("status_bits"),
                }
                for atom in atoms
            ],
            "bonds": bonds,
            "duplicate_bonds": metadata.get("duplicate_bonds"),
            "unity_atom_attributes": metadata.get("unity_atom_attributes"),
            "opaque_sections": metadata.get("opaque_sections"),
        }
    return {}


def compare(fmt: str, reference: dict[str, Any], candidate: dict[str, Any]) -> dict[str, Any]:
    ref_records = reference.get("records", [])
    candidate_records = candidate.get("records", [])
    errors: list[str] = []
    max_delta = 0.0
    coordinate_scales: list[float] = []
    if len(ref_records) != len(candidate_records):
        errors.append(f"record count {len(candidate_records)} != {len(ref_records)}")
    for index, (ref, got) in enumerate(zip(ref_records, candidate_records)):
        for field in (
            "graph_key",
            "atom_count",
            "bond_count",
            "total_formal_charge",
            "isotopes",
        ):
            if got.get(field) != ref.get(field):
                errors.append(f"record {index}: {field} differs")
        delta, scale = coordinate_comparison(
            ref.get("pairwise_distances", []),
            got.get("pairwise_distances", []),
            allow_uniform_scale=fmt == "cdxml",
        )
        if delta is None:
            errors.append(f"record {index}: coordinate pair count differs")
        else:
            max_delta = max(max_delta, delta)
            if scale is not None:
                coordinate_scales.append(scale)
            if delta > 1e-3:
                errors.append(f"record {index}: coordinate distance delta {delta:.6g} > 0.001")
        if selected_metadata(fmt, got.get("metadata", {})) != selected_metadata(
            fmt, ref.get("metadata", {})
        ):
            errors.append(f"record {index}: format metadata differs")
    return {
        "status": "exact_semantic_match" if not errors else "semantic_mismatch",
        "records": len(candidate_records),
        "max_pairwise_distance_delta": max_delta,
        "coordinate_scales": coordinate_scales,
        "errors": errors,
    }


def executable_version(executable: str) -> str:
    completed = run([executable, "-V"])
    output = (completed.stdout or completed.stderr).strip()
    return output.splitlines()[0] if output else "unknown"


def source_provenance() -> dict[str, Any]:
    revision = run(["git", "rev-parse", "HEAD"])
    status = run(["git", "status", "--porcelain", "--untracked-files=normal"])
    if revision.returncode != 0 or status.returncode != 0:
        raise RuntimeError("cannot determine source provenance")
    return {"commit": revision.stdout.strip(), "dirty": bool(status.stdout.strip())}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--openbabel", default="obabel")
    parser.add_argument(
        "--probe",
        type=Path,
        default=ROOT / "target/release/examples/file_io_semantic_probe",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    if not args.probe.exists():
        raise SystemExit(
            f"semantic probe not found: {args.probe}; run "
            "cargo build -p chematic-mol --example file_io_semantic_probe --release"
        )

    formats: dict[str, Any] = {}
    chematic_errors: list[str] = []
    wins = 0
    ties = 0
    losses = 0
    with tempfile.TemporaryDirectory(prefix="chematic-openbabel-semantic-") as temp:
        temp_dir = Path(temp)
        for fmt, fixture in FIXTURES.items():
            baseline = snapshot(args.probe, fmt, fixture)
            chematic_path = temp_dir / f"chematic.{fmt}"
            openbabel_path = temp_dir / f"openbabel.{fmt}"
            rewrite_with_chematic(args.probe, fmt, fixture, chematic_path)
            rewrite_with_openbabel(args.openbabel, fmt, fixture, openbabel_path)
            chematic_result = compare(fmt, baseline, snapshot(args.probe, fmt, chematic_path))
            openbabel_result = compare(fmt, baseline, snapshot(args.probe, fmt, openbabel_path))
            formats[fmt] = {
                "fixture": str(fixture.relative_to(ROOT)),
                "fixture_sha256": hashlib.sha256(fixture.read_bytes()).hexdigest(),
                "chematic": chematic_result,
                "openbabel": openbabel_result,
            }
            chematic_exact = chematic_result["status"] == "exact_semantic_match"
            openbabel_exact = openbabel_result["status"] == "exact_semantic_match"
            if not chematic_exact:
                chematic_errors.append(f"{fmt}/chematic: {chematic_result['errors']}")
            if chematic_exact and not openbabel_exact:
                wins += 1
                formats[fmt]["result"] = "chematic"
            elif chematic_exact == openbabel_exact:
                ties += 1
                formats[fmt]["result"] = "tie"
            else:
                losses += 1
                formats[fmt]["result"] = "openbabel"

    report = {
        "schema_version": 1,
        "target_version": workspace_version(ROOT),
        "source": source_provenance(),
        "gate": "openbabel_file_io_fixture_semantics_v1",
        "status": "local-verified" if not chematic_errors else "failed",
        "observer": "the same chematic Rust semantic probe parses both rewritten outputs",
        "coordinate_tolerance": 1e-3,
        "score": {"chematic_wins": wins, "ties": ties, "openbabel_wins": losses},
        "formats": formats,
        "tool_versions": {"openbabel": executable_version(args.openbabel)},
        "claim_boundary": (
            "checked-in single-molecule fixture graph, charge, isotope, bond/stereo key, "
            "coordinates, and bounded format metadata only; not broad corpus parity"
        ),
        "errors": chematic_errors,
    }
    encoded = json.dumps(report, indent=2) + "\n"
    if args.output:
        target = args.output if args.output.is_absolute() else ROOT / args.output
        target.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    return 1 if chematic_errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
