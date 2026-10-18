#!/usr/bin/env python3
"""Measure CheMatic Cleanup against RDKit on a frozen stress-corpus slice.

The input-compatibility filter is identical to bench_rdkit_workflows.py. The
historical ``standardize()`` result is retained as a separate baseline because
that API performs transformations outside RDKit Cleanup's contract.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import subprocess
from pathlib import Path

import chematic
import rdkit
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "validation/cosmolkit_comparison/corpora/stress-v1.smi"
METALS = {
    3, 4, 11, 12, 13, 19, 20, *range(21, 32), *range(37, 51),
    *range(55, 84), *range(87, 117),
}


def smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def chematic_output(operation) -> tuple[str | None, str | None]:
    try:
        return operation().rdkit_smiles, None
    except Exception as exc:
        return None, f"{type(exc).__name__}: {exc}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sample-size", type=int, default=256)
    args = parser.parse_args()
    if args.sample_size < 2:
        parser.error("--sample-size must be at least 2")
    if rdkit.__version__ != "2026.09.1":
        raise RuntimeError(f"expected RDKit 2026.09.1; got {rdkit.__version__}")
    RDLogger.DisableLog("rdApp.*")

    all_rows = [
        line.split()[0]
        for line in CORPUS.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    sampled = (
        all_rows
        if args.sample_size >= len(all_rows)
        else [
            all_rows[round(i * (len(all_rows) - 1) / (args.sample_size - 1))]
            for i in range(args.sample_size)
        ]
    )
    native_path = Path(importlib.import_module("chematic.chematic").__file__)
    result: dict[str, object] = {
        "source_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "source_sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (
                ROOT / "crates/chematic-chem/src/standardize.rs",
                ROOT / "crates/chematic-py/src/mol_methods/transformations.rs",
            )
        },
        "corpus": str(CORPUS.relative_to(ROOT)),
        "corpus_sha256": hashlib.sha256(CORPUS.read_bytes()).hexdigest(),
        "selection": (
            f"{len(sampled)} evenly spaced rows, then identical parsed "
            "RDKit-compatible SMILES"
        ),
        "versions": {"chematic": chematic.__version__, "rdkit": rdkit.__version__},
        "native_extension_sha256": hashlib.sha256(native_path.read_bytes()).hexdigest(),
        "source_rows": len(all_rows),
        "sampled_rows": len(sampled),
    }

    compatible = 0
    input_compatible = 0
    legacy_exact = 0
    cleanup_exact = 0
    normalizer_only_exact = 0
    normalize_stage_exact = 0
    rdkit_changed = 0
    excluded = []
    oracle_errors = []
    operation_errors = []
    mismatches = []
    for index, text in enumerate(sampled):
        try:
            cmol = chematic.from_smiles(text)
            c_input = cmol.rdkit_smiles
        except Exception:
            cmol = None
            c_input = None
        rmol = Chem.MolFromSmiles(text)
        r_input = smiles(rmol) if rmol is not None else None
        if cmol is None or c_input != r_input:
            excluded.append({"sample_index": index, "input": text,
                             "chematic": c_input, "rdkit": r_input})
            continue

        input_compatible += 1
        try:
            expected = smiles(rdMolStandardize.Cleanup(rmol))
        except Exception as exc:
            oracle_errors.append({"sample_index": index, "input": text,
                                  "error": f"{type(exc).__name__}: {exc}"})
            continue
        compatible += 1
        try:
            normalized = smiles(rdMolStandardize.Normalize(rmol))
        except Exception:
            normalized = None
        legacy, legacy_error = chematic_output(cmol.standardize)
        normalizer_only, normalize_error = chematic_output(cmol.normalize_groups)
        actual, cleanup_error = chematic_output(cmol.cleanup)
        for operation, error in (
            ("standardize", legacy_error),
            ("normalize_groups", normalize_error),
            ("cleanup", cleanup_error),
        ):
            if error is not None:
                operation_errors.append({"sample_index": index,
                                         "operation": operation, "error": error})
        legacy_exact += legacy == expected
        cleanup_exact += actual == expected
        normalizer_only_exact += normalizer_only == expected
        normalize_stage_exact += normalizer_only == normalized
        rdkit_changed += expected != r_input
        if actual != expected:
            metal = any(atom.GetAtomicNum() in METALS for atom in rmol.GetAtoms())
            mismatches.append({
                "sample_index": index,
                "input": text,
                "rdkit_input": r_input,
                "rdkit_cleanup": expected,
                "chematic_cleanup": actual,
                "chematic_error": cleanup_error,
                "category": "metal" if metal else (
                    "reionization" if normalized != expected else "other"
                ),
            })

    result.update({
        "input_compatible_rows": input_compatible,
        "compatible_rows": compatible,
        "excluded": excluded,
        "oracle_errors": oracle_errors,
        "operation_errors": operation_errors,
        "rdkit_cleanup_changed_rows": rdkit_changed,
        "legacy_standardize_exact": legacy_exact,
        "normalizer_only_vs_cleanup_exact": normalizer_only_exact,
        "normalize_stage_exact": normalize_stage_exact,
        "cleanup_exact": cleanup_exact,
        "mismatches": mismatches,
    })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(f"legacy standardize: {legacy_exact}/{compatible}")
    print(f"normalizer only: {normalizer_only_exact}/{compatible}")
    print(f"normalize vs normalize: {normalize_stage_exact}/{compatible}")
    print(f"cleanup: {cleanup_exact}/{compatible}")
    print(f"mismatches: {len(mismatches)}")


if __name__ == "__main__":
    main()
