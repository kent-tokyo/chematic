#!/usr/bin/env python3
"""Measure CheMatic Cleanup against RDKit on a frozen stress-corpus slice.

The input-compatibility filter is identical to bench_rdkit_workflows.py. The
historical ``standardize()`` result is retained as a separate baseline because
that API performs transformations outside RDKit Cleanup's contract.
"""

from __future__ import annotations

import argparse
from collections import Counter
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

NON_TETRAHEDRAL_TAGS = {
    Chem.ChiralType.CHI_SQUAREPLANAR,
    Chem.ChiralType.CHI_TRIGONALBIPYRAMIDAL,
    Chem.ChiralType.CHI_OCTAHEDRAL,
}


def smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def chematic_output(operation) -> tuple[str | None, str | None]:
    try:
        return operation().rdkit_smiles, None
    except Exception as exc:
        return None, f"{type(exc).__name__}: {exc}"


def graph_key(text: str | None, *, keep_charge: bool, keep_bond_order: bool) -> str | None:
    """Return a deterministic RDKit graph key for residual classification.

    This is deliberately not an equality oracle for chemistry. It only lets
    the report distinguish charge assignment from connectivity/bond-order
    policy without hiding the original exact strings.
    """
    if text is None:
        return None
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return None
    copy = Chem.RWMol(mol)
    for atom in copy.GetAtoms():
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        atom.SetIsotope(0)
        if not keep_charge:
            atom.SetFormalCharge(0)
    for bond in copy.GetBonds():
        bond.SetStereo(Chem.BondStereo.STEREONONE)
        if not keep_bond_order:
            bond.SetBondType(Chem.BondType.SINGLE)
            bond.SetIsAromatic(False)
    for atom in copy.GetAtoms():
        atom.SetIsAromatic(False)
    return Chem.MolToSmiles(copy, canonical=True, isomericSmiles=False)


def classify_mismatch(
    source: Chem.Mol,
    normalized: str | None,
    expected: str | None,
    actual: str | None,
    cleanup_error: str | None,
) -> tuple[str, str]:
    """Classify one Cleanup residual without treating it as an exact match."""
    if cleanup_error is not None or actual is None:
        if cleanup_error and "kekulize" in cleanup_error.lower():
            return "aromatic_n_oxide_kekulization", "not_implemented"
        return "operation_error", "unsupported_or_invalid_output"

    has_metal = any(atom.GetAtomicNum() in METALS for atom in source.GetAtoms())
    if has_metal:
        if any(atom.GetChiralTag() in NON_TETRAHEDRAL_TAGS for atom in source.GetAtoms()):
            return "metal_stereo_preservation", "policy_difference"
        exact_graph = graph_key(actual, keep_charge=True, keep_bond_order=True)
        expected_graph = graph_key(expected, keep_charge=True, keep_bond_order=True)
        neutral_graph = graph_key(actual, keep_charge=False, keep_bond_order=True)
        expected_neutral_graph = graph_key(expected, keep_charge=False, keep_bond_order=True)
        if exact_graph is not None and exact_graph == expected_graph:
            return "canonical_representation", "representation_only"
        if neutral_graph is not None and neutral_graph == expected_neutral_graph:
            return "metal_charge_assignment", "not_implemented"
        topology = graph_key(actual, keep_charge=False, keep_bond_order=False)
        expected_topology = graph_key(expected, keep_charge=False, keep_bond_order=False)
        if topology is not None and topology == expected_topology:
            return "metal_bond_order_assignment", "not_implemented"
        return "metal_disconnection_policy", "not_implemented"

    if normalized != expected:
        return "acid_base_reionization", "not_implemented"

    if actual is not None and expected is not None:
        if "[S+]" in actual and "[O-]" in actual and "[S+]" not in expected:
            return "sulfur_oxide_representation", "representation_only"
        if "[H]=C" in smiles(source) or "[CH]" in actual:
            return "explicit_hydrogen_normalization", "not_implemented"

    topology = graph_key(actual, keep_charge=False, keep_bond_order=False)
    expected_topology = graph_key(expected, keep_charge=False, keep_bond_order=False)
    if topology is not None and topology == expected_topology:
        return "functional_group_normalization", "not_implemented"
    return "unclassified_semantic_difference", "needs_adjudication"


def classify_reionize_mismatch(
    source: str | None, expected: str | None, actual: str | None, error: str | None
) -> tuple[str, str]:
    if error is not None or actual is None:
        return "operation_error", "unsupported_or_invalid_output"
    if expected == source:
        return "independent_ionization_policy", "policy_difference"
    topology = graph_key(actual, keep_charge=False, keep_bond_order=False)
    expected_topology = graph_key(expected, keep_charge=False, keep_bond_order=False)
    if topology is not None and topology == expected_topology:
        return "acid_base_charge_assignment", "not_implemented"
    return "hydrogen_or_resonance_assignment", "not_implemented"


def classify_normalize_mismatch(
    expected: str | None, actual: str | None, error: str | None
) -> tuple[str, str]:
    if error is not None or actual is None:
        if error and "kekulize" in error.lower():
            return "aromatic_n_oxide_kekulization", "not_implemented"
        return "operation_error", "unsupported_or_invalid_output"
    topology = graph_key(actual, keep_charge=False, keep_bond_order=False)
    expected_topology = graph_key(expected, keep_charge=False, keep_bond_order=False)
    if topology is not None and topology == expected_topology:
        return "functional_group_normalization", "not_implemented"
    return "normalization_graph_assignment", "not_implemented"


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
                ROOT / "crates/chematic-chem/src/standardize/cleanup.rs",
                ROOT / "crates/chematic-chem/src/standardize/groups.rs",
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
    reionize_exact = 0
    rdkit_changed = 0
    excluded = []
    oracle_errors = []
    operation_errors = []
    mismatches = []
    classifications: Counter[str] = Counter()
    support_statuses: Counter[str] = Counter()
    reionize_classifications: Counter[str] = Counter()
    reionize_support_statuses: Counter[str] = Counter()
    reionize_mismatches = []
    normalize_classifications: Counter[str] = Counter()
    normalize_support_statuses: Counter[str] = Counter()
    normalize_mismatches = []
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
        try:
            expected_reionized = smiles(rdMolStandardize.Reionize(rmol))
        except Exception:
            expected_reionized = None
        legacy, legacy_error = chematic_output(cmol.standardize)
        normalizer_only, normalize_error = chematic_output(cmol.normalize_groups)
        actual, cleanup_error = chematic_output(cmol.cleanup)
        actual_reionized, reionize_error = chematic_output(cmol.reionize)
        for operation, error in (
            ("standardize", legacy_error),
            ("normalize_groups", normalize_error),
            ("cleanup", cleanup_error),
            ("reionize", reionize_error),
        ):
            if error is not None:
                operation_errors.append({"sample_index": index,
                                         "operation": operation, "error": error})
        legacy_exact += legacy == expected
        cleanup_exact += actual == expected
        normalizer_only_exact += normalizer_only == expected
        normalize_stage_exact += normalizer_only == normalized
        reionize_exact += actual_reionized == expected_reionized
        rdkit_changed += expected != r_input
        if normalizer_only != normalized:
            classification, support_status = classify_normalize_mismatch(
                normalized, normalizer_only, normalize_error
            )
            normalize_classifications[classification] += 1
            normalize_support_statuses[support_status] += 1
            normalize_mismatches.append({
                "sample_index": index,
                "input": text,
                "rdkit_normalize": normalized,
                "chematic_normalize": normalizer_only,
                "chematic_error": normalize_error,
                "classification": classification,
                "support_status": support_status,
            })
        if actual != expected:
            classification, support_status = classify_mismatch(
                rmol, normalized, expected, actual, cleanup_error
            )
            classifications[classification] += 1
            support_statuses[support_status] += 1
            mismatches.append({
                "sample_index": index,
                "input": text,
                "rdkit_input": r_input,
                "rdkit_cleanup": expected,
                "chematic_cleanup": actual,
                "chematic_error": cleanup_error,
                "classification": classification,
                "support_status": support_status,
            })
        if actual_reionized != expected_reionized:
            classification, support_status = classify_reionize_mismatch(
                r_input, expected_reionized, actual_reionized, reionize_error
            )
            reionize_classifications[classification] += 1
            reionize_support_statuses[support_status] += 1
            reionize_mismatches.append({
                "sample_index": index,
                "input": text,
                "rdkit_reionize": expected_reionized,
                "chematic_reionize": actual_reionized,
                "chematic_error": reionize_error,
                "classification": classification,
                "support_status": support_status,
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
        "normalize_classification_counts": dict(sorted(normalize_classifications.items())),
        "normalize_support_status_counts": dict(sorted(normalize_support_statuses.items())),
        "normalize_mismatches": normalize_mismatches,
        "cleanup_exact": cleanup_exact,
        "classification_counts": dict(sorted(classifications.items())),
        "support_status_counts": dict(sorted(support_statuses.items())),
        "mismatches": mismatches,
        "reionize_exact": reionize_exact,
        "reionize_classification_counts": dict(sorted(reionize_classifications.items())),
        "reionize_support_status_counts": dict(sorted(reionize_support_statuses.items())),
        "reionize_mismatches": reionize_mismatches,
    })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(f"legacy standardize: {legacy_exact}/{compatible}")
    print(f"normalizer only: {normalizer_only_exact}/{compatible}")
    print(f"normalize vs normalize: {normalize_stage_exact}/{compatible}")
    print(f"cleanup: {cleanup_exact}/{compatible}")
    print(f"mismatches: {len(mismatches)}")
    print(f"reionize: {reionize_exact}/{compatible}")
    print(f"reionize mismatches: {len(reionize_mismatches)}")


if __name__ == "__main__":
    main()
