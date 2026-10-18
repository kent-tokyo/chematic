#!/usr/bin/env python3
"""Classify CTfile readback residuals against a pinned RDKit release."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import subprocess
from collections import Counter
from pathlib import Path

import chematic
import rdkit
from rdkit import Chem, RDLogger

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "validation/cosmolkit_comparison/corpora/stress-v1.smi"


def canonical(mol: Chem.Mol, *, stereo: bool = True) -> str:
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=stereo)


def topology_key(mol: Chem.Mol) -> str:
    copy = Chem.RWMol(mol)
    for atom in copy.GetAtoms():
        atom.SetFormalCharge(0)
        atom.SetIsotope(0)
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        atom.SetIsAromatic(False)
    for bond in copy.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
        bond.SetStereo(Chem.BondStereo.STEREONONE)
    return Chem.MolToSmiles(copy, canonical=True, isomericSmiles=False)


def classify(expected: Chem.Mol, actual: Chem.Mol) -> str:
    if canonical(expected, stereo=False) == canonical(actual, stereo=False):
        return "stereo_only"
    if topology_key(expected) == topology_key(actual):
        expected_atoms = [(a.GetAtomicNum(), a.GetFormalCharge()) for a in expected.GetAtoms()]
        actual_atoms = [(a.GetAtomicNum(), a.GetFormalCharge()) for a in actual.GetAtoms()]
        return "charge_or_bond_order" if expected_atoms != actual_atoms else "bond_order_or_aromaticity"
    return "connectivity_or_atom_identity"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-rdkit", default="2026.09.1")
    args = parser.parse_args()
    if rdkit.__version__ != args.expected_rdkit:
        raise RuntimeError(f"expected RDKit {args.expected_rdkit}; got {rdkit.__version__}")
    RDLogger.DisableLog("rdApp.*")
    rows = [
        line.split()[0] for line in CORPUS.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    counts: Counter[str] = Counter()
    residuals = []
    for index, text in enumerate(rows):
        reference = Chem.MolFromSmiles(text)
        try:
            source = chematic.from_smiles(text)
            source_smiles = source.rdkit_smiles
        except Exception as exc:
            counts["input_refused"] += 1
            residuals.append({"row": index, "input": text, "classification": "input_refused",
                              "error": f"{type(exc).__name__}: {exc}"})
            continue
        if reference is None or source_smiles != canonical(reference):
            counts["input_incompatible"] += 1
            residuals.append({"row": index, "input": text, "classification": "input_incompatible",
                              "chematic_input": source_smiles,
                              "rdkit_input": None if reference is None else canonical(reference)})
            continue
        try:
            reference_block = Chem.MolToMolBlock(reference)
        except Exception as exc:
            counts["oracle_write_failed"] += 1
            residuals.append({
                "row": index,
                "input": text,
                "classification": "oracle_write_failed",
                "error": f"{type(exc).__name__}: {exc}",
            })
            continue
        expected = Chem.MolFromMolBlock(reference_block)
        if expected is None:
            counts["oracle_readback_failed"] += 1
            residuals.append({"row": index, "input": text, "classification": "oracle_readback_failed"})
            continue
        try:
            try:
                block = source.rdkit_mol_block_2d()
                writer = "rdkit_mol_block_2d"
            except (AttributeError, ValueError):
                block = source.to_mol_block()
                writer = "to_mol_block_fallback"
        except Exception as exc:
            counts["writer_refused"] += 1
            residuals.append({"row": index, "input": text, "classification": "writer_refused",
                              "error": f"{type(exc).__name__}: {exc}"})
            continue
        actual = Chem.MolFromMolBlock(block)
        if actual is None:
            counts["readback_failed"] += 1
            residuals.append({"row": index, "input": text, "classification": "readback_failed",
                              "writer": writer, "block_sha256": hashlib.sha256(block.encode()).hexdigest()})
            continue
        expected_smiles = canonical(expected)
        actual_smiles = canonical(actual)
        if expected_smiles == actual_smiles:
            counts["exact"] += 1
            continue
        category = classify(expected, actual)
        counts[category] += 1
        residuals.append({"row": index, "input": text, "classification": category,
                          "writer": writer, "chematic_readback": actual_smiles,
                          "rdkit_readback": expected_smiles})

    extension = Path(importlib.import_module("chematic.chematic").__file__)
    report = {
        "schema": "chematic.rdkit-molblock-readback-parity.v1",
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "versions": {"chematic": chematic.__version__, "rdkit": rdkit.__version__},
        "corpus": {"path": str(CORPUS.relative_to(ROOT)), "rows": len(rows),
                   "sha256": hashlib.sha256(CORPUS.read_bytes()).hexdigest()},
        "native_extension_sha256": hashlib.sha256(extension.read_bytes()).hexdigest(),
        "classification_counts": dict(sorted(counts.items())),
        "residuals": residuals,
        "contract": (
            "Exact RDKit canonical-isomeric readback is agreement. Every other row is retained and "
            "classified; input/oracle failures are not charged to the writer and no row is replaced."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(report["classification_counts"], indent=2))


if __name__ == "__main__":
    main()
