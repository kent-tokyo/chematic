#!/usr/bin/env python3
"""Score native SMILES read-back against pinned RDKit and COSMolKit.

Build the current source first:
  cargo build -p chematic-smiles --example canonical_roundtrip_dump --offline

Then run with a pinned RDKit version and COSMolKit 0.5.0:
  python check_native_smiles_roundtrip.py --corpus corpora/stress-v1.smi \
      --chematic ../../target/debug/examples/canonical_roundtrip_dump \
      --expected-rdkit 2025.09.3

The oracle is RDKit's canonical isomeric SMILES of each source molecule.
Each engine's native output is parsed by RDKit before comparison. Report both
strict read-back and verified equivalence: exact stereochemical SMILES after
renumbering up to 10,000 graph matches, then identical unsanitized aromatic
input graphs with matching formulas and stereo-aware canonical SMILES, followed
by a bounded resonance search for remaining molecules without specified
stereochemistry. Other differences remain unresolved.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import cosmolkit
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors


def has_specified_stereo(mol: Chem.Mol) -> bool:
    atom_stereo = any(
        atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for atom in mol.GetAtoms()
    )
    bond_stereo = any(
        bond.GetStereo() != Chem.BondStereo.STEREONONE for bond in mol.GetBonds()
    )
    return atom_stereo or bond_stereo


def resonance_smiles(mol: Chem.Mol) -> set[str]:
    flags = (
        Chem.ResonanceFlags.ALLOW_CHARGE_SEPARATION
        | Chem.ResonanceFlags.ALLOW_INCOMPLETE_OCTETS
        | Chem.ResonanceFlags.UNCONSTRAINED_CATIONS
        | Chem.ResonanceFlags.UNCONSTRAINED_ANIONS
        | Chem.ResonanceFlags.KEKULE_ALL
    )
    supplier = Chem.ResonanceMolSupplier(mol, flags=flags, maxStructs=256)
    return {Chem.MolToSmiles(candidate) for candidate in supplier if candidate is not None}


def stereo_isomorphic(expected: Chem.Mol, actual: Chem.Mol) -> bool:
    if expected.GetNumAtoms() != actual.GetNumAtoms() or expected.GetNumBonds() != actual.GetNumBonds():
        return False
    reference = Chem.MolToSmiles(expected, canonical=False)
    # RDKit's useChirality matching does not distinguish every non-tetrahedral
    # permutation (for example @SP1 versus @SP2). Compare each graph mapping
    # with the same atom order instead. Noncanonical SMILES then encodes the
    # local neighbor order and its tetrahedral/non-tetrahedral stereo tag.
    for match in expected.GetSubstructMatches(
        actual, uniquify=False, useChirality=False, maxMatches=10000
    ):
        order = [0] * len(match)
        for actual_index, expected_index in enumerate(match):
            order[expected_index] = actual_index
        aligned = Chem.RenumberAtoms(actual, order)
        if Chem.MolToSmiles(aligned, canonical=False) == reference:
            return True
    return False


def same_raw_aromatic_graph(
    expected_text: str, actual_text: str, expected: Chem.Mol, actual: Chem.Mol
) -> bool:
    """Check the full stereo-aware aromatic graph before Kekulization."""
    if rdMolDescriptors.CalcMolFormula(expected) != rdMolDescriptors.CalcMolFormula(actual):
        return False
    raw_expected = Chem.MolFromSmiles(expected_text, sanitize=False)
    raw_actual = Chem.MolFromSmiles(actual_text, sanitize=False)
    if raw_expected is None or raw_actual is None:
        return False
    if not any(atom.GetIsAromatic() for atom in raw_expected.GetAtoms()):
        return False
    # Unsanitized isomeric SMILES includes atom and bond stereo, so equality
    # cannot waive an inversion elsewhere in the same molecule.
    return Chem.MolToSmiles(raw_expected, isomericSmiles=True) == Chem.MolToSmiles(
        raw_actual, isomericSmiles=True
    )


def classify(
    expected: Chem.Mol,
    actual: Chem.Mol | None,
    expected_text: str | None = None,
    actual_text: str | None = None,
) -> str:
    if actual is None:
        return "invalid"
    if Chem.MolToSmiles(expected) == Chem.MolToSmiles(actual):
        return "strict"
    if stereo_isomorphic(expected, actual):
        return "stereo_isomorphic"
    if expected_text is not None and actual_text is not None:
        if same_raw_aromatic_graph(expected_text, actual_text, expected, actual):
            return "raw_aromatic_graph"
    # The resonance supplier may change bond orders and charges. Do not use it
    # to waive a mismatch involving specified atom or bond stereochemistry.
    if not has_specified_stereo(expected) and not has_specified_stereo(actual):
        if resonance_smiles(expected) & resonance_smiles(actual):
            return "resonance"
    return "unresolved"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--corpus", required=True, type=Path)
    parser.add_argument("--chematic", required=True, type=Path)
    parser.add_argument("--expected-rdkit", default="2025.09.3")
    args = parser.parse_args()
    if Chem.rdBase.rdkitVersion != args.expected_rdkit or cosmolkit.__version__ != "0.5.0":
        raise SystemExit(f"requires RDKit {args.expected_rdkit} and COSMolKit 0.5.0")
    RDLogger.DisableLog("rdApp.*")
    corpus = args.corpus.resolve()
    binary = args.chematic.resolve()
    lines = corpus.read_text().splitlines()
    completed = subprocess.run([str(binary), str(corpus)], text=True, capture_output=True, check=True)
    output = {}
    for line in completed.stdout.splitlines():
        index, smiles = line.split("\t", 1)
        output[int(index)] = smiles

    eligible = 0
    skipped = 0
    matched = {"chematic": 0, "cosmolkit": 0}
    labels = ("strict", "stereo_isomorphic", "raw_aromatic_graph", "resonance", "unresolved", "invalid")
    classifications = {engine: {label: 0 for label in labels} for engine in matched}
    misses = {"chematic": [], "cosmolkit": []}
    for index, line in enumerate(lines):
        tokens = line.split()
        if not tokens:
            continue
        source = tokens[0]
        expected_mol = Chem.MolFromSmiles(source)
        if expected_mol is None:
            skipped += 1
            continue
        expected = Chem.MolToSmiles(expected_mol)
        eligible += 1
        try:
            cosmolkit_smiles = cosmolkit.Molecule.from_smiles(source).to_smiles()
        except Exception:
            cosmolkit_smiles = None
        for engine, smiles in (
            ("chematic", output.get(index)),
            ("cosmolkit", cosmolkit_smiles),
        ):
            actual_mol = Chem.MolFromSmiles(smiles) if smiles and not smiles.startswith("ERROR:") else None
            actual = Chem.MolToSmiles(actual_mol) if actual_mol is not None else None
            classification = classify(expected_mol, actual_mol, source, smiles)
            classifications[engine][classification] += 1
            if classification == "strict":
                matched[engine] += 1
            else:
                misses[engine].append({
                    "index": index, "classification": classification,
                    "expected": expected, "actual": actual,
                })
    print(json.dumps({
        "corpus_sha256": hashlib.sha256(corpus.read_bytes()).hexdigest(),
        "chematic_binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
        "rdkit_version": Chem.rdBase.rdkitVersion,
        "cosmolkit_version": cosmolkit.__version__,
        "eligible": eligible,
        "skipped": skipped,
        "matched": matched,
        "classifications": classifications,
        "verified_equivalent": {
            engine: counts["strict"] + counts["stereo_isomorphic"]
            + counts["raw_aromatic_graph"] + counts["resonance"]
            for engine, counts in classifications.items()
        },
        "unresolved_indices": {
            engine: [row["index"] for row in rows if row["classification"] == "unresolved"]
            for engine, rows in misses.items()
        },
        "mismatches": {engine: rows[:20] for engine, rows in misses.items()},
    }, sort_keys=True))


if __name__ == "__main__":
    main()
