#!/usr/bin/env python3
"""Compare the potential tetrahedral stereocenter count with RDKit.

The chematic API intentionally documents ``num_stereocenters`` as the
potential-center count, not the count of assigned R/S labels.  RDKit's
``CalcNumAtomStereoCenters`` is the matching oracle.  Unspecified centers and
assigned centers are therefore not silently mixed with ``CalcNumUnspecified``
or CIP-label counts.
"""

import hashlib
import argparse
import json
import sys
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
CORPUS = ROOT / "scripts" / "descriptor_census_corpus.smi"
DEFAULT_OUT = ROOT / "validation" / "results" / "descriptor-stereocenter-rdkit-parity-v1.0.13.json"
# The comparator pin used by the current accuracy manifest.  A caller may
# override this explicitly, but the default must agree with the documented
# and recorded local evidence rather than silently probing whatever happens
# to be installed.
EXPECTED_RDKIT = "2025.09.3"


def rdkit_atom_stereocenter_indices(mol):
    """Recover the atom-level set behind RDKit's count API.

    ``FindPotentialStereo(flagPossible=True)`` also reports planar or
    otherwise non-assignable candidates.  Probe each candidate independently
    with a temporary tetrahedral tag and retain only atoms for which RDKit can
    assign a CIP code.  The resulting set is checked against
    ``CalcNumAtomStereoCenters`` rather than being used as a replacement for
    that count oracle.
    """
    indices = set()
    for info in Chem.FindPotentialStereo(mol, cleanIt=False, flagPossible=True):
        if info.type != Chem.StereoType.Atom_Tetrahedral:
            continue
        probe = Chem.Mol(mol)
        probe.GetAtomWithIdx(int(info.centeredOn)).SetChiralTag(
            Chem.ChiralType.CHI_TETRAHEDRAL_CW
        )
        Chem.AssignStereochemistry(probe, force=True, cleanIt=False)
        atom = probe.GetAtomWithIdx(int(info.centeredOn))
        if atom.HasProp("_CIPCode"):
            indices.add(int(info.centeredOn))
    return indices


def rdkit_potential_tetrahedral_indices(mol):
    """Return RDKit's explicit potential tetrahedral atom set.

    This is a fallback for structures where assigning one arbitrary temporary
    tetrahedral tag does not produce a CIP code, even though RDKit's count
    oracle still includes the atom.  The caller accepts it only when its size
    agrees with ``CalcNumAtomStereoCenters``; otherwise the row remains
    unresolved instead of being silently scored.
    """
    return {
        int(info.centeredOn)
        for info in Chem.FindPotentialStereo(mol, cleanIt=False, flagPossible=True)
        if info.type == Chem.StereoType.Atom_Tetrahedral
    }


def python_provenance(chematic):
    """Record the exact Python extension used by the parity run."""
    package_init = Path(chematic.__file__).resolve()
    extension = next(package_init.parent.glob("chematic*.so"), None)
    return {
        "executable": str(Path(sys.executable).resolve()),
        "package_version": getattr(chematic, "__version__", None),
        "package_init": str(package_init),
        "extension": str(extension) if extension else None,
        "extension_sha256": hashlib.sha256(extension.read_bytes()).hexdigest()
        if extension
        else None,
        "pythonpath": sys.path,
        "workspace_path_visible": str(ROOT) in str(package_init),
        "import_origin": (
            "temporary_extracted_wheel"
            if "site-packages" not in str(package_init)
            else "site_package"
        ),
    }


def load_corpus_rows(corpus: Path) -> list[str]:
    """Load plain SMILES or JSONL records carrying a ``smiles`` field.

    CIP investigations keep provenance and bucket labels in JSONL.  Treating
    those records as literal SMILES silently turns every row into a parse
    failure, so the runner makes the input shape explicit and skips only the
    optional manifest record.
    """
    rows: list[str] = []
    for line_number, line in enumerate(corpus.read_text(encoding="utf-8").splitlines(), 1):
        value = line.strip()
        if not value:
            continue
        if corpus.suffix.lower() == ".jsonl":
            try:
                record = json.loads(value)
            except json.JSONDecodeError as exc:
                raise ValueError(f"{corpus}:{line_number}: invalid JSONL: {exc}") from exc
            if not isinstance(record, dict):
                raise ValueError(f"{corpus}:{line_number}: JSONL row must be an object")
            if record.get("_manifest") is True:
                continue
            smiles = record.get("smiles")
            if not isinstance(smiles, str) or not smiles.strip():
                raise ValueError(f"{corpus}:{line_number}: row is missing a non-empty smiles field")
            rows.append(smiles.strip())
        else:
            rows.append(value)
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument("--expected-rdkit", default=EXPECTED_RDKIT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    if rdBase.rdkitVersion != args.expected_rdkit:
        print(f"FATAL: expected rdkit=={args.expected_rdkit}, got {rdBase.rdkitVersion}", file=sys.stderr)
        return 2
    corpus = args.corpus.resolve()
    if not corpus.exists():
        print(f"FATAL: missing corpus {corpus}", file=sys.stderr)
        return 2

    raw = corpus.read_bytes()
    try:
        rows = load_corpus_rows(corpus)
    except ValueError as exc:
        print(f"FATAL: {exc}", file=sys.stderr)
        return 2
    compared = 0
    rdkit_parse_failures = 0
    chematic_failures = 0
    mismatch_examples = []
    mismatch_count = 0
    atom_false_positive_count = 0
    atom_false_negative_count = 0
    atom_mismatch_examples = []
    error_examples = []
    error_count = 0
    atom_oracle_fallback_count = 0
    atom_oracle_unresolved_count = 0
    # Import after the environment check so this gate cannot accidentally use
    # a stale global installation when invoked from the source-built venv.
    try:
        import chematic
    except Exception as exc:  # pragma: no cover - environment diagnostic
        print(f"FATAL: could not import source-built chematic: {exc}", file=sys.stderr)
        return 2
    provenance = python_provenance(chematic)

    for index, smi in enumerate(rows):
        rd_mol = Chem.MolFromSmiles(smi)
        if rd_mol is None:
            rdkit_parse_failures += 1
            continue
        try:
            cm = chematic.from_smiles(smi)
            chematic_value = int(cm.num_stereocenters)
            chematic_indices = set(int(i) for i in cm.potential_stereocenter_indices)
        except Exception as exc:
            chematic_failures += 1
            error_count += 1
            if len(error_examples) < 20:
                error_examples.append({"index": index, "smiles": smi, "error": repr(exc)})
            continue
        rdkit_value = int(rdMolDescriptors.CalcNumAtomStereoCenters(rd_mol))
        rdkit_indices = rdkit_atom_stereocenter_indices(rd_mol)
        if len(rdkit_indices) != rdkit_value:
            direct_indices = rdkit_potential_tetrahedral_indices(rd_mol)
            if len(direct_indices) == rdkit_value:
                rdkit_indices = direct_indices
                atom_oracle_fallback_count += 1
            else:
                atom_oracle_unresolved_count += 1
        false_positives = sorted(chematic_indices - rdkit_indices)
        false_negatives = sorted(rdkit_indices - chematic_indices)
        atom_false_positive_count += len(false_positives)
        atom_false_negative_count += len(false_negatives)
        if (false_positives or false_negatives) and len(atom_mismatch_examples) < 50:
            atom_mismatch_examples.append(
                {
                    "index": index,
                    "smiles": smi,
                    "chematic_indices": sorted(chematic_indices),
                    "rdkit_indices": sorted(rdkit_indices),
                    "false_positives": false_positives,
                    "false_negatives": false_negatives,
                }
            )
        compared += 1
        if chematic_value != rdkit_value:
            mismatch_count += 1
            if len(mismatch_examples) < 50:
                mismatch_examples.append(
                    {
                        "index": index,
                        "smiles": smi,
                        "chematic": chematic_value,
                        "rdkit": rdkit_value,
                        "assigned_rdkit": sum(
                            1 for atom in rd_mol.GetAtoms() if atom.HasProp("_CIPCode")
                        ),
                    }
                )

    result = {
        "profile": "A1-stereocenter-potential",
        "version": "1.0.13",
        "rdkit_version": rdBase.rdkitVersion,
        "chematic_version": getattr(chematic, "__version__", None),
        "python_provenance": provenance,
        "corpus": str(corpus.relative_to(ROOT)) if corpus.is_relative_to(ROOT) else str(corpus),
        "corpus_sha256": hashlib.sha256(raw).hexdigest(),
        "input_rows": len(rows),
        "compared": compared,
        "rdkit_parse_failures": rdkit_parse_failures,
        "chematic_failures": chematic_failures,
        "mismatch_count": mismatch_count,
        "mismatch_examples": mismatch_examples,
        "mismatches": mismatch_examples,
        "mismatch_examples_truncated": mismatch_count > len(mismatch_examples),
        "atom_false_positive_count": atom_false_positive_count,
        "atom_false_negative_count": atom_false_negative_count,
        "atom_oracle_fallback_count": atom_oracle_fallback_count,
        "atom_oracle_unresolved_count": atom_oracle_unresolved_count,
        "atom_mismatch_examples": atom_mismatch_examples,
        "atom_mismatch_examples_truncated": (
            atom_false_positive_count + atom_false_negative_count
            > sum(
                len(row["false_positives"]) + len(row["false_negatives"])
                for row in atom_mismatch_examples
            )
        ),
        "error_count": error_count,
        "error_examples": error_examples,
        "error_examples_truncated": error_count > len(error_examples),
        "gate_passed": (
            compared == len(rows)
            and rdkit_parse_failures == 0
            and chematic_failures == 0
            and mismatch_count == 0
            and atom_false_positive_count == 0
            and atom_false_negative_count == 0
            and atom_oracle_unresolved_count == 0
            and error_count == 0
        ),
        "definition": "potential tetrahedral stereocenters; assigned and unspecified centers included",
        "oracle": "rdMolDescriptors.CalcNumAtomStereoCenters",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"stereocenter RDKit parity: {compared}/{len(rows)} compared, "
        f"mismatches={mismatch_count}, gate_passed={result['gate_passed']}"
    )
    return 0 if result["gate_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
