#!/usr/bin/env python3
"""Record bounded writer-stability checks and optional RDKit reproduction."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import chematic as ch
import chematic.chematic as native
import chematic._rdkit_smiles_checked as implementation

ISSUE_CASES = [
    "C[N@+]1(F)CCC2(CC[C@H](Cl)CC2)CC1",
    "C[N@+]1(F)CC[C@]2(CC[C@H](Cl)CC2)CC1",
    "C[C@]1(F)CCC2(CC[C@H](Cl)CC2)CC1",
    "C[C@]1(F)CCC2(CC[C@](I)(Cl)CC2)CC1",
]


def inspect(source):
    row = {"source": source}
    try:
        mol = ch.from_smiles(source)
        before = mol.to_moljson()
        try:
            result = ch.rdkit_smiles_checked(mol)
        except ch.RdkitSmilesStabilityError as error:
            row.update(
                status=error.code,
                first=error.first_smiles,
                roundtrip=error.roundtrip_smiles,
                atom_order=error.atom_order,
            )
        except (ValueError, RuntimeError) as error:
            row.update(status="writer_refused", error=str(error))
        else:
            assert ch.from_smiles(result).rdkit_smiles == result
            row.update(status="stable", output=result)
        assert mol.to_moljson() == before, "input mutation"
    except (ValueError, RuntimeError) as error:
        row.update(status="input_refused", error=str(error))
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--reference", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    corpus = root / "crates/chematic-smiles/tests/data/rdkit_canonical_smiles.tsv"
    sources = [
        line.split("\t")[0]
        for line in corpus.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    rows = [inspect(source) for source in sources]
    evidence = {
        "scope": (
            "Python fixed-point plus 32-order bounded atom-permutation guard; "
            "not a stereo identity proof or exhaustive canonicality proof"
        ),
        "native_version": ch.__version__,
        "native_sha256": hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
        "guard_sha256": hashlib.sha256(
            Path(implementation.__file__).read_bytes()
        ).hexdigest(),
        "corpus_sha256": hashlib.sha256(corpus.read_bytes()).hexdigest(),
        "counts": dict(Counter(row["status"] for row in rows)),
        "issue_cases": [inspect(source) for source in ISSUE_CASES],
        "corpus": rows,
    }
    if args.reference:
        from rdkit import Chem, rdBase

        previous = Chem.GetUseLegacyStereoPerception()
        reference = []
        try:
            Chem.SetUseLegacyStereoPerception(True)
            for source in ISSUE_CASES:
                first = Chem.MolToSmiles(Chem.MolFromSmiles(source))
                second = Chem.MolToSmiles(Chem.MolFromSmiles(first))
                reference.append(
                    {
                        "source": source,
                        "first": first,
                        "roundtrip": second,
                        "stable": first == second,
                    }
                )
        finally:
            Chem.SetUseLegacyStereoPerception(previous)
        evidence["reference"] = {
            "version": rdBase.rdkitVersion,
            "legacy_stereo": True,
            "cases": reference,
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(evidence, indent=2) + "\n")
    print(
        json.dumps(
            {
                "counts": evidence["counts"],
                "issue_statuses": [row["status"] for row in evidence["issue_cases"]],
            }
        )
    )


if __name__ == "__main__":
    main()
