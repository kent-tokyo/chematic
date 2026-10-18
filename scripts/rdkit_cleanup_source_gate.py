#!/usr/bin/env python3
"""Compare source-built Cleanup JSONL with RDKit on the stress corpus.

Generate the JSONL with ``cargo run -p chematic-chem --example
cleanup_rdkit_dump --offline < validation/cosmolkit_comparison/corpora/stress-v1.smi``.
Only rows whose parsed input SMILES agree are scored.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "validation/cosmolkit_comparison/corpora/stress-v1.smi"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def value(result: dict) -> tuple[str | None, str | None]:
    if "Ok" in result:
        return result["Ok"], None
    return None, result.get("Err", "missing result")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dump", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--previous-source", type=Path)
    args = parser.parse_args()
    if args.baseline and args.previous_source:
        parser.error("choose one comparison report")
    if rdBase.rdkitVersion != "2026.09.1":
        parser.error(f"expected RDKit 2026.09.1; found {rdBase.rdkitVersion}")
    RDLogger.DisableLog("rdApp.*")

    corpus = [
        line.split()[0] if line.split() else None
        for line in CORPUS.read_text(encoding="utf-8").splitlines()
    ]
    records = [json.loads(line) for line in args.dump.read_text().splitlines()]
    if len(records) != sum(value is not None for value in corpus):
        raise ValueError("source dump does not cover every nonblank corpus row")
    rows = {}
    for record in records:
        index = record["index"]
        if index in rows or index >= len(corpus) or record["input"] != corpus[index]:
            raise ValueError(f"source dump index/input mismatch at {index}")
        rows[index] = record

    compatible, exact, changed = set(), set(), set()
    excluded, errors, mismatches = [], [], []
    for index, source in sorted(rows.items()):
        input_text = source["input"]
        rdkit_mol = Chem.MolFromSmiles(input_text)
        rdkit_input = (
            Chem.MolToSmiles(rdkit_mol, canonical=True, isomericSmiles=True)
            if rdkit_mol is not None else None
        )
        source_input, parse_error = value(source.get("rdkit_input", {}))
        if parse_error is not None or source_input != rdkit_input:
            excluded.append({"index": index, "input": input_text,
                             "source_input": source_input, "rdkit_input": rdkit_input})
            continue
        try:
            expected = Chem.MolToSmiles(
                rdMolStandardize.Cleanup(rdkit_mol), canonical=True, isomericSmiles=True
            )
        except Exception as exc:
            errors.append({"index": index, "input": input_text,
                           "oracle_error": f"{type(exc).__name__}: {exc}"})
            continue
        compatible.add(index)
        actual, source_error = value(source.get("cleanup", {}))
        if actual == expected:
            exact.add(index)
        else:
            mismatches.append({"index": index, "input": input_text,
                               "rdkit_cleanup": expected, "chematic_cleanup": actual,
                               "source_error": source_error})
        if expected != rdkit_input:
            changed.add(index)

    comparison = None
    if args.baseline:
        baseline = json.loads(args.baseline.read_text())
        if baseline["sampled_rows"] != len(corpus):
            raise ValueError("baseline and source sample counts differ")
        old_scored = set(range(len(corpus))) - {
            row["sample_index"] for field in ("excluded", "oracle_errors")
            for row in baseline[field]
        }
        if old_scored != compatible:
            raise ValueError("baseline and source input compatibility sets differ")
        old_mismatches = {row["sample_index"] for row in baseline["mismatches"]}
        comparison = {
            "baseline_exact": baseline["cleanup_exact"],
            "fixed_indices": sorted(old_mismatches & exact),
            "regressed_indices": sorted((old_scored - old_mismatches) - exact),
        }
    if args.previous_source:
        previous = json.loads(args.previous_source.read_text())
        if (previous["source_rows"] != len(corpus)
                or previous["compatible_rows"] != len(compatible)
                or previous["corpus_sha256"] != sha256(CORPUS)
                or previous["rdkit_version"] != rdBase.rdkitVersion
                or {row["index"] for row in previous["oracle_errors"]}
                    != {row["index"] for row in errors}
                or {row["index"] for row in previous["excluded"]}
                    != {row["index"] for row in excluded}):
            raise ValueError("previous source report uses a different corpus or oracle")
        previous_mismatches = {row["index"] for row in previous["mismatches"]}
        current_mismatches = {row["index"] for row in mismatches}
        comparison = {
            "previous_source": str(args.previous_source),
            "previous_exact": previous["exact_rows"],
            "fixed_indices": sorted(previous_mismatches - current_mismatches),
            "regressed_indices": sorted(current_mismatches - previous_mismatches),
        }

    report = {
        "profile": "source_cleanup_vs_rdkit_2026_09_1",
        "source_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "source_sha256": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (
                ROOT / "crates/chematic-chem/src/standardize.rs",
                ROOT / "crates/chematic-chem/src/standardize/cleanup.rs",
                ROOT / "crates/chematic-chem/src/standardize/groups.rs",
                ROOT / "crates/chematic-core/src/molecule/builder.rs",
                ROOT / "crates/chematic-chem/examples/cleanup_rdkit_dump.rs",
                ROOT / "crates/chematic-smiles/src/lib.rs",
                ROOT / "crates/chematic-smiles/src/rdkit/mod.rs",
            )
        },
        "corpus": str(CORPUS.relative_to(ROOT)),
        "corpus_sha256": sha256(CORPUS),
        "source_dump_sha256": sha256(args.dump),
        "rdkit_version": rdBase.rdkitVersion,
        "source_rows": len(corpus),
        "compatible_rows": len(compatible),
        "exact_rows": len(exact),
        "rdkit_changed_rows": len(changed),
        "excluded": excluded,
        "oracle_errors": errors,
        "mismatches": mismatches,
        "baseline_comparison": comparison,
        "limits": ["same parsed inputs only", "RDKit agreement is not independent chemical truth"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(f"Cleanup: {len(exact)}/{len(compatible)} exact")
    if comparison:
        print(f"Fixed: {len(comparison['fixed_indices'])}; "
              f"regressed: {len(comparison['regressed_indices'])}")


if __name__ == "__main__":
    main()
