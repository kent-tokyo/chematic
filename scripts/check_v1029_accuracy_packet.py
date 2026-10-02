#!/usr/bin/env python3
"""Check v1.0.29 published-wheel and later source-candidate evidence without RDKit.

The archived rows are exposed data. This validates the published artifact's
outcomes and their correspondence to the separately classified #634/#635
residuals; it does not bless a later source candidate or a sealed cohort.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

from check_rdkit_rebaseline_evidence import validate_smarts_cell_correspondence

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "validation/results"
PREFIX = "2026.03.6-2026-10-02"
ROW_PATH = RESULTS / f"rdkit-python-chemistry-rows-v1.0.29-vs-{PREFIX}.jsonl.gz"
ROW_SHA256 = "e9b8592b69a0663608ab51d0d83d3db644a515a7cbfedefe5ad5aff61780b968"
GZIP_SHA256 = "95655e9159ee4d5ea59049723b5e2c352616826ade2ccb1500ceb47ab5f91012"
WHEEL_SHA256 = "a9ceb3c685e6cf4abbbdb1ae1cdbaa76b9ba13c8015d53750ba617cba3f7c1d8"
SOURCE_DEV_WHEEL_SHA256 = "b03642d52f4312c9f55221a70bdb7403cd6373fc4b0095b87ea2df7a04dc83e1"


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def check(condition: bool, message: str, errors: list[str]) -> None:
    if not condition:
        errors.append(message)


def main() -> int:
    try:
        current = read_json(RESULTS / f"rdkit-python-chemistry-v1.0.29-vs-{PREFIX}.json")
        previous = read_json(RESULTS / f"rdkit-python-chemistry-v1.0.28-wheel-vs-{PREFIX}.json")
        reaction = read_json(RESULTS / f"reaction-product-parity-v1.0.29-vs-{PREFIX}.json")
        matrix = read_json(ROOT / "benchmarks/2026-10-02-v1.0.29-python-op-matrix.json")
        source_hba = read_json(RESULTS / "rdkit-hba-source-dev-wheel-5k-2026-10-02.json")
        classified = read_json(
            RESULTS / "rdkit-rebaseline-residual-classification-v1.0.27-issue634-vs-2026.03.6-2026-09-28.json"
        )
        exposed = read_json(
            RESULTS / "rdkit-rebaseline-python-chemistry-summary-v1.0.27-issue634-vs-2026.03.6-2026-09-28.json"
        )
        compressed = ROW_PATH.read_bytes()
        raw = gzip.decompress(compressed)
        rows = [json.loads(line) for line in raw.splitlines()]
    except (OSError, ValueError, gzip.BadGzipFile, json.JSONDecodeError) as exc:
        print(f"v1.0.29 accuracy packet invalid: {exc}", file=sys.stderr)
        return 1

    errors: list[str] = []
    check(hashlib.sha256(compressed).hexdigest() == GZIP_SHA256, "compressed rows hash", errors)
    check(hashlib.sha256(raw).hexdigest() == ROW_SHA256, "raw rows hash", errors)
    check(current.get("rows", {}).get("archived_path") == str(ROW_PATH.relative_to(ROOT))
          and current.get("rows", {}).get("archived_sha256") == GZIP_SHA256,
          "archived row reference", errors)
    check(len(rows) == 10000, "row count", errors)
    check(all(row.get("input_index") == i for i, row in enumerate(rows)), "row indices", errors)
    check(all(row.get("status") == "completed" for row in rows), "terminal row status", errors)

    for version, summary in (("1.0.28", previous), ("1.0.29", current)):
        check(summary.get("chematic_version") == version, f"{version}: version", errors)
        check(summary.get("rdkit_version") == "2026.03.6", f"{version}: oracle", errors)
        check(summary.get("gate_passed") is True, f"{version}: runner gate", errors)
        check(summary.get("sealed_accuracy_cohort_reused") is False, f"{version}: sealed status", errors)
        check(
            summary.get("row_accounting")
            == {"input_count": 10000, "completed_count": 10000, "parse_failure_count": 0},
            f"{version}: row accounting", errors,
        )
        check(summary.get("rows", {}).get("sha256") == ROW_SHA256, f"{version}: row digest", errors)
        check(summary.get("counts") == current.get("counts"), f"{version}: operation counts", errors)
    check(exposed.get("rows", {}).get("sha256") == ROW_SHA256,
          "classified historical source and published rows diverge", errors)

    counts = Counter()
    abstentions = Counter()
    for row in rows:
        smiles = row.get("smiles_parse_write", {})
        check(smiles.get("semantic_roundtrip") is True,
              f"row {row['input_index']}: SMILES semantic change", errors)
        cip = row.get("cip", {})
        counts["cip_exact" if cip.get("exact") else "cip_difference"] += 1
        if not cip.get("exact"):
            candidate = cip.get("chematic_atoms", {})
            oracle = cip.get("rdkit_atoms", {})
            missing = set(oracle) - set(candidate)
            unresolved = cip.get("chematic_unresolved", {})
            check(len(missing) == 1 and set(unresolved) == missing,
                  f"row {row['input_index']}: unaccounted CIP label", errors)
            check(all(oracle.get(atom) == label for atom, label in candidate.items())
                  and cip.get("chematic_bonds") == cip.get("rdkit_bonds")
                  and all(cip.get("index_correspondence", {}).values()),
                  f"row {row['input_index']}: wrong confident CIP label", errors)
            abstentions.update(unresolved.values())
        counts["morgan_exact" if row.get("morgan", {}).get("exact") else "morgan_difference"] += 1
        if not row.get("morgan", {}).get("exact"):
            check("unsupported" in str(row.get("morgan", {}).get("chematic_error", "")).lower(),
                  f"row {row['input_index']}: untyped Morgan failure", errors)
        smarts = row.get("smarts", {})
        counts["smarts_cells"] += smarts.get("query_count", 0)
        counts["smarts_differences"] += smarts.get("difference_count", 0)

    expected = {"cip_exact": 9995, "cip_difference": 5, "morgan_exact": 9999,
                "morgan_difference": 1, "smarts_cells": 310000, "smarts_differences": 200}
    check(dict(counts) == expected, "operation totals", errors)
    check(dict(abstentions) == {"oracle_unstable": 4, "lone_pair_center": 1},
          "CIP abstention reasons", errors)
    validate_smarts_cell_correspondence(rows, classified.get("smarts", {}), errors)

    check(reaction.get("rdkit_version") == "2026.03.6", "reaction oracle", errors)
    check(reaction.get("chematic", {}).get("version") == "1.0.29"
          and reaction.get("chematic", {}).get("artifact_sha256") == WHEEL_SHA256,
          "reaction published artifact", errors)
    check(reaction.get("case_file_sha256") == hashlib.sha256(
        (ROOT / "validation/reaction_product_parity_cases.json").read_bytes()).hexdigest(),
        "reaction fixture hash", errors)
    check(reaction.get("agreement") == {"matched": 57, "total": 57},
          "reaction agreement", errors)
    reaction_cases = reaction.get("cases", [])
    check(len(reaction_cases) == 57
          and len({case.get("id") for case in reaction_cases}) == 57,
          "reaction row correspondence", errors)
    check(all(case.get("agree") is True and not case.get("missing_in_chematic")
              and not case.get("extra_in_chematic") for case in reaction_cases),
          "reaction product-set accounting", errors)

    hba = next((op for op in matrix.get("operations", []) if op.get("op") == "hba"), {})
    check(matrix.get("chematic", {}).get("version") == "1.0.29"
          and matrix.get("chematic", {}).get("artifact_sha256") == WHEEL_SHA256
          and matrix.get("rdkit", {}).get("version") == "2026.03.6",
          "matrix artifact identity", errors)
    check(matrix.get("corpus", {}).get("common_valid_rows") == 5000
          and hba.get("output_agreement") == {"compared": 5000, "agree": 3641},
          "HBA baseline", errors)
    check(source_hba.get("build_profile") == "dev"
          and source_hba.get("wheel_sha256") == SOURCE_DEV_WHEEL_SHA256
          and source_hba.get("rdkit_version") == "2026.03.6"
          and source_hba.get("corpus", {}).get("sha256")
          == "1c47371dcbe37f4e0a141bf545b72bf238de2761fa3894fa251a552d84728d3e"
          and source_hba.get("counts") == {
              "input": 5000, "compared": 5000, "exact_named": 5000,
              "exact_native": 5000, "parse_failed": 0}
          and source_hba.get("differences") == []
          and source_hba.get("gate_passed") is True,
          "source dev-wheel HBA result", errors)
    if errors:
        print("v1.0.29 accuracy packet invalid:\n" + "\n".join(f"- {e}" for e in errors[:30])
              + (f"\n... {len(errors) - 30} more" if len(errors) > 30 else ""),
              file=sys.stderr)
        return 1
    print("v1.0.29 accuracy packet OK: 10,000 rows; CIP 9,995 + 5 typed; "
          "Morgan 9,999 + 1 typed; SMARTS 200/310,000 classified; "
          "reactions 57/57; HBA published 3,641/5,000; source dev-wheel 5,000/5,000")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
