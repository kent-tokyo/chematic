#!/usr/bin/env python3
"""Classify public functional-group descriptor differences against RDKit.

This development probe deliberately uses the readable, non-sealed cases from
``tpsa_functional_group_probe.py``.  It records every descriptor independently
instead of converting an aggregate percentage into an adoption claim.  In
particular, CheMatic's CLI exposes its native molecular-weight table; an
observed mass difference is classified as a contract difference until the
RDKit-compatible mass API is invoked through an equivalent binding.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

from tpsa_functional_group_probe import CASES

# This Kekule tautomer is deliberately a classification-only residual.  The
# pass/fail TPSA atom-type probe above contains only established regressions;
# keeping a known mismatch here preserves its public evidence without turning
# the narrow green gate into an aggregate accuracy claim.
CLASSIFICATION_ONLY_CASES = {
    "2_pyridone": "O=C1C=CC=CN1",
}


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CLI = ROOT / "target" / "debug" / "chematic"

CASE_STRATA = {
    "carbon13_ethanol": ["isotope"],
    "carbon11_ethanol": ["isotope"],
    "deuterated_ethanol": ["isotope"],
    "2_hydroxypyridine": ["tautomer", "aromatic"],
    "2_pyridone": ["tautomer", "kekule"],
    "phenoxide": ["charge", "aromatic"],
}

FIELDS = {
    "molecular_weight": {
        "rdkit": lambda mol, modules: modules["Descriptors"].MolWt(mol),
        "tolerance": 0.01,
        "comparison": "native_mass_vs_rdkit_molwt_contract_difference",
    },
    "rdkit_molecular_weight": {
        "rdkit": lambda mol, modules: modules["Descriptors"].MolWt(mol),
        "tolerance": 0.01,
        "comparison": "numeric_rdkit_mass_profile",
    },
    "exact_mass": {
        "rdkit": lambda mol, modules: modules["Descriptors"].ExactMolWt(mol),
        "tolerance": 1e-6,
        "comparison": "numeric_exact_mass",
    },
    "heavy_atoms": {
        "rdkit": lambda mol, modules: mol.GetNumHeavyAtoms(),
        "tolerance": 0.0,
        "comparison": "integer",
    },
    "logp": {
        "rdkit": lambda mol, modules: modules["Crippen"].MolLogP(mol),
        "tolerance": 1e-2,
        "comparison": "numeric",
    },
    "molar_refractivity": {
        "rdkit": lambda mol, modules: modules["Crippen"].MolMR(mol),
        "tolerance": 1e-2,
        "comparison": "numeric",
    },
    "tpsa": {
        "rdkit": lambda mol, modules: modules["rdMolDescriptors"].CalcTPSA(
            mol, includeSandP=True
        ),
        "tolerance": 1e-6,
        "comparison": "numeric",
    },
    "hbd": {
        "rdkit": lambda mol, modules: modules["rdMolDescriptors"].CalcNumHBD(mol),
        "tolerance": 0.0,
        "comparison": "integer",
    },
    "hba": {
        "rdkit": lambda mol, modules: modules["rdMolDescriptors"].CalcNumHBA(mol),
        "tolerance": 0.0,
        "comparison": "integer",
    },
    "fsp3": {
        "rdkit": lambda mol, modules: modules["rdMolDescriptors"].CalcFractionCSP3(mol),
        "tolerance": 1e-12,
        "comparison": "numeric",
    },
    "aromatic_ring_count": {
        "rdkit": lambda mol, modules: modules["rdMolDescriptors"].CalcNumAromaticRings(mol),
        "tolerance": 0.0,
        "comparison": "integer",
    },
    "rotatable_bonds": {
        "rdkit": lambda mol, modules: modules["rdMolDescriptors"].CalcNumRotatableBonds(mol),
        "tolerance": 0.0,
        "comparison": "integer_rdkit_default",
    },
}


def chematic_descriptors(cli: Path, smiles: str) -> dict[str, float | int]:
    run = subprocess.run(
        [str(cli), "descriptors", smiles], capture_output=True, text=True
    )
    if run.returncode:
        raise RuntimeError(f"CheMatic CLI rejected '{smiles}': {run.stderr.strip()}")
    return json.loads(run.stdout)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, default=DEFAULT_CLI)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cli = args.cli.resolve()
    if not cli.is_file():
        parser.error(f"CheMatic CLI not found: {cli}")

    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

    modules = {
        "Crippen": Crippen,
        "Descriptors": Descriptors,
        "rdMolDescriptors": rdMolDescriptors,
    }
    rows = []
    summary = {field: {"strict_matches": 0, "mismatches": []} for field in FIELDS}
    for case_id, smiles in {**CASES, **CLASSIFICATION_ONLY_CASES}.items():
        rd_mol = Chem.MolFromSmiles(smiles)
        if rd_mol is None:
            raise RuntimeError(f"RDKit rejected public probe {case_id}: {smiles}")
        actual = chematic_descriptors(cli, smiles)
        values = {}
        for field, spec in FIELDS.items():
            expected = spec["rdkit"](rd_mol, modules)
            observed = actual[field]
            error = abs(float(observed) - float(expected))
            strict = error <= spec["tolerance"]
            values[field] = {
                "rdkit": expected,
                "chematic": observed,
                "absolute_error": error,
                "strict_match": strict,
            }
            if strict:
                summary[field]["strict_matches"] += 1
            else:
                summary[field]["mismatches"].append(case_id)
        rows.append(
            {
                "id": case_id,
                "smiles": smiles,
                "strata": CASE_STRATA.get(case_id, ["functional_group"]),
                "descriptors": values,
            }
        )

    with cli.open("rb") as handle:
        cli_sha256 = hashlib.file_digest(handle, "sha256").hexdigest()
    result = {
        "schema_version": 1,
        "profile": "descriptor_public_functional_group_classification_v1",
        "status": "development_classification",
        "scope": (
            "Public, non-sealed functional-group classification; not an unused "
            "evaluation cohort or adoption decision."
        ),
        "oracle": {
            "engine": "RDKit",
            "version": Chem.rdBase.rdkitVersion,
            "operations": {
                field: {key: value for key, value in spec.items() if key != "rdkit"}
                for field, spec in FIELDS.items()
            },
        },
        "implementation": {
            "engine": "CheMatic",
            "operation": "CLI descriptors <SMILES> JSON fields",
            "cli_sha256": cli_sha256,
        },
        "rows": len(rows),
        "summary": summary,
        "cases": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"rows": result["rows"], "summary": summary}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
