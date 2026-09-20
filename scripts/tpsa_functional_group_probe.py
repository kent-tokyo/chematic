#!/usr/bin/env python3
"""Run a public, non-sealed RDKit TPSA functional-group probe.

The cases are deliberately small, named SMILES covering the N/O/S/P branches
implemented by ``chematic_chem::tpsa``.  This complements corpus-scale
evidence with readable atom-type regressions; it is neither a replacement for
an unused evaluation cohort nor input to one.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path


CASES = {
    "water": "O",
    "alcohol": "CCO",
    "phenol": "Oc1ccccc1",
    "ether": "COC",
    "epoxide": "C1CO1",
    "carbonyl": "CC=O",
    "carboxylic_acid": "CC(=O)O",
    "carboxylate": "CC(=O)[O-]",
    "ester": "CC(=O)OC",
    "amide": "CC(=O)N",
    "primary_amine": "CCN",
    "tertiary_amine": "CN(C)C",
    "quaternary_ammonium": "C[N+](C)(C)C",
    "nitrile": "CC#N",
    "imine": "CC=N",
    "azo": "CN=NC",
    "nitro": "C[N+](=O)[O-]",
    "nitroso": "CN=O",
    "n_oxide": "C[n+]1ccccc1[O-]",
    "azide": "CN=[N+]=[N-]",
    "pyridine": "n1ccccc1",
    "pyrrole": "c1cc[nH]c1",
    "imidazole": "c1ncc[nH]1",
    "sulfide": "CSC",
    "thiol": "CS",
    "sulfoxide": "CS(=O)C",
    "sulfone": "CS(=O)(=O)C",
    "sulfonamide": "CS(=O)(=O)N",
    "sulfonate": "CS(=O)(=O)[O-]",
    "thiophene": "c1ccsc1",
    "phosphate": "COP(=O)(O)O",
    "phosphonate": "CP(=O)(O)O",
    "phosphorothioate": "COP(=S)(O)O",
    "phosphine": "CP(C)C",
    "phosphazene_implicit_h": "CP(=N)C",
    "phosphazene_fully_substituted": "CP(=N)(C)C",
    "phosphazene_explicit_h": "[PH](=N)C",
    "phosphazene_amine": "CP(=N)N",
    "seleninic_acid": "C[Se](=O)O",
    "boronic_acid": "B(O)O",
    "hydrazine": "CNN",
    "guanidine": "NC(=N)N",
    "urea": "NC(=O)N",
    "isocyanate": "CN=C=O",
    "isothiocyanate": "CN=C=S",
    "triazine": "n1cncnc1",
    # Non-sealed representation strata: isotope labels, tautomer spellings,
    # and an anionic aromatic oxygen. These are public, hand-auditable inputs
    # rather than any row or derivative from a sealed candidate.
    "carbon13_ethanol": "[13CH3]CO",
    "carbon11_ethanol": "[11CH3]CO",
    "deuterated_ethanol": "CCO[2H]",
    "2_hydroxypyridine": "Oc1ccccn1",
    "phenoxide": "[O-]c1ccccc1",
}

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CLI = ROOT / "target" / "debug" / "chematic"


def chematic_tpsa(cli: Path, smiles: str) -> float:
    run = subprocess.run(
        [str(cli), "descriptors", smiles], capture_output=True, text=True
    )
    if run.returncode:
        raise RuntimeError(f"chematic CLI rejected '{smiles}': {run.stderr.strip()}")
    return float(json.loads(run.stdout)["tpsa"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, default=DEFAULT_CLI)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cli = args.cli.resolve()
    if not cli.is_file():
        parser.error(f"chematic CLI not found: {cli}")

    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    rows = []
    for case_id, smiles in CASES.items():
        rd_mol = Chem.MolFromSmiles(smiles)
        if rd_mol is None:
            raise RuntimeError(f"RDKit rejected public probe {case_id}: {smiles}")
        expected = rdMolDescriptors.CalcTPSA(rd_mol, includeSandP=True)
        actual = chematic_tpsa(cli, smiles)
        rows.append(
            {
                "id": case_id,
                "smiles": smiles,
                "rdkit_tpsa": expected,
                "chematic_tpsa": actual,
                "absolute_error": abs(actual - expected),
                "strict_match": abs(actual - expected) <= 1e-6,
            }
        )

    failures = [row["id"] for row in rows if not row["strict_match"]]
    with cli.open("rb") as handle:
        cli_sha256 = hashlib.file_digest(handle, "sha256").hexdigest()

    result = {
        "schema_version": 1,
        "profile": "tpsa_public_functional_group_probe_v1",
        "status": "development_regression",
        "scope": "Public, non-sealed functional-group probe; not an unused evaluation cohort or adoption decision.",
        "oracle": {
            "engine": "RDKit",
            "version": Chem.rdBase.rdkitVersion,
            "operation": "rdMolDescriptors.CalcTPSA(includeSandP=True)",
        },
        "implementation": {
            "engine": "chematic",
            "operation": "CLI descriptors <SMILES> JSON tpsa",
            "cli_name": cli.name,
            "cli_sha256": cli_sha256,
        },
        "rows": len(rows),
        "strict_tolerance": 1e-6,
        "strict_matches": len(rows) - len(failures),
        "failures": failures,
        "cases": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in ("rows", "strict_matches", "failures")}))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
