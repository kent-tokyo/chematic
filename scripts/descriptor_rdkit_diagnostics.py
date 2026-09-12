#!/usr/bin/env python3
"""Measure and classify chematic/RDKit descriptor differences.

The script intentionally requires a source-built chematic package exposing
``Mol.rdkit_mw``.  It does not substitute the native ``mw`` value when that
profile is absent, because doing so would make the compatibility result
unverifiable.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


FIELDS = {
    "molecular_weight": ("rdkit_mw", "mw", 0.01),
    "hba": ("rdkit_hba", "hba", 0.0),
    "hbd": ("hbd", "hbd", 0.0),
    "tpsa": ("tpsa", "tpsa", 0.1),
    "logp": ("logp", "logp", 0.01),
    "molar_refractivity": ("molar_refractivity", "mr", 0.01),
    "fsp3": ("fsp3", "fsp3", 0.001),
    "aromatic_ring_count": ("rdkit_aromatic_ring_count", "aromatic_ring_count", 0.0),
}
LEGACY_FIELDS = {
    "molecular_weight": ("rdkit_mw", "mw", 0.01),
    "hba": ("hba", "hba", 0.0),
    "hbd": ("hbd", "hbd", 0.0),
    "tpsa": ("tpsa", "tpsa", 0.1),
    "logp": ("logp", "logp", 0.01),
    "molar_refractivity": ("molar_refractivity", "mr", 0.01),
    "fsp3": ("fsp3", "fsp3", 0.001),
    "aromatic_ring_count": ("aromatic_ring_count", "aromatic_ring_count", 0.0),
}
STRICT_TOLERANCES = {
    "molecular_weight": 1e-6,
    "hba": 0.0,
    "hbd": 0.0,
    "tpsa": 1e-6,
    "logp": 1e-6,
    "molar_refractivity": 1e-6,
    "fsp3": 1e-6,
    "aromatic_ring_count": 0.0,
}


def classify(name: str, rd_mol, delta: float) -> str:
    if name == "molecular_weight":
        if any(atom.GetIsotope() for atom in rd_mol.GetAtoms()):
            return "isotope_table"
        if any(atom.GetAtomicNum() in (5, 16, 34) for atom in rd_mol.GetAtoms()):
            return "atomic_mass_table"
        return "implicit_hydrogen_or_valence"
    if name in {"aromatic_ring_count", "fsp3"}:
        return "ring_perception_or_aromaticity"
    if name in {"logp", "molar_refractivity"}:
        return "crippen_contribution_or_aromaticity"
    if name in {"hba", "hbd", "tpsa"}:
        return "charge_or_valence_rule"
    return "unclassified"


def load_smiles(path: Path) -> list[str]:
    """Load either a headerless .smi file or a CSV with a SMILES column."""
    if path.suffix.lower() in {".smi", ".smiles"}:
        values = []
        for line in path.read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                values.append(line.split()[0])
        return values
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return []
    column = "SMILES" if "SMILES" in rows[0] else next(iter(rows[0]))
    return [row[column].strip() for row in rows if row.get(column, "").strip()]


def python_provenance(chematic) -> dict:
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
        "workspace_path_visible": str(ROOT) in str(package_init),
        "import_origin": (
            "temporary_extracted_wheel"
            if "site-packages" not in str(package_init)
            else "site_package"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("smiles_csv", type=Path)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument(
        "--native-regression",
        type=Path,
        default=Path("validation/results/native-descriptor-regression-v1.0.13.json"),
        help="native-profile regression artifact used to prove the default API was preserved",
    )
    parser.add_argument(
        "--legacy-api",
        action="store_true",
        help="use the v1 baseline attribute map and record it explicitly",
    )
    args = parser.parse_args()

    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
    import rdkit
    import chematic

    smiles_values = load_smiles(args.smiles_csv)
    field_map = LEGACY_FIELDS if args.legacy_api else FIELDS
    try:
        probe = chematic.from_smiles("CS")
        if not hasattr(probe, "rdkit_mw"):
            raise RuntimeError("source-built chematic must expose Mol.rdkit_mw")
    except Exception as exc:
        raise SystemExit(f"compatibility profile unavailable: {exc}") from exc

    stats = {
        name: {"tolerance": tol, "strict_tolerance": STRICT_TOLERANCES[name],
               "parsed": 0, "matches": 0, "strict_matches": 0, "mismatches": 0,
               "mae": 0.0, "median_abs_error": None, "p95_abs_error": None,
               "max_abs_error": 0.0, "errors": [], "_deltas": []}
        for name, (_, _, tol) in FIELDS.items()
    }
    causes: Counter[str] = Counter()
    parsed = failures = unsupported_values = 0
    raw_rows = []
    for index, smiles in enumerate(smiles_values):
        if not smiles:
            continue
        rd_mol = Chem.MolFromSmiles(smiles)
        try:
            ch_mol = chematic.from_smiles(smiles)
        except (TypeError, ValueError, RuntimeError):
            ch_mol = None
        if rd_mol is None or ch_mol is None:
            failures += 1
            raw_rows.append({"index": index, "smiles": smiles, "status": "invalid_input"})
            continue
        parsed += 1
        expected = {
            "molecular_weight": Descriptors.MolWt(rd_mol),
            "hba": Lipinski.NumHAcceptors(rd_mol),
            "hbd": Lipinski.NumHDonors(rd_mol),
            "tpsa": rdMolDescriptors.CalcTPSA(rd_mol, includeSandP=True),
            "logp": Crippen.MolLogP(rd_mol),
            "molar_refractivity": Crippen.MolMR(rd_mol),
            "fsp3": rdMolDescriptors.CalcFractionCSP3(rd_mol),
            "aromatic_ring_count": rdMolDescriptors.CalcNumAromaticRings(rd_mol),
        }
        actual = {name: getattr(ch_mol, attr) for name, (attr, _, _) in field_map.items()}
        raw_fields = {}
        for name, values in expected.items():
            value = actual[name]
            if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(float(value)):
                unsupported_values += 1
                raw_fields[name] = {"status": "unsupported", "reason": "non_finite_or_non_numeric"}
                continue
            delta = abs(float(value) - float(values))
            item = stats[name]
            item["parsed"] += 1
            item["mae"] += delta
            item["_deltas"].append(delta)
            item["max_abs_error"] = max(item["max_abs_error"], delta)
            if delta <= item["tolerance"]:
                item["matches"] += 1
            if delta <= item["strict_tolerance"]:
                item["strict_matches"] += 1
                strict_passed = True
            else:
                item["mismatches"] += 1
                strict_passed = False
                cause = classify(name, rd_mol, delta)
                causes[cause] += 1
                if len(item["errors"]) < 20:
                    item["errors"].append({"smiles": smiles, "delta": delta, "cause": cause,
                                            "chematic": value, "rdkit": values})
            raw_fields[name] = {
                "status": "ok",
                "chematic": value,
                "rdkit": values,
                "absolute_error": delta,
                "matched": delta <= item["tolerance"],
                "strict_passed": strict_passed,
            }
        raw_rows.append({"index": index, "smiles": smiles, "status": "ok", "fields": raw_fields})
    for item in stats.values():
        deltas = sorted(item.pop("_deltas"))
        item["mae"] = item["mae"] / item["parsed"] if item["parsed"] else None
        if deltas:
            item["median_abs_error"] = deltas[(len(deltas) - 1) // 2]
            item["p95_abs_error"] = deltas[min(len(deltas) - 1, math.ceil(len(deltas) * 0.95) - 1)]
    native_regression = json.loads(args.native_regression.read_text())
    raw_payload = json.dumps(raw_rows, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
    result = {"schema_version": 2, "corpus": str(args.smiles_csv),
              "corpus_sha256": hashlib.sha256(args.smiles_csv.read_bytes()).hexdigest(),
              "rdkit_version": rdkit.__version__,
              "chematic_version": getattr(chematic, "__version__", "unknown"),
              "python_provenance": python_provenance(chematic),
              "rows": len(smiles_values), "parsed": parsed, "parse_failures": failures,
              "unsupported_values": unsupported_values,
              "native_default_preserved": native_regression.get("gate_passed") is True
              and native_regression.get("failed") == 0,
              "native_regression_artifact": str(args.native_regression),
              "profile": "rdkit_compat_v1_legacy" if args.legacy_api else "rdkit_compat_v2",
              "contract": "rdkit_descriptor_semantics_v1",
              "api_profile": "legacy_native_attributes" if args.legacy_api else "rdkit_compat_attributes",
              "field_attribute_map": {name: attrs[0] for name, attrs in field_map.items()},
              "raw_rows": raw_rows,
              "raw_rows_sha256": hashlib.sha256(raw_payload).hexdigest(),
              "fields": stats,
              "mismatch_causes": dict(sorted(causes.items()))}
    args.json.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
