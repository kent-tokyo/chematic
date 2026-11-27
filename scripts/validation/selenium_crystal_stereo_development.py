#!/usr/bin/env python3
"""Measure declared stereo on mapped selenium-containing development crystals."""

import argparse
import hashlib
import importlib.util
import json
import platform
import sys
from pathlib import Path

import chematic
from rdkit import Chem, rdBase

SDF_HASH = "e39322d1eba60bc041fa380f043146e95247d4209d000374591582e6e2acd117"
SPLIT_HASH = "e13f3e18f11eeb75ccd9dd3e96a28426e47ef1d8bc56fa28c53b02ae9bfc8da8"
INPUT_HASH = "e3c2cf4aa5c7c5740fd3d23c670e26535515b97f16a02ee2bf387661513726a1"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def indexed_graph_matches(native, reference):
    atoms = [a[1] for a in native.atom_table] == [a.GetAtomicNum() for a in reference.GetAtoms()]
    native_bonds = {tuple(sorted((b[0], b[1]))) for b in native.bond_table}
    ref_bonds = {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in reference.GetBonds()}
    return atoms and native_bonds == ref_bonds


def assign_3d(query, coordinates):
    molecule = Chem.Mol(query)
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    for i, xyz in enumerate(coordinates):
        conformer.SetAtomPosition(i, xyz)
    molecule.AddConformer(conformer, assignId=True)
    Chem.AssignStereochemistryFrom3D(molecule, replaceExistingTags=True)
    return molecule


def target_selenium_centers(query):
    centers = []
    for atom in query.GetAtoms():
        if (atom.GetAtomicNum() != 34 or atom.GetFormalCharge() != 0 or atom.GetDegree() != 3
                or atom.GetTotalNumHs() != 0 or atom.GetChiralTag() not in
                (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW)):
            continue
        double = [b for b in atom.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE]
        if (len(double) == 1 and double[0].GetOtherAtom(atom).GetAtomicNum() in (7, 8, 16, 34)
                and all(b.GetBondType() in (Chem.BondType.SINGLE, Chem.BondType.DOUBLE)
                        for b in atom.GetBonds())):
            centers.append(atom.GetIdx())
    return centers


def process(row, supplier):
    query = Chem.MolFromSmiles(row["smiles"])
    record = {"index": row["index"], "name": row["name"], "smiles": row["smiles"]}
    try:
        native = chematic.from_smiles(row["smiles"])
        if not indexed_graph_matches(native, query):
            raise ValueError("SMILES parser indexed graph mismatch")
        raw = supplier[row["index"]]
        if raw is None:
            raise ValueError("source SDF parse failed")
        Chem.SanitizeMol(raw)
        ref = Chem.RemoveHs(raw)
        if query.GetNumAtoms() != ref.GetNumAtoms() or query.GetNumBonds() != ref.GetNumBonds():
            raise ValueError("source and SMILES graph sizes differ")
        matches = ref.GetSubstructMatches(query, uniquify=False, useChirality=True, maxMatches=2)
        if not matches:
            raise ValueError("source and SMILES chirality-aware mapping failed")
        mapping = matches[0]
        coordinates = [list(ref.GetConformer().GetAtomPosition(i)) for i in mapping]
        expected = {a.GetIdx(): a.GetProp("_CIPCode") for a in query.GetAtoms() if a.HasProp("_CIPCode")}
        record["mapping"] = list(mapping)
        record["mapping_count_capped_at_two"] = len(matches)
        record["neutral_selenium_targets"] = target_selenium_centers(query)
        record["coordinates"] = coordinates
        record["checks"] = []
        for mirror in (False, True):
            xyz = [[-x, y, z] if mirror else [x, y, z] for x, y, z in coordinates]
            rdkit = assign_3d(query, xyz)
            native_report = native.stereo_3d_diagnostics(xyz)
            rdkit_atoms = []
            for atom, expected_cip in expected.items():
                a = rdkit.GetAtomWithIdx(atom)
                observed = a.GetProp("_CIPCode") if a.HasProp("_CIPCode") else None
                rdkit_atoms.append({"atom": atom, "expected": expected_cip, "observed": observed,
                                    "matches": None if observed is None else observed == expected_cip})
            record["checks"].append({"mirror": mirror, "native": native_report, "rdkit_cip": rdkit_atoms})
        record["status"] = "ok"
    except (ValueError, RuntimeError, IndexError) as error:
        record.update({"status": "error", "error": f"{type(error).__name__}: {error}"})
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--development", required=True)
    parser.add_argument("--sdf", required=True)
    parser.add_argument("--split", required=True)
    parser.add_argument("--wheel", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    if (digest(args.sdf), digest(args.split), digest(args.development)) != (SDF_HASH, SPLIT_HASH, INPUT_HASH):
        raise ValueError("input checksum mismatch")
    rdBase.DisableLog("rdApp.*")
    development = json.loads(Path(args.development).read_text(encoding="utf-8"))["results"]
    split = json.loads(Path(args.split).read_text(encoding="utf-8"))["rows"]
    split_by_index = {r["index"]: r for r in split}
    selected, parse_errors = [], []
    for row in development:
        membership = split_by_index.get(row["index"])
        if membership is None or membership["split"] != "development" or membership["name"] != row["name"]:
            raise ValueError(f"not a verified development record: {row['index']}")
        molecule = Chem.MolFromSmiles(row["smiles"])
        if molecule is None:
            parse_errors.append(row["index"])
        elif any(atom.GetAtomicNum() == 34 for atom in molecule.GetAtoms()):
            selected.append(row)
    supplier = Chem.SDMolSupplier(args.sdf, sanitize=False, removeHs=False)
    results = [process(row, supplier) for row in selected]
    report = {"selection": "all selenium-containing parsed records in the exposed COD development partition",
              "development_count": len(development), "parse_failures": parse_errors,
              "selected_count": len(selected), "rdkit_version": rdBase.rdkitVersion,
              "platform": platform.platform(), "python": sys.version,
              "development_sha256": digest(args.development), "sdf_sha256": digest(args.sdf),
              "split_sha256": digest(args.split), "wheel_sha256": digest(args.wheel),
              "extension_sha256": digest(importlib.util.find_spec("chematic.chematic").origin),
              "script_sha256": digest(__file__), "results": results}
    Path(args.output).write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
