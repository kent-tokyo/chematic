#!/usr/bin/env python3
"""Compare generated conformers with Platinum development structures."""

import argparse
import hashlib
import importlib.metadata
import json
import platform
import time
from pathlib import Path

import chematic
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdMolAlign
from platinum_geometry import check_geometry, prepare_geometry_checks


def graph_matches(other, mol):
    atoms = [a[1] for a in other.atom_table]
    reference_atoms = [a.GetAtomicNum() for a in mol.GetAtoms()]
    bonds = sorted((min(a, b), max(a, b), order) for a, b, order, _ in other.bond_table)
    reference_bonds = sorted(
        (
            min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
            max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
            str(b.GetBondType()).upper(),
        )
        for b in mol.GetBonds()
    )
    return atoms == reference_atoms and bonds == reference_bonds


def rms_values(mol, reference, conformers):
    result = []
    for coords in conformers:
        probe = Chem.Mol(mol)
        probe.RemoveAllConformers()
        conformer = Chem.Conformer(probe.GetNumAtoms())
        for i, point in enumerate(coords):
            conformer.SetAtomPosition(i, point)
        probe.AddConformer(conformer)
        heavy = Chem.RemoveHs(probe)
        result.append(rdMolAlign.GetBestRMS(heavy, reference, maxMatches=100000))
    return result


def stereo_violations(mol, conformers):
    expected = Chem.Mol(mol)
    Chem.AssignStereochemistry(expected, cleanIt=True, force=True)
    atoms = {
        a.GetIdx(): a.GetProp("_CIPCode")
        for a in expected.GetAtoms()
        if a.HasProp("_CIPCode")
    }
    bonds = {
        b.GetIdx(): str(b.GetStereo())
        for b in expected.GetBonds()
        if b.GetStereo() in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ)
    }
    violations = []
    for coords in conformers:
        probe = Chem.Mol(mol)
        probe.RemoveAllConformers()
        conf = Chem.Conformer(probe.GetNumAtoms())
        for index, point in enumerate(coords):
            conf.SetAtomPosition(index, point)
        probe.AddConformer(conf)
        Chem.AssignStereochemistryFrom3D(probe)
        Chem.AssignStereochemistry(probe, cleanIt=True, force=True)
        bad_atoms = [
            index
            for index, code in atoms.items()
            if not probe.GetAtomWithIdx(index).HasProp("_CIPCode")
            or probe.GetAtomWithIdx(index).GetProp("_CIPCode") != code
        ]
        bad_bonds = [
            index
            for index, code in bonds.items()
            if str(probe.GetBondWithIdx(index).GetStereo()) != code
        ]
        violations.append({"atoms": bad_atoms, "bonds": bad_bonds})
    return violations


def score(reference, count, seed, profile, geometry, repair_iters, independent_mmff=False):
    heavy_ref = Chem.RemoveHs(reference)
    smiles = Chem.MolToSmiles(heavy_ref, isomericSmiles=True)
    graph = Chem.MolFromSmiles(smiles)
    if graph is None or len(Chem.GetMolFrags(graph)) != 1:
        return {"status": "invalid_or_multifragment"}
    mol = Chem.AddHs(graph)
    if not independent_mmff and not AllChem.MMFFHasAllMoleculeParams(mol):
        return {"status": "rdkit_mmff_unsupported", "smiles": smiles}
    if not independent_mmff:
        other = chematic.from_smiles(smiles).add_hydrogens()
        if not graph_matches(other, mol):
            return {"status": "atom_graph_mismatch", "smiles": smiles}
    result = {"status": "ok", "smiles": smiles}
    checks = prepare_geometry_checks(mol) if geometry else None
    if checks is not None:
        result["geometry_reference"] = {
            "bonds_supported": len(checks[0]),
            "bonds_total": mol.GetNumBonds(),
            "angles_supported": len(checks[1]),
            "angles_total": sum(
                a.GetDegree() * (a.GetDegree() - 1) // 2 for a in mol.GetAtoms()
            ),
        }
    start = time.perf_counter()
    try:
        if independent_mmff:
            other = chematic.from_smiles(smiles).add_hydrogens()
            if not graph_matches(other, mol):
                raise ValueError("indexed atom or bond graph mismatch")
        if profile == "compatible":
            raw = other.rdkit_embed_multiple_confs(count, random_seed=seed)
            result["chematic_raw_seconds"] = time.perf_counter() - start
            result["chematic_raw_rmsd"] = rms_values(mol, heavy_ref, raw)
            result["chematic_raw_stereo"] = stereo_violations(mol, raw)
            if checks is not None:
                result["chematic_raw_geometry"] = [
                    check_geometry(conf, checks) for conf in raw
                ]
        start = time.perf_counter()
        method = (
            other.rdkit_embed_multiple_confs_mmff
            if profile == "compatible"
            else other.etkdg_embed_multiple_confs_mmff_corrected
        )
        kwargs = {"random_seed": seed, "max_mmff_iters": 200}
        if profile == "corrected":
            kwargs["repair_mmff_iters"] = repair_iters
        candidates = method(count, **kwargs)
        result["chematic_seconds"] = time.perf_counter() - start
        result["chematic_count"] = len(candidates)
        result["chematic_statuses"] = [status for status, _, _ in candidates]
        result["chematic_converged"] = sum(status == 0 for status, _, _ in candidates)
        result["chematic_rmsd"] = rms_values(
            mol, heavy_ref, [coords for _, _, coords in candidates]
        )
        result["chematic_stereo"] = stereo_violations(
            mol, [coords for _, _, coords in candidates]
        )
        if checks is not None:
            result["chematic_geometry"] = [
                check_geometry(coords, checks) for _, _, coords in candidates
            ]
        result["chematic_energies"] = [energy for _, energy, _ in candidates]
    except Exception as exc:
        result["chematic_error"] = repr(exc)
    generated = Chem.Mol(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.numThreads = 1
    start = time.perf_counter()
    try:
        ids = list(AllChem.EmbedMultipleConfs(generated, numConfs=count, params=params))
        raw_time = time.perf_counter() - start
        raw_coords = [
            [
                list(generated.GetConformer(cid).GetAtomPosition(i))
                for i in range(generated.GetNumAtoms())
            ]
            for cid in ids
        ]
        minimized = (
            AllChem.MMFFOptimizeMoleculeConfs(generated, numThreads=1, maxIters=200)
            if ids else []
        )
        result["rdkit_seconds"] = time.perf_counter() - start
        result["rdkit_raw_seconds"] = raw_time
        result["rdkit_raw_rmsd"] = rms_values(mol, heavy_ref, raw_coords)
        result["rdkit_raw_stereo"] = stereo_violations(mol, raw_coords)
        if checks is not None:
            result["rdkit_raw_geometry"] = [
                check_geometry(conf, checks) for conf in raw_coords
            ]
        assert len(ids) == len(minimized)
        coords = [
            [
                list(generated.GetConformer(cid).GetAtomPosition(i))
                for i in range(generated.GetNumAtoms())
            ]
            for cid in ids
        ]
        result["rdkit_count"] = len(ids)
        result["rdkit_statuses"] = [status for status, _ in minimized]
        result["rdkit_converged"] = sum(status == 0 for status, _ in minimized)
        result["rdkit_rmsd"] = rms_values(mol, heavy_ref, coords)
        result["rdkit_stereo"] = stereo_violations(mol, coords)
        if checks is not None:
            result["rdkit_geometry"] = [check_geometry(conf, checks) for conf in coords]
        result["rdkit_energies"] = [energy for _, energy in minimized]
    except Exception as exc:
        result["rdkit_error"] = repr(exc)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sdf", type=Path, required=True)
    parser.add_argument("--split", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--partition", choices=("development", "holdout"), default="development")
    parser.add_argument("--invalid-exclusion-audit", type=Path)
    parser.add_argument("--count", type=int, default=20)
    parser.add_argument("--limit", type=int, default=20)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--profile", choices=("compatible", "corrected"), default="compatible"
    )
    parser.add_argument("--geometry", action="store_true")
    parser.add_argument("--repair-iters", type=int, default=0)
    parser.add_argument("--independent-mmff", action="store_true",
                        help="attempt and record each engine even when RDKit lacks MMFF parameters")
    args = parser.parse_args()
    if args.count < 1 or args.limit < 1 or args.repair_iters < 0:
        parser.error("count and limit must be positive and repair-iters nonnegative")
    if args.profile == "compatible" and args.repair_iters:
        parser.error("repair-iters requires the corrected profile")
    split = json.loads(args.split.read_text())
    assert hashlib.sha256(args.sdf.read_bytes()).hexdigest() == split["source_sha256"]
    if args.partition == "holdout":
        if args.invalid_exclusion_audit is None:
            parser.error("holdout requires an invalid-exclusion graph audit")
        audit = json.loads(args.invalid_exclusion_audit.read_text())
        if (audit["split_sha256"] != hashlib.sha256(args.split.read_bytes()).hexdigest()
            or audit["source_sha256"] != split["source_sha256"]
            or audit["invalid_exclusion_count"] != len(split["invalid_exclusions"])
            or audit["overlaps"]):
            parser.error("invalid-exclusion audit does not establish a clean split")
    rows = [r for r in split["rows"] if r.get("split") == args.partition]
    rows.sort(key=lambda r: hashlib.sha256(r["name"].encode()).digest())
    selected = {r["index"]: r for r in rows[: args.limit]}
    supplier = Chem.SDMolSupplier(str(args.sdf), removeHs=False)
    output = {
        "rdkit": rdBase.rdkitVersion,
        "chematic": importlib.metadata.version("chematic"),
        "platform": platform.platform(),
        "split_sha256": hashlib.sha256(args.split.read_bytes()).hexdigest(),
        "source_sha256": split["source_sha256"],
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "profile": args.profile,
        "geometry_checked": args.geometry,
        "repair_iters": args.repair_iters,
        "candidate_limit": args.count,
        "mmff_iterations": 200,
        "independent_mmff": args.independent_mmff,
        "seed": args.seed,
        "selection": f"{args.partition} only sorted by SHA256(name)",
        "results": [],
    }
    for index in sorted(selected):
        row = selected[index]
        try:
            reference = supplier[index]
            assert reference is not None and reference.GetProp("_Name") == row["name"]
            result = score(
                reference,
                args.count,
                args.seed,
                args.profile,
                args.geometry,
                args.repair_iters,
                args.independent_mmff,
            )
        except Exception as exc:
            result = {"status": "error", "error": repr(exc)}
        result.update(index=index, name=row["name"])
        output["results"].append(result)
        checkpoint = args.output.with_suffix(args.output.suffix + ".tmp")
        checkpoint.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
        checkpoint.replace(args.output)
        print(index, row["name"], result["status"], flush=True)


if __name__ == "__main__":
    main()
