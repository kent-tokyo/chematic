#!/usr/bin/env python3
"""Compare timed checked ensembles with and without provisional embeddings."""
import argparse
import hashlib
import json
import platform
from pathlib import Path

import chematic
from chematic import chematic as extension
from rdkit import Chem, rdBase

from platinum_development_probe import graph_matches, rms_values
from platinum_geometry import check_geometry, prepare_geometry_checks
from platinum_stereo_strict import stereo_violations_strict


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def screen(molecule, reference, conformers, checks):
    if not conformers:
        return []
    distances = rms_values(molecule, reference, conformers)
    stereo = stereo_violations_strict(molecule, conformers)
    return [{"rmsd": rms, "stereo": state, "geometry": check_geometry(coords, checks)}
            for rms, state, coords in zip(distances, stereo, conformers)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("sdf", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seconds", type=float, default=0.2)
    args = parser.parse_args()
    source = json.loads(args.input.read_text(encoding="utf-8"))
    if digest(args.sdf) != source["source_sha256"]:
        parser.error("reference SDF hash differs")
    supplier = Chem.SDMolSupplier(str(args.sdf), removeHs=False)
    results = []
    for index, row in enumerate(source["results"]):
        result = {"index": row["index"], "name": row["name"]}
        try:
            reference = supplier[row["index"]]
            if reference is None or reference.GetProp("_Name") != row["name"]:
                raise ValueError("reference parse or identity mismatch")
            reference = Chem.RemoveHs(reference)
            molecule = Chem.AddHs(Chem.MolFromSmiles(row["smiles"]))
            other = chematic.from_smiles(row["smiles"]).add_hydrogens()
            if not graph_matches(other, molecule):
                raise ValueError("atom or bond graph mismatch")
            checks = prepare_geometry_checks(molecule)
            result["geometry_coverage"] = {
                "bonds_supported": len(checks[0]), "bonds_total": molecule.GetNumBonds(),
                "angles_supported": len(checks[1]),
                "angles_total": sum(a.GetDegree() * (a.GetDegree() - 1) // 2 for a in molecule.GetAtoms())}
            order = (False, True) if index % 2 == 0 else (True, False)
            result["order"] = ["preview" if enabled else "partial" for enabled in order]
            for enabled in order:
                response = chematic.generate_conformers_timed(
                    row["smiles"], 20, timeout_seconds=args.seconds, return_partial=True,
                    return_preview=enabled, max_mmff_iters=200, repair_mmff_iters=0)
                conformers = [candidate[2] for candidate in response.get("conformers", [])]
                preview = response.get("in_progress")
                result["preview" if enabled else "partial"] = {
                    "status": response["status"], "elapsed_seconds": response["elapsed_seconds"],
                    "attempts_completed": response.get("attempts_completed"),
                    "conformers": response.get("conformers", []), "in_progress": preview,
                    "completed_screen": screen(molecule, reference, conformers, checks),
                    "preview_screen": screen(molecule, reference, [preview["coordinates"]], checks) if preview else [],
                    "error": response.get("error"), "checkpoint_error": response.get("checkpoint_error")}
        except Exception as error:
            result["error"] = repr(error)
        results.append(result)
        print(f"{len(results)}/{len(source['results'])}", flush=True)
    module = Path(chematic.__file__).parent / "_timed_conformers.py"
    output = {"input_sha256": digest(args.input), "sdf_sha256": digest(args.sdf),
              "wrapper_sha256": digest(module), "extension_sha256": digest(Path(extension.__file__)),
              "script_sha256": digest(Path(__file__)), "platform": platform.platform(),
              "python": platform.python_version(), "rdkit": rdBase.rdkitVersion,
              "seconds": args.seconds, "num_confs": 20, "results": results}
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
