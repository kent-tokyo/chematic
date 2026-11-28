"""Compare bounded stereo repair alternatives on the exposed development partition."""
import hashlib
import json
import math
import platform
from pathlib import Path

import chematic
from rdkit import Chem, rdBase

root = Path(__file__).resolve().parents[3]
source_path = root / "validation/results/2026-11-27-selenium-crystal-stereo/source-development-input.json"
selection_path = root / "validation/results/2026-11-27-stereo-repair-clashes/results.json"
source = json.loads(source_path.read_text())
selection = json.loads(selection_path.read_text())
selected = {row["index"]: row for row in selection["results"] if row["status"] == "evaluated"}
rows = [row for row in source["results"] if row["index"] in selected]
assert len(rows) == len(selected) == 576


def pair_clashes(coords, bonds):
    return {(i, j) for i in range(len(coords)) for j in range(i+1, len(coords))
            if (i,j) not in bonds and math.dist(coords[i],coords[j]) < 0.8}


results = []
for row in rows:
    try:
        mol = chematic.from_smiles(row["smiles"])
        before = mol.generate_3d_etkdg()
        baseline = mol.repair_stereo_3d(before, clash_cutoff=0.8)
        candidate = mol.repair_stereo_3d(before, clash_cutoff=0.8, max_alternative_attempts=3)
        record = {"index": row["index"], "baseline_ok": baseline["ok"],
                  "candidate_ok": candidate["ok"], "retry": candidate["retry"],
                  "baseline_matches_previous": baseline["ok"] == selected[row["index"]]["guarded_ok"],
                  "failures": candidate["failures"]}
        if candidate["ok"]:
            after = candidate["coordinates"]
            bonds = {tuple(sorted((a,b))) for a,b,*_ in mol.bond_table}
            new_pairs = sorted(pair_clashes(after,bonds) - pair_clashes(before,bonds))
            stereo = mol.stereo_3d_diagnostics(after)
            record.update(new_pairs=new_pairs, stereo=stereo,
                          max_bond_length_change=max((abs(math.dist(before[a],before[b]) - math.dist(after[a],after[b]))
                                                     for a,b in bonds), default=0.0),
                          coordinates_unchanged_from_baseline=(after == baseline["coordinates"] if baseline["ok"] else None))
            if not baseline["ok"]:
                reference = Chem.MolFromSmiles(row["smiles"])
                assert reference is not None and reference.GetNumAtoms() == len(after)
                assert [(a.GetSymbol(), a.GetAtomicNum()) for a in reference.GetAtoms()] == [(a[0], a[1]) for a in mol.atom_table]
                assert {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in reference.GetBonds()} == bonds
                Chem.AssignStereochemistry(reference, force=True, cleanIt=True)
                expected = {a.GetIdx(): a.GetProp("_CIPCode") for a in reference.GetAtoms() if a.HasProp("_CIPCode")}
                conformer = Chem.Conformer(reference.GetNumAtoms())
                for i, point in enumerate(after):
                    conformer.SetAtomPosition(i, point)
                reference.AddConformer(conformer)
                Chem.AssignStereochemistryFrom3D(reference, replaceExistingTags=True)
                Chem.AssignStereochemistry(reference, force=True, cleanIt=True)
                observed = {i: reference.GetAtomWithIdx(i).GetProp("_CIPCode") if reference.GetAtomWithIdx(i).HasProp("_CIPCode") else None for i in expected}
                record["rdkit_cip_expected"] = expected
                record["rdkit_cip_observed"] = observed
                record["rdkit_cip_agreement"] = bool(expected) and expected == observed
        results.append(record)
    except Exception as exc:
        results.append({"index": row["index"], "error": f"{type(exc).__name__}: {exc}"})

print(json.dumps({"input_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
                  "selection_sha256": hashlib.sha256(selection_path.read_bytes()).hexdigest(),
                  "cutoff_angstrom": 0.8, "max_alternative_attempts": 3,
                  "platform": platform.platform(), "python": platform.python_version(), "rdkit": rdBase.rdkitVersion,
                  "results": results}, indent=2, allow_nan=False))
