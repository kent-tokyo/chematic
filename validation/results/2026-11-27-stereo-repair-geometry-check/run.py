"""Compare opt-in post-repair geometry checks on exposed development cases."""
import hashlib
import json
import math
import platform
from pathlib import Path

import chematic
from rdkit import Chem, rdBase

rdBase.DisableLog("rdApp.*")
root = Path(__file__).resolve().parents[3]
source_path = root / "validation/results/2026-11-27-selenium-crystal-stereo/source-development-input.json"
selection_path = root / "validation/results/2026-11-27-implicit-h-planarity/alternatives-after.json"
source = json.loads(source_path.read_text())
selection = json.loads(selection_path.read_text())
selected = {row["index"]: row for row in selection["results"]}
rows = [row for row in source["results"] if row["index"] in selected]
assert len(rows) == len(selected) == 576

tetra = {Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW}
results = []
for row in rows:
    try:
        mol = chematic.from_smiles(row["smiles"])
        before = mol.generate_3d_etkdg()
        kwargs = {"clash_cutoff": 0.8, "max_alternative_attempts": 3}
        plain = mol.repair_stereo_3d(before, **kwargs)
        checked = mol.repair_stereo_3d(before, require_regular_geometry=True, **kwargs)
        result = {"index": row["index"], "baseline_ok": plain["ok"],
                  "checked_ok": checked["ok"], "fully_verified": checked["fully_verified"],
                  "baseline_matches_previous": plain["ok"] == selected[row["index"]]["candidate_ok"],
                  "retry": checked["retry"], "clash_check": checked["clash_check"],
                  "regular_geometry_check": checked["regular_geometry_check"],
                  "failures": checked["failures"]}
        if plain["ok"] and checked["ok"]:
            result["unchanged_from_baseline"] = plain["coordinates"] == checked["coordinates"]
        if checked["ok"]:
            coords = checked["coordinates"]
            result["post_screen"] = mol.stereo_3d_diagnostics(coords, require_regular_geometry=True)
            bonds = {tuple(sorted((a,b))) for a,b,*_ in mol.bond_table}
            result["max_bond_length_change"] = max(
                (abs(math.dist(before[a], before[b]) - math.dist(coords[a], coords[b]))
                 for a,b in bonds), default=0.0)
            def pairs(x):
                return {(i,j) for i in range(len(x)) for j in range(i+1,len(x))
                        if (i,j) not in bonds and math.dist(x[i],x[j]) < 0.8}
            result["new_clashes"] = sorted(pairs(coords) - pairs(before))
            ref = Chem.MolFromSmiles(row["smiles"])
            assert ref is not None
            assert [(a.GetSymbol(), a.GetAtomicNum(), a.GetFormalCharge()) for a in ref.GetAtoms()] == [(a[0], a[1], a[2]) for a in mol.atom_table], "atom mapping"
            assert {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in ref.GetBonds()} == bonds, "bond mapping"
            declared = {a.GetIdx(): a.GetChiralTag() for a in ref.GetAtoms() if a.GetChiralTag() in tetra}
            conf = Chem.Conformer(ref.GetNumAtoms())
            for i, point in enumerate(coords):
                conf.SetAtomPosition(i, point)
            ref.AddConformer(conf)
            Chem.AssignStereochemistryFrom3D(ref, replaceExistingTags=True)
            reports = {r["atom"]: r for r in result["post_screen"]["atom_stereo"]}
            centers = []
            for idx, original in declared.items():
                observed = ref.GetAtomWithIdx(idx).GetChiralTag()
                item = {"atom": idx, "declared": str(original), "rdkit": str(observed), "screened": reports.get(idx)}
                if observed in tetra and idx in reports and reports[idx]["status"] in ("satisfied", "violated"):
                    item["agrees"] = (reports[idx]["status"] == "satisfied") == (observed == original)
                centers.append(item)
            result["rdkit_centers"] = centers
        results.append(result)
    except Exception as exc:
        results.append({"index": row["index"], "error": f"{type(exc).__name__}: {exc}"})

print(json.dumps({"source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
                  "selection_sha256": hashlib.sha256(selection_path.read_bytes()).hexdigest(),
                  "platform": platform.platform(), "python": platform.python_version(),
                  "rdkit": rdBase.rdkitVersion, "cutoff": 0.8, "max_alternative_attempts": 3,
                  "results": results}, indent=2, allow_nan=False))
