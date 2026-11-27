"""Compare declared stereo diagnostics against RDKit's 3D chiral tags."""
import hashlib
import json
import platform
from pathlib import Path

import chematic
from rdkit import Chem, rdBase

rdBase.DisableLog("rdApp.*")
root = Path(__file__).resolve().parents[3]
source_path = root / "validation/results/2026-11-27-selenium-crystal-stereo/source-development-input.json"
source = json.loads(source_path.read_text())
tetra = {Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW}
results = []
for row in source["results"]:
    if "@" not in row["smiles"]:
        results.append({"index": row["index"], "status": "no_atom_stereo_marker"})
        continue
    try:
        smiles = row["smiles"]
        mol = chematic.from_smiles(smiles)
        ref = Chem.MolFromSmiles(smiles)
        if ref is None:
            results.append({"index": row["index"], "status": "rdkit_parse_failed"})
            continue
        assert [(a.GetSymbol(), a.GetAtomicNum(), a.GetFormalCharge()) for a in ref.GetAtoms()] == [(a[0], a[1], a[2]) for a in mol.atom_table], "atom mapping"
        assert {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in ref.GetBonds()} == {tuple(sorted((a, b))) for a, b, *_ in mol.bond_table}, "bond mapping"
        expected = {a.GetIdx(): a.GetChiralTag() for a in ref.GetAtoms() if a.GetChiralTag() in tetra}
        coords = mol.generate_3d_etkdg()
        raw = mol.stereo_3d_diagnostics(coords)
        screened = mol.stereo_3d_diagnostics(coords, require_regular_geometry=True)
        conformer = Chem.Conformer(ref.GetNumAtoms())
        for i, point in enumerate(coords):
            conformer.SetAtomPosition(i, point)
        ref.AddConformer(conformer)
        Chem.AssignStereochemistryFrom3D(ref, replaceExistingTags=True)
        r = {entry["atom"]: entry for entry in raw["atom_stereo"]}
        t = {entry["atom"]: entry for entry in screened["atom_stereo"]}
        centers = []
        for idx, declared in expected.items():
            observed = ref.GetAtomWithIdx(idx).GetChiralTag()
            item = {"atom": idx, "declared_tag": str(declared), "rdkit_tag": str(observed),
                    "raw": r.get(idx), "screened": t.get(idx)}
            if observed in tetra and idx in r and r[idx]["status"] in ("satisfied", "violated"):
                item["rdkit_agrees_with_declaration"] = observed == declared
                item["raw_agrees_with_rdkit"] = (r[idx]["status"] == "satisfied") == (observed == declared)
            centers.append(item)
        results.append({"index": row["index"], "status": "evaluated", "centers": centers})
    except Exception as exc:
        results.append({"index": row["index"], "status": "error", "error": f"{type(exc).__name__}: {exc}"})

print(json.dumps({"input_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
                  "platform": platform.platform(), "python": platform.python_version(),
                  "rdkit": rdBase.rdkitVersion, "results": results}, indent=2, allow_nan=False))
