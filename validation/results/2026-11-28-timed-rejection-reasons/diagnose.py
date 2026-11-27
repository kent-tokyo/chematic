"""Exercise rejection classification using development rule-based geometries.

These are simulated candidate rows, not timed-ensemble results.
"""
import hashlib
import json
import platform
from pathlib import Path

import chematic
from chematic._timed_conformers import _filter_regular_geometry

root = Path(__file__).resolve().parents[3]
source_path = root / "validation/results/2026-11-28-timed-preview-regular/results.json"
fixture_path = root / "crates/chematic-py/tests/stereo_distorted_prototype_719.json"
source = json.loads(source_path.read_text(encoding="utf-8"))
fixture = json.loads(fixture_path.read_text(encoding="utf-8"))
results = []
for row in source["results"]:
    try:
        mol = chematic.from_smiles(row["smiles"])
        coords = mol.generate_3d_etkdg()
        retained, rejected, reasons = _filter_regular_geometry(mol, [[1, 0.0, coords, None]])
        report = mol.stereo_3d_diagnostics(coords, require_regular_geometry=True)
        results.append({"index": row["index"], "retained": len(retained), "rejected": rejected,
                        "reasons": reasons, "report": report})
    except Exception as exc:
        results.append({"index": row["index"], "error": f"{type(exc).__name__}: {exc}"})

controls = [
    ("distorted", fixture["smiles"], fixture["coords"]),
    ("violated", "F/C=C/F", [[0.,1.,0.],[0.,0.,0.],[1.,0.,0.],[1.,1.,0.]]),
    ("unsupported", "[Pt@SP1](F)(Cl)(Br)I", [[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[-1.,0.,0.],[0.,-1.,0.]]),
]
control_results = []
for name, smiles, coords in controls:
    try:
        mol = chematic.from_smiles(smiles)
        retained, rejected, reasons = _filter_regular_geometry(mol, [[1, 0.0, coords, None]])
        control_results.append({"name": name, "retained": len(retained), "rejected": rejected,
                                "reasons": reasons})
    except Exception as exc:
        control_results.append({"name": name, "error": f"{type(exc).__name__}: {exc}"})

print(json.dumps({"source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
                  "fixture_sha256": hashlib.sha256(fixture_path.read_bytes()).hexdigest(),
                  "platform": platform.platform(), "python": platform.python_version(),
                  "results": results, "controls": control_results}, indent=2, allow_nan=False))
