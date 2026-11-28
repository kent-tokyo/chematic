"""Compare legacy and opt-in preview checks on fixed exposed development inputs."""
import hashlib
import json
import platform
from pathlib import Path

import chematic

root = Path(__file__).resolve().parents[3]
source_path = root / "validation/results/2026-11-27-selenium-crystal-stereo/source-development-input.json"
source = json.loads(source_path.read_text(encoding="utf-8"))
selected = []
selection_errors = []
for row in sorted(source["results"], key=lambda r: r["index"]):
    if "@" not in row["smiles"]:
        continue
    try:
        mol = chematic.from_smiles(row["smiles"])
        heavy_atoms = sum(a[1] != 1 for a in mol.atom_table)
        if heavy_atoms <= 22:
            selected.append({"index": row["index"], "smiles": row["smiles"],
                             "heavy_atoms": heavy_atoms})
            if len(selected) == 30:
                break
    except Exception as exc:
        selection_errors.append({"index": row["index"], "error": f"{type(exc).__name__}: {exc}"})
assert len(selected) == 30
results = []
for row in selected:
    try:
        mol = chematic.from_smiles(row["smiles"]).add_hydrogens()
        coords = mol._timed_prepare_corrected(max_iterations=0).embed(42)
        if coords is None:
            results.append({**row, "status": "embedding_failed"})
            continue
        baseline = mol._timed_stereo_valid(coords)
        strict = mol._timed_stereo_valid(coords, require_regular_geometry=True)
        raw_report = mol.stereo_3d_diagnostics(coords)
        strict_report = mol.stereo_3d_diagnostics(coords, require_regular_geometry=True)
        assert baseline == (raw_report["violated"] == 0 and raw_report["unevaluable"] == 0)
        assert strict == (strict_report["violated"] == 0 and strict_report["unevaluable"] == 0)
        results.append({**row, "status": "embedded", "baseline_preview_valid": baseline,
                        "strict_preview_valid": strict, "raw": raw_report,
                        "strict": strict_report, "coordinates": coords})
    except Exception as exc:
        results.append({**row, "status": "error", "error": f"{type(exc).__name__}: {exc}"})

print(json.dumps({"source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
                  "selection": "first 30 ascending source indices with @ and <=22 parsed heavy atoms",
                  "selection_errors": selection_errors, "seed": 42,
                  "platform": platform.platform(), "python": platform.python_version(),
                  "results": results}, indent=2, allow_nan=False))
