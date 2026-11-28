"""Evaluate stereo repair on the previously selected exposed COD development subset."""
import hashlib
import json
import math
import platform
from pathlib import Path

import chematic

root = Path(__file__).resolve().parents[3]
source_path = root / "validation/results/2026-11-27-selenium-crystal-stereo/source-development-input.json"
selection_path = root / "validation/results/2026-11-27-boron-initial-geometry/before.json"
source = json.loads(source_path.read_text())
selection = json.loads(selection_path.read_text())
indices = {row["index"] for row in selection["results"]}
rows = [row for row in source["results"] if row["index"] in indices]
assert len(rows) == len(indices) == 24


def distance(a, b):
    return math.dist(a, b)


def clash_count(coords, bonds):
    return sum(distance(coords[i], coords[j]) < 0.8
               for i in range(len(coords)) for j in range(i + 1, len(coords))
               if (i, j) not in bonds)


results = []
for row in rows:
    try:
        mol = chematic.from_smiles(row["smiles"])
        before = mol.generate_3d_etkdg()
        bonds = {tuple(sorted((a, b))) for a, b, *_ in mol.bond_table}
        diagnostic_before = mol.stereo_3d_diagnostics(before)
        repair = mol.repair_stereo_3d(before)
        record = {"index": row["index"], "before": diagnostic_before,
                  "ok": repair["ok"], "fully_verified": repair["fully_verified"],
                  "repaired": repair["repaired"], "unevaluable": repair["unevaluable"],
                  "failures": repair["failures"], "clashes_before": clash_count(before, bonds)}
        if repair["ok"]:
            after = repair["coordinates"]
            record["after"] = mol.stereo_3d_diagnostics(after)
            record["clashes_after"] = clash_count(after, bonds)
            record["max_bond_length_change"] = max(
                (abs(distance(before[a], before[b]) - distance(after[a], after[b]))
                 for a, b in bonds), default=0.0)
        results.append(record)
    except Exception as exc:
        results.append({"index": row["index"], "error": f"{type(exc).__name__}: {exc}"})

output = {"input_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
          "selection_sha256": hashlib.sha256(selection_path.read_bytes()).hexdigest(),
          "platform": platform.platform(), "python": platform.python_version(),
          "selection_count": len(indices), "results": results}
print(json.dumps(output, indent=2, allow_nan=False))
