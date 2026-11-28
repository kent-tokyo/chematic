"""Compare opt-in clash rejection with an independent pair-distance screen."""
import hashlib
import json
import math
import platform
from pathlib import Path

import chematic

root = Path(__file__).resolve().parents[3]
input_path = root / "validation/results/2026-11-27-selenium-crystal-stereo/source-development-input.json"
source = json.loads(input_path.read_text())
results = []
for row in source["results"]:
    smiles = row["smiles"]
    if not any(marker in smiles for marker in ("@", "/", "\\")):
        results.append({"index": row["index"], "status": "no_stereo_marker"})
        continue
    try:
        mol = chematic.from_smiles(smiles)
        before = mol.generate_3d_etkdg()
        initial = mol.stereo_3d_diagnostics(before)
        if not initial["violated"]:
            results.append({"index": row["index"], "status": "no_reported_violation",
                            "reported": initial["reported"], "unevaluable": initial["unevaluable"]})
            continue
        plain = mol.repair_stereo_3d(before)
        guarded = mol.repair_stereo_3d(before, clash_cutoff=0.8)
        record = {"index": row["index"], "status": "evaluated", "plain_ok": plain["ok"],
                  "guarded_ok": guarded["ok"], "check": guarded["clash_check"],
                  "failures": guarded["failures"]}
        if plain["ok"]:
            bonds = {tuple(sorted((a, b))) for a, b, *_ in mol.bond_table}
            def pairs(coords):
                return {(i, j) for i in range(len(coords)) for j in range(i + 1, len(coords))
                        if (i, j) not in bonds and math.dist(coords[i], coords[j]) < 0.8}
            before_pairs, after_pairs = pairs(before), pairs(plain["coordinates"])
            independent_new = sorted(after_pairs - before_pairs)
            reported_new = sorted(tuple(item["atoms"]) for item in guarded["clash_check"]["new"])
            record.update(independent_before_count=len(before_pairs),
                          independent_after_count=len(after_pairs),
                          independent_new=independent_new,
                          pair_agreement=independent_new == reported_new,
                          decision_agreement=guarded["ok"] == (not independent_new))
        results.append(record)
    except Exception as exc:
        results.append({"index": row["index"], "status": "error",
                        "error": f"{type(exc).__name__}: {exc}"})

print(json.dumps({"input_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
                  "cutoff_angstrom": 0.8, "platform": platform.platform(),
                  "python": platform.python_version(), "results": results},
                 indent=2, allow_nan=False))
