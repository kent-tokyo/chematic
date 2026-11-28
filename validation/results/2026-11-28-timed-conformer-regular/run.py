"""Same-input timed conformer check on a fixed exposed development subset."""
import argparse
import hashlib
import json
import platform
from pathlib import Path

import chematic

root = Path(__file__).resolve().parents[3]
source_path = root / "validation/results/2026-11-28-timed-preview-regular/results.json"
source = json.loads(source_path.read_text(encoding="utf-8"))
rows = source["results"][:12]
assert len(rows) == 12 and all(row["status"] == "embedded" for row in rows)
parser = argparse.ArgumentParser()
parser.add_argument("mode", choices=("before", "after"))
mode = parser.parse_args().mode
results = []
for row in rows:
    for partial in (False, True):
        options = {"timeout_seconds": 20.0, "random_seed": 42,
                   "max_iterations": 0, "max_mmff_iters": 30,
                   "repair_mmff_iters": 0, "return_partial": partial}
        try:
            plain = chematic.generate_conformers_timed(row["smiles"], 1, **options)
            item = {"index": row["index"], "partial": partial, "plain": plain}
            if mode == "after":
                item["screened"] = chematic.generate_conformers_timed(
                    row["smiles"], 1, conformer_require_regular_geometry=True, **options)
            results.append(item)
        except Exception as exc:
            results.append({"index": row["index"], "partial": partial,
                            "error": f"{type(exc).__name__}: {exc}"})

print(json.dumps({"mode": mode, "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
                  "selection": "first 12 from preceding fixed 30-molecule development sample",
                  "platform": platform.platform(), "python": platform.python_version(),
                  "options": {"num_confs": 1, "random_seed": 42, "max_iterations": 0,
                              "max_mmff_iters": 30, "repair_mmff_iters": 0,
                              "timeout_seconds": 20.0},
                  "results": results}, indent=2, allow_nan=False))
