"""Compare implicit and reconstructed-explicit H on exposed development inputs."""
import hashlib
import json
import math
import platform
from pathlib import Path
import chematic

root = Path(__file__).resolve().parents[3]
source_path = root / "validation/results/2026-11-27-selenium-crystal-stereo/source-development-input.json"
source = json.loads(source_path.read_text())
results = []
for row in source["results"]:
    if "@" not in row["smiles"]:
        results.append({"index": row["index"], "status": "no_stereo_marker"})
        continue
    try:
        mol = chematic.from_smiles(row["smiles"])
        coords = mol.generate_3d_etkdg()
        n = mol.atom_count
        original = mol.atom_table
        neighbors = [[] for _ in range(n)]
        for a, b, *_ in mol.bond_table:
            neighbors[a].append(b)
            neighbors[b].append(a)
        targets = [i for i, atom in enumerate(original)
                   if atom[1] == 6 and atom[4] == 1 and len(neighbors[i]) == 3
                   and all(original[j][1] != 1 for j in neighbors[i])]
        implicit = mol.stereo_3d_diagnostics(coords, require_regular_geometry=True)
        reported = {x["atom"]: x for x in implicit["atom_stereo"]}
        targets = [i for i in targets if i in reported]
        if not targets:
            results.append({"index": row["index"], "status": "no_target"})
            continue
        expanded = mol.add_hydrogens()
        atoms = expanded.atom_table
        assert [(a[0], a[1], a[2]) for a in atoms[:n]] == [(a[0], a[1], a[2]) for a in original], "original atom order"
        added_parent = {}
        for a, b, *_ in expanded.bond_table:
            if a >= n and b < n:
                added_parent[a] = b
            elif b >= n and a < n:
                added_parent[b] = a
        extra = []
        skipped = {}
        target_set = set(targets)
        for idx in range(n, expanded.atom_count):
            assert atoms[idx][1] == 1 and idx in added_parent, "added hydrogen mapping"
            parent = added_parent[idx]
            p = coords[parent]
            direction = [0., 0., 1.]
            if parent in target_set:
                unit = []
                for child in neighbors[parent]:
                    vector = [coords[child][i] - p[i] for i in range(3)]
                    length = math.sqrt(sum(x*x for x in vector))
                    if length < 1e-10:
                        skipped[parent] = "zero_bond_direction"
                        break
                    unit.append([x / length for x in vector])
                if len(unit) == 3:
                    direction = [-sum(v[i] for v in unit) for i in range(3)]
                    length = math.sqrt(sum(x*x for x in direction))
                    if length < 1e-3:
                        skipped[parent] = "undefined_phantom_direction"
                        direction = [0., 0., 1.]
                    else:
                        direction = [x / length for x in direction]
            extra.append([p[i] + direction[i] for i in range(3)])
        explicit = expanded.stereo_3d_diagnostics(coords + extra, require_regular_geometry=True)
        explicit_report = {x["atom"]: x for x in explicit["atom_stereo"]}
        centers = []
        for idx in targets:
            item = {"atom": idx, "implicit": reported[idx], "skip": skipped.get(idx)}
            if idx not in skipped:
                assert list(added_parent.values()).count(idx) == 1, "target hydrogen count"
                item["explicit"] = explicit_report.get(idx)
            centers.append(item)
        results.append({"index": row["index"], "status": "evaluated", "centers": centers})
    except Exception as exc:
        results.append({"index": row["index"], "status": "error", "error": f"{type(exc).__name__}: {exc}"})

print(json.dumps({"input_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
                  "platform": platform.platform(), "python": platform.python_version(),
                  "results": results}, indent=2, allow_nan=False))
