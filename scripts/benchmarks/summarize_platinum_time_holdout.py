#!/usr/bin/env python3
"""Report paired scaffold-cluster intervals for the full time-budget holdout."""
import argparse
import json
import random
from collections import defaultdict
from pathlib import Path

from summarize_platinum_time_budget import recovered


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input", type=Path)
    p.add_argument("split", type=Path)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--bootstrap", type=int, default=10000)
    p.add_argument("--support", type=Path,
                   help="independent MMFF support audit for a common-domain convergence analysis")
    p.add_argument("--groups", type=Path,
                   help="precomputed, outcome-independent holdout group assignments")
    p.add_argument("--require-strict", action="store_true",
                   help="require geometry and complete strict-stereo scoring")
    args = p.parse_args()
    data = json.loads(args.input.read_text())
    split = json.loads(args.split.read_text())
    if data["source_sha256"] != split["source_sha256"]:
        raise ValueError("source hash mismatch")
    expected = {r["index"]: r["scaffold_sha256"] for r in split["rows"] if r.get("split") == "holdout"}
    rows = data["results"]
    indices = [row["index"] for row in rows]
    if set(indices) != set(expected) or len(indices) != len(expected) or args.bootstrap < 1:
        raise ValueError("holdout is incomplete, duplicated, or bootstrap count is invalid")
    if args.require_strict:
        if not data.get("geometry_checked") or data.get("stereo_policy") not in (
            "strict", "rdkit-declared-stereo-strict-v1"
        ):
            raise ValueError("strict stereo and geometry scoring are required")
        for row in rows:
            if row.get("error"):
                continue
            for engine in ("chematic", "rdkit"):
                if engine in row.get("engine_errors", {}):
                    continue
                record = row.get(engine)
                if record is None:
                    raise ValueError("engine result is missing")
                count = record.get("accepted")
                if not isinstance(count, int) or count < 0 or any(
                    len(record.get(key, [])) != count for key in ("rmsd", "stereo", "geometry", "statuses")
                ):
                    raise ValueError("accepted conformer scoring is incomplete")
                required = {"atoms", "bonds", "violated_atoms", "violated_bonds",
                            "unevaluable_atoms", "unevaluable_bonds"}
                if any(not required.issubset(stereo) for stereo in record["stereo"]):
                    raise ValueError("strict stereo details are missing")
    geometry = data.get("geometry_checked", False)
    result = {"rows": len(rows), "seconds": data["seconds"], "bootstrap_samples": args.bootstrap, "bootstrap_seed": 42, "stages": {}}
    rng = random.Random(42)
    analyses = [("recovery", False, rows), ("converged_recovery", True, rows)]
    common = None
    if args.support:
        import hashlib
        support = json.loads(args.support.read_text())
        if (support.get("source_sha256") != data["source_sha256"]
                or support.get("split_sha256") != hashlib.sha256(args.split.read_bytes()).hexdigest()
                or support.get("partition") != "holdout"):
            raise ValueError("MMFF support audit provenance differs")
        support_rows = support.get("results", [])
        by_index = {row["index"]: row for row in support_rows}
        if len(by_index) != len(support_rows) or set(by_index) != set(expected):
            raise ValueError("MMFF support audit is incomplete or duplicated")
        if any(by_index[row["index"]].get("name") != row.get("name") for row in rows):
            raise ValueError("MMFF support audit identity differs")
        common = [row for row in rows if by_index[row["index"]].get("rdkit_mmff") is True
                  and by_index[row["index"]].get("chematic_mmff") is True
                  and by_index[row["index"]].get("graph_matches") is True]
        if not common:
            raise ValueError("common MMFF domain is empty")
        result["common_mmff_records"] = len(common)
        result["mmff_support_unknown_records"] = sum(
            any(not isinstance(row.get(key), bool) for key in
                ("rdkit_mmff", "chematic_mmff", "graph_matches"))
            for row in support_rows
        )
        result["support_sha256"] = hashlib.sha256(args.support.read_bytes()).hexdigest()
        analyses.append(("common_mmff_converged_recovery", True, common))
    if args.groups:
        import hashlib
        groups = json.loads(args.groups.read_text())
        if (groups.get("source_sha256") != data["source_sha256"]
                or groups.get("split_sha256") != hashlib.sha256(args.split.read_bytes()).hexdigest()):
            raise ValueError("group manifest provenance differs")
        group_rows = groups.get("rows", [])
        by_index = {row["index"]: row for row in group_rows}
        if len(by_index) != len(group_rows) or set(by_index) != set(expected):
            raise ValueError("group manifest is incomplete or duplicated")
        for row in rows:
            group = by_index[row["index"]]
            if (group.get("name") != row.get("name")
                    or group.get("scaffold_sha256") != expected[row["index"]]
                    or group.get("group") not in ("ordinary", "flexible", "macrocycle", "unclassified")):
                raise ValueError("group manifest identity or label differs")
        result["groups_sha256"] = hashlib.sha256(args.groups.read_bytes()).hexdigest()
        for name in ("ordinary", "flexible", "macrocycle", "unclassified"):
            selected = [row for row in rows if by_index[row["index"]]["group"] == name]
            if selected:
                analyses.append((f"{name}_recovery", False, selected))
                if common is not None:
                    common_indices = {row["index"] for row in common}
                    group_common = [row for row in selected if row["index"] in common_indices]
                    if group_common:
                        analyses.append((f"{name}_common_mmff_converged_recovery", True, group_common))
    for label, converged, selected in analyses:
        clusters = defaultdict(list)
        for row in selected:
            delta = int(recovered(row, "chematic", geometry, converged)) - int(recovered(row, "rdkit", geometry, converged))
            clusters[expected[row["index"]]].append(delta)
        totals = [(sum(values), len(values)) for values in clusters.values()]
        samples = []
        for _ in range(args.bootstrap):
            draw = [totals[rng.randrange(len(totals))] for _ in totals]
            samples.append(sum(value for value, _ in draw) / sum(size for _, size in draw))
        samples.sort()
        lo = int(0.025 * (args.bootstrap - 1))
        hi = int(0.975 * (args.bootstrap - 1))
        result["stages"][label] = {
            "records": len(selected),
            "chematic_recovered": sum(recovered(row, "chematic", geometry, converged) for row in selected),
            "rdkit_recovered": sum(recovered(row, "rdkit", geometry, converged) for row in selected),
            "difference": sum(v for v, _ in totals) / len(selected),
            "scaffold_clusters": len(totals),
            "paired_scaffold_cluster_95_interval": [samples[lo], samples[hi]],
        }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
