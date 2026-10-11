"""Audit source-derived OOP lookup models against the fixed live RDKit oracle.

This models the preceding exact/wildcard resolver and the equivalence resolver;
it does not call Rust's parameter API. The energy run/tests exercise production.
Run from the repository root in the source diagnostic Python environment.
"""

import hashlib
import itertools
import json
import re
from collections import Counter
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

ROOT = Path("validation/results/mmff94-correctness-2026-10-11")


def read(path):
    raw = Path(path).read_bytes()
    return raw, {"path": str(path), "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}


def main():
    assert rdBase.rdkitVersion == "2026.03.6"
    raw, corpus = read(ROOT / "torsion/candidate.jsonl")
    rows = [json.loads(line) for line in raw.splitlines()]
    raw, parameters = read("crates/chematic-ff/src/mmff94_energy/oop_stbn.rs")
    data = raw.decode().split("pub static MMFF94_OOP:")[1].split("];")[0]
    table = {
        (int(i), int(j), int(k), int(l)): float(value)
        for i, j, k, l, value in re.findall(r"\((\d+), (\d+), (\d+), (\d+), ([-\d.]+)\)", data)
    }
    raw, metadata = read("crates/chematic-ff/src/mmff94_numeric_type_registry.rs")
    equivalence = {
        int(i): list(map(int, re.findall(r"\d+", levels)))
        for i, levels in re.findall(r"id: (\d+),.*?equivalence_levels: \[([^]]+)\]", raw.decode(), re.S)
    }

    def previous(center, outer):
        for i, k, l in itertools.permutations(outer):
            if (i, center, k, l) in table:
                return table[i, center, k, l]
        for i, k, l in [(0, 0, 0), (outer[0], 0, 0), (0, outer[1], 0), (0, 0, outer[2])]:
            if (i, center, k, l) in table:
                return table[i, center, k, l]
        return None

    def resolved(center, outer):
        for level in range(4):
            i, k, l = sorted(equivalence[t][level] for t in outer)
            if (i, center, k, l) in table:
                return table[i, center, k, l]
        return None

    queries = 0
    before = Counter()
    after = []
    centers = Counter()
    for row in rows:
        if row["status"] != "ok":
            continue
        mol = Chem.AddHs(Chem.MolFromSmiles(row["smiles"]))
        properties = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
        assert properties is not None
        for atom in mol.GetAtoms():
            if atom.GetDegree() != 3:
                continue
            i, k, l = [nb.GetIdx() for nb in atom.GetNeighbors()]
            j = atom.GetIdx()
            center = properties.GetMMFFAtomType(j)
            outer = [properties.GetMMFFAtomType(x) for x in (i, k, l)]
            oracle = properties.GetMMFFOopBendParams(mol, i, j, k, l)
            queries += 1
            if oracle is not None:
                centers[center] += 1
            old = previous(center, outer)
            if old != oracle:
                before[center, tuple(sorted(outer)), old, oracle] += 1
            new = resolved(center, outer)
            if new != oracle:
                after.append({"input_index": row["input_index"], "center": center, "outer": outer, "model": new, "oracle": oracle})
    result = {
        "kind": "source_derived_lookup_models_not_a_direct_rust_api_measurement",
        "rdkit_runtime": rdBase.rdkitVersion,
        "corpus": corpus,
        "parameter_source": parameters,
        "metadata_source": metadata,
        "comparable_inputs": sum(row["status"] == "ok" for row in rows),
        "degree_three_centers_queried": queries,
        "oracle_parameterized_centers_by_type": dict(centers),
        "before_mismatch_instances": sum(before.values()),
        "before_mismatch_tuples": [
            {"center": center, "outer": outer, "before": old, "oracle": oracle, "instances": count}
            for (center, outer, old, oracle), count in sorted(before.items())
        ],
        "after_mismatches": after,
    }
    (ROOT / "oop-vdw/oop-parameter-audit.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    assert not after
    print(f"{queries} centers queried; {sum(before.values())} previous mismatches; 0 equivalence-model mismatches")


if __name__ == "__main__":
    main()
