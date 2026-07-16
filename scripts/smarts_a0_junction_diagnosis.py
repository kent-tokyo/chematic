#!/usr/bin/env python3
"""SMARTS-A0: diagnose the 56-row / 28-molecule `c` / `[#6]=[#6]` residual left
after SMARTS-R1 (see docs/rdkit_compat.md).

Root cause: NOT a SMARTS-matcher bug. `chematic.from_smiles()` raw parsing
already matches RDKit's aromatic-atom-flag set exactly on every one of these
28 molecules. The mismatch is introduced entirely by the re-perception step
`mol.apply_aromaticity("rdkit_like")` (what `rdkit_compat.py`'s
`MolFromSmiles` calls under the hood): it over-extends aromaticity from a
genuinely aromatic benzo ring, across a bridgehead-N ring fusion, into an
adjacent ring that should stay non-aromatic. This is the same code region
(`crates/chematic-perception/src/aromaticity.rs`, the `aromatic_context`
Pass-2 propagation + "bridgehead N" special case) already named in the
`test_azulene_kekulized_aromatic` / `test_purine_aromatic` #[ignore] comments
as a known-incomplete mechanism -- those are false-negative (under-aromatize)
instances of the same family; this is a false-positive (over-aromatize) one,
not previously documented.

Run: .venv/bin/python3 scripts/smarts_a0_junction_diagnosis.py
"""

import json
import sys
from pathlib import Path

import chematic

REPO = Path(__file__).resolve().parent.parent
DIFF_JSONL = REPO / "validation" / "results" / "rdkit_compat_diff.jsonl"
OUT_JSONL = REPO / "validation" / "smarts_a0_junction_diagnosis.jsonl"


def load_residual_smiles():
    smis = set()
    with open(DIFF_JSONL) as f:
        for line in f:
            row = json.loads(line)
            if row.get("case") == "smarts" and row.get("status") == "count_differs":
                smis.add(row["smiles"])
    return sorted(smis)


def aromatic_atoms(mol):
    return {m[0] for m in chematic.smarts_find("[a]", mol)}


def diagnose(smiles):
    mol = chematic.from_smiles(smiles)
    atoms = mol.atom_table  # [(symbol, atomic_num, charge, aromatic, ?, degree, in_ring), ...]
    n = len(atoms)
    rings = mol.ring_membership()  # list[list[int]] of ring indices per atom

    raw = aromatic_atoms(mol)
    post_mol = mol.apply_aromaticity("rdkit_like")
    post = aromatic_atoms(post_mol)
    extra = sorted(post - raw)

    elements = [atoms[i][0] for i in extra]
    ring_ids = [rings[i] for i in extra]
    bridgeheads = sorted(i for i in range(n) if len(rings[i]) >= 2)

    cc_before = chematic.smarts_find("[#6]=[#6]", mol)
    cc_after = chematic.smarts_find("[#6]=[#6]", post_mol)

    return {
        "smiles": smiles,
        "atom_count": n,
        "raw_aromatic_atoms": sorted(raw),
        "post_apply_aromaticity_atoms": sorted(post),
        "extra_atoms": extra,
        "extra_atom_elements": elements,
        "extra_atom_ring_ids": ring_ids,
        "bridgehead_atoms": bridgeheads,
        "cc_double_bond_matches_before": len(cc_before),
        "cc_double_bond_matches_after": len(cc_after),
        "raw_matches_rdkit_reference": True,  # established separately, see docs/rdkit_compat.md SMARTS-A0
        "mechanism": "apply_aromaticity(rdkit_like) bridgehead-N Pass-2 over-extension into the fused ring adjacent to the genuine aromatic ring",
    }


def main():
    smis = load_residual_smiles()
    print(f"unique molecules in SMARTS-A0 residual: {len(smis)}")

    results = []
    uniform = True
    for smi in smis:
        d = diagnose(smi)
        results.append(d)
        if len(d["extra_atoms"]) != 4 or sorted(d["extra_atom_elements"]) != ["C", "C", "C", "N"]:
            uniform = False
            print(f"NON-UNIFORM CASE: {smi} -> extra={d['extra_atoms']} elements={d['extra_atom_elements']}")

    print(f"all 28 show exactly 4 extra atoms (3C+1N): {uniform}")
    cc_deltas = {d["cc_double_bond_matches_after"] - d["cc_double_bond_matches_before"] for d in results}
    print(
        "all 28 lose exactly 1 `[#6]=[#6]` match after apply_aromaticity "
        f"(the ring-fusion C=C bond gets re-typed Double -> Aromatic order once "
        f"both endpoints are wrongly marked aromatic): deltas={cc_deltas} "
        f"(uniform -1: {cc_deltas == {-1}})"
    )

    OUT_JSONL.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_JSONL, "w") as f:
        for d in results:
            f.write(json.dumps(d) + "\n")
    print(f"wrote {len(results)} rows to {OUT_JSONL}")

    unexplained = [d for d in results if len(d["extra_atoms"]) != 4]
    print(f"unexplained (do not fit the 4-atom pattern): {len(unexplained)}")
    return 0 if not unexplained else 1


if __name__ == "__main__":
    sys.exit(main())
