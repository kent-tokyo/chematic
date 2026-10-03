#!/usr/bin/env python3
"""Compare six source-only SMARTS ring models with pinned RDKit ring sets.

Consumes JSONL from `cargo run -p chematic-smarts --example ring_refusal_probe`.
This is a diagnostic; it never turns a typed refusal into an accepted match.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import sys
from pathlib import Path

from rdkit import Chem, rdBase

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"
PROFILE = ROOT / "validation/results/v1.0.31-source-smarts-310k-profiles.json"
CORPUS_SHA256 = "f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f"
ROWS = [9, 23, 28, 29, 30, 34]


def ring_sets(rings: object) -> set[tuple[int, ...]]:
    return {tuple(sorted(ring)) for ring in rings}


def bond_sets(mol: Chem.Mol, remap: list[int] | None = None) -> set[tuple[tuple[int, int], ...]]:
    result = set()
    for ring in mol.GetRingInfo().BondRings():
        pairs = []
        for bond_index in ring:
            bond = mol.GetBondWithIdx(bond_index)
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if remap is not None:
                a, b = remap[a], remap[b]
            pairs.append(tuple(sorted((a, b))))
        result.add(tuple(sorted(pairs)))
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument("--profile", type=Path, default=PROFILE)
    parser.add_argument("--permutations", type=int, default=16)
    args = parser.parse_args()
    if args.permutations < 0:
        parser.error("permutations must be nonnegative")
    if rdBase.rdkitVersion != "2026.03.6":
        parser.error(f"RDKit {rdBase.rdkitVersion} is not the pinned 2026.03.6 oracle")
    raw = args.corpus.read_bytes()
    if hashlib.sha256(raw).hexdigest() != CORPUS_SHA256:
        parser.error("corpus bytes differ from the pinned 10,000-molecule input")
    smiles = raw.decode("utf-8").splitlines()
    profile = json.loads(args.profile.read_text(encoding="utf-8"))
    parity = profile["profiles"]["parity"]
    if (parity["refused_rows"] != ROWS or parity["typed_refusal_cells"] != 18
            or parity["refused_by_query"] != {"[R1]": 6, "[R2]": 6, "[R3]": 6}):
        parser.error("source profile no longer identifies exactly the pinned refusals")
    records = [json.loads(line) for line in sys.stdin if line.strip()]
    if [record["input_index"] for record in records] != ROWS:
        parser.error("Rust probe row IDs or order differ from the pinned six")

    for record in records:
        index = record["input_index"]
        if record["smiles"] != smiles[index]:
            raise ValueError(f"row {index}: Rust and pinned corpus SMILES disagree")
        rd = Chem.MolFromSmiles(smiles[index])
        if rd is None:
            raise ValueError(f"row {index}: RDKit could not parse pinned SMILES")
        diag = record["diagnostics"]
        if (diag["atom_count"] != rd.GetNumAtoms()
                or diag["elements"] != [atom.GetSymbol() for atom in rd.GetAtoms()]
                or diag["charges"] != [atom.GetFormalCharge() for atom in rd.GetAtoms()]
                or diag["symmetrized_status"] != "Complete"):
            raise ValueError(f"row {index}: atom mapping or bounded ring search differs")
        oracle = ring_sets(rd.GetRingInfo().AtomRings())
        base = ring_sets(diag["base_ring_atom_sets"])
        sym = ring_sets(diag["symmetrized_ring_atom_sets"])
        oracle_bonds = bond_sets(rd)
        base_bonds = {tuple(map(tuple, ring)) for ring in diag["base_ring_bond_sets"]}
        sym_bonds = {tuple(map(tuple, ring)) for ring in diag["symmetrized_ring_bond_sets"]}
        if (len(base_bonds) != len(base) or len(sym_bonds) != len(sym)
                or len(oracle_bonds) != len(oracle)):
            raise ValueError(f"row {index}: duplicate or incomplete ring edge set")
        counts = [rd.GetRingInfo().NumAtomRings(i) for i in range(rd.GetNumAtoms())]
        approximate_mismatches = [i for i, (actual, expected) in enumerate(
            zip(diag["approximate_ring_counts"], counts, strict=True)) if actual != expected]
        shared_mismatches = [i for i, (actual, expected) in enumerate(
            zip(diag["shared_ring_counts"], counts, strict=True)) if actual != expected]
        rng = random.Random(20261003 + index)
        ring_families = {frozenset(oracle)}
        bond_families = {frozenset(oracle_bonds)}
        canonical = Chem.MolToSmiles(rd, isomericSmiles=True)
        for _ in range(args.permutations):
            order = list(range(rd.GetNumAtoms()))
            rng.shuffle(order)
            renumbered = Chem.RenumberAtoms(rd, order)
            if Chem.MolToSmiles(renumbered, isomericSmiles=True) != canonical:
                raise ValueError(f"row {index}: atom permutation changed the molecule")
            Chem.GetSymmSSSR(renumbered)
            remapped = ring_sets([[order[atom] for atom in ring]
                                   for ring in renumbered.GetRingInfo().AtomRings()])
            ring_families.add(frozenset(remapped))
            bond_families.add(frozenset(bond_sets(renumbered, order)))
        def membership(family: frozenset[tuple[int, ...]]) -> tuple[int, ...]:
            return tuple(sum(atom in ring for ring in family) for atom in range(rd.GetNumAtoms()))

        count_families = {membership(family) for family in ring_families}
        result = {
            "input_index": index,
            "rdkit_ring_count": len(oracle),
            "base_ring_count": len(base),
            "symmetrized_ring_count": len(sym),
            "approximate_extra_ring_count": diag["approximate_extra_ring_count"],
            "shared_extra_ring_count": diag["shared_extra_ring_count"],
            "base_missing_ring_sets": sorted(oracle - base),
            "base_extra_ring_sets": sorted(base - oracle),
            "sym_missing_ring_sets": sorted(oracle - sym),
            "sym_extra_ring_sets": sorted(sym - oracle),
            "sym_missing_bond_sets": sorted(oracle_bonds - sym_bonds),
            "sym_extra_bond_sets": sorted(sym_bonds - oracle_bonds),
            "approximate_count_mismatch_atoms": approximate_mismatches,
            "shared_count_mismatch_atoms": shared_mismatches,
            "rdkit_distinct_ring_families_under_permutation": len(ring_families),
            "rdkit_distinct_ring_bond_families_under_permutation": len(bond_families),
            "rdkit_distinct_ring_count_vectors_under_permutation": len(count_families),
            "source_sym_is_rdkit_permutation_family": frozenset(sym) in ring_families,
            "source_sym_bonds_are_rdkit_permutation_family":
                frozenset(sym_bonds) in bond_families,
            "source_approximate_counts_in_rdkit_permutations":
                tuple(diag["approximate_ring_counts"]) in count_families,
            "source_shared_counts_in_rdkit_permutations":
                tuple(diag["shared_ring_counts"]) in count_families,
        }
        print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
