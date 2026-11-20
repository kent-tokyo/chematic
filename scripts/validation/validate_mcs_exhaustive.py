#!/usr/bin/env python3
"""Independently enumerate the MCS of small, nonaromatic development graphs."""

import argparse
from collections import Counter
from itertools import combinations, product
from pathlib import Path
import json

import chematic
from rdkit import Chem, rdBase
from rdkit.Chem import rdFMCS

FIXTURES = (
    "C", "CC", "CCC", "CCCCC", "CC(C)C", "C=C", "CC=O", "COC",
    "C1CC1", "C1CCC1", "C1CCOC1", "C1CC2CC12", "C12C3C1C23.CCCCC",
    "[13CH3]CN", "[12CH3]CN", "[13CH3][NH3+]", "[NH3+]C", "NCC",
    "[*]O", "[*]C",
)
MULTI = (
    ("C1CC1", "CCC", "C1CCC1"),
    ("C12C3C1C23.CCCCC", "C12C3C1C23.CCCCC", "C12C3C1C23.CCCCC"),
    ("CC=O", "CC=N", "CCC"),
    ("[13CH3][NH3+]", "[NH3+][13CH3]", "[NH3+][12CH3]"),
    ("C1CCOC1", "COC", "CCOCC"),
)


def graph(molecule, charge, isotope):
    atoms = [
        (a.GetAtomicNum(), a.GetFormalCharge() if charge else 0,
         a.GetIsotope() if isotope else 0)
        for a in molecule.GetAtoms()
    ]
    edges = {(min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
              max(b.GetBeginAtomIdx(), b.GetEndAtomIdx())): b.GetBondType()
             for b in molecule.GetBonds()}
    return atoms, edges


def induced_vertices(edges):
    """Return the vertices of a connected edge set, or an empty set."""
    vertices = set().union(*map(set, edges))
    reached = {next(iter(edges))[0]}
    while True:
        expanded = reached.union(*(edge for edge in edges if reached.intersection(edge)))
        if expanded == reached:
            return reached if reached == vertices else set()
        reached = expanded


def embeds(labels, edges, target):
    nodes, bonds = target
    neighbors = {i: [] for i in labels}
    for (a, b), order in edges.items():
        neighbors[a].append((b, order))
        neighbors[b].append((a, order))
    possible = {i: [j for j, atom in enumerate(nodes) if atom == label]
                for i, label in labels.items()}
    order = sorted(labels, key=lambda i: (len(possible[i]), -len(neighbors[i]), i))
    assigned = {}
    occupied = set()

    def visit(depth):
        if depth == len(order):
            return True
        node = order[depth]
        for candidate in possible[node]:
            if candidate in occupied or any(
                bonds.get(tuple(sorted((candidate, assigned[neighbor])))) != kind
                for neighbor, kind in neighbors[node] if neighbor in assigned
            ):
                continue
            assigned[node] = candidate
            occupied.add(candidate)
            if visit(depth + 1):
                return True
            occupied.remove(candidate)
            del assigned[node]
        return False

    return visit(0)


def exhaustive_score(graphs, objective):
    source, *targets = sorted(graphs, key=lambda g: len(g[1]))
    atoms, edges = source
    priority = lambda pair: pair[::-1] if objective == "bonds" else pair
    best = (0, 0)
    for idx, label in enumerate(atoms):
        if all(embeds({idx: label}, {}, target) for target in targets):
            best = (1, 0)
            break
    for size in range(1, len(edges) + 1):
        for subset in combinations(edges, size):
            vertices = induced_vertices(subset)
            score = (len(vertices), size)
            if not vertices or priority(score) <= priority(best):
                continue
            query = {i: atoms[i] for i in vertices}
            if all(embeds(query, {e: edges[e] for e in subset}, target) for target in targets):
                best = score
    return best


def returned_score(smarts, molecules):
    if smarts is None:
        return [0, 0], True
    query = Chem.MolFromSmarts(smarts)
    if query is None:
        return None, False
    return [query.GetNumAtoms(), query.GetNumBonds()], all(
        mol.HasSubstructMatch(query) for mol in molecules
    )


def run_case(inputs, charge, isotope, objective):
    oracle_mols = [Chem.MolFromSmiles(s) for s in inputs]
    expected = exhaustive_score(
        [graph(mol, charge, isotope) for mol in oracle_mols], objective
    )
    candidate_mols = [chematic.from_smiles(s) for s in inputs]
    smarts, timed_out = chematic.find_mcs_smarts_checked(
        candidate_mols, timeout_ms=2000, match_charge=charge, match_isotope=isotope,
        bond_count_first=objective == "bonds",
    )
    actual, rematched = returned_score(smarts, oracle_mols)
    p = rdFMCS.MCSParameters()
    p.AtomCompareParameters.MatchFormalCharge = charge
    p.AtomCompareParameters.MatchIsotope = isotope
    p.MaximizeBonds = True
    p.Timeout = 2
    rdkit = rdFMCS.FindMCS(oracle_mols, p)
    return {
        "inputs": inputs, "match_charge": charge, "match_isotope": isotope,
        "objective": objective,
        "expected_atoms_bonds": expected, "candidate_atoms_bonds": actual,
        "candidate_smarts": smarts, "rematched": rematched, "timed_out": timed_out,
        "passed": tuple(actual or []) == expected and rematched and not timed_out,
        "rdkit_bond_first_atoms_bonds": [rdkit.numAtoms, rdkit.numBonds],
        "rdkit_canceled": rdkit.canceled,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--objective", choices=("atoms", "bonds", "both"), default="both")
    args = parser.parse_args()
    inputs = list(combinations(FIXTURES, 2)) + list(MULTI)
    objectives = ("atoms", "bonds") if args.objective == "both" else (args.objective,)
    rows = [run_case(mols, charge, isotope, objective)
            for items in inputs for mols in (items, items[::-1])
            for charge, isotope in product((False, True), repeat=2)
            for objective in objectives]
    totals = Counter("pass" if row["passed"] else "fail" for row in rows)
    args.output.write_text(json.dumps({
        "cohort": "exposed small graphs with enumerable connected edge subsets",
        "objectives": objectives,
        "rdkit_version": rdBase.rdkitVersion,
        "chematic_import": chematic.__file__,
        "totals": totals, "rows": rows,
    }, indent=2) + "\n")
    print(dict(totals))
    if totals["fail"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
