#!/usr/bin/env python3
"""Issue #227 Priority 1A-3, Phase 1.2: exhaustive sub-bucket evidence for the
37-atom `terminal_oxygen_o2cm_umbrella_gap` residual left after PR #239.

For each atom, mirrors RDKit's real `AtomTyper.cpp` `case 8:` (aliphatic
oxygen), 1-neighbor ("terminal") branch logic directly against the parsed
RDKit molecule (not re-derived heuristically) to determine which named
condition (isCarboxylateO / isNitroO / isNOxideO / isThioSulfinateO /
isSulfateO / isPhosphateOrPerchlorateO) actually fires, then maps that to
one of the exclusive reporting buckets. `unclassified` must be 0.

Usage: .venv/bin/python scripts/mmff94_o2cm_subbucket_evidence_227.py
Reads: validation/results/mmff94_atom_type_bucket_evidence_227.jsonl
Writes: validation/results/mmff94_o2cm_subbucket_evidence_227.jsonl (stdout)
        summary to stderr
"""
import json
import sys

from rdkit import Chem

SRC = "validation/results/mmff94_atom_type_bucket_evidence_227.jsonl"


def total_degree(atom):
    return atom.GetDegree() + atom.GetTotalNumHs()


def total_valence(atom):
    explicit = sum(b.GetBondTypeAsDouble() for b in atom.GetBonds())
    return explicit + atom.GetTotalNumHs()


def is_terminal_o(atom):
    return (
        atom.GetSymbol() == "O"
        and atom.GetDegree() <= 1
        and atom.GetTotalNumHs() == 0
    )


def count_terminal_o_neighbors(central):
    return sum(1 for n in central.GetNeighbors() if is_terminal_o(n))


def count_deg2_n_neighbors(central):
    return sum(
        1
        for n in central.GetNeighbors()
        if n.GetSymbol() == "N" and total_degree(n) == 2
    )


def count_terminal_s_neighbors(central):
    return sum(
        1
        for n in central.GetNeighbors()
        if n.GetSymbol() == "S" and total_degree(n) == 1
    )


def classify(mol, o_idx):
    """Returns (fired_condition, subbucket, detail) mirroring AtomTyper.cpp
    case 8 lines 1554-1737 at the pinned commit."""
    o_atom = mol.GetAtomWithIdx(o_idx)
    nbrs = list(o_atom.GetNeighbors())
    if len(nbrs) != 1:
        return ("no_single_neighbor", "unclassified", f"O has {len(nbrs)} neighbors")
    central = nbrs[0]
    bond = mol.GetBondBetweenAtoms(o_idx, central.GetIdx())
    bond_double = bond.GetBondType() == Chem.BondType.DOUBLE
    bond_single = bond.GetBondType() == Chem.BondType.SINGLE
    elem = central.GetSymbol()

    n_o = count_terminal_o_neighbors(central)
    n_n = count_deg2_n_neighbors(central)
    n_s = count_terminal_s_neighbors(central)

    detail = {
        "central_element": elem,
        "central_idx": central.GetIdx(),
        "central_charge": central.GetFormalCharge(),
        "central_degree": central.GetDegree(),
        "central_total_degree": total_degree(central),
        "central_total_valence": total_valence(central),
        "bond_order": "DOUBLE" if bond_double else ("SINGLE" if bond_single else str(bond.GetBondType())),
        "n_terminal_o_on_central": n_o,
        "n_deg2_n_on_central": n_n,
        "n_terminal_s_on_central": n_s,
    }

    if elem == "C":
        if n_o == 2:
            return ("isCarboxylateO", "carboxylate_terminal_o", detail)
        # bond_double -> isCarbonylO (routes to 7, not part of the O2CM
        # residual by construction -- if we get here for a 37-atom residual
        # row, it's a genuine surprise, not expected).
        return ("none_fired_C", "unclassified", detail)

    if elem == "N":
        if n_o >= 2:
            return ("isNitroO", "nitro_terminal_o", detail)
        if bond_single and n_o == 1:
            deg = total_degree(central)
            val = total_valence(central)
            if val == 4:
                return ("isNOxideO", "n_oxide_terminal_o", detail)
            if deg == 2 or val == 3:
                return ("isOxideOBondedToN", "other_classified", {**detail, "note": "routes to OM(35), not O2CM -- should not appear in this residual"})
        return ("none_fired_N", "unclassified", detail)

    if elem == "S":
        if n_s == 1:
            return ("isThioSulfinateO", "sulfoxide_terminal_o", {**detail, "note": "thiosulfinate, reported under sulfoxide_terminal_o (S=O family)"})
        if bond_single or (bond_double and (n_o + n_n) > 1):
            # Sub-distinguish sulfone / sulfonate / sulfate / sulfonamide
            # for reporting only -- RDKit's own code does not distinguish
            # these, they all share the single isSulfateO condition.
            c_neighbors = sum(1 for n in central.GetNeighbors() if n.GetSymbol() == "C")
            bridging_o = sum(
                1
                for n in central.GetNeighbors()
                if n.GetSymbol() == "O" and not is_terminal_o(n)
            )
            if c_neighbors >= 1 and n_o >= 2:
                return ("isSulfateO", "sulfone_terminal_o", {**detail, "note": f"C-bonded S with {n_o} terminal O -- sulfone/sulfonate/sulfonamide family"})
            if bridging_o >= 1 or c_neighbors == 0:
                return ("isSulfateO", "sulfate_terminal_o", {**detail, "note": "no direct S-C bond -- sulfate family"})
            return ("isSulfateO", "sulfonate_terminal_o", detail)
        if bond_double and (n_o + n_n) == 1:
            return ("isSulfoxideO", "sulfoxide_terminal_o", {**detail, "note": "routes to O=C(7), not O2CM -- should not appear in this residual"})
        return ("none_fired_S", "unclassified", detail)

    if elem == "P":
        c_neighbors = sum(1 for n in central.GetNeighbors() if n.GetSymbol() == "C")
        sub = "phosphonate_terminal_o" if c_neighbors >= 1 else "phosphate_terminal_o"
        return ("isPhosphateOrPerchlorateO", sub, detail)

    if elem == "Cl":
        return ("isPhosphateOrPerchlorateO", "perchlorate_terminal_o", detail)

    return ("no_matching_central_element", "unclassified", detail)


def main():
    rows = [json.loads(l) for l in open(SRC)]
    sub = [r for r in rows if r.get("bucket") == "terminal_oxygen_o2cm_umbrella_gap"]

    out = []
    bucket_counts = {}
    unclassified = 0
    for r in sub:
        mol = Chem.MolFromSmiles(r["smiles"])
        if mol is None:
            out.append({**r, "subbucket": "unclassified", "fired_condition": "smiles_parse_failed"})
            unclassified += 1
            continue
        fired, subbucket, detail = classify(mol, r["atom_index"])
        if subbucket == "unclassified":
            unclassified += 1
        bucket_counts[subbucket] = bucket_counts.get(subbucket, 0) + 1
        out.append({**r, "fired_rdkit_condition": fired, "subbucket": subbucket, "rdkit_side_detail": detail})
        print(json.dumps(out[-1]))

    summary = {
        "total": len(sub),
        "unclassified": unclassified,
        "bucket_counts": bucket_counts,
    }
    print(json.dumps(summary, indent=2), file=sys.stderr)


if __name__ == "__main__":
    main()
