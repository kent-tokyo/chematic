#!/usr/bin/env python3
"""
Phase 1/2 of the ECFP4/canonical_smiles/InChI shared-residual mechanism
discrimination (see plan: identify why ~13% of molecules still diverge
across all three after apply_aromaticity(), even though per-atom and
per-bond order-independent MULTISETS are identical).

The multiset check used in scripts/ecfp4_agreement.py (tier 5) has a known
blind spot: it cannot detect a symmetric swap where two atoms of identical
(atomic_number, degree) trade `aromatic` flags between the aromatic-origin
and Kekule-origin-then-perceived representations -- same aggregate profile,
different physical assignment. This script builds a REAL per-atom/per-bond
correspondence between the two representations (not index-based, which is
meaningless since the Kekule respelling has a different atom traversal
order; not multiset-based, which has the blind spot above) using RDKit's
canonical atom ranking as an independent oracle: since both RDKit parses
represent the identical molecule, `Chem.CanonicalRankAtoms(rd, breakTies=
True)` assigns the SAME rank to the SAME physical atom in both parses.
Chematic parses each SMILES string in the same atom order RDKit does for
that string (verified directly, not assumed -- see the module docstring's
sibling verification in ecfp4_agreement.py's methodology), so chematic
atom index == RDKit atom index for a *given* string, giving a full
physical-atom correspondence path: rank -> RDKit index (per string) ->
chematic index (per string).

Source-grounded context (see the plan file / memory for full citations):
apply_aromaticity_ex() (crates/chematic-perception/src/aromaticity.rs)
normalizes every aromaticity-model bond to BondOrder::Aromatic uniformly --
so a "leftover Kekule single/double pattern" hypothesis is already refuted
by source reading. What's left to test here is whether the aromaticity
MODEL itself (AromaticityModel::aromatic_atoms/aromatic_bonds) assigns a
genuinely different SET of atoms/bonds as aromatic between the two
representations -- which this script checks directly, per physical atom
and bond, for the first time.

Usage:
    .venv/bin/python scripts/aromaticity_mechanism_probe.py [SMILES.csv]
        --molecules smi1,smi2,...       # explicit SMILES to probe
        --positive-control              # run against a known-differing
                                         # (41-set) example first
"""
import argparse
import os
import sys


def rdkit_correspondence(rd_a, rd_b):
    """Map RDKit atom index in rd_a -> RDKit atom index in rd_b, for the
    same physical atom, via canonical rank. Returns None if the two mols
    don't have a valid rank bijection (shouldn't happen for the same
    molecule, but don't assume)."""
    from rdkit import Chem

    ranks_a = list(Chem.CanonicalRankAtoms(rd_a, breakTies=True))
    ranks_b = list(Chem.CanonicalRankAtoms(rd_b, breakTies=True))
    if set(ranks_a) != set(ranks_b) or len(set(ranks_a)) != rd_a.GetNumAtoms():
        return None
    rank_to_idx_b = {r: i for i, r in enumerate(ranks_b)}
    return {i: rank_to_idx_b[r] for i, r in enumerate(ranks_a)}


def probe_pair(smi, kek_smi, chematic, Chem):
    """Compare cm_arom (parsed from `smi`) against cm_kek_perceived (parsed
    from `kek_smi`, then apply_aromaticity()'d), per PHYSICAL atom/bond via
    RDKit canonical-rank correspondence. Returns a dict of findings, or
    None if any precondition (parse, correspondence, atom-order match)
    fails."""
    rd_a = Chem.MolFromSmiles(smi)
    rd_b = Chem.MolFromSmiles(kek_smi)
    if rd_a is None or rd_b is None:
        return None

    corr = rdkit_correspondence(rd_a, rd_b)
    if corr is None:
        return None

    try:
        cm_a = chematic.from_smiles(smi)
        cm_b = chematic.from_smiles(kek_smi).apply_aromaticity()
    except Exception:
        return None

    at_a, at_b = cm_a.atom_table, cm_b.atom_table
    if len(at_a) != rd_a.GetNumAtoms() or len(at_b) != rd_b.GetNumAtoms():
        return None
    # Verify the same-order assumption directly rather than trusting it --
    # chematic atom i's element must match RDKit atom i's element for BOTH
    # strings, or the index correspondence built above is meaningless.
    if [row[0] for row in at_a] != [a.GetSymbol() for a in rd_a.GetAtoms()]:
        return None
    if [row[0] for row in at_b] != [a.GetSymbol() for a in rd_b.GetAtoms()]:
        return None

    atom_flag_diffs = []
    for ia, ib in corr.items():
        arom_a = at_a[ia][3]
        arom_b = at_b[ib][3]
        if arom_a != arom_b:
            atom_flag_diffs.append(
                {"rdkit_idx_a": ia, "rdkit_idx_b": ib, "element": at_a[ia][0],
                 "aromatic_a": arom_a, "aromatic_b": arom_b}
            )

    # Build a (a1,a2) -> (btype, aromatic) map for cm_a's bonds, keyed by
    # RDKit-index pairs (chematic index == rdkit index per string).
    bonds_a = {}
    for a1, a2, btype, arom in cm_a.bond_table:
        bonds_a[frozenset((a1, a2))] = (btype, arom)
    bonds_b = {}
    for a1, a2, btype, arom in cm_b.bond_table:
        bonds_b[frozenset((a1, a2))] = (btype, arom)

    bond_diffs = []
    for key_a, (btype_a, arom_a) in bonds_a.items():
        ia1, ia2 = tuple(key_a)
        key_b = frozenset((corr[ia1], corr[ia2]))
        if key_b not in bonds_b:
            bond_diffs.append({"rdkit_pair_a": sorted(key_a), "issue": "no corresponding bond in b"})
            continue
        btype_b, arom_b = bonds_b[key_b]
        if btype_a != btype_b or arom_a != arom_b:
            bond_diffs.append(
                {"rdkit_pair_a": sorted(key_a), "btype_a": btype_a, "btype_b": btype_b,
                 "aromatic_a": arom_a, "aromatic_b": arom_b}
            )

    ecfp4_differs = cm_a.ecfp4() != cm_b.ecfp4()
    csmi_differs = cm_a.canonical_smiles_mode("normal") != cm_b.canonical_smiles_mode("normal")

    return {
        "smiles": smi,
        "kekule_smiles": kek_smi,
        "n_atom_flag_diffs": len(atom_flag_diffs),
        "atom_flag_diffs": atom_flag_diffs[:5],
        "n_bond_diffs": len(bond_diffs),
        "bond_diffs": bond_diffs[:5],
        "ecfp4_differs": ecfp4_differs,
        "canonical_smiles_differs": csmi_differs,
    }


def build_kek_smi(smi, Chem):
    rd_kek = Chem.MolFromSmiles(smi)
    if rd_kek is None:
        return None
    try:
        Chem.Kekulize(rd_kek, clearAromaticFlags=True)
    except Exception:
        return None
    return Chem.MolToSmiles(rd_kek, kekuleSmiles=True, canonical=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("smiles_csv", nargs="?", default="~/Downloads/SMILES.csv")
    parser.add_argument("--molecules", default=None, help="Comma-separated explicit SMILES to probe")
    parser.add_argument("--sample-n", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
    except ImportError:
        sys.exit("rdkit not installed. pip install rdkit")
    import chematic
    import random

    if args.molecules:
        smis = [s.strip() for s in args.molecules.split(",") if s.strip()]
        for smi in smis:
            kek_smi = build_kek_smi(smi, Chem)
            if kek_smi is None:
                print(f"{smi}: could not build Kekule respelling, skipped")
                continue
            result = probe_pair(smi, kek_smi, chematic, Chem)
            if result is None:
                print(f"{smi}: probe preconditions failed (parse/correspondence/order-check), skipped")
                continue
            print(f"SMILES: {smi}")
            print(f"  atom aromatic-flag diffs (true correspondence): {result['n_atom_flag_diffs']}")
            for d in result["atom_flag_diffs"]:
                print(f"    {d}")
            print(f"  bond diffs (true correspondence): {result['n_bond_diffs']}")
            for d in result["bond_diffs"]:
                print(f"    {d}")
            print(f"  ecfp4 differs: {result['ecfp4_differs']}, canonical_smiles differs: {result['canonical_smiles_differs']}")
            print()
        return

    path = os.path.expanduser(args.smiles_csv)
    with open(path) as f:
        lines = [l.strip() for l in f if l.strip()]
    smis = [l.split(",")[0].strip() for l in lines if l.split(",")[0].strip().lower() != "smiles"]

    rng = random.Random(args.seed)
    sample = smis if len(smis) <= args.sample_n else rng.sample(smis, args.sample_n)

    # Reproduce tier 5's split: which molecules have identical assignment
    # multisets (the "89" set) vs disagreeing ones (the "41" set)?
    from collections import Counter

    def atom_multiset(mol):
        return Counter((row[1], row[3], row[5]) for row in mol.atom_table)

    def bond_multiset(mol):
        at = mol.atom_table
        ms = Counter()
        for a1, a2, btype, _arom in mol.bond_table:
            e1, e2 = at[a1][1], at[a2][1]
            ms[(min(e1, e2), max(e1, e2), btype)] += 1
        return ms

    set_89 = []  # multiset-identical, ecfp4 still differs
    set_41 = []  # multiset-differs
    n_checked = 0
    for smi in sample:
        kek_smi = build_kek_smi(smi, Chem)
        if kek_smi is None:
            continue
        try:
            cm_a = chematic.from_smiles(smi)
            cm_b = chematic.from_smiles(kek_smi).apply_aromaticity()
        except Exception:
            continue
        n_checked += 1
        if cm_a.ecfp4() == cm_b.ecfp4():
            continue
        same_assignment = (
            atom_multiset(cm_a) == atom_multiset(cm_b)
            and bond_multiset(cm_a) == bond_multiset(cm_b)
        )
        if same_assignment:
            set_89.append((smi, kek_smi))
        else:
            set_41.append((smi, kek_smi))

    print(f"n_checked={n_checked}, |89-set|={len(set_89)}, |41-set|={len(set_41)}")
    print()

    print("=== POSITIVE CONTROL: probe against the 41-set (known to differ) ===")
    n_control_detected = 0
    for smi, kek_smi in set_41[:10]:
        result = probe_pair(smi, kek_smi, chematic, Chem)
        if result is None:
            print(f"  {smi}: probe preconditions failed, skipped")
            continue
        detected = result["n_atom_flag_diffs"] > 0 or result["n_bond_diffs"] > 0
        if detected:
            n_control_detected += 1
        print(f"  {smi[:60]}: atom_diffs={result['n_atom_flag_diffs']}, "
              f"bond_diffs={result['n_bond_diffs']}, detected={detected}")
    print(f"  Positive control: {n_control_detected}/{min(10, len(set_41))} correctly detected a difference")
    print()

    if n_control_detected == 0 and set_41:
        print("POSITIVE CONTROL FAILED -- the correspondence tool did not detect ANY "
              "difference on molecules already known to differ. Do not trust results "
              "on the 89-set below. Stopping.")
        return

    print("=== Probing the 89-set (multiset-identical, ecfp4 still differs) ===")
    n_atom_diff = 0
    n_bond_diff = 0
    n_neither = 0
    examples_neither = []
    for smi, kek_smi in set_89:
        result = probe_pair(smi, kek_smi, chematic, Chem)
        if result is None:
            continue
        if result["n_atom_flag_diffs"] > 0:
            n_atom_diff += 1
        elif result["n_bond_diffs"] > 0:
            n_bond_diff += 1
        else:
            n_neither += 1
            if len(examples_neither) < 5:
                examples_neither.append(result["smiles"])

    print(f"  |89-set| = {len(set_89)}")
    print(f"  true per-atom aromatic-flag correspondence differs: {n_atom_diff}")
    print(f"  true per-bond correspondence differs (flags matched): {n_bond_diff}")
    print(f"  neither differs under true correspondence (mystery deepens): {n_neither}")
    if examples_neither:
        print(f"  examples where neither differs: {examples_neither}")


if __name__ == "__main__":
    main()
