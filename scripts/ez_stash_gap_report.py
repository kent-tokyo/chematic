#!/usr/bin/env python3
"""
EZ-S1 full-corpus before/after accounting for `chematic_chem::assign_cip`
(legacy engine) E/Z output.

Three fixes landed together in `crates/chematic-chem/src/cip.rs`:

1. `substituent_is_up` now reads `Molecule::bond_direction` -- the side
   channel used when a SMILES `/`/`\\` marker lands on a bond between two
   aromatic-flagged atoms (e.g. a ring bond flanking an exocyclic C=N).
   Previously it read only the bond's own `order`, so these bonds produced
   no E/Z at all.
2. That read exposed a second, pre-existing, unrelated bug in
   `highest_stereo_sub`: it returned the highest-priority substituent
   *among those carrying an explicit marker*, not the true highest-priority
   substituent overall -- silently using a lower-priority marked
   substituent's raw side instead of deriving the true-highest one's side
   as its geometric complement. Affects every trisubstituted alkene end
   where the marked substituent isn't the CIP-highest one, aromatic or not.
3. A missing tie guard in `highest_stereo_sub`: when the only two
   substituents at an alkene end are a genuine CIP priority tie (e.g. the
   two ring branches of an *unsubstituted*, symmetric ring's ipso carbon),
   there is no stereogenic bond to report -- swapping them maps the
   molecule onto itself. Without the guard, a stable sort silently picked
   whichever substituent happened to come first in adjacency order, which
   is not stable across atom renumbering (e.g. a `canonical_smiles` round
   trip), flipping the reported side arbitrarily.

This is an ORACLE-FIRST accounting: the RDKit E/Z oracle is built by
scanning every SMILES in the corpus directly (not just the ones chematic
happened to assign something for), so a bond chematic silently never
assigns anything to -- not just one it assigns the wrong code to -- is
counted (`*_missing`). Axial (allene) chirality is excluded from the
oracle: RDKit doesn't label individual bonds of a cumulated diene the same
way, and it isn't what this accounting claims to measure.

Consumes two TSV snapshots produced by
`cargo run -p chematic-chem --release --example ez_stash_gap_snapshot --
<SMILES.csv> <out.tsv>` -- one row per assignment
(`smiles\\tkind\\tatom_idx\\tpartner_idx\\tcode`), `kind` is `tetra`/`ez`/
`allene`/`parse_fail`; `partner_idx` is the double bond's other atom for
`ez`/`allene` rows (used to key against RDKit's bond-level `_CIPCode`).

Usage:
    .venv/bin/python scripts/ez_stash_gap_report.py baseline.tsv candidate.tsv SMILES.csv
"""

import sys

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler


def load_snapshot(path):
    """-> ({(smiles, kind, atom_idx): (partner_idx, code)}, {parse_fail smiles})"""
    rows = {}
    parse_fails = set()
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            smi, kind = parts[0], parts[1]
            if kind == "parse_fail":
                parse_fails.add(smi)
                continue
            atom_idx, partner, code = int(parts[2]), parts[3], parts[4]
            rows[(smi, kind, atom_idx)] = (partner, code)
    return rows, parse_fails


def ez_index(rows):
    """{(kind='ez' only)} -> {smiles: {frozenset(bond atoms): code}}"""
    idx = {}
    for (smi, kind, atom_idx), (partner, code) in rows.items():
        if kind != "ez" or not partner:
            continue
        idx.setdefault(smi, {})[frozenset((atom_idx, int(partner)))] = code
    return idx


def diff_by_kind(baseline, candidate, kind):
    b = {(smi, idx): code for (smi, k, idx), (_, code) in baseline.items() if k == kind}
    c = {(smi, idx): code for (smi, k, idx), (_, code) in candidate.items() if k == kind}
    keys = set(b) | set(c)
    newly = [k for k in keys if k not in b and k in c]
    lost = [k for k in keys if k in b and k not in c]
    flipped = [k for k in keys if k in b and k in c and b[k] != c[k]]
    return newly, lost, flipped


def is_allene_central(atom):
    dbl = sum(1 for b in atom.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE)
    return dbl == 2 and atom.GetDegree() == 2


def rdkit_ez_oracle(smi):
    """{frozenset(bond atom idxs): code} for smi's plain (non-allene) E/Z
    bonds, or None if RDKit can't parse/label it."""
    rd = Chem.MolFromSmiles(smi)
    if rd is None:
        return None
    try:
        rdCIPLabeler.AssignCIPLabels(rd)
    except Exception:
        return None
    out = {}
    for b in rd.GetBonds():
        if not b.HasProp("_CIPCode"):
            continue
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        if is_allene_central(a1) or is_allene_central(a2):
            continue
        out[frozenset((a1.GetIdx(), a2.GetIdx()))] = b.GetProp("_CIPCode")
    return out


def main():
    args = sys.argv[1:]
    if len(args) < 3:
        print(__doc__)
        sys.exit(1)
    baseline_path, candidate_path, csv_path = args[0], args[1], args[2]

    baseline, baseline_parse_fail = load_snapshot(baseline_path)
    candidate, candidate_parse_fail = load_snapshot(candidate_path)
    base_idx = ez_index(baseline)
    cand_idx = ez_index(candidate)

    with open(csv_path) as f:
        smis = [line.strip() for line in f if line.strip()]

    rdkit_parse_fail = 0
    rdkit_ez_total = 0
    cand_correct = cand_wrong = cand_missing = 0
    base_correct = base_wrong = base_missing = 0
    cand_extra = 0
    base_extra = 0
    wrong_examples = []
    missing_examples = []
    extra_examples = []

    for smi in smis:
        oracle = rdkit_ez_oracle(smi)
        if oracle is None:
            rdkit_parse_fail += 1
            continue
        rdkit_ez_total += len(oracle)

        if smi not in candidate_parse_fail:
            c_map = cand_idx.get(smi, {})
            for key, ocode in oracle.items():
                c_code = c_map.get(key)
                if c_code is None:
                    cand_missing += 1
                    missing_examples.append((smi, sorted(key), ocode))
                elif c_code == ocode:
                    cand_correct += 1
                else:
                    cand_wrong += 1
                    wrong_examples.append((smi, sorted(key), c_code, ocode))
            for key in c_map:
                if key not in oracle:
                    cand_extra += 1
                    extra_examples.append((smi, sorted(key), c_map[key]))

        if smi not in baseline_parse_fail:
            b_map = base_idx.get(smi, {})
            for key, ocode in oracle.items():
                b_code = b_map.get(key)
                if b_code is None:
                    base_missing += 1
                elif b_code == ocode:
                    base_correct += 1
                else:
                    base_wrong += 1
            for key in b_map:
                if key not in oracle:
                    base_extra += 1

    ez_newly, ez_lost, ez_flipped = diff_by_kind(baseline, candidate, "ez")
    allene_newly, allene_lost, allene_flipped = diff_by_kind(baseline, candidate, "allene")
    tetra_newly, tetra_lost, tetra_flipped = diff_by_kind(baseline, candidate, "tetra")

    print("=== EZ-S1 oracle-first E/Z completeness (RDKit rdCIPLabeler) ===")
    print(f"corpus size: {len(smis)}")
    print(f"RDKit parse failures: {rdkit_parse_fail}")
    print(f"chematic parse failures (baseline): {len(baseline_parse_fail)}")
    print(f"chematic parse failures (candidate): {len(candidate_parse_fail)}")
    print(f"rdkit_ez_total: {rdkit_ez_total}")
    print()
    print(f"candidate_correct: {cand_correct} (gate: == rdkit_ez_total)")
    print(f"candidate_wrong:   {cand_wrong} (gate: == 0)")
    print(f"candidate_missing: {cand_missing} (gate: == 0)")
    print(f"candidate_extra:   {cand_extra} (gate: == 0)")
    print()
    print(f"baseline_correct: {base_correct}")
    print(f"baseline_wrong:   {base_wrong}")
    print(f"baseline_missing: {base_missing}")
    print(f"baseline_extra:   {base_extra}")
    print()
    for ex in wrong_examples[:20]:
        print("  WRONG", ex)
    for ex in missing_examples[:20]:
        print("  MISSING", ex)
    for ex in extra_examples[:20]:
        print("  EXTRA", ex)

    print()
    print("=== kind-separated before/after diff (all rows, not just oracle-covered) ===")
    print(f"ez newly assigned: {len(ez_newly)}")
    print(f"ez lost:           {len(ez_lost)} (gate: == 0)")
    print(f"ez flipped:        {len(ez_flipped)}")
    print(f"allene newly assigned: {len(allene_newly)} (gate: == 0)")
    print(f"allene lost:           {len(allene_lost)} (gate: == 0)")
    print(f"allene flipped:        {len(allene_flipped)} (gate: == 0)")
    print(f"R/S (tetra) newly assigned: {len(tetra_newly)}")
    print(f"R/S (tetra) lost:           {len(tetra_lost)}")
    print(f"R/S (tetra) changed:        {len(tetra_flipped)} (gate: == 0)")


if __name__ == "__main__":
    main()
