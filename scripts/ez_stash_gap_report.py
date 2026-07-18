#!/usr/bin/env python3
"""
EZ-S1 full-corpus before/after accounting for `chematic_chem::assign_cip`
(legacy engine) E/Z output.

Two fixes landed together in `crates/chematic-chem/src/cip.rs`:

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
   as its geometric complement. This affects every trisubstituted alkene
   end where the marked substituent isn't the CIP-highest one, aromatic or
   not; fixing it flips some *already-assigned* E/Z labels.

Consumes two TSV snapshots produced by
`cargo run -p chematic-chem --release --example ez_stash_gap_snapshot --
<SMILES.csv> <out.tsv>` -- one row per assignment
(`smiles\\tatom_idx\\tpartner_idx\\tcode`), `partner_idx` empty for R/S/r/s
rows, the double-bond's other atom for E/Z rows (used to key against
RDKit's bond-level `_CIPCode`).

Usage:
    .venv/bin/python scripts/ez_stash_gap_report.py baseline.tsv candidate.tsv [SMILES.csv]
"""

import sys

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler


def load_snapshot(path):
    rows = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            smi, idx, partner, code = parts[0], int(parts[1]), parts[2], parts[3]
            rows[(smi, idx)] = (partner, code)
    return rows


def main():
    args = sys.argv[1:]
    if len(args) < 2:
        print(__doc__)
        sys.exit(1)
    baseline_path, candidate_path = args[0], args[1]

    baseline = load_snapshot(baseline_path)
    candidate = load_snapshot(candidate_path)
    all_keys = set(baseline) | set(candidate)

    newly_assigned = [k for k in all_keys if k not in baseline and k in candidate]
    lost = [k for k in all_keys if k in baseline and k not in candidate]
    changed = [
        k for k in all_keys if k in baseline and k in candidate and baseline[k] != candidate[k]
    ]
    rs_codes = {"R", "S", "r", "s"}
    rs_changed = [
        k for k in changed if candidate[k][1] in rs_codes or baseline[k][1] in rs_codes
    ]
    ez_flipped = [
        k for k in changed if baseline[k][1] in ("E", "Z") and candidate[k][1] in ("E", "Z")
    ]
    allene_flipped = [k for k in ez_flipped if "=C=" in k[0]]
    allene_new = [k for k in newly_assigned if "=C=" in k[0] and candidate[k][1] in ("E", "Z")]

    oracle_cache = {}

    def oracle_for(smi):
        if smi not in oracle_cache:
            rd = Chem.MolFromSmiles(smi)
            bond_codes = {}
            if rd is not None:
                try:
                    rdCIPLabeler.AssignCIPLabels(rd)
                    for b in rd.GetBonds():
                        if b.HasProp("_CIPCode"):
                            bond_codes[frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))] = (
                                b.GetProp("_CIPCode")
                            )
                except Exception:
                    pass
            oracle_cache[smi] = bond_codes
        return oracle_cache[smi]

    def oracle_code_for(key, partner):
        smi, idx = key
        if not partner:
            return None
        return oracle_for(smi).get(frozenset((idx, int(partner))))

    # Every E/Z row present on either side, checked against a freshly
    # regenerated RDKit oracle -- this is the authoritative accuracy number
    # (newly-assigned agreement is a subset of it).
    ez_keys = [
        k
        for k in all_keys
        if (k in baseline and baseline[k][1] in ("E", "Z"))
        or (k in candidate and candidate[k][1] in ("E", "Z"))
    ]
    baseline_correct = candidate_correct = newly_correct = both_wrong = oracle_missing = 0
    regressions = []
    for k in ez_keys:
        partner = (candidate.get(k) or baseline.get(k))[0]
        oracle_code = oracle_code_for(k, partner)
        if oracle_code is None:
            oracle_missing += 1
            continue
        b_ok = k in baseline and baseline[k][1] == oracle_code
        c_ok = k in candidate and candidate[k][1] == oracle_code
        baseline_correct += b_ok
        candidate_correct += c_ok
        newly_correct += c_ok and not b_ok
        both_wrong += not b_ok and not c_ok
        if b_ok and not c_ok:
            regressions.append((k, baseline[k], candidate[k], oracle_code))

    new_agree = new_disagree = 0
    new_disagree_examples = []
    for k in newly_assigned:
        partner, code = candidate[k]
        if code not in ("E", "Z"):
            continue
        oracle_code = oracle_code_for(k, partner)
        if oracle_code is None:
            continue
        if oracle_code == code:
            new_agree += 1
        else:
            new_disagree += 1
            new_disagree_examples.append((*k, code, oracle_code))

    print("=== EZ-S1 full-corpus before/after ===")
    print(f"E/Z newly assigned:        {len(newly_assigned)}")
    print(f"E/Z lost:                  {len(lost)} (gate: must be 0)")
    print(f"R/S changed:               {len(rs_changed)} (gate: must be 0)")
    print(f"E<->Z flipped (existing):  {len(ez_flipped)}")
    print(f"  of which allene-tagged (heuristic '=C=' substring): {len(allene_flipped)}")
    print(f"allene newly assigned (heuristic): {len(allene_new)}")
    print()
    print(f"newly-assigned E/Z vs RDKit: agree={new_agree} disagree={new_disagree}")
    for ex in new_disagree_examples[:20]:
        print("  DISAGREE", ex)
    print()
    print(
        f"all E/Z rows vs RDKit: baseline_correct={baseline_correct} "
        f"candidate_correct={candidate_correct} newly_correct={newly_correct} "
        f"both_wrong={both_wrong} oracle_missing={oracle_missing}"
    )
    print(f"regressions (baseline correct -> candidate incorrect): {len(regressions)} (gate: must be 0)")
    for r in regressions[:20]:
        print("  REGRESSION", r)


if __name__ == "__main__":
    main()
