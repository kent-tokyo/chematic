#!/usr/bin/env python3
"""
Before/after safety audit for the issue #92 fix (benzene substituent
misclassification in `chematic-iupac`'s disubstituted/trisubstituted
naming paths).

Consumes two TSVs from `cargo run -p chematic-iupac --release --example
iupac_name_snapshot -- <SMILES.csv> <out.tsv>` (`smiles\\tstatus\\tname`,
status in PARSE_FAIL/OK/NOT_SUPPORTED), one generated before the fix and one
after, over the SAME corpus.

Every input SMILES is accounted for in exactly one bucket (no silent skips).
Every molecule downgraded from OK to NOT_SUPPORTED is written out AND
classified by structural root cause via RDKit -- an "unclassified" residual
means the fix rejected something for a reason this audit doesn't yet
recognize, which must be explained before merge, not shrugged off.

Usage:
    .venv/bin/python scripts/iupac_snapshot_diff.py <before.tsv> <after.tsv> \\
        [--downgraded-out downgraded.jsonl]
"""

import argparse
import json
import sys
from collections import defaultdict

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

SUPPORTED_SUBSTITUENT_ELEMENTS = {"C": "methyl", "O": "hydroxy", "N": "amino", "F": "fluoro", "Cl": "chloro", "Br": "bromo", "I": "iodo"}
EXPECTED_H = {"C": 3, "O": 1, "N": 2, "F": 0, "Cl": 0, "Br": 0, "I": 0}


def load(path):
    rows = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            smi, status = parts[0], parts[1]
            name = parts[2] if len(parts) > 2 else ""
            rows[smi] = (status, name)
    return rows


def find_pure_benzene_rings(mol):
    """6-membered aromatic all-carbon rings (the only rings the di/tri-
    substituted naming path in rings.rs ever dispatches into)."""
    ri = mol.GetRingInfo()
    rings = []
    for ring in ri.AtomRings():
        if len(ring) != 6:
            continue
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if all(a.GetIsAromatic() and a.GetSymbol() == "C" for a in atoms):
            rings.append(set(ring))
    return rings


def classify_downgrade(smi):
    """Best-effort single-tag classification of why a molecule that used to
    get a (wrong) name from the di/tri-substituted path now correctly fails.
    Mirrors classify_simple_benzene_substituent's rejection order."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return ["rdkit_parse_fail"]

    rings = find_pure_benzene_rings(mol)
    if not rings:
        return ["no_pure_benzene_ring_in_rdkit_view"]

    tags = set()
    for ring in rings:
        attach_atoms = [
            a
            for a in mol.GetAtoms()
            if a.GetIdx() in ring
            and any(nb.GetIdx() not in ring and nb.GetSymbol() != "H" for nb in a.GetNeighbors())
        ]
        if len(attach_atoms) not in (2, 3):
            continue  # not a di/tri-substituted benzene case at all
        for attach in attach_atoms:
            heavy_ext = [nb for nb in attach.GetNeighbors() if nb.GetIdx() not in ring and nb.GetSymbol() != "H"]
            if len(heavy_ext) != 1:
                tags.add("more_than_one_external_heavy_neighbor")
                continue
            first = heavy_ext[0]
            bond = mol.GetBondBetweenAtoms(attach.GetIdx(), first.GetIdx())
            if bond.GetBondType() != Chem.BondType.SINGLE:
                tags.add("non_single_attachment_bond")
                continue
            if first.GetIsAromatic():
                tags.add("aromatic_substituent")
                continue
            if first.GetFormalCharge() != 0:
                tags.add("charged_substituent")
                continue
            if first.GetIsotope() != 0:
                tags.add("isotopic_substituent")
                continue
            heavy_neighbors_of_first = [
                nb for nb in first.GetNeighbors() if nb.GetSymbol() != "H" and nb.GetIdx() != attach.GetIdx()
            ]
            if heavy_neighbors_of_first:
                tags.add("extended_substituent_beyond_first_atom")
                continue
            sym = first.GetSymbol()
            if sym not in SUPPORTED_SUBSTITUENT_ELEMENTS:
                tags.add("unsupported_substituent_element")
                continue
            total_h = first.GetTotalNumHs()
            if total_h != EXPECTED_H[sym]:
                tags.add("wrong_hydrogen_count_for_element")
                continue
            # This attach point actually matches an accepted shape -- the
            # downgrade must come from a SIBLING attach point on this ring.
    if not tags:
        tags.add("unclassified")
    return sorted(tags)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("before_path")
    p.add_argument("after_path")
    p.add_argument("--downgraded-out", default=None)
    args = p.parse_args()

    before = load(args.before_path)
    after = load(args.after_path)

    if set(before) != set(after):
        print("FATAL: before/after TSVs do not cover the same SMILES set", file=sys.stderr)
        sys.exit(1)

    counts = defaultdict(int)
    tag_counts = defaultdict(int)
    downgraded_rows = []
    changed_name_examples = []
    newly_ok_examples = []

    for smi, (b_status, b_name) in before.items():
        a_status, a_name = after[smi]
        counts["input"] += 1

        if b_status == "PARSE_FAIL" or a_status == "PARSE_FAIL":
            counts["parse_fail"] += 1
            continue
        if b_status == "OK":
            counts["before_ok"] += 1
        if a_status == "OK":
            counts["after_ok"] += 1

        if b_status == "OK" and a_status == "NOT_SUPPORTED":
            counts["downgraded_to_not_supported"] += 1
            tags = classify_downgrade(smi)
            for t in tags:
                tag_counts[t] += 1
            downgraded_rows.append({"smiles": smi, "before_name": b_name, "tags": tags})
        elif b_status != "OK" and a_status == "OK":
            counts["newly_ok"] += 1
            if len(newly_ok_examples) < 20:
                newly_ok_examples.append((smi, a_name))
        elif b_status == "OK" and a_status == "OK":
            if b_name == a_name:
                counts["still_ok_same_name"] += 1
            else:
                counts["still_ok_changed_name"] += 1
                if len(changed_name_examples) < 20:
                    changed_name_examples.append((smi, b_name, a_name))

    accounted = (
        counts["parse_fail"]
        + counts["downgraded_to_not_supported"]
        + counts["newly_ok"]
        + counts["still_ok_same_name"]
        + counts["still_ok_changed_name"]
        + sum(
            1
            for smi, (b_status, _) in before.items()
            if b_status != "PARSE_FAIL"
            and after[smi][0] != "PARSE_FAIL"
            and b_status != "OK"
            and after[smi][0] != "OK"
        )
    )

    print(f"input:                        {counts['input']}")
    print(f"parse_fail:                   {counts['parse_fail']} (gate: == 0)")
    print(f"before_ok:                    {counts['before_ok']}")
    print(f"after_ok:                     {counts['after_ok']}")
    print(f"downgraded_to_not_supported:  {counts['downgraded_to_not_supported']} (intended safety fix, gate: not required to be 0)")
    print(f"newly_ok:                     {counts['newly_ok']} (gate: == 0)")
    print(f"still_ok_same_name:           {counts['still_ok_same_name']}")
    print(f"still_ok_changed_name:        {counts['still_ok_changed_name']} (gate: == 0)")
    assert accounted == counts["input"], f"every input SMILES must land in exactly one bucket: accounted={accounted} != input={counts['input']}"
    print()

    if tag_counts:
        print("downgraded_to_not_supported root-cause breakdown (a row can carry multiple tags):")
        for tag, count in sorted(tag_counts.items(), key=lambda kv: -kv[1]):
            print(f"  {tag}: {count}")
        unclassified = tag_counts.get("unclassified", 0)
        print()
        print(f"unclassified downgrade residual: {unclassified} (gate: must be 0)")
    else:
        unclassified = 0
    print()

    if newly_ok_examples:
        print("newly_ok examples (should be empty):")
        for smi, n in newly_ok_examples:
            print(f"  {smi} -> {n!r}")
    if changed_name_examples:
        print("still_ok_changed_name examples (should be empty):")
        for smi, b, a in changed_name_examples:
            print(f"  {smi}: {b!r} -> {a!r}")

    if args.downgraded_out:
        with open(args.downgraded_out, "w") as f:
            for row in downgraded_rows:
                f.write(json.dumps(row) + "\n")
        print()
        print(f"wrote {len(downgraded_rows)} downgraded rows to {args.downgraded_out}")

    failed = counts["parse_fail"] != 0 or counts["newly_ok"] != 0 or counts["still_ok_changed_name"] != 0 or unclassified != 0
    if failed:
        print()
        print("GATE FAILED")
        sys.exit(1)


if __name__ == "__main__":
    main()
