#!/usr/bin/env python3
"""
Atom-invariant PARTITION agreement between chematic's `EcfpInvariantMode::
RdkitMorgan` and RDKit's default Morgan atom invariant
(`rdMolDescriptors.GetConnectivityInvariants`).

This is hash-VALUE-independent: chematic uses FNV-1a, RDKit uses its own
hash, so raw integer agreement is meaningless and not what's tested. What's
tested is whether the two implementations partition a molecule's atoms into
the SAME equivalence classes -- do the atoms RDKit considers
invariant-identical also come out invariant-identical in chematic, and vice
versa. Same methodology as this project's existing "invariant partition
agreement" measurement (docs/verification_coverage.md).

Every input SMILES is accounted for in exactly one bucket -- chematic parse
failures and RDKit parse failures are counted explicitly (gated on ==0), not
silently skipped, so a "100%" figure can't hide an incomplete comparison.

Consumes the TSV from `cargo run -p chematic-fp --release --example
rdkit_invariant_snapshot -- <SMILES.csv> <out.tsv>` (`smiles\\tatom_idx\\t
invariant` per atom row, or `smiles\\tPARSE_FAIL\\t<error>` for a chematic
parse failure).

Residual classification (for any partition_mismatch, which should be zero):
checked against a short list of known-shape causes (degree semantics,
explicit H, isotope delta, formal charge, ring membership, aromatic/Kekule)
and reported bucketed; anything left over is "unclassified".

Usage:
    .venv/bin/python scripts/ecfp_rdkit_invariant_parity.py <chematic.tsv> [SMILES.csv]
"""

import os
import sys
from collections import defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

_PARSE_PARAMS = Chem.SmilesParserParams()
_PARSE_PARAMS.removeHs = False


def parse_smiles(smi):
    """RDKit's default MolFromSmiles silently strips 'trivial' explicit H
    atoms, which chematic keeps as real atoms -- disable that so atom counts
    (and hence atom_idx correspondence) stay comparable."""
    return Chem.MolFromSmiles(smi, _PARSE_PARAMS)


def load_chematic(path):
    """smiles -> [invariant per atom, in atom_idx order] | "PARSE_FAIL" """
    by_smi = defaultdict(dict)
    parse_fail = set()
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            smi, tag = parts[0], parts[1]
            if tag == "PARSE_FAIL":
                parse_fail.add(smi)
                continue
            by_smi[smi][int(tag)] = int(parts[2])
    out = {}
    for smi, idx_map in by_smi.items():
        n = max(idx_map) + 1
        out[smi] = [idx_map.get(i) for i in range(n)]
    return out, parse_fail


def partition(values):
    """[v0, v1, ...] -> sorted tuple of sorted-index-tuples grouped by equal value."""
    groups = defaultdict(list)
    for i, v in enumerate(values):
        groups[v].append(i)
    return tuple(sorted(tuple(sorted(g)) for g in groups.values()))


def rdkit_invariants(rd):
    return list(rdMolDescriptors.GetConnectivityInvariants(rd, True))


def classify_mismatch(rd):
    """Best-effort single-tag classification of why the partitions differ."""
    has_explicit_h = any(a.GetSymbol() == "H" for a in rd.GetAtoms())
    has_isotope = any(a.GetIsotope() != 0 for a in rd.GetAtoms())
    has_charge = any(a.GetFormalCharge() != 0 for a in rd.GetAtoms())
    has_aromatic = any(a.GetIsAromatic() for a in rd.GetAtoms())
    has_ring = rd.GetRingInfo().NumRings() > 0
    tags = []
    if has_explicit_h:
        tags.append("explicit_hydrogen")
    if has_isotope:
        tags.append("isotope_delta")
    if has_charge:
        tags.append("formal_charge")
    if has_aromatic:
        tags.append("aromatic_kekule_representation")
    if has_ring:
        tags.append("ring_membership")
    return tags if tags else ["unclassified"]


def main():
    args = sys.argv[1:]
    if len(args) < 1:
        print(__doc__)
        sys.exit(1)
    chematic_path = args[0]
    csv_path = args[1] if len(args) > 1 else "~/Downloads/SMILES.csv"

    chematic, chematic_parse_fail = load_chematic(chematic_path)

    with open(os.path.expanduser(csv_path)) as f:
        smis = [line.strip() for line in f if line.strip()]

    counts = {
        "input_smiles": len(smis),
        "chematic_parse_fail": 0,
        "rdkit_parse_fail": 0,
        "atom_count_mismatch": 0,
        "partition_match": 0,
        "partition_mismatch": 0,
    }
    tag_counts = defaultdict(int)
    mismatch_examples = []

    for smi in smis:
        if smi in chematic_parse_fail:
            counts["chematic_parse_fail"] += 1
            continue
        if smi not in chematic:
            # Molecule chematic assigned zero atom invariants for (should not
            # happen -- assign a row for every heavy atom) -- treat as a
            # chematic-side gap, not silently dropped.
            counts["chematic_parse_fail"] += 1
            continue

        rd = parse_smiles(smi)
        if rd is None:
            counts["rdkit_parse_fail"] += 1
            continue

        chem_vals = chematic[smi]
        rd_vals = rdkit_invariants(rd)
        if len(chem_vals) != len(rd_vals):
            counts["atom_count_mismatch"] += 1
            if len(mismatch_examples) < 10:
                mismatch_examples.append((smi, ["atom_count_mismatch"]))
            continue

        if partition(chem_vals) == partition(rd_vals):
            counts["partition_match"] += 1
        else:
            counts["partition_mismatch"] += 1
            tags = classify_mismatch(rd)
            for t in tags:
                tag_counts[t] += 1
            if len(mismatch_examples) < 20:
                mismatch_examples.append((smi, tags))

    accounted = sum(
        counts[k]
        for k in (
            "chematic_parse_fail",
            "rdkit_parse_fail",
            "atom_count_mismatch",
            "partition_match",
            "partition_mismatch",
        )
    )

    print(f"input_smiles:         {counts['input_smiles']}")
    print(f"chematic_parse_fail:  {counts['chematic_parse_fail']} (gate: == 0)")
    print(f"rdkit_parse_fail:     {counts['rdkit_parse_fail']} (gate: == 0)")
    print(f"atom_count_mismatch:  {counts['atom_count_mismatch']} (gate: == 0)")
    print(f"partition_match:      {counts['partition_match']}")
    print(f"partition_mismatch:   {counts['partition_mismatch']} (gate: == 0)")
    assert accounted == counts["input_smiles"], (
        f"every input SMILES must land in exactly one bucket: "
        f"accounted={accounted} != input_smiles={counts['input_smiles']}"
    )
    denom = counts["partition_match"] + counts["partition_mismatch"]
    pct = 100 * counts["partition_match"] / denom if denom else 0.0
    print(f"partition agreement (of comparable pairs): {pct:.4f}% (gate: 100%)")
    print()
    if tag_counts:
        print("mismatch tag breakdown (a mismatch can carry multiple tags):")
        for tag, count in sorted(tag_counts.items(), key=lambda kv: -kv[1]):
            print(f"  {tag}: {count}")
        unclassified = tag_counts.get("unclassified", 0)
        print()
        print(f"unclassified residual: {unclassified} (gate: must be 0 to ship as 'RdkitMorgan')")
    print()
    for smi, tags in mismatch_examples:
        print(f"  MISMATCH [{tags}] {smi}")

    failed = (
        counts["chematic_parse_fail"] != 0
        or counts["rdkit_parse_fail"] != 0
        or counts["atom_count_mismatch"] != 0
        or counts["partition_mismatch"] != 0
        or counts["partition_match"] != counts["input_smiles"]
    )
    if failed:
        print()
        print("GATE FAILED")
        sys.exit(1)


if __name__ == "__main__":
    main()
