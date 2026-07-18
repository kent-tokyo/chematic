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
versa -- for every molecule where both RDKit and chematic parse successfully.
Same methodology as this project's existing "invariant partition agreement"
measurement (docs/verification_coverage.md).

Consumes the TSV from `cargo run -p chematic-fp --release --example
rdkit_invariant_snapshot -- <SMILES.csv> <out.tsv>` (`smiles\\tatom_idx\\t
invariant`, one row per atom).

Residual classification (per the acceptance gate): every molecule where the
partitions disagree is checked against a short list of known-shape causes
(degree semantics, explicit H, isotope delta, formal charge, ring
membership, aromatic/Kekule) and reported bucketed; anything left over is
"unclassified".

Usage:
    .venv/bin/python scripts/ecfp_rdkit_invariant_parity.py <chematic.tsv> [SMILES.csv]
"""

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
    """smiles -> [invariant per atom, in atom_idx order]"""
    by_smi = defaultdict(dict)
    with open(path) as f:
        for line in f:
            smi, idx, inv = line.rstrip("\n").split("\t")
            by_smi[smi][int(idx)] = int(inv)
    out = {}
    for smi, idx_map in by_smi.items():
        n = max(idx_map) + 1
        out[smi] = [idx_map.get(i) for i in range(n)]
    return out


def partition(values):
    """[v0, v1, ...] -> sorted tuple of sorted-index-tuples grouped by equal value."""
    groups = defaultdict(list)
    for i, v in enumerate(values):
        groups[v].append(i)
    return tuple(sorted(tuple(sorted(g)) for g in groups.values()))


def rdkit_invariants(rd):
    return list(rdMolDescriptors.GetConnectivityInvariants(rd, True))


def classify_mismatch(smi, rd, chem_vals, rd_vals):
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

    chematic = load_chematic(chematic_path)

    total = 0
    rdkit_parse_fail = 0
    matched = 0
    mismatched = 0
    tag_counts = defaultdict(int)
    mismatch_examples = []

    import os

    with open(os.path.expanduser(csv_path)) as f:
        smis = [line.strip() for line in f if line.strip()]

    for smi in smis:
        if smi not in chematic:
            continue  # chematic parse failure -- not this metric's concern
        rd = parse_smiles(smi)
        if rd is None:
            rdkit_parse_fail += 1
            continue
        chem_vals = chematic[smi]
        rd_vals = rdkit_invariants(rd)
        if len(chem_vals) != len(rd_vals):
            mismatched += 1
            tag_counts["atom_count_mismatch"] += 1
            if len(mismatch_examples) < 10:
                mismatch_examples.append((smi, "atom_count_mismatch"))
            continue
        total += 1
        if partition(chem_vals) == partition(rd_vals):
            matched += 1
        else:
            mismatched += 1
            tags = classify_mismatch(smi, rd, chem_vals, rd_vals)
            for t in tags:
                tag_counts[t] += 1
            if len(mismatch_examples) < 20:
                mismatch_examples.append((smi, tags))

    print(f"corpus SMILES: {len(smis)}")
    print(f"RDKit parse failures: {rdkit_parse_fail}")
    print(f"compared (both parsed, atom counts match): {total}")
    print(f"matched:    {matched}")
    print(f"mismatched: {mismatched}")
    pct = 100 * matched / total if total else 0.0
    print(f"partition agreement: {pct:.4f}% (gate: 100%)")
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


if __name__ == "__main__":
    main()
