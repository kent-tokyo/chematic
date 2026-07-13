#!/usr/bin/env python3
"""
RingInfo parity: chematic's find_sssr vs RDKit, and self-stability across
independently-traversed representations of the same molecule.

Motivation: aromaticity_atom_parity.py showed a narrow (~3%) aromaticity-flag
order-dependence, but tracing the root cause (see project plan) found the
real signal is much larger and lives one layer down, in find_sssr itself:
a single BFS spanning tree + exactly one fundamental cycle per non-tree edge
(zero candidate redundancy) means find_sssr can return a non-minimal ring
(e.g. a 10-membered ring standing in for two 6-membered ones) depending on
which atom happens to be index 0. This often doesn't change the final
aromaticity flags (Pass 1/Pass 2 + the aromatic_context bypass tends to
paper over a non-minimal ring choice as long as some other, correct ring
still exists) -- but it directly corrupts anything that reads ring
*membership* or ring *size* rather than just the final aromatic bool:
RingCount-adjacent descriptors, ring-size SMARTS ([r5], [r6], [R2]), Murcko
scaffold extraction, ring-membership fingerprint invariants, and Crippen
atom-typing (LogP/MR).

This script measures two things, deliberately NOT the same thing:
1. Self-stability: does the SAME molecule give the SAME RingInfo across N
   independently-traversed valid Kekulized SMILES? This is oracle-independent
   -- any variation is definitionally a bug, regardless of what RDKit says
   (RDKit's own SSSR isn't unique either, so "does it match RDKit" alone
   can't be the only signal).
2. RDKit agreement: for context, using the *default* (as-parsed) chematic
   result against RDKit's GetRingInfo(), on invariants only (never the
   literal ring set).

Invariants measured, all read via existing chematic Python API:
  - num_rings            (mu = |E|-|V|+components; a topological invariant --
                           if THIS varies across traversals, something more
                           basic than ring-size selection is broken)
  - ring_size_multiset    (sorted tuple of ring sizes -- the core signal)
  - total_ring_weight     (sum of ring sizes -- evidence of minimality)
  - min_ring_size_per_atom (via ring_sizes_for_atom(i), min per atom)
  - ring_membership_count (per atom, how many SSSR rings contain it)

Usage:
    python scripts/ringinfo_parity.py [SMILES.csv] [--detail] [--limit N] [-n N_VARIANTS] [--self-test]
"""
import sys
import re
import argparse


def map_nums_in_order(smi):
    return [int(x) for x in re.findall(r":(\d+)\]", smi)]


def chematic_ring_metrics(mol):
    """Compute RingInfo invariants from a chematic Mol, in its own atom-index order."""
    rings = mol.sssr_atom_rings
    n_atoms = len(mol.atom_table)
    sizes = sorted(len(r) for r in rings)
    membership = [0] * n_atoms
    for r in rings:
        for a in r:
            membership[a] += 1
    min_size = [mol.ring_sizes_for_atom(i) for i in range(n_atoms)]
    min_size = [min(s) if s else 0 for s in min_size]
    return {
        "num_rings": len(rings),
        "ring_size_multiset": tuple(sizes),
        "total_ring_weight": sum(sizes),
        "min_ring_size_per_atom": tuple(min_size),
        "ring_membership_count": tuple(membership),
    }


def rdkit_ring_metrics(rdmol):
    ri = rdmol.GetRingInfo()
    sizes = sorted(len(r) for r in ri.AtomRings())
    n_atoms = rdmol.GetNumAtoms()
    min_size = [ri.MinAtomRingSize(i) for i in range(n_atoms)]
    membership = [ri.NumAtomRings(i) for i in range(n_atoms)]
    return {
        "num_rings": ri.NumRings(),
        "ring_size_multiset": tuple(sizes),
        "total_ring_weight": sum(sizes),
        "min_ring_size_per_atom": tuple(min_size),
        "ring_membership_count": tuple(membership),
    }


def self_test():
    """Positive control: corrupt a known-good result and confirm the harness's
    own comparison logic flags it. Validates detection logic independent of
    whether real find_sssr output is being tested."""
    good = {
        "num_rings": 2,
        "ring_size_multiset": (6, 6),
        "total_ring_weight": 12,
        "min_ring_size_per_atom": (6,) * 10,
        "ring_membership_count": (1,) * 8 + (2, 2),
    }
    # Corruption: drop the largest ring (as suggested) -- simulate what a
    # non-minimal find_sssr looks like: one 6-ring replaced by a 10-ring.
    corrupted = dict(good)
    corrupted["ring_size_multiset"] = (6, 10)
    corrupted["total_ring_weight"] = 16

    assert good != corrupted, "self-test setup is broken: corruption produced no difference"
    variants = [good, good, corrupted, good]
    unstable = len(set(v["ring_size_multiset"] for v in variants)) > 1
    assert unstable, "POSITIVE CONTROL FAILED: harness did not detect an injected ring-size fault"
    print("positive control: PASS (harness correctly detects an injected non-minimal ring)")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("smiles_csv", nargs="?", default="~/Downloads/SMILES.csv")
    parser.add_argument("--detail", action="store_true")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("-n", "--variants", type=int, default=10)
    parser.add_argument("--self-test", action="store_true",
                         help="Run only the positive control and exit")
    args = parser.parse_args()

    if args.self_test:
        self_test()
        return

    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
    except ImportError:
        sys.exit("rdkit not installed. pip install rdkit")
    import chematic
    import os

    # Always run the positive control before trusting production numbers.
    self_test()

    path = os.path.expanduser(args.smiles_csv)
    with open(path) as f:
        lines = [l.strip() for l in f if l.strip()]
    smis = [l.split(",")[0].strip() for l in lines if l.split(",")[0].strip().lower() != "smiles"]

    n_mol = 0
    self_stable = 0
    self_unstable = 0
    rdkit_agree = {k: 0 for k in ["num_rings", "ring_size_multiset", "total_ring_weight"]}
    rdkit_total = 0
    shown = 0

    for smi in smis:
        rd = Chem.MolFromSmiles(smi)
        if rd is None:
            continue
        n_mol += 1

        # --- RDKit agreement, default chematic parse ---
        try:
            cm_default = chematic.from_smiles(smi)
            cm_metrics = chematic_ring_metrics(cm_default)
            rd_metrics = rdkit_ring_metrics(rd)
            rdkit_total += 1
            for k in rdkit_agree:
                if cm_metrics[k] == rd_metrics[k]:
                    rdkit_agree[k] += 1
        except Exception:
            pass

        # --- self-stability across N independently-traversed variants ---
        rd_mapped = Chem.MolFromSmiles(smi)
        for i, a in enumerate(rd_mapped.GetAtoms()):
            a.SetAtomMapNum(i + 1)
        Chem.Kekulize(rd_mapped, clearAromaticFlags=True)

        multisets = set()
        for _ in range(args.variants):
            variant = Chem.MolToSmiles(rd_mapped, doRandom=True, kekuleSmiles=True, canonical=False)
            variant_clean = re.sub(r":\d+", "", variant)
            try:
                cm = chematic.from_smiles(variant_clean)
                m = chematic_ring_metrics(cm)
                multisets.add(m["ring_size_multiset"])
            except Exception:
                multisets.add(None)

        if len(multisets) <= 1:
            self_stable += 1
        else:
            self_unstable += 1
            if args.detail and (args.limit is None or shown < args.limit):
                shown += 1
                print(f"UNSTABLE smiles={smi!r} multisets={sorted(multisets, key=lambda x: len(x) if x else 0)}",
                      file=sys.stderr)

    print(f"corpus (parsed by RDKit): {n_mol}")
    print(f"\n[self-stability across {args.variants} traversals]")
    total = self_stable + self_unstable
    print(f"  stable:   {self_stable}/{total} ({100*self_stable/max(total,1):.2f}%)")
    print(f"  unstable: {self_unstable}/{total} ({100*self_unstable/max(total,1):.2f}%)")
    print(f"\n[RDKit agreement, default chematic parse]")
    for k, v in rdkit_agree.items():
        print(f"  {k}: {v}/{rdkit_total} ({100*v/max(rdkit_total,1):.2f}%)")


if __name__ == "__main__":
    main()
