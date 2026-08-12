#!/usr/bin/env python3
"""RDKit-oracle derivation of the @SP1/@SP2/@SP3 permutation-remap rule.

The OpenSMILES spec defers the precise geometric rule for square-planar
stereo (which pair of the 4 SMILES-order neighbor positions sit *trans* to
each other, for each of SP1/SP2/SP3) to a diagram that isn't reproduced in
the machine-readable spec text. This script derives the rule empirically
against chematic's own pinned RDKit oracle instead of hand-deriving it from
an unavailable diagram -- matching this project's established
oracle-verification convention (see validation/platinum/, scripts/
platinum_rdkit_oracle.py, etc.).

Method: for each of 2 differently-shaped test molecules (a simple
4-distinct-ligand center, and a dative-bond mix matching the cisplatin/
transplatin corpus's own donor/acceptor grammar), enumerate all 4! = 24
neighbor-order permutations x 3 tags = 72 SMILES variants each. For every
variant, predict (via the pair-of-pairs rule below) which tag the SAME
physical molecule needs when re-spelled with neighbors in the fixed
reference order 0,1,2,3, then confirm RDKit's own canonical SMILES agrees
that the two spellings describe the same molecule.

The derived rule (also implemented independently in Rust,
crates/chematic-smiles/src/canonical.rs's `remap_square_planar`): each tag
names a partition of positions {0,1,2,3} into its two trans-pairs
(SP1={0,2}|{1,3}, SP2={0,1}|{2,3}, SP3={0,3}|{1,2}). Applying a neighbor
permutation to that pair-of-pairs and matching the result against the 3
templates predicts which tag a reordered SMILES needs to describe the same
molecule -- this script's job is to confirm that prediction holds for every
one of the 24 x 3 x 2 = 144 cases, not to assume it.

Scope note: this script deliberately covers branch-only neighbor ordering
(the general case `remap_square_planar` is written against -- it operates on
whatever `original`/`canonical` neighbor-id sequences it's given, with no
assumption about *how* those sequences were produced). Ring-closure-based
neighbor ordering (chelate rings, e.g. carboplatin/oxaliplatin's grammar) is
NOT re-derived here: an earlier version of this script modeled ring-closure
SMILES-encounter order incorrectly (assumed the digit's *closing* occurrence
sets the encounter position; RDKit's actual convention disagreed), which
would have silently tested a different permutation than the one intended.
Rather than reverse-engineer RDKit's internal ring-closure-encounter-order
convention in Python, that case is verified directly and more reliably at
the Rust level instead -- see
`crates/chematic-smiles/tests/square_planar_stereo.rs`'s chelate-shaped
fixture, which round-trips a real ring-closure SMILES through chematic's
own parser and canonical writer and checks the result against RDKit,
with no intermediate Python re-modeling step to get wrong.

Usage: python scripts/square_planar_permutation_oracle.py
"""

import itertools
import sys

from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

TAGS = ["SP1", "SP2", "SP3"]

# Each shape names 4 distinct ligand *elements* and which of them are dative
# donors. A plain (covalent) ligand's SMILES token is the same regardless of
# position. A dative ligand's token is position-dependent: written BEFORE Pt
# it's a suffix arrow ("N->"), written AFTER Pt (branch or trailing chain)
# it's a prefix arrow ("<-N") -- the arrow always points donor -> acceptor,
# so it flips depending on which side of "[Pt...]" the donor ends up on
# after permutation.
SHAPES = {
    "simple_4_distinct_ligands": {
        "ligands": ["C", "F", "Cl", "[H]"],
        "dative": set(),
    },
    "dative_bond_mix": {
        # matches the cisplatin/transplatin corpus's own donor/acceptor
        # grammar (N/O donors, halide acceptors).
        "ligands": ["N", "O", "Cl", "Br"],
        "dative": {"N", "O"},
    },
}

TEMPLATE = "{a0}[Pt@{tag}]({a1})({a2}){a3}"


def ligand_token(element, slot_index, dative_set):
    """Render `element` as it should appear at SMILES slot 0 (before Pt) or
    slots 1-3 (after Pt, branch or trailing chain)."""
    if element not in dative_set:
        return element
    return f"{element}->" if slot_index == 0 else f"<-{element}"


def build_smiles(shape_key, order, tag):
    shape = SHAPES[shape_key]
    ligands = shape["ligands"]
    dative_set = shape["dative"]
    tokens = [
        ligand_token(ligands[ligand_idx], slot, dative_set)
        for slot, ligand_idx in enumerate(order)
    ]
    return TEMPLATE.format(a0=tokens[0], a1=tokens[1], a2=tokens[2], a3=tokens[3], tag=tag)


def trans_pairs(tag):
    return {
        "SP1": [(0, 2), (1, 3)],
        "SP2": [(0, 1), (2, 3)],
        "SP3": [(0, 3), (1, 2)],
    }[tag]


def predict_remap(tag, original_order, canonical_order):
    """Predict the new tag when neighbors move from `original_order` to
    `canonical_order` (both permutations of range(4), naming which physical
    ligand-index sits at each SMILES position)."""
    pos_in_canonical = {v: i for i, v in enumerate(canonical_order)}
    new_pairs = []
    for i, j in trans_pairs(tag):
        a, b = pos_in_canonical[original_order[i]], pos_in_canonical[original_order[j]]
        new_pairs.append(tuple(sorted((a, b))))
    new_pairs = tuple(sorted(new_pairs))
    for candidate in TAGS:
        if tuple(sorted(trans_pairs(candidate))) == new_pairs:
            return candidate
    return None


def main():
    total = 0
    mismatches = []
    reference_order = (0, 1, 2, 3)

    for shape_key in SHAPES:
        for order in itertools.permutations(range(4)):
            for tag in TAGS:
                total += 1
                smi = build_smiles(shape_key, order, tag)
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    mismatches.append((shape_key, order, tag, "PARSE_FAILED", smi))
                    continue
                canon = Chem.MolToSmiles(mol)

                predicted_tag = predict_remap(tag, list(order), list(reference_order))
                if predicted_tag is None:
                    mismatches.append((shape_key, order, tag, "NO_PREDICTION", smi))
                    continue
                reference_smi = build_smiles(shape_key, reference_order, predicted_tag)
                ref_mol = Chem.MolFromSmiles(reference_smi)
                if ref_mol is None:
                    mismatches.append(
                        (shape_key, order, tag, "REFERENCE_PARSE_FAILED", reference_smi)
                    )
                    continue
                ref_canon = Chem.MolToSmiles(ref_mol)

                if canon != ref_canon:
                    mismatches.append(
                        (
                            shape_key,
                            order,
                            tag,
                            f"MISMATCH predicted={predicted_tag} got_canon={canon} ref_canon={ref_canon}",
                            smi,
                        )
                    )

    print(
        f"Checked {total} (shape x permutation x tag) cases across {len(SHAPES)} molecule shapes."
    )
    if mismatches:
        print(f"FAILED: {len(mismatches)} mismatches")
        for m in mismatches[:20]:
            print(f"  {m}")
        sys.exit(1)
    print("All cases matched the predicted remap rule. 0 mismatches.")


if __name__ == "__main__":
    main()
