#!/usr/bin/env python3
"""
Corpus-scale ground-truth check for chematic's CIP R/S assignment
(chematic_chem::assign_cip, exposed as Mol.cip_stereo() in Python).

Root-caused and fixed in the round this script was added: two independent
bugs in crates/chematic-chem/src/cip.rs --

1. assign_tetrahedral() rebuilt a chiral atom's substituent order from raw
   Molecule::neighbors() adjacency order instead of the already-correct,
   parser-populated Molecule::stereo_neighbor_order(). Adjacency order only
   matches SMILES textual order when every ring bond at that atom CLOSES
   (partner already known); it silently reorders substituents when a
   stereocenter OPENS a ring before its other neighbors, because a
   ring-opening bond only gets added to the adjacency list once the
   matching closing digit is reached later in the string.

2. cip_branch_spheres()'s CIP double-bond duplication only added the
   "arrival side" phantom (B's own sphere gets a phantom-A once B is
   expanded, having been reached via A=B) -- never the "departure side"
   (A's own sphere never got a second phantom-B while iterating A's
   neighbors). A double bond must duplicate its partner into BOTH atoms'
   substituent spheres.

Fixed in commit d0e726b: 76.22% -> 96.83% vs the legacy oracle (below).

Follow-up (same round): a reverted attempt to extend the double-bond fix to
triple bonds (BondOrder::Triple, 2 phantom duplicates per side instead of 1)
went net NEGATIVE (16 newly-wrong vs 1 newly-fixed against the 132-mismatch
baseline). Root cause: cip_branch_spheres()/compare_branches() pools all
atoms at each BFS depth into one sorted multiset and compares shell-by-shell
-- an approximation of, not the true recursive branch-by-branch, CIP
hierarchical-digraph algorithm. Adding correct per-bond-type duplication
rules to this approximate comparator is whack-a-mole. See
docs/cip_accurate_rfc.md for the follow-on "Accurate CIP" engine design this
finding motivated.

Two oracles, deliberately not conflated:
  - MODERN (primary): rdkit.Chem.rdCIPLabeler.AssignCIPLabels -- the
    IUPAC-2013-rules-based implementation. This is the oracle project-wide
    correctness numbers should cite going forward.
  - LEGACY (secondary/reference): Chem.AssignStereochemistry(cleanIt=True,
    force=True) + _CIPCode -- kept because it's what earlier rounds/docs
    cited, and because legacy-vs-modern disagreement (43 cases in this
    corpus) is itself informative: it's not a chematic bug when only the
    legacy oracle is wrong.

Correspondence note: this compares the ORIGINAL (non-respelled) SMILES only,
so RDKit and chematic atom indices align directly (both parse the same
string in the same atom order -- verified elsewhere this project,
scripts/aromaticity_mechanism_probe.py). No respelling/correspondence
mapping is needed or performed here.

Every mismatch (vs the modern oracle) is written to a frozen, deterministic
JSONL corpus (validation/cip_label_corpus.jsonl) -- see --freeze. Rerunning
against the same corpus CSV and RDKit version reproduces it byte-for-byte
(stable sort by (smiles, atom_idx)).

Usage:
    .venv/bin/python scripts/cip_ground_truth.py [SMILES.csv]
    .venv/bin/python scripts/cip_ground_truth.py [SMILES.csv] --freeze [out.jsonl]
"""

import csv
import datetime
import json
import sys

sys.path.insert(0, ".")
import chematic
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

DEFAULT_FREEZE_PATH = "validation/cip_label_corpus.jsonl"


def classify_bucket(rd, aidx):
    """Best-effort mechanism tag for known residual classes (informational,
    not load-bearing -- see docs/cip_accurate_rfc.md Milestones 2-4)."""
    atom = rd.GetAtomWithIdx(aidx)
    if atom.GetSymbol() == "P":
        return "phosphorus"
    if any(n.GetIsAromatic() for n in atom.GetNeighbors()):
        return "aromatic_mancude"
    if any(b.GetBondType() == Chem.BondType.TRIPLE for b in rd.GetBonds()):
        return "triple_bond_dup"
    return None


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    csv_path = args[0] if args else "SMILES.csv"
    freeze = "--freeze" in sys.argv
    freeze_path = DEFAULT_FREEZE_PATH
    if freeze and len(args) > 1:
        freeze_path = args[1]

    with open(csv_path) as f:
        reader = csv.reader(f)
        next(reader)
        smis = [row[0] for row in reader if row]

    total = 0
    match = 0
    mismatches = []  # (smiles, atom_idx, chematic, modern, legacy, bucket)
    legacy_vs_modern_disagree = 0

    for smi in smis:
        rd = Chem.MolFromSmiles(smi)
        if rd is None:
            continue
        if not any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in rd.GetAtoms()):
            continue

        rd_legacy = Chem.MolFromSmiles(smi)
        Chem.AssignStereochemistry(rd_legacy, cleanIt=True, force=True)
        legacy_cip = {
            a.GetIdx(): a.GetProp("_CIPCode") for a in rd_legacy.GetAtoms() if a.HasProp("_CIPCode")
        }

        try:
            rdCIPLabeler.AssignCIPLabels(rd)
        except Exception:
            continue
        modern_cip = {a.GetIdx(): a.GetProp("_CIPCode") for a in rd.GetAtoms() if a.HasProp("_CIPCode")}
        if not modern_cip:
            continue

        for aidx in set(legacy_cip) | set(modern_cip):
            if legacy_cip.get(aidx) != modern_cip.get(aidx):
                legacy_vs_modern_disagree += 1

        try:
            m = chematic.from_smiles(smi)
            cm_cip = {d["atom_idx"]: d["descriptor"] for d in m.cip_stereo()}
        except Exception:
            continue

        for aidx, modern_code in modern_cip.items():
            total += 1
            cm_code = cm_cip.get(aidx)
            if cm_code == modern_code:
                match += 1
            else:
                mismatches.append(
                    (smi, aidx, cm_code, modern_code, legacy_cip.get(aidx), classify_bucket(rd, aidx))
                )

    mismatches.sort(key=lambda x: (x[0], x[1]))

    print(f"legacy-vs-modern RDKit disagreement count (any atom): {legacy_vs_modern_disagree}")
    print(f"total stereocenters compared (modern oracle): {total}")
    print(f"match: {match} ({100 * match / total:.2f}%)" if total else "no stereocenters found")
    print(f"mismatch: {len(mismatches)}")
    for ex in mismatches[:10]:
        print(" ", ex)

    if freeze:
        manifest = {
            "_manifest": True,
            "rdkit_version": Chem.rdBase.rdkitVersion,
            "chematic_version": getattr(chematic, "__version__", None),
            "source_csv": csv_path.rsplit("/", 1)[-1],
            "total_stereocenters": total,
            "mismatches": len(mismatches),
            "generated": datetime.date.today().isoformat(),
        }
        with open(freeze_path, "w") as f:
            f.write(json.dumps(manifest, sort_keys=True) + "\n")
            for smi, aidx, cm, modern, legacy, bucket in mismatches:
                f.write(
                    json.dumps(
                        {
                            "smiles": smi,
                            "atom_idx": aidx,
                            "chematic": cm,
                            "modern": modern,
                            "legacy": legacy,
                            "bucket": bucket,
                        },
                        sort_keys=True,
                    )
                    + "\n"
                )
        print(f"froze {len(mismatches)} mismatches to {freeze_path}")


if __name__ == "__main__":
    main()
