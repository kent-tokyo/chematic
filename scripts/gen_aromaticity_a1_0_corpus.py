#!/usr/bin/env python3
"""Generate validation/aromaticity_a1_0_corpus.jsonl for Aromaticity-A1-0.

Three buckets, per the user's spec:
- false_positive: chematic marks atoms aromatic that RDKit does not. All 32
  cases are `KNOWN_BRIDGEHEAD_N_FALSE_POSITIVES` in
  crates/chematic-perception/src/aromaticity.rs's `test_known_regressions_
  from_bridgehead_n_fix` (already pinned pre-existing regression corpus --
  root-caused in SMARTS-A0 / docs/rdkit_compat.md), plus PR #86's minimal
  14-atom reproducer.
- false_negative: chematic fails to mark atoms aromatic that RDKit does.
  azulene + purine (the two #[ignore]d tests) plus the 3
  `KNOWN_ORDER_DEPENDENT_FALSE_NEGATIVES` cases.
- negative_control: molecules where chematic and RDKit already agree, used
  to catch "stopped the over-propagation, also stopped correct propagation"
  regressions in any future A1-1 fix. Includes 4-way isolation of the FP
  mechanism's two necessary conditions (bridgehead N, exocyclic C=C) and
  representative carbonyl / non-aromatic-conjugated / plain-fused cases.

Run: python3 scripts/gen_aromaticity_a1_0_corpus.py
"""

import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "validation" / "aromaticity_a1_0_corpus.jsonl"

# Kept in sync by hand with crates/chematic-perception/src/aromaticity.rs's
# KNOWN_BRIDGEHEAD_N_FALSE_POSITIVES const -- both are frozen pinned data, a
# corpus regeneration script isn't warranted for something that only changes
# if that Rust const changes (see trace_matches_ring_pi_electrons_on_corpus
# for the guard that keeps the *trace logic* honest against this same data).
FALSE_POSITIVE = [
    "C[Si](C)(C)C1=CC=C(C2=CC3=CC=CC=C3C3=NCCCN23)C=C1",
    "C1=C(C2=CC=C(CCC3=CC=CC=C3)C=C2)N2CCCN=C2C2=CC=CC=C12",
    "ClC1=CC=C(OCC2=CC3=CC=CC=C3C3=NCCCN23)C=C1",
    "N[C@@H](CC1=CC=CC=C1)C1=CC2=CC=CC=C2C2=NCCCN12",
    "CC(C)(C)C1=CC=C(C2=C(CC3=CC=CC=C3)C3=CC=CC=C3C3=NCCCN32)C=C1",
    "C[Si](C)(C)C1=CC=C(C2=C(CC3=CC=CC=C3)C3=CC=CC=C3C3=NCCCN32)C=C1",
    "C1=C(C2=CC=C(C3=CC=CC=C3)C=C2)N2CCCN=C2C2=CC=CC=C12",
    "C1=C(C2=CC=C(OCC3=CC=CC=C3)C=C2)N2CCCN=C2C2=CC=CC=C12",
    "COC1=C(OC)C(OC)=CC(C2=CC3=CC=CC=C3C3=NCCCN23)=C1",
    "CC1=CC2=CC=CC=C2C2=NCCCN12",
    "CC(C)(C)C1=CC=C(C2=CC3=C(C=C(NC(=O)NC4CCCCC4)C=C3)C3=NCCCN23)C=C1",
    "C1=CC=C(CCC2=CC=C(C3=C(CC4=CC=CC=C4)C4=CC=CC=C4C4=NCCCN43)C=C2)C=C1",
    "CCCCC1=C(C2=CC=C(CCC3=CC=CC=C3)C=C2)N2CCCN=C2C2=CC=CC=C12",
    "CCCCC1=C(C2=CC=C(C(C)(C)C)C=C2)N2CCCN=C2C2=CC=CC=C12",
    "CCCCCCC1=CC2=CC=CC=C2C2=NCCCN12",
    "CCOC1=CC=C(CC2=C(CCCC3=CC=CC4=CC=CC=C34)N3CCCN=C3C3=CC=CC=C23)C=C1",
    "CCOC1=CC=C(CC2=C(C3=CC=C(CCC4=CC=CC=C4)C=C3)N3CCCN=C3C3=CC=CC=C23)C=C1",
    "CN(C)CCC1=C(C2=CC=C(C(C)(C)C)C=C2)N2CCCN=C2C2=CC=CC=C12",
    "CC(C)(C)C1=CC=C(C2=CC3=C(C=C(N/C(S)=N/C4CCCCC4)C=C3)C3=NCCCN23)C=C1",
    "C1=C(/C=C/C2=CC=CC=C2)N2CCCN=C2C2=CC=CC=C12",
    "CC(C)(C)C1=CC=C(C2=CC3=CC=CC=C3C3=NCCCN23)C=C1",
    "CC(C)(C)C1=CC=C(C2=CC3=C(C=C(NC(=O)CC4=CC=CC=N4)C=C3)C3=NCCCN23)C=C1",
    "CC(C)(C)C1=CC=C(C2=CC3=C(C=C(NC(=O)NC4=C(Cl)C=C(Cl)C=C4)C=C3)C3=NCCCN23)C=C1",
    "C1=C(CC2=CC=CC=C2)C2=CC=CC=C2C2=NCCCN12",
    "ClC1=CC=C(C2=CC3=CC=CC=C3C3=NCCCN23)C=C1",
    "C1=C(C2=CC=CC=C2)N2CCCN=C2C2=CC=CC=C12",
    "CC(C)(C)C1=CC=C(C2=CC3=C(C=C(N(CC4=CC=CC=C4)CC4=CC=CC=C4)C=C3)C3=NCCCN23)C=C1",
    "CC(C)(C)C1=CC=C(C2=CC3=C(C=C(N)C=C3)C3=NCCCN23)C=C1",
    "CC1=C2C(=NC=C1)N(C1CC1)C1=NC=CC=C1C(=O)N2C",
    "CC(=O)N1C2=NC=CC=C2C(=O)N(C)C2=CC=CN=C21",
    "CN1C(=O)C2=CC=CN=C2N(C(C)(C)C)C2=NC=CC=C21",
    "CCCN1C2=NC=CC=C2C(=O)N(C)C2=CC=CN=C21",
]

FALSE_POSITIVE_EXTRA = [
    ("C1=Cc2ccccc2C2=NCCCN12", "PR #86 minimal 14-atom bare-core reproducer"),
]

# Kept in sync by hand with KNOWN_ORDER_DEPENDENT_FALSE_NEGATIVES.
FALSE_NEGATIVE_ORDER_DEPENDENT = [
    "N1=C2C(N(CC(O)=O)C(=O)N=C2N(C2C=C(C(F)(F)F)C=C(C=2)C(F)(F)F)C2C1=CC=CC=2)=O",
    "[C@H]12N(C([C@H](NC(=O)[C@H]([C@H](OC(=O)[C@@H](N(C)C(CN(C)C1=O)=O)C(C)C)C)NC(=O)C1C=C(OC)C(C)=C3OC4=C(C)C(=O)C(=C(C4=NC=13)C(=O)N[C@H]1C(=O)N[C@@H](C(C)C)C(N3[C@H](C(=O)N(CC(N([C@H](C(C)C)C(O[C@H]1C)=O)C)=O)C)CCC3)=O)N)C(C)C)=O)CCC2",
    "C12N(C3C=CC=CC=3)C3=NC(=O)N(C)C(C3=NC1=CC=CC=2)=O",
]

FALSE_NEGATIVE_EXTRA = [
    ("C1=CC2=CC=CC=CC2=C1", "azulene, Kekulized -- #[ignore]d test_azulene_kekulized_aromatic"),
    ("c1cnc2[nH]cnc2n1", "purine -- #[ignore]d test_purine_aromatic"),
]

NEGATIVE_CONTROL = [
    ("C1=Cc2ccccc2C2=CCCC12", "bridgehead N removed (all-carbon 3rd ring) -- FP mechanism needs both conditions"),
    ("C1=Cc2ccccc2C2=CCNC12", "N present but not at the bridgehead -- FP mechanism needs both conditions"),
    ("C1Cc2ccccc2C2=NCCCN12", "bridgehead N kept, exocyclic C=C saturated -- FP mechanism needs both conditions"),
    ("c1ccc2ccccc2c1", "naphthalene -- plain linear-fused benzenoid, must stay correct"),
    ("C1CCc2ccccc2C1", "tetralin -- benzo-fused saturated ring, must stay correct"),
    ("c1ccc2[nH]ccc2c1", "indole -- fused 5+6 heteroaromatic, must stay correct"),
    ("c1ccc2ncccc2c1", "quinoline -- fused 6+6 heteroaromatic, must stay correct"),
    ("c1ccccc1C(C)=O", "acetophenone -- benzene + exocyclic C=O, ruled out as the FP mechanism"),
    ("O=C1CCc2ccccc21", "indanone -- fused bicyclic + exocyclic C=O, ruled out as the FP mechanism"),
    ("O=C1C=Nc2ccccc21", "heteroaromatic-junction carbonyl, ruled out as the FP mechanism"),
    ("C1=CC=CC=C1", "1,3-cyclohexadiene -- non-aromatic conjugated ring, must stay non-aromatic"),
    ("C1=CC=CC=CC1", "cycloheptatriene -- non-aromatic (4n) conjugated ring, must stay non-aromatic"),
    # Added for A1-1a: does the atom-8-type exocyclic-carbonyl-in-ring rule
    # (CarbonExocyclicHeteroatomDouble -> 0pi, the OTHER atom in the false-
    # positive scaffold's 6pi sum) match RDKit at all, or is it itself the
    # thing to fix? Advisor-flagged blocking check for A1-1a's fix-route
    # decision -- confirmed against real RDKit (see docs/rfcs/aromaticity_a1_rfc.md):
    # RDKit marks the WHOLE ring (including the carbonyl carbon) aromatic in
    # all three, so this rule is correct and NOT the fix target.
    ("O=c1cccccc1", "tropone -- exocyclic-carbonyl-IN-RING rule, RDKit agrees, must stay aromatic"),
    ("O=c1cccc[nH]1", "2-pyridone -- exocyclic-carbonyl-IN-RING rule, RDKit agrees, must stay aromatic"),
    ("O=c1ccocc1", "4-pyranone -- exocyclic-carbonyl-IN-RING rule, RDKit agrees, must stay aromatic"),
    (
        "c1ccn2ccccc12",
        "indolizine -- TRUE bridgehead N shared by two BOTH-valid rings (contrast with the FP scaffold, where the bridgehead N's other ring is sp3-broken)",
    ),
    ("c1ccc2cc3ccccc3cc2c1", "anthracene -- 3 linearly fused benzo rings, must stay correct"),
]


def main():
    rows = []
    for smi in FALSE_POSITIVE:
        rows.append({"bucket": "false_positive", "smiles": smi, "note": "test_known_regressions_from_bridgehead_n_fix"})
    for smi, note in FALSE_POSITIVE_EXTRA:
        rows.append({"bucket": "false_positive", "smiles": smi, "note": note})
    for smi in FALSE_NEGATIVE_ORDER_DEPENDENT:
        rows.append({"bucket": "false_negative", "smiles": smi, "note": "test_known_order_dependent_regressions"})
    for smi, note in FALSE_NEGATIVE_EXTRA:
        rows.append({"bucket": "false_negative", "smiles": smi, "note": note})
    for smi, note in NEGATIVE_CONTROL:
        rows.append({"bucket": "negative_control", "smiles": smi, "note": note})

    for i, row in enumerate(rows):
        row["case_id"] = f"{row['bucket']}-{i:03d}"

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        for row in rows:
            f.write(json.dumps(row) + "\n")

    by_bucket = {}
    for row in rows:
        by_bucket[row["bucket"]] = by_bucket.get(row["bucket"], 0) + 1
    print(f"wrote {len(rows)} rows to {OUT}")
    for bucket, n in by_bucket.items():
        print(f"  {bucket}: {n}")


if __name__ == "__main__":
    main()
