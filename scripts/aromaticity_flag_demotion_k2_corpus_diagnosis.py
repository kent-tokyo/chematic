#!/usr/bin/env python3
"""fix/aromaticity-flag-demotion-k2: before/after diff of `apply_aromaticity`'s
per-atom/per-bond aromatic flags across the 5000-molecule descriptor-census
corpus (`scripts/descriptor_census_corpus.smi`, reused from PR #137 -- not
rebuilt for this PR), for both documented calling conventions
(`apply_aromaticity`'s own doc comment: "may be kekulized... or may retain
Aromatic bond orders from the SMILES parser").

Reads two JSONL dumps produced by
`crates/chematic-perception/examples/aromaticity_flag_demotion_k2_corpus.rs`
-- one run against the pre-fix crate, one against the post-fix crate (same
example file both times; it only calls public API, so it is
behavior-neutral with respect to the fix itself) -- and reports, per
pathway ("raw": `apply_aromaticity(&raw)` directly; "kekulized":
`chematic_core::kekulize`+`apply_kekule` first, then `apply_aromaticity`):

- how many atoms/bonds/molecules flip, and in which direction (demotion vs
  promotion);
- the self-consistency invariant (checked in BOTH directions: no
  `Aromatic`-order bond may have a non-aromatic endpoint atom, and no
  `aromatic: true` ring atom may have every one of its ring bonds be a
  non-aromatic order) before and after;
- explicitly, RDKit-agreement counted separately from self-consistency:
  atom-level demotions on the "kekulized" pathway make some molecules that
  used to coincidentally AGREE with RDKit (stale flag survived, matched by
  luck) now DISAGREE with RDKit (correctly demoted, but RDKit itself still
  says aromatic) -- this is a genuine, expected, RDKit-agreement regression
  on this specific pathway/corpus, not swept into a net "0 regressions"
  framing. See the K2 PR description for why (RFC §6's own explicit
  warning: fixing the promote-only bug is "not a pure strictly-additive
  improvement").

Run:
    .venv/bin/python scripts/aromaticity_flag_demotion_k2_corpus_diagnosis.py \\
        validation/results/aromaticity_flag_demotion_k2_corpus_before.jsonl \\
        validation/results/aromaticity_flag_demotion_k2_corpus_after.jsonl
"""

import json
import sys
from pathlib import Path

from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
SUMMARY_PATH = ROOT / "validation" / "results" / "aromaticity_flag_demotion_k2_corpus_summary.json"


def load(path):
    out = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            row = json.loads(line)
            out[row["smiles"]] = row
    return out


def analyze(before, after, pathway):
    n_atom_changed_mols = 0
    n_bond_changed_mols = 0
    atom_flip_count = 0
    bond_flip_count = 0
    demotions = 0
    promotions = 0
    before_inconsistent = 0
    after_inconsistent = 0
    became_consistent = 0
    became_inconsistent = 0
    changed_smiles = []

    for smi in before:
        b = before[smi][pathway]
        a = after[smi][pathway]
        b_atoms = {x["idx"]: x["aromatic"] for x in b["atoms"]}
        a_atoms = {x["idx"]: x["aromatic"] for x in a["atoms"]}
        b_bonds = {(x["a1"], x["a2"]): x["aromatic"] for x in b["bonds"]}
        a_bonds = {(x["a1"], x["a2"]): x["aromatic"] for x in a["bonds"]}
        atom_diff = [i for i in b_atoms if b_atoms[i] != a_atoms.get(i)]
        bond_diff = [k for k in b_bonds if b_bonds[k] != a_bonds.get(k)]
        if atom_diff:
            n_atom_changed_mols += 1
            atom_flip_count += len(atom_diff)
            for i in atom_diff:
                if b_atoms[i] and not a_atoms[i]:
                    demotions += 1
                elif not b_atoms[i] and a_atoms[i]:
                    promotions += 1
            changed_smiles.append(smi)
        if bond_diff:
            n_bond_changed_mols += 1
            bond_flip_count += len(bond_diff)
        if not b["consistent"]:
            before_inconsistent += 1
        if not a["consistent"]:
            after_inconsistent += 1
        if not b["consistent"] and a["consistent"]:
            became_consistent += 1
        if b["consistent"] and not a["consistent"]:
            became_inconsistent += 1

    return {
        "n_molecules": len(before),
        "molecules_with_atom_flag_change": n_atom_changed_mols,
        "molecules_with_bond_flag_change": n_bond_changed_mols,
        "atom_flips": atom_flip_count,
        "atom_demotions": demotions,
        "atom_promotions": promotions,
        "bond_flips": bond_flip_count,
        "inconsistent_before": before_inconsistent,
        "inconsistent_after": after_inconsistent,
        "became_consistent": became_consistent,
        "became_inconsistent_REGRESSION": became_inconsistent,
        "changed_smiles": changed_smiles,
    }


def rdkit_atom_aromatic(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return [a.GetIsAromatic() for a in mol.GetAtoms()]


def rdkit_agreement_check(before, after, pathway, changed_smiles):
    """For every molecule whose atom flags changed on this pathway, check
    RDKit-agreement before vs after directly (not inferred from
    self-consistency) -- classifies each changed molecule as
    newly-agrees-with-rdkit, newly-disagrees-with-rdkit (a real regression
    in RDKit-agreement terms even though self-consistency improved), or
    other (agreement unchanged, e.g. RDKit itself couldn't be compared)."""
    newly_agrees = []
    newly_disagrees = []
    other = []
    for smi in changed_smiles:
        rd = rdkit_atom_aromatic(smi)
        if rd is None:
            other.append(smi)
            continue
        b_atoms = before[smi][pathway]["atoms"]
        a_atoms = after[smi][pathway]["atoms"]
        if len(rd) != len(b_atoms):
            other.append(smi)
            continue
        b_agree = all(b_atoms[i]["aromatic"] == rd[i] for i in range(len(rd)))
        a_agree = all(a_atoms[i]["aromatic"] == rd[i] for i in range(len(rd)))
        if not b_agree and a_agree:
            newly_agrees.append(smi)
        elif b_agree and not a_agree:
            newly_disagrees.append(smi)
        else:
            other.append(smi)
    return {
        "checked": len(changed_smiles),
        "newly_agrees_with_rdkit": len(newly_agrees),
        "newly_disagrees_with_rdkit_REGRESSION": len(newly_disagrees),
        "other_or_unchanged_agreement": len(other),
        "newly_disagrees_examples": newly_disagrees[:10],
    }


def main():
    before_path = sys.argv[1] if len(sys.argv) > 1 else "validation/results/aromaticity_flag_demotion_k2_corpus_before.jsonl"
    after_path = sys.argv[2] if len(sys.argv) > 2 else "validation/results/aromaticity_flag_demotion_k2_corpus_after.jsonl"

    before = load(before_path)
    after = load(after_path)
    assert set(before) == set(after), "before/after dumps must cover the identical SMILES set"

    result = {"n_molecules": len(before), "pathways": {}}
    for pathway in ("raw", "kekulized"):
        r = analyze(before, after, pathway)
        print(f"--- pathway: {pathway} ---")
        for k, v in r.items():
            if k == "changed_smiles":
                continue
            print(f"  {k}: {v}")
        if r["atom_demotions"] > 0:
            print("  cross-checking RDKit agreement on every changed molecule (this may take a moment)...")
            rd_check = rdkit_agreement_check(before, after, pathway, r["changed_smiles"])
            print(f"  rdkit_agreement_check: {rd_check}")
            r["rdkit_agreement_check"] = rd_check
        result["pathways"][pathway] = r

    SUMMARY_PATH.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_PATH.write_text(json.dumps(result, indent=2))
    print(f"\nwrote {SUMMARY_PATH}")


if __name__ == "__main__":
    main()
