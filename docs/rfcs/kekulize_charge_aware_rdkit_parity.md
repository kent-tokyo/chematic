# K1: charge-aware `atom_must_be_matched` — RDKit source-grounded rationale

**Status:** implemented (`fix/kekulize-charge-aware-k1`, forked from `main` at
`b07ac653f79eeed03d60303e8849bcd71c1f6ce0`).

**Scope:** `atom_must_be_matched` in
`crates/chematic-core/src/kekulization.rs` only. Does not touch
`build_molecule_from_model` (`crates/chematic-perception/src/aromaticity.rs`,
the atom/bond aromatic-flag rebuild loop — that promote-only-never-demote
behavior is a separate, deliberately deferred fix, tracked as "K2"), does not
change `AromaticityAlgorithm`'s default, does not add a new aromaticity
model, and does not touch `chematic-smiles`.

## Background

PR #134's diagnosis (`docs/rfcs/aromaticity_rdkit_parity_rfc.md`) found that
`chematic_core::kekulize()` hard-fails for 6 fixture classes — tropylium
cation, imidazolium, pyridinium, pyrylium, tellurophene, phosphole — because
`atom_must_be_matched`'s lone-pair-donor rules are charge-blind (and Te(52)
is missing from the donor list entirely). This document is the RDKit
source audit that RFC section 6 called for before implementing the fix.

## RDKit's actual algorithm

RDKit's kekulization (`Code/GraphMol/Kekulize.cpp`, pinned commit
`8afba32ec539dcb2369bc84549d802aca3f7eb39`, fetched directly from
`https://raw.githubusercontent.com/rdkit/rdkit/8afba32ec539dcb2369bc84549d802aca3f7eb39/Code/GraphMol/Kekulize.cpp`)
does **not** classify atoms by a hardcoded per-element donor/acceptor table
the way chematic's `atom_must_be_matched` does. Instead, `markDbondCands`
computes a per-atom **target valence** and checks whether the atom is
exactly one bond order short of it. The relevant block (lines 153–213):

```cpp
sbo += at->getTotalNumHs();
auto dv =
    PeriodicTable::getTable()->getDefaultValence(at->getAtomicNum());
auto chrg = at->getFormalCharge();
if (isEarlyAtom(at->getAtomicNum())) {
  chrg = -chrg;  // fix for GitHub #65
}
// special case for carbon - see GitHub #539
if (at->getAtomicNum() == 6 && chrg > 0) {
  chrg = -chrg;
}
dv += chrg;
int tbo = at->getTotalValence();
int nRadicals = at->getNumRadicalElectrons();
...
if (totalDegree + nRadicals >= dv) {
  // if our degree + nRadicals exceeds the default valence,
  // there's no way we can take a double bond, just continue.
  continue;
}

// we're a candidate if our total current bond order + nRadicals + 1
// matches the valence state
if (dv == (sbo + 1 + nRadicals)) {
  dBndCands[allAtm] = 1;
} else if (!nRadicals && at->getNoImplicit() && dv == (sbo + 2)) {
  dBndCands[allAtm] = 1;
}
```

Here `sbo` is the atom's current bond-order sum (each aromatic ring bond
counts as 1, plus `getTotalNumHs()`). An atom becomes a double-bond
candidate — chematic's `must_match = true` — iff `dv == sbo + 1`: it is
exactly one bond order short of its target valence `dv`, and a double bond
(which adds +1 over a single bond) would close that gap exactly.

Two mechanisms make this formula charge-aware in *opposite* directions
depending on element:

1. **For most elements (including N, O, P, S, Se, Te — none of these are
   "early" atoms):** `chrg` is added to `dv` directly, unflipped.
   `isEarlyAtom` (`Code/GraphMol/Atom.cpp`, "Determine whether or not an
   element is to the left of carbon") is a hardcoded 119-entry table; the
   entries actually fetched and checked for this audit:
   `false, // #6 C`, `false, // #7 N`, `false, // #8 O`, `false, // #14 Si`,
   `false, // #15 P`, `false, // #16 S`, `false, // #34 Se`,
   `false, // #52 Te` — every element this fix touches reads `false`. A
   **positive** charge therefore *raises* `dv`, meaning the atom needs
   *more* bond order to reach its target — which a double bond can supply.
2. **For carbon specifically**, an explicit second special case
   (`Kekulize.cpp` lines 164–166, their own comment: `// special case for
   carbon - see GitHub #539`) flips the sign *again* when
   `atomicNum == 6 && chrg > 0`: a cationic carbon's `dv` *drops* instead of
   rising, because a carbocation is electron-deficient (empty p-orbital, one
   fewer bond needed), not electron-rich like a protonated heteroatom.

This is the exact, source-cited explanation for why RFC root causes A and B
require **opposite-direction** fixes: heteroatom cations (N+, O+) must be
*added* to the must-match set on protonation, while carbon cations must be
*removed* from it.

## Empirical verification against the pinned RDKit build

`rdkit==2026.03.3` is installed in this repo's `.venv` and pinned to the same
commit. Rather than trust a hand-traced simulation of the C++ above, every
claim in this document was independently confirmed by running
`Chem.MolFromSmiles` + `Chem.Kekulize(mol, clearAromaticFlags=True)` and
inspecting which atoms end up incident to a `DOUBLE` bond (i.e., RDKit's own
real `dBndCands` verdict, observed rather than re-derived):

| Fixture | SMILES | Atom | charge | numH | `has_double_bond` (RDKit) | chematic verdict needed |
|---|---|---|---|---|---|---|
| tropylium_cation | `c1ccc[cH+]cc1` | C (idx 4) | +1 | 1 | **False** | must NOT match (Root Cause B) |
| imidazolium | `c1c[nH+]c[nH]1` | N (idx 2) | +1 | 1 | **True** | must match (Root Cause A) |
| imidazolium | `c1c[nH+]c[nH]1` | N (idx 4) | 0 | 1 | False | must NOT match (existing rule, unchanged) |
| pyridinium | `c1cc[nH+]cc1` | N (idx 3) | +1 | 1 | **True** | must match (Root Cause A) |
| pyrylium | `c1cc[o+]cc1` | O (idx 3) | +1 | 0 | **True** | must match (Root Cause A) |
| tellurophene | `c1cc[te]c1` | Te (idx 3) | 0 | 0 | False | must NOT match (Root Cause C: Te missing from donor list) |
| phosphole | `c1cc[pH]c1` | P (idx 3) | 0 | 1 | False | must NOT match (Root Cause D, kekulize-layer only: P missing entirely) |
| pyridine (regression) | `c1ccncc1` | N (idx 3) | 0 | 0 | True | must match (unchanged) |
| pyrrole (regression) | `c1cc[nH]c1` | N (idx 3) | 0 | 1 | False | must NOT match (unchanged) |
| cyclopentadienyl_anion (regression) | `c1cc[cH-]c1` | C (idx 3) | -1 | 1 | False | must NOT match (unchanged, existing anion rule) |
| boron_azine (regression) | `b1ccccn1` | B (idx 0) | 0 | 0 | True | must match (unchanged, existing boron rule) |

(Reproduction script: parse each SMILES with RDKit, kekulize, and check
`any(b.GetBondType() == Chem.BondType.DOUBLE for b in atom.GetBonds())` per
atom — run against the pinned `2026.03.3` build in this repo's `.venv`.)

## The fix

`atom_must_be_matched` now:

- Requires `atom.charge <= 0` for the O/S/Se/Te lone-pair-donor exemption
  (was unconditional), and adds Te(52) to that element list (Root Causes A
  + C). A charged chalcogen (`[o+]`) falls through to the final `_ => true`
  catch-all, matching RDKit's `dv += chrg` (unflipped for O/S/Se/Te).
- Requires `atom.charge <= 0` for the N-with-H exemption, and extends the
  same rule (plus the anion rule and the neutral-substituent rule) to
  P(15) (Root Causes A + D-at-the-kekulize-layer — P's Hückel/aromaticity
  support in `chematic-perception` remains a separate, already-documented,
  opt-in-shaped sprint; this fix only makes `kekulize()` itself handle P
  correctly, since `kekulize()` has no dependency on the Hückel model).
- Adds a new `6 if atom.charge > 0 => false` arm for cationic carbon (Root
  Cause B), mirroring RDKit's own carbon-specific sign flip.

Every one of the 6 previously-failing fixtures now kekulizes successfully,
and — checked directly against `scripts/aromaticity_rdkit_parity_diagnosis.py`'s
own evidence fields, not assumed — every one produces a bond-by-bond Kekulé
structure **byte-identical** to RDKit's own choice
(`kekule_bond_mismatch_pairs == []` for all 6).
