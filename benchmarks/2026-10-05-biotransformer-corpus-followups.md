# BioTransformer rule corpus: follow-ups after the first record

Second pass over the BioTransformer public rule corpus of
[the first record](2026-10-05-biotransformer-rule-corpus.md) (same rules,
reactants, RDKit 2026.03.6, harness and host). It closes the differences
that record left open, except the classified policies. Source lane: code
revision `f64fff42` (the evidence commit that follows changes documents and
result files only).

## Results

| Pairs where RDKit has a product | v1.0.35 | `e2b188a8` | Source |
|---|---:|---:|---:|
| Implicit H: same product sets | 775 / 1,051 | 1,045 | **1,051** |
| Implicit H: chematic refuses | 241 | 2 | 0 |
| Explicit H: same product sets | 3,807 / 6,652 | 6,453 | **6,601** |
| Explicit H: chematic refuses | 2,478 | 95 | 0 |
| Explicit H: different sets or RDKit truncated | 367 | 104 | 51 |

Rules refused on every reactant: 241 → 5 → 0. Run time 858 s.
[Summary](../validation/results/biotransformer-corpus-followups-2026-10-05.json),
[rows](../validation/results/biotransformer-corpus-followups-rows-2026-10-05.jsonl)
(rule IDs, reactants and products; no SMIRKS).

### Remaining differences (explicit H only)

| Class | Rows | Rules |
|---|---:|---|
| RDKit raw products all fail its sanitize; chematic gives the implicit-H product (policy, #734 item 5) | 185 | 36 |
| RDKit keeps H atoms over-valent (`[SH2]`) | 15 | BTMR0215, two thiourea standardizations |
| RDKit keeps H atoms as radicals (`[CH]`) | 20 | 8 |
| RDKit stopped at 1,000 raw products (H permutations); its set is a subset | 12 | BTMR1063 |
| RDKit's `RemoveHs` drops an E/Z whose marker chematic writes on an H atom (a macrolide's ring-closure alkene); chematic's own `remove_hydrogens` keeps it | 3 | BTMR0368, 0665, 1171 |
| RDKit cannot kekulize an aromatic ring joined by unspecified bonds to the ring the template wrote in Kekulé form (naphthalene); chematic keeps the products | 1 | BTMR1031 |

## Fixes

| First record's open item | Rows | Source |
|---|---:|---|
| E/Z lost when an H atom carrying a conjugated double bond's marker is removed | 36 | `remove_hydrogens` moves the marker to the end's other substituent and rewrites conjugated markers together (new `chematic_core::ez_markers`); the canonical writer carries markers on heavy atoms in explicit-H molecules |
| E/Z of double bonds a template creates or copies | 4 | template markers next to carried ones are reconciled; the RDKit profile follows RDKit's stereo atoms (legacy CIP rank) into the product: replaced by the end's other reactant substituent when left out, dropped when still present but detached, dropped for a template double bond with template neighbours at both ends that spells no stereo |
| `reactant_too_large` | 92 | Python reaction functions accept 1,000 atoms, explicit H included (was 300) |
| `R<n>` ring basis (BTMR1099; BTMR0774 `X4R1` on quinuclidine and adamantane) | 5 | RDKit profile matches ring counts and sizes on RDKit's symmetrized rings; typed `ring_model_unresolved` when the ring model cannot be settled |
| unmapped product atom naming no element | 5 rules | RDKit profile builds a dummy atom `*`; bare `*` now parses in SMILES |
| `ambiguous_stereo_bond_order` (BTMR0822) | 1 | an aromatic product atom keeps no tetrahedral tag |
| product valence refusals | 2 | a template-aliphatic atom keeping an aromatic bond keeps its H count (BTMR0254); thiophene S-oxidation (BTMR0747) |
| other product-set differences | 6 | indole [nH] beside a zero-order bond is a kekulization candidate (BTMR1031 on indole); imine E/Z of the two S-oxidation matches agree (canonical writer fix below) |
| product `[C;+10]`, `[C;H15]` | | already parsed; RDKit's results confirmed in tests |

Found on the way:

- The canonical writer mirrored the E/Z of a stashed direction read from the
  bond's atom2: its carrier resolution wrote atom1-relative markers that
  emission read against the stash's anchor. With the markers of 572 E/Z
  molecules (exposed 10k) moved into atom2-anchored stashes, 92 of 272
  checked canonicalized differently; now 0 of 572. `remove_hydrogens`,
  `Molecule::remove_atom`/`remove_bond`, metal disconnection and tautomer
  atom reordering dropped the anchor of such stashes; they keep it now.
- Explicit-H SMILES written by `canonical_smiles(add_hydrogens(m))` for the
  572 E/Z molecules: RDKit 2026.03.6 read 22 back with an E/Z dropped or
  inverted (its `RemoveHs` moves a marker off an H atom onto a conjugated
  bond); now 1 (a macrocycle). chematic reads all 572 back unchanged, and
  random SMILES of the explicit-H forms canonicalize identically (0 of
  1,716 differ; implicit-H stereo molecules 0 of 6,748).
- A product `@`/`~` bond is never chosen as an E/Z marker carrier (it was
  written as a single bond).

## Gates

| Gate | Result |
|---|---|
| 83 reaction fixtures, checked provenance (source wheel, CPython 3.11) | 80 graph/origin/map, 3 jointly invalid |
| xsmarts-autoconf lane (`6ce94c9`) | 75/83, the same 8 expected differences, 0 unexpected |
