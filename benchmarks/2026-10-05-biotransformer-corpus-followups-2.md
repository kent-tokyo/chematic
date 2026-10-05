# BioTransformer rule corpus: second follow-up batch

Third pass over the BioTransformer public rule corpus of
[the first record](2026-10-05-biotransformer-rule-corpus.md), after
[the first follow-ups](2026-10-05-biotransformer-corpus-followups.md) (same
rules, reactants, RDKit 2026.03.6, harness and host). It closes most of what
that record left open and the canonical-writer, MOL-reader and kekulization
issues found while checking the fixes. Source lane: code revision `e070bd43`
(the evidence commit that follows changes documents and result files only).

## Results

| Pairs where RDKit has a product | `e2b188a8` | `f64fff42` | Source |
|---|---:|---:|---:|
| Implicit H: same product sets | 1,045 / 1,051 | 1,051 | **1,051** |
| Implicit H: chematic refuses | 2 | 0 | 0 |
| Explicit H: same product sets | 6,453 / 6,652 | 6,601 | **6,604** |
| Explicit H: chematic refuses | 95 | 0 | 0 |
| Explicit H: different sets or RDKit truncated | 104 | 51 | 48 |

Run time 651 s.
[Summary](../validation/results/biotransformer-corpus-followups-2-2026-10-05.json),
[rows](../validation/results/biotransformer-corpus-followups-2-rows-2026-10-05.jsonl)
(rule IDs, reactants and products; no SMIRKS).

### Remaining differences (explicit H only)

| Class | Rows | Rules |
|---|---:|---|
| RDKit raw products all fail its sanitize; chematic gives the implicit-H product (policy, #734 item 5) | 185 | 36 |
| RDKit keeps H atoms over-valent (`[SH2]`) | 15 | BTMR0215, two thiourea standardizations |
| RDKit keeps H atoms as radicals (`[CH]`) | 20 | 8 |
| RDKit stopped at 1,000 raw products (H permutations); its set is a subset | 12 | BTMR1063 |
| RDKit cannot kekulize an aromatic ring joined by unspecified bonds to the ring the template wrote in Kekulé form (naphthalene); chematic keeps the products | 1 | BTMR1031 |

The macrolide rows of the previous record (BTMR0368, 0665, 1171) now agree:
the canonical writer puts the ring-closure alkene's marker on the ring bond
(written at the ring opening) instead of an H atom, so RDKit's `RemoveHs`
keeps the E/Z. One raw-unsanitizable row (BTMR1032, a furocoumarin whose
product keeps zero-order bonds) depends on the order in which chematic adds
carried bonds, which follows a hash map: pre-sizing that map during this
batch gave a Kekulé form RDKit cannot sanitize, so the pre-sizing was
dropped; an order-independent product is open. RDKit 2026.09.1 is not on PyPI (2026.3.6 is the latest
Python release), so the corpus was not rerun against it.

## Fixes

| Issue | Evidence | Change |
|---|---|---|
| Explicit-H E/Z written on an H atom at a ring closure | 572 explicit-H E/Z molecules (exposed 10k): RDKit reads 572/572 back unchanged (was 571); 2,288 random explicit-H spellings, 0 differ | an end whose other substituent is H may carry its marker on the ring-closure bond |
| Imine E/Z beside a ring-fusion N (pyrido[1,2-a]pyrimidine) | E/Z of 1,492 molecules of the exposed 10k and a 15,000-row check set against RDKit's CIPLabeler: 1,949/1,949 double bonds agree (6 differed) | kekulization tries neutral three-aromatic-bond N/P as lone-pair donors first, as RDKit; 17 molecules that did not kekulize now do |
| `[s+]([O-])` of a 1,2,5-thiadiazole S-oxide did not kekulize | 5 molecules of the same 25,000 | a three-bonded aromatic `[s+]` is a donor; only the two-bonded cation (thiopyrylium) needs a double bond |
| Kekulé squaraine diimines: one imine written inverted or unmarked | Kekulé random spellings of the 1,492 E/Z molecules: 0 order-variant, 0 E/Z changed; 12 of 25,000 canonical strings change, all such molecules; RDKit reads 24,999 outputs back as the input (1 bridgehead N) | the carrier-respelling pass leaves coupled ends to the joint resolver; a double bond in a ring of fewer than eight atoms no longer ties its neighbours' marker orientation |
| MOL reader: a centre with one wedge/hash and two flat bonds a few degrees apart with different lengths read inverted | RDKit-written stereo blocks: exposed 10k 1,666 → **1,678** / 1,687, ChEMBL 5k 1,647 → **1,663** / 1,670 | parity from the hash's sign and the angular order of the flat bonds (unit bond directions) |

Remaining MOL reader rows: 15 that RDKit reads the same way from its own
block (its 2D layout draws a cis macrocycle diene bond trans: 8 trichothecene
macrolides in the exposed 10k and 7 in ChEMBL 5k; one hydrazone), and one
bridgehead amine N (below).
chematic-written blocks read by RDKit are unchanged (1,681/1,687 and
1,663/1,670).

## CIP abstentions

ChEMBL 5k stereo rows: 20 centres are `OracleUnstable` (cyclophosphazene P)
and 2 `Tied` (an adamantane bridgehead, which RDKit's CIPLabeler leaves
unlabelled too); exposed 10k: one `LonePairCenter`.

- With the phosphorus fence lifted, the accurate engine gives RDKit's
  CIPLabeler label for all 20 P centres on both Kekulé spellings of the P=N
  ring and under 10 random atom orders each (220/220); both engines flip
  the label with the Kekulé spelling. The fence stays: lifting it is a
  policy decision (the label is not a property of the molecule).
- The bridgehead N stays unlabelled: chematic's reading of `[N@]`/`[N@@]`
  with one ring-closure digit on the centre disagrees with RDKit's (22 of
  40 random spellings), while spellings with none or two agree.

## Performance

| Operation (Linux, CPython 3.11) | Before | Source |
|---|---:|---:|
| 974 rules × 100 corpus molecules, `run_smirks_checked(..., rdkit_compat=True)` | 15.7 s (`e2b188a8`) | **1.15 s** |
| canonical SMILES of the products (7,887 / 7,882) | 1.41 s (`e2b188a8`); 2.49 s with the pre-pass below | 1.49 s |
| SMARTS `[C^2]`, 5,000 molecules, Python `smarts_find` (best of 3, two runs) | 32–37 ms (v1.0.34), 42–44 ms (v1.0.35) | 35 ms |
| `[C^2]` matching, 1,500 molecules (instructions) | 121M (`108106bc`) | 66M |
| 13 reaction cases × 100 (`run_reactants` instructions) | 271M (v1.0.34), 346M (`108106bc`) | 317M |

The SMIRKS-text entry points share prepared templates through a bounded
process-wide cache; RDKit's parse clean-up and ring model are memoized per
molecule, and the RDKit-parity ring matcher runs only where RDKit's rings
differ from the SSSR. The first follow-ups added a canonical-writer pre-pass
for explicit-H E/Z molecules that changed no output on the 572 molecules or
the 1,064 E/Z corpus products and made writing them 2.5x slower; it is
removed. `^n` looks up conjugation only for atoms with fewer than four
neighbours, and the product valence check reads Kekulé orders without
copying the product.

## Gates

| Gate | Result |
|---|---|
| 83 reaction fixtures, checked provenance (source wheel, CPython 3.11) | 80 graph/origin/map, 3 jointly invalid |
| xsmarts-autoconf lane (`6ce94c9`) | 75/83, the same 8 expected differences, 0 unexpected |
| Rust workspace tests, Python tests (1,021), clippy `-D warnings` | pass |
