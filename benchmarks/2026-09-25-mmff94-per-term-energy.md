# MMFF94 per-term same-coordinate energy — 2026-09-25 (#637)

This record adds RDKit's per-term MMFF94 energies to the fixed 265-molecule
Tier A+B same-coordinate lane and uses them to locate and fix the remaining
CheMatic residuals. RDKit 2026.03.6 generates one explicit-H conformer per
row. RDKit and source-built CheMatic 1.0.25 then evaluate those identical
coordinates. For each term (bond, angle, stretch-bend, out-of-plane, torsion,
van der Waals, electrostatic), RDKit's energy is obtained with an
`MMFFMolProperties` that has only that term enabled. On every row the seven
values sum to RDKit's total within 2.3e-12 kcal/mol.

The before and after runs use byte-identical coordinates (the checker
compares every coordinate hash).

## Before and after

| Metric | Before (v1.0.25 MMFF94, `049c12e3`) | After (`13d70a2e`) |
|---|---:|---:|
| Median absolute total delta | 0.215678 kcal/mol | 0.004010 kcal/mol |
| p90 absolute total delta | 1.211657 kcal/mol | 0.030229 kcal/mol |
| Maximum absolute total delta | 9.874153 kcal/mol | 0.321211 kcal/mol |
| Rows within 1 kcal/mol | 230/262 | 262/262 |
| Rows within 5 kcal/mol | 260/262 | 262/262 |

Per-term absolute deltas (rows above 0.1 / rows above 1 kcal/mol / maximum):

| Term | Before | After |
|---|---|---|
| Bond | 7 / 0 / 0.225 | 0 / 0 / 0.000 |
| Angle | 35 / 22 / 5.660 | 7 / 0 / 0.358 |
| Stretch-bend | 29 / 13 / 1.745 | 0 / 0 / 0.003 |
| Out-of-plane | 1 / 0 / 0.173 | 1 / 0 / 0.173 |
| Torsion | 28 / 12 / 2.715 | 0 / 0 / 0.048 |
| van der Waals | 138 / 24 / 5.583 | 0 / 0 / 0.080 |
| Electrostatic | 18 / 18 / 12.069 | 0 / 0 / 0.000 |

## Root causes of the former >5 kcal/mol rows

- **Row 166** (`chembl_tier_b_0101`, a catechol hydroxylamine): 9.874 kcal/mol,
  mostly electrostatic (+12.07) partly offset by angle (−5.12).
  - Cause: hydrogen typing. Phenolic O–H must be HOCC (29), not HOR (21). The
    hydroxylamine O–H must be HOR (21), not HOS (33).
  - These types feed the bond-charge increments, vdW, and angle parameters.
- **Row 231** (`chembl_tier_b_0166`, a peptide): 5.528 kcal/mol, almost
  entirely vdW (+5.58).
  - Cause: MMFF omits the R* asymmetry correction for pairs that include a
    hydrogen-bond donor. CheMatic applied it to every pair, for example CR/HNCO
    R* 3.530 vs RDKit 3.276 Å.

## Fixes (each checked against RDKit 2026.03.6 values)

1. **Hydrogen types come from the parent's MMFF type** (MMFFHDEF semantics):
   - HOCC 29 for phenol, enol, and imidic O–H;
   - HOP 24 for O–H on phosphorus;
   - HOS 33 for O–H on sulfur;
   - HOR 21 for N–OH and thioacid O–H;
   - HNSO/HNCS 28 on NSO2, NC%N, and NC=S nitrogen;
   - HNR+ (36) family on NR+, N+=C, NCN+, NGD+, NPD+, and NIM+ nitrogen;
   - HP 71 on phosphorus;
   - HOM 21 for hydroxide.
2. **vdW:** no asymmetry term when either atom is a donor.
3. **Stretch-bend:**
   - lookup is the exact (stretch-bend type, i, j, k) row in either
     orientation, then the periodic-row default (RDKit order);
   - the former retry at stretch-bend type 0 gave, for example, −0.411 instead
     of 0.30 for a type-1 37-37-37 biaryl term;
   - no term is added at a linear central atom.
4. **Out-of-plane:** lookup steps down the equivalence levels before the
   `(0, j, 0, 0)` wildcard. For example, an aryl amide carbonyl `(7, 3, 10, 37)`
   now resolves to 0.116, not 0.130.
5. **Angle:** the equivalence ladder covers types beyond 55 through the
   numeric-type registry. For example, pyridinium N+ (58) resolves
   `(0, 37, 1, 58)` through `(2, 1, 10)`.
6. **Torsion:** no term is added for i-j-k-i in 3-membered rings.

## Atom-type census (10k exposed corpus)

`scripts/mmff94_atom_type_census.py` compares every atom's MMFF94 type with
RDKit on `validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi`. That is
9,774 comparable molecules: 214,990 heavy atoms and 190,737 hydrogens. 22
CheMatic typed refusals and 204 molecules RDKit cannot set up for MMFF are
counted separately. The before census uses the v1.0.25 release wheel
(SHA-256 `03e3ee5f11bd04a1677e71fa3b0f4544eb131b96a4e947f904189e6801720208`).

| | Before (v1.0.25) | After |
|---|---:|---:|
| Hydrogens with a different type | 3,101 | 56 |
| …of which the parent's type agrees | 3,045 | 0 |
| Heavy atoms with a different type | 2,908 | 2,908 |

The remaining 56 hydrogen differences all follow a heavy-atom difference. The
heavy-atom differences are unchanged by this work. They concentrate in
aromaticity perception of Kekulé and charged inputs, for example fullerenes,
tropone-like rings, and pyrrolo-fused systems. They remain an open A6 gate.

## Remaining residuals

- **Angle:** 7 rows are between 0.13 and 0.36 kcal/mol. Every atom type agrees
  on these rows: five trimethoxybenzoyl piperazines, one cinnamoyl piperazine
  amide, and one bis-aminoquinolinium.
- **Out-of-plane:** 0.173 kcal/mol on one aminoquinolinium macrocycle, whose
  atom types also agree.
- **vdW:** at most 0.080 kcal/mol.

These rows are recorded but not assigned a cause here.

## Gradient consistency

On the identical coordinates, CheMatic's analytic gradient was compared with a
central difference of its own total energy.

- Rows 166 and 231, at δ = 1e-5 Å: maximum scaled error 4.25e-8 and 4.11e-8.
- Row 178 (largest remaining residual):
  - at δ = 1e-5 Å, error 1.57e-5 kcal mol⁻¹ Å⁻¹ (scaled 1.57e-6);
  - the error scales with δ²: 1.57e-3 at 1e-4, and 1.73e-7 (scaled 1.06e-7)
    at δ = 1e-6, which is the separate packet;
  - this is finite-difference truncation at the embedded geometry's
    near-linear C–O–C angles (177.5°), not an analytic-gradient error.

## Provenance and boundary

- Before:
  - source `049c12e3`: v1.0.25 plus the #634/#635 fixes and this runner; MMFF94
    code identical to v1.0.25
  - CheMatic wheel SHA-256
    `74e68d4dddf5f959373a838a028e5beba96d1a08fc680ee9927a15ecf6778469`
- After:
  - source `13d70a2e`
  - wheel SHA-256
    `3976e4d8aa176696869aa3fafa40f53cf98e1b9b9926b2771a9bbeab738d91eb`
- Both runs:
  - clean tracked tree
  - Linux x86_64, Python 3.11.15
  - `rdkit==2026.3.6` manylinux wheel SHA-256
    `dcf09347e1a0dc1d7e5fcc34f5da79d1f4bb0e04ecadee669449b82ae2492799`
- Rows: 265/265 terminal; 262 comparable, 2 RDKit embedding failures, and 1
  declared-unsupported probe. No rows are dropped.
- The RDKit coordinates here were generated on Linux x86_64. They differ from
  the 2026-09-23 macOS arm64 packet, so compare totals across the two records
  only at the summary level.

This is source evidence, not a published-package result. It does not cover
conformer quality, minimization convergence or timeouts, stereo preservation,
speed, or heavy-atom typing. Those remain separate A6 gates.

Evidence (checked by `scripts/check_mmff94_current_energy_evidence.py`):

- `validation/results/mmff94-same-explicit-h-energy-per-term-v1.0.25-before-issue637-2026-09-25.json`
  and `.jsonl`
- `validation/results/mmff94-same-explicit-h-energy-per-term-v1.0.25-issue637-2026-09-25.json`
  and `.jsonl`
- `validation/results/mmff94-same-explicit-h-gradient-delta1e-6-v1.0.25-issue637-2026-09-25.json`
  and `.jsonl`
- `validation/results/mmff94-atom-type-census-v1.0.25-before-issue637-2026-09-25.json`
- `validation/results/mmff94-atom-type-census-v1.0.25-issue637-2026-09-25.json`

Commands:

```
git checkout 13d70a2e
maturin build --release -i python3 -m crates/chematic-py/Cargo.toml -o <wheels>
pip install <wheels>/chematic-1.0.25-cp311-cp311-manylinux_2_35_x86_64.whl
python3 scripts/mmff94_same_explicit_h_energy.py \
  --output <rows.jsonl> --summary <summary.json> \
  --rdkit-artifact <rdkit wheel> --schematic-artifact <chematic wheel> \
  --expected-rdkit 2026.03.6 --expected-schematic 1.0.25 \
  --gradient-input-index 166 --gradient-input-index 231 --gradient-input-index 178
python3 scripts/mmff94_same_explicit_h_energy.py ... \
  --gradient-input-index 178 --gradient-delta 1e-6
python3 scripts/mmff94_atom_type_census.py \
  --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi --output <census.json>
```

The before packet repeats these commands at `049c12e3` with the `11a4ea27`
wheel.
