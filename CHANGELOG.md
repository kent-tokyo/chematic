# Changelog

This file keeps concise, user-visible `chematic` release summaries. Detailed release
history through v1.0.25 is preserved in the
[changelog archive](docs/archive/changelog-through-v1.0.25.md). Validation and
benchmark claims remain scoped to their dated records.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and public releases follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

- SMIRKS reactant templates are parsed as SMARTS (#734): query features
  such as `[CD1:1]`, `[CX4:1]`, `[C;H3:1]`, `[C,N:1]`, `[!O:1]` and
  `[$([OH]):1]` now match instead of failing with a SMILES parse error.
  Templates that are valid SMILES take the previous path unchanged. A
  SMARTS-only template with `@`/`@@` or `/`/`\` stereo is refused with a
  typed error rather than ignoring the stereo.
- SMIRKS product templates follow RDKit's product semantics: query-only
  features (`X`, `D`, `R`, `$()`, `,`/`!` alternatives) are dropped and the
  element, aromaticity, isotope, chirality, H count (`H3` or `;H3`), charge
  and atom map are applied. A mapped atom whose query names no single
  element keeps the reactant's element; an unmapped one is a typed error.
  `[#N;H<n>]` accepts any single-digit H count.
- SMARTS `[H]`, `[2H]`, `[H+]` and `[H:1]` are hydrogen atoms, per
  Daylight/OpenSMARTS and RDKit; `[*H]`, `[CH]`, `[H1]` and `[C;H]` keep
  the H-count meaning.
- A mapped atom takes a changed product element (#734): `[C:1]>>[N:1]` on
  ethane gives `CN` instead of returning the reactant unchanged. As in
  RDKit, the product element applies when it differs from the reactant
  template's; a `[*:1]` product atom keeps the matched element.
- Reaction products are also checked against RDKit's sanitize valence rules
  (#734): neutral four-bonded N (`N(=C)(C)C`), uncharged N-oxide spellings
  (`N(=O)(C)C`) and ipso-substituted aromatic carbons are dropped, as RDKit
  drops them; aromatic products must kekulize. `[#7+:1]`-style charged
  atomic-number atoms are supported, so `[#7:1]>>[#7+:1]C` on pyridine gives
  aromatic `C[n+]1ccccc1`.
- Implicit hydrogens on electron-poor charged atoms follow the isoelectronic
  rule (#734): N2+ and C+ take valence 3, B- valence 4, so
  `[N:1]>>[N+2:1]` on methylamine gives `C[NH2+2]` (was `C[NH4+2]`). Atoms
  with four or more valence electrons (N+, O-, S+, C-...) are unchanged.
- Reactants with explicit hydrogen atoms (`add_hydrogens`) give the same
  products as their implicit-H forms (#734); the edited atoms keep explicit
  H atoms. Templates that match `[H]` still see them.
- SMARTS gains RDKit's `z`/`Z` (heteroatom / aliphatic heteroatom neighbour
  count), `d` (non-hydrogen degree), ranges such as `[CD{1-2}]`, `[R{1-}]`,
  `[+{1-2}]`, and dative bonds `->`/`<-` (#734). Bare `D`, `X` and `v` mean
  `D1`, `X1`, `v1`. New `AtomPrimitive` and `BondPrimitive` variants are
  added for these.
- On 45 reaction-rule cases (the issue's table plus BioTransformer-style
  phase I/II rules) `run_smirks` gives the same product sets as RDKit
  2026.03.6 on 44; the remaining one is a chiral reactant template, already
  documented as outside the pinned RDKit comparison.

## [1.0.32] - 2026-10-03

- Added atom-origin and product-template-map arrays to the opt-in checked
  reaction API in Python and WASM/Node. Ordinary product molecules and the
  existing reaction API are unchanged. Pinned RDKit 2026.03.6 comparisons on
  Linux and macOS release-profile source wheels classify 76/83 rows as exact
  on graph, origin and map; three are typed unsupported, one a diagnosed
  refusal and three invalid in both. WASM release-mode Node tests pass.
  Published-package verification remains open.
- Added opt-in Python `Mol.find_matches_rdkit_parity()`. Source-wheel CI
  matches 309,982/310,000 pinned SMARTS match sets; the other 18 are typed
  unsupported, not exact matches. Native SMARTS is unchanged and published
  package parity is open.
- Corrected MMFF94 aromatic-state propagation on compact fused rings and
  added typed optimizer termination plus stage-level failure diagnostics.
  Source-wheel atom-type differences are 451 heavy and eight H on 10,000
  molecules, versus 2,908 and 56 on published v1.0.31. This does not clear
  geometry, convergence, conformer-quality or cross-platform gates; Linux
  stereo failures on rows 53/246 remain under #739.
- Added versioned reaction, SMARTS and A6 evidence records and kept source,
  published and historical results distinct. See the
  [benchmark index](benchmarks/README.md) for row-level records.

## [1.0.31] - 2026-10-03

- Fixed WASM `MolHandle.formula()` to use the shared Hill-order formula.
  Published v1.0.30 npm and Python differed in spelling on 1,575/5,000
  exposed ChEMBL rows (element counts agreed); the v1.0.31 package still
  requires a published-artifact rerun before claiming cross-binding parity.
- Fixed WASM `run_reactants` and `enumerate_library_2way` JSON escaping for
  product SMILES containing E/Z backslashes. Rust and Node/WASM regression
  tests cover the source fix; published-package verification is still needed.
- Updated the MCP schema-conformance test dependency to `jsonschema` 0.58.2.
  This is a dev-dependency update, not a runtime chemistry change.

## [1.0.30] - 2026-10-02

- Corrected the named `rdkit_hba`, Python `CalcNumHBA`, and WASM RDKit-profile
  `hba` path for substituted aromatic N. The hash-verified published macOS
  arm64 v1.0.30 wheel matches RDKit 2026.03.6 on 5,000/5,000 exposed ChEMBL
  rows (published v1.0.29: 3,641/5,000). This is HBA-specific evidence.
- Added a reproducible, exposed-cohort v1.0.29 Python-wheel accuracy packet:
  10,000 chemistry rows, 310,000 SMARTS cells, 57 reaction fixtures, and a
  5,000-row operation matrix. The reaction gate now rejects unsanitizable
  oracle products and reports missing/extra product sets; it does not claim
  broad SMIRKS equivalence. The operation timing is diagnostic, not a paired
  speed claim.

## [1.0.29] - 2026-10-02

- SMIRKS product templates treat atomic-number atoms as literals (#679):
  an unmapped `[#6](=[#8])[#6]` now adds an acetyl group exactly like
  `C(=O)C` (organic-subset symbols get implicit hydrogens; `[#14]` stays
  `[Si]`; `;H1` is kept), and a mapped `[#7:1]` keeps the matched reactant
  atom's aromaticity instead of producing an extra product with the flag
  flipped. Bare bracket product atoms such as `[C]` likewise leave the
  hydrogen count to valence rules, as mapped atoms already did; an explicit
  `[CH]`/`[NH2]` still pins it. Reactant-side expansion is unchanged.
- SMARTS bracket charges written with repeated signs (`[++]`, `[--]`,
  `[Fe+++]`) parse as one total charge (+2, -2, +3) per Daylight/OpenSMARTS,
  not as several ±1 primitives joined by AND (#680).

## [1.0.28] - 2026-09-28

- Fixed accurate-mode CIP labels decided by Rule 4b when an embedded
  stereocentre's back-to-root ligand is ranked: that Rule 1a comparison pooled
  each sphere into one multiset instead of exploring branch by branch. The
  RDKit 2026.03.6 rebaseline row 4480 (`CO[C@@H]1[C@@H](N)[C@@H](OC)[C@@H](O)[C@H]1O`,
  atom 3) is now S, adjudicated by hand (#634); no other label changed across
  2,539 stereo-tagged test molecules.
- Accurate-mode CIP reports stereo-tagged centres with three explicit ligands
  (bridgehead amines, sulfoxides) as unresolved with the reason
  `lone_pair_center` (Rust `CipUnresolvedReason::LonePairCenter`, Python
  `"lone_pair_center"`, WASM `"lonePairCenter"`) instead of omitting them
  silently. They still get no R/S label.
- On the pinned, exposed 10,000-row RDKit 2026.03.6 comparison, 9,995 rows
  agree exactly and the other five are explicit abstentions (four
  `oracle_unstable`, one `lone_pair_center`). This is source-level, cohort-bound
  evidence, not complete CIP parity or a published-package speed result.

## [1.0.27] - 2026-09-26

- Added Rust `Molecule::set_tag` / `atom_tag` for caller-managed atom labels
  preserved by molecule clone, core graph edits, reaction apply, fragments,
  and aromaticity perception. Private, lazily allocated tag storage preserves
  existing `Atom` struct literals and equality. Tags do not affect SMILES or
  canonicalization; write/parse requires explicit atom-order remapping. Labels
  are `1..=u16::MAX`, need not be unique, and `None` / `Some(0)` clear them.
  Tag changes invalidate cached molecular views.
- Restored the WASM torsion-scan demo API and added browser smoke coverage for
  its typed response; this does not expand the supported 3D chemistry domain.
- Strengthened the RDKit 2026.03.6 SMARTS residual checker to account for
  every classified query/target cell. The remaining 200 of 310,000 cells are
  documented compatibility boundaries, not exact parity.
- Published a separately pinned v1.0.26-wheel MMFF94 quality packet with all
  265 input rows retained and independent geometry/stereo scoring. This is
  historical evidence for v1.0.26, not a v1.0.27 remeasurement; broader A6
  typing, convergence, and conformer-quality gates remain open.

## [1.0.26] - 2026-09-25

- Improved RDKit-compatible SMARTS matching, aromatic/ring perception, and
  selected fingerprint hot paths. The dated source differential preserves the
  measured output boundary; shared-VM timing is not a package or universal
  performance claim.
- **Behavior change (Python/WASM SMARTS):** public SMARTS APIs now match a
  perceived RDKit-parity aromatic view. A Kekulé benzene matches `c`, and no
  longer matches `[#6]=[#6]`; returned atom indices still refer to the input.
  Rust's default matcher is unchanged.
- Made accurate E/Z CIP output bond-keyed, including bond endpoints, and
  improved same-coordinate MMFF94 per-term agreement with RDKit 2026.03.6.
  These are scoped source-evidence improvements, not complete parity,
  convergence, conformer-quality, or published-package claims.

Detailed inputs, results, and limits are in [validation](docs/validation.md)
and the [benchmark index](benchmarks/README.md).

## [1.0.25] - 2026-09-25

- Added atom-output and source-atom provenance for SMILES, fragments, and Rust
  reaction products.
- Corrected the named RDKit-compatible Python HBA profile and plain SMILES
  aromatic/non-aromatic ring-closure spelling.

## [1.0.24] - 2026-09-24

- Optimized ring perception, RDKit-parity aromatic preparation, and
  SMARTS-existence checks without changing the checked source outputs.

## [1.0.23] - 2026-09-24

- Improved RDKit-defined output agreement for compatible fingerprints, MACCS,
  QED, Murcko scaffolds, and selected descriptors.
- Aligned SMARTS implicit-bond semantics and added `rdkit_tpsa`.

## [1.0.22] - 2026-09-24

- Added derived caches and a versioned RDKit operation matrix with explicit
  equivalence and environment boundaries.

## [1.0.21] - 2026-09-23

- Fixed aromatic-stash E/Z canonical SMILES and added strict documentation-site
  checks.

## Earlier releases

See the [v1.0.0–v1.0.20 archive](docs/archive/changelog-through-v1.0.25.md)
for release summaries and historical details.
