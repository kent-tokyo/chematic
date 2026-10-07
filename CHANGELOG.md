# Changelog

This file keeps concise, user-visible `chematic` release summaries. Detailed release
history through v1.0.25 is preserved in the
[changelog archive](docs/archive/changelog-through-v1.0.25.md). Validation and
benchmark claims remain scoped to their dated records.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and public releases follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Published v1.0.36 reruns (PyPI Linux wheel, sdist, npm, crate) give the
source's results on the BioTransformer corpus, the 83 reaction fixtures,
xsmarts-autoconf, the 310k SMARTS gate, CIP, hybridization, MMFF typing, the
MOL writer and the WASM Node tests. Record:
`benchmarks/2026-10-06-v1036-published-reruns-and-followups.md`.

- 2D layout: the clash relief stops counting a move once it cannot win
  (identical drawings, `write_mol` 652M → 589M instructions on 2,000 rows),
  so components of more than 60 atoms now also get their ±45° forks
  relieved: narrow branch atoms 1,757 → 916 of 15,000 rows; three-atom
  groups take the wide turns. Record:
  `benchmarks/2026-10-08-734-754-followups-batch22.md`.
- Python: the parsers (`from_smiles`, `from_inchi`, `from_mol_block`,
  `from_mol_block_with_coords`, `parse_mmcif`, SMARTS patterns) raise
  `chematic.ChematicInputError`, a `ValueError` subclass with `category`
  (`malformed` / `resource_limit`), `code` (`smiles_parse`, `input_too_large`,
  …) and `format`. `except ValueError` keeps working. Nucleic-acid edit
  sequences are replayed with the same envelopes in Rust, Python and WASM.
- 2D layout: porphyrins are drawn from a template (pentagons round the
  macrocycle, a central metal at the centre); perfluoroalkyl and other chains
  run straight through CF2/CMe2 atoms; a ring-fusion atom's substituent goes
  into the exterior gap (an angle comparison past 2π sent it into a ring).
  As drawn, 232 of 15,000 rows have a crossing (was 242) and 23 a clash.
  CI: the published-wheel chemistry gates record chematic's outputs on macOS
  and Windows without RDKit and compare them on Linux against a hash-pinned
  RDKit wheel (the first runs failed on a missing macOS x86-64 RDKit wheel
  and on wrong rule-table paths). Published v1.0.37: Linux wheel, npm, crate
  and sdist give the gate counts (`validation/published-wheel-chemistry-gates-expected-v1.0.37.json`).
  A6 harness: a lowest-of-10 stereo-safe MMFF94 arm. Record:
  `benchmarks/2026-10-08-734-754-followups-batch20.md`.
- 2D layout: adamantane-type cages (adamantane, hexamine) are drawn as
  RDKit's projection instead of with a clash; a one- or two-atom substituent
  may swing up to 150° round its hinge in the clash relief (a tropane's
  N-methyl, a pinane's methyls go into a free face); crowded forks open to
  90° where that does as well; the relief no longer stacks two atoms. As
  drawn, 242 of 15,000 rows have a crossing (was 277; RDKit 402) and 25 a
  clash (was 28); 14,757 are clean. SVG text and canonical SMILES are
  faster (identical output). Record:
  `benchmarks/2026-10-08-734-754-followups-batch19.md`.
- 2D layout: a ring atom with two ring bonds and two substituents draws
  them 30° either side of the ring's outward bisector (phenytoin and
  4,4-disubstituted glutarimides drew a crossing); as drawn, 277 of 15,000
  rows have a crossing (was 300; RDKit 402). CI: the binding surface
  inventory and the MMFF94 source-wheel gate are refreshed for this branch
  (the gate pinned v1.0.36's 81 differing heavy atoms). Record:
  `benchmarks/2026-10-07-734-754-followups-batch18.md`.
- 2D layout of crowded molecules: the layout's clash relief may turn a
  branch about an atom whose double bonds are all in rings (a Kekulé-written
  aromatic ring, a ring P=N), and fused or spiro systems of three or more
  rings that overlap are redrawn by stress majorization. As drawn, 26 of
  15,000 exposed-10k and ChEMBL-5k rows have a clash (batch 16: 53; RDKit
  276), 14,699 are clean (RDKit 14,504). `write_mol` layout work is
  cheaper (byte-identical: flat grids, fewer square roots and allocations).
  Record: `benchmarks/2026-10-07-734-754-followups-batch17.md`.
- 2D layout: a branch point's three bonds are drawn 120° apart (they were
  150°/150°/60°, so an ester's C=O ran 60° from its C-C and the two aryl
  rings of a tetrasubstituted alkene such as tamoxifen overlapped); where
  that crowds a component, the narrower fork is kept. A stereocentre's
  terminal neighbour is no longer moved onto another atom when its wedge is
  placed. As drawn, 53 of 15,000 exposed-10k and ChEMBL-5k rows have a clash
  (was 155; RDKit 276) and 14,669 are clean (was 14,541; RDKit 14,504).
  `scripts/layout_2d_quality_vs_rdkit.py` counts narrow branch points.
  Record: `benchmarks/2026-10-07-734-754-followups-batch16.md`.
- Python `to_mol_block()` and `to_mol_block_with_report()` write 2D
  coordinates for every molecule (a molecule without stereo was written with
  every atom at the origin), as the WASM `to_mol_block` and RDKit do.
- The reaction SVG (Python `reaction_svg`, WASM `depict_reaction_svg`),
  similarity maps and the WASM depiction preflight draw the declared stereo
  too. Rust: `chematic_depict::depict_reaction_svg_prepared`,
  `similarity_map_options`.
- Canonical SMILES: an E/Z mark pins only its double bond and that bond's
  substituents for the symmetry search (it pinned the substituents'
  neighbours too, so a tert-butyl's methyls on an imine were branched on):
  explicit-H exposed-10k rows 4111 and 4384 take 2.3 ms and 0.4 ms (about
  60 and 40 ms), strings unchanged on all 30,305 audit inputs.
- Depictions draw the declared stereo: Python `svg()`, `depict_data()`,
  `to_eps()`, `to_pdf()`, highlighted/alert SVGs, `depict_grid()` and reports,
  and the WASM SVG, depict-data and grid functions draw E/Z double bonds with
  their declared geometry and one wedge or hash per stereocentre (the MOL
  writer's depiction), and no longer draw SMILES `/` `\` marks as wedges.
  Read back by RDKit, the drawings carry the input stereo on 1,686/1,687
  exposed-10k and 1,670/1,670 ChEMBL-5k stereo rows (PyPI v1.0.36:
  286 and 317). Rust: `chematic_mol::stereo_depiction::depiction_with_stereo`
  and `with_stereo_depiction`; `chematic_depict::depict_svg_grid_with_layouts`,
  `svg_to_pdf`, `relieve_layout_clashes`. Record:
  `benchmarks/2026-10-07-734-754-followups-batch15.md`.
- Python `SDWriter` writes Å coordinates with stereo (it wrote depiction
  units, 40 Å bonds with y flipped, and E/Z as "either").
- 2D layout: batch 14's clash relief made `write_mol` up to 20× slower; it is
  now 1.8× batch 13 in instructions. Bridged ring systems that still
  overlap are redrawn by stress majorization (bridged rows with a clash in
  the layout: 212 → 5 of 531; RDKit 102), and the clash relief runs again
  after the MOL writer's E/Z reflections. Two collinear bonds no longer
  count as crossing.
- 3D pipeline: a minimization stopped by a stereo constraint at its first
  steps re-embeds from up to three other seeds (A6 row 0161 converges; A6
  264/265 at the default budget).
- 3D pipeline: the MMFF94 iteration budget defaults to 1,000 (Rust
  `PipelineV2Config::minimal()` was 300; Python `safe`/`stereo_safe`
  `force_field_max_iterations` was 200). A6 rows converged within the budget
  263/265 (135 at 200), 265/265 sound, stereo-clean, clash-free; about 20%
  more time. Pass the old value to keep the old behaviour. Record:
  `benchmarks/2026-10-07-734-754-followups-batch14.md`.
- 2D layout (depiction and MOL coordinates): three or more chain branches
  spread evenly, bridged rings (norbornane, tropane, quinuclidine,
  bicyclo[2.2.2]) drawn around their bridge, and acyclic branches mirrored
  or turned to clear clashes and crossings. Exposed 10k + ChEMBL 5k rows
  without a clash or crossing: 10,619 → 14,592 of 15,000 (RDKit
  `Compute2DCoords` 14,504). Coordinates change; RDKit reads the written
  stereo back as before.
- Accurate CIP: a carbon whose rule-5 comparison runs through a lone-pair
  centre (quinuclidinone C6 beside a stereo bridgehead N) gets RDKit's
  label instead of none.
- Performance: canonical SMILES with ring E/Z markers (worst exposed-10k row
  63.5 → 6.3 ms; corpus −40% instructions, identical strings); `write_mol`
  −2.2%; SMARTS `^n` matching −17%; `run_reactants` on Kekulé-written
  reactants −10%.
- New scripts: `scripts/layout_2d_quality_vs_rdkit.py` (clashes and
  crossings against RDKit's layout) and
  `scripts/a6_conformer_quality_vs_rdkit.py` (A6 conformers relaxed and
  ranked by RDKit's MMFF94 against RDKit's best of 10 ETKDGv3 conformers);
  `public_package_3d_chematic.py --force-field-max-iterations`.
- MMFF94 energies match RDKit 2026.03.6 term by term on every
  single-fragment row of the exposed 10k and ChEMBL 5k (14,684): RDKit's
  formal-charge sharing for partial charges (carboxylate/sulfonate heads,
  N5M rings, amidinium/imidazolium N, N=N+ next to diazonium), out-of-plane
  terms on every RDKit-covered trigonal type, no angle-type-0 retry before
  the empirical rule, and RDKit's 100 Å non-bonded threshold. New
  `Mmff94Options { ignore_interfragment_interactions }` (Python/WASM
  `ignore_interfrag_interactions`) leaves out pairs between fragments as
  RDKit's `MMFFGetMoleculeForceField` does; the default keeps them. Record:
  `benchmarks/2026-10-07-734-754-followups-batch13.md`.
- `Mmff94EnergyModel::minimize_bfgs`: a port of RDKit's BFGS minimizer
  (`Mmff94Convergence::Rdkit` follows RDKit's trajectory from the same
  start; `MaxGradient` runs to chematic's absolute residual-force test). The
  3D pipeline's MMFF94 step uses it (A6: 135 of 265 rows converged within
  200 iterations, was 100; 265/265 sound, stereo-clean, clash-free) and
  retries with the earlier L-BFGS when the result crosses declared stereo.
- `find_matches_rdkit_parity(smarts, profile="2026.09.1")` (Rust
  `RdkitRingCountModel::RelevantCycles`): `[R<n>]` counted over relevant
  cycles as in RDKit 2026.09.1 (310,000/310,000 cells of the native
  2026.09.1 grid; no `ring_model_ambiguous` refusals). New
  `chematic_perception::relevant_cycles`.
- Accurate CIP gives RDKit's `r`/`s` to pseudo-asymmetric saturated ring
  lone-pair centres.
- Performance: canonical SMILES of rigid symmetric graphs (a 60-atom random
  cubic graph took seconds, now 3 ms; identical output on 30,305 inputs);
  MMFF94 typing −22% instructions vs v1.0.36 on aromatic input, −0.6% on
  Kekulé input (was +11%).
- Harness: `published_wheel_chemistry_gates.py --resume --corpus-shards N`.
- MMFF94 energies match RDKit 2026.03.6 term by term (within 1e-12 kcal/mol
  on the 262 same-coordinate rows; the largest difference was 0.32): three
  out-of-plane terms per trigonal centre (one was added), RDKit's unrounded
  angle and stretch-bend constants and its linear-angle form, no 10 Å cutoff
  on van der Waals in the default energy and gradient (the opt-in cutoff
  path keeps it), RDKit's torsion lookup through equivalence levels with
  Halgren's empirical rule as fallback (`mmff94_torsion_term_params`), `/`
  `\` bonds read as single for torsion types, and RDKit's stretch-bend
  orientation for equal end types. Minimized geometries change slightly.
  Record: `benchmarks/2026-10-06-734-754-followups-batch12.md`.
- `cip_stereo(mode="accurate")` labels saturated ring lone-pair centres
  (bridgehead amines, cyclic phosphines, sulfonium) as RDKit reads the
  parsed spelling (RDKit inverts `@`/`@@` when one ring-closure digit is on
  the centre); such labels carry `spelling_dependent` (Python),
  `spellingDependent` (WASM) and `cip_label_depends_on_smiles_spelling`
  (Rust), and InChI identity checks refuse them. SMILES parsing is
  unchanged. Exposed 10k: 4,395 agree, no abstention.
- Python: `Mol.to_mol_block`, `to_mol_block_2d` and `to_mol_v3000` warn
  with the new `chematic.StereoLossWarning` when the block loses stereo; the
  returned block is unchanged.
- `chematic_chem::rdkit_canonical_atom_ranks` and
  `chematic_perception::rdkit_canonical_atom_ranks_with_bond_stereo`: RDKit's
  canonical ranks including `STEREOE`/`STEREOZ` bonds (all E/Z rows of the
  three corpora identical to `CanonicalRankAtoms`).
- Performance: canonical SMILES of large explicit-H molecules in reaction
  product atom order (the automorphism check extends from the mapped
  region: 34 ms to 24 µs on a 117-atom product; output unchanged); MMFF94
  typing reads atom adjacency instead of scanning bonds and skips RDKit's
  canonical Kekulé ranking for isolated aromatic rings (−11% instructions
  on aromatic input).
- Harness: `published_wheel_chemistry_gates.py` resolves its path
  arguments; `biotransformer_rule_corpus.py --rules-slice/--reactants-slice`.
  The published v1.0.36 Linux aarch64 wheel gives every Linux x86-64 count.
- MMFF94 typing matches RDKit 2026.03.6 on every heavy atom of the exposed
  10k, ChEMBL 5k and Kekulé-written 15k corpora (81, 5 and 39 atoms differed
  in v1.0.36): the aromaticity pass reads RDKit's canonical Kekulé structure
  (new `chematic_perception::rdkit_canonical_kekule`, a port of RDKit's
  canonical ranking and `kekulizeWorker`; fullerene cages), RDKit's ring
  list, and RDKit's aromatic view of Kekulé-written input; non-aromatic
  nitrogen is a literal port of RDKit's rules (NGD+, NCN+, NM, NSO, imine/azo
  N next to sulfonyl) and imidazolium C2 is CIM+. Record:
  `benchmarks/2026-10-06-734-754-followups-batch11.md`.
- `rdkit_sssr_ring_order` gives RDKit's ring order on all 8,830 ring rows of
  the exposed 10k (8,825): the SMILES parser records ring closures so
  `Molecule::rdkit_bond_order` lists bonds as RDKit's parser creates them,
  and metal bonds RDKit makes dative are skipped.
- RDKit-parity aromaticity takes RDKit's symmetrized rings on cages (a
  fullerene hexagon was missed); `^n` 135,000/135,000 cells.
- MOL writer: a declared ring E/Z bond the first layout cannot draw is drawn
  on the reversed-order layout (bryostatin; ChEMBL 5k 1,670/1,670 read back).
- Docs: `cip_stereo(mode="accurate")` labels E/Z with the accurate engine
  (it agrees with rdCIPLabeler on 780/780 and 678/678 bonds); the default
  legacy mode is unchanged.
- CI: `published-wheel-chemistry-gates.yml` runs the chemistry gates on the
  published macOS and Windows wheels (`scripts/published_wheel_chemistry_gates.py`).
- Dependencies: cc 1.5.1, jsonschema 0.58.4, pyo3 0.29.3, smallvec 1.16.2.
- Performance: the SMIRKS template cache (`PreparedReaction::shared`) is 16
  least-recently-used shards per reading instead of one map cleared whole
  when full (a working set used to be prepared again after every 2,048 new
  templates), and a hit no longer allocates.
- Performance: `write_mol` on stereo molecules (12% fewer instructions,
  same bytes): E/Z perception skips the CIP ranking and small-ring search
  for double bonds without `/`/`\` markers, and the 2D layout reuses the
  molecule's memoized SSSR.
- `scripts/check_release_docs_consistency.py` fails on a `docs/` link that
  leaves the MkDocs tree (the v1.0.36 Pages deploy failed on one).
- Harness: `biotransformer_rule_corpus.py --chematic-python` runs chematic
  under another interpreter (published CPython 3.9 wheels);
  `chematic_chemistry_dump.py` / `compare_chemistry_dump_rdkit.py` do the
  same for CIP, hybridization, MMFF and the MOL writer;
  `rdkitjs_reaction_corpus.mjs` runs the corpus with RDKit.js.

### Added

- Added the versioned `chematic.nucleic-acid.v1` document model for bounded,
  explicit DNA/RNA strand metadata. It preserves ordered residues, standard or
  known modified-base identity, sugar identity, atom ownership, linear
  phosphodiester linkage references, and annotations without inferring a
  sequence or flattening the document into a molecule. Rust, Python, and
  WASM/Node share validation and metadata-edit fixtures plus stable typed error
  categories for ambiguous mapping, unsupported topology, unknown
  modifications, and resource limits (#715).

### Fixed

- Aromaticity perception no longer treats a neutral O/S/Se/Te atom with more
  than two coordinated neighbours (hydrogens included) as a lone-pair donor.
  This corrects the hypercoordinate sulfur cases from #767 for both Kekulé
  and explicit aromatic input while preserving thiophene, furan, pyridine
  N-oxide, and the other reported controls.

## [1.0.36] - 2026-10-06

Remaining #734 work: BioTransformer's public rule tables (983 rules, about
the 1,200 the reporter used) were run through RDKit 2026.03.6 and the RDKit
profile on 400 corpus molecules, as parsed and with explicit hydrogens.
Where RDKit has a product, the source gives the same product sets for
1,051/1,051 implicit-H pairs (v1.0.35: 775) and 6,604/6,652 explicit-H
pairs (v1.0.35: 3,807), and refuses none (v1.0.35: 241 rules refused on
every molecule). Records: `benchmarks/2026-10-05-biotransformer-rule-corpus.md`,
`benchmarks/2026-10-05-biotransformer-corpus-followups.md`,
`benchmarks/2026-10-05-biotransformer-corpus-followups-2.md`,
`benchmarks/2026-10-06-734-754-followups-batch9.md` and
`benchmarks/2026-10-06-734-754-followups-batch10.md`; the apply
model is written down in `docs/smirks-apply-model.md`.

- **Canonical SMILES E/Z fix:** a direction stashed beside a bond and read
  from the bond's atom2 (reaction products, reactions on aromatic E/Z
  carriers) was written mirrored when the writer chose that bond as marker
  carrier. `remove_hydrogens`, atom/bond removal, metal disconnection and
  tautomer atom reordering no longer drop the stash's anchor.
- `remove_hydrogens` keeps an E/Z whose only marker at one end was on a
  removed H atom: the marker moves to the end's other substituent, and the
  markers of conjugated double bonds are rewritten together
  (`chematic_core::ez_markers`).
- **Behaviour change:** canonical SMILES of explicit-H molecules carry E/Z
  markers on heavy atoms where possible, a ring-closure bond included (its
  marker is written at the ring opening). RDKit's `RemoveHs` read 22 of 572
  explicit-H E/Z molecules written by chematic with an E/Z dropped or
  inverted; now none.
- Canonical SMILES of two E/Z double bonds coupled through a shared bond
  (squaraine diimines in Kekulé form) kept both: one was written inverted
  or without a marker, depending on the input spelling. A double bond in a
  ring of fewer than eight atoms no longer ties the marker orientation of
  its neighbours. RDKit-written MOL blocks of such molecules (14 rows of
  two ChEMBL sets) now read back as RDKit's molecule.
- MOL reader: a centre with one wedge or hash takes its parity from the
  angular order of the flat bonds; two flat bonds a few degrees apart with
  different lengths gave the opposite centre (artemisinin peroxide-bridge
  carbons and a morphinan in RDKit's layout, 14 rows). RDKit-written stereo
  blocks read back as RDKit's molecule for 1,678/1,687 exposed-10k and
  1,663/1,670 ChEMBL-5k rows (were 1,666 and 1,647); the remaining 15 of
  the 16 are read by RDKit itself the same way from its own block.
- Kekulization: a neutral ring-fusion N or P with three aromatic bonds is
  tried as a lone-pair donor first, as RDKit does (a pyrido[1,2-a]pyrimidine
  beside `c(=O)`/`c(=S)`/`c(=N)` did not kekulize, and the imine E/Z of six
  BioTransformer products was ranked without the ring's duplicate atoms; it
  now agrees with RDKit's CIPLabeler on 1,949/1,949 double bonds of 1,492
  E/Z molecules); a three-bonded `[s+]` (`[s+]([O-])` of a thiadiazole
  S-oxide) is a donor too. 22 molecules of a 25,000-molecule check set
  kekulize that did not.
- RDKit profile: E/Z of double bonds copied from the reactant follows
  RDKit's stereo atoms (dropped when the stereo atom stays in the product
  but is no longer bonded to its end, e.g. a Baeyer–Villiger ring
  expansion; dropped for a template double bond with template neighbours at
  both ends that spells no stereo); template markers next to carried ones
  no longer conflict (an enyne reduction gave both E and Z).
- RDKit profile: `R<n>` and ring sizes in reactant templates count RDKit's
  symmetrized rings (a quinuclidine N is `R3`); a ring model that cannot be
  settled is a typed `ring_model_unresolved`.
- RDKit profile: an unmapped product atom naming no single element
  (`P([!#1!#6;O,$([O-])])`) is a dummy atom `*`, as in RDKit (it was a
  parse error). SMILES accepts an unbracketed `*`.
- An aromatic product atom keeps no tetrahedral tag (a steroid A-ring
  aromatization was refused as `ambiguous_stereo_bond_order`); an
  explicit-H ring CH beside a rewritten C=C keeps its H count; an indole
  [nH] beside a zero-order bond is a kekulization candidate as in RDKit;
  a zero-order bond is never written as an E/Z carrier. Surrounding
  whitespace in a SMIRKS is ignored.
- Python reaction functions accept reactants of up to 1,000 atoms, explicit
  H atoms included (was 300).
- Reaction products no longer depend on hash-map iteration order: carried
  atoms are added from template atoms in reactant order and carried bonds in
  product atom order (a mixed aromatic product's Kekulé form, and so whether
  it passes RDKit's sanitize, changed with a map's capacity, BTMR1032).
- Accurate CIP labels ring sulfoxides and selenoxides (RDKit reads their
  spellings as chematic does); ring centres with three single bonds
  (bridgehead amines, cyclic phosphines, sulfonium ions) stay
  `LonePairCenter`, since RDKit's reading of their `@`/`@@` changes with the
  spelling while chematic follows OpenSMILES.
- **Behaviour change:** RDKit-model hybridization (`hybridization_per_atom`,
  SMARTS `^n`): a multiple bond conjugates only toward a conjugation
  candidate (the amine N of metal dithiocarbamates is sp3), and the `=O` of
  a neutral Cl/Br/I oxo acid is sp3 as RDKit's clean-up reads it. Exposed
  10k 220,014/220,015 atoms, ChEMBL 5k all 138,655, `^n` 134,999/135,000
  cells (were 220,003, 138,654 and 134,990).
- **Behaviour change:** MMFF94 typing: an aromatic S is thiophene S5 (44)
  only in a five-membered ring, and a two-connected N double bonded to
  neither C nor N is NM (62), as in RDKit. Heavy atoms typed differently
  from RDKit: exposed 10k 89 → 88, ChEMBL 5k 103 → 53.
- **Behaviour change:** accurate CIP labels phosphorus on a ring carrying a
  double bond (cyclophosphazenes) with RDKit's CIPLabeler label instead of
  `OracleUnstable` (ChEMBL 5k: 4,186 labels agree, no abstention; 520/520
  spellings of 12 P molecules). The label flips between the ring's Kekulé
  spellings in both libraries: Python marks it `"kekule_dependent": True`,
  WASM `"kekuleDependent": true`, Rust
  `chematic_chem::cip_label_depends_on_kekule_spelling`; identity
  deduplication still fails closed on it.
- New `chematic_perception::rdkit_sssr_ring_order`: RDKit's ring list, in
  RDKit's order (port of `findSSSR`/`symmetrizeSSSR`; identical lists on
  8,825/8,830 exposed-10k and all 4,929 ChEMBL-5k ring rows).
- **Behaviour change:** MMFF94 aromaticity takes rings in RDKit's order and,
  like RDKit, counts only rings it has accepted toward an exocyclic double
  bond (staurosporine aglycones, porphyrins). Heavy atoms typed differently
  from RDKit: exposed 10k 88 → 81 (78 in three fullerenes), ChEMBL 5k
  53 → 5.
- **Behaviour change:** RDKit-model hybridization counts a bracket atom's
  radical electrons below an octet, as RDKit's `assignRadicals` does (the
  Mg of a magnesium acetylacetonate is SP3; `C[CH2]` both SP3): all 220,015
  exposed-10k atoms agree.
- MOL writer: a centre whose layout leaves two neighbours 5–15 degrees apart
  gets a wedge when no other drawing works (RDKit reads 1,686/1,687
  exposed-10k and 1,669/1,670 ChEMBL-5k stereo rows back, was 1,681 and
  1,663).
- Performance: SMIRKS entry points that take the template as text
  (`run_reactants*`, Python `run_smirks_checked`, WASM checked reactions)
  share prepared templates through a bounded process-wide cache
  (`PreparedReaction::shared`); RDKit's parse clean-up and its ring model
  are memoized per molecule, and the RDKit-parity ring matcher runs only
  where RDKit's rings differ from the SSSR. The BioTransformer corpus loop
  (974 rules × 100 molecules) reacts in 1.15 s instead of 15.7 s
  (`e2b188a8`). A pre-pass of the canonical writer over explicit-H E/Z
  molecules that changed no output (and made writing them 2.5x slower) is
  gone. SMARTS `^n` looks up conjugation only for atoms with fewer than
  four neighbours (`[C^2]` matching: 122M → 66M instructions on 1,500
  molecules), and the product valence check reads Kekulé orders in place
  (13 reaction cases: 1.28x → 1.16x v1.0.34's instructions). Product
  construction reserves the product once and keeps its atom tables as
  vectors (13 reaction cases: 0.94x v1.0.34's instructions, same products;
  re-added hydrogens follow product atom order). `write_mol` writes laid-out
  atom lines without the formatter (5% fewer instructions, same bytes).

- SMIRKS with more `[#6:n]`-style atoms than the aromatic/aliphatic
  expansion enumerates (256 combinations) are applied instead of refused
  (136 BioTransformer rules failed with "atomic-number SMARTS expansion
  exceeds limit"). The RDKit profile no longer enumerates at all: reactant
  atoms are SMARTS and product atoms aliphatic, as in RDKit, whose sanitize
  then perceives aromaticity again.
- **Behaviour change (RDKit profile):** a product `[#6:1]` mapped to a
  reactant `[c:1]`, and a product `[#6;a:1]`, are built aliphatic, as RDKit
  builds them (an aromatic ring cleavage now gives RDKit's product; a
  dioxygenation that leaves one such atom in an aromatic ring gives RDKit's
  empty set). The native profile keeps pairing aromaticity (#679).
- Product bond expressions are read as RDKit's SMARTS parser reads them
  (about 80 rules failed to parse): the first order primitive (`=,:` is
  double, `-;!@` single, `!=` double), and no order after a leading ring
  primitive (`!@-`, `@=`, `@`), which gives a zero-order bond.
  **Behaviour change:** a product `~` keeps the reactant bond between two
  atoms already bonded and is otherwise a zero-order bond that adds nothing
  to valence (`[N+:1]>>[N+:1]~C` gives RDKit's `C[NH3+]~C`; it counted as a
  single bond).
- The RDKit profile accepts `/` and `\` in SMARTS-only reactant templates
  (they match a single or aromatic bond, as RDKit matches them).
- A mapped product atom naming no single element (`[F,Cl,Br,I;-:7]`) keeps
  the matched element and applies a spelled charge, H count or isotope; an
  H count or charge inside a list takes the first alternative written, as
  RDKit does.
- A product whose template writes Kekulé bonds into an aromatic ring is
  judged after RDKit's sanitize steps (kekulize what is still aromatic,
  perceive aromaticity again) instead of being dropped.
- Product kekulization uses RDKit's candidate rule for neutral N and P and
  for `[o+]`/`[s+]`: benzothiophene S-oxidation and methyl hydroxylation of
  an imidazo[1,2-a]pyrimidinone (a ring-fusion N beside `c(=O)`) gave no
  product in either profile.
- E/Z next to the reaction centre is kept when the template rewrites the
  carrier single bond, when the carrier was an explicit H atom, and when the
  direction sits on an aromatic ring bond (`C/N=c1\c(O)c(O)c1`, dropped by
  every template, even `[Cl:1]>>[Cl:1]`).
- An explicit-H ring the template dearomatizes gets RDKit's H counts (ring
  carbons are refilled instead of left as radicals).
- SMARTS `$` is a quadruple bond, as in RDKit. It was a parse error in
  SMARTS, and in a SMIRKS reactant template it matched any bond
  (`[C:1]$[O:2]>>[C:1].[O:2]` split ethanol; found by `xsmarts-autoconf
  tweak`).
- **Behaviour change:** SMILES dummy atoms keep a charge, isotope or H count
  (`[*-]`, `[13*]`, `[*H2]`); they were dropped.
- Evidence: the published v1.0.35 PyPI Linux wheel (CPython 3.9), npm
  package and crates.io crate each give 80 graph/origin/map matches and
  three jointly invalid rows on the 83 reaction fixtures; the PyPI sdist
  gives 75/83 xsmarts-autoconf flags with the same 8 expected differences.
  The source gives the same on both. `tools/published_rust_checked_gate`
  pins the published crate version.
- Scripts: `biotransformer_rule_corpus.py` runs BioTransformer rule tables
  through RDKit and chematic (rules are not vendored; the tables are pinned
  by SHA-256 in the output).

## [1.0.35] - 2026-10-05

SMARTS/SMIRKS dialect fixes from the xsmarts-autoconf report (#734, #754).
On the tool's 83 behaviour flags the source gives RDKit 2026.03.6's value for
75 (published v1.0.34: 65); the remaining 8 are documented policies or an
adapter difference. Record: `benchmarks/2026-10-04-xsmarts-autoconf-v1034.md`.

- **Behaviour change:** SMARTS isotope primitives are enforced by default,
  as in RDKit: `[12C]` matches only a carbon labelled 12, and an unlabelled
  atom has isotope 0 (`[0C]` matches it). `MatchConfig { use_isotopes:
  false, .. }` keeps the old behaviour. MACCS and MCS keep ignoring isotopes.
- SMARTS `v` counts an aromatic atom's Kekulé bond orders (benzene `c` is
  `v4`, pyrrole `n` `v3`), as in RDKit.
- SMARTS `/` and `\` match a single or aromatic bond and do not constrain
  E/Z, as in RDKit (`C/C` matched nothing before, and `F/C=C/F` depended on
  the input spelling).
- The SMARTS parser follows RDKit for chirality classes (`@TH1`–`@TH2`,
  `@AL1`–`@AL2`, `@SP1`–`@SP3`, `@TB1`–`@TB20`, `@OH1`–`@OH30`; other
  numbers are errors) and rejects `i`, `^6`/`^7` and one-letter lowercase
  symbols other than `b c n o p s` (they parsed and matched nothing).
- `run_smirks` / `run_reactants` match Kekulé-written aromatic reactants as
  their aromatic form, as RDKit does after sanitizing its input
  (`[C:1]>>[C:1]` no longer matches `C1=CC=CC=C1`). Atom and bond indices,
  and stereo, are unchanged.
- Explicit-H reactants (`add_hydrogens`): an atom that loses a bond or bond
  order gets its hydrogens back, as RDKit refills them, instead of being
  left a radical (`[C:1][C:2]>>[C:1].[C:2]` on an explicit-H spiro ring,
  `[C:1]=[O:2]>>[C:1]-[O:2]` on acetic acid, `[O:1]>>[O+:1]` on dimethyl
  ether gives `C[OH+]C`).
- An H-free charged aromatic carbon (`[c+]1ccccc1`, `[c-]1ccccc1`) kekulizes,
  as in RDKit; `[c:1]>>[c+:1]` on benzene gives products instead of none.
- SMIRKS product bonds spelled as alternatives (`=,:`) take the first one,
  and product component grouping (`([C:1].[O:2])`) gives one product object,
  as in RDKit (both were parse errors). Under the checked RDKit profile a
  grouped template is `typed_unsupported`.
- A product `:` bond flags its atoms aromatic, as in RDKit; outside a ring
  the product fails sanitize and is dropped (it was returned as `C:C`).
- A mapped atom whose element changes re-derives its H count, as in RDKit
  (`[C:1]>>[N:1]` on `[13CH3]C` gives `CN` for the labelled carbon instead
  of a rejected five-valent N).
- An explicit-H aromatic atom keeps its H count through an edit that leaves
  its degree, element and charge alone (`[n:1]>>[n:1]` on explicit-H
  pyrrole returned nothing).
- `remove_hydrogens` keeps the H atoms RDKit's `RemoveHs` keeps (isolated
  `[H+]`/`[H]`, `[H][H]`, hydrides, H on `*`), turns a removed H on a
  non-organic atom into its H count (`[Na][H]` → `[NaH]`), and keeps each
  atom's total H count: explicit-H pyrrole came back as `c1ccnc1`. On the
  exposed 10k corpus, canonical SMILES after `add_hydrogens` then
  `remove_hydrogens` differ from the input for 45 rows (was 652; the 45 are
  exocyclic imine E/Z spellings).
- **Fixed stereo in `random_smiles` and the non-canonical `write`:**
  `random_smiles` dropped the tetrahedral neighbour order and E/Z bond
  directions when it renumbered atoms, and `write` printed stored `@`/`@@`
  without re-expressing them for its own neighbour order (and dropped the
  tag of an unbracketed stereocentre). On the 1,687 stereo rows of the
  exposed 10k corpus, 1,999 of 5,061 random SMILES (three seeds) encoded a
  different molecule according to RDKit; none do now. Reaction products
  written with `write` are affected the same way.
- **3D (#739):** chematic-3d and chematic-ff take transcendental functions
  from the pure-Rust `libm` crate (as wasm32 builds already did), so
  conformers no longer depend on the host C library. Coordinates change in
  the last bits; on Linux the stereo-safe MMFF94 A6 arm now gives 265/265
  (rows 53 and 246 failed with `FinalStereoViolation`).
- **Behaviour change:** `uncharge` follows RDKit's `Uncharger`: cations
  without H keep their charge, charge-separated groups (nitro, N-oxide) stay,
  and enough anions stay charged to balance them (it set every charge to 0,
  turning nitro groups into radicals). `CanonicalMode::Backbone` keeps the
  old force-neutral key.
- `reionize` deprotonates only the O–H of carboxylic acids and phenols and
  protonates only aliphatic primary/secondary amines (an ether or ester O
  was given a negative charge with two bonds, depending on atom order).
- `canonical_tautomer` / `standardize` no longer turn amides, ureas and
  carbamates into iminols when a methoxyarene is present (216 of the exposed
  10k rows), and drop tetrahedral tags on atoms the shift made sp2 (the
  result depended on input atom order).
- `add_hydrogens` keeps E/Z bond directions held in the side table (aromatic
  exocyclic imines), stereo groups, R-group labels and atom tags;
  `remove_hydrogens` keeps stereo groups, R-group labels and tags. The
  add/remove round trip now leaves canonical SMILES unchanged on all 10,000
  exposed rows.
- MOL V2000/V3000: atoms of type-4 (aromatic) bonds are read as aromatic;
  the writer uses Kekulé bonds when type 4 would lose an aromatic atom's H
  (`c1cc[nH]c1` came back as `c1ccnc1`).
- SMARTS `^n` uses RDKit's hybridization model (conjugation-aware: amide N,
  ester O and carboxylate O⁻ are SP2; plain H atoms have none); per-atom
  agreement with RDKit on the exposed 10k corpus rose from 90.7% to 99.94%.
  `k` alone is "in a ring" and `k0` matches every atom, as in RDKit.
- SMILES bracket atoms read two-digit charges and H counts (`[C+10]`,
  `[CH10]`), as RDKit does.
- Kekulization: an aromatic N/P cation already at valence 4
  (`C[n+]1(C)cccc1`) needs no ring double bond, as in RDKit. Reaction
  products with a plain double bond between two aromatic ring atoms are
  dropped (RDKit's sanitize rejects them).
- Scripts: `reaction_python_checked_candidates.py` and
  `run_published_npm_checked_candidates.mjs` emit checked 83-row candidates
  from an interpreter or npm package that cannot load RDKit, and
  `reaction_python_checked_provenance_gate.py --candidates` classifies them;
  `xsmarts_autoconf_lane.py` (CI, Linux) pins the autoconf commit and fails
  on any flag difference not listed in `validation/xsmarts_autoconf_expected.json`.
  `tools/published_rust_checked_gate` runs the checked RDKit profile of the
  published crate on the 83 rows; `check_published_npm_formula_ez_json.mjs`
  reruns the v1.0.31 WASM formula and E/Z JSON fixes on a published tarball.
  `check_smarts_independent_corpus.py` runs the 31 SMARTS queries on corpus
  rows outside the exposed 10k lane (chematic and RDKit in separate
  interpreters if needed); `check_rdkit_ring_order_dependence.py` lists ring
  queries whose RDKit answer depends on atom order.
- Published v1.0.34 evidence: the PyPI Linux wheel (CPython 3.9), the npm
  package and the crates.io crate each give 80 graph/origin/map matches and
  three jointly invalid rows on the 83 reaction fixtures; npm `formula()`
  spells all 5,000 exposed ChEMBL rows as published Python does (v1.0.30
  npm: 1,575 differed), and E/Z product SMILES parse as JSON.
- **MOL writer stereo:** `write_mol` / `write_mol_v3000` (and Python
  `to_mol_block`, WASM, CLI) write a molecule with stereo on 2D coordinates
  that express it: one wedge or hash per tetrahedral centre, drawn from the
  centre and checked by re-reading it (also at the far end of a wedge that
  meets another centre), and E/Z set by the geometry (a side is reflected;
  in a macrocycle one end of the double bond is moved across its ring
  neighbours). Before, a molecule without coordinates lost every
  tetrahedral centre, and SMILES `/`/`\` markers were written as wedge codes
  only chematic read back. RDKit 2026.03.6 reads 1,681 of the 1,687
  exposed-10k stereo rows and 1,663 of 1,670 ChEMBL-5k rows back as the
  input (v1.0.34: none); no row reads back with an inverted centre. The rest
  lose a centre in a cage layout (6 rows each) or write one ring E/Z bond as
  "either" (stereo 3 / `CFG=2`; 1 row). A centre whose bonds the layout
  draws in nearly one direction gets a second layout (atoms in reverse
  order) or has a two-bond bridge neighbour moved into the widest gap.
  Stereogenic double bonds without declared E/Z are written "either", so
  the drawing does not invent stereo. Kekulé bonds written for an aromatic
  H that type 4 would lose are chosen in canonical atom order. The
  coordinates path wrote the atom-map column one character late (RDKit read
  every atom as mapped 0); fixed.
- **Stereo loss is reported:** `write_mol_with_stereo_report` /
  `write_mol_v3000_with_stereo_report` return the block with a
  `MolStereoLoss` (centres written without a wedge, E/Z bonds written
  "either", square-planar centres, enhanced stereo groups V2000 cannot
  hold). Python: `Mol.to_mol_block_with_report()`, and `strict=True` on
  `to_mol_block` / `to_mol_block_2d` / `to_mol_v3000` raises `ValueError`
  instead of returning a lossy block. WASM: `to_mol_block_strict`,
  `mol_block_stereo_loss_json`; `to_mol_block` now uses the writer's stereo
  layout (it passed plain layout coordinates, so E/Z was not set).
- 2D layout: a fused ring system stopped placing rings after a bound that
  shrank as rings were placed, so the last ring of a pentacyclic system fell
  to a fallback with stretched bonds; a ring fused on two atoms goes to the
  side where it lands on fewer placed atoms. Molecules with a bond outside
  0.9–1.1 bond lengths: 233 of 10,000 exposed rows (v1.0.34: 272), ChEMBL 5k
  242 (297). This changes some SVG depictions.
- `canonical_tautomer` no longer turns peptide amides into enols when a
  molecule has a phenol: a change to the aromatic O/C candidate set must
  stay within the aromatic system and the heteroatoms on it. On the exposed
  10k corpus, rows with more enol groups than RDKit's canonical tautomer
  113 → 4, amide-derived enamines 44 → 0, S ylides 4 → 0; same tautomer as
  RDKit 7,137 → 7,218 of 10,000.
- **MOL reader E/Z:** conjugated and branched double bonds are read from
  coordinates (two double bonds sharing a carrier bond are reconciled by
  flipping one, as RDKit does for polyenes; the reader used to reject
  retinoids, amidines and dienes with a branch), a wedge bond can carry the
  direction of an adjacent double bond, and double bonds in rings of fewer
  than eight atoms get no E/Z. RDKit-written MOL blocks of the exposed-10k
  stereo rows read back as RDKit's molecule for 1,666 / 1,687 rows (v1.0.34:
  1,548); ChEMBL 5k: 1,647 / 1,670 (1,596).
- Canonical SMILES drop a `/`/`\` marker that specifies no E/Z (its double
  bond's other end has no marker, or the bond is in a ring of fewer than
  eight atoms), as RDKit does, and E/Z carrier choice no longer depends on
  a marker next to a non-stereogenic double bond (`=C(C#N)2`, `=CD2`).
  Canonical output is unchanged for 14,992 of 15,000 corpus rows; the 8
  changes move a carrier, and RDKit reads every row's output as the input
  molecule. `assign_ez_bonds` gives no E/Z to a double bond in a ring of
  fewer than eight atoms.
- **Behaviour change:** `hybridization_per_atom` (Python
  `Mol.hybridization_per_atom`) follows RDKit's model on RDKit's aromaticity
  view (amide/enamine N and aryl ether O are sp2; five- and six-orbital
  centres report 0); it used bond orders only. Agreement with RDKit
  `GetHybridization` on the exposed 10k corpus: 220,003 / 220,015 atoms (was
  199,617). `hall_kier_alpha` and `carbon_types` use it.
  `chematic_smarts::rdkit_hybridization` is public.
- MMFF94 atom typing follows RDKit's `setMMFFHeavyAtomType` for aliphatic
  carbon (by total degree; CR4E 30, CO2M 41 and CNN+ 57 were never
  assigned, and a carbon of an MMFF-aromatic seven-ring was typed alkyl),
  nitroso N (46), NSO N (48), N-oxides (67/68), N=N-N / N=N-C amide-like N
  (10), =N= (53) and isonitrile N (61), phosphorus (by degree: thiophosphate
  P was 26), oxonium / pyrylium / water O (49/51/70) and perchlorate Cl
  (77). Heavy atoms typed differently from RDKit on the exposed 10k corpus:
  89 (was 451; 78 of them in one fullerene cage); ChEMBL 5k: 103.
- **CIP (Accurate):** a multiple bond at the stereocentre itself is not
  expanded, as in RDKit's CIPLabeler (a P=O oxygen ranks like the O⁻ of
  `[P+][O-]`); acyclic phosphorus centres are labelled (only P on an
  unsaturated ring, where Kekulé respellings change RDKit's label, stays
  `OracleUnstable`); a lone pair outside rings is a lowest-priority phantom
  ligand (sulfoxides, sulfinamides, sulfoximines, phosphines, selenoxides:
  540 / 540 centres over random atom orders agree with RDKit). Ring
  lone-pair centres (bridgehead amines) stay `LonePairCenter`. Exposed 10k:
  4,394 labels agree with RDKit, 1 abstention (was 5), none differ.
- Native SMIRKS: a stereocentre that loses its implicit H and gains one new
  neighbour keeps its configuration with the new group in the H's place. The
  RDKit profile (`rdkit_compat=True`) gives RDKit's opposite arrangement
  (raw tag copy onto an appended bond).
- Performance (Linux, exposed 10k, time relative to v1.0.34, >1 is
  slower): canonical SMILES 1.06x, add/remove H 1.09x, SMARTS `^n` 1.5x
  (RDKit model), MMFF94 typing 0.96x, UFF 3D 1.0x, reactions 1.07x (Kekulé
  aromatic reactants are aromatized first). **`write_mol` is about 9.5x
  slower** (8.2 → 77.9 ms for 10,000 molecules): v1.0.34 wrote zero
  coordinates and dropped stereo; the source lays out each molecule with
  stereo (≈21 ms), places and checks wedges and E/Z (≈21 ms), and kekulizes
  in canonical order where type 4 would lose an aromatic H (most of ≈20 ms
  for the molecules without stereo). An earlier draft of this change took
  190.1 ms; a molecule rebuild for every wedge check, angles recomputed in a
  sort comparator and the general float formatter were removed. Writing
  given coordinates (`write_mol_with_coords`) skips the layout.
- Evidence: published PyPI v1.0.34 (CPython 3.9) passes the opt-in SMARTS
  310k gate (309,982 exact, 18 typed refusals, none unexpected) with the
  archived oracle (`check_python_smarts_parity_310k.py --archived-oracle`);
  a clean source commit (`b014d970`, empty diff; libm math, new MMFF
  typing) gives 265 / 265 A6 successes, 102 gradient-converged, all sound,
  stereo-clean and clash-free on the external scorer. (An earlier 265 / 265
  record was built from an uncommitted tree.)

## [1.0.34] - 2026-10-04

- Pinned the published RDKit.js 2026.09.1 rebaseline against the v1.0.33
  browser package on the exposed 10,000-row corpus, with separate Morgan,
  CIP, SMARTS and runtime records. The 12 old/new `[R2]`/`[R3]` SMARTS
  differences were reproduced with a same-binary native C++ legacy-ring
  probe. A corrected RDKit Python 2026.03.6 baseline and a separate
  2026.09.1 Python execution plan were recorded; no distributed 2026.09.1
  Python result or v1.0.34 package result is inferred from them.
- `run_smirks_checked(..., rdkit_compat=True)` and the WASM
  `run_reactants_checked(..., true)` follow RDKit 2026.03.6 reaction
  stereochemistry instead of declining tetrahedral reactant templates
  (#734): reactant `@`/`@@` do not restrict matching, and product tags
  follow RDKit's inversion flags (retain, invert, create, remove) and its
  bond orders, including its raw tag copy when a product atom's neighbours
  cannot be traced to the reactant — so the reordered L-alanine spelling
  under the identity template gives D-alanine, as RDKit does. Reactant tags
  RDKit's SMILES parser would drop (no stereocentre) are dropped first.
  SMARTS-only reactant templates may carry `@` under this profile. When a
  product's stereo depends on the order of an atom's ring-closure bonds,
  which RDKit takes from ring-closure numbers that a parsed molecule does
  not keep, the result is `typed_unsupported` /
  `ambiguous_stereo_bond_order`; a match rejected only by the native E/Z
  check gives `ez_reactant_template_semantics`. New Rust API:
  `run_reactants_traced_rdkit_2026_03_6`,
  `PreparedReaction::run_reactants_traced_rdkit_2026_03_6`,
  `RdkitProfileOutcome`.
- Reaction products keep the configuration of stereocentres carried
  through unchanged. v1.0.33 could invert one when the product's bond order
  around it differed from the reactant's: esterifying `N[C@@H](C)C(O)=O`
  with `[C:1](=[O:2])[OH:3]>>[C:1](=[O:2])OC` gave the D-alanine ester.
- A bond between two mapped atoms that the reactant template does not
  spell is kept, as in RDKit (github #1387): `[N:1][C:2][C:3]>>...` on
  cyclopropylamine keeps the ring instead of opening it.
- On the 83 pinned reaction fixtures, Linux and macOS release-profile source
  wheels match graph, atom origins and template maps on 80 (was 77); the
  remaining 3 are invalid inputs refused by both RDKit and chematic. The
  WASM/Node 83-row test also passes. These are source-build checks, not
  published-package results or general SMIRKS parity. The contributor's
  broader BioTransformer rule collection remains untested.

## [1.0.33] - 2026-10-04

- SMIRKS reactant templates are parsed as SMARTS (#734): query features
  such as `[CD1:1]`, `[CX4:1]`, `[C;H3:1]`, `[C,N:1]`, `[!O:1]` and
  `[$([OH]):1]` now match instead of failing with a SMILES parse error.
  Every stereo-free template matches by its SMARTS reading, as in RDKit:
  `[*:1]` is any atom (was carbon only), `[N+0:1]` and `[CH0:1]` keep their
  zero charge / H count, and an implicit bond between aromatic atoms is
  single or aromatic. Templates with `@`/`@@` or `/`/`\` keep the SMILES
  reading for the stereo checks; a SMARTS-only template with stereo is
  refused with a typed error rather than ignoring the stereo.
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
- Reaction products are checked with RDKit's sanitize rules instead of the
  native valence model (#734): RDKit's valence table and isoelectronic
  charge rule (Li/Na/K/Mg/Ca and transition metals unrestricted, `Cl+` like
  S, `P-2` three), its clean-up spellings (neutral nitro `N(=O)=O`, azide
  `N=N#N`, perchlorate `OCl(=O)(=O)=O` are valid), a dative bond counting
  for its acceptor only, aromatic atoms only in rings, and aromatic products
  must kekulize. Neutral four-bonded N (`N(=C)(C)C`), uncharged N-oxide
  spellings (`N(=O)(C)C`) and ipso-substituted aromatic carbons are dropped,
  as RDKit drops them. `[#7+:1]`-style charged
  atomic-number atoms are supported, so `[#7:1]>>[#7+:1]C` on pyridine gives
  aromatic `C[n+]1ccccc1`.
- Implicit hydrogens on electron-poor charged atoms follow the isoelectronic
  rule (#734): N2+ and C+ take valence 3, B- valence 4, so
  `[N:1]>>[N+2:1]` on methylamine gives `C[NH2+2]` (was `C[NH4+2]`). Atoms
  with four or more valence electrons (N+, O-, S+, C-...) are unchanged.
- Mapped product atoms follow RDKit for hydrogens, charge and isotope
  (#734): a template that does not spell the charge keeps the reactant's
  (`[O-:1]>>[O:1]` leaves `C[O-]`, `[N:1]>>[N:1]C` on `C[NH3+]` gives
  `C[NH2+]C`); an explicit `H0`/`+0` applies; an isotope in the template
  applies, none keeps the reactant's. An atom whose template degree is
  unchanged keeps its explicit H count (`[n:1]>>[n:1]` keeps pyrrole's
  `[nH]`), otherwise its hydrogens follow RDKit's valence lists for every
  element (`[C:1]>>[Se:1]` gives `C[SeH]`, `[C:1]>>[Si:1]` `C[SiH3]`).
- An explicit `H0` on a product-only atom uses ordinary valence inference,
  as RDKit does: `[C:1]>>[C:1][CH0]` on methane gives ethane, while `H0`
  on a mapped reactant atom remains an explicit override.
- A dative bond `->` in a SMIRKS or reaction SMARTS is no longer read as a
  reaction arrow (#734); SMARTS agents (`>[O;X2]>`) are ignored instead of
  failing the template, and `ReactionMatch::atom_map_positions` accepts a
  SMIRKS whose product side has atomic-number primitives.
- Reactants with explicit hydrogen atoms (`add_hydrogens`) give the same
  products as their implicit-H forms (#734); the edited atoms keep explicit
  H atoms. Templates that match `[H]` still see them.
- SMARTS gains RDKit's `z`/`Z` (heteroatom / aliphatic heteroatom neighbour
  count), `d` (non-hydrogen degree), ranges such as `[CD{1-2}]`, `[R{1-}]`,
  `[+{1-2}]`, and dative bonds `->`/`<-` (#734). Bare `D`, `X` and `v` mean
  `D1`, `X1`, `v1`; bare `x` and `h` mean "at least one". New
  `AtomPrimitive` and `BondPrimitive` variants are added for these.
- SMARTS parsing fixes found in review (#734): `C<C` and similar inputs are
  parse errors instead of a panic; out-of-range range bounds
  (`[+{120-130}]`, `[D{999-1}]`) match nothing instead of overflowing or
  being misread; `[r0]` and `r` ranges treat acyclic atoms as ring size 0
  and `k` ranges never match them, as RDKit does; counts and charges read
  multiple digits (`[D12]`, `[+10]`; previously `[D12]` was `D1` and
  isotope 2); `[Xe]` parses as xenon.
- On 89 local reaction-rule cases (the issue's table, BioTransformer-style
  phase I/II rules and the review's edge cases) `run_smirks` gives the same
  product sets as RDKit 2026.03.6 after re-canonicalizing with RDKit. On the
  83 pinned fixtures the checked Python gate matches graph, atom origins and
  template maps on 77 (was 76); 3 chiral templates stay typed unsupported
  and 3 invalid inputs are refused. These are source-build results, not
  published-package evidence.

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
