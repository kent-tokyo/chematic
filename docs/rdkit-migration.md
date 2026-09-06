# RDKit → chematic Feature Support Matrix

Honest, feature-by-feature capability-transfer guidance for teams coming
from RDKit. This page classifies each area as **Supported**, **Partially
supported**, or **Not currently supported**, names the real chematic API
alongside RDKit's, and states known residuals — it is not a marketing
document, and it does not claim general superiority over RDKit.

See also: [`rdkit_cheatsheet.md`](rdkit_cheatsheet.md) (side-by-side API
snippets for common tasks) and [`rdkit-comparison.md`](rdkit-comparison.md)
(prose comparison of chematic vs. RDKit for teams evaluating which library
to use). This page's added value over both: a per-feature-area
Supported/Partial/Not-supported classification with explicit RDKit-API ↔
chematic-API rows, sourced against this repository's own code and
CHANGELOG rather than general impressions.

**What "verified" means on this page**: every chematic function name below
was checked against `crates/chematic-py/python/chematic/__init__.pyi` or
`crates/chematic-py/src/*.rs` directly (grep + read, not memory). RDKit
function names are its well-documented, stable public API (`Chem.*`,
`AllChem.*`, `rdMolDescriptors.*`, ...) — this page does not assert
anything about RDKit's internal implementation it has not verified against
its public interface.

---

## Legend

- **Supported** — chematic has a direct, tested equivalent.
- **Partially supported** — chematic has *an* equivalent, but with a
  narrower scope, a different default, or a documented residual gap.
- **Not currently supported** — no chematic equivalent exists today.

---

## SMILES parsing

**Supported.**

| | RDKit | chematic |
|---|---|---|
| Parse | `Chem.MolFromSmiles(smiles)` | `chematic.from_smiles(smiles)` (Python); `chematic_smiles::parse` (Rust); `parse_smiles` (WASM) |
| Validate only | `Chem.MolFromSmiles(smiles) is not None` | `chematic.is_valid_smiles(smiles)` |
| CXSMILES | `Chem.MolFromSmiles(s, params)` with `Chem.SmilesParserParams` | `chematic.from_cxsmiles(s)` / `chematic_smiles::parse_cxsmiles` |

Migration note: chematic's parser is a from-scratch Rust implementation, not
a port of RDKit's — expect the same molecule for valid SMILES, but do not
assume identical error messages or identical handling of edge-case/
malformed input.

## Canonical SMILES

**Partially supported — known residual, do not treat as a safe dedup key today.**

| | RDKit | chematic |
|---|---|---|
| Canonicalize | `Chem.MolToSmiles(mol, canonical=True)` | `mol.smiles` / `chematic.canonical_smiles(...)` (`chematic_smiles::canonical_smiles`) |

Structural correctness of canonical SMILES round-tripping is measured at
**100% (0/5000)** on a 5,000-mol ChEMBL corpus and **0/33** on a dedicated
acyclic-polyene corpus (tretinoin/β-carotene/lycopene-class molecules),
per README's "Known Limitations" section — an earlier corruption bug (two
independent parser bugs, both ring-closure-specific) is fixed.

The residual gap is narrower but real: `canonical_smiles()` can still emit
two different, individually-valid `/`/`\` spellings for the same E/Z
double-bond system depending on input traversal order, in roughly **1 in
18** stereo-bearing molecules (measured worst-of-10 on the same 5,000-mol
corpus — 275/5000 unstable, confirmed 100% cosmetic, not structural
corruption). **Do not use `canonical_smiles()` as a dedup or cache key**
until this closes; use `apply_aromaticity()`-normalized output as your own
dedup key in the meantime if this matters for your use case. Full figures
and root-cause detail: README.md's "Known Limitations" section.

For a fail-closed identity key, use `canonical_smiles_stable_key()`. It
reparses and re-canonicalizes the candidate and returns no key when the
canonical spelling is not self-stable or when multiple independent E/Z
systems remain coupled. This is a safety boundary, not a claim that the
historical 275/5000 E/Z-only residual has been eliminated. The
recovered Issue #11 corpus and the current diagnostic parameters/results are
pinned in `validation/canonical_original_corpus_manifest.json`.

## SMARTS matching

**Supported.**

| | RDKit | chematic |
|---|---|---|
| Match | `mol.GetSubstructMatch(patt)` / `GetSubstructMatches(patt)` | `mol.has_substructure(smarts)` / `mol.find_matches(smarts)` (Python); `chematic_smarts::find_matches` (Rust); `smarts_match_atoms` (WASM) |
| Validate | `Chem.MolFromSmarts(smarts) is not None` | `chematic.is_valid_smarts(smarts)` |
| Bulk | `[mol.HasSubstructMatch(patt) for mol in mols]` | `chematic.bulk.substructure_search(smarts, smiles_list)` / `chematic.bulk.substructure_match(smarts, mols)` (Rayon-parallel) |

## Fingerprints

**Supported** for the common bit-vector families; **partially supported**
for exact-parity claims against RDKit's implementation (not independently
re-verified on this pass — see the validation and limitations sections for what has
been measured).

| | RDKit | chematic |
|---|---|---|
| ECFP4 (Morgan, radius 2) | `AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)` | `mol.ecfp4()` (bytes) / `mol.ecfp4_numpy()` (NumPy `(2048,)` uint8) |
| ECFP6 | `AllChem.GetMorganFingerprintAsBitVect(mol, 3, nBits=2048)` | `mol.ecfp6()` |
| Chiral ECFP4 | `useChirality=True` | `mol.ecfp4_chiral()` / `mol.rdkit_ecfp_config(include_chirality=True)` |
| FCFP4 | `useFeatures=True` | `mol.fcfp4()` |
| MACCS 166-bit | `MACCSkeys.GenMACCSKeys(mol)` | `mol.maccs()` / `mol.maccs_numpy()` |
| Atom-pair | `Pairs.GetAtomPairFingerprintAsBitVect(mol)` | `mol.atom_pair_fp()` (bytes) |
| Topological torsion | `Torsions.GetTopologicalTorsionFingerprintAsIntVect(mol)` | `mol.torsion_fp()` |
| Layered | `Chem.LayeredFingerprint(mol, layerFlags=0x7F)` | `mol.layered_fp_layers()` — 7-layer list, documented as equivalent to RDKit's `layerFlags=0x7F` |
| MAP4 (not in RDKit core) | — | `mol.map4()` / `mol.map4_numpy()` |
| Tanimoto similarity | `DataStructs.TanimotoSimilarity(fp1, fp2)` | `chematic.tanimoto(fp1, fp2)` |
| Bulk fingerprints | list comprehension over `mols` | `chematic.bulk.ecfp4(smiles_list)` (NumPy `(N, 2048)` uint8, Rayon-parallel) |

Performance note: batch ECFP4 was measured at **~54.7 µs/mol vs. RDKit's
~94.3 µs/mol** (1.7×) on a diverse 5,000-mol ChEMBL corpus, and
**6.76 µs/mol vs. ~44.5 µs/mol** (6.6×) on a small repeated-fixture
10,000-mol batch — both numbers from `docs/benchmark.md` (measured Python
3.13.6, Apple M4, chematic v0.18.0, RDKit 2026.03.4, 2026-08-23). Cite these
exact, already-measured figures if you cite speed at all — do not
extrapolate to other batch sizes or fingerprint types not in that table.

## Descriptors

**Supported** for the commonly-used physicochemical set.

| | RDKit | chematic |
|---|---|---|
| One descriptor | `Descriptors.MolWt(mol)`, `rdMolDescriptors.CalcTPSA(mol)`, etc. | `mol.mw`, `mol.tpsa`, ... (property access) |
| All at once | manual loop over `Descriptors._descList` | `mol.descriptors()` — dict of 70+ descriptor *functions* in one call (190+ individual *values*, since a few functions such as MQN/BCUT2D/autocorr2d return multi-value arrays — see `docs/rdkit-comparison.md`'s descriptor accuracy caveat) |
| Bulk / DataFrame | manual loop + `pd.DataFrame(...)` | `chematic.bulk.descriptors(smiles_list)` / `chematic.descriptors_df(smiles_list)` (Rayon-parallel, returns a list of dicts / DataFrame directly) |

Accuracy vs. RDKit, per `docs/benchmark.md` and README's badge comment
(4,999-mol ChEMBL subset, chematic v0.18.0 vs. RDKit 2026.03.4, measured
2026-08-23): HBA/HBD/ARC **100%**, MW **99.82%** (±0.01 Da — a genuine
corpus-wide check, added this release), TPSA **100% within ±0.1 Å²**
(README's "TPSA edge cases" bullet notes a residual 0.3%/16-molecule gap in
exotic phosphazene/S=N=P chemistry), LogP (Crippen) **100%*** (max Δ =
1.1×10⁻¹³). These are the only descriptor-accuracy figures this page
cites, and only because README/CHANGELOG already document them.

## Aromaticity

**Partially supported — documented, root-caused gap.**

chematic applies Hückel 4n+2 per SSSR ring independently; RDKit uses
fused-ring electron delocalization. Per README's "Aromaticity model"
bullet: aromaticity-flag parity on Kekulized input is measured at **96.3%**
worst-of-10 (5,000-mol ChEMBL); visible differences concentrate in
N-heterocycles (pyridone, quinolone, indolizine) and non-alternant/
bridgehead-heavy structures (azulene, purine). Root cause is documented as
an `aromatic_context` bypass mechanism in the default compatibility path.

The opt-in `AromaticityAlgorithm::RdkitLike` model is covered by a public
regression gate for three fused macrocycle holdouts, purine (9 aromatic atoms),
and azulene (10). The default Hückel model now has a deliberately narrow
all-carbon odd/odd fused-envelope fallback (including azulene), but remains
model-distinct from the broader RDKit-like path; select and record the
RDKit-like model when parity is needed.
This gate is a supported regression boundary, not a claim of universal RDKit
parity: bridgehead-N and other fused/non-alternant topologies remain known
residuals in both model comparisons.

## Stereocenter / CIP assignment

**Partially supported.**

Per README's badge comment: stereocenter count **99.96%** (legacy) /
**98.6%** (new CIP `FindPotentialStereo`-equivalent path); CIP R/S label
**96.30%** vs. modern RDKit `rdCIPLabeler`, **96.83%** vs. legacy RDKit CIP
assignment. Square-planar (`@SP1`/`@SP2`/`@SP3`-equivalent) stereo is read
automatically from MOL/SDF; write is opt-in only via 3 specific `_checked`
functions — see [`format-capabilities.md`](format-capabilities.md#molsdf).

The default assignment remains the fast legacy CIP path. The hierarchical
accurate engine is opt-in through `CipMode::Accurate` (and the named Python
and WASM binding endpoints). It improves agreement on many structures, but
is not a universal guarantee: symmetric cages, tied Rule 4b/pseudoasymmetric
cases, and unsupported or budget-limited structures may remain unresolved.
The current frozen 155-case residual report resolves 140/140 non-phosphorus
cases with no wrong confident labels or regressions. Fifteen phosphorus rows
remain a separate representation-unstable oracle class and are now surfaced
as typed `OracleUnstable` results; unresolved CIP is represented as unresolved,
never as a guessed R/S label.

## Conformer generation

**Partially supported — real gap in generation quality, not just parity measurement.**

| | RDKit | chematic |
|---|---|---|
| Single conformer | `AllChem.EmbedMolecule(mol)` (ETKDGv3) | `mol.generate_3d()` — distance geometry + DREIDING minimization |
| Multiple conformers | `AllChem.EmbedMultipleConfs(mol, numConfs=N)` | `mol.conformer_ensemble(n, rmsd_threshold=0.5)` — RMSD-based pruning |
| Torsion-knowledge-aware pipeline | ETKDGv3's built-in torsion preferences | `mol.embed_pipeline_v2(config)` — opt-in v2 pipeline: torsion-knowledge-aware distance geometry + stereo verification/repair + policy-gated force field |

README's "Use RDKit if" section states this directly: RDKit's ETKDGv3
includes ML-assisted torsion corrections chematic does not have. The
feature-maturity table in README.md marks 3D conformer generation
(distance geometry + MMFF94) as **Experimental**.

The legacy `generate_3d()`/`generate_coords` path is kept for compatibility
and may produce a usable starting geometry without matching ETKDGv3 quality.
Use the opt-in `embed_pipeline_v2` when bounded work, stage evidence, and
typed failure outcomes are required; even a successful pipeline result is not
a claim of conformational-quality or RDKit parity.

**Correction (2026-08-23):** an earlier version of this page said no
RDKit-comparison figure for conformer RMSD/TFD existed in this
repository's docs. That was wrong — `validation/results/mmff94_bci_gap_
227_phase2_report.md` already measures RMSD (mean 1.685 Å) and TFD (mean
0.2228) against RDKit's ETKDGv3+MMFF94 on the project's 265-molecule
corpus (`pipeline_v2_mmff94_strict`, last re-measured v0.17.0, 241/265
success). This page still does not attempt a deeper quantitative
characterization of the conformer-quality gap than that one summary
figure — see the limitations in this guide for the fuller picture,
including a live defect found in the public `Mol.conformer_ensemble()`
API (distinct from the `embed_pipeline_v2` path measured above).

## Force-field optimization

**Partially supported.**

| | RDKit | chematic |
|---|---|---|
| MMFF94 | `AllChem.MMFFOptimizeMolecule(mol)` | `mol.minimize_mmff94(coords)`; energy via `mol.mmff94_total_energy(coords)` / `mol.mmff94_energy_breakdown(coords)` |
| UFF | `AllChem.UFFOptimizeMolecule(mol)` | `mol.minimize_uff(coords)` — per its own docstring, UFF covers all elements including metals, unlike chematic's MMFF94, which is limited |
| DREIDING (not in RDKit core) | — | `mol.minimize_dreiding(coords)` |

Known residual: MMFF94 atom-typing issue #337 — in the archived v0.18.0
development record,
one sub-bug (aryl isothiocyanate cumulated-double-bond CSP carbon) is
fixed; 6 of the original 8 affected molecules remain an honestly-disclosed
residual, root-caused to a genuine RDKit Kekulization/MMFF-aromaticity-
perception artifact for a specific fused, macrocyclic ring topology (a
pyridinium-conjugated exocyclic-amine scaffold) rather than a locally-
statable atom-typing rule gap — 32/6,693 type-mismatched and 56/6,693
charge-mismatched atoms remain on the 264-molecule reference corpus. See
the [archived roadmap and audit notes](archive/README.md)
and the public force-field documentation
for the full writeup; this page does not re-derive it.

The MMFF94 implementation contains all seven energy-term families, but this
describes implemented terms, not complete chemical coverage. Missing typing or
parameters and non-convergent minimization remain valid outcomes, especially
for charged, metal-containing, fused, or otherwise difficult structures.

## Molecular depiction (2D)

**Supported** for SVG; **not currently supported** for PNG-only workflows
without the optional `png` feature.

| | RDKit | chematic |
|---|---|---|
| Single molecule SVG | `Draw.MolToImage(mol)` / `rdMolDraw2D.MolDraw2DSVG` | `mol._repr_svg_()` (Jupyter auto-render) / `chematic-depict::depict_svg(mol)` (Rust) |
| Grid of molecules | `Draw.MolsToGridImage(mols)` | `chematic.depict_grid(mols, cols)` |
| Reaction depiction | `Draw.ReactionToImage(rxn)` | `chematic.reaction_svg(reaction_smiles)` |

## InChI

**Partially supported — two distinct code paths with different accuracy characteristics.**

| | RDKit | chematic |
|---|---|---|
| Standard InChI | `Chem.MolToInchi(mol)` (vendored InChI C library) | `mol.standard_inchi` / `mol.standard_inchikey` — bit-exact, via the vendored InChI C library (v1.07.5), **requires the `native-inchi` Cargo feature** |
| Default (no C dependency) | — | `mol.inchi` / `mol.inchikey` — pure-Rust approximation, not bit-exact |

Migration note: if your pipeline depends on bit-exact standard InChI
without enabling `native-inchi`, chematic's default path is an
approximation, not a drop-in replacement — README's "Use RDKit if" section
states this explicitly as a reason to prefer RDKit.

## Substructure search

**Supported** — see the SMARTS matching row above; this row exists
separately only because RDKit users often look up "substructure search" as
its own topic. Same functions apply: `mol.has_substructure`/
`mol.find_matches` (single), `chematic.bulk.substructure_search`/
`substructure_match` (bulk, Rayon-parallel).

## Reactions / SMIRKS

**Supported** for template-based reaction application; **not currently
supported** for RDKit's full reaction-standardization/validation toolkit
beyond what's listed.

| | RDKit | chematic |
|---|---|---|
| Apply a reaction | `rxn = AllChem.ReactionFromSmarts(smirks); rxn.RunReactants(reactants)` | `chematic.run_smirks(smirks, reactants)` |
| Combinatorial library enumeration | manual loop over `RunReactants` | `chematic.enumerate_library(smirks, fragment_sets)` |
| Reaction fingerprint similarity | (community recipes, not core RDKit) | `chematic.tanimoto_reaction_fp(rxn1, rxn2)` |
| Query a reaction by SMARTS | manual matching | `chematic.query_reaction(reaction_smiles, smarts)` / `chematic.batch_query_reactions(reactions, smarts)` |
| Atom economy / mass balance | manual computation | `chematic.atom_economy(reaction_smiles)` / `chematic.balance_check(reaction_smiles)` |

## Maximum Common Substructure (MCS)

**Supported** with explicit configuration and timeout reporting.

| | RDKit | chematic |
|---|---|---|
| MCS of a set of molecules | `rdFMCS.FindMCS(mols)` | `chematic.find_mcs(mols, ...)` → `Optional[Mol]` |
| Distinguish timeout from no match | result status | `chematic.find_mcs_checked(mols, ...)` → `(Optional[Mol], bool)` |

Migration note: chematic exposes atom/bond comparison modes, timeout,
ring-matching strictness, charge/isotope/chirality matching, and related
options as keyword arguments. The option model and search semantics are not a
drop-in copy of RDKit's `MCSParameters`; verify application-specific results
before switching.

## CIF / materials & simulation formats

**Supported**, broader than RDKit's out-of-the-box coverage for this
category — to our knowledge, RDKit's core Python package does not ship
native mmCIF/PQR/QCSchema/ORCA/Cube/OpenDX/LAMMPS readers (some may be
reachable indirectly via other tools in the wider ecosystem, which this
page does not attempt to catalogue or verify).
See [`format-capabilities.md`](format-capabilities.md) for the full format
matrix, including which of these formats has **zero WASM
exposure** (plain CIF) and which Python bindings are read-only (plain CIF).
This is a real area where chematic covers formats RDKit's core package
does not — stated narrowly, without a general "chematic beats RDKit on
formats" claim.

## Python batch processing

**Supported.**

| | RDKit | chematic |
|---|---|---|
| Bulk parse | `[Chem.MolFromSmiles(s) for s in smiles_list]` | `chematic.bulk.parse(smiles_list)` (Rayon-parallel) |
| Bulk fingerprints / descriptors | manual `for` loop, optionally with `multiprocessing` | `chematic.bulk.ecfp4(...)`, `chematic.bulk.maccs(...)`, `chematic.bulk.descriptors(...)`, `chematic.descriptors_df(...)` — Rayon-parallel inside a single PyO3 call, no `multiprocessing` boilerplate needed |
| Bulk Tanimoto | `DataStructs.BulkTanimotoSimilarity` | `chematic.bulk.tanimoto(...)` / `chematic.bulk.tanimoto_matrix(...)` / `chematic.bulk.tanimoto_search(query, library)` |

## WASM / browser usage

**Not applicable to RDKit's core package** — RDKit has no first-party
Python-style WASM bindings; RDKit.js is a separate community project.
chematic ships `chematic-wasm` directly from the same Rust source as the
Python bindings. Measured 2026-09-06 from the v1.0.8 release candidate (see the
[artifact record](../benchmarks/2026-09-06-wasm-size-v1.0.8.md)): chematic's WASM bundle is
**3.58 MB raw / 1.31 MB gzip**, versus RDKit.js's `RDKit_minimal.wasm` at
**6.91 MB raw** (gzip not independently measured) — about 2.0× smaller on a
raw-to-raw basis. See [`format-capabilities.md`](format-capabilities.md)
for exactly which formats are and are not exposed at the WASM layer (plain
CIF, notably, is not).

---

## What this page deliberately does not claim

- No "full compatibility" claim anywhere on this page — every row is
  scoped to what was actually verified against source.
- No performance claim beyond the specific, already-measured figures in
  `docs/benchmark.md`/`CHANGELOG.md`, cited with their measurement context
  (environment, corpus, version) intact.
- No claim of general superiority over RDKit. Where chematic has a real,
  narrow advantage (WASM bundle size, native mmCIF/PQR/QCSchema/ORCA/Cube/
  OpenDX/LAMMPS I/O), it is stated as exactly that — narrow and specific —
  not generalized.
- Known residuals (canonical-SMILES E/Z direction-normalization gap,
  aromaticity `aromatic_context` gap, MMFF94 issue #337 residual, InChI
  approximation without `native-inchi`) are stated here, not hidden, and
  are sourced from README.md's "Known Limitations" section and
  CHANGELOG.md rather than invented for this page.
