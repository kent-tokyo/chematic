# Migrate from RDKit

CheMatic is not a drop-in RDKit replacement. Choose a declared operation
profile, preserve unsupported inputs, and compare results on your own data
before changing a production workflow. The current release is **v1.0.33**;
the latest published-package comparison packet is pinned to **v1.0.30**.

For runnable calls, see the [RDKit cheat sheet](rdkit_cheatsheet.md). For
operation-level measurement and unsupported states, see the
[compatibility dashboard](compatibility-dashboard.md),
[validation report](validation.md), and [compatibility scope](compatibility-scope.md).

## API map

| Task | RDKit | chematic | Boundary |
|---|---|---|---|
| Parse SMILES | `Chem.MolFromSmiles(s)` | `chematic.from_smiles(s)`; Rust `chematic_smiles::parse` | Invalid inputs and error text need not match. |
| Write SMILES | `Chem.MolToSmiles(m)` | `mol.smiles`; Rust `canonical_smiles` | Text equality is not a general identity test. |
| Match SMARTS | `mol.HasSubstructMatch(q)` | `mol.has_substructure(q)` / `mol.find_matches(q)`; opt-in `mol.find_matches_rdkit_parity(q)` | Native query interpretation has recorded residuals. The opt-in result distinguishes no-match from typed refusal; match atom order can differ. |
| Descriptors | `Descriptors.MolWt(m)`, `rdMolDescriptors.CalcTPSA(m)` | `mol.mw`, `mol.tpsa`, `mol.descriptors()` | Choose native or named RDKit-compatible profiles explicitly. |
| Fingerprints | `rdFingerprintGenerator.GetMorganGenerator(...)` | `mol.ecfp4()` or named RDKit-compatible Morgan API | Native ECFP and RDKit Morgan bit positions are different profiles. |
| Similarity | `DataStructs.TanimotoSimilarity(a, b)` | `chematic.tanimoto(a, b)` | Use fingerprints from the same declared profile. |
| 2D depiction | `rdMolDraw2D.MolDraw2DSVG` | `mol._repr_svg_()`; Rust `chematic_depict::depict_svg` | SVG supported; PNG requires the optional feature. |
| Reactions | `AllChem.ReactionFromSmarts(...).RunReactants(...)` | `chematic.run_smirks(...)` | Bounded template application, not full SMIRKS parity. Reactant templates are SMARTS queries (`X`, `D`, `;`, `,`, `!`, `$()`); product templates follow RDKit's product semantics (query-only features dropped, element changes applied, unspelled charge and H kept from the reactant) and products must pass RDKit's sanitize valence rules. Canonical SMILES spelling differs from RDKit's; compare by InChIKey or re-canonicalize. |
| MCS | `rdFMCS.FindMCS(mols)` | `chematic.find_mcs(...)` / `find_mcs_checked(...)` | Check timeout separately from no match. |
| Batch processing | Python loops / `BulkTanimotoSimilarity` | `chematic.bulk.*`, `descriptors_df(...)` | Preserve original input index and all failure outcomes. |
| Browser | `@rdkit/rdkit` | `@kent-tokyo/chematic` | Compare the published packages on the same operation and host. |

The Python import surface is typed, but similarly named calls may have
different defaults, result shapes, or supported chemistry. The
[language-binding guide](language-bindings.md) maps Rust, Python and WASM
errors and ownership. The [format matrix](format-capabilities.md) lists
read/write coverage rather than assuming every RDKit format is supported.

## Chemistry boundaries to check first

- **Identity and stereo:** `canonical_smiles()` is a representation, not a
  universal deduplication key. Use the fail-closed
  `canonical_smiles_stable_key()` only within its documented domain. The
  accurate CIP mode is opt-in; unresolved centers remain unresolved. On a
  pinned exposed 10,000-row source comparison, 9,995 rows agree and five
  return typed abstentions. This is not complete CIP parity.
- **Aromaticity and SMARTS:** the default per-SSSR model and opt-in
  `AromaticityAlgorithm::RdkitLike` are distinct. The v1.0.30 published
  audit retains 200 disagreements in 310,000 SMARTS query/target cells.
  A locally built source wheel's opt-in Python profile matches 309,982 cells
  and marks 18 as typed unsupported; it is not yet published-package evidence.
  Check `result["status"]` before using `result["matches"]`: `[]` under `ok`
  means no match, while `None` under a refusal is not a negative match.
  Do not assume `R`/`r` ring-count semantics or a Kekulé/aromatic input
  spelling is interchangeable across engines.
- **Fingerprints:** native `ecfp4()` uses its own hashing. Use the named
  RDKit-compatible Morgan API when comparing RDKit bits or rebuilding an
  existing search index; never mix the two profiles in one index.
- **InChI:** the default pure-Rust `mol.inchi`/`mol.inchikey` path is not
  Standard InChI. `mol.standard_inchi`/`standard_inchikey` require the
  optional `native-inchi` feature and its C library.
- **3D:** `generate_3d`, `embed_pipeline_v2`, MMFF94 and related force-field
  paths remain Experimental. A typed successful result does not prove
  ETKDGv3 conformer quality, universal parameter coverage or RDKit speed.
- **Rich formats and reactions:** CDXML presentation editing and
  Markush/polymer expansion are bounded. V3000 coordination/haptic
  semantics are not generally interchangeable. The published v1.0.30
  reaction audit matches the original 57 fixtures but retains differences
  in an 83-case extension; it does not predict yield or selectivity.
  `run_smirks_checked(..., rdkit_compat=True)` reproduces RDKit 2026.03.6
  reaction stereochemistry, including its bond-order-dependent tag copy;
  it returns `typed_unsupported` (`ambiguous_stereo_bond_order`) when the
  answer depends on ring-closure numbers the input molecule does not keep.
  The native `run_smirks` keeps its own stereo semantics (reactant `@`/`@@`
  filter matches; product tags are geometric).

## Move a persisted workflow safely

1. Pin both package versions, operation options, corpus hash and failure
   policy. Do not substitute a source-build result for a published package.
2. Choose a fingerprint and aromaticity/CIP profile. Rebuild persisted
   fingerprints and search indexes; do not import RDKit bit vectors as native
   chematic ECFP.
3. Compare semantic graph, stereo, numeric tolerances and failed/refused rows
   on representative held-out inputs. Do not use canonical-SMILES string
   equality as the only test.
4. In browser code, release `MolHandle` objects with `.free()` and keep
   unsupported or cancelled rows visible. The [browser guide](use-cases/browser-app.md)
   has working examples.

The [v1.0.30 artifact packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-02-v1.0.30-published-artifact-gates.md)
and [benchmark index](https://github.com/kent-tokyo/chematic/tree/main/benchmarks)
provide versions, denominators and raw records. Older measurements remain
historical; none imply blanket parity or superiority over RDKit.
