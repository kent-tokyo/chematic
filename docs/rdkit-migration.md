# Migrate from RDKit

CheMatic is not a drop-in RDKit replacement. Choose a declared operation
profile, preserve unsupported inputs, and compare results on your own data
before changing a production workflow. The current release is **v1.2.0**;
the latest complete hosted Python chemistry run is **v1.1.7**, the latest
durable checked-in cross-platform packet is **v1.1.3**, and the
broadest published-package RDKit performance packet is **v1.0.42**.

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
| Path enumeration | `Chem.FindAllPathsOfLengthN` | `mol.rdkit_paths_of_length(length)` | Ordered atom/bond indices, rooted/H/shortest options and bounded search. See [scope](rdkit-path-enumeration.md). |
| Atom environments | `Chem.FindAtomEnvironmentOfRadiusN` | `mol.rdkit_atom_environment(radius, rooted_at_atom)` | Ordered model bonds and atom distances; hydrogen/size flags and bounded search. See [scope](rdkit-atom-environment.md). |
| NumPy graph matrices | `Chem.GetAdjacencyMatrix`, `Chem.GetDistanceMatrix` | Same function names on `chematic.rdkit_compat` | Bounded to 1024 input atoms, empty cache prefix and zero adjacency fill; refuses prepared-index mismatch. See [NumPy facade](compat-graph-matrices.md). |
| Graph queries | `Chem.GetAdjacencyMatrix`, `Chem.GetShortestPath` | `mol.rdkit_adjacency_matrix()`, `mol.rdkit_shortest_path(start, end)` | Sanitized RDKit atom order; adjacency emptyVal=0 only; 1024 atom / 16384 bond input caps. See [graph-query scope](rdkit-graph-queries.md). |
| Similarity | `DataStructs.TanimotoSimilarity(a, b)` | `chematic.tanimoto(a, b)` | Use fingerprints from the same declared profile. |
| 2D depiction | `rdMolDraw2D.MolDraw2DSVG` | `mol._repr_svg_()`; Rust `chematic_depict::depict_svg` | SVG supported; PNG requires the optional feature. |
| Conformer transforms | `rdMolTransforms.Set*`, `ComputeCanonicalTransform`, `CanonicalizeConformer` | `mol.set_*`, `mol.compute_canonical_transform(...)`, `mol.canonicalize_conformer(...)`; Rust `chematic_3d::mol_transforms` | chematic returns new coordinates instead of mutating a conformer. Degenerate principal-axis frames are not byte-identical across solvers. |
| Batch conformer RMSD | `rdMolAlign.GetAllConformerBestRMS`, `GetAllConformerBestRMSToRef` | `mol.rdkit_all_conformer_best_rms(...)`, `mol.rdkit_all_conformer_best_rms_to_reference(...)`; Rust and WASM equivalents | RDKit-compatible atom matching and result order. Coordinates remain caller-owned; the batch APIs do not mutate conformers. |
| Molecule properties | `mol.GetProp(key, autoConvert=False, default=...)`, typed `Get*Prop` | Same names on native `Mol` and `chematic.rdkit_compat.Mol` | Omitting `default` raises `KeyError`; an explicit fallback, including `None`, is returned unchanged. Typed setters preserve `int`, unsigned, `float`, and `bool` identity under `autoConvert=True`, including across pickle round trips. |
| Molecule name | `mol.GetName()`, `mol.SetName(name)` | Same names on native and `rdkit_compat` `Mol`/`RWMol` | Uses the conventional `_Name` property. An unset name returns an empty string; SDF writers place it in the MOL header rather than an SD data field. |
| Reactions | `AllChem.ReactionFromSmarts(...).RunReactants(...)` | `chematic.run_smirks(...)` | Bounded template application, not full SMIRKS parity. Reactant templates are SMARTS queries (`X`, `D`, `;`, `,`, `!`, `$()`); product templates follow RDKit's product semantics (query-only features dropped, element changes applied, unspelled charge and H kept from the reactant) and products must pass RDKit's sanitize valence rules. Canonical SMILES spelling differs from RDKit's; compare by InChIKey or re-canonicalize. |
| Tautomers | `rdMolStandardize.TautomerEnumerator` | `mol.rdkit_canonical_tautomer()`, `mol.rdkit_tautomers()` | Named RDKit 2026.03.1-compatible profile; preserve refusals and compare on your corpus. |
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
  paths remain Experimental. v1.0.40 adds a named fixed-seed RDKit-compatible
  ETKDGv3 path: the recorded Linux x86-64 source comparison matches all
  successful RDKit rows on the exposed ChEMBL 5k and RDKit.js 10k corpora.
  This does not prove other seeds, topologies, platforms, ensemble quality,
  universal parameter coverage, published-package parity or RDKit speed.
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
