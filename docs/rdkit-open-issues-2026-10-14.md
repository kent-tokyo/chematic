# RDKit open issue implementation audit — 2026-10-14

Inspected 45 open issues returned by GitHub's
`repo:rdkit/rdkit is:issue is:open` search and compared relevant reports with
chematic source at `ffe682c3` plus the existing working-tree changes.
Issue state is a snapshot; an open issue can already have a proposed upstream fix.

| RDKit issue | Source finding | Decision |
| --- | --- | --- |
| [#9719: traversal-order-dependent 1–4 stereo](https://github.com/rdkit/rdkit/issues/9719) | The native legacy bounds predicate reproduced order dependence when called with an arbitrary reversed outer pair; ordinary builder paths were already oriented, so this does not establish a default matrix defect. | [Added a checked four-atom path query](alkene-path-stereo.md) in Rust/Python with endpoint normalization, relative cis/trans/unspecified results, 34,608 permutation/direction decisions and unchanged bounds/embedding defaults. |
| [#9701: pseudoasymmetric rings and duplicate stereo enumeration](https://github.com/rdkit/rdkit/issues/9701) | Native inositol keys were already stable; aminoadamantane enumeration hit the cap because CH2 atoms were counted. | [Corrected native candidate selection](native-stereo-candidate-budget.md): four aminoadamantane assignments, nine inositols, metadata preservation and 60 InChI-set comparisons with RDKit 2026.09.1. |
| [#9629: atom-renumbering changes stereo counts](https://github.com/rdkit/rdkit/issues/9629) | Native baseline counted one center but lacked a public permutation API. | [Added checked atom permutations](atom-renumbering-stereo.md) with remapped stereo metadata and cold caches; 838 RDKit-reference permutations, including the identity-ring regression. |
| [#9675: query H folding loses/inverts tetrahedral configuration](https://github.com/rdkit/rdkit/issues/9675) | No query-H folding API existed; the main native matcher compared raw chirality tags. | [Local Rust/Python candidate](query-hydrogen-folding.md): conservative H folding, complete-frame parity, lexical ring slots and 198 live RDKit-reference cases; recursive folding is rejected. |
| [#7981: enrichment factor returns an empty list](https://github.com/rdkit/rdkit/issues/7981) | No ranked screening enrichment API existed. RDKit 2026.09.1 reproduced the top-90% and full-library defects. | [Implemented native Rust/Python EF](screening-enrichment.md): ceiling cutoffs, full-library endpoints, order/duplicates and explicit undefined-input errors; 3,048 ordinary live-RDKit comparisons. |
| [#9647: oversized dense matrix allocation](https://github.com/rdkit/rdkit/issues/9647) | Existing RDKit-compatible distance APIs used unchecked sizes and infallible buffers. | [Implemented allocation boundary](rdkit-distance-matrix-limits.md): checked dimensions, 1024-input/prepared-atom cap, fallible buffers and Python errors; applies to the two named distance APIs. |
| [#6716: MOL2 UNITY formal charge](https://github.com/rdkit/rdkit/issues/6716) | The editable reader understood `charge`, but the immutable source projection ignored it; its graph writer also omitted nonzero charges. | Implement the missing projection and graph round trip in this change. |
| [#9458: sulfur resonance enumeration](https://github.com/rdkit/rdkit/issues/9458) | A bounded local N/O/S–C–N/O/S cation-shift API was added in a follow-up. | [Implemented scope and source validation](local-cation-resonance.md): both oxygen/sulfur pairs, indexed-state deduplication, charge/H conservation, limits and explicit status; this is not general resonance enumeration. |
| [#9688: EEM charge conservation](https://github.com/rdkit/rdkit/issues/9688) | `eem_charges` and ionic conservation tests already exist. | Keep the regression; no duplicate implementation. |
| [#9678: allene stereo](https://github.com/rdkit/rdkit/issues/9678) | Allene stereo round-trip support already exists. | Keep the regression; no duplicate implementation. |
| [#9691: reaction-template order and E/Z](https://github.com/rdkit/rdkit/issues/9691) | Product stereo cleanup and an order-invariance regression already exist. | Keep the regression; no duplicate implementation. |
| [#9705: V3000 bond stereo collections](https://github.com/rdkit/rdkit/issues/9705) | Bond stereo-group storage and STEB collection support already exist. | [Added checked Rust/Python group construction and inspection](checked-stereo-groups.md), including mixed atom/bond membership and 216 V3000/geometry checks. Existing STE/STEB parsing/writing is reused. [Added checked full alternatives](stereo-group-inversion.md) for relative tetrahedral/geometric groups, with 528 RDKit-reference comparisons. [Added bounded correlated alternatives](correlated-stereo-group-alternatives.md) across independent AND/OR groups, with 1,056 direct-configuration reference comparisons and explicit mask ordering. [Added conservative meso-group reduction](meso-stereo-group-reduction.md) with group/label-preserving graph proofs in every independent context: 672 decisions, 768 verified witnesses and 672 V3000 round trips. General canonicalization, axes and CXSMILES extensions remain outside these candidates. |
| [#9343: internal-coordinate bounds](https://github.com/rdkit/rdkit/issues/9343) | Implementation and Python tests are present as pre-existing uncommitted work. | Preserve that work; not part of this change or its validation claim. |
| [#9247: recursive SMARTS stereo](https://github.com/rdkit/rdkit/issues/9247) | Matcher changes and Python tests are present as pre-existing uncommitted work. | Preserve that work; not part of this change or its validation claim. |
| [#9524: SA-score smoothing discontinuity](https://github.com/rdkit/rdkit/issues/9524) | chematic uses its own fragment table and score mapping, without the reported logarithmic smoothing function. | Not the same defect. SA-score calibration would need an independent evaluation. |

## Implemented MOL2 contract

- Explicit `UNITY_ATOM_ATTR` `charge` takes precedence, including zero.
  Without it, `N.4` supplies +1 in both native readers.
- Attributes refer to source atom IDs, including nonsequential IDs; partial
  charges remain separate and never supply formal charge.
- Unknown referenced atoms in **charge** attributes, invalid integers,
  overflow, missing values, truncated records, duplicate charges/groups, and
  duplicate UNITY sections return `Mol2Error`.
- Unknown attributes remain source metadata. They are not interpreted as
  charges; detached unknown metadata is retained. Immutable record output
  preserves the exact input, including CRLF.
- Graph-only MOL2 output emits nonzero formal charges in UNITY attributes.
  `NO_CHARGES` describes its partial-charge column.
- MOL2-to-CJSON conversion preserves projected formal charges and separately
  verifies partial charges through source-ID correspondence. Conversion reports
  continue to disclose unverified chemistry and omitted metadata.

This is charge projection and preservation, not a new MOL2 valence sanitizer
or general charge-perception algorithm. The RDKit-compatible MOL2 parser under
`chematic_smiles::rdkit` retains its separate compatibility behavior.

## Evidence

The original five-atom ammonium fixture is retained at
`validation/rdkit_issues/formats/ammonium_unity_charge.mol2`.
Locally installed RDKit **2026.09.1** returned `None` with a nitrogen-valence
error when UNITY was present, and parsed the same graph as charge +1 when
UNITY was removed. RDKit 2025.09.3 showed the same behavior.

Before this change, four of the six initial tests failed: negative charge
projection, explicit-zero precedence, invalid attribute rejection, and
MOL2-to-CJSON charge preservation. The original ammonium fixture already
passed chematic through its N.4 fallback; that pre-existing behavior is not
claimed as a new fix.

Validation commands:

```sh
cargo test -p chematic-mol --offline --quiet
cargo clippy -p chematic-mol --all-targets --offline --no-deps -- -D warnings
```

The new regression suite is `crates/chematic-mol/tests/mol2_unity_charge.rs`.
Local validation passed 975 tests (including seven new regressions), with
zero failures or ignored tests. Clippy with warnings denied and formatting
checks passed for the changed Rust files.
These are local source checks; no package publication or upstream RDKit change
is included.

## October 16 follow-up: checked 3D matrix numerics

The allocation guard motivated by [#9647](https://github.com/rdkit/rdkit/issues/9647)
now has an [opt-in numeric extension](rdkit-checked-3d-distance-matrix.md). It
refuses nonfinite coordinates and genuinely unrepresentable distances, while
recovering finite norms when intermediate squares overflow or underflow.
These are locally investigated numeric boundaries; the upstream issue remains
a report of matrix-allocation safety, and its large-input crash was not run.
