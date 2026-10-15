# Proof-backed reduction of a redundant relative stereo group

RDKit [#9705](https://github.com/rdkit/rdkit/issues/9705) proposes treating a
mixed relative group's full alternative as the simultaneous inversion of its
atom centers and geometric double bonds. Sometimes that full alternative is
the same drawn compound. The relative ambiguity can then be represented as ABS
without changing the allowed configurations.

This candidate adds conservative reduction of **one selected AND/OR group**.
It builds on [full group inversion](stereo-group-inversion.md) and
[correlated alternatives](correlated-stereo-group-alternatives.md).

```python
import chematic

source = chematic.from_smiles(r"C/C=C/[C@H](C)/C=C\C")
source = source.with_stereo_groups([("and", 1, [3], [1, 5])])
report = source.reduce_meso_stereo_group(0)
assert report["reduced"]
assert report["molecule"].stereo_groups == [("abs", 0, [3], [1, 5])]
assert report["molecule"].smiles == source.smiles
mapping = report["witnesses"][0][1]  # source atom index -> full-alternative atom index
```

Here both double bonds and the center participate in the proof. Putting only
the center in the group gives a different alternative, so that group is
retained. The same API reduces a meso diol group while retaining its true
racemic counterpart.

The term meso here follows #9705's full-alternative redundancy proposal.
Swapping E/Z is part of that group operation; the API does not identify every
IUPAC meso structure or equate a geometric swap with a general mirror operation.

## Exact mapping contract

No canonical-SMILES equality or hash collision can justify a promotion. For
each assignment of all other relative groups, the implementation searches for
a complete bijection from the selected group's drawn state to its full
alternative. A successful witness must preserve:

- atomic element, isotope, charge, aromaticity, resolved H count and degree;
- atom maps and nonchemical atom tags;
- each stored group's atom and bond membership independently, including ABS;
- adjacency and chemical bond order;
- every specified tetrahedral configuration through its complete four-ligand
  reference frame and the permutation's parity;
- every specified geometric double-bond fact, accounting for endpoint reversal
  and choosing the opposite substituent as a reference.

The search prunes only by exact attributes, edges and configured parity. It
uses a bounded recursive graph search with at most 64 atoms. Direction tokens
are interpreted through the existing E/Z fact model, so shared markers and
alternate reference choices do not become raw-token comparisons.

If every other-group assignment has a witness, the selected members join ABS.
If ABS already exists, its members come first and the selected atoms/bonds are
appended; the selected relative entry is removed. Otherwise its stored position
becomes ABS. Other groups retain their kind, number and members. Remaining
stored group indices may shift after removal; no canonical renumbering occurs.

The drawn configuration, indexed atom/bond graph, tags and reference frames
remain unchanged. The output has cold derived caches and cleared stored atom
CIP annotations. Python typed molecule properties are copied independently.
Input metadata/caches are unchanged. Original drawing coordinates can still
be used because this operation changes grouping without changing configuration.

## Results and conservative retention

Python `Mol.reduce_meso_stereo_group(index, *, max_visits=1000000)` returns:

| Field | Meaning |
| --- | --- |
| `molecule` | Independent output molecule |
| `reduced` | True only after all required pairs have verified mappings |
| `checked_pairs` | Configuration pairs examined |
| `search_visits` | Candidate feasibility probes across the whole call |
| `witnesses` | `(lower_mask, atom_permutation)` tuples for successful pairs |

Masks follow the relative-group order used by the family API. Each witness maps
a configuration whose selected bit is zero to the configuration with that bit
set. When a pair has no mapping, the output retains the original groups. The
report can include proofs for earlier pairs; no partial promotion is performed.

**Retained is not a general non-meso classification.** Preserving source labels
and exact group identities can deliberately block an otherwise chemically
valid equivalence. The mapped diol fixture demonstrates this: the unlabelled
drawn compounds coincide, but its distinct source atom-map numbers block a
label-preserving proof. The API keeps that record's relative group.

In Rust, `reduce_meso_stereo_group(mol, index, max_visits)` returns
`Result<MesoStereoGroupReduction, MesoStereoGroupError>`. Its witness list holds
`MesoStereoWitness { mask, atom_permutation }`, using core `AtomIdx` values.

## Bounds and unsupported inputs

Hard graph bounds are **64 atoms** and **128 bonds**. The reused family boundary
allows at most six relative groups / 64 representatives, giving at most 32
selected-group pairs. A single shared probe budget covers all pairs, including
candidate selection/pruning. Default Python budget is 1,000,000 probes; accepted
explicit budgets are integers from 1 through 10,000,000. Bool/float/text indices
and budgets are refused.

Budget exhaustion returns an error with no output reduction. It is never
reported as a failed equivalence proof. Invalid group indices/metadata and ABS
selection also return errors before an unsafe graph lookup.

The proof boundary rejects unspecified locally tetrahedral sites (one or zero
H and four single-bond ligand slots), including achiral cases it has not proved
safe to ignore. It also rejects wildcard/R groups, incomplete or nonunique
stereo frames, allene/non-tetrahedral/axial assignments, query/dative/zero-order
bonds, loops/parallel edges, stashed directions, and conflicting or unaccounted
literal directions. These checks apply globally, including unselected centers,
so a proof cannot silently ignore another configuration. Relative small-ring
bonds below eight atoms fail through the existing inversion boundary.

No 3D feasibility, stereogenicity cleanup, full stereo canonicalization,
canonical representative selection, group sorting, compound counting or
bond-indexed CXSMILES extension is added. Same-number mixed atom/bond grouping
follows the #9705 proposal and existing storage contract; BIOVIA interoperability
has not been independently confirmed.

## Validation

The native proof is checked independently with RDKit **2026.09.1**, legacy
stereo perception and `removeHs=False`. The reference directly flips selected
atom tags and bond E/Z states, then regenerates neighbor directions. Each native
witness must appear in RDKit's complete chirality-aware substructure mappings;
explicit attribute/label/group checks further constrain that mapping. Direct
references have no conformer: a separate molecule supplies drawing coordinates
for V3000 tests, avoiding stale coordinates influencing inverted configurations.

Twelve original synthetic fixtures cover atom-only/mixed/bond-only redundancy,
true racemates, center-only flips, isotope/map distinctions, existing ABS,
independent relative groups, explicit H and substituted alkene references.
Twenty-four numberings, AND/OR and stored/reverse group order produce:

- **672 decisions:** 384 promotions and 288 conservative retentions;
- **768 independently verified graph witnesses** across 1,056 context pairs;
- **672 V3000 round trips**, keeping drawing, ABS members and other groups;
- zero decision, witness, allowed-configuration-set, metadata or source-mutation
  mismatches.

Before/after allowed sets are checked using drawn InChI, not as enhanced-record
identity keys. A separate diagnostic reproduces RDKit canonicalization dropping
mixed bond members; its atom-only meso reduction is a passing control. Neither
mixed canonicalization nor canonical strings are used as the reduction oracle.

Local validation passed 948 Rust library tests, 16 new proof contracts, 14
family contracts and 15 inversion contracts: **993 passed, one pre-existing
library test ignored**. Python passed **142 tests** across proof reduction,
family expansion, inversion and checked groups. Strict Clippy, formatting,
static binding inventory, product naming, the refactor budget and strict MkDocs
also passed; no refactor allowances were added.

Fixture provenance, MIT OR Apache-2.0 licensing, minimization, expected decisions,
pinned settings and payload hash are in
`validation/rdkit_issues/stereo/meso-group-proof.json`. No upstream attachment
or proprietary coordinates were copied. Regenerate the hashed evidence with:

```sh
python scripts/validate_meso_stereo_groups.py
```

This is a local, uncommitted, unpublished source candidate in
`/private/tmp/chematic-meso-groups`, preserving previous candidates through
RDFile. Its patch is `validation/rdkit_issues/meso-stereo-groups-2026-10-15.patch`
and the record is `validation/results/rdkit-meso-stereo-groups-2026-10-15.json`.
The patch requires the preceding group construction, inversion and family APIs.
