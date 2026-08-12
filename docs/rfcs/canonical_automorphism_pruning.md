# Canonical SMILES: automorphism-orbit-aware branch pruning

Date: 2026-07-29. Branch: `fix/canonical-automorphism-pruning`. Trigger: RENKIN
(external sibling project) reported `apply_retro`/`run_reactants` at up to a
45x slowdown between chematic 0.4.25 and later versions, traced to
`canonical_smiles()`'s individualize-refine branch-and-minimize search
introduced by `be5dbb1` (2026-07-10). PR #189 (merged `97c87e3`) fixed a
narrower, real bug (a redundant second `write_all()` call, ~13.4% measured
improvement) but explicitly left the combinatorial branch explosion itself
unfixed (`docs/rfcs/reaction_transform_perf.md`, "Remaining gap"). This is that
follow-up.

## Root cause, precisely

`be5dbb1` fixed a genuine correctness bug (`canonical_smiles(parse(x))` could
disagree with `canonical_smiles(parse(y))` for two spellings of the same
molecule) by adding individualize-refine branch-and-minimize: whenever plain
Morgan-rank refinement (`morgan_ranks`) plateaus with ties, every possible
tie-resolution is tried (`enumerate_discrete_ranks`) and the
lexicographically smallest resulting string kept. This is *necessary* for
correctness, but the enumeration is combinatorially blind to *why* a cell is
tied: a 1-WL / Morgan-rank refinement cell can be tied because it is a
genuine automorphism orbit (every choice of representative yields an
identical result -- provably redundant to explore more than one) or because
it merely *contains* an orbit alongside genuinely different candidates (every
choice must be explored). The old code always assumed the latter,
unconditionally, for every tied cell -- hence the documented up-to-168,219
branches on real ChEMBL molecules with independent 3-fold-symmetric Boc/tBu
groups.

**Why Morgan-rank-only pruning would be unsound** (the naive "fix" this PR
does NOT take): merging two atoms into one orbit merely because they share
the same Morgan rank conflates "same 1-WL cell" with "same automorphism
orbit". These are NOT equivalent -- a 1-WL refinement cell can contain
multiple distinct orbits (a standard, well-known limitation of
color-refinement algorithms; see McKay & Piperno 2014 and Piperno
arXiv:0804.4881, used here for algorithmic design inspiration only, not
ported). `canonical_automorphism.rs`'s own exhaustive test suite pins a
concrete witness of this: the disjoint union of a triangle (C3) and a square
(C4) is 2-regular throughout, so 1-WL/Morgan-rank refinement can never split
it into more than one cell (both components look identical to purely local
neighbor-hash iteration) -- yet no automorphism can map a triangle vertex to
a square vertex (components of different order can never correspond under a
graph automorphism). A rank-hash-only pruning scheme would wrongly treat all
7 vertices as one interchangeable orbit.

## Algorithm overview

Four new/changed files, matching the required layer split:

- `crates/chematic-smiles/src/canonical_partition.rs` -- the writer-visible
  vertex/edge coloring (`VertexColor`/`EdgeColor`/`CanonicalColoredGraph`)
  and an *exact* partition refinement (`exact_refine`) that never merges two
  cells on hash-collision alone (full signature equality, only using a
  hash/sort for lookup speed -- standard `HashMap`/`sort`/`dedup` behavior,
  not the forbidden pattern).
- `crates/chematic-smiles/src/canonical_automorphism.rs` -- an exact,
  backtracking bijection search (`has_colored_automorphism_mapping`):
  verifies a full graph automorphism (never a partial or subgraph match),
  forces `from -> to`, keeps already-individualized singleton cells fixed
  (via cell-restricted candidate search), and independently re-verifies the
  complete mapping from scratch before returning `true`.
- `crates/chematic-smiles/src/canonical_search.rs` -- `exact_orbit_representatives`
  (pairwise union-find over one target cell's members, gated on an actual
  verified automorphism test per pair, never Morgan-rank equality) and
  `search_canonical` (streaming DFS: exact refinement, singleton check,
  orbit partition of the target cell, individualize-refine exactly one
  representative per orbit, recurse). Also the fallible bounded API
  (`CanonicalizationLimits`/`CanonicalizationError`/
  `canonical_smiles_with_limits`) and feature-gated instrumentation
  counters.
- `crates/chematic-smiles/src/canonical.rs` -- unchanged in its public
  meaning (`morgan_ranks`, `equivalent_atom_classes`, `are_atoms_equivalent`
  all still the old hash-based refinement, per this PR's own constraint).
  `winning_individualized_ranks` (the shared internal helper used by both
  `canonical_smiles` and `canonical_atom_order`) now calls the new
  orbit-pruned search by default, with the original exhaustive enumeration
  (`legacy_winning_individualized_ranks`, unchanged) kept as (a) a
  last-resort fallback for the new engine's own internal-invariant error and
  (b) the implementation behind the test-only exhaustive oracle
  (`canonical_smiles_exhaustive_oracle`).

## Why pruning cannot change surviving output

The rank vectors handed to `CanonicalWriter` (and therefore every byte of the
output string) are produced by the *exact same* `individualize` +
`refine_ranks` primitives the legacy enumeration used, unchanged. Orbit
pruning only ever decides *which atoms in a target cell get individualized
at all* -- it selects a subset of the same candidate branches the legacy
enumeration would visit; it never perturbs what a surviving branch computes.
So any leaf this engine reaches is byte-for-byte identical to what the
unpruned exhaustive enumeration would have produced for that same
individualization sequence. This is why the engine can be directly
cross-checked against `canonical_smiles_exhaustive_oracle` per molecule
(section "Small-graph and chemical fixture results" below) rather than
merely "look plausible".

**Why skipping proven-redundant branches cannot change the *minimum*
either**: if atoms `x` and `y` are in the same automorphism orbit (a real,
verified graph automorphism `phi` fixes everything already individualized
and maps `x -> y`), then the subtree reachable by individualizing `x` and
the subtree reachable by individualizing `y` are related by `phi`: every leaf
string reachable from one subtree is reachable from the other (composing
with `phi`). So the *set* of leaf strings under `x` equals the *set* of leaf
strings under `y`, and in particular their minima are equal to each other
and to the true global minimum from that cell. Exploring only one
representative therefore never removes the global minimum from
consideration -- it only removes duplicate copies of it.

## Why Morgan-rank equality is not proof, and what replaces it

`initial_invariant` (in `canonical.rs`, unchanged) packs atomic number,
degree, charge, isotope, aromaticity, and explicit-H flag into a `u64`, then
`refine_ranks` iterates FNV-hash-based neighbor aggregation to a fixpoint.
Two atoms sharing a Morgan rank means no *local, hash-summarized* feature
distinguishes them within the number of refinement rounds run -- it does
**not** mean a genuine graph automorphism relates them (the triangle/square
witness above proves this can fail even for the *whole graph*, not just one
cell). This PR never treats Morgan-rank equality as license to prune; every
orbit union in `exact_orbit_representatives` is gated on an independently
re-verified `has_colored_automorphism_mapping` call.

## The writer-visible graph coloring

Audited directly from `CanonicalWriter`'s `emit_atom`/`write_chain`/
`find_ring_closures`/`corrected_chirality` and the plain writer's
`crate::writer` (same crate, same rules per `raw_bond_direction`'s own doc
comment: the two writers are kept from silently disagreeing).

**Vertex color** (`VertexColor` in `canonical_partition.rs`) includes:
- `wildcard` (the SMILES `*` atom always emits `[*]`, element is irrelevant)
- `atomic_number` (0 when wildcard)
- `isotope` (bracket-forcing; written verbatim when present)
- `charge` (bracket-forcing; written verbatim)
- `aromatic` (lowercase vs uppercase element symbol)
- `h_state` (`atom.hydrogen_count.is_some()`/value -- `needs_bracket` is
  keyed on this being `Some` at all, independent of the numeric H count
  `implicit_hcount()` would separately compute, so `Some(0)` and `None` are
  writer-distinguishable even when the *effective* H count matches)
- `atom_map` (always emitted verbatim when present -- see "atom-map policy"
  below)
- `chirality` (`@`/`@@`/none discriminant)
- `stereo_unique` (see "judgment call" below -- not itself a literal writer
  attribute, but the mechanism this PR uses to stay sound around chirality
  and E/Z without deriving full parity-aware automorphism math)

**Edge color** (`EdgeColor`) includes:
- `order_class`: all 13 `BondOrder` variants mapped to distinct classes
  (Single/Double/Triple/Quadruple/Aromatic/Up/Down/Zero/Dative/QueryAny/
  QuerySingleOrDouble/QuerySingleOrAromatic/QueryDoubleOrAromatic) -- never
  collapsed the way `canonical.rs`'s *separate* `bond_order_value` (used
  only for Morgan-rank neighbor hashing, unchanged) collapses Single/Up/
  Down/Dative into one value; this module needs the finer distinction to
  avoid false automorphisms between e.g. a plain single bond and a dative
  bond.
- `from_is_donor`: for a `Dative` bond only, `true` when the atom queried
  *from* is `bond.atom1` (`chematic_core::BondOrder::Dative`'s own doc:
  `atom1 -> atom2` is donor -> acceptor). Always computed relative to the
  query side (`edge_color(from, bidx)`), so comparing `edge_color(u, ..)`
  against `edge_color(phi(u), ..)` for the *same* relative direction
  correctly preserves donor/acceptor identity under a candidate mapping,
  without needing a separate direction-normalization step.

**Excluded attributes, and why each is provably safe to exclude**:
- `cip_code` (assigned CIP R/S/E/Z label): never read by `CanonicalWriter`
  or `crate::writer` (grep-verified -- `chematic-chem`'s `assign_cip` writes
  it, but this crate never emits it in `canonical_smiles`/`write` output). Not
  writer-visible, safe to exclude.
- `Molecule::stereo_groups` (enhanced-stereo `|&1:0,3|`-style annotations):
  never referenced anywhere in `crate::canonical` or `crate::writer`
  (grep-verified). Enhanced-stereo group output is a CXSMILES-specific
  concern (`crate::cx`), out of scope for plain `canonical_smiles`. Not
  writer-visible for this module's purposes, safe to exclude.
- Raw stored `/`/`\` position/orientation as a literal attribute of one
  specific bond: per this PR's task spec, "the raw stored position of `/`/
  `\` is a spelling detail, not a chemical difference" -- `CanonicalWriter`
  itself re-derives/re-orients direction per write-order via
  `effective_order`/`normalize_ez`/`resolve_ez_markers`, so encoding the
  *raw* literal Up-vs-Down value as an immutable per-bond color (rather
  than its real cis/trans meaning) would be over-fitting to spelling. This
  PR does not attempt to derive the fully resolved, write-order-independent
  E/Z meaning ahead of the search (see the judgment call below for why);
  instead of building a shortcut that could be wrong, it takes the strictly
  conservative route of never pruning around such bonds at all, which is
  safe (never a false merge) at the cost of some missed optimization
  opportunity in molecules that mix E/Z stereo with unrelated symmetric
  groups (verified NOT to disable pruning of unrelated groups --
  `ez_pinning_does_not_disable_unrelated_cf3_pruning` in
  `canonical_search.rs`).

## Atom-map policy

`CanonicalWriter::emit_atom` always writes `:{map}` when `atom.atom_map` is
`Some`, unconditionally -- canonical SMILES **preserves** atom maps verbatim,
it never strips or renumbers them. `VertexColor::atom_map` therefore
participates in vertex color like every other writer-visible attribute: two
atoms with different (or present-vs-absent) atom-map numbers are never
considered automorphic.

## Judgment call: stereo-bearing atoms and their direct neighbors are pinned

**The call**: any atom with `chirality != Chirality::None`, or incident to a
bond carrying real-or-potential direction information (`BondOrder::Up`/
`Down`, or a non-`None` `Molecule::bond_direction` stash), together with
*all of that atom's direct graph neighbors*, gets a globally unique vertex
color (`stereo_unique: Some(atom_index)`) -- i.e. such an atom can only ever
be mapped to itself by any accepted automorphism.

**Why this is necessary, not just convenient** (found empirically, not
predicted up front): an earlier version of this PR pinned *only* the
stereo-bearing atom itself, reasoning "the tag is on that one atom, so
forcing it to map to itself is enough". This is wrong. A tetrahedral
center's `@`/`@@` tag is defined relative to the *order* of its direct
neighbors (`Molecule::stereo_neighbor_order`); transposing two of those
neighbors with each other inverts the encoded configuration even though the
stereocenter atom itself never moves anywhere. This was caught by an
*existing* regression test,
`ring_digit_reuse_inside_stereocenter_branch_minimal` (a stereocenter with
two topologically-swappable ring-neighbor `CH2` atoms): the atom-only pin
passed every per-atom color/edge check yet a candidate automorphism swapping
the two ring neighbors silently inverted the written chirality tag between
the first and second canonicalization of the same molecule (an idempotence
break, i.e. a real false-prune bug). Pinning every direct neighbor too
closes this: chirality parity depends only on the stereocenter's direct
neighbor list, never on anything further away, so pinning that whole list
pointwise is sufficient (not merely convenient) to prevent any accepted
automorphism from touching it. The same 1-hop rule is applied uniformly to
E/Z-direction-bearing bonds (a stereogenic alkene end's *unmarked*
substituent -- implied via substituent count rather than an explicit bond
marker -- is one hop from the marked bond's endpoint and must stay pinned
too, for the analogous reason).

**Cost**: false negatives only (section 8's sanctioned trade) -- a
stereocenter's own neighbors, and the stereo-bearing bond's endpoints, are
never pruned even in cases where a more careful (full parity-aware
automorphism check) analysis might prove it safe. Verified this does not
regress the required performance targets: none of benzene/adamantane/
coronene/CF3/tBu/Boc/pivaloyl carry any chirality or E/Z bonds, so the pin is
a complete no-op for the high-symmetry performance corpus. Verified the pin
does not *disable* pruning of an unrelated symmetric group sharing a
molecule with an E/Z system (`ez_pinning_does_not_disable_unrelated_cf3_pruning`).

**Alternative considered and rejected**: deriving a full parity-aware
automorphism check (permutation-parity of the mapped neighbor order via the
existing `permutation_is_odd` machinery, mirroring `corrected_chirality`'s
own trick) would recover the strictly-necessary rather than
strictly-sufficient pruning boundary. Not implemented in this PR: the
1-hop-pin is simple, easy to verify by exhaustive small-graph and chemical
fixture testing, has zero cost on every performance-gated fixture, and the
false-negative cost (per section 8) is explicitly sanctioned. Left as a
well-scoped future refinement if a real molecule is found where it matters.

## Orbit-pruning safety argument (summary)

1. `initial_partition` combines the current node's Morgan-derived `ranks`
   (already encoding every individualization done so far) with the full,
   static `VertexColor` -- so an already-individualized singleton, or an
   atom distinguished only by an attribute `initial_invariant` doesn't carry
   (isotope/charge/atom-map/stereo), is never merged into a larger cell by
   this step.
2. `exact_refine` refines that composite partition to a fixpoint using full
   signature equality (own color + sorted multiset of (edge color, neighbor
   cell)) -- never a hash-collision shortcut.
3. `exact_orbit_representatives` only merges two members of the *current*
   target cell when `has_colored_automorphism_mapping` returns `true` for
   that specific pair, under the *current node's* coloring (not a
   global/root-only symmetry group) -- so previously-individualized
   singletons are automatically respected (the cell-restricted candidate
   search in `extend_mapping` can only ever map a singleton to itself).
4. `has_colored_automorphism_mapping` builds a complete bijection (never a
   partial/subgraph match), independently re-verifies it in full
   (`verify_full_bijection`) before returning `true`, and forces `from ->
   to`.
5. A false negative (missing a real orbit) only costs performance --
   `exact_orbit_representatives` still explores that branch. A false
   positive is prevented **given step 1's `ranks` correctly reflects every
   individualization already committed** -- every union is gated on an
   independently re-verified bijection over `VertexColor`/`EdgeColor`, never
   rank/hash equality alone.

   That qualification is not vacuous (independent Round-2 false-prune
   audit, PR #193): `ranks` is produced by the pre-existing, unchanged
   `individualize` + `refine_ranks`. `individualize` itself is exact integer
   arithmetic (zero collision risk), but `refine_ranks`'s subsequent
   `fnv_hash_sequence`/`normalize_ranks` groups by raw 64-bit hash-value
   *equality* -- so a genuine hash collision there could in principle
   re-merge an already-individualized atom with a formerly-tied sibling,
   and step 1's `VertexColor` alone cannot catch that (it deliberately
   carries no search-time individualization history, only intrinsic atom
   attributes). This exposure is not new: `refine_ranks` rank-equality is
   *already* the sole basis for this crate's pre-existing
   `equivalent_atom_classes`/`are_atoms_equivalent` public APIs, with
   identical collision risk, unrelated to orbit pruning. What this PR
   changes is the *consequence*: the legacy exhaustive engine would turn a
   hypothetical collision into redundant-but-still-correct
   over-exploration (every member of a wrongly-merged cell still gets
   individualized); this PR's pruned engine could instead silently skip a
   genuinely distinct branch. Not observed on any fixture, the `n<=5`
   exhaustive suite, hundreds of randomized fuzz trials, or the
   5,000-molecule corpus; would require a correlated 64-bit FNV-1a collision
   reconstructing an entire real symmetry's cell structure to manifest.
   Closing it fully would mean threading a parallel, hash-free
   individualization-state vector through the search's hot path --
   judged out of this PR's scope (a bigger change, defending against a risk
   already implicitly accepted crate-wide), so disclosed here instead of
   left as an unqualified "structurally impossible" claim. See
   `canonical_search::exact_orbit_representatives`'s doc comment for the
   same account in the code itself.

## Old vs. new budget semantics

Old: `enumerate_discrete_ranks` capped at `MAX_INDIVIDUALIZE_BRANCHES =
10_000`; on exhaustion, the child loop silently `break`s and the caller
takes `.min_by`/`.min_by_key` over whatever prefix of branches had been
explored, presented as if it were the true answer. A sufficiently small
budget (unreachable at 10_000 for any real molecule so far, but directly
reachable once a caller can configure the budget) can make the *very first*
node return an empty `Vec`, which `.unwrap_or_default()` turns into an empty
canonical SMILES string -- silently wrong, not merely slow.

New: `canonical_smiles`/`canonical_atom_order` (the existing infallible
signatures, kept source-compatible) use
`CanonicalizationLimits::unbounded()` -- no cap on search nodes or
automorphism tests. This is safe specifically *because* orbit pruning
removes the combinatorial blowup that made the old fixed cap necessary:
recursion depth is bounded by atom count regardless of branching factor
(mathematical termination was never actually at risk, only wall-clock cost
was), and with orbit-redundant branches pruned away, the *necessary*
branch count for every fixture and the 5,000-molecule corpus measured here
stays small. `CanonicalizationError::InvalidInternalMapping` (should never
happen) falls back to the original, unchanged exhaustive enumeration
(`legacy_winning_individualized_ranks`) so a hypothetical internal-invariant
bug in the new engine degrades to "slow but correct" rather than a panic or
silently wrong output. New callers that want a hard ceiling instead should
use the new fallible `canonical_smiles_with_limits`, which returns
`Err(CanonicalizationError::SearchBudgetExceeded { .. })` -- never `Ok` with
a truncated/wrong string, never an empty string, never an input-order-
dependent silent fallback.

## A verified side effect: coronene idempotence

`tests/canonical_robustness.rs`'s `coronene_canonical_known_bug` was
`#[ignore]`d as a known limitation ("coronene canonical SMILES oscillates").
It now passes under this engine. Before un-ignoring it, this was verified
(not assumed): the engine's output matches
`canonical_smiles_exhaustive_oracle` for coronene in *both* directions
(parsing the original aromatic SMILES, and re-parsing the resulting
canonical string), the two canonicalizations agree with each other
(idempotence), and the result is stable across 32 rotated/reflected
relabelings of the same molecule (all in
`canonical_search::tests::unbounded_matches_exhaustive_oracle_on_symmetric_molecules`
and the `coronene_canonical_known_bug` update itself). This is reported
individually and explicitly, per this PR's own requirement for any
previously-broken case whose output changes.
