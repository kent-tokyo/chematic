# Square-planar coordination stereochemistry (`@SP1`/`@SP2`/`@SP3`)

Date: 2026-08-12. Branch: `feat/square-planar-stereo`. Trigger: user's own
competitive-scoring exercise named this the single highest-leverage next
step -- the one P0 gap `validation/platinum/FEASIBILITY.md` measured but did
not fix: cisplatin and transplatin collapsed to the same canonical SMILES
identity, because `Chirality` was tetrahedral-only and the parser hard-erred
on `@SP1`-style syntax.

Scope: SMILES parsing, canonicalization, and writing only. MOL/SDF I/O,
`chematic-3d` embedding, and CIP/R-S-analog labeling *for* square-planar
centers are explicitly out of scope (tracked as follow-ups, not started).

## Why the OpenSMILES spec text wasn't enough

The OpenSMILES spec defines `@SP1`/`@SP2`/`@SP3` by reference to a diagram
that isn't reproduced in the text spec, and doesn't give a formula. Rather
than guess, the geometric semantics were derived empirically against
chematic's local RDKit oracle (`.venv`, 2026.03.3), matching this project's
established oracle-verification convention (CIP, aromaticity, Morgan
fingerprints all followed the same pattern).

## The trans-pairing rule

Given the 4 SMILES-order neighbor positions 0,1,2,3 (same convention as
`stereo_neighbor_order`, chematic's existing tetrahedral mechanism), each
tag names which pair of positions sits trans (~180°) to each other:

- **SP1**: (0,2) trans, (1,3) trans
- **SP2**: (0,1) trans, (2,3) trans
- **SP3**: (0,3) trans, (1,2) trans

Verified two ways:

1. RDKit 3D ETKDGv3 embedding of tagged test molecules, measuring the actual
   Pt-neighbor bond angles.
2. Cross-checked against RDKit's own documented example
   (`Docs/Book/RDKit_Book.rst`): cisplatin = `Cl[Pt@SP1](Cl)(<-[NH3])<-[NH3]`,
   transplatin = the same with `@SP2`. The derived rule reproduces both:
   SP1 puts the two `Cl` cis, SP2 puts them trans.

## The permutation-remap rule (for canonicalization)

Canonicalizing a molecule reorders neighbors relative to the SMILES they were
parsed from, so an `@SPn` tag written against the original neighbor order
must be remapped to whichever tag describes the same physical arrangement
against the *canonical* neighbor order. This piece is not obvious from the
3 base cases by hand-reasoning -- one manual derivation attempt during this
work was wrong when checked against RDKit.

Resolved by exhaustive enumeration instead of hand-derivation: all 4! = 24
neighbor permutations x 3 tags, canonicalized via RDKit, grouped by
resulting identity. Re-verified on 2 independent molecule shapes (a plain
4-distinct-ligand molecule, and a dative-bond mix using `->`/`<-` donor
syntax) -- 144/144 cases, 0 mismatches. Script:
`scripts/square_planar_permutation_oracle.py` (re-runnable independently).

The resulting rule: treat each tag as a partition of `{0,1,2,3}` into its
two trans-pairs, apply the neighbor-id permutation (original SMILES order ->
canonical DFS order) to that pair-of-pairs, and match the result against the
3 templates above to find the new tag. Implemented as
`remap_square_planar()` in `crates/chematic-smiles/src/canonical.rs`.

**Scope decision -- ring-closure/chelate shapes dropped from the oracle
script.** Two additional shapes modeling ring-closure-based chelates (e.g.
carboplatin/oxaliplatin's grammar) were tried in the Python script, but
RDKit's actual stereo-perception convention for *which* neighbor a
ring-closure digit counts as, encounter-order-wise, didn't match either of
two reasonable placements tried by hand, producing 90/288 mismatches in both
attempts (0 mismatches in the 2 branch-only shapes throughout). Rather than
keep reverse-engineering RDKit's internal ring-closure encounter-order
convention in Python, ring-closure coverage was moved to a direct Rust-level
round-trip test instead
(`chelate_ring_closure_shaped_fixture_round_trips` in
`crates/chematic-smiles/tests/square_planar_stereo.rs`), which exercises
chematic's own real parser and writer rather than a re-modeled Python
approximation of RDKit's internal convention.

## Design: new `Chirality` variant, not a parallel field

`Chirality::SquarePlanar(SquarePlanarPermutation)` was added as a 4th enum
variant rather than a separate field alongside the existing 3-variant enum.
The alternative (an `Option<SquarePlanarPermutation>` field next to
`chirality: Chirality`) would let the two fields silently disagree, and
critically would NOT force every exhaustive `match` on `Chirality` in the
codebase to be revisited.

That forcing turned out to matter: `chematic-chem/src/cip.rs` and
`chematic-cip/src/assign.rs` (plus a third site found only by reading the
code directly, `chematic-cip/src/resolver.rs`) gated "is this atom a real
tetrahedral stereocenter" on `chirality == Chirality::None`, an *equality*
check rather than an exhaustive match. Equality checks are invisible to the
compiler when a new variant is added -- an `@SP1`-tagged Pt atom would have
silently fallen through into the tetrahedral CIP algorithm and produced a
bogus R/S-shaped code with no error, no panic, nothing to catch it in
review. All three sites were changed to `!chirality.is_tetrahedral()`
(`is_tetrahedral()` is the new shared helper,
`matches!(self, CounterClockwise | Clockwise)`).

One residual, deliberately left as documented rather than fixed:
`chematic-cip/src/rule4b.rs`'s `nearest_embedded`/`embedded_chain` (an
auxiliary "find a nearby chirality-bearing atom to break a Rule-4b tie"
search) still checks `!= Chirality::None`, so it can "find" a square-planar
atom that isn't a real tetrahedral center. This is safe, not a correctness
bug: every caller eventually routes the found atom through
`resolve_chirality`, which does gate on `is_tetrahedral()` and returns
`None` for a square-planar atom -- so the search just fails closed to
`SkipReason::Tied` one level later than it ideally would, a conservative
false-negative. No fixture in this codebase reaches this path (it needs a
genuine tetrahedral Rule-4b tie *and* a reachable square-planar center in
the same digraph). See the doc comment on `nearest_embedded` itself.

## Fail-closed on data-integrity problems, unlike the tetrahedral fallback

The existing tetrahedral canonicalization fallback (no recorded
`stereo_neighbor_order`, or a length mismatch) passes the original 2-state
tag through unchanged -- safe, because "unchanged" for a 2-state tag against
an unverified neighbor order is still a valid state, just possibly not
provably correct.

That safety argument does not hold for a 3-state tag: passing an `@SPn` tag
through unchanged against a reordered-but-unverified neighbor list can
silently describe a *different, plausible-but-wrong* stereoisomer, with no
error signal -- exactly the failure class this feature exists to eliminate.
So `corrected_chirality`'s fallback paths and `remap_square_planar()` itself
drop to `Chirality::None` (unspecified) on any data-integrity problem
(missing neighbor order, wrong length, duplicate atom ids) for square-planar
centers specifically, rather than reusing the tetrahedral pass-through.
Covered by `malformed_duplicate_neighbor_order_drops_to_unspecified` in the
test suite.

## Explicitly out of scope

- MOL/SDF I/O (`mol2000.rs`/`mol3000.rs`, `stereo2d_local.rs`).
- `chematic-3d` embedding / distance geometry.
- CIP/R-S-analog labeling *for* square-planar centers (only the gate fix
  above, preventing a *wrong* tetrahedral code, is in scope -- no new
  square-planar-specific CIP rule was added).
- `@TH`/`@AL`/`@TB`/`@OH` semantic support. All are now syntactically
  recognized (via `peek_chirality_class`) and produce a dedicated
  `SmilesError::UnsupportedChiralityClass` diagnostic instead of a generic
  parse error, but none are implemented.
- `chematic-py`/`chematic-wasm` bindings: verified neither exposes
  `Chirality` as a typed value (SMILES text only), so no binding changes
  were needed.
- `chematic-smarts`: independent chirality matcher, untouched.

## Verification

- `crates/chematic-smiles/tests/square_planar_stereo.rs`: 8 tests, including
  the literal killer-benchmark condition (cisplatin != transplatin once
  `@SP`-tagged), the full 24x3 oracle-parity table re-verified end to end
  through chematic's own parser + canonical writer, idempotence, the
  malformed-input fail-closed case, and the ring-closure/chelate round-trip.
- New CIP regression tests in `chematic-chem/src/cip.rs` and
  `chematic-cip/src/tests.rs` assert an `@SP1`-tagged Pt center gets no CIP
  code at all, not a wrong one.
- Non-regression: all 18 entries of `validation/platinum/pt_corpus.jsonl`
  (none of which carry an `@SP` tag) canonicalize byte-identically
  before/after this change (`git stash` before/after diff).
- `cargo test --workspace`, `cargo clippy --workspace --all-targets -- -D
  warnings`, `cargo fmt --all -- --check` all clean.
