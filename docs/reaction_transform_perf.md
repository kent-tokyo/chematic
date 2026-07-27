# `run_reactants`/`apply_retro` performance regression: root cause and fix

Date: 2026-07-27. Trigger: an external consumer, RENKIN, measured `apply_retro`
(its own thin wrapper around `run_reactants` + fragment splitting) at
**247.5µs/call on chematic 0.6.0 / master, vs 17.3µs/call on chematic 0.4.25**
— a ~14.3x per-call cost increase, with call *counts* flat-to-lower on the
newer version (ruling out "more products/matches per call" as the mechanism).
RENKIN's own investigation (their `artifacts/perf_root_cause/` directory, not
reproduced here) bisected the regression to a single dependency-version bump
in their own repo, `chematic 0.4.25 -> 0.4.30` (a **9.15x** jump in isolation,
`7.46µs -> 68.25µs` per call on a 6-target sentinel), and stopped there,
explicitly leaving "which chematic-internal change" as future work. This
document is that follow-up: a bisect *inside* chematic's own history between
those two versions, the root cause, and the fix.

## Corpus and methodology note

`v0.4.25` has **no git tag** in this repository (tags jump `v0.4.22` ->
`v0.4.28`; 0.4.23-0.4.27 were published to crates.io without a corresponding
git tag — see `f8c9961 chore(release): fix version-sync drift...`). The
release commits `3dd8db9` (`release: v0.4.25`) and `19e410b` (`chore: bump
version 0.4.29 -> 0.4.30`, matching the existing `v0.4.30` tag) were used in
their place.

RENKIN's own Stage 2 microbenchmark (`chematic_version_benchmark.json`, 30
root target molecules × a 50-template random sample) found **no** chematic-side
regression (`run_reactants` warm: 0.87x, i.e. *faster*, on 0.6.0 vs 0.4.25).
This is real and reproduced here too (see below) — and it is also exactly why
that benchmark missed the actual regression. The mechanism (below) is specific
to **highly symmetric molecules** (plain rings, cages, `CF3`/`tBu`-style
substituents), which drug-like *root* targets rarely are, but which
retrosynthesis *intermediates* and common building blocks very often are. A
30-root-molecule sample essentially never hits it; a benchmark that also
canonicalizes fragments/intermediates does.

## Bisect: the regression is not in `chematic-rxn`

`chematic-rxn`'s own dependencies are `chematic-core`, `chematic-smiles`,
`chematic-smarts` (no `chematic-perception`, no `chematic-chem`). Scoping
`git log 3dd8db9..19e410b` to those crates plus `chematic-perception`
(a transitive dependency of `chematic-smarts`) narrows ~167 total commits in
the 0.4.25→0.4.30 range down to ~28. Only one of those touches
`chematic-rxn` itself: `16db7f1 fix(rxn): preserve E/Z geometry in
run_reactants products (closes #50)` — and its diff is a single `if`
branch selecting bond endpoint order; it cannot plausibly cause a 9x
regression, and a standalone A/B build confirmed it does not (see below).

RENKIN's own `apply_retro` calls `chematic::rxn::run_reactants`, then for
every product: `canonical_smiles` → split on `.` → `parse` each fragment →
`chematic::chem::standardize` → `canonical_smiles` again. That pulls in
`chematic-chem` (via `standardize`), which depends on `chematic-fp`, which
depends on `chematic-rxn` — i.e. the *actual* regression surface for a real
consumer is wider than `chematic-rxn`'s own Cargo.toml suggests. The real
culprit turned out to be in neither `chematic-rxn` nor `chematic-chem`: it's
in `chematic-smiles`, a dependency of both.

### Standalone A/B build (this investigation's own bisect)

A throwaway `rxnbench` crate (path-dependent on a `git worktree` checked out
to each candidate commit, not part of this repository) ran
`templates_extracted_5000.smi` (RENKIN's real, USPTO-50k-derived template
corpus, read only — never copied into this repo) against
`benchmark_targets.smi` (RENKIN's 76-molecule representative root-target
set) at chematic 0.4.25 vs 0.4.30:

| version | templates × targets | run_reactants calls | ns/call |
|---|---|---|---|
| 0.4.25 (`3dd8db9`) | 5000 × 42 parsed | 210,000 | 10,944.6 |
| 0.4.30 (`19e410b`) | 5000 × 42 parsed | 210,000 | 10,480.9 |

**No regression** — reproducing RENKIN's own Stage 2 finding exactly (root
targets alone don't trigger it). Work counts (`product_molecules`,
`fragments_total`, `canonical_smiles_calls`, `standardize_calls`) were
byte-identical between versions, confirming the two versions do exactly the
same amount of *work* here.

### The actual reproducer: symmetric molecules

Directly timing `canonical_smiles()` (repeated calls on the same parsed
`Molecule`, 200 reps, release build) at 0.4.25 vs 0.4.30:

| molecule | atoms | 0.4.25 mean | 0.4.30 mean | ratio |
|---|---|---|---|---|
| adamantane (cage) | 10 | 6.06µs | 273.09µs | **45.1x** |
| coronene (7 fused rings) | 24 | 18.60µs | 902.31µs† | **48.5x** |
| plain benzene | 6 | 3.58µs | 169.36µs | **47.3x** |
| cholesterol (steroid, chiral, low symmetry) | 28 | 21.08µs | 113.39µs | 5.4x |
| morphine (chiral, no automorphism ties) | 21 | 27.06µs | 41.47µs | 1.5x |
| fused 7-ring, steroid-like (chiral) | 22 | 21.48µs | 31.74µs | 1.5x |

† high variance on this one run; a repeat with per-call timing gave a stable
~340-460µs mean for the same molecule — still a large, real regression, just
noisier at this magnitude on a shared machine.

The pattern is unambiguous: **the slowdown correlates with molecular
symmetry** (automorphism-rich structures), not with molecule size, ring
fusion complexity, or SSSR cost (a plain 6-atom benzene ring shows the same
~47x hit as a 24-atom fused system; a 28-atom *chiral* steroid does not).
This ruled out the initially-suspected SSSR/Horton's-algorithm and
`[rN]`/`[kN]` SMARTS-predicate changes in the same commit range (both real
commits in this window, but neither correlates with the observed pattern —
confirmed by branch-count instrumentation, not just plausibility).

### The causal commit

`be5dbb1 fix(smiles): individualize-refine + branch-and-minimize in
morgan_ranks` (2026-07-10, between `3dd8db9` and `19e410b`). This commit
fixed a real correctness bug (`canonical_smiles(parse(x))` could disagree
with `canonical_smiles(parse(y))` for two spellings of the same molecule, due
to a silent input-order-dependent tie-break) by adding an
individualize-refine branch-and-minimize step
(`enumerate_discrete_ranks`/`winning_individualized_ranks` in
`crates/chematic-smiles/src/canonical.rs`): when plain Morgan-rank refinement
plateaus with ties still present, every possible tie-resolution is tried and
the lexicographically smallest resulting string is kept. This is *necessary*
for correctness (proven in the commit's own message, and re-derived
independently here) — but for highly symmetric molecules, "every possible
tie-resolution" is a real, deliberately un-optimized combinatorial
enumeration: the commit's own code comment on `MAX_INDIVIDUALIZE_BRANCHES`
documents up to **168,219 branches** observed on real ChEMBL molecules with
multiple Boc/pivaloyl (`tBu`) protecting groups, and explicitly defers
"automorphism-aware branch pruning (nauty-style)" as future work, not
attempted there.

Instrumenting `enumerate_discrete_ranks` confirmed the mechanism directly:
adamantane produces 24 individualize-refine branches, coronene and benzene
12 each (vs. 1 branch — no ties at all — for the chiral, low-symmetry
molecules above). Critically: **every branch for a given molecule writes the
same final canonical string** (`distinct_strings == 1` in all cases measured,
even though the raw discrete rank vectors are all pairwise distinct) — the
branching is 100% correct but, for pure automorphism orbits, 100% redundant
work, exactly as the original commit's own docstring predicts
("if a cell IS an orbit, every choice yields an automorphic result... correct
but redundant").

## Fix implemented here

`winning_individualized_ranks` (the internal helper both `canonical_smiles`
and `canonical_atom_order` share) already had to write *every* branch's full
string via `CanonicalWriter::write_all()` to find the minimum — but its
caller, `canonical_smiles`, then called `write_all()` a **second, fully
redundant** time on the already-known winning ranks, on every single call
(tied or not). The fix (`crates/chematic-smiles/src/canonical.rs`) changes
`winning_individualized_ranks` to return `(ranks, winning_string)` instead of
just `ranks`, so `canonical_smiles` reuses the already-computed string
instead of re-deriving it.

This is a pure, behavior-preserving refactor: the exact same set of
candidate strings is computed, in the exact same order, and the exact same
minimum is selected (`Iterator::min_by`, like the `min_by_key` it replaces,
returns the *first* minimal element on ties — verified from the stdlib
docs, not assumed). **Every existing test passes unchanged, including the
five golden-string tests `be5dbb1` itself had to update**
(`crates/chematic-chem/src/mmp.rs`, `crates/chematic-wasm/src/tests.rs`) and
all `issue50_*` E/Z regression tests in `crates/chematic-rxn/src/transform.rs`
— those are the two families of test this specific change could plausibly
have broken, and both were checked explicitly, not just "the suite is green."

### What this fix does *not* do

It does **not** implement automorphism-aware branch pruning. That is a
substantial, separate undertaking (the original commit's own author
explicitly deferred it as "nauty-style... not attempted here", and it is not
attempted here either) that would meaningfully reduce or eliminate the
45-48x regression on genuinely symmetric molecules (plain rings, cages,
protecting groups). What's shipped here removes exactly one fully-redundant
write per call — a real, measurable, and completely safe win, but a partial
one. See "Remaining gap" below.

## Internal work counters (Phase 2)

`chematic-rxn` gained a `perf-instrumentation` Cargo feature (off by default,
zero cost when disabled — see `crates/chematic-rxn/src/perf_counters.rs`)
exposing process-global counters: `run_reactants_calls`,
`reaction_parse_calls`, `reactant_query_match_calls`, `vf2_match_count`,
`match_combination_count`, `build_product_calls`,
`product_sets_before_dedup`, `product_molecules_built`,
`atoms_copied_to_products`, `bonds_copied_to_products`. These cover
`run_reactants`'s own internals; the RENKIN-side concepts of
`split_fragments`/dedup/fragment-filtering are apply_retro-level (not part of
`run_reactants`'s public API), so the benchmark harness below tracks its own
equivalents (`canonical_smiles_calls`, `standardize_calls`,
`fragments_before_filter`, `fragments_parsed_ok`,
`errors_swallowed_by_unwrap_or_default`) rather than adding RENKIN-specific
concepts to the library itself. `fragments_parsed_ok` counts fragments that
parsed successfully, not post-filter survivors — this harness does not
reimplement RENKIN's own aromatic-atom-without-a-ring-closure fragment
filter, so the name doesn't promise one.

## Corpus-weighted benchmark (Phase 1)

`crates/chematic-rxn/examples/reaction_transform_perf_report.rs`. Accepts
external corpora via `RENKIN_TEMPLATES` / `RENKIN_PROBE` env vars (never
requires copying RENKIN's data into this repo); falls back to small,
hand-authored fixtures committed at `crates/chematic-rxn/fixtures/` when
unset. Mirrors RENKIN's real `apply_retro` → `split_fragments` → `standardize`
call pattern (not `run_reactants` in isolation — that's exactly what missed
the regression the first time), and **feeds round-1 fragments back in as
round-2 probe molecules** (`RXN_PERF_DEPTH`, default 2) so symmetric
intermediates — not just root targets — enter the measured population,
matching how a real multi-step retrosynthesis search actually behaves.

```
RAYON_NUM_THREADS=2 cargo run --release -p chematic-rxn \
    --features perf-instrumentation --example reaction_transform_perf_report
```

(`RAYON_NUM_THREADS=2` is set for environment comparability with RENKIN's own
methodology; this harness itself is single-threaded — see the `ponytail:`
note in the file.)

### Before/after, hand-authored witness fixtures (7 templates × 11 molecules × depth 2, 336 calls)

| | before (main) | after (this fix) | improvement |
|---|---|---|---|
| total elapsed | 702-758ms | 596-598ms | ~15-22% |
| mean/call | 2.09-2.26ms | 1.78ms | ~15-21% |
| p50 | 109.5-114.7µs | 101.7-103.0µs | ~7-11% |
| p95 | 3.01ms | 2.77-2.90ms | ~4-8% |
| p99 | 67.8-68.1ms | 57.3-57.6ms | ~15% |
| max | 72.9-86.2ms | 65.3-66.2ms | ~10-24% |

**The two arms were not instrumentation-symmetric**, disclosed rather than
re-measured away: the "before" run used a `main` checkout with
`perf_counters.rs` copied in but `transform.rs` left uninstrumented (its
counters are permanently zero, i.e. dead code on that arm), while "after" was
built with `--features perf-instrumentation` live, paying real
`AtomicU64::fetch_add` overhead on every `run_reactants`/`build_product` call
(840 `build_product` calls in this run alone). That bias is conservative —
the fixed arm carried extra cost the baseline didn't — so **the improvement
above is understated, not inflated**.

(Ranges are two repeated runs per side, release build, `RAYON_NUM_THREADS=2`,
same machine.) This is a real, honest, consistent improvement across the
mean *and* the tail (p95/p99/max) — but a modest one, not a restoration of
RENKIN's full ~9x. **This is expected and disclosed, not a shortfall being
hidden**: the witness corpus mixes symmetric (large-regression) and
asymmetric (small/no-regression) molecules, and this fix only removes the
universal single-redundant-write cost, not the underlying branch-enumeration
cost that dominates the symmetric cases. RENKIN's own full 30-target
integration re-measurement (out of scope here per the task's explicit
instruction) is the right place to see the aggregate, corpus-realistic
effect at RENKIN's actual scale and template/molecule mix.

## Remaining gap (explicitly not fixed here)

The dominant cost for genuinely symmetric molecules — up to ~48x on plain
rings/cages in isolation — is **not** eliminated by this fix. Closing it
requires detecting when a tied refinement cell is a genuine automorphism
orbit (in which case every branch is provably redundant, per `be5dbb1`'s own
docstring) and exploring only one representative, which is real,
correctness-sensitive graph-automorphism work, not a small patch. Two
approaches were evaluated and rejected for this PR:

- **Byte-identical rank-vector dedup before writing**: measured directly
  (adamantane/coronene/benzene) — 0% duplicate raw rank vectors even though
  100% of branches for a given molecule produce the *same final string*.
  This confirms the redundancy is real but means naive vector-level dedup
  buys nothing; the equivalence only shows up after writing.
- **True-twin detection** (atoms with identical neighbor sets, e.g. a `CF3`
  group's three fluorines): a rigorous, cheap, and correct pruning rule for
  that *specific* symmetry pattern — but it does not fire on any of the
  ring/cage witnesses measured above (their symmetry is rotational, not
  twin-based), so it would not move this PR's own positive witness and was
  not implemented here. Left as a well-scoped, believed-correct follow-up.

Filed as future work, not started: general automorphism-orbit-aware pruning
in `enumerate_discrete_ranks`/`winning_individualized_ranks`.

## Why fix forward, not pin to 0.4.25

RENKIN's own `next_gate.json` lists "pin chematic back to 0.4.25" as one
candidate. That was not the direction taken: 0.4.25 predates the E/Z
geometry fix (`16db7f1`, issue #50) and the individualize-refine correctness
fix (`be5dbb1`) themselves — pinning back would silently reintroduce two
confirmed correctness bugs (nondeterministic/wrong E/Z geometry in
`run_reactants` products; input-order-dependent canonical SMILES) in
exchange for the performance this fix restores only partially anyway. This
is a chematic-internal regression, so the chematic-internal fix (forward,
not backward) is the correct direction, per explicit user instruction for
this task.
