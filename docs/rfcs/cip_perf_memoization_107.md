# CIP rank_children pairwise memoization: issue #107

Companion to `docs/rfcs/rank_children_heavy_tail_diagnosis_107.md` (diagnostic pass,
PR #208): implements the optimization that diagnosis identified as most
promising, with the richer cache key that diagnosis's own caveat called for.

## Rebased onto latest main (post-#208/#211 merge)

Rebased cleanly (no conflicts) onto `main` after both PR #211 (UFF rescue
fix) and PR #208 (this issue's own diagnosis pass) merged. `compare.rs`'s
overlap with PR #208's `DecisionStep` NodeId-field addition was checked
semantically, not just for a textual merge conflict: the memoization's cache
check/insert and PR #208's trace-recording both live in the same
`compare_ligands` function, and the cache path was written to still push
exactly one `DecisionStep` per logical call regardless of hit/miss — so
nothing needed to be deduplicated or reconciled; both changes compose as
originally designed. Re-verified fresh, not carried over from the pre-rebase
commit: `cargo build -p chematic-cip --lib` clean, `cargo test -p
chematic-cip --lib` 63/63 (61 pre-existing + 2 new), the 12 specific
integration tests this change is gated on (`rule4b_production_parity` 9 +
`uncharacterized_diagnosis` 1 + `mancude_decision_regression` 1 +
`mancude_corpus_classification` 1) 12/12, `cargo test -p chematic-chem --lib
-- cip` 45/45, `bash scripts/check.sh` full pass, and a fresh byte-identical
diff (1,286/1,286 stereo-bearing molecules, assignments + skip reasons,
0 differences) against the rebased, pre-memoization `main`. See the
Performance section below for the post-rebase timing re-measurement.

## Recap of the diagnosis

PR #208 measured that 95.4% of all `compare_ligands` comparisons corpus-wide
share an isomorphic-subtree signature with another comparison, and that every
one of 52,908 repeating signature-pair buckets measured (full 5,000-molecule
corpus + both of the issue's own named worst-case fixtures) produced the same
outcome on every repeat. But that measurement used `branch_signature` — a
hash of atomic number + isotope only, which does **not** encode MANCUDE
fractional atomic numbers or ring-closure ancestor identity — as an explicit,
disclosed **measurement convenience**, not a key a shipped cache should use.

## What shipped: memoization keyed on the literal NodeId pair, not a content signature

`PairwiseCacheKey` (`crates/chematic-cip/src/compare.rs`), scoped to one
`CompareContext` (one stereocenter resolution — never global, never shared
across resolutions or molecules):

- `left`/`right`: the literal `NodeId`s, canonicalized to `left.0 <= right.0`
  so `(a, b)` and `(b, a)` share one entry — the stored outcome is inverted
  on lookup when the query is in the opposite direction (issue #107's
  "outcome-direction normalization" requirement).
- `rule_mode`: a forward-compatible discriminator (currently one variant,
  `Rule1a2` — this crate implements Rules 1a/1b/2 only).
- `mancude_identity`: the attached `MancudeContext`'s pointer address, or
  `None`.
- `budget`: the digraph's own `CipBudget` (now `Hash`-derived).

Keying on the literal `NodeId` (not a content signature) makes ring-closure/
duplicate-node provenance a non-issue by construction: a `NodeId` already
uniquely identifies one exact, immutable digraph node (including a
`RingDuplicate`'s `closure_atom` or a `MultipleBondDuplicate`'s
`duplicated_atom`/`bond_order`) for the whole digraph's lifetime — two nodes
that merely *look* alike get different `NodeId`s and therefore independently
cached entries, never conflated.

`mancude_identity`/`budget` are in fact constant for one `CompareContext`'s
entire lifetime (one digraph, one resolution) — within a single cache
instance these fields never actually vary. They're included anyway so the key
type stays correct by construction if a future change ever widens the
cache's scope, rather than relying on an unenforced "this cache only ever
sees one digraph" invariant.

## Why this is safe: only successful comparisons are cached, and reused answers are provably identical

The comparison result for a given `(left, right)` pair is a pure function of
the digraph, the `MancudeContext`, and the two nodes — it does not depend on
`CompareContext::recursive_calls`'s current value, only on whether that
counter would *exceed* `max_recursive_calls` mid-recursion. A cache hit never
recomputes, so it can never itself trigger `BudgetExceeded` — meaning
memoization can only ever *reduce* spurious budget exhaustion (a comparison
that would have re-paid its own recursive cost a second time and tipped the
shared counter over the ceiling instead gets served for free), never
increase it. This is exactly the asymmetric guarantee issue #107's absolute
gate asks for ("`BudgetExceeded`が増加しないこと").

The trace (`ComparisonTrace`) still records one `DecisionStep` per logical
`compare_ligands` call regardless of hit or miss (labeled `"1a/2 (cached)"`
on a hit) — anything reading the trace sees a complete, faithful log of
every comparison a caller actually asked for.

## Memory bound

The cache can hold at most as many entries as `compare_ligands` calls made
in one resolution — itself already bounded by `CompareContext::max_recursive_calls`
(default 1,000,000). In practice, orders of magnitude smaller: the largest
observed cache size in this measurement pass was 581 entries (the
`penicillin_core`-family fixture's tied stereocenters), and the full corpus's
actual worst-case molecule (an 11-stereocenter polyphenolic tannin, below)
never exceeds a few hundred entries per stereocenter.

## Absolute gate — measured, not assumed

All measurements below compare the exact same 5,000-molecule corpus
(`SMILES.csv`) run through this crate's own real `assign_cip_accurate_experimental`
entry point, with vs. without the memoization patch (`git stash` toggling
only `crates/chematic-cip/{budget,compare,digraph}.rs` — a genuine
before/after on identical code otherwise).

- **Assignments byte-identical**: 1,286/1,286 stereo-bearing molecules'
  `(assignments, skipped)` pairs diffed exactly (elapsed-time column
  excluded) — **0 differences**.
- **Skip reasons byte-identical**: covered by the same diff (skip reasons
  are part of the same dumped tuple).
- **MANCUDE D=36/E=0 unchanged**: the 3 fixtures in
  `tests/mancude_decision_regression.rs` (`mancude_decision_a0_fixtures`)
  pass unchanged. The full D/E classification is fundamentally a statement
  about whether a fractional atomic number ever changes a *final* resolved
  label — already covered by the byte-identical assignment diff above, which
  spans the entire corpus, not just these 3 fixtures.
- **`rule4b_production_parity` (9 fixtures) and `uncharacterized_diagnosis`
  (1 fixture) unchanged**: full chematic-cip test suite passes, 61 lib +
  9 + 1 + 1 + 1 integration tests, all unchanged.
- **`chematic-chem`'s production consumer unaffected**: `chematic-chem` is
  the one real (non-dev) downstream crate depending on `chematic-cip`
  (`crates/chematic-chem/src/cip.rs`, via `assign_cip_accurate_experimental`
  directly). Its own 45 CIP-specific tests — including
  `cip_mode_legacy_fast_matches_assign_cip_byte_for_byte` and
  `cip_mode_accurate_does_not_hide_oracle_unstable_answers` — pass unchanged.
- **`BudgetExceeded` count non-increasing**: 0 molecules in this corpus hit
  `BudgetExceeded` either before or after (see reasoning above for why this
  is structurally guaranteed, not just empirically true on this corpus).
- **Rust/Python/WASM public results unchanged**: `chematic-cip` has no
  Python or WASM binding surface today (confirmed in PR #208's own
  investigation — zero references anywhere in `chematic-py`/`chematic-wasm`,
  and no re-export from the `chematic` umbrella crate); the Rust-level
  contract is the byte-identical assignment diff above.
- **`bash scripts/check.sh`**: full pass (fmt, clippy, full workspace test
  suite, deny, publish-graph, version).

## Performance — measured, three separate numbers, never conflated

**95.4% is never reported as a speedup.** It was PR #208's *comparisons-saved
upper bound*, not a timing claim. The numbers below are wall-clock,
independently measured on the same hardware, same corpus, same code except
for this one patch.

**Two measurement passes exist**, taken in two different sessions on the same
machine, and they disagree substantially on absolute magnitude (though not on
correctness — both passes' assignments are byte-identical). The first pass
(directly below, kept for continuity) was taken before this branch was
rebased onto `main` post-PR #208/#211 merge; the second (further below) is
the **authoritative, current** measurement, taken fresh after that rebase, in
one back-to-back before/after run so both sides share identical system load.
**Do not average or reconcile the two — treat the post-rebase numbers as
current, the pre-rebase numbers as historical reference only.** The
before-side absolute times differ by roughly 2–3x between the two passes
(the after/memoized side is far more stable, since it's dominated by
per-call overhead rather than redundant recursive work) — almost certainly
system-load variance between sessions on a machine that was running many
concurrent builds across other worktrees, not a code change, since the two
passes measure the identical pre-memoization commit for their "before" side.

### Pre-rebase measurement (historical reference only)

#### Corpus-wide (1,286 stereo-bearing molecules)

| | Before | After | Ratio |
|---|---|---|---|
| Total wall-clock | 2,415.81ms | 288.39ms | 8.4x |
| p50 | 0.1103ms | 0.0749ms | 1.5x |
| p95 | 4.2482ms | 0.4451ms | 9.5x |
| p99 | 31.4959ms | 1.0989ms | 28.6x |
| max | 370.0361ms | 16.7311ms | 22.1x |

#### Today's actual worst-case molecule in the corpus (pre-rebase pass)

- Before: **370.04ms**. After: **3.73ms**. **99.2x**.
- Assignments identical both times: `[(9,R),(18,R),(20,R),(30,R),(47,R),(48,R),(51,R),(60,R),(62,S),(72,R),(81,R)]`, `skipped=[]`.

### Post-rebase measurement (current, authoritative — re-run fresh on `main` post-#208/#211-merge)

#### Corpus-wide (5,000-molecule corpus total; percentiles over the 1,286 stereo-bearing molecules)

| | Before | After | Ratio |
|---|---|---|---|
| Total wall-clock (median of 3 runs) | 1,050.47ms | 286.43ms | 3.7x |
| p50 | 0.0737ms | 0.0702ms | ~1.0x (no regression) |
| p95 | 1.7411ms | 0.4376ms | 4.0x |
| p99 | 15.4362ms | 0.9327ms | 16.5x |
| max (see below — this is two different molecules, before vs. after) | 134.31ms | 16.33ms | not a single ratio, see below |

p50 barely moves — most stereocenters were already cheap (PR #208's own
99.3%-trivially-resolved finding) and pay a small, constant cache-lookup
overhead for no benefit. The win is concentrated exactly where the diagnosis
said the cost was: the tail.

**The corpus's worst-case molecule is not the same molecule before vs. after
memoization** — a fixture that benefits enormously from caching can drop out
of the "worst" slot entirely, letting a fixture that benefits *less*
(because it has a lower cache-hit rate) take its place:

- The **pre-memoization worst case** (134.31ms) is the same 11-stereocenter
  polyphenolic tannin from the pre-rebase pass above (`Oc1cc(O)c2c(c1)O[C@H]
  (c1ccc(O)c(O)c1)[C@H](O)[C@H]2c1c(O)cc(O)c2c1O[C@@]1(c3ccc(O)c(O)c3)Oc3cc(O)
  c4c(c3[C@@H]2[C@H]1O)O[C@H](c1ccc(O)c(O)c1)[C@H](O)[C@H]4c1c(O)cc(O)c2c1O
  [C@H](c1ccc(O)c(O)c1)[C@H](O)C2`). After memoization it drops to **3.63ms**
  (**~37.0x**) — identical assignments. This molecule has a high internal
  cache-hit rate, so memoization helps it the most.
- The **post-memoization worst case** is a different molecule, an
  adamantane-amide/piperazine fixture (`O=C(NCCCCN1CCN(c2cccc3ccccc23)CC1)
  C12C[C@H]3C[C@@H](C1)C[C@@H](C2)C3`) — 123.96ms before, **16.33ms** after
  (**~7.6x**), identical `[(25,LowerS),(27,LowerS),(30,LowerS)]` assignments.
  This molecule's cache hit rate (54.68%, measured directly on its
  root-level `rank_children` call — 701 hits / 581 misses / 1282 total) is
  lower than the tannin's, so it benefits less and becomes the new tail.

Neither number alone is "the" worst-case speedup — both are reported so
neither is cherry-picked.

#### The originally-named issue #107 worst-case fixture — NOT this work's achievement

Issue #107 originally named a tannin/ellagitannin-like fixture with a
reported ~89,250 comparisons. That same fixture, measured today with **zero
caching at all** (i.e. even a hypothetical build with every cache lookup
forced to miss), needs only **77 total `compare_ligands` calls** (9 of which
this memoization serves from cache; the other 68 are misses that a
from-scratch computation would need regardless of caching). Since the
"total distinct calls needed" figure (77) is a property of the comparator's
own algorithmic structure — not affected by whether a cache exists — the
~89,250 → 77 reduction happened via **unrelated intervening commits** to the
comparator itself, before this memoization work started. This memoization
PR's actual contribution on this specific fixture is small: 9 saved calls
out of 77 (11.69% hit rate), not the full historical reduction.

#### Cache construction cost and memory bound

Per-stereocenter cache sizes (final entry count = `cache_misses`) range from
12 to 581 entries across the fixtures measured directly. `PairwiseCacheKey`
is 48 bytes (`std::mem::size_of`, measured); paired with its 1-byte
`BranchComparison` value the in-memory `(key, value)` pair is 56 bytes
before hashmap bucket overhead. At the largest observed single-resolution
cache size in this corpus (581 entries), that's on the order of **35–60 KB**
peak, freed immediately when that stereocenter's `CompareContext` drops —
never global, never persisted across stereocenters or molecules. Cache
construction itself (the `HashMap` insert on a miss) is not separately
measurable from useful computation — a miss's cost is dominated by the
recursive comparison it performs, not the O(1) insert after it.

## Tests

- `test_pairwise_cache_hit_matches_fresh_computation` — a repeated `(a, b)`
  call within one context is a cache hit, and that hit's value matches an
  independent, uncached (fresh context, fresh digraph) computation of the
  same comparison.
- `test_pairwise_cache_respects_direction` — `compare_ligands(a, b)` then
  `compare_ligands(b, a)` in the same context: the second call is a cache
  hit (not a miss/recompute), and its value is the correctly-inverted
  opposite of the first — proving the direction-normalization logic, not
  just that *some* value gets reused.
- Full existing `chematic-cip` (61 lib + 12 integration) and `chematic-chem`
  CIP (45) test suites pass unchanged — the strongest evidence, since many
  of these assert exact R/S/E-Z labels on real molecules, not just "compiles
  and doesn't panic."

## Not done here (explicit scope)

- No change to `minimize_uff`/any force-field code (unrelated to this
  issue).
- No widening of the cache beyond one `CompareContext`/resolution — issue
  #107 explicitly asked for this to stay bounded, not global.
- No change to `resolver.rs`'s own, separate, higher-level embedded-
  chirality-sign cache (`assign.rs`'s `resolve_one_pseudoasymmetric`-style
  functions already memoize *resolved signs* per embedded atom across a
  Rule 4b/5 pass) — this PR's cache operates one level lower (pairwise
  branch comparisons within one `rank_children`/`compare_ligands` call
  tree) and composes with, rather than replaces, that existing mechanism.
