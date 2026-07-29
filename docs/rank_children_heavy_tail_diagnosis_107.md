# CIP `rank_children` heavy-tail diagnosis: issue #107 (CIP-Perf-A1)

Diagnostic pass only — no optimization in this document or its companion PR,
matching the issue's own explicit framing ("this issue is for
tracking/discussion only"). The instrumentation tool this diagnosis is built
on lives at `crates/chematic-cip/examples/rank_children_heavy_tail_diagnosis.rs`
and ships no production behavior change.

## Question

CIP-Perf-A0 (PR #106) found that even among the 99.3% of stereocenters that
resolve "trivially" at Rules 1a/1b/2 (no Rule 4b/5 needed), a small tail does
40–3700x the median comparator work. Issue #107 asks five specific questions
about *why*, in a suggested investigation order, before any optimization is
attempted:

1. Same-`NodeId`-pair re-comparison count
2. Isomorphic-subtree-pair re-comparison count
3. Whether `rank_children`'s comparison matrix is reusable across sibling groups
4. Whether duplicate leaves can be collapsed into equivalence classes before comparison
5. Cost split: digraph construction vs. comparison itself

## Method

`crates/chematic-cip/src/trace.rs`'s `DecisionStep` already records one entry
per `compare_ligands` call (`left_kind`/`right_kind`/`outcome`/`rule`/
`ranking_parent`), but not the compared nodes' own `NodeId`s. The only
production edit this diagnosis required was adding `left_node`/
`right_node: NodeId` to that struct (populated at its one construction site in
`compare.rs`) — a pure data-carrying addition to an already debug-only trace
type, no behavior change.

With real `NodeId`s available, the diagnostic tool re-traces the two
worst-case fixtures issue #107 names directly (the 89,250-comparison fused
polyphenolic macrocycle and the 36,198-comparison adamantane-cage), plus a
full-corpus aggregate over the same 5,000-molecule `SMILES.csv` CIP-Perf-A0
used, computing, per stereocenter:

- Exact repeated `(NodeId, NodeId)` pairs (finding 1).
- Isomorphic-subtree pairs, via `CipDigraph::branch_signature` — an
  already-existing, already-correct, order/numbering-invariant structural
  hash (finding 2).
- **A discriminating safety check on finding 2**, added after an advisor
  review caught that raw signature-repeat counts only bound *candidate*
  savings, not *safe* savings (see Caveats below): for every signature pair
  that repeats, whether every occurrence produced the same canonical
  comparison outcome (`Higher`/`Lower`/`Equal`/`Unresolved`, normalized for
  which side is which), or whether the same signature pair produced more than
  one distinct outcome somewhere in the corpus.
- Whether a `rank_children` call's whole sibling-signature *set* recurs
  elsewhere in the same resolution, i.e. whether the entire pairwise matrix —
  not just individual pairs — is redundant (finding 3).
- What fraction of comparisons involve a leaf node (`MultipleBondDuplicate`,
  `RingDuplicate`, or `ImplicitHydrogen` — childless by construction, decidable
  by atomic number alone) (finding 4).
- Construction cost (digraph build + full `expand_all`) vs. comparison cost
  (`rank_children`'s pairwise fill on an already-fully-materialized graph)
  (finding 5).

## Results

### Named worst-case fixtures (from the issue directly)

| Fixture | total comparisons | [1] distinct NodeId pairs / max repeats | [2] distinct signature pairs / max repeats | [2] safety check | [4] leaf-involved / both-leaf |
|---|---|---|---|---|---|
| Fused polyphenolic macrocycle (atom 65) | 89,250 | 453 / 640 | 239 / 8,695 | **233/233 buckets homogeneous, 0 mixed** | 77.2% / 10.6% |
| Adamantane cage (atoms 25, 27, 30 — 3 equivalent centers) | 36,198 each | 1,093 / 118 | 84 / 10,448 | **82/82 buckets homogeneous, 0 mixed** | 98.5% / 50.1% |

### Full-corpus aggregate (4,188 stereocenters examined)

- Total comparisons: 2,096,066.
- **Upper-bound** savings if every repeated isomorphic-subtree-pair comparison
  were memoized: 2,000,261 (95.4%).
- **Discriminating safety check, corpus-wide**: 52,908 repeating
  signature-pair buckets total, **52,908 homogeneous / 0 mixed** — every
  center's within-resolution signature-pair repeats produced the same
  outcome every time, on this corpus. 0/4,188 centers had any mixed bucket.
- Leaf-involved comparisons: 89.1% of all comparisons touch at least one
  leaf node; 33.4% are leaf-vs-leaf (decidable by atomic number alone, no
  recursion into either side needed).
- Construction (digraph build + eager `expand_all`): 1,300.1ms total.
  Comparison (pairwise fill on an already-fully-materialized graph): 759.7ms
  total. Ratio 0.58x.

## Answering the issue's five questions

1. **Same-`NodeId`-pair re-comparison is real (up to 640 repeats of one pair)
   but not itself a bug.** `rank_children`'s own pairwise fill visits each
   `(i, j)` exactly once per call; a repeated `NodeId` pair means the *same*
   descendant pair was reached again from a *different* ancestor comparison
   elsewhere in the digraph — the same phenomenon finding 2 measures with a
   coarser (isomorphism, not identity) key.
2. **Isomorphic-subtree-pair re-comparison is the dominant cost pattern, and —
   critically — it was empirically safe to collapse everywhere tested.**
   95.4% of all comparisons corpus-wide share a signature pair with at least
   one other comparison, and every one of the 52,908 repeating buckets
   (across the full corpus and both named worst-case fixtures) produced the
   same outcome on every repeat. This is the strongest, most directly
   actionable finding in this diagnosis.
3. **`rank_children`'s comparison matrix is reusable across sibling groups
   most of the time.** In the polyphenolic fixture, 70/150 `rank_children`
   calls share their exact sibling-signature set with another call in the
   same trace; in the adamantane fixture, 367/374 do. A cache keyed on the
   sibling set (not just individual pairs) would avoid rebuilding whole
   matrices, not just individual cells.
4. **Duplicate leaves dominate comparison volume and are cheap to special-case
   without touching the general recursive path.** 89.1% of all comparisons
   corpus-wide involve at least one leaf node; 33.4% involve two. A leaf vs.
   leaf comparison is decidable by atomic number alone — no digraph
   recursion, no `rank_children` call, no `compare_by_level` machinery needed
   — making this the cheapest and least risky candidate for a narrow,
   separately-gated fast path.
5. **Construction (1,300.1ms) numerically exceeds comparison (759.7ms) in
   this measurement, but this split does not represent two disjoint slices of
   one real production run** (see Caveats).

## Caveats — read before treating any number above as a target

- **Finding 2's "upper bound" vs. "safe" distinction is the load-bearing
  correction in this diagnosis.** `branch_signature` hashes atomic number and
  isotope only. It does **not** encode MANCUDE fractional atomic numbers
  (`CipNode::atomic_number`, a separate field the signature never reads) or
  which ancestor a `RingDuplicate` closes back to (it hashes `closure_atom`'s
  *element*, not its identity). Two subtrees with equal signatures are
  therefore isomorphic *up to atom identity*, not proven identical for every
  input this comparator can see — a signature-only cache is a **candidate**
  key, not a **proven-safe** one. The 0-mixed-buckets result above is real,
  corpus-measured evidence that the candidate key happens to be safe on every
  case this diagnosis ran — it is not a proof that it is safe in general.
  Issue #107's own "First implementation candidate" section already
  anticipated this: it specifies the cache key must include comparison
  rule/mode, `MancudeContext` identity, and budget semantics in addition to
  the two `NodeId`s — a real implementation should use that richer key, not
  the bare `branch_signature` this diagnosis used as a measurement
  convenience.
- **Finding 5's construction number is not production cost.** The tool calls
  `expand_all(root)` to fully materialize the reachable digraph *before*
  timing comparison, specifically so lazy node materialization isn't charged
  to the comparison timer. Production's actual `expand_children` is lazy and
  only ever materializes nodes the comparator truly visits during its
  recursive descent — for centers that resolve early (most of them, per
  CIP-Perf-A0's own 99.3% figure), the real construction cost is almost
  certainly far below 1,300.1ms's implied per-center share. The 0.58x ratio
  measured here is "comparison cost given an eagerly-built graph" vs. "cost
  to eagerly build that graph" — not two disjoint slices of one real,
  lazily-executed run. A future measurement of true lazy production cost
  would need separate instrumentation of `expand_children` itself, not
  attempted here.
- Adamantane's three listed centers (atoms 25, 27, 30) are almost certainly
  symmetry-equivalent (same molecule, same comparator-size numbers to the
  last digit) — reported as three rows for completeness against the trace
  output, not as three independent data points.
- This diagnosis reuses CIP-Perf-A0's own corpus (`~/Downloads/SMILES.csv`,
  5,000 molecules, same file referenced in `docs/cip_accurate_rfc.md`'s
  MANCUDE-Decision-A0 entry) rather than a fresh draw, to keep this result
  directly comparable to the issue's own baseline numbers.

## What this does and doesn't authorize

This diagnosis identifies isomorphic-subtree-pair memoization (finding 2) and
leaf-vs-leaf fast-pathing (finding 4) as the two most promising, most
directly evidenced optimization candidates, and rank_children-matrix reuse
(finding 3) as a secondary candidate once pair-level memoization exists. It
does **not** implement any of them. Per issue #107's own absolute gate for any
future optimization PR (assignments byte-identical, skip_reasons
byte-identical, MANCUDE-Decision-A0's D/E classification unchanged, the 3
`mancude_decision_regression.rs` fixtures unchanged, RDKit oracle agreement
unchanged, `BudgetExceeded` count non-increasing), any real implementation
still needs its own correctness verification against a cache key that
includes `MancudeContext` identity and budget semantics — the bare
`branch_signature` key measured here is a diagnostic convenience, not the
key a shipped cache should use unmodified.
