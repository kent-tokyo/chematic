# Typed binary bit-pair co-occurrence

`chematic.rdkit_compat.rdInfoTheory.BitCorrMatGenerator` now counts how often
selected binary fingerprint bits occur together. This adds a missing chematic
capability alongside the typed InfoTheory input work associated with open
[RDKit #9601](https://github.com/rdkit/rdkit/issues/9601). RDKit already provides
the generator; the issue does not report a co-occurrence algorithm defect.
This is an independently implemented Python/NumPy API without an RDKit runtime
dependency, not a complete solution to the wrapper tracking issue.

```python
from chematic.rdkit_compat import ExplicitBitVect, rdInfoTheory as I

counter = I.BitCorrMatGenerator()
counter.SetBitList([5, 1, 3])  # Caller order is significant.
for bits in [[1, 3], [1, 5], [1, 3, 5], [3]]:
    fp = ExplicitBitVect(8)
    for bit in bits:
        fp.SetBit(bit)
    counter.CollectVotes(fp)
assert counter.GetCorrMatrix().tolist() == [2.0, 1.0, 2.0]
```

For selected IDs `[b0, b1, b2, ...]`, the result orders pairs
`(b1,b0), (b2,b0), (b2,b1), ...`. Each entry is the number of collected records
with **both** bits on. No diagonal is returned. This follows the
[RDKit InfoTheory API](https://www.rdkit.org/docs/source/rdkit.ML.InfoTheory.rdInfoTheory.html).
Despite the historical `CorrMatrix` name, these are unnormalized counts, not
Pearson correlation coefficients, probabilities or evidence of predictive value.
Empty vectors and records with only one selected bit add no pair counts.

## Typed source and state contract

`SetBitList(bitList)` accepts 2..1,024 unique integer IDs in `0..4294967294`,
including bounded iterables/generators and NumPy integer arrays. Python/NumPy
bool values, floats, text, duplicate IDs and out-of-range IDs refuse. Selected
order is retained; the caller container is not stored.

`CollectVotes(bitVect)` uses the structural `GetNumBits()` / `GetOnBits()`
binary-vector interface. chematic and RDKit explicit vectors, RDKit sparse
vectors and caller-defined lazy sources work. Plain ID rows, packed strings,
bytes and count vectors are not alternate source types. Width must be in
`1..4294967295`, contain all selected IDs and remain consistent until reset.
IDs must be unique integers inside that width, including unselected IDs.
No dense allocation proportional to the declared width occurs.

The inspected RDKit sparse getter exposes some high unsigned IDs as negative
signed values. Such views refuse rather than silently reinterpreting them;
caller sources can supply verified unsigned positions instead. The API's
configuration/width/duplicate restrictions are explicit profiles, not all-input
RDKit compatibility claims. No unsafe upstream invalid-input probes were run.

Before configuration, voting and matrix access raise `ValueError`. A configured
generator returns zeros before any votes. Successful `SetBitList`, including the
same list, resets counts, source width and resource counters; this lifecycle
matches the inspected RDKit behavior. Failed configuration/vote collection
leaves previous state unchanged. Caller getter/iteration exceptions propagate;
`AttributeError` during method lookup is treated as a missing method. Sources
stay alive while lazy views are consumed and are released afterward.

`GetCorrMatrix()` returns an independent owned one-dimensional float64 array.
Counts use uint64 internally and remain exact at the record cap. Changing input
containers or returned arrays cannot change retained counts. This describes the
generator's own operations, not side effects in caller-defined methods or
concurrent mutations of the same generator instance.

## Bounds and interpretation

- At most 1,024 selected bits, giving 523,776 pair counts (about 4 MiB for the
  internal uint64 array, excluding temporary buffers and Python objects).
- At most 65,536 raw on-bit IDs per record. Reading consumes at most one extra
  ID to detect overflow; bit-list configuration similarly consumes at most
  1,025 IDs. Duplicates refuse rather than being deduplicated.
- At most one million votes and 50 million actual pair increments per
  configuration. A vote with `k` selected bits on costs `k*(k-1)/2` increments.
  Limits are checked before changing counters; reset starts a new budget.

Wrong types raise `TypeError`; invalid values, state and bounds raise
`ValueError`; allocation failures propagate as `MemoryError`. Individual caller
methods have no timeout. There is no labeling, entropy ranking, normalization,
fingerprint generation, chemical interpretation or Rust/WASM API in this addition.
The existing [feature ranker](info-bit-ranker.md) remains a separate operation.

## Validation

262 Python tests pass: 37 new co-occurrence cases plus existing entropy,
chi-square and compatibility regressions. Tests cover exact small binary sets,
selected-bit permutations, empty/unselected bits, failed input state, ownership,
actual 1,024-selected/65,536-raw-ID boundaries and the production 50-million
pair-update boundary. The record-limit lifecycle uses a lowered test budget;
one million Python voting calls were not run.

A live RDKit 2026.09.1 oracle checks 128 datasets through explicit/sparse,
lazy NumPy-ID and chematic explicit sources: 512 final matrix comparisons,
22,240 pair counts and 12,056 paired votes. Counts agree exactly with RDKit
and independent integer membership sums. All 1,536 prefix checks and 1,536
order-reversal/duplication/reset checks pass. Two source-profile diagnostics
(width changes without reset and high signed sparse views) are separate from
ordinary comparisons. No ordinary failures were discarded.

Evidence and the executable oracle are under
`validation/results/bit-cooccurrence-2026-10-18/`. Python/NumPy performs the new
calculation; a verified existing native extension is reused to import the
package and run legacy regressions. No fresh current-main native build,
published-wheel or performance validation is claimed.

## Python reference

::: chematic._bit_correlation
