# RDKit-compatible conjugation assignment

The native RDKit-compatible sanitizer now computes each atom's conjugation
eligibility once per pass. This is an incremental performance improvement
related to open [RDKit #9475](https://github.com/rdkit/rdkit/issues/9475), verified
on 2026-10-19. It changes neither the chemical rules nor any public API.

The former nested loop repeatedly evaluated a neighbor's valence, hydrogen
count and electron contribution. The new loop keeps the resulting eligibility
and degree thresholds in a local temporary table. These values depend on
topology, valence, H counts, charges and radicals; they do not depend on the
conjugated flags assigned during the pass. The temporary state is discarded
after every call. Up to 64 atoms use inline storage; larger models spill to
heap storage through the existing SmallVec dependency.

For models containing only single/dative bond types, the function returns
after resetting conjugated flags from aromatic flags. Such bonds contribute
at most one unit of valence to either endpoint, so they cannot initiate the
multiple-bond pair required by the old loop. The reset still preserves any
existing aromatic flags, including models where aromatic flags and stored bond
types differ.

## Correctness checks

The pre-change loop is retained only in tests. The existing 1,131-input
compatibility fixture prepares 1,125 models and explicitly records six parser
or sanitization refusals. Across both all-false and all-true initial conjugated
flags, the new and old implementations produce identical complete internal
models: 2,250 fixture comparisons. A second call produces the same model in
each case. Further controls cover empty models, 64/65/256-atom storage
boundaries, long polyenes, radicals, charged groups, explicit H, silicon,
hypervalent sulfur, metal coordination and independent aromatic flags.

This is native Rust source validation, including a fresh candidate build.
The measurements do not use a published Python extension as an implementation
fixture. The six refused inputs are retained in the differential log; they
are not counted as successful chemical comparisons.

## Measured stage performance

On macOS arm64 with Rust 1.97.0, nine alternating-order paired rounds compare
the frozen loop and the candidate in the same release test binary. Each
round performs 300 passes over prepared models. Parsing, sanitization setup,
model cloning and test comparison formatting are outside the timed region.

| Prepared models | Count | Old median | New median | Time reduction |
|---|---:|---:|---:|---:|
| Complete accepted fixture | 1,125 | 52.71 ms | 23.11 ms | 56.2% |
| Only single/dative bond types | 445 | 5.06 ms | 0.69 ms | 86.3% |
| Includes another bond type | 680 | 47.21 ms | 22.64 ms | 52.0% |

These are timings of `set_conjugation` on this prepared fixture, not of the
complete SMILES parser/sanitizer or an application workload. They do not
establish a package-wide speedup, a published release improvement, or RDKit
performance superiority. The upstream Enamine REAL benchmark packet was not
used. LTO was disabled and codegen units set to 16 for this local comparison;
other compiler profiles, inputs and machines may differ.

Raw rounds, refusals, compiler settings, source/fixture/lockfile hashes and the
reproduction command are under
`validation/results/conjugation-cache-2026-10-19/`. The opt-in timing test is
ignored by normal test runs and is not a CI performance threshold.
