# MCS non-induced and status follow-up

This record documents the source fix at `0059d41869be6287cac67f43970dea11e3e2238a`
relative to published chematic 1.2.0. It is a focused regression and API-status
check, not a broad claim of RDKit FMCS parity.

## Scope

- Host: macOS 27.0.1, Apple silicon (`arm64`)
- Python: CPython 3.13.6
- Candidate wheel:
  `chematic-1.2.0-cp313-cp313-macosx_11_0_arm64.whl`
- Candidate wheel SHA-256:
  `185009fed9a0d61e4bd09acb6ef3e2a8f9ee08382d885ec3b1e85823ac0e19ea`
- Comparator used for the 1,980-row profile: RDKit 2025.09.3

The candidate version remains 1.2.0 because this is unreleased source. The
commit and wheel hash, rather than the package version alone, identify it.

## Correctness result

Published 1.2.0 treats a closure edge present in only one input as mandatory,
which makes the result depend on input order:

| Inputs | Published 1.2.0 | Candidate |
|---|---:|---:|
| `C1CC1`, `CCC` | 2 atoms / 1 bond | 3 atoms / 2 bonds |
| `CCC`, `C1CC1` | 3 atoms / 2 bonds | 3 atoms / 2 bonds |

The candidate constructs a common non-induced subgraph: bonds absent from one
input are omitted from the common result. Nine mixed ring, chain, aromatic and
multiple-bond molecules were checked in all 36 unordered pairs and both input
orders; atom and bond counts were invariant.

With `complete_rings_only=True`, a ring-to-chain comparison still returns an
empty result. The fix therefore does not silently weaken the complete-ring
contract.

The versioned comparator profile passed 1,980/1,980 rows under its native
profile definition. The raw RDKit comparison is 1,976/1,980 because the
existing `AnyHeavyAtom` explicit-hydrogen policy differs; that known policy
difference is not changed or hidden by this fix.

## Resource-limit status

Python `find_mcs_checked` now accepts `max_search_nodes`. Exhaustive completion
returns `was_timed_out=False`; a time limit returns `True`; a native search-node
limit raises `ChematicInputError` with:

- category: `resource_limit`
- code: `mcs_resource_limit`
- format: `mcs`

WASM returns separate `wasTimedOut` and `wasResourceLimited` fields. A resource
limit is therefore no longer reported as a successful non-timeout result.

## Verification

- `chematic-smarts` release MCS tests: 42 passed
- `chematic-wasm` release MCS tests: 13 passed
- Python focused MCS tests: 17 passed
- `cargo clippy` for `chematic-smarts`, `chematic-py` and `chematic-wasm`: clean
- `cargo fmt --check` and `git diff --check`: clean

The historical frozen holdout was inspected only as a regression diagnostic;
it is no longer sealed and was not used as a release score. Its MCS exact rows
increase from 31/64 to 34/64, while five rows now end in an explicit search
resource limit. A new frozen holdout must be selected independently after the
candidate is fixed.

## Focused timing

The table shows medians from nine runs, 200 evaluations per pair. It measures
one host and three cases only.

| Pair | Published 1.2.0 | Candidate | Change |
|---|---:|---:|---:|
| ring / chain | 1.633 ms | 0.485 ms | 3.37x faster |
| benzene / toluene | 1.739 ms | 1.814 ms | 4.3% slower |
| aspirin / paracetamol | 186.317 ms | 200.759 ms | 7.7% slower |

The correctness repair is not presented as a general speed improvement. The
larger search cases require follow-up profiling, alternating blocks and a
second host before any performance claim.

## Remaining boundaries

- This does not establish full RDKit FMCS feature or output parity.
- The MCS representation is still projected through the first input rather
  than emitted as a fully generalized SMARTS query.
- The explicit-hydrogen `AnyHeavyAtom` policy difference remains documented.
- Five inspected historical rows reach the native resource limit.
- The next independent MCS holdout must not reuse the inspected v1 rows.
