# Supplied tautomer selection with a typed score callback

`chematic.rdkit_pick_tautomer_scored(candidates, score_func, *, max_candidates=1000)`
adds custom integer scoring to supplied-candidate selection. The existing
native `rdkit_pick_tautomer` keeps its default score and behavior.
This is an adjacent capability motivated by the typed iterable/callback and
ownership work in [RDKit #9601](https://github.com/rdkit/rdkit/issues/9601).
RDKit already supports custom scores in `TautomerEnumerator.PickCanonical`;
this does not implement nanobind or resolve the whole tracking issue.

```python
import chematic

candidates = (chematic.from_smiles(s) for s in ["CC(C)=O", "C=C(C)O"])
selected = chematic.rdkit_pick_tautomer_scored(
    candidates, lambda mol: -mol.rdkit_tautomer_score()
)
assert selected.rdkit_smiles == "C=C(C)O"
```

The callback is `Callable[[Mol], SupportsIndex]`. Python/NumPy integer return
values in signed 32-bit range are supported; bool and floating-point values
are rejected. Higher scores win. Equal scores select the lexicographically
smallest native RDKit-port canonical SMILES; identical keys keep the first
source occurrence. No candidates are generated or family membership checked.
A singleton stays that source graph. No stereo cleanup is performed, and
selection does not establish physical stability or an experimental energy.

## Ownership and callback order

Each yielded native `Mol` is copied using the existing versioned native state
codec before the input iterator advances. Its graph, indexed stereo metadata,
atom maps/tags, title and typed molecule properties follow that codec's
contract. No external pickle stream or caller-supplied binary payload is read.
Derived caches and source spelling are not part of a source-text preservation
promise.

All candidates and canonical keys are prepared before the first scoring call.
The callback then runs once per candidate in input order, receiving another
independent snapshot. It can retain or edit that view without changing the
retained source snapshot or the result. The returned winner is an independent
copy of the selected source. Changes to source objects while advancing a
reused generator cannot alter already prepared snapshots.

Invalid input, iterator errors or limits produce no scoring calls and no
partial winner. Callback exceptions propagate unchanged and produce no
winner. A later callback error does not roll back earlier callback effects in
caller code. Captured external objects, files or other side effects remain
under the caller's control.

## Limits

| Input | Bound |
| --- | --- |
| `max_candidates` | 1–4096 integer, excluding bool; default 1000 |
| Per-candidate graph | 1024 atoms and 4096 bonds |
| Combined copied graph atoms | 16,384 |
| Per-candidate encoded native state | 1 MiB |
| Combined retained encoded state | 32 MiB |
| Callback score | -2,147,483,648 through 2,147,483,647 |

Duplicate candidates count toward all limits. At most one additional input
item is requested to detect count overflow. The native codec applies its own
bounds during encoding; wrapper byte limits are checked after encoding, and
combined limits can transiently include one extra candidate. Payload limits
are not a total-process-memory bound. Callbacks can retain views or perform
arbitrary work, so no execution deadline or memory bound for caller code is
promised. Empty collections raise `ValueError`; a zero-atom native `Mol` is a
valid candidate if its canonical representation is supported.

This is a distinct chematic function, not a complete `rdMolStandardize`
namespace or `PickCanonical` replacement. Newer RDKit prefers specified
stereo centers before lexical ties. For a constant score on `CC(F)O` and
`C[C@H](F)O`, this API chooses the unspecified graph by lexical key while
RDKit 2026.09.1 chooses the specified graph. That difference is documented
and excluded from parity counts. Default native selection retains its
[existing policy](tautomer-candidate-selection.md).

## Verification

The candidate uses Python for staging/selection and existing native methods
for state encoding, restoration and canonical keys. Local checks use a
hash-verified older native extension. The codec/property implementation is
unchanged from that fixture; current-main native aromaticity has a separate
nine-membered-ring change, so the reference packet deliberately covers
acyclic families. No fresh current-main binary or whole-extension equivalence
is claimed.

Tests include reuse/lifetime handling, typed properties, callback isolation,
preflight failures, score/size limits and stereo/isotope/map preservation.
The live RDKit 2026.09.1 packet compares pure callbacks on eight supplied
acyclic families; caller mutation is tested locally and is not a cross-library
mutation-equivalence claim. Full evidence and reproduced differences are in
`validation/results/tautomer-score-callback-2026-10-18/`.
