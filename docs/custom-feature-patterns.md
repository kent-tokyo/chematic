# Custom SMARTS feature patterns

`FeaturePatternSet` compiles and owns an ordered iterable of SMARTS strings,
then reuses it for atom-feature masks and native circular fingerprints:

```python
import chematic

features = chematic.FeaturePatternSet(p for p in ["[#8]", "C=O", "[#6]"])
mol = chematic.from_smiles("CC(=O)O")
assert features.invariants(mol) == [4, 6, 3, 1]
packed = features.fingerprint(mol, radius=2, nbits=2048)
```

Pattern `i` sets bit `i` on **every** atom participating in any match. Pattern
order and duplicates matter. Up to 32 patterns are supported; all 32 bits are
preserved, so 32 copies of `O` assign `0xffffffff` to matching oxygens. Empty
definitions produce zero masks rather than silently selecting default features.
An unmatched atom also receives zero. `patterns` returns a separate list of the
original definitions; mutating it cannot change the compiled set.

The Python object owns Rust query graphs and source strings. Lists, tuples,
generators and other iterables work; no Python query pointers are retained.
Construction consumes at most 33 items, refusing the 33rd. Iteration exceptions
propagate unchanged, and strings/bytes passed as the collection or non-string
elements raise `TypeError`. Invalid SMARTS and size violations raise `ValueError`.
Rust callers use `chematic_fp::feature_patterns::FeaturePatternSet::new`.

Each pattern is limited to 4096 UTF-8 bytes and 32 top-level query atoms.
The parser also applies its existing recursive-SMARTS depth limit. Targets
are limited to 256 atoms / 512 bonds. Matching uses the RDKit-parity aromatic
view with source atom indices, native SSSR SMARTS ring semantics, isotope
constraints and chirality ignored. This is not a claim that every native SMARTS
predicate matches every RDKit profile.

`max_visits` defaults to 100,000 **per pattern**, including recursive matching,
and must be 1..1,000,000. Exhaustion raises without returning partial masks or
a fingerprint; the same object can subsequently be used with a larger budget.
Arithmetic always uses unsigned 32-bit feature masks; more than 32 definitions
are refused before any bit shift.

The packed fingerprint uses chematic's FNV circular hashing, with all four
little-endian bytes of each mask as the initial seed. It emits every atom at
every radius, including zero-mask atoms. Radius 0..3 and widths 128, 256, 512,
1024 and 2048 are supported. Aromatic/Kekulé forms use the same perceived
bond view. The result is LSB-first immutable `bytes` of length `nbits / 8`.
These are **not RDKit Morgan bit identifiers**, nor the existing default FCFP
algorithm's one-byte seeds. Existing FCFP/Morgan methods and compatibility
wrapper options retain their current behavior.

## Upstream report and validation

[RDKit #9601](https://github.com/rdkit/rdkit/issues/9601) tracks wrapper input
contracts. Its proposed [PR #9615](https://github.com/rdkit/rdkit/pull/9615)
reports custom Morgan feature patterns being ignored and unsafe handling above
32 patterns. The feature-mask semantics are documented by the official
[FingerprintUtil implementation](https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/Fingerprints/FingerprintUtil.cpp).

An original probe on RDKit 2025.09.3 reproduced the ignored definitions:
both `O` and `C` custom patterns yielded the default atom masks `[0, 0, 3]`
for ethanol, rather than the requested `[0, 0, 1]` and `[1, 1, 0]`. The installed
Boost wrapper also accepted 33 definitions while ignoring them; no C++ overflow
or unsafe native call was executed. Before this candidate, chematic had no
custom-pattern set API. This adds an explicit native API, not a complete
`GetMorganFeatureAtomInvGen` compatibility implementation or upstream repair.

`validation/rdkit_issues/custom-feature-patterns/probe.py` compares atom masks
with an independent OR of RDKit SMARTS matches across 150 cases: 15 molecules,
nine prefix sets (including empty) and a 32-pattern set. Controls include
multi-atom matches, isotope and recursive patterns, monocyclic rings, and
aromatic/Kekulé spellings. Fingerprint tests check caller-pattern influence,
upper-bit preservation, supported shapes, representation invariance and errors;
no RDKit fingerprint-bit parity is claimed.

Dated source/artifact hashes and logs are under
`validation/results/custom-feature-patterns-2026-10-17/`. These are local source
checks, not published-wheel or completed hosted-CI evidence.
Local validation passed 352 Rust tests (six new), 211 Python tests (49 new),
and all 150 direct SMARTS-mask comparisons. Strict Clippy, format, source stubs,
refactor budgets, binding/manifests, naming, unsafe-surface and strict MkDocs
checks passed.
