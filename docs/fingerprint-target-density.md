# Fingerprint folding to a target density

This isolated, unreleased candidate implements the correction proposed in
[RDKit #9715](https://github.com/rdkit/rdkit/issues/9715), verified OPEN on
2026-10-15 (updated 2026-10-15T08:57:28Z). RDKit 2026.09.1 still folds when
the fingerprint density is **above** its target, and refuses a next width
exactly equal to `minLength`. The reported correction is to fold while
density is **below** the target and the next width is **at least** the minimum.
chematic had no target-density folding API; this adds one in Rust, Python and
the Python RDKit adapter.

## Example

```python
import chematic

# 16 on bits in a 128-bit packed fingerprint: density 0.125.
raw = b"\xff\xff" + bytes(14)
result = chematic.fold_to_target_density(raw, density=0.25, min_length=64)
assert result["bit_width"] == 64
assert result["fingerprint"] == b"\xff\xff" + bytes(6)
assert result["density"] == 0.25
assert result["folds"] == 1
assert result["target_reached"]

# Molecular fingerprints use the same packed-byte entry point.
folded = chematic.fold_to_target_density(chematic.from_smiles("CCO").ecfp4())
```

Python returns a `FoldDensityResult` dictionary:

| Key | Meaning |
| --- | --- |
| `fingerprint` | Folded LSB-first packed `bytes` |
| `bit_width` | Actual final width, including non-byte-aligned widths |
| `original_bits` | Validated input width |
| `on_bit_count` | Number of set bits after folding |
| `density` | `on_bit_count / bit_width` |
| `folds` | Number of halving steps |
| `target_reached` | Whether final density is at least the requested target |
| `stop_reason` | `target_reached` or `minimum_length` |

Rust exports `fold_to_target_density(&BitVecN, density, min_length)` and
`fold_packed_to_target_density(&[u8], bit_width, density, min_length)` from
`chematic_fp`. Both return `Result<FoldDensityResult, FoldDensityError>`.
The result contains the native `BitVecN`, original width, fold count and target
status; its `density()` method computes the achieved density.

## Folding and validation contract

Every step uses `new_width = floor(old_width / 2)` and sets output bit
`old_index % new_width` for every input on bit. Colliding bits combine with OR.
The existing `BitVec2048::fold` retains its published XOR behavior; the new
operation has a separate entry point.

- Fold only while density is strictly below the target and the next width is
  at least `min_length`. Equality at the target stops; equality at the minimum
  allows the fold.
- Odd widths are supported, including tail bits. With odd widths, collisions
  can reduce density: bits 0/2/4 at width 5 fold to bit 0 at width 2.
- All-zero fingerprints can stop at the minimum with an unmet positive target.
  A target of zero needs no folding. An input already shorter than the minimum
  is returned at its original width. No operation enlarges a fingerprint.
- The source is unchanged. Every result, including a no-op, is an independent
  copy. The returned status discloses whether the target was achieved.
- The width must be positive and at most **1048576 bits** (128 KiB packed).
  The cap is checked before input decoding or copying. Each step reduces the
  width, so at most 20 folds are possible. No unlimited-width option is exposed.
- `density` must be finite in `[0, 1]`; `min_length` must be positive. A minimum
  larger than the input is valid and prevents folding.
- Python takes `bytes`. `bit_width` is keyword-only and defaults to
  `len(fp) * 8`. Packed input must have exactly `ceil(bit_width / 8)` bytes and
  zero unused high padding bits. This permits exact widths such as 5 or 65.
- Rust returns typed errors; Python raises `ValueError` for contract failures,
  `OverflowError` for negative/out-of-range unsigned integer options, and
  `TypeError` for non-`bytes` packed input.

The width cap bounds this operation's work and buffer sizes. Existing
`BitVecN` allocations are used; allocator failure is not converted into a
new fallible-allocation contract. Fingerprint generation, database-wide
folding policies, similarity recalibration and WASM bindings are outside this
candidate. Differently sized packed results require an explicit comparison
policy before similarity calculation.

## RDKit adapter

```python
from chematic.rdkit_compat import DataStructs, ExplicitBitVect

fp = ExplicitBitVect(128)
for bit in range(16):
    fp.SetBit(bit)
folded = DataStructs.FoldToTargetDensity(fp, density=0.25, minLength=64)
assert len(folded) == 64
assert folded.GetNumOnBits() == 16
assert len(fp) == 128
```

The adapter accepts its own `ExplicitBitVect` and returns a new one. It
delegates to the native operation; `__len__` and `GetNumOnBits` are added to
the adapter vector. `SparseBitVect` and RDKit-native vector objects are outside
the adapter contract. Values use corrected #9715 semantics, so they deliberately
differ from the installed buggy RDKit wrapper. RDKit can return the original
object when no folding occurs; this adapter always returns an independent copy.

There is an additional upstream primitive boundary: RDKit 2026.09.1
`FoldFingerprint(fp, 2)` rejects a **2-bit** input with
`ValueError: invalid fold factor`, for both explicit and sparse vectors.
Its 3-to-1 fold works. chematic supports both boundaries with modulo-OR
semantics. The validator records this distinction separately from #9715.

## Evidence

- Rust: **356 passed** — 346 existing library tests and 10 new integration
  tests, zero failures or ignored tests.
- Python: **182 passed** — 20 new folding tests and 162 existing adapter tests.
- Independent modular-set oracle: **20259 exact cases**, including every bit
  pattern for widths 1–8, six density thresholds, six minimum lengths, plus
  byte/word/odd/fingerprint boundaries through width 4096.
- RDKit **2026.09.1** explicit and sparse vector primitives: **40262 complete
  exact comparisons** and **14744 exact successful fold steps**. Another
  **256 comparisons** stop at the known 2-to-1 primitive exception; those
  terminal cases are covered by the independent oracle. Unexpected mismatches:
  **zero**. The RDKit oracle uses the corrected loop around `FoldFingerprint`,
  not the installed buggy `FoldToTargetDensity` wrapper.
- Adapter results match all **20259 cases**, with input immutability checked.
- Seven of nine original, minimal regression recipes differ from the installed
  RDKit wrapper, reproducing the reported direction/minimum issues. Two recipes
  are unchanged controls. Exact maximum width and first rejected width are
  exercised in Rust/Python boundary tests.

The recipe manifest is `validation/rdkit_issues/fold-density/fixtures.json`.
Each small fixture records provenance, license, redistribution, minimization,
pinned oracle, expected units/tolerance and normalized payload SHA-256.
Generated case IDs and their compact normalized inputs/outputs feed a packet
SHA-256; the checked-in recipe regenerates that packet without external data.
`scripts/validate_fold_density.py` saves counts, the installed wrapper's source
hash, upstream repro results and packet hash in
`validation/results/rdkit-fold-density-2026-10-15.json`.

```sh
cargo test -p chematic-fp --lib --test fold_density --offline --quiet
cargo clippy -p chematic-fp -p chematic-py --all-targets --offline --no-deps -- -D warnings
python3 -m pytest crates/chematic-py/tests/test_fold_density.py crates/chematic-py/tests/test_rdkit_compat.py -q
python3 scripts/validate_fold_density.py
```

Formatting, static binding inventory, refactor budget, Python stub syntax,
product naming and strict MkDocs checks also pass. These are local source
checks with a debug Python extension, without native InChI; no publication
or upstream modification is claimed.

The source candidate is `/private/tmp/chematic-fold-density`, inherited from
`/private/tmp/chematic-native-matrix`. The incremental patch is
`validation/rdkit_issues/fold-density-2026-10-15.patch`, applicable to that base.
Only this document, the patch and validation record are copied to the primary
checkout. Source changes remain isolated, uncommitted and unpushed. The sparse
local baseline commit anchors the diff and is not an integration commit.
