# Aromatic ring triple and same-coordinate USR follow-up

Source candidate based on `origin/main` `7704ee2030dda8f9ca54a07eb3fc22e94623fb49`.
Comparator: RDKit 2026.09.1 in Python 3.13.6 on macOS 27.0.1 arm64.
The [machine fingerprint packet](2026-10-25-rdkit-core-usr-accuracy.json)
records corpus hashes, exact per-algorithm counts, and every mismatch (none).
The measured CLI binary SHA-256 was
`61cc65c0164d4bf8e3e7c34f7f12f7c66534840fb60d2265e16978a9e6ce81a5`;
the comparator script SHA-256 was
`efa8850f2b324a7bad6a4fb1c74fe8344370b572c0dc2dfa868fcb183ca99cd0`.

## Chemical graph, stereo search, and fingerprints

The one remaining stress-row mismatch from the previous packet contained a
triple bond in an aromatic ring. RDKit accepted the graph with explicit
valence 4 at the affected carbon. CheMatic's RDKit-compatible sanitizer added
an extra floor-sum check and rejected it. Removing that non-RDKit check lets
the common sanitized graph reach all five fingerprint implementations.
The separate non-aromatic valence check and a six-bond aromatic rejection
remain in the regression tests.

| Corpus | RDKit-accepted rows | Exact per fingerprint |
| --- | ---: | ---: |
| Stress v1 | 4,072 | 4,072/4,072 |
| Browser 10k | 10,000 | 10,000/10,000 |
| Descriptor census | 5,000 | 5,000/5,000 |

Each row was compared for Morgan R2, Pattern, RDK, Topological Torsion, and
Atom Pair 2048-bit fingerprints: 95,360 exact result cells. The same gate
found zero mismatches on 324 stereo-search pairs and 1,653 stereo-bearing
self-searches. These are fixed corpora, not a proof of all RDKit graph or
fingerprint behavior. The previous one-row gap was resolved; no RDKit
refusals were removed from the denominator.

## 3D USR on identical coordinates

The existing `usr_descriptors` returns unstandardized second and third
central moments and uses coordinate-value tie breaks. RDKit's `GetUSR`
uses standard deviation, the signed cube root of standardized skewness,
and first-index tie breaks. A separate `rdkit_usr_descriptors` API now
implements those semantics on caller-supplied coordinates, including explicit
hydrogens. It rejects missing or nonfinite geometry with typed errors.
The Python binding exposes the same calculation as
`chematic.rdkit_usr_descriptors(coords)`.

The [RDKit oracle](../validation/rdkit_usr_3d_2026_09_1.json) was generated
from the previously frozen explicit-H coordinates in
`validation/rdkit_3d_g1_fixtures.json` without reembedding. It contains 24
conformers from 22 molecules, including a macrocycle, halogens, sulfur and
phosphorus. All 288 moments agree within `1e-9`; the maximum absolute
difference is `1.352e-14`. The legacy profile differs by more than `1e-6`
on all 24 conformers. The oracle records the source fixture SHA-256 and RDKit
version. The generator is
[`gen_rdkit_usr_3d_reference.py`](../scripts/generators/gen_rdkit_usr_3d_reference.py).
The Python wheel passed both binding tests; a direct four-point comparison
against RDKit 2026.09.1 had zero absolute difference across all 12 values.

## Reproduce

```sh
cargo build --locked -p chematic-cli
python scripts/validation/rdkit_cli_stereo_fingerprint_gate.py \
  --binary target/debug/chematic \
  --stereo-self-corpus scripts/fixtures/descriptor_census_corpus.smi \
  --fingerprint-corpus validation/cosmolkit_comparison/corpora/stress-v1.smi \
  --fingerprint-corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --fingerprint-corpus scripts/fixtures/descriptor_census_corpus.smi
python scripts/generators/gen_rdkit_usr_3d_reference.py \
  --output /tmp/rdkit-usr-oracle.json
cargo test --locked -p chematic-3d --test rdkit_usr_parity
# With chematic installed from this source into the active Python environment:
python -m pytest -q crates/chematic-py/tests/test_rdkit_usr.py
```

The fingerprint comparison concerns bit exactness, not molecule recognition
outside these inputs. The USR comparison fixes the coordinates on both sides;
it does not compare conformer generation, force-field physics, convergence,
diversity, or experimental geometry. Neither result establishes broad RDKit
superiority or performance superiority.
