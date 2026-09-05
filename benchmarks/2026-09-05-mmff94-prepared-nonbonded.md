# MMFF94 prepared nonbonded terms — 2026-09-05

Status: local release-mode evidence for the current `1.0.6` source tree.
No version, tag, or published-artifact claim is made here.

The prepared MMFF94 model now resolves van der Waals combination parameters
and electrostatic charge products once during topology preparation. Repeated
energy evaluations, including the finite-difference gradient probes used by
L-BFGS, no longer repeat those topology-only lookups. Pair exclusion and 1-4
scaling semantics are unchanged.

Environment: Apple arm64 host, `cargo bench -p chematic-ff --bench
mmff94_bench --offline -- --quick`, release profile, six fixed SMILES fixtures.
Criterion used the Plotters backend because Gnuplot was unavailable.

| Benchmark | Median-ish Criterion interval |
| --- | ---: |
| Prepared energy, 6 molecules | 105.69–108.64 µs |
| One-shot energy, 6 molecules | 2.9327–3.3038 ms |
| L-BFGS, 8 iterations, 6 molecules | 158.98–179.25 ms |

The prepared/one-shot values are not a like-for-like speedup claim: one-shot
includes model construction while prepared energy reuses the model. The
purpose of this lane is to pin the optimized path and guard output parity.

Validation:

```text
CARGO_TARGET_DIR=/private/tmp/chematic-binding-target \
  cargo test -p chematic-ff --lib --offline
```

Result: 202 passed, 0 failed. The broader finite-difference-to-analytic
gradient replacement and coordinate-dependent bounded neighbor list remain
separate experimental work and are not marked complete by this change.
