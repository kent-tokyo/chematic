# MMFF94 finite-difference gradient parallel lane — 2026-09-05

Status: local release-mode evidence for the current `1.0.6` source tree.
No version, tag, or published-artifact claim is made here.

Independent atom coordinate probes in the prepared L-BFGS gradient are now
parallel for molecules with at least 16 atoms. Molecules below that threshold
retain the sequential, low-allocation path. Each atom writes only its own
three gradient components, so no floating-point reduction order changes.

Environment: Apple arm64 host, release profile, Criterion `--quick`, Plotters
backend. The new large-molecule lane uses a 24-carbon chain and two L-BFGS
iterations:

| Benchmark | Criterion interval |
| --- | ---: |
| `mmff94_lbfgs_large_24atom_2iter` | 5.4491–9.1009 ms |
| Existing `mmff94_lbfgs_6mol_8iter` | 46.346–46.420 ms |

Criterion also reported the existing lane's cached-baseline change as
non-significant (`p = 0.05`); this is not treated as a causal A/B speedup
claim because the baseline was not captured in a fresh paired run. The result
is retained as a reproducible smoke/performance lane, not a headline
multiplier.

Validation:

```text
CARGO_TARGET_DIR=/private/tmp/chematic-binding-target \
  cargo test -p chematic-ff --lib --offline
```

Result: 202 passed, 0 failed. Analytic gradients and coordinate-dependent
bounded neighbor lists remain separate experimental gates.
