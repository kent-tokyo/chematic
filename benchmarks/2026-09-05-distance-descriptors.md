# Shared distance descriptor context

Implemented 2026-09-05 in the local v1.0.6 working tree. This record covers
semantic reuse and parity, not a published-artifact performance claim.

`DistanceDescriptorBundle` computes AutoCorr2D, Moran, and Geary from one
topological distance matrix. The existing `autocorr_2d`, `moran_autocorr`, and
`geary_autocorr` functions remain lazy and continue to compute only their own
family when called independently.

Parity was checked for ethane, benzene, and aspirin. The `chematic-chem` suite
passed **842 tests** with one pre-existing ignored test. Python binding
`cargo check --locked`, clippy, formatting, and diff checks also passed.

BCUT2D is intentionally not routed through this bundle: it uses a Burden
connectivity matrix and eigenvalue calculation rather than topological shortest
path distances. It remains a separate optimization target.
