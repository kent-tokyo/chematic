# Descriptor topology context

Implemented 2026-09-05 in the local v1.0.6 working tree. This record covers
semantic reuse and parity, not a published-artifact performance claim.

`chematic_chem::TopologyBundle` now shares the heavy-atom list and membership
set across Wiener, Kappa 1–3, and Chi 0–4/0v–4v calculations. The full Python
bulk descriptor path uses the bundle. `descriptors_array` uses it when at least
two of those topology groups are requested, while preserving lazy scalar
calculation for a single selected group.

The bundle was checked against the existing scalar APIs for ethane, propane,
benzene, and aspirin. The `chematic-chem` library suite passed **841 tests**
with one pre-existing ignored test; Python binding compilation passed with
`--locked`, and clippy/fmt/diff checks passed.

This slice deliberately does not claim that all distance-matrix descriptors
share one matrix yet. AutoCorr, Moran/Geary, BCUT, and related families remain
separate follow-up work because their atom weighting and disconnected-graph
contracts differ.
