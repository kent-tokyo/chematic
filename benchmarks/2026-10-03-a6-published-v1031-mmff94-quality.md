# Published v1.0.31 MMFF94 stereo-safe quality rerun

The published PyPI `chematic==1.0.31` CPython 3.13/macOS arm64 wheel
(SHA-256 `b481316bb5a49cb572b035c805316e38a381f4f0e6fa141a05d1da3dc6f05c1f`)
was installed in an isolated environment on macOS 27.0.1 arm64. The existing
Tier A+B 265-molecule manifest and `scripts/public_package_3d_chematic.py`
were used with seed 20260801, eight attempts, 200 force-field iterations,
20-second per-row limit, diagnostic-only ring torsions, and the
`chematic_pipeline_v2_mmff94_strict_stereo_safe` arm. The runner completed
265/265 rows; no failures or fallback occurred.

| Dimension | v1.0.31 result |
|---|---:|
| Bond/angle/OOP/stretch-bend terms missing | 0 each |
| Torsion terms missing | 32, across six rows (permitted by this strict bond-angle arm) |
| Force-field convergence | 100/265; **165/265 not converged** |
| Independent geometry soundness | 265/265 |
| Independent stereo clean | 265/265 |
| Independent gross-clash count | 0 |

The independent scorer used the same source bytes as the v1.0.26 packet
(SHA-256 `2f631a1368cade9c588fc9f19490bfa26d326519c780b37ead15d92616c6bbb9`).
It scored all 265 v1.0.31 successes and the 264 successful rows from the
archived, pinned RDKit 2026.03.6 run (archive SHA-256
`a856235a39f37e9f8385ff4d60fc4efc8fdf116eedb43effaf4e3893b041dd1f`).
There were zero integrity-error rows; RDKit's archived row 8 remains an
`internal_error`. Reusing the archived RDKit coordinates is an output-quality
comparison, **not** a new paired speed run.

The v1.0.31 coordinates, convergence flags and per-row parameter coverage
are exactly equal to all 265 corresponding v1.0.26 published rows. Canonical
coordinate-array SHA-256 is
`f5fb982a8159b5a3625592bf4194bcf4969e78ecadccaf48ccfc00cbe1685a8e`
for both releases. The v1.0.31 raw runner JSONL SHA-256 is
`c2525cf7a12f7fd77b1881ee868b49785b7b62b44e006eb708d2c630b6e0d850`;
the scored JSONL SHA-256 is
`73495316864f2f924dde0475adacb27945269cb3888b5595ff23f1a136fa6189`.
Those raw files were kept outside Git; they can be regenerated with the
version-pinned wheel and the repository runners. The [machine summary](../validation/results/a6-published-v1031-mmff94-quality.json)
retains the digests and counts.

## Same explicit-H coordinates, separate from 3D generation

`scripts/mmff94_same_explicit_h_energy.py` was also run with the same
published v1.0.31 wheel and pinned RDKit 2026.03.6. RDKit generated a single
explicit-H coordinate set per row; both engines evaluated that exact set.
Statuses were 262 comparable, two RDKit embedding failures and one declared
unsupported row. All 262 coordinate hashes and statuses matched the
historical v1.0.19 packet, so the energy change is not explained by a changed
input geometry. The previous >5 kcal/mol rows 166 and 231 are now 0.004 and
0.048 kcal/mol apart; all 262 comparable total energies are within 1 kcal/mol
(maximum 0.262; median 0.0043). The raw output SHA-256 is
`6ede7c998827d97bca1436246d0b5290e8e1ecf6e7b25991a44d2f954e4d7054`.

Per-term residuals are still visible: three angle rows and one OOP row exceed
0.1 kcal/mol. On the two formerly large-residual rows, the current wheel's
analytic-gradient/central-difference check has maximum scaled error
4.37×10⁻⁸. This closes the old **total-energy residual** on the comparable
cohort, not full atom typing, every per-term tolerance, the three uncomparable
rows, convergence, or conformer quality.

The energy command used the isolated published-wheel Python environment with
RDKit 2026.03.6 on its import path:

```sh
python3 scripts/mmff94_same_explicit_h_energy.py \
  --output /tmp/mmff94-same-coordinates.jsonl --seed 20260913 \
  --gradient-input-index 166 --gradient-input-index 231 \
  --gradient-delta 1e-5
```

The runner invocation (in an isolated environment containing the wheel) was:

```sh
python3 scripts/public_package_3d_chematic.py \
  --wheel /path/to/chematic-1.0.31-cp313-cp313-macosx_11_0_arm64.whl \
  --tiers AB --arms chematic_pipeline_v2_mmff94_strict_stereo_safe \
  --output /tmp/chematic-3d.jsonl \
  --metadata-output /tmp/chematic-3d.meta.json \
  --wall-budget-seconds 1800
cargo run -p chematic-3d --example pipeline_v2_vs_rdkit_common_scorer -- \
  --chematic-rows /tmp/chematic-3d.jsonl \
  --rdkit-rows validation/results/mmff94-public-v1026-20260926/rdkit-3d.jsonl \
  --pair chematic_pipeline_v2_mmff94_strict_stereo_safe rdkit_etkdgv3_mmff94 \
  --output /tmp/chematic-common-scored.jsonl
```

This closes only the **published v1.0.31 geometry/stereo/clash** recheck.
Convergence, full force-field typing, per-term energy boundaries,
and independent conformer-quality non-inferiority remain open. It does not
establish MMFF94 parity or a speed advantage.
