# Zero-seed ensemble regression

RDKit [#9738](https://github.com/rdkit/rdkit/issues/9738) was OPEN when checked
on 2026-10-24 JST. The multiplicative ensemble schedule also made Chematic
attempt the same seed zero repeatedly. `native-before.txt` records the failing
internal-distance regression on base `38b4c92c430d9b3bfe70e47e97293b9017ced990`.

The fix uses seeds `0, 2, 3, ...` for a zero base seed. The first attempt
preserves existing single embedding, and one is skipped because the underlying
generator aliases zero and one. Positive seeds keep their old multiplicative
schedule. Corrected-gradient, MMFF and extra quality-profile batches use the
same seed derivation. Signed seed/index overflow is checked before constructing
the ensemble; failure omission and candidate order remain deterministic.

## Fixed diagnostic

Four authored flexible inputs: `CCCCCCO`, `CCCCCC`, `NCCCCCC`, `CC(C)CCCO`.
Each runs seeds 0/1/42, five conformers, no pruning or MMFF, repeats the full
run, and compares the three-candidate prefix. Each graph has explicit Hs.

| Runtime | Seed 0: varying internal geometry | Seeds 1 and 42: varying internal geometry | Repetition / prefix |
|---|---|---|---|
| Prior local Chematic runtime | 0/4 | 8/8 | 12/12 each |
| RDKit 2026.09.1 | 0/4 | 8/8 | 12/12 each |
| Patched Chematic runtime | 4/4 | 8/8 | 12/12 each |

Variation means that at least one conformer pair has a heavy-atom interatomic
distance change greater than 0.1 angstrom. Translation, rotation and reflection
alone cannot produce this change. It is not a statement about optimized minima,
energy, ensemble coverage, pruning quality, performance or toolkit-wide quality.
Distinct seeds do not guarantee distinct conformations for arbitrary molecules.

The prior runtime is the previously built `86afa09f` candidate (base main plus
the unrelated isotope-formula addition), not a separately rebuilt main/release
binary. All 12 first conformers and all eight positive-seed ensembles compare
exactly with that runtime on this host. The native red/green regression provides
the direct main-source reproduction. RDKit coordinates are retained as an
independent within-engine diagnostic; cross-engine coordinate equality is not
asserted. `summary.json` pins source, runtime-packet and observer hashes; the
gzip packets include every coordinate array and per-pair distance metric.

## Reproduce

Build/install the candidate Python extension, then run with its import root:

```sh
python3 scripts/validation/validate_zero_seed_ensemble.py \
  --engine chematic --python-root /path/to/candidate/python-root \
  --label candidate --output /tmp/zero-seed-candidate.json.gz
/path/to/rdkit-2026.09.1/bin/python scripts/validation/validate_zero_seed_ensemble.py \
  --engine rdkit --label rdkit-2026.09.1 --output /tmp/zero-seed-rdkit.json.gz
cargo test --locked -p chematic-3d --test zero_seed_ensemble --test rdkit_embed
cargo test --locked -p chematic-3d --lib rdkit_embed::
PYTHONPATH=/path/to/candidate/python-root PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 \
  python3 -m pytest crates/chematic-py/tests/test_zero_seed_ensemble.py \
  crates/chematic-py/tests/test_rdkit_embed.py
```

The diagnostic has fixed inputs and candidate counts. It has no timing or
process-memory claim. WASM single embedding is unchanged; the shared Rust crate
was also compile-checked for `wasm32-unknown-unknown`.
