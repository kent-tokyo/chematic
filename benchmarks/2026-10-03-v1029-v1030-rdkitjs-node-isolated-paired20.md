# Published npm/WASM paired Node comparison — 2026-10-03

This is an exposed, single-host Node measurement of the **published**
`@kent-tokyo/chematic` v1.0.29 and v1.0.30 tarballs and official
`@rdkit/rdkit` 2026.03.6. It is not a browser benchmark or an all-operation
speed claim. [Raw 20-block record](2026-10-03-v1029-v1030-rdkitjs-node-isolated-paired20.json)
SHA-256:
`edb2b289698dc022db9791b0307f6cb001b141d961045324a524b60e0882ae37`.
The [runner](../scripts/bench_published_wasm_paired.mjs) and
[offline checker](../scripts/check_published_wasm_paired.py) pin the artifact
hashes and verify every paired block, outcome digest, interval and RSS value.

## Contract

- Input: the first 250 of 5,000 rows in
  [`descriptor_census_corpus.smi`](../scripts/descriptor_census_corpus.smi)
  (full-file SHA-256
  `d6f2ba3f128296f935007f0b0813aa97b6ebcc2457e014ddca2213ddd655276c`).
  Every row parsed in all three implementations; atom counts and radius-2,
  2,048-bit Morgan fingerprints were **250/250 bit-identical** before timing.
  The timed fingerprint outputs are checked against the same digest.
- Published tarball SHA-256: CheMatic v1.0.29
  `d0af78a8d6b711a13b63985d4b36079e245441e99f1712532555dc69f11569fa`;
  v1.0.30 `fd427a161c11b14e6cca51e78e0c35ff05bbc52faf0539a2ddc2b7090e534201`;
  RDKit.js `3b86b72775394ae997fb96ca854c6e7c5840927d3e899c1ef117d42bca573c87`.
  The extracted JS and WASM file hashes and loaded versions are checked too.
- Environment: Apple M4, macOS arm64, Node v24.5.0. Each observation runs in
  its own Node process after 20 warm-up rows. Each of 20 blocks uses ABBA or
  BAAB order (alternating), two observations per side. The 95% interval is
  10,000 deterministic bootstrap resamples of paired block log-ratios.
  A practical win requires the **lower** interval bound to exceed 1.10.
- `parse_only` times SMILES parse and object release. `parse_morgan`
  includes parse, Morgan generation and release. `prepared_first_use`
  parses before timing and then computes each fingerprint once.
  `prepared_reused` also primes each molecule's fingerprint before timing.
  Import/WASM initialization is measured separately. Peak RSS includes the
  entire fresh Node process; it does not isolate library allocations.

Ratio above 1 means the left-hand implementation is faster than the right.

| Timed lane | v1.0.29 / v1.0.30 (95% CI) | v1.0.30 / RDKit.js (95% CI) |
|---|---:|---:|
| Parse only | 1.01 (0.99–1.03) | 28.59 (28.17–29.01) |
| Parse + Morgan | 1.00 (0.98–1.02) | 6.09 (5.97–6.18) |
| Prepared first use | 0.99 (0.97–1.01) | 1.35 (1.31–1.39) |
| Prepared reused | 0.99 (0.97–1.01) | 1.76 (1.71–1.80) |

The two CheMatic releases are indistinguishable here under the declared
10% rule. CheMatic v1.0.30 clears that rule against RDKit.js on all four
**Node** lanes for these 250 fingerprint-equivalent rows. The result does not
establish faster canonical writing, all 63 operations, browser startup, or
general RDKit chemistry parity.

As separate whole-process observations, parse-only median peak RSS was
67.0 MiB (CheMatic) versus 145.3 MiB (RDKit.js); prepared-reused was
79.1 versus 151.2 MiB. Median import-plus-WASM initialization was about
14 versus 83 ms. These figures include runtime, module loading, and allocator
effects, so they are not a matched per-molecule memory or browser metric.
The raw JSON retains all 640 child-process timings, starts and RSS values.

## Reproduce

Fetch and extract the three named npm tarballs into separate directories,
verify the hashes above, then run:

```sh
node scripts/bench_published_wasm_paired.mjs \
  --corpus scripts/descriptor_census_corpus.smi --rows 250 \
  --v29-dir OUT/v29/package --v29-tarball OUT/v29/kent-tokyo-chematic-1.0.29.tgz \
  --v30-dir OUT/v30/package --v30-tarball OUT/v30/kent-tokyo-chematic-1.0.30.tgz \
  --rdkit-dir OUT/rdkit/package --rdkit-tarball OUT/rdkit/rdkit-rdkit-2026.3.6.tgz \
  --output OUT/paired.json
python3 scripts/check_published_wasm_paired.py
```

The checker validates the checked-in record, not a newly generated
`OUT/paired.json`. A different machine should retain its own raw output and
must not overwrite this M4 record.
