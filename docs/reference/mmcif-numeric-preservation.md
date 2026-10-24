# mmCIF numeric preservation

`write_mmcif` now retains the finite binary64 values of Cartesian coordinates,
occupancy, isotropic B factors and all six cell parameters. It emits enough
decimal digits for an exact read/write/read round trip, using scientific
notation when a plain decimal would exceed 24 characters.

Previously, the writer rounded coordinates to three decimal places,
occupancy/B factors to two, cell lengths to four and cell angles to three.
For example, `x=1.123456789` became `1.123` and `occupancy=0.3333333333333333`
became `0.33`. Rust, Python and WASM/Node now share the corrected writer.
Signatures and record shapes are unchanged. Output text has different decimal
spelling, so consumers should parse CIF numbers rather than require fixed
decimal widths.

The wwPDB dictionary defines
[Cartesian x](https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_atom_site.Cartn_x.html)
as a floating point value in Angstroms; mmCIF is not constrained by PDB's
fixed-width coordinate fields. This implementation retains the *parsed numeric
value*. Original decimal spelling, significant-figure intent, uncertainty
suffixes, comments and complete source bytes are outside this writer's contract.

## Development comparison

The packet has three controls, twelve reproduction cases, and four separate
arithmetic storage cases. All are original synthetic MIT fixtures. Source fields
are checked against an authored numeric oracle by an independent Python observer
before any engine is scored. The observer checks ordered atom IDs/elements,
Cartesian coordinates, occupancy, B factors and cell parameters. It requires
explicit numeric columns and compares binary64 bits; absent columns are losses
even when a parser's default would happen to match the source value.

The twelve reproduction cases have the following exact-field results. Python
and WASM have identical counts; native Rust tests independently verify the
same packet and the existing modeled atom records.

| Observed field | Baseline chematic | Candidate chematic | Open Babel 3.2.1 |
|---|---:|---:|---:|
| Ordered atom IDs/elements | 12/12 | 12/12 | 12/12 |
| Cartesian coordinates | 4/12 | 12/12 | 4/12 |
| Occupancy | 9/12 | 12/12 | 0/12 |
| Isotropic B factor | 9/12 | 12/12 | 0/12 |
| Cell parameters, when present | 0/3 | 3/3 | 0/3 |
| All declared numeric fields together | 0/12 | 12/12 | 0/12 |

The three controls retain all declared numeric fields in baseline and candidate
chematic. Open Babel retains their coordinates and the one control cell, but
omits occupancy/B-factor columns, so its combined result is 0/3. In the four
separate arithmetic storage cases, combined numeric preservation is 1/4 for the
baseline, 4/4 for the candidate and 0/4 for Open Babel. Every engine returns an
output for every case under the final process budget; refusals or resource
failures are not being counted as numeric losses in these counts.

Open Babel's coordinate/cell failures here contain rounded values; its occupancy
and B-factor failures are missing columns. This is a measured advantage for
these explicit numeric fields in this development packet, not a general
mmCIF or Open Babel ranking.

The Open Babel lane uses `obabel -immcif -ommcif -ab` on the same original bytes.
Its documented [mmCIF read option](https://openbabel.org/docs/FileFormats/Macromolecular_Crystallographic_Info.html)
disables bond inference, matching chematic's atom-record workflow. Each process
is bounded to 60 seconds and 1 MiB combined stdout/stderr; resource limits,
refusals and no-output results remain explicit. No timing comparison is made.

The baseline Python/WASM artifacts were built from
`342d1d5f9ad09aa6ade4e3079a902cba856209e2`. The relevant mmCIF reader/writer and
binding sources match main base `fcbf6d0a`. Records retain source/runtime hashes,
every emitted CIF, errors and field observations. Their recorded HEAD identifies
the pre-commit checkout; candidate source hashes identify the exact changed
files. The independent observer supports this packet's single-block, scalar
and flat-loop profile; it is not a complete CIF or dictionary validator.

Native tests also cover signed zero and finite binary64 extrema, checking that
numeric magnitude cannot expand one small atom row beyond 256 bytes. The public
Python/WASM packet includes four storage cases separately from the twelve
reproductions. Those extremes do not establish physical or deposition validity.
The JavaScript packet contains no negative zero: native Rust tests cover that
case, while a caller's `JSON.stringify` can normalize it before a WASM call.

## Scope and reproduction

This result concerns finite numeric fields in the existing parsed atom-record
profile. The existing `label_*`/`auth_*` collapse, first-data-block boundary,
unhandled categories and unknown-column behavior are unchanged. Use the
[mmCIF chain record API](mmcif-chain-records.md) when retaining source and complete
identity namespaces matters. No result here establishes complete mmCIF metadata,
uncertainties, bonds, stereo, schema validation, CLI support or speed superiority.

Build a Python package and a Node-target WASM module from the selected checkout,
then run:

```sh
PYTHONPATH=/absolute/python-package-root python3 scripts/validation/check_mmcif_numeric_preservation.py \
  --wasm-module /absolute/pkg-node/chematic_wasm.js --output /tmp/mmcif-numeric.json
python3 -m unittest scripts.tests.test_mmcif_numeric_observer -q
cargo test -p chematic-mol --test mmcif_numeric_preservation
```

Full fixtures, before/after outputs and identities are in the
[numeric packet](https://github.com/kent-tokyo/chematic/tree/main/validation/mmcif_numeric_preservation).
Native, Python and Node regression suites run in CI, with seven independent
observer negative-control tests that need no chemistry package.
