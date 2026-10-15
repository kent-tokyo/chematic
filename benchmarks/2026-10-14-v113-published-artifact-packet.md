# v1.1.3 published Python artifact packet

This record closes the missing public-artifact validation boundary for
CheMatic 1.1.3. Hosted run
[38474540633](https://github.com/kent-tokyo/chematic/actions/runs/38474540633)
downloads exact public PyPI wheels and runs the version-pinned chemistry gates
on Linux x86-64, macOS arm64, macOS x86-64, and Windows x86-64.

## Result

All four wheel lanes match the checked-in expected manifest. The bounded
headline counts are:

- reactions: 80 exact graph/origin/map rows and three inputs invalid in both
  engines;
- SMARTS: 310,000/310,000 cells with no unexpected result;
- Morgan radius 2, 2,048 bits: 10,000/10,000 exposed rows and 5,000/5,000
  ChEMBL rows bit-identical to RDKit 2026.03.6;
- unusual-input descriptors: HBD 4,072/4,072 and QED 4,071/4,071 comparable
  rows, with one RDKit QED error retained separately;
- A6: 265/265 independently sound and stereo-clean, with no gross clash.

CIP, hybridization, MMFF94 typing, MOL stereo, and the BioTransformer
classification counts also match the expected manifest. Existing MMFF
differences and typed boundaries remain visible; a passing gate does not turn
them into parity claims.

## Wheels

| Lane | Public wheel SHA-256 |
|---|---|
| Linux x86-64 CPython 3.9 | `1ebeaf3ab2c8bc125d984e999f139e3c7cc195d735efaa71cf8a472c17d809a8` |
| macOS arm64 CPython 3.12 | `4c9ad6d18cc933f4fa42f4761c97c080d7ca14161906a05103c617c8f1818416` |
| macOS x86-64 CPython 3.12 | `bb2a32556d3f942686ffccbbbf5fc010b33cc429c7be01e9e4c57b79d020cad9` |
| Windows x86-64 CPython 3.12 | `15664700d9c83065110b14b14174bb014dd572594d54e79d97e01282ae7e06fd` |

The machine-readable packet is
[`validation/results/published-wheel-chemistry-gates-v1.1.3.json`](../validation/results/published-wheel-chemistry-gates-v1.1.3.json).
The complete mismatch, refusal, and unsupported accounting remains in
[`validation/published-wheel-chemistry-gates-expected-v1.1.3.json`](../validation/published-wheel-chemistry-gates-expected-v1.1.3.json)
and the retained workflow artifacts.

## Boundary

The platform jobs execute CheMatic's public wheel. Their serialized outputs
are compared on Linux against the SHA-pinned RDKit 2026.03.6 wheel when the
same RDKit build is unavailable on the originating platform. The packet is
Python-only: npm/WASM and crates.io evidence remain separate. Source-built
timings and post-release source fixes are not promoted to v1.1.3 claims here.
