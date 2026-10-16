# Inline model-gate evidence packet

Read [the report](../2026-10-16-python-model-hook-inline.md) with the predeclared measurement plan and `decision.json`. Complete single-epoch matrices are compressed without dropping any raw lane. `observer-*` files contain fresh post-timing values; their elapsed fields were discarded. Source output dumps and wheel hashes are in `parity.json`; full binaries and 14,072-input streams remain local at `/private/tmp/chematic-model-hook-inline-evidence-20261016`.

`candidate.patch` reconstructs the two changed/added Rust files from base `3e7f04b20593af4b07569d886b56afc1184468fe`. Compiler flags and build logs, the complete source/input manifest, package-local native library hashes, frozen-function extraction proof and validation commands are included. The cleanup record lists only unused regenerable compiler metadata/static libraries removed after process/open-file checks. No timing worker was restarted.
