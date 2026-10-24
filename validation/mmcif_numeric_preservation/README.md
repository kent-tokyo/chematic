# mmCIF numeric preservation packet

Original MIT fixtures authored on 2026-10-24 from wwPDB `atom_site`/`cell`
definitions: three controls, twelve reproductions and four arithmetic storage
cases. The authored oracle is finite binary64 numeric values and ordered atom
IDs/elements, not original numeric spelling or uncertainty suffixes.

`baseline.json` and `candidate.json` retain every Python/WASM/Open Babel output,
error, resource outcome, source/runtime identity and independent observation.
Baseline runtime commit: `342d1d5f9ad09aa6ade4e3079a902cba856209e2`; relevant
sources match main base `fcbf6d0a`. Candidate source hashes identify the changed
source snapshot separately from its pre-commit HEAD.

See [measurement scope and reproduction](../../docs/reference/mmcif-numeric-preservation.md).
The packet does not establish physical/deposition validity, complete CIF/schema
or namespace preservation, connectivity, source-byte retention or throughput.
