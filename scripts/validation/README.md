# Validation entry points

`validate_*` scripts check a specific chemistry or artifact contract. Their
inputs, pinned oracle versions, and expected output records are described in
the linked `docs/` or repository-root `validation/` material. Run them from
the repository root. This directory holds executable checks; the top-level
`validation/` directory holds evidence and fixtures.

`check_fingerprint_cross_binding.py` replaces nine older one-operation
fingerprint parity drivers. Select `--profile` and a new `--output` path for
Rust/Python/Node-WASM comparison. It retains each operation's data field and
comparison boundary; existing versioned results are historical records.
The runner checks the Python and generated Node package versions against the
workspace version before producing a new report, and records their artifact
hashes.
