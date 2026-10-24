# Validation entry points

`validate_*` scripts check a specific chemistry or artifact contract. Their
inputs, pinned oracle versions, and expected output records are described in
the linked `docs/` or repository-root `validation/` material. Run them from
the repository root. This directory holds executable checks; the top-level
`validation/` directory holds evidence and fixtures.

The `check_same_process_*_contract.py` runners validate individual format
contracts. `check_same_process_contracts.py` checks their retained eight-format
bundle and is called by `scripts/check.sh`. Their historical result JSON files
remain in `benchmarks/`; run the checks from the repository root with the
version and inputs specified by the dated benchmark records.

`check_fingerprint_cross_binding.py` replaces nine older one-operation
fingerprint parity drivers. Select `--profile` and a new `--output` path for
Rust/Python/Node-WASM comparison. It retains each operation's data field and
comparison boundary; existing versioned results are historical records.
The runner checks the Python and generated Node package versions against the
workspace version before producing a new report, and records their artifact
hashes.

`score_trust_holdouts.py` retains the one-shot scorer for frozen trust
holdouts. Its published row-level score is historical evidence; use the exact
candidate artifact and manifest described in the dated benchmark record.

Manual parity and runtime probes also live here: `check_conversion_report_parity.py`,
`check_descriptor_rdkit_holdout.py`, `check_mcp_runtime_inventory.py`,
`check_native_descriptor_regression.py`, and `check_streaming_gzip_openbabel.py`.
They have no current CI caller; supply the source-built package, external tools,
and input paths required by each command's `--help` before using its output as
new evidence.

`compare_rdkitjs_reaction_runs.py` compares two pinned RDKit.js reaction runs.
Its existing differential report remains under `validation/results/`.
