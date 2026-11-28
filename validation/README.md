# Validation evidence

This directory holds the inputs, expected values and recorded results behind
chematic's compatibility and accuracy claims. Executable checks live in
[`scripts/validation/`](../scripts/validation/README.md); the current summary of
what has been verified is the
[validation report](../docs/reference/validation.md).

Existing paths are referenced by CI, benchmark records and documentation, so
add new material beside them rather than moving files.

## Where things are

| Content | Location |
|---|---|
| Per-format and per-contract packets (cases with input and expected values, plus a recorded comparison) | One directory per contract, such as [`bgf/`](bgf/), [`cml_3d/`](cml_3d/), [`v3000_source_roundtrip/`](v3000_source_roundtrip/) |
| Shared molecule corpora and their manifests | [`benchmark_corpora/`](benchmark_corpora/) |
| Dataset, holdout and release-gate manifests | [`manifests/`](manifests/) |
| Frozen independent holdouts | [`trust_holdouts/`](trust_holdouts/) |
| Comparisons with other toolkits | [`rdkit/`](rdkit/), [`openbabel_comparison/`](openbabel_comparison/), [`cdk_indigo_comparison/`](cdk_indigo_comparison/), [`cosmolkit_comparison/`](cosmolkit_comparison/) |
| Upstream issue reproductions | [`rdkit_issues/`](rdkit_issues/), [`openbabel_issues/`](openbabel_issues/), dated `rdkit-*`/`openbabel-*` directories |
| Versioned result records (release channels, gates, parity runs) | [`results/`](results/) |
| Maintainer attestations (for example, unused-data statements for holdouts) | [`attestations/`](attestations/) |

Files directly in this directory, both records and a few runners, are kept at
the paths that other material cites. Put new material in a subdirectory.

## Worked example: BGF round-trip preservation

A per-format packet keeps input, expected values and the recorded result
together.

1. **Input and expected values:** [`bgf/cases.json`](bgf/cases.json) lists 15
   valid and 26 invalid BIOGRF 200 records. Each valid case has an `input`
   string and the `expected` fields an independent fixed-column observer
   should read back.
2. **Command:** build the CLI, then run the observer from the repository root:

   ```bash
   cargo build -p chematic-cli
   python3 scripts/validation/check_bgf_preservation.py \
     --cli target/debug/chematic \
     --obabel "$(command -v obabel)" \
     --output /tmp/bgf-comparison.json
   ```

3. **Recorded result:**
   [`bgf/openbabel-3.2.1-comparison.json`](bgf/openbabel-3.2.1-comparison.json)
   stores per-field preserved counts for chematic and Open Babel, the tool
   versions, and hashes of the packet, source and binaries.
4. **Documentation:** [`docs/reference/bgf-records.md`](../docs/reference/bgf-records.md)
   states the supported scope and links the record.

Write a new run to a new path. Replace a checked-in record only when the
packet, tool versions and hashes it records still describe what you ran.

## CI-produced results

Some gates run in [CI](../.github/workflows/ci.yml) and write to
`validation/results/` inside the runner. For example, the SMIRKS property fuzz
compares against
[`smirks_property_fuzz_expected.json`](smirks_property_fuzz_expected.json):

```bash
python3 scripts/validation/smirks_property_fuzz.py --templates 2000 --seed 754 \
  --expected validation/smirks_property_fuzz_expected.json \
  --output validation/results/smirks-property-fuzz-linux.json \
  --rows validation/results/smirks-property-fuzz-linux-probes.jsonl
```

These outputs are uploaded as workflow artifacts and are not committed. Only
records that a document cites as evidence are checked in.

## Adding a validation

1. Put the cases (input with expected values) in a new subdirectory named for
   the contract.
2. Put the runner in `scripts/validation/` and give it an explicit `--output`.
3. Record the tool versions, corpus, options and hashes in the result, and keep
   refusals and failures in its counts.
4. Link the result from the relevant page under `docs/reference/`.
