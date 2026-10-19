# Reproducible diagnostics

These scripts reproduce a named residual, probe a published contract, or
classify an existing evidence corpus. They are kept because an associated
result, CI check, test, or source example still uses them. Run them from the
repository root with the inputs and dependency versions named in each script
or its linked record. Outputs belong under `validation/` or a temporary
directory, never beside the script.

The `stereo2d_*` scripts form one fixture family; the `aromaticity_*`,
`canonical_*`, and `ez_*` scripts cover separate mechanisms. A one-time
experiment with no retained result or current reference should be removed
after its finding is recorded.
