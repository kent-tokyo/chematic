# Reproducible diagnostics

These scripts reproduce a named residual, probe a published contract, or
classify an existing evidence corpus. They are kept because an associated
result, CI check, test, or source example still uses them. Run them from the
repository root with the inputs and dependency versions named in each script
or its linked record. Outputs belong under `validation/` or a temporary
directory, never beside the script.

The `stereo2d_*` scripts form one fixture family; the `aromaticity_*`,
`canonical_*`, `ez_*`, and `mmff94_*` scripts cover separate mechanisms. A one-time
experiment with no retained result or current reference should be removed
after its finding is recorded.

`v0_14_0_ez_bound_release_gate_report.py` replays the checked-in v0.13.0 and
v0.14.0 pipeline rows. Its historical result reports failed release gates;
moving the script does not reinterpret that result.
