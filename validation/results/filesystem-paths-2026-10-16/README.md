# Python filesystem-path validation

This is source validation for the filesystem subset of RDKit #9601.
See `local-validation.json` for exact commits, source/artifact hashes, command
lines, build profile and limitations. Compressed logs retain the actual runs.

The old retained module passed 30 string acceptance cases and failed 75 new
path cases. The four affected I/O/report source files match the current pinned
base. The candidate passed all 105 acceptance cases and 390 selected checks
in total, with 30 skipped raw-filename file round trips. Python's own file
creation returns EILSEQ on this macOS filesystem; those controls separately
verify OSError at the chematic boundary and remain executable on supporting
POSIX filesystems. This is not Linux/Windows runtime evidence.

Rebuild the candidate using the normal native Python tooling, set PYTHONPATH
to that package, and run the pytest command recorded in the JSON. The original
fixtures are created in pytest's temporary directory by
`crates/chematic-py/tests/test_filesystem_paths.py`; no upstream dataset or code
is copied. Paths containing Japanese characters and spaces are exercised.

An AST test checks all 15 path stub annotations. This is not a complete static
type checker or a wheel-content test. No chemical algorithms or native formats
change, and this record makes no speed, release or general RDKit parity claim.
