# Complete molecule reset — 2026-10-16

The original native C++ probe reproduces Open Babel 3.2.1 keeping the title
across four Clear/untitled-read cycles. `openbabel-probe.json` records the actual
output. Build it against the installed Open Babel using the compiler/linker
flags from `pkg-config --cflags --libs openbabel-3`, then run the executable.
No upstream code or fixture was copied.

The baseline retained Python module is tied to main commit a89446d3. Its native
RWMol and compatibility wrapper sources match the current base. Ten new reset
controls failed because Clear did not exist. This is API-absence evidence,
not an existing chematic title-leak reproduction.

The candidate passes 162 core unit tests and 173 Python tests, including three
new native lifecycle checks, eleven new Python checks and two 100-cycle reuse
loops. Existing snapshots/source properties remain intact. Rebuilt methanol
MOL output equals fresh output. Native metadata fixtures test storage/reset
behavior without claiming chemically valid enhanced-stereo assignments.

`local-validation.json` contains commands, exact source/artifact hashes,
versions, build profile and scope. Compressed logs retain actual runs. Rebuild
the native module using normal Python tooling, set PYTHONPATH to that package,
and execute the recorded pytest command. The native Rust tests use
`cargo test -p chematic-core --lib --offline`.

These are local source/lifecycle controls, not publication, speed, complete
Open Babel state compatibility or cross-platform proof. Remaining hosted
checks are tracked on the PR. The compatibility wrapper's Clear is a chematic
extension.
