# Native V3000 3D numeric regeneration

`cases.json` contains original MIT synthetic inputs: nine paired records and
three native-only finite extremes. Expected coordinate tokens are authored;
no engine supplies the expected numeric values. Two paired controls have
mirrored tetrahedral geometry; other cases are numeric stress/control inputs,
not physically qualified molecular geometries. Repeated writers/round trips
are not extra corpus samples or an independent holdout.

Build and run from the repository root with RDKit available to Python:

```sh
cargo build --locked -p chematic-mol --example v3000_3d_numeric_dump
python scripts/validation/check_v3000_3d_numeric.py \
  --native-executable target/debug/examples/v3000_3d_numeric_dump \
  --output /tmp/v3000-3d-candidate.json.gz
cargo test --locked -p chematic-mol --lib --test v3000_3d_numeric
python -m unittest scripts.tests.test_v3000_3d_numeric
```

The source and each output are independently read by RDKit 2026.09.1 with
`sanitize=False`, `removeHs=False`, `strictParsing=True`. CIP assignment uses
`AssignStereochemistry(cleanIt=True, force=True)`. Source qualification first
checks the authored coordinate bits. Output coordinates, indexed atomic/bond
fields and atom-indexed source CIP codes are separate observations. The raw
numeric observer handles continuation lines without either engine's reader.
Six observer controls detect coordinate rounding, signed-zero loss, element/
row changes, changed stereo and unreadable/incomplete output.

The native plain/checked paths share the existing 3D conformer writer. The
checked wrapper continues to validate square-planar consistency rather than
arbitrary edited metadata or all invalid numbers. Open Babel uses precisely
`obabel -imol -omol -x3`; output is V3000, not default V2000. Each subprocess
has a 10-second and 1-MiB combined-output bound.

Saved baseline: debug executable built from main
`aeae3f53deb34579be60c054cdf22cd686d9254d` with only the diagnostic example
added, copied before production edits. The candidate is a working-tree build
on that parent. `--baseline --source-revision <SHA> --probe-source <saved.rs>`
records the pinned baseline source and suppresses candidate-only acceptance.
A revision flag is operator provenance, not a proof of arbitrary executable
origin. Reports save executable/probe/production-source/packet/harness/Open
Babel/RDKit-parser hashes, versions, commands, raw output/errors, source and
per-case observations. Gzip timestamps are zero.

Reports: `../results/v3000-3d-numeric-{baseline,candidate}-2026-10-24.json.gz`.
Candidate exact coordinates: 9/9 in both native outputs. Older native and Open
Babel 3.2.1: 2/9 each. All three paths retain the nine indexed graph profiles
and source CIP arrays. Native-only min-subnormal, min-normal and max-finite
cases are excluded from external scores. No speed/release/general superiority
claim follows from this packet. 2D and source-retention paths are separate.
