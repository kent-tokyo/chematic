# MOPAC Cartesian keyword continuations

[Open Babel #1231](https://github.com/openbabel/openbabel/issues/1231) requests
reading more than one MOPAC keyword line. Native `parse_mopac_cartesian` and
`write_mopac_cartesian` now support the standalone `+` continuation profile.
The existing Python functions delegate to these same Rust entry points.

```rust
use chematic_mol::{parse_mopac_cartesian, write_mopac_cartesian};

let input = "PM7 XYZ +\nCHARGE=+1\nexample\ncomment\nC 0 0 0\n";
let geometry = parse_mopac_cartesian(input).unwrap();
assert_eq!(geometry.keywords, "PM7 XYZ +\nCHARGE=+1");
assert_eq!(geometry.title, "example");
assert_eq!(geometry.comment, "comment");
assert_eq!(geometry.source_text, input);
let output = write_mopac_cartesian(&geometry).unwrap();
assert_eq!(parse_mopac_cartesian(&output).unwrap().keywords, geometry.keywords);
```

## Physical header lines

A standalone `+` outside double-quoted values announces another keyword line.
It can appear at the beginning, middle or end of the current line. Every
nonfinal keyword line must contain it; the final line must not. The following
two physical lines are always the title and comment, including when blank.
This follows the official [OpenMOPAC + specification](https://openmopac.net/Manual/plus.html).

`CHARGE=+1` is an ordinary token and does not continue a line. A `+`, `++`, `&`
or `A0` inside a double-quoted value stays opaque text. Double quotes must
balance on each physical keyword line; escape processing and multiline quoted
values are outside this profile. Other calculation keywords remain opaque.
No Hamiltonian, charge, spin, calculation validity or bond topology is inferred.

`MopacCartesianGeometry.keywords` now represents an LF-separated physical
keyword block rather than always a single line. It retains the continuation
markers, spacing, tabs and empty keyword lines. CRLF input normalizes this
field to LF while `source_text` retains the exact original input. The writer
accepts LF blocks, validates their continuation structure and emits their
current contents. It retains the existing numeric serialization and geometry
semantics; it does not flatten keywords or return stale source text.

## Limits and boundaries

- At most 64 physical keyword lines, including an empty final keyword line.
- At most 65,536 bytes in the LF-separated keyword block, counting separators.
- At most 16,384 bytes per keyword line; caller parsing limits may tighten this.
- The existing input/atom/description/serialized-output limits still apply.
- Controls other than tab are refused in header fields. Invalid source
  diagnostics use physical line numbers, including after continuations.

The writer rejects extra physical keyword lines without an announcing `+`,
and a final line that still announces a continuation. Atom flags, labels,
translation vectors and original-source behavior retain the existing
[Cartesian geometry contract](openbabel-mopac-cartesian.md).

`++` and `&` have different semantics in the
[OpenMOPAC ++](https://openmopac.net/Manual/plusplus.html) and
[OpenMOPAC &](https://openmopac.net/Manual/ampersand.html) specifications.
They remain unsupported outside quoted values. Atomic-unit `A0` input,
internal/hybrid coordinates, global/interleaved comment syntax, multiple jobs
and running MOPAC calculations are not added. The separate archive reader's
restart-deck selection is unchanged.

## Verification

Original native fixtures cover continuation placement, 2–64 lines, signed
keyword values, quoted marker text, blank descriptions, CRLF/source retention,
coordinate flags/labels, writer validation, byte limits and physical diagnostics.
Native library and existing Cartesian/archive regression tests are run directly.
The new Python regression cases are supplied for a freshly rebuilt extension;
they are not claimed as locally executed against an older binary.

Open Babel 3.2.1 reads the single-line controls but returns no molecule for
the authored continuation decks. The reference run therefore separates direct
continuation diagnostics from explicit normalization controls: after collapsing
these unquoted synthetic keyword headers to one line, both original inputs
and native writer outputs are cross-read for elements/coordinates. This does
not claim that Open Babel directly supports the new multiline output, retains
all metadata, or runs the specified calculation. Reproducible records are in
`validation/results/mopac-keyword-continuation-2026-10-18/`.
