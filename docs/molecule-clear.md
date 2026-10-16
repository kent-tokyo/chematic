# Reset an editable molecule completely

[Open Babel #3032](https://github.com/openbabel/openbabel/issues/3032) describes
a reused molecule retaining its previous title after `Clear()`. chematic's
baseline had no complete reset entry point. It now provides Rust
`Molecule::clear()` and Python `RWMol.Clear()` with a fresh-object contract.
This is an additive API, rather than a correction to an existing chematic
`Clear` implementation.

```python
import chematic

source = chematic.from_smiles('CCO')
source.SetName('old title')
editor = chematic.RWMol(source)
snapshot = editor.GetMol()
editor.Clear()
assert editor.GetName() == ''
assert editor.GetNumAtoms() == 0
assert snapshot.GetName() == 'old title'
assert editor.AddAtom(7) == 0
```

## Contract

Rust `clear(&mut self)` replaces the graph with a freshly built empty molecule.
It releases atoms, bonds and adjacency, atom tags, R labels, stereo groups and
bond members, stereo-neighbor order, stashed bond directions/anchors,
SMILES ring-order metadata and every derived cache. Previously cloned molecules
and retained derived views remain valid. It releases storage and invalidates
old indices; subsequent atom insertion starts at index zero.

Python `RWMol.Clear()` also resets every typed property and the conventional
`_Name` title. Source `Mol` objects and snapshots returned by `GetMol()` remain
independent. Calling Clear repeatedly is supported. Empty graph observables
match a newly constructed editor, including existing empty-formula behavior.
No default parser/writer behavior changes.

The `chematic.rdkit_compat.RWMol` wrapper exposes the same operation as a
chematic extension. The native stub includes the complete current editable
subset and the zero-argument `Clear() -> None` signature. The core graph has no
wrapper-owned title/property store; those are reset by the Python editor.
There are no keep flags, conversion-stream resets or WASM entry points here.

## Bounded validation — 2026-10-16

Base: `faa8092e6050ad6250c68498970245cfbc6579e9`.
The original four-cycle native probe reproduced Open Babel 3.2.1 retaining
`original title 0` through `original title 3` after a successful clear and an
untitled read. Fresh Open Babel molecules had empty titles, with the same
formula as the reused molecules.

The retained pre-change Python module failed all ten reset controls with
`AttributeError` because both native and compatibility editors lacked `Clear`.
The candidate passed 162 core Rust unit tests (three new lifecycle controls)
and 173 Python tests (11 new and 162 existing compatibility tests).

Controls cover all native metadata stores, three scalar cache slots, the ring
count cache, a retained typed derived view, source/snapshot independence,
string/integer/double/boolean properties, tetrahedral groups, isotopes/maps,
200 Python reuse cycles and index-zero rebuilding. Rebuilt methanol MOL output
matches a freshly parsed methanol graph byte-for-byte. These are lifecycle and
serialization controls, not broad chemical-equivalence or speed evidence.

The native probe is
`validation/openbabel_issues/molecule-clear/openbabel_probe.cpp`; source/module
hashes, raw outputs, commands and compressed logs are retained in
`validation/results/molecule-clear-2026-10-16/`. Strict Clippy, formatting,
refactor budget, static binding inventory, shared manifest, unsafe/naming
checks and strict MkDocs passed. This is local macOS source validation, not a
published wheel or a completed Linux/Windows/WASM validation claim.
