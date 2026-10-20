# Explicit PDB CONECT graph

RDKit [#9581](https://github.com/rdkit/rdkit/issues/9581) reports a ligand
whose sulfur is rejected with valence seven when proximity bonding is enabled.
The original input reproduces this with RDKit 2026.09.1. Disabling proximity
bonding already works in RDKit: the explicit sulfur bond-order sum is six.

Native `from_pdb` / `from_pdb_strict` previously offered only geometric bond
inference. On the original ADX input, the unchanged native `from_pdb` yields
44 bonds and five single sulfur edges; the explicit graph has 43 bonds with
four sulfur neighbours and order sum six. The new API reads an explicit graph without adding distance-based
edges. This is a native alternative with a deliberately narrow input contract,
not a replacement for RDKit's sanitizer or a general PDB chemistry reader.

## API

```python
frame = chematic.from_pdb_conect(pdb_text)
mol = frame["mol"]
coords = frame["coords"]       # source XYZ in Å, including explicit H
serials = frame["serials"]     # source atom serials in molecule/coordinate order

# Only when every extra CONECT slot really denotes a covalent target:
frame = chematic.from_pdb_conect(pdb_text, allow_extended_targets=True)
```

Rust exports `parse_pdb_conect`, `parse_pdb_conect_with_options`, and
`parse_pdb_conect_with_options_and_limits` from `chematic_mol`. The returned
`PdbConectFrame` contains `molecule`, `coords` and `serials`.
`PdbConectOptions` defaults to refusing extended target slots.

WASM exports `mol_from_pdb_conect(text, allowExtendedTargets)` and
`pdb_conect_frame_json(text, allowExtendedTargets)`. The JSON frame contains
`atoms` (`atomic_number`, `charge`), zero-based `bonds` (`[a,b,order]`),
`coords` and `serials`. Pass `false` for the conservative profile.

## Accepted profile

- Unmodelled ATOM/HETATM records with blank alternate-location columns.
  MODEL/ENDMDL input is refused, including a single explicit model. Select
  a model/location upstream before using this API.
- Positive unique decimal atom serials, explicit recognized element fields,
  finite fixed-column XYZ, and a valid residue-sequence field. Atom names
  never determine the element. Canonical element symbols and uppercase
  two-letter symbols are accepted; D/T isotope aliases are refused.
- Formal charge comes only from blank or `1+` through `9+` / `1-` through `9-`
  in columns 79–80. No residue charge inference is performed.
- CONECT targets must refer to listed atoms. Self-links, incomplete nonblank
  target slots, unknown references and conflicting reciprocal multiplicities
  are refused. Nonsequential source IDs and continuation rows are supported.
- This **multiplicity profile** interprets repeated targets in one direction,
  including across continuation rows, as single/double/triple bond order.
  If both directions are given, they must agree and are counted once.
  Multiplicity above three is refused. This convention extends the connectivity
  fields; the wwPDB specification itself does not establish these bond orders.
- Missing CONECT edges stay absent. With no CONECT records the graph has no
  edges, even if atoms overlap. Every listed H atom remains in the graph.

The [wwPDB 3.3 CONECT format](https://www.wwpdb.org/documentation/file-format-content/format33/sect10.html)
provides four targets per physical row, continuation rows, and redundant
opposite declarations. The default rejects nonblank target slots after the
fourth. Legacy contact slots and the six-target row in the original #9581
input cannot be safely interpreted as covalent targets automatically;
`allow_extended_targets=True` makes that interpretation explicit.

## Limits and losses

The Rust limit-aware API combines `PdbParseLimits` (input bytes, line bytes,
atoms and models) with `max_bonds` and `max_targets` for total directed target
occurrences. Defaults are 256 MiB input, 1 MiB lines, 10,000 atoms, 30,000
unique bonds and 180,000 targets. Python uses these defaults; WASM also applies
its existing tighter input/atom caps. Errors are typed `PdbConectError` in Rust
and exceptions in the bindings. There is no allocation from an unchecked
serial number or target count.

This is a graph projection: residue/name/chain fields, occupancy, B-factors,
ANISOU, LINK, SSBOND and other source metadata are not returned or interpreted.
Use source-preserving PDB selection separately when those fields matter.
No valence sanitization, bond-order inference, aromaticity, isotope or stereo
perception is performed. Implicit H follows the existing native molecule model
and is outside the explicit-graph comparison. Existing PDB APIs are unchanged.

## Validation

The authored fixture is `validation/rdkit_issues/pdb-conect/sulfur.pdb`: four
sulfur neighbours with order sum six, charged oxygen atoms, nonsequential IDs
and a nearby isolated explicit H. Rust, Python and real Node/WASM contract
tests cover graph/coordinate order, reciprocal/continuation declarations,
extended-slot opt-in, malformed inputs and resource fences.

`scripts/validation/validate_pdb_conect.py` compares indexed atoms, formal charges,
single/double/triple edges, source coordinates and serials against pinned
RDKit 2026.09.1 with `sanitize=False, removeHs=False, proximityBonding=False`.
It checks Python, WASM frame JSON, and the WASM molecule handle through a
V3000 graph writer. Permutations, coordinate translation/scaling, record ordering, CRLF,
one-way/reciprocal declarations and continuation multiplicities are varied.
The optional unedited #9581 PDB block supplies a separate defect reproduction.
Raw inputs/outputs and source/runtime hashes are retained with the result.

The local [2026-10-16 source packet](https://github.com/kent-tokyo/chematic/tree/codex/rdkit-pdb-conect-20261016/validation/results/pdb-conect-2026-10-16)
records **102/102** matches in each of the Python, WASM frame and WASM handle
graph lanes (101 authored cases plus the original ADX issue input). The
original input is accepted only with extended targets explicitly enabled;
its default RDKit sanitization failure is recorded separately. This is local
source validation, not a published package or completed hosted-CI result.

This evidence concerns explicit graph preservation only. It does not establish
general PDB compatibility, implicit-H equivalence, stereo correctness, speed
or superiority to RDKit.
