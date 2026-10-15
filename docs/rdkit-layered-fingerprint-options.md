# Configurable RDKit Layered fingerprints

This candidate fills the gap between chematic's existing fixed-default Layered
fingerprint and the parameter controls exposed by
[RDKit `LayeredFingerprint`](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html#rdkit.Chem.rdmolops.LayeredFingerprint).
The reference is **RDKit 2026.09.1**, tag `Release_2026_09_1`, commit
`fece8caa860bdf6c9c82bdb9253f22a0cf87150c`. RDKit classifies this family as
experimental, so parameter/result compatibility is version-specific.

The existing `rdkit_layered_fp()` / `Mol.rdkit_layered_fp()` profile remains
2048 bits with all six active layers and paths of 1–7 bonds. It and the new
API share the same per-path hashing code. chematic's separate native
`layered_fp()` / `layered_fp_layers()` family uses a different seven-layer
scheme and is not interchangeable with these RDKit bit positions.

## API

Rust exports `RdkitLayeredConfig`, `RdkitLayeredResult`,
`RdkitLayeredError` and `rdkit_layered_config(mol, &config)` from
`chematic_fp`. The result includes packed LSB-first `bits`, `fp_size`,
`atom_counts`, `selected_paths` and `search_visits`.

Python exposes `Mol.rdkit_layered_config()` with the corresponding keyword
arguments and a detached result dictionary:

```python
import chematic

mol = chematic.from_smiles("CC(C)O")
initial_counts = [0] * mol.atom_count
result = mol.rdkit_layered_config(
    layer_flags=0x25,
    min_path=1, max_path=5, nbits=1024,
    from_atoms=[0], atom_counts=initial_counts,
)
packed = result["bits"]                 # bytes, LSB-first
counts = result["atom_counts"]           # a new list
assert initial_counts == [0] * mol.atom_count

linear = mol.rdkit_layered_config(
    branched_paths=False, from_atoms=[1],
)
```

Python result keys are `bits`, `nbits`, `atom_counts`, `selected_paths`
and `search_visits`. The bytes use `ceil(nbits/8)` storage with zero unused
high bits in the last byte. A packed mask uses the same convention and must
have exactly that length and valid zero padding. A mask is provided as
`set_only_bits=bytes(...)`; integer lists are not silently interpreted as masks.

| Option | Default | Checked contract |
| --- | --- | --- |
| `layer_flags` | `0xffffffff` | u32; six active low bits, unknown high bits ignored as in the reference |
| `min_path`, `max_path` | 1, 7 | 1 ≤ min ≤ max ≤ 12 bonds |
| `nbits` / Rust `fp_size` | 2048 | 1–65,536, including odd/non-power-of-two sizes |
| `branched_paths` | True | False requires an explicit nonempty root list |
| `from_atoms` | None / empty | stored zero-based atom indices; empty means the ordinary all-subgraphs calculation |
| `set_only_bits` | None | packed bytes selecting which output bits may be set |
| `atom_counts` | None | optional initial u32 counters, at least atom_count entries and at most 65,536 |
| `max_visits` | 1,000,000 | positive, maximum 10,000,000 enumeration visits |
| `max_paths` | 100,000 | positive, maximum 1,000,000 selected path occurrences; also caps the linear frontier |

Flags select topology (`0x01`), bond order (`0x02`), atom types (`0x04`),
ring membership (`0x08`), minimum ring size (`0x10`) and aromaticity (`0x20`).
Counts increase once per accepted path containing an atom, not once per layer,
new output bit or distinct bit. A collision with an already-set bit still
contributes a path count. If the mask rejects every selected layer for a path,
that path does not increase any counter. Flags with no active layers likewise
leave the initial counts unchanged.

Initial counter lists are copied instead of mutated in place. Values add with
u32 wrapping, matching the reference's unsigned arithmetic. Extra trailing
entries remain unchanged. With no initial list, the returned list contains
one zero-initialized counter per stored atom. An explicitly supplied short or
empty list is rejected rather than interpreted as an absent output request.
Python integer arguments/list entries exclude bool, float and text; the
branched flag must be a bool. Bounded integer iterators are accepted for
roots/counters and rejected after their size cap.

## Root and path semantics

Branched subgraphs use only bonds whose endpoints are not explicit H atoms.
A rooted branched subgraph is any connected subgraph containing a requested
root. Duplicate roots retain their multiplicity: bits are unchanged while
counters can increase more than once for the same subgraph. The checked API
retains the caller's stored atom indices, including H slots. Explicit H atoms
have no eligible bonds and receive no added counts. Comparing with RDKit
therefore requires the corresponding indexed graph, such as parsing with
`SmilesParserParams.removeHs=False`.

Rooted linear paths start at the requested atom. They retain RDKit's specific
range-wide closure convention: revisiting an earlier atom is allowed only in
the final extension at the requested maximum length, except immediate
backtracking. Deduplication uses the set of bond indices within each length
and root. This can make a shorter closed-ring contribution disappear when
`max_path` increases; it is preserved reference behavior, not replaced with
an assumed monotonic simple-path definition.

### Unrooted linear refusal

The pinned RDKit C++ unrooted `branchedPaths=False` branch passes atom paths
to code that treats their indices as bond indices. In isolated subprocess
probes during this implementation, `CC` and `CCO` both terminated with signal
11; some ring graphs instead returned numerical fingerprints. The checked
chematic API rejects this entire combination with
`unrooted_linear_unsupported`, avoiding a molecule-dependent unsafe branch.
It does not silently substitute another definition. Supply `from_atoms` when
requesting linear paths. An all-atoms root list is a rooted calculation with
its own multiplicities and is not an alias for the unrooted branch.

## Checked graph and resource boundaries

The current profile accepts at most 512 atoms and 1024 bonds. Complex query
bonds, zero/dative bonds, loops and parallel edges are refused. Ordinary SMILES
dummy atoms use atomic number zero. The existing RDKit aromaticity view is
reused, with atom/bond correspondence verified before applying caller indices.
If sanitization changes that correspondence or introduces an unsupported bond
type, calculation returns a typed error.

Branched enumeration is bounded connected-subgraph search. Root filtering
uses the global unique subgraphs and applies root multiplicities; roots do not
promise reduced search work. Rooted linear enumeration bounds its intermediate
frontier as well as accepted paths. Budgets cover enumeration, with graph caps
also bounding preparation; they are not wall-time guarantees. Exhaustion
returns no partial fingerprint or partially updated caller counters.
`selected_paths` and `search_visits` are chematic diagnostics, not RDKit API
return values or cross-engine performance measurements.

Rust errors carry `code` and `detail`. Python calculation errors subclass `ValueError`
through `ChematicInputError` and expose `category`, `code` and
`format="rdkit_layered_config"`. Category distinguishes malformed options,
unsupported branches and resource limits. Invalid Python argument types are caught before core calculation and may
raise TypeError or plain ValueError without those attributes. No unsupported
option is rounded or coerced to another profile.

Ring hashing retains the existing SSSR-based implementation. The measured
ring/bridged/spiro examples below are reference checks, not proof of arbitrary
ring-perception or chemistry parity. This is not a SMARTS-query fingerprint
API, an unrestricted RDKit wrapper, or a complete RDKit feature implementation.

## Source and independent validation

Behavioral details derive from the pinned RDKit
[Fingerprint source](https://github.com/rdkit/rdkit/blob/fece8caa860bdf6c9c82bdb9253f22a0cf87150c/Code/GraphMol/Fingerprints/Fingerprints.cpp)
and
[Subgraph source](https://github.com/rdkit/rdkit/blob/fece8caa860bdf6c9c82bdb9253f22a0cf87150c/Code/GraphMol/Subgraphs/Subgraphs.cpp).
The BSD-3-Clause copyright/license inventory is updated in
`THIRD_PARTY_NOTICES.md`; source comments retain the applicable notice.
No RDKit binary is bundled or called by the implementation.

The original minimized corpus is described in
`validation/rdkit_issues/fingerprints/layered-config.json`, including license,
provenance, normalized payload hash, exact oracle and settings. It includes
empty/singleton graphs, branches, rings, fused/bridged/spiro/cage examples,
aromatic/Kekulé inputs, charges, explicit/isotopic H and dummy atoms. The empty
native graph is constructed with `RWMol().GetMol()` because empty SMILES input
is rejected by chematic's parser.

`python scripts/validate_layered_config.py` compares all output bits and
accumulated counters exactly. It covers every active-layer subset, all 78
supported path ranges, odd and maximum bit lengths, masks, counter seeds and
wrapping, repeated/endpoint roots, range-dependent closures and atom
permutations. The oracle's bits are packed independently. Separate controls
compare the existing fixed default with the new default and check input graph,
properties and seed-list nonmutation. Exception/diagnostic statistics are not
represented as bit/counter parity.

The run passed **13,378** cases across 20 original fixtures, with exact bits
and **75,874** accumulated counter values. This includes 2,790 masked cases,
2,850 rooted-linear cases and 360 permuted-atom cases, covering all 64 active
layer subsets and all 78 supported min/max path ranges. There were no bit,
counter or input-mutation mismatches. **110** controls preserved the existing
fixed-profile output.

Rust validation passed **346 existing library + 12 new integration tests**;
Python passed **64** tests. Clippy with warnings denied, formatting, strict
MkDocs, stub syntax, static binding inventory and refactoring budget checks
passed. These are targeted local checks, not a full-workspace hosted CI run.

Evidence is saved to
`validation/results/rdkit-layered-config-2026-10-15.json`, with source,
fixture and compiled-extension hashes. It is a local source validation,
not a release or hosted CI result.

```sh
cargo test -p chematic-fp --lib --test rdkit_layered_config --offline --quiet
cargo clippy -p chematic-fp -p chematic-py --all-targets --offline --no-deps -- -D warnings
cargo fmt --all -- --check
python -m pytest crates/chematic-py/tests/test_layered_config.py -q
python scripts/validate_layered_config.py
python scripts/check_refactor_budget.py
python scripts/check_binding_surface.py
python scripts/check_product_names.py
mkdocs build --strict
```

The source candidate preserves the previous work through indexed SDF selections.
Its review patch is `validation/rdkit_issues/layered-config-2026-10-15.patch`;
applying it requires that inherited context. Main integration, feature commit,
push and publication remain separate steps. Shared hashing removes one resolved
refactoring function allowance; no new allowance was added.
