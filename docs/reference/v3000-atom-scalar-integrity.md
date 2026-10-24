# V3000 atom scalar integrity

V3000 atom parsing now refuses malformed, duplicate or unrepresentable `CHG`,
`MASS` and positional atom mapping values. Previously an invalid value could
silently become charge zero, an unspecified isotope or an absent mapping.
Attribute-looking text inside quoted strings or parenthesized values no longer
changes the atom's charge or isotope.

## Contract

| Field | Accepted storage values | When absent |
|---|---|---|
| `CHG` | Integer in -128..127 | Charge 0 |
| `MASS` | Integer in 1..65535 | Unspecified isotope; `D` defaults to mass 2 |
| Positional atom map | Integer in 0..65535 | Required positional field; 0 means unspecified |

These are the core model's storage limits, not a declaration of chemical validity
or full CTfile conformance. `MASS=0`, a bare `CHG`/`MASS`, fractions, exponents,
overflow and repeated attributes (including identical repeats) return
`InvalidAtomLine` with the field name and physical source line. An invalid
`MASS` on `D` also errors instead of taking the deuterium fallback.
No public function signatures change.

The [CTfile specification](https://mdli.com/downloads/literature/ctfile.pdf)
describes integer charge, positive isotope mass and positional atom mapping.
The scope-aware scanner considers `CHG` and `MASS` only at attribute starts
outside quoted strings and parenthesized values. Doubled quotes are handled;
retained opaque text survives editable V3000 regeneration. This is not a complete
opaque-attribute grammar validator. `HCOUNT`, R-group and other attribute
interpretation retain their existing behavior. Query, valence, bond, hydrogen,
stereo and coordinate validation are outside this change.

## Fixed packet, 2026-10-24

The baseline is commit `342d1d5f9ad09aa6ade4e3079a902cba856209e2`, built from a
clean detached worktree. Baseline and candidate use the same authored packet,
observer, RDKit 2025.09.3 reader and Open Babel 3.2.1 executable. Python and Node
WASM runtimes were rebuilt from their respective source trees. Evidence records
hash the production sources, fixture packet, harness, extensions, WASM artifacts,
RDKit parser and Open Babel binary. Paths describe the local measurement run and
are not promises that temporary build artifacts remain available.

| Positive scope | Baseline Python / WASM | Candidate Python / WASM | Open Babel default / atom classes |
|---|---|---|---|
| 8 ordinary scalar cases | 8/8 / 8/8 | 8/8 / 8/8 | 6/8 / 6/8 |
| 2 local quoted `LABEL` compatibility controls | 0/2 / 0/2 | 2/2 / 2/2 | 0/2 / 0/2 |

Ordinary cases cover charge, isotope, the `D` alias and atom mapping. The two
Open Babel ordinary mismatches lose atom mapping. `LABEL` controls are a local
compatibility extension; they are kept separate from ordinary scalar scores.
The new positive improvement is the quoted-text interpretation fix. Ordinary
scalar preservation already passed on the baseline.

RDKit observes **emitted SMILES**, with sanitization disabled and hydrogens not
removed. Each result must contain exactly one unbonded atom and match four
authored fields: atomic number, formal charge, isotope and atom mapping. Both
direct and common conversion routes must agree with that oracle. Local MolJSON
additionally checks element, charge and isotope; its schema has no atom-map
field. No hydrogen-count, multi-atom graph or stereochemical score is inferred.

Open Babel runs `obabel -imol -osmi -xn`, and separately the same command with
`-xa`. These are its documented
[SMILES output options](https://openbabel.org/docs/FileFormats/SMILES_format.html):
omit names and include atom classes when available. Both profiles are retained
even though `-xa` does not recover the maps in this installed build.

## Refusal and source retention are separate checks

For 25 malformed, duplicate or storage-overflow cases, baseline Python and WASM
produce field-specific errors on 0/25 cases; the candidate does so on 25/25,
across direct, common-to-SMILES and common-to-V3000 routes. This is a refusal
diagnostic, not a positive preservation score. Open Babel's `-xa` profile
converts 23 and emits no output for 2; it records no process errors, signals or
resource-limit events. Empty output is not credited as safe refusal.

Source-preserving Rust, Python and WASM routes retain the 10 accepted source
blocks exactly, on both baseline and candidate. Baseline byte retention did not
prevent the incorrect quoted-field interpretation. Byte equality alone therefore
does not prove the parsed chemistry is correct. CLI default conversion regenerates
V3000 rather than retaining its exact source bytes; the native and CLI tests
check the regenerated scalar fields and quoted text separately.

Native tests also cover continuation-line diagnostics, storage endpoints,
parenthesized text and doubled quotes. Python, WASM and CLI tests cover accepted
and refused inputs. CLI errors emit no partial stdout and leave an existing
output file untouched. Observer negative controls test field changes, invalid
types and distinct process outcomes without requiring RDKit in CI.

## Reproduction and limits

The [packet and commands](https://github.com/kent-tokyo/chematic/tree/main/validation/v3000_atom_scalars)
include the original synthetic inputs and both compressed observations. The
candidate gate requires every positive field/source observation and every
field-specific refusal. Each external process is bounded to 10 seconds and 1 MiB
of combined output. This is a small source-build regression packet, not a
published-package result, throughput benchmark or Open Babel-wide superiority
claim.
