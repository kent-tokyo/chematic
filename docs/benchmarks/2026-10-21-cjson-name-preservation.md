# CJSON molecule-name preservation, 2026-10-21

The checked CJSON reader previously treated root `name` as opaque and selected
only `properties.name`. A standard root title became the placeholder `chematic`
on CJSON → MOL2. The [ChemicalJSON schema](https://github.com/OpenChemistry/chemicaljson/blob/main/cjson.schema)
defines the optional molecule title at root `name`.

The source candidate recognizes both fields. One name, or two identical names,
resolves verbatim. Conflicting declarations remain in the immutable source,
with `name_conflict: true` and no selected name. MOL2 conversion refuses the
conflict even under `allow_loss`, while same-format conversion retains the
exact source. CJSON → SMILES explicitly reports name loss, including conflicts.
MOL2 → CJSON emits a root name and the matching legacy alias and verifies both.
The graph-only permissive CJSON reader remains compatible.

## Development packet

The packet consists of eleven original inputs, each containing one helium
atom to isolate title handling: root ASCII/Unicode names, a legacy alias,
matching dual declarations, an absent-name control, conflicting declarations,
empty/newline/boundary-space names and two reverse MOL2 → CJSON cases. Both
engines receive identical bytes. Expected titles come from the authored input;
neither engine supplies expected values. Integer graph/coordinate checks and
verbatim title checks have no tolerance. Source geometry is supplied, not
generated or optimized.

| Terminal outcome, all 11 inputs retained | Before name fix | Fixed candidate | Open Babel 3.2.1 default |
|---|---:|---:|---:|
| Required packet fields preserved | 5 | 7 | 3 |
| Output with required field loss | 2 | 0 | 4 |
| Explicit refusal | 0 | 4 | 0 |
| Output on conflicting/out-of-profile title input | 4 | 0 | 4 |
| Process failure / unresolved output | 0 | 0 | 0 |

The seven supported cases include the unnamed control, for which no title
preservation is claimed. Four declared conflicts/title-profile exclusions are
refusals, not successful conversions. Both root and legacy names count as
retained meaning in CJSON output; legacy-only storage is not counted as data
loss. The two reverse cases therefore pass for both engines before this fix.
Emitting the standard root field improves compatibility with root-only readers,
but is not scored as a separate semantic win over legacy output.

The pre-edit production baseline is `bb682cbd`, after incorporating current
main and repairing the charge-packet script location. Raw input text/hashes,
stdout/stderr, terminal outcomes and field scores are archived in the
[baseline](https://github.com/kent-tokyo/chematic/blob/main/validation/results/cjson-names-baseline-2026-10-21.json),
[candidate](https://github.com/kent-tokyo/chematic/blob/main/validation/results/cjson-names-candidate-2026-10-21.json)
and [validation index](https://github.com/kent-tokyo/chematic/blob/main/validation/results/cjson-names-validation-2026-10-21.json).

This is an authored development packet including the motivating failure, not
an independent holdout, published-package proof or toolkit-wide superiority.
General atom typing, isotope/stereo interpretation, implicit hydrogens and
arbitrary metadata are outside the packet. Strict cross-format conversion
continues to refuse because other source losses remain. Malformed name types
are rejected by the checked reader; valid conflicting strings remain available
through exact source retention.

## Reproduce

```sh
cargo test --locked -p chematic-mol --test cjson_names
cargo build --locked -p chematic-cli
python3 -m scripts.validation.check_cjson_names \
  --cli target/debug/chematic --obabel /opt/homebrew/bin/obabel \
  --output /tmp/cjson-names-new.json
python3 -m pytest -q scripts/tests/test_cjson_names.py
```

The comparator is pinned to Open Babel 3.2.1 and new result paths are required.
The source fixtures and recipe are original project contributions under
`MIT OR Apache-2.0`; provenance, minimization and normalized input hashes are
included in the report. The previous charge-packet script now lives under
`scripts/validation/`; its original archived measurements remain unchanged.
