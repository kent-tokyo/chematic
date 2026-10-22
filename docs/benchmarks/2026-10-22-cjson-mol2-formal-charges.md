# Indexed formal charges in checked CJSON → MOL2 conversion

The checked CJSON → MOL2 adapter retained partial-charge columns but represented
formal charge only for +1 nitrogen (`N.4`). Anions, metals and other charges
were dropped, disabling verified graph correspondence even when the numeric
source formal charges were explicit. The adapter now writes every nonzero
projected formal charge as a `UNITY_ATOM_ATTR` charge attribute, using the
candidate atom's sequential ID.

```python
report = chematic.convert_with_report(source, "cjson", "mol2", policy="allow_loss")
# The entire cross-format result still has losses; inspect all observations.
assert report["candidate_graph_verified"]
formal = next(o for o in report["observations"] if o["field"] == "formal_charges")
assert formal["status"] == "preserved"
```

Rust, CLI, Python and WASM use the same adapter. The new observation group is
additive within the existing `chematic.conversion-report.v1` envelope. No API
signature, error enum or parser is added. Existing +1 nitrogen remains `N.4`
and also receives an explicit UNITY charge attribute. Neutral output has no
new section. ATOM partial charges and their declared method stay independent:
no charges are calculated or converted between formal and partial charge.
When there is no partial dataset, `NO_CHARGES` and absent partial-charge columns
remain unchanged even if the formal charge is nonzero.

## Verification and scope

The candidate is reparsed with the immutable MOL2 reader. The existing graph
check requires indexed elements, aromatic flags, formal charges, isotopes and
bond endpoints/orders. `formal_charges: preserved` requires this full projected
graph correspondence. It does not merely compare net charge. If, for example,
a source isotope cannot be emitted, the isotope is still `dropped`, formal
charges are `unverified` and mappings remain empty, even though the numeric
formal-charge attribute was written.

The overall result remains `completed_with_loss` under `allow_loss`, because
source representation, generated MOL2 identity/typing, opaque metadata and
stereo verification retain their existing limitations. `strict` still refuses
cross-format output with `output: null`. Source-preserving CJSON retention is
unchanged. Formal charges defaulted from an absent source array are still
separately reported as inferred defaults. The supported i8 storage range is
not a claim of chemically plausible charges throughout that range.

UNITY appears **before BOND**, matching the default MOL2 writer's compatibility
profile. The existing Open Babel [section-order issue #2350](https://github.com/openbabel/openbabel/issues/2350)
describes readers that instead require ATOM/BOND adjacency; this adapter does
not claim compatibility with those readers or with Antechamber. Choosing an
alternate order can change charge recognition and is outside this adapter's
current contract. Original atom typing, isotope retention, general MOL2
chemistry, conformers, protonation and docking preparation are not established.

## Paired packet

`validation/cjson_mol2_formal_charges/cases.json` contains eight original MIT
inputs: an oxygen anion, oxide, nitrogen cation, zinc ion, sodium/chloride pair,
carboxylate, explicit-H ammonium and neutral control. Each declares expected
indexed elements, formal charges, partial charges, coordinates and bonds, plus
provenance, units, minimization and SHA-256. LF and CRLF produce 16 variants.
Geometry is synthetic. The packet was used during implementation and is not
an independent holdout.

The baseline native Python extension was reused from the attested build in
PR #1133. That build's native sources are main `b9412a06`; its other changes
concern ordinary Python/WASM CJSON source retention, not this native checked
adapter. A fresh baseline Python package was copied from main and its extension
hash verified. Candidate bindings were built from this patch on the same base.
Source and binary identities are recorded in the validation index.

A separately compiled C++ probe linked to Open Babel **3.2.1** reads every
output as MOL2. It sets the partial-charges-perceived flag before querying
stored partial charges, preventing lazy charge calculation in that observer.
The probe reports indexed elements, per-atom formal/partial charge, coordinates
and bond endpoints/orders. Hand-authored source values are the oracle; neither
engine supplies expected values.

| Route, decoded by the external MOL2 probe | Indexed graph | Indexed formal charge | Partial charge bits | Coordinate bits |
| --- | ---: | ---: | ---: | ---: |
| Baseline checked CJSON → MOL2 | 16/16 | 6/16 | 16/16 | 16/16 |
| Candidate checked CJSON → MOL2 | 16/16 | 16/16 | 16/16 | 16/16 |
| Open Babel 3.2.1 `-icjson -omol2` | 14/16 | 4/16 | 0/16 | 14/16 |

The indexed-graph metric here compares elements and source-indexed bond
endpoints/orders independently of formal charge. Formal, partial and coordinate
metrics require that indexed-graph check first. Thus a failed coordinate metric
need not mean coordinates changed. All converter/probe invocations completed.
The Open Babel carboxylate output changed the source-indexed resonance
assignment (double bond and anionic oxygen exchanged), accounting for its two
indexed-graph failures; its raw coordinates still match. This does not establish
loss of chemical equivalence under resonance normalization. Other charged-ion
failures concern exact source-assigned formal charges. The partial-charge score
measures retention of the declared local CJSON dataset, not charge-model quality.

Native full graph verification, which also includes formal charges, improves
from 6/16 to 16/16. Python and real Node/WASM produce identical full reports for
all 16 candidates. Strict refusal, isotope-loss reporting and absence of invented
partial-charge columns are tested separately. This evidence is specific to this
local CJSON field profile and an explicit indexed-preservation workflow; it is
not toolkit-wide superiority, a general format guarantee, timing evidence or
a published release claim.

## Reproduction

Build a Python extension and Node-target WASM module from this source, and build
the original reader probe against Open Babel development headers/libraries:

```sh
c++ -std=c++17 $(pkg-config --cflags openbabel-3) \
  validation/cjson_mol2_formal_charges/openbabel_read_probe.cpp \
  $(pkg-config --libs openbabel-3) -o /tmp/cjson-mol2-ob-probe
PYTHONPATH=/path/to/python-package python3 scripts/validation/check_cjson_mol2_formal_charges.py \
  --wasm-module /path/to/chematic_wasm.js --openbabel-probe /tmp/cjson-mol2-ob-probe \
  --obabel /path/to/obabel --output /tmp/cjson-mol2-formal-charges.json
```

`--baseline` records earlier native output without enforcing candidate success;
it can omit the optional WASM module. All subprocesses have time/output bounds.
The observer tests reject equal-net-charge substitutions, atom reordering,
formal/partial-charge confusion, signed-zero loss and malformed probe rows.
Native tests additionally exercise every i8 formal-charge value without a
partial dataset and the output byte cap. CI discovers native/Python/WASM
regressions and runs the independent observer tests.

Full per-input reports, stdout/stderr, observations and hashes are retained in
`validation/results/cjson-mol2-formal-charges-{baseline,candidate}-2026-10-22.json.gz`.
The accompanying validation index identifies source/artifacts and compressed
build/test logs under `validation/results/cjson-mol2-formal-charges-2026-10-22/`.
