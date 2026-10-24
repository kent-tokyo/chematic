# CDXML bond class and direction

The CDXML graph reader and writer now retain quadruple bonds and directed
dative bonds. Previously, both could become ordinary single bonds. The reader
also refuses declared bond orders outside its supported profile instead of
silently treating them as single.

The profile accepts an absent `Order` as single, or explicit `1`, `2`, `3`, `4`,
`1.5` (aromatic) and `dative`. These names follow the
[archived ChemDraw SDK Bond_Order reference](https://iupac.github.io/IUPAC-FAIRSpec/cdx_sdk/properties/Bond_Order.htm).
A dative edge retains its `B` → `E` direction. `WedgeEnd`/`WedgedHashEnd` can
reverse a single bond's stereo reference endpoint; they do not reverse an
explicit dative edge. The native writer emits `Order="4"` and
`Order="dative"` while keeping the stored dative endpoints.

Fractional orders other than `1.5`, orders five/six and ionic, hydrogen or
three-center declarations are outside this graph profile. They return
`CdxmlError::UnsupportedBondOrder` in Rust and an error in Python, WASM and CLI
common-format conversion. An empty or unrecognized declaration also refuses;
only absence invokes the single-bond default. This is not full CDXML validation.

## Bounded development comparison

The original MIT packet has four single/double/triple controls, three quadruple
cases, five directed-dative cases and twelve required refusals. Every input is
well-formed XML checked with ElementTree. The graphs are authored examples, not
engine-generated chemical or physical-validity oracles. Node-order and endpoint
reversal cases check atom correspondence; display variants check that drawing
decoration cannot replace the explicit bond class or dative direction.

| Route | Quadruple, baseline → candidate | Dative, baseline → candidate |
|---|---:|---:|
| Python direct reader, stereo option false/true (each) | 0/3 → 3/3 | 0/5 → 5/5 |
| Python common CDXML → SMILES | 0/3 → 3/3 | 0/5 → 5/5 |
| Python CDXML → CDXML graph round trip | 0/3 → 3/3 | 0/5 → 5/5 |
| Python writer from authored SMILES | 0/3 → 3/3 | 0/5 → 5/5 |
| WASM equivalents of these routes (each) | 0/3 → 3/3 | 0/5 → 5/5 |
| Open Babel 3.2.1 CDXML → SMILES | 3/3 | 0/5 |

All four controls preserve the observed graph in every measured route, before
and after. Candidate reader/common routes refuse all twelve unsupported or
unrecognized declarations. The baseline and Open Babel emit output for all
twelve; those outputs are recorded as outside this profile, not assigned an
invented chemical oracle. Ten cases use SDK-defined orders that the bounded
core cannot represent; two use empty/unrecognized declarations.

The Open Babel command is `obabel -icdxml -osmi -xn -xa -l1`. Stopping after the
first fragment matches the direct reader's scope and avoids a subsequent EOF
parse diagnostic in this single-fragment packet. Every process is bounded to
10 seconds and 1 MiB combined output. All measured Open Babel invocations
complete successfully. Its documented [CDXML support](https://openbabel.org/docs/FileFormats/ChemDraw_CDXML_format.html)
is limited to chemical structure information.

SMILES outputs are observed with RDKit 2025.09.3, with sanitization and H-node
removal disabled. RDKit supplies no CDXML interpretation or chemical oracle.
The observer explicitly distinguishes `DATIVE` from `SINGLE`, even though both
have numeric order one, and checks donor/acceptor orientation. Unique elements
identify this packet's two atoms; charge, isotope and non-node H counts must
also match. Writer and CDXML round-trip outputs are inspected independently
with ElementTree, without reusing chematic's CDXML reader.

The frozen baseline runtime is
`342d1d5f9ad09aa6ade4e3079a902cba856209e2`; its CDXML decoder and shared XML
attribute implementation match main base `08c2da71`. The full before/after
outputs, errors, observations and runtime/source hashes are retained in the
[packet](https://github.com/kent-tokyo/chematic/tree/main/validation/cdxml_bond_meaning).
The recorded HEAD is the pre-commit checkout; candidate hashes identify the
changed source. These are source-candidate results, not a published-package or
general Open Babel ranking.

## Scope and reproduction

Native tests additionally inspect `BondOrder`, endpoints and input-indexed 2D
coordinates, and test writing from independently parsed authored SMILES. CLI
tests check successful graph conversion and empty stdout on unsupported orders.
Python and Node regressions cover both stereo-option settings. Seven observer
negative controls run without RDKit in CI.

This change covers the listed graph classes. It does not establish general
aromaticity/stereo, query, coordinate-conversion or source-byte/layout retention.
The legacy writer's behavior for other unrepresentable core bond classes is
unchanged and outside this packet; this is not a new general checked writer.
ChemDraw-application compatibility and speed are not measured.

```sh
PYTHONPATH=/absolute/python-package-root python3 scripts/validation/check_cdxml_bond_meaning.py \
  --wasm-module /absolute/pkg-node/chematic_wasm.js --output /tmp/cdxml-bonds.json
python3 -m unittest scripts.tests.test_cdxml_bond_observer -q
cargo test -p chematic-mol --test cdxml_bond_meaning
cargo test -p chematic-cli --test cdxml_bond_meaning
```
