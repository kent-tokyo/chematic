# MOL2 numeric preservation

The native editable/generated MOL2 writers retain every supplied finite
coordinate and partial-charge value as binary64, including negative zero.
They keep the existing four-decimal spelling when it already reparses to the
same bits; otherwise they emit enough decimal digits to recover the value.
For example, `1.2345678901234567` stays that value instead of becoming `1.2346`.
V2000 fixed-width formatting is unchanged.

```rust
let mut record = chematic_mol::parse_mol2_record(source)?;
record.coords[0].0 = 1.2345678901234567;
record.atoms[0].partial_charge = Some(-0.12345678901234566);
let output = chematic_mol::write_mol2_record(&record);
```

This applies to `write_mol2_record` and its section-order variant, `write_mol2`
and its section-order variant, and `mol2_source::write_mol2_source`.
Graph-only writers still generate `NO_CHARGES` and placeholder zero partial
charges. They do not transfer a source charge dataset, names or SYBYL types.
Missing coordinate rows still default to zero in these existing APIs.
This change adds no validation of arbitrary edited metadata or non-finite values.

Immutable source retention (`write_mol2_source_record`, Python/WASM
`roundtrip_mol2_record`, ordinary MOL2 self-conversion) already retains original
text and remains separate. The new measurement exercises the native regenerating
writers; it does not claim a new Python/WASM coordinate-output API or CLI route.

## Measured scope — 2026-10-24

Eight original synthetic inputs isolate exact four-decimal controls, longer
coordinate mantissas, tiny displacements, charge precision, signed zero, large
coordinates, rounding boundaries and a two-atom record. Both engines rewrite
MOL2 to MOL2. The native candidate uses parse → editable-record write;
Open Babel 3.2.1 uses `obabel -imol2 -omol2` with default options.
All eight sources and outputs are also read with RDKit 2026.09.1, with
`sanitize=False`, `removeHs=False`, `cleanupSubstructures=False`.
Expected values come from authored source tokens, independently of either
engine. The observer checks exact bits and indexed elements/bond endpoints/orders.

| Regenerating path | Indexed graph | Coordinates, exact bits | Partial charges, exact bits | All three |
|---|---:|---:|---:|---:|
| Native candidate, either section order | 8/8 | 8/8 | 8/8 | 8/8 |
| Previous native `aeae3f53` | 8/8 | 3/8 | 5/8 | 2/8 |
| Open Babel 3.2.1, default MOL2 rewrite | 8/8 | 3/8 | 5/8 | 2/8 |

All eight source records and the candidate/Open Babel outputs are readable by RDKit.
Both native graph-only writers also retain coordinates in 8/8 cases; source
partial charges are outside their contract and are not credited as preserved.
Three additional native-only cases cover the smallest subnormal, smallest normal
and largest finite binary64 values. They are excluded from the external score.
Native tests also cover edited values, repeated regeneration and stratified
finite bit patterns, while retaining the older wide-name/Unicode layout checks.

The packet is a small numerical regression diagnostic, not an independent
holdout, chemical-validity test, performance result or general Open Babel ranking.
The external reader check validates these tokens in RDKit, not every downstream
consumer; fixed-column consumers may impose additional width restrictions.
No atom typing, stereochemistry or unsupported metadata guarantee is added.
See the [Open Babel MOL2 format options](https://openbabel.org/docs/FileFormats/Sybyl_Mol2_format.html)
for the comparator's documented options.

Fixtures, runner, provenance and compressed raw records are in
`validation/mol2_numeric_preservation/`,
`scripts/validation/check_mol2_numeric_preservation.py` and
`validation/results/mol2-numeric-preservation-{baseline,candidate}-2026-10-24.json.gz`.
