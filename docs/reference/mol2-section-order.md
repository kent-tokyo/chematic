# MOL2 section-order compatibility

The editable/generated MOL2 writers offer an explicit `AtomBondAdjacent`
policy for readers expecting ATOM immediately followed by BOND, as requested
in [Open Babel #2350](https://github.com/openbabel/openbabel/issues/2350).
The default keeps the existing ATOM → UNITY_ATOM_ATTR → BOND order and bytes.

This is a reader-specific choice. On the nine original controls, Open Babel
3.2.1 preserved formal charges with the default order. With adjacent ATOM/BOND,
it ignored the trailing UNITY formal charges in all six charged graph controls.
The rich `N.4` record still read as +1, consistent with its atom-type annotation.
Changing the default would therefore regress a measured interchange route.

## Usage

Python exposes a keyword-only option:

```python
mol = chematic.from_smiles('[O-2]')
normal = mol.to_mol2()
adjacent = mol.to_mol2(atom_bond_adjacent=True)
```

JavaScript/WASM accepts an optional second boolean:

```javascript
const normal = wasm.smiles_to_mol2('[O-2]');
const adjacent = wasm.smiles_to_mol2('[O-2]', true);
```

Rust provides `write_mol2_with_section_order` and
`write_mol2_record_with_section_order` with `Mol2SectionOrder::AtomBondAdjacent`.
The existing `write_mol2` and `write_mol2_record` retain their behavior.
No CLI option is added. No exports are added to Python or WASM.

The option moves UNITY after BOND without dropping its attributes. It preserves
formal-charge text, partial-charge values, coordinates, source atom IDs/types,
residue/status fields and opaque sections in the editable record profile.
Generated graph output continues to use element-only types, zero coordinates
and placeholder partial charges; this does not assign SYBYL types or calculate
charges. It does not make every consumer understand trailing UNITY attributes.

Immutable source writers, `roundtrip_mol2_record`, and same-format MOL2
conversion retain original text and section order, including LF/CRLF.

## Bounded validation — 2026-10-16

Baseline source: `3e7f04b20593af4b07569d886b56afc1184468fe`.
The original fixtures cover six charged graphs, two neutral graphs and one
editable record with nonsequential IDs, a partial charge, unknown attributes
and an opaque section.

- 27 Rust integration checks and 15 MOL2 unit checks passed.
- 57 Python checks passed, including 11 section-order/preservation checks.
- 11 WASM section-order checks and the existing MOL2 preservation script passed.
- All nine default native outputs equal the baseline bytes. For the opt-in
  outputs, every section's contents equal the baseline; seven outputs only
  change section order, and the two neutral outputs remain byte-identical.
- Eight graph outputs match byte-for-byte across native Rust, Python and WASM.
- Open Babel 3.2.1 read both profiles with the expected atom/bond counts in
  nine cases. All nine default formal-charge vectors match. Six adjacent
  profile charge vectors differ, as disclosed above.

The runnable evidence is in
`validation/openbabel_issues/mol2-section-order/validate.py`; saved native,
Python, Node and external-read outputs are in
`validation/results/mol2-section-order-2026-10-16/`.
Antechamber is not installed and was not executed. This does not establish a
complete Amber workflow, general chemical equivalence, a speed gain or a
published package change.
