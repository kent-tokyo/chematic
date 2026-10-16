# Explicit ChemDraw block formats

RDKit [#9604](https://github.com/rdkit/rdkit/issues/9604) requests explicit
AUTO/CDX/CDXML reader selection. Its linked
[upstream proposal](https://github.com/rdkit/rdkit/pull/9597) describes binary
CDX being routed through XML parsing. RDKit 2026.09.1 reproduces this on the
authored ethanol CDX fixture: the block reader returns no molecules, the file
reader raises an XML parsing error, and the explicit `format` keyword is absent.

The new native API adds a bounded binary CDX molecular profile and explicit
dispatch. Existing `from_cdxml` and other CDXML entry points are unchanged.

```python
result = chematic.from_chemdraw_block(binary_bytes, format="AUTO")
molecules = result["molecules"]
coordinates = result["coordinates"]
assert result["source_bytes"] == binary_bytes

# Consent to listed presentation omissions, when appropriate:
result = chematic.from_chemdraw_block(
    binary_bytes, format="CDX", allow_metadata_loss=True
)
print(result["omissions"])  # [{"offset": byte_offset, "tag": numeric_tag}, ...]

# Explicit XML uses UTF-8 bytes too:
result = chematic.from_chemdraw_block(xml_text.encode("utf-8"), format="CDXML")
```

Rust exports `parse_chemdraw_block` / `parse_chemdraw_block_with_limits` and
`ChemDrawFormat`, `ChemDrawOptions`, `ChemDrawLimits`, `ChemDrawBlock`,
`ChemDrawError` from `chematic_mol`. Formats are `Auto`, `Cdx`, `Cdxml`;
`ChemDrawFormat::from_name` accepts their case-insensitive names.
Errors use a new non-exhaustive enum, preserving existing public error enums.

WASM takes `Uint8Array` input and exposes:

- `mols_from_chemdraw_block(bytes, format, allowMetadataLoss)`: source-ordered
  molecule handles for nonempty fragments.
- `chemdraw_block_report_json(bytes, format, allowMetadataLoss)`: detected
  format, exact `source_bytes` as a byte array, `projected_cdxml`,
  `molecule_count`, corresponding `coordinates`, and CDX `omissions`.

Both operations use the same native decoder. Reading handles and a report
separately decodes twice; this is not a throughput API. Python returns both
in one dict, using `molecules` instead of `molecule_count`.

## Binary CDX profile

The implementation uses the vendor specification preserved by
[IUPAC FAIRSpec](https://iupac.github.io/IUPAC-FAIRSpec/cdx_sdk/):
[tagged binary framing](https://iupac.github.io/IUPAC-FAIRSpec/cdx_sdk/IntroCDX.htm),
[nodes](https://iupac.github.io/IUPAC-FAIRSpec/cdx_sdk/Node.htm) and
[bond orders](https://iupac.github.io/IUPAC-FAIRSpec/cdx_sdk/properties/Bond_Order.htm).
It supports little-endian fields and short/extended property lengths. The
22-byte header seen in RDKit's `structure_1.cdx` and the documented 28-byte
header layout are accepted with their reserved bytes checked. A document
object, balanced object endings, final file terminator and no trailing bytes
are required. AUTO detects only the exact eight-byte CDX signature; everything
else takes the checked UTF-8 CDXML path. Forced formats never fall back.

Supported CDX objects are document, page, group, fragment, element node and
bond. Nested molecular abbreviations, reaction schemes, queries, radicals,
enhanced stereo, custom properties and unknown objects are refused even when
metadata loss is enabled. Atom/bond IDs must be positive and child IDs unique
within their container. References must resolve inside their fragment.

Atoms require 2D positions. Source vertical/horizontal signed 16.16 values
become `[x,y]` in ChemDraw point units with Y increasing downward. Coordinates
and atoms retain source order; they are not converted to Å. Element nodes
default to carbon when Element is absent. Explicit elements 1–118, signed
8-bit charges, isotope masses and 0–255 explicit hydrogen counts are supported.
Absent H counts use the existing native valence model, not a new CDX sanitizer.

Single, double, triple and 1.5 aromatic bond orders are supported. Aromatic
edges mark their incident atoms aromatic; there is no general aromaticity
perception. Query combinations, dative/fractional orders outside this profile,
and non-solid `Display`/`Display2` values are refused. In particular, CDX
wedge/hash stereo is outside this first profile. Geometry-derived E/Z uses
the existing checked 2D direction assignment, so SMILES/MOL output can retain
supported double-bond configurations; tied/ambiguous cases abstain.

Known presentation properties and text objects require
`allow_metadata_loss=True`. Every omitted CDX tag is listed with its byte
position and the exact source remains available. Unknown chemical fields
cannot be bypassed by that flag. `projected_cdxml` is an intermediate molecular
projection, not a vendor document round trip or a replacement for source bytes.

## XML policy and limits

CDXML is checked for one balanced root and bounded nesting/items before
delegating to the existing molecular reader. The new API also completes E/Z
direction metadata from the source geometry. It retains that reader's
graph-only treatment of presentation and unhandled attributes. `omissions`
lists **CDX tags only**; an empty list on XML is not proof that XML metadata
was interpreted or preserved in the molecules. Original XML bytes are retained.

Rust/Python defaults are 16 MiB input, 200,000 binary/XML items, depth 32,
1 MiB per binary property and 32 MiB projected XML, plus existing
`CdxmlParseLimits` for atoms/bonds/fragments/lines/attributes. Node budgets are
per fragment; the global item/input caps bound a multi-fragment document.
WASM additionally caps input at 1,000,000 bytes and JSON output at 16 MiB.
Limits and malformed/unsupported input return errors without a partial result.

## Evidence

`scripts/generate_chemdraw_fixtures.py` creates 56 original binary recipes with
independently assembled CDXML counterparts: seven molecular profiles, two
header lengths, reversed property ordering and optional presentation data.
They include charged/isotopic atoms, explicit H counts, aromaticity, E/Z,
multiple fragments and reused local IDs across fragments.

Rust/Python/real Node-WASM tests cover dispatch, source bytes, coordinates,
metadata consent, unsupported-chemistry refusal, every truncated prefix,
extended lengths, references, duplicate IDs/properties and resource boundaries.
`scripts/validate_chemdraw_block.py` uses RDKit 2026.09.1's independent XML
reader as the graph oracle and retains all inputs/outputs plus source/runtime
hashes. It compares binary and XML AUTO/explicit modes, including indexed atom
identities, charge, isotope, H count, bonds and canonical isomeric identity.
The original binary-dispatch failure is recorded separately.

The local packet in `validation/results/chemdraw-block-2026-10-16` records
**224/224** Python and **224/224** WASM graph/identity comparisons, with exact
source-byte, coordinate and projection-report agreement. These expand the
56 binary recipes across forced/AUTO binary and forced/AUTO XML modes.

This is a local concrete-molecule candidate. Full CDX/ChemDraw support,
reaction/query reading, vendor write interoperability, speed and public-package
publication are not established by these fixtures.
