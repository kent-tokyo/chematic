# CML namespace-aware molecular projection

The Cartesian reader and CML→CJSON bridge select chemistry by XML namespace
and molecular scope. Changing a CML prefix does not change the graph. A foreign
element named `atom` does not add an atom, and a foreign attribute named `x3`
does not override the source coordinate.

Previously, `from_cml_with_3d` refused a normal `c:molecule` document but accepted
an unprefixed foreign `atom` inside the same file as an extra chemical atom.
The prefix and annotation spellings therefore changed the projected molecule.
The shared native reader now uses the existing bounded XML dependency to resolve
expanded names before graph construction.

```python
import chematic

source = '''<c:molecule xmlns:c="http://www.xml-cml.org/schema" xmlns:ui="urn:ui">
  <c:atomArray>
    <c:atom id="oxygen" elementType="O" x3="1" y3="2" z3="3"/>
    <ui:atom id="label" elementType="C" x3="9" y3="8" z3="7"/>
  </c:atomArray>
</c:molecule>'''
mol, xyz = chematic.from_cml_with_3d(source)
assert mol.atom_count == 1
assert xyz == [[1.0, 2.0, 3.0]]
```

The same selection is used by Rust `parse_cml_with_3d[_with_limits]`, Python
`convert_format(source, "cml", "cjson")`, and WASM/Node
`convert_common_format(source, "cml", "cjson")`. No new public symbol or return
shape is introduced. Legacy graph/2D readers keep their existing profile.

## Local projection rules

- Accept one molecule root, or one direct molecule inside a CML wrapper.
  CML elements use `http://www.xml-cml.org/schema`, with any bound prefix or
  a default namespace. Legacy unnamespaced documents are also accepted.
- Chemical atom/bond records belong to direct `atomArray`/`bondArray` children
  in that molecule's namespace. Direct atom/bond children remain a legacy
  extension. Formula composition and other metadata records do not become
  molecular graph atoms. Foreign subtrees remain excluded even if they contain
  descendants that reintroduce the CML namespace.
- Only unqualified chemical attributes supply atom IDs, elements, coordinates
  and bond data. A qualified metadata attribute cannot replace a chemical field.
  XML character references are decoded once by the XML parser.
- Atom IDs are unique/nonempty; `elementType` is required. Bond `atomRefs2`
  contains exactly two different existing atom IDs and `order` is required.
  Duplicate graph bonds are rejected. Bond arrays may precede atom arrays;
  source atom order and declared endpoint order are retained.
- Complete finite x3/y3/z3 coordinates are returned in angstroms in source
  order; absent geometry stays absent. Partial/mixed triplets, packed and
  fractional coordinates remain outside this profile.
- XML must be well formed. DTD processing is disabled. All elements, including
  excluded annotations, count toward XML budgets. Existing input, physical and
  normalized line, atom and bond limits apply; depth is capped at 64 and the
  derived XML-node budget at four million.

These are bounded molecular-projection rules, not complete CML schema validation.
The [CML molecular convention](https://www.xml-cml.org/convention/molecular)
defines namespace identity, required atom elements and two-atom bond references;
the [namespace examples](https://www.xml-cml.org/documentation/basic-concepts.html)
show equivalent prefixed elements and mixed foreign metadata.

Opaque properties, namespaces as source text, 2D layout and stereochemical
declarations are not retained in the CJSON projection. Use the checked
same-format source-retention route when the original XML document must survive.
Its exact-source contract is separate from this interpreted graph subset.
The namespace fix does not extend CJSON's representable bond types; query-bond
and stereochemical round trips remain outside this comparison packet.

## Evidence

The original synthetic packet `validation/cml_namespace_projection/cases.json`
contains 20 accepted profiles and 28 refusal profiles, each with compact, LF
and CRLF layout. This gives 60 positive and 84 negative cases, not 144 distinct
molecules. Positives cover prefix changes, namespace resets/rebinding, foreign
atom/bond/attribute annotations, formula scopes, forward references, character
references and labelled ions. Negatives cover XML/DTD/depth failures, duplicate
expanded attributes or decoded IDs, missing required fields, invalid endpoints
and unsupported coordinate declarations.

Python's independent ElementTree parser checks authored source graph fields and
coordinate bits. A separate CJSON observer then checks the emitted indexed graph,
formal charges, isotope annotations and Cartesian bits without relying solely
on chematic readback. Python and real Node/WASM pass every positive and refusal
case. Native tests also verify the graph, endpoint order, coordinates and bounds.

The dated packet `validation/cml-namespace-projection-2026-10-24/` retains local
baseline/candidate binary identities, source identities, comparator commands, all outputs
and refusals. Open Babel 3.2.1 is measured on the same input and CJSON target,
with both default output and `-xs`. Comparator acceptance of a negative case
means acceptance outside this local profile; it does not automatically establish
a CML-specification violation. This packet does not establish full CML, stereo,
opaque-metadata, performance or toolkit-wide superiority.
