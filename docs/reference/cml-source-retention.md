# Keep CML source through ordinary conversion

`convert_format(text, "cml", "cml")` in Python, `convert_common_format` in WASM
and `chematic convert --input-format cml --output-format cml` now return the
exact input text after bounded XML/document checks. This avoids projecting a
document into a molecule and rewriting it, which discarded source geometry,
IDs and extensions. No final newline is added or removed.

```python
import chematic

source = open("input.cml", encoding="utf-8", newline="").read()
assert chematic.convert_format(source, "cml", "cml") == source
assert chematic.roundtrip_cml_source(source.encode("utf-8")) == source
```

The dedicated `roundtrip_cml_source` API is available in Rust, Python and WASM.
Python's dedicated API accepts UTF-8 `str` or `bytes`; the common conversion
wrapper retains its existing `str` input contract. Rust also provides
`roundtrip_cml_source_with_limits` with `ConversionLimits`.

This reuses the retention profile already available through
`convert_with_report(..., policy="preserve_source")`. That report continues to
state `verification_scope="exact_source_bytes"`,
`candidate_graph_verified=false` and `source_semantics="opaque_retained"`.
The new ordinary path returns the source without allocating/encoding a JSON
report, and checks its output budget against the raw UTF-8 text.

The profile requires one well-formed XML 1.0 document whose root is a
`molecule` in `http://www.xml-cml.org/schema`. Namespace prefixes are allowed.
Nested CML molecules, DTDs, unsupported XML declarations/encodings, malformed
later XML and resource-limit violations are refused before returning output.
The XML parser has a 200,000-node and 64-element-depth cap; shared conversion
limits also bound bytes, lines, line length and CML atom/bond element counts.
WASM has its usual 1,000,000-byte input and 10,000-CML-atom limits. Python rejects
coordinate, charge, name and comment overrides on this source-retention route.
CLI failures leave stdout empty and existing output files untouched.

Coordinates, number spellings, IDs, namespaces, comments, CDATA, property order,
unknown attributes and stereo declarations are retained in their source form.
The operation does **not** verify chemical validity, schema conformance, atom
references or stereo meaning. For example, an unresolved bond reference can be
retained opaquely. A successful retention is not a verified molecular parse.
Unqualified molecule roots and CML wrapper/hierarchy documents accepted by
other graph readers are outside this profile and now cause refusal in ordinary
CML self-conversion. Graph parsing/editing and cross-format conversion keep
their separate contracts.

On six original development documents × three newline forms, Python, WASM and
CLI each retained exact source and all independently observed selected fields
in 18/18 evaluations. Open Babel 3.2.1 retained those selected fields in 9/18
evaluations in both default and `-xp -xm` modes. Losses involved simultaneous
2D/3D coordinates, high-precision numbers and opaque attributes/properties;
this is a bounded source-retention comparison, not toolkit-wide superiority.
See the [packet and complete results](https://github.com/kent-tokyo/chematic/tree/main/validation/cml_common_retention).

Open Babel's [CML documentation](https://openbabel.org/docs/FileFormats/Chemical_Markup_Language.html)
describes its recognized field subset, dimensionality selection and metadata/
property options. This retention profile serves the workflow of carrying a
document intact without interpreting all its contents.
