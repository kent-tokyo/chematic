# CML to CJSON preservation packet

Original synthetic single-molecule fixtures: 11 accepted profiles and 25 refusal
profiles, each with compact/LF/CRLF representation (33 + 75 evaluations). Neither
engine is the oracle. ElementTree decodes source fields and a separate typed JSON
observer checks emitted fields without using chematic readback.

Run `scripts/validation/check_cml_cjson_projection.py` with fresh Python/WASM/CLI
builds. `--baseline` records an earlier runtime without imposing candidate gates;
`--runtime-revision` labels its provenance. The earlier runtime is local
`ed050821`, not a release or main build. Gzip records keep complete input/output,
process diagnostics, commands, bounds, hashes and version. Observer controls live
in `scripts/tests/test_cml_cjson_projection.py`.

See [contract and measured limits](../../docs/reference/cml-cjson-projection.md).
Refusals are separate from selected-field preservation. General CML schema,
chemical validity, geometry-derived stereo, ID spelling and opaque source XML
are outside this comparison.
