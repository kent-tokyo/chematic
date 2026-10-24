# CJSON to CML with checked graph, geometry and name

`project_cjson_to_cml` in Rust, Python and WASM exports the local ChemicalJSON
graph/Cartesian/name profile as CML. Python `convert_format`, WASM
`convert_common_format` and the CLI's ordinary `convert` command use this path
for CJSON → CML. The old route discarded source XYZ and names.

```python
import chematic

text = '{"name":"carbonyl","atoms":{"elements":{"number":[6,8]},"coords":{"3d":[0,0,0,1.25,0,1]}},"bonds":{"connections":{"index":[0,1]},"order":[2]}}'
cml = chematic.convert_format(text, "cjson", "cml")
assert cml == chematic.project_cjson_to_cml(text.encode("utf-8"))
```

The projector retains source-indexed elements, formal charges, isotope masses,
bond endpoints/orders, optional Cartesian coordinates in Angstrom, and a
consistent `name` / `properties.name`. Generated CML IDs `a1` through `aN`
follow the source array order. Finite binary64 coordinates retain their bits,
including signed zero; missing 3D remains missing. XML name text is escaped,
including carriage returns that XML would otherwise normalize. An explicit
empty name remains present.

Before returning any text, the projector reparses the candidate CML and checks
those indexed graph fields, coordinate presence/bits and name. It refuses
malformed CJSON, empty graphs, conflicting names, name characters outside XML
1.0, partial-charge datasets and every opaque field. This includes stereo
extensions, unit cells, 2D arrays and vendor metadata. Source information that
cannot be carried is rejected instead of dropped. Python rejects coordinate,
charge, name and comment overrides on this route. The dedicated Python API
accepts UTF-8 `str` or `bytes`; the common wrapper remains str-only. Rust also
offers `project_cjson_to_cml_with_limits` with `ConversionLimits`. Input, atom,
bond and output limits apply, including WASM's 1 MB input / 10,000-atom caps.

This is a selected-field projection. JSON spellings and the distinction between
the two name aliases are not preserved. The CML writer adds total H counts
using the existing valence inference. CML derives aromatic atom flags from
aromatic bonds; those flags are not explicit fields in the CJSON reader profile
and are outside the equivalence check. Chemical validity, inferred H validity
and coordinate-derived stereo are not verified. The report-based conversion
API retains its existing routes; it does not yet expose CJSON → CML reports.

In the original development packet, all three bindings retain all selected
fields on 24/24 common input variants. The three additional aromatic-4 cases
are separate local compatibility controls, since chematic interprets order 4
as aromatic and that convention is not shared universally. Open Babel 3.2.1
retains all selected common fields on 3/24 variants; individual graph/coordinate
results, missing outputs and name/isotope/charge coverage are reported
separately. These exposed synthetic tests do not establish toolkit-wide
superiority or chemical equivalence. See the
[complete packet, outputs and reproduction instructions](https://github.com/kent-tokyo/chematic/tree/main/validation/cjson_cml_projection).

The primary [ChemicalJSON description](https://github.com/OpenChemistry/chemicaljson/blob/main/chemicaljson.md)
and [Open Babel CML documentation](https://openbabel.org/docs/FileFormats/Chemical_Markup_Language.html)
describe the interchange formats and CML options. The accepted fields, local
aliases and refusals above are this projector's explicit contract.
