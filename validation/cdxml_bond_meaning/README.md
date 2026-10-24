# CDXML bond meaning packet

Original synthetic MIT fixtures authored on 2026-10-24 from the archived
ChemDraw SDK: 4 graph controls, 3 quadruple cases, 5 directed-dative cases and
12 explicit profile refusals. XML well-formedness is checked; no schema,
ChemDraw-application or physical-validity claim is made.

`cases.json` supplies the authored atom/bond/direction oracle. `baseline.json`
and `candidate.json` retain every public-route output or error, runtime/source
identity and independent observation. SMILES is observed with RDKit, CDXML
output with ElementTree. Numeric order one alone cannot identify a dative bond.
Unsupported-input output is classified separately from preservation/refusal.

Frozen baseline runtime: `342d1d5f9ad09aa6ade4e3079a902cba856209e2`; its CDXML
decoder/attribute implementation match main base `08c2da71`. Candidate source
hashes identify the changed snapshot independently of the pre-commit HEAD.
See [scope, counts and reproduction](../../docs/reference/cdxml-bond-meaning.md).
