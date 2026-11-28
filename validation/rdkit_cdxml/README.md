# RDKit ChemDraw CDXML corpus

ChemDraw-saved `.cdxml` files copied unchanged from RDKit's
`Code/GraphMol/test_data/CDXML` (commit `0618b7055a6a42522b84d1acb1fd9f19f20b2ada`),
under RDKit's BSD 3-Clause license (see `LICENSE-RDKIT`).

`chemdraw_template2.cdxml` to `chemdraw_template9.cdxml` (about 1.2 MB) are
left out; `chemdraw_template1.cdxml` covers the same template-library shapes.
`bad-*.cdxml` and `malformed.cdxml` are intentionally invalid in RDKit's suite.

`expected_smiles.json` holds the expected SMILES RDKit's own
`cdxml_parser_catch.cpp` asserts for a subset of these files. The tests compare
them after canonicalizing both sides with chematic, so they check structure and
stereo, not RDKit's SMILES spelling.

Used by `crates/chematic-mol/tests/cdxml_rdkit_corpus.rs` (issues #1316,
#1317 and #1319).
