# Checked molecular CIF bond output

Rust and Python expose `write_cif_with_bonds_checked(mol, coords, cell=None)`
alongside the existing atomic-position writer. In Rust, `cell` is
`Option<&UnitCell>` and the return type is `Result<String, CifBondWriteError>`.
Python raises `ValueError` on invalid inputs.

```python
import chematic

mol = chematic.from_smiles("[H]O[H]")
text = chematic.write_cif_with_bonds_checked(
    mol, [[0.9, 0, 0], [0, 0, 0], [0, 1.1, 0]]
)
```

## Output contract

- One finite Cartesian Å triple per graph atom, in graph atom order. Optional
  cell validation and fractional conversion match `write_cif_checked`.
- Atom labels use the existing writer's element symbol and one-based graph
  index. Bond rows retain graph edge order and endpoint orientation and refer
  to those exact labels.
- `_geom_bond_distance` uses the supplied Cartesian positions, with stable
  `hypot` evaluation and round-trip float output. Finite strictly positive
  distances are required, including for subnormal and very large inputs.
- Both `_geom_bond_site_symmetry_1` and `_geom_bond_site_symmetry_2` are `.`.
  No minimum-image wrapping, periodic connectivity or symmetry expansion is
  performed. A 9.8 Å bond stays 9.8 Å in a 10 Å cell.
- `_ccdc_geom_bond_type` is S, D, T or A. Dative, zero-order, quadruple,
  directional and query orders reject instead of being converted to single
  bonds. Edgeless graphs produce no bond loop; bonds are never guessed.

This exports positions and explicit connectivity. Formal charges, isotopes,
implicit hydrogens, stereo annotations, source names, occupancies and other
metadata are outside this projection. It is not a full chemical round trip.
`parse_cif` still returns positions with **no bonds**. No plain-CIF WASM or MCP
API is added. The periodic `crystal` adapter has its own contract.

## Open Babel issue boundary

[Open Babel #2726](https://github.com/openbabel/openbabel/issues/2726) was still
open when inspected on 2026-10-16. Its comments clarify that write options use
`-xg`, rather than the report's `-ag` read option. The
[Open Babel CIF documentation](https://openbabel.org/docs/FileFormats/Crystallographic_Information_File.html)
already describes bond writing with `-g` and bond reading with `-B`.
This is a previously missing **chematic** feature, not evidence of unresolved
Open Babel bond-output support or an upstream bug fix.

## Validation

Hermetic Rust and Python tests cover labels, explicit H, S/D/T/A distinctions,
edgeless/disconnected graphs, cell-independent distances, extreme representable
norms, invalid inputs and unsupported orders. The v1.1.3 exhaustive public
error contract is retained by using a separate non-exhaustive writer error.

`scripts/validation/validate_cif_bonds.py` parses the output independently with Gemmi
0.7.5, checking atom references, coordinate conversion, edge sets, bond types,
symmetry markers and Cartesian distances. It also records bounded Open Babel
CLI probes with both `-ag` and `-xg`; CIF read timeouts are reported separately.
The original fixture manifest and dated source-candidate receipt are under
`validation/openbabel_issues/cif-bonds/` and `validation/results/`.

Local source-candidate results on base `52328da3` (2026-10-16): 22 Rust tests
passed (8 new bond tests, 8 existing checked-writer tests, 5 coordinate tests,
1 v1.1.3 API contract); 63 Python tests passed (18 new and 45 existing).
Gemmi checked 372 items across 21 cases with zero failures. Open Babel 3.2.1
wrote no bond loop with `-ag` and did write one with `-xg`. Its CIF reading
probe exceeded five seconds, so external chemical read-back is inconclusive.
These are local source results, not a published package claim.
