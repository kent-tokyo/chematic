"""Tests for 2D wedge/hash stereo perception + diagnostics wired into the
MOL/SDF readers (from_mol_block_with_diagnostics, from_mol_v3000_with_diagnostics,
SdfRecord.stereo_diagnostics()).

See docs/rfcs/stereo2d_reader_integration_rfc.md for the design background and
crates/chematic-mol/tests/stereo_reader_integration.rs for the Rust-side
equivalents of these same fixtures.
"""
import chematic


def _atom_line(x, y, sym):
    return f"{x:>10.4f}{y:>10.4f}{0.0:>10.4f} {sym:<3} 0  0  0  0  0  0  0  0  0  0  0  0"


def _bond_line(a1, a2, btype, stereo):
    return f"{a1:>3}{a2:>3}{btype:>3}{stereo:>3}"


def _chfclbr_v2000(wedge_bonds):
    """A V2000 CHFClBrI block. `wedge_bonds` maps 1-based substituent index
    (2=F, 3=Cl, 4=Br, 5=I) to an MDL stereo code (1=up, 6=down)."""
    lines = [
        "test",
        "  chematic",
        "",
        "  5  4  0  0  0  0  0  0  0  0999 V2000",
        _atom_line(0.0, 0.0, "C"),
        _atom_line(-1.0, 0.4, "F"),
        _atom_line(0.9, 0.7, "Cl"),
        _atom_line(-0.5, -1.1, "Br"),
        _atom_line(0.8, -0.6, "I"),
    ]
    for sub in (2, 3, 4, 5):
        lines.append(_bond_line(1, sub, 1, wedge_bonds.get(sub, 0)))
    lines.append("M  END")
    lines.append("")
    return "\n".join(lines)


VALID_WEDGE_BLOCK = _chfclbr_v2000({2: 1})
CONTRADICTORY_WEDGE_BLOCK = _chfclbr_v2000({2: 1, 3: 1})

V3000_EITHER_BLOCK = """either
  chematic

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0 0
M  V30 2 F -1.0 0.4 0 0
M  V30 3 Cl 0.9 0.7 0 0
M  V30 4 Br -0.5 -1.1 0 0
M  V30 5 I 0.8 -0.6 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2 CFG=2
M  V30 2 1 1 3
M  V30 3 1 1 4
M  V30 4 1 1 5
M  V30 END BOND
M  V30 END CTAB
M  END
"""

V3000_WEDGE_BLOCK = """wedge
  chematic

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0 0
M  V30 2 F -1.0 0.4 0 0
M  V30 3 Cl 0.9 0.7 0 0
M  V30 4 Br -0.5 -1.1 0 0
M  V30 5 I 0.8 -0.6 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2 CFG=1
M  V30 2 1 1 3
M  V30 3 1 1 4
M  V30 4 1 1 5
M  V30 END BOND
M  V30 END CTAB
M  END
"""


def test_valid_wedge_produces_no_diagnostics_and_defined_chirality():
    mol, name, coords, diagnostics = chematic.from_mol_block_with_diagnostics(
        VALID_WEDGE_BLOCK
    )
    assert diagnostics == []
    assert "@" in mol.smiles


def test_contradictory_wedge_produces_diagnostic_and_no_chirality():
    mol, name, coords, diagnostics = chematic.from_mol_block_with_diagnostics(
        CONTRADICTORY_WEDGE_BLOCK
    )
    assert diagnostics == [{"atom_idx": 0, "reason": "contradictory_wedges"}]
    assert "@" not in mol.smiles


def test_v3000_wedge_matches_v2000():
    mol, _, _, diagnostics = chematic.from_mol_v3000_with_diagnostics(V3000_WEDGE_BLOCK)
    assert diagnostics == []
    mol_v2000, _, _, _ = chematic.from_mol_block_with_diagnostics(VALID_WEDGE_BLOCK)
    assert mol.smiles == mol_v2000.smiles


def test_v3000_either_cfg_does_not_produce_a_wedge():
    mol, _, _, diagnostics = chematic.from_mol_v3000_with_diagnostics(V3000_EITHER_BLOCK)
    assert diagnostics == []
    assert "@" not in mol.smiles


def test_direct_parse_and_sdf_supplier_return_identical_diagnostics():
    direct_mol, _, _, direct_diag = chematic.from_mol_block_with_diagnostics(
        CONTRADICTORY_WEDGE_BLOCK
    )

    sdf = CONTRADICTORY_WEDGE_BLOCK + "$$$$\n"
    records = list(chematic.iter_sdf_str(sdf))
    assert len(records) == 1
    assert records[0].stereo_diagnostics() == direct_diag
    assert records[0].smiles == direct_mol.smiles
