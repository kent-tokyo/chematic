"""Maestro writer boundaries and exact V2000 coordinate extraction."""
from dataclasses import FrozenInstanceError
import pytest
import chematic as ch


def water():
    atoms = [("O", (0., 0., .125)), ("H", (.96, 0., .125)), ("H", (-.24, .93, .5))]
    text = "water\n  authored         3D\ncomment\n  3  2  0  0  0  0  0  0  0  0999 V2000\n"
    for symbol, xyz in atoms:
        text += ''.join(f'{v:10.4f}' for v in xyz) + f' {symbol:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n'
    return text + "  1  2  1  0  0  0  0\n  1  3  1  0  0  0  0\nM  END\n"


def write(smiles="[H]O[H]", **kwargs):
    mol = ch.from_smiles(smiles)
    return ch.write_mae_checked(mol, [[i, i+.25, -i-.5] for i in range(mol.atom_count)],
                                allow_metadata_loss=True, **kwargs)


def test_graph_geometry_and_source_are_unchanged():
    mol = ch.from_smiles("[H]O[H]")
    mol.SetName('water "quoted" \\ tail')
    coords = [[0., .25, -.5], [1., 1.25, -1.5], [2., 2.25, -2.5]]
    snapshot = repr(coords), mol.to_moljson(), mol.GetPropsAsDict()
    result = ch.write_mae_checked(mol, coords, allow_metadata_loss=True)
    assert result.title == mol.GetName() and result.source_text is None
    assert result.atom_count == 3 and result.bond_count == 2
    assert ' m_atom[3] {' in result.text and ' m_bond[2] {' in result.text
    assert '1 0 0.25 -0.5 1 0' in result.text
    assert '\\"quoted\\"' in result.text and '\\\\ tail' in result.text
    assert snapshot == (repr(coords), mol.to_moljson(), mol.GetPropsAsDict())
    assert 'sd_properties' in result.omitted_features
    with pytest.raises(FrozenInstanceError):
        result.text = "changed"


@pytest.mark.parametrize("suffix", ["", ">  <Cat. No.>\nA-123\n\n$$$$\n"])
@pytest.mark.parametrize("crlf", [False, True])
def test_single_v2000_or_sd_record_keeps_source_and_z(suffix, crlf):
    source = water() + suffix
    if crlf:
        source = source.replace("\n", "\r\n")
    result = ch.v2000_to_mae_checked(source, allow_metadata_loss=True)
    assert result.source_text == source and result.title == "water"
    assert result.atom_count == 3 and result.bond_count == 2
    assert '1 0 0 0.125 8 0' in result.text
    assert '3 -0.23999999999999999 0.93000000000000005 0.5 1 0' in result.text


@pytest.mark.parametrize("smiles", ["[Na+].[Cl-]", "[H+]", "[H][H]", "[H]N([H])[H]",
    "[H]C([H])([H])[H]", "[H]C(=O)[H]", "[H]C#N", "[H]P([H])[H]", "[H]S[H]"])
def test_supported_explicit_h_and_ions(smiles):
    assert write(smiles).atom_count == ch.from_smiles(smiles).atom_count


@pytest.mark.parametrize("smiles", ["C", "[C]", "[H]c1ccccc1", "[2H]O[H]", "[H:1]O[H]",
    "*", "[He]", "[Fe+2]", "F/C=C/F", "F[C@](Cl)(Br)I"])
def test_unsupported_chemistry_refuses(smiles):
    with pytest.raises(ValueError):
        write(smiles)


@pytest.mark.parametrize("title,error", [(1, TypeError), (b"x", TypeError),
    ("bad\nname", ValueError), ("bad\rname", ValueError), ("bad\tname", ValueError),
    ("λ", ValueError), ("x"*257, ValueError)])
def test_bad_titles(title, error):
    with pytest.raises(error):
        write(title=title)


def test_empty_and_boundary_quoted_titles():
    assert write(title="").title == ""
    assert write(title="x"*256).title == "x"*256
    assert write(title='{} ::: # " \\').title == '{} ::: # " \\'


@pytest.mark.parametrize("coordinates,error", [(None, TypeError), ("abc", TypeError),
    ([], ValueError), ([[0,0]]*3, ValueError), ([[0,0,0,0]]*3, ValueError),
    ([[0,0,"0"]]*3, TypeError), ([[0,0,True]]*3, TypeError),
    ([[0,0,float('nan')]]*3, ValueError), ([[0,0,float('inf')]]*3, ValueError),
    ([[0,0,1e9]]*3, ValueError), ([[0,0,1e-101]]*3, ValueError)])
def test_coordinates(coordinates, error):
    with pytest.raises(error):
        ch.write_mae_checked(ch.from_smiles('[H]O[H]'), coordinates, allow_metadata_loss=True)


def test_consent_and_molecule_types():
    with pytest.raises(TypeError):
        ch.write_mae_checked(None, [], allow_metadata_loss=True)
    with pytest.raises(ValueError, match="consent"):
        ch.write_mae_checked(ch.from_smiles('[H]O[H]'), [])
    with pytest.raises(TypeError):
        ch.write_mae_checked(ch.from_smiles('[H]O[H]'), [], allow_metadata_loss=1)


@pytest.mark.parametrize("change", [lambda s: s.replace("V2000", "V3000"),
    lambda s: s.replace("M  END", "M  ISO  1   1  18\nM  END"),
    lambda s: s.replace("M  END", "M  RAD  1   1   1\nM  END"),
    lambda s: s.replace("M  END", "M  ALS  query\nM  END"),
    lambda s: s.replace("  1  2  1  0", "  1  2  1  1"),
    lambda s: s.replace("M  END", "M  STY  1   1 SUP\nM  END"),
    lambda s: s.replace("M  END", ""), lambda s: s + "arbitrary\n",
    lambda s: s + "$$$$\n" + s, lambda s: s + "$$$$\n$$$$\n",
    lambda s: s[:100], lambda s: s.replace("  3  2", "999999")])
def test_v2000_unsupported_or_multiple_records(change):
    with pytest.raises(ValueError):
        ch.v2000_to_mae_checked(change(water()), allow_metadata_loss=True)


def test_v2000_types_and_caps(monkeypatch):
    import chematic._maestro as implementation
    with pytest.raises(TypeError):
        ch.v2000_to_mae_checked(water().encode(), allow_metadata_loss=True)
    for name, value in [("_MAX_ATOMS", 2), ("_MAX_BONDS", 1), ("_MAX_OUTPUT", 8)]:
        with monkeypatch.context() as local:
            local.setattr(implementation, name, value)
            with pytest.raises(ValueError):
                write()
    monkeypatch.setattr(implementation, "_MAX_INPUT", 8)
    with pytest.raises(ValueError):
        ch.v2000_to_mae_checked(water(), allow_metadata_loss=True)


@pytest.mark.parametrize("column,value", [(34," 1"), (39,"  1"), (42,"  1"),
                                         (45,"  1"), (48,"  2"), (60,"  1")])
def test_raw_atom_flags_are_never_silently_projected(column, value):
    rows = water().splitlines()
    rows[4] = rows[4][:column] + value + rows[4][column+len(value):]
    with pytest.raises(ValueError):
        ch.v2000_to_mae_checked('\n'.join(rows)+'\n', allow_metadata_loss=True)


def test_raw_duplicate_bond_refuses():
    source = water().replace('  3  2', '  3  3').replace('M  END', '  1  2  1  0  0  0  0\nM  END')
    with pytest.raises(ValueError):
        ch.v2000_to_mae_checked(source, allow_metadata_loss=True)


def test_v2000_charge_records_survive_without_bonds():
    rows = water().splitlines()
    rows[3] = rows[3].replace('  3  2', '  2  0')
    a = rows[4][:31] + 'Na ' + rows[4][34:]
    b = rows[5][:31] + 'Cl ' + rows[5][34:]
    source = '\n'.join(rows[:4]+[a,b,'M  CHG  2   1   1   2  -1','M  END'])+'\n'
    result = ch.v2000_to_mae_checked(source,allow_metadata_loss=True)
    assert result.atom_count == 2 and result.bond_count == 0
    assert '11 1\n' in result.text and '17 -1\n' in result.text


@pytest.mark.parametrize("token", ["1_0", "nan", "1e-999999", "1e-101"])
def test_raw_coordinates_refuse_invalid_or_lost_values(token):
    rows = water().splitlines()
    rows[4] = f'{token:>10}' + rows[4][10:]
    with pytest.raises(ValueError):
        ch.v2000_to_mae_checked('\n'.join(rows)+'\n',allow_metadata_loss=True)


@pytest.mark.parametrize("record", ["M  CHG  1   9   0", "M  CHG  1   1 257",
    "M  CHG  2   1   0   1   0", "M  CHG  0", "M  CHG  2   1   0",
    "M CHG  1   1   0", "M  CHG  1 0_1   0", "M  CHG  1   1 0_0"])
def test_malformed_charge_data_is_not_ignored(record):
    with pytest.raises(ValueError):
        ch.v2000_to_mae_checked(water().replace('M  END',record+'\nM  END'),allow_metadata_loss=True)
