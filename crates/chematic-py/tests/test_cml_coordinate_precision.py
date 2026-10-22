"""Original binary64 writer checks; MIT OR Apache-2.0, no geometry generation."""
import struct
import xml.etree.ElementTree as ET

import chematic
import pytest


@pytest.mark.parametrize('value', [
    0.0, -0.0, 0.123456789012345, -1.23456789012345,
    1.23456789012345e-8, 1000000.123456789, 0.000049999999999,
    0.3333333333333333, 5e-324, 2.2250738585072014e-308,
    1.7976931348623157e308, 1.0000000000000002,
])
def test_to_cml_preserves_supplied_coordinate_values(value):
    mol = chematic.from_smiles('CO')
    supplied = [[value, -value], [1.523456789012345, -0.334567890123456]]
    text = mol.to_cml(supplied)
    atoms = list(ET.fromstring(text).iter('{http://www.xml-cml.org/schema}atom'))
    assert [a.attrib['elementType'] for a in atoms] == ['C', 'O']
    for atom, xy in zip(atoms, supplied):
        for axis, expected in zip(('x2', 'y2'), xy):
            assert struct.pack('>d', float(atom.attrib[axis])) == struct.pack('>d', expected)


def test_to_cml_without_coordinates_does_not_generate_positions():
    text = chematic.from_smiles('CO').to_cml()
    atoms = list(ET.fromstring(text).iter('{http://www.xml-cml.org/schema}atom'))
    assert atoms
    assert all('x2' not in a.attrib and 'y2' not in a.attrib for a in atoms)
