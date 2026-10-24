"""Authored two-atom graph observer; dative class is not numeric order 1.

Unique elements identify atoms in this development packet. CDXML is observed
with ElementTree; emitted SMILES with RDKit (never the chemical oracle).
No schema, physical-validity, stereo, aromaticity or general isomorphism claim.
"""
from __future__ import annotations
import xml.etree.ElementTree as ET

ORDER = {'1': 'single', '2': 'double', '3': 'triple', '4': 'quadruple', 'dative': 'dative'}


def compare(expected, atoms, bond):
    fields = dict(atoms=False, bond_class=False, endpoints=False,
                  dative_direction=False if expected['bond']['kind'] == 'dative' else None)
    wanted = expected['atoms']
    if len(atoms) != len(wanted) or len({a['number'] for a in atoms}) != len(atoms):
        return fields
    numbers = [a['number'] for a in wanted]
    if len(set(numbers)) != len(numbers) or set(numbers) != {a['number'] for a in atoms}:
        return fields
    positions = {i: numbers.index(a['number']) for i, a in enumerate(atoms)}
    fields['atoms'] = sorted(atoms, key=lambda a: numbers.index(a['number'])) == wanted
    if bond is None:
        return fields
    fields['bond_class'] = bond['kind'] == expected['bond']['kind']
    try:
        actual = [positions[bond[key]] for key in ('begin', 'end')]
    except KeyError:
        return fields
    wanted_ends = [expected['bond'][key] for key in ('begin', 'end')]
    fields['endpoints'] = sorted(actual) == sorted(wanted_ends)
    if fields['dative_direction'] is not None:
        fields['dative_direction'] = fields['bond_class'] and actual == wanted_ends
    return fields


def cdxml_fields(text):
    root = ET.fromstring(text)
    fragments = root.findall('.//fragment')
    if len(fragments) != 1:
        raise ValueError('one fragment required by this observer')
    fragment = fragments[0]
    nodes = fragment.findall('n')
    if len({node.attrib['id'] for node in nodes}) != len(nodes):
        raise ValueError('duplicate node id')
    ids = {node.attrib['id']: i for i, node in enumerate(nodes)}
    atoms = [dict(number=int(n.attrib.get('Element', '6')), isotope=int(n.attrib.get('Isotope', '0')),
                  charge=int(n.attrib.get('Charge', '0')), hydrogens=int(n.attrib['NumHydrogens'])) for n in nodes]
    bonds = fragment.findall('b')
    if len(bonds) != 1:
        raise ValueError('one bond required by this observer')
    a = bonds[0].attrib
    return atoms, dict(kind=ORDER[a.get('Order', '1').strip()], begin=ids[a['B']], end=ids[a['E']])


def smiles_fields(text):
    from rdkit import Chem
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if len(lines) != 1:
        raise ValueError('one SMILES record required')
    params = Chem.SmilesParserParams()
    params.sanitize = False
    params.removeHs = False
    mol = Chem.MolFromSmiles(lines[0], params)
    if mol is None:
        raise ValueError('invalid SMILES')
    mol.UpdatePropertyCache(strict=False)
    atoms = [dict(number=a.GetAtomicNum(), isotope=a.GetIsotope(), charge=a.GetFormalCharge(),
                  hydrogens=a.GetTotalNumHs(includeNeighbors=True)) for a in mol.GetAtoms()]
    if mol.GetNumBonds() != 1:
        raise ValueError('one bond required')
    b = mol.GetBondWithIdx(0)
    kinds = {'SINGLE': 'single', 'DOUBLE': 'double', 'TRIPLE': 'triple', 'QUADRUPLE': 'quadruple', 'DATIVE': 'dative'}
    return atoms, dict(kind=kinds[str(b.GetBondType())], begin=b.GetBeginAtomIdx(), end=b.GetEndAtomIdx())


def observe(case, result, format):
    output = result.get('output')
    if not isinstance(output, str) or not output.strip():
        outcome = 'error' if result.get('error') else 'process_error' if result.get('exit_code', 0) else 'no_output'
        return dict(outcome=outcome, preserved=False,
                    profile_refused=case['group'] == 'refusal' and 'unsupported CDXML bond order' in result.get('error', ''))
    if case['group'] == 'refusal':
        return dict(outcome='output_outside_profile', preserved=False, profile_refused=False)
    try:
        atoms, bond = (cdxml_fields if format == 'cdxml' else smiles_fields)(output)
        fields = compare(case['expected'], atoms, bond)
        return dict(outcome='graph_output', fields=fields, actual_atoms=atoms, actual_bond=bond,
                    preserved=all(value for value in fields.values() if value is not None))
    except (ValueError, KeyError, TypeError, ET.ParseError, RuntimeError) as error:
        return dict(outcome='invalid_output', preserved=False, observer_error=str(error))
