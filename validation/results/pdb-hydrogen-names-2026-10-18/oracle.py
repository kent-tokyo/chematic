"""Independent PDB naming recipe, exact byte checks and RDKit cross-reading."""
import hashlib
import itertools
import json
from pathlib import Path

import rdkit
from rdkit import Chem
import chematic as ch

ROOT = Path(__file__).resolve().parent
report = dict(rdkit=rdkit.__version__, cases=0, atom_cross_reads=0,
              selected_hydrogens=0, changed_anisou=0, exact_lines=0,
              idempotence_checks=0, default_selection_checks=0,
              implicit_model_bond_checks=0, explicit_model_coordinate_checks=0)


def source(order, crlf, model):
    rows, meta = [], []
    for logical in order:
        group, local = divmod(logical, 4)
        serial = (logical + 1) * 11
        name = ['H1  ', ' H  ', ' H2 ', ' H  '][local]
        element = 'C' if local == 0 else 'H'
        kind = 'HETATM' if local in (0, 2) else 'ATOM  '
        chain, residue, insertion = ('A', -2, 'I') if group == 0 else (' ', 7, ' ')
        line = (f'{kind}{serial:5d} {name} LIG {chain}{residue:4d}{insertion}   '
                f'{logical * .137:8.3f}{group * 2.:8.3f}{local * .231:8.3f}'
                f'{.75:6.2f}{23.45:6.2f}          {element:>2}  ')
        rows.append(line)
        rows.append('ANISOU' + line[6:28] + ''.join(f'{n:7d}' for n in (101, 202, 303, -11, 22, -33)) + line[70:])
        meta.append(dict(group=group, local=local, serial=serial, name=name))
    if model:
        rows = ['MODEL        7'] + rows + ['ENDMDL']
    # Explicit, symmetric ethane connectivity; no distance-based bond inference.
    for logical in range(8):
        group, local = divmod(logical, 4)
        targets = [group * 4 + i for i in (1, 2, 3)] + [(1-group)*4] if local == 0 else [group*4]
        rows.append('CONECT' + f'{(logical+1)*11:5d}' + ''.join(f'{(i+1)*11:5d}' for i in targets))
    rows.append('END')
    sep = '\r\n' if crlf else '\n'
    return sep.join(rows) + sep, meta


def read(text):
    mol = Chem.MolFromPDBBlock(text, sanitize=False, removeHs=False, proximityBonding=False)
    assert mol is not None and mol.GetNumAtoms() == 8
    atoms = []
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        atoms.append((atom.GetAtomicNum(), info.GetSerialNumber(), info.GetName(),
                      info.GetResidueName(), info.GetChainId(), info.GetResidueNumber(),
                      info.GetInsertionCode(), info.GetAltLoc(), info.GetIsHeteroAtom(),
                      info.GetOccupancy(), info.GetTempFactor(),
                      tuple(mol.GetConformer().GetAtomPosition(atom.GetIdx()))))
    bonds = sorted((min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                    max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), str(b.GetBondType()))
                   for b in mol.GetBonds())
    return atoms, bonds


permutations = [(0, 1, 2, 3), (3, 2, 1, 0), (1, 3, 0, 2), (2, 0, 3, 1)]
for permutation, groups, crlf, model in itertools.product(permutations, [(0, 1), (1, 0)], [False, True], [False, True]):
    order = [g * 4 + local for g in groups for local in permutation]
    text, meta = source(order, crlf, model)
    baseline_atoms, baseline_bonds = read(text)
    assert len(baseline_bonds) == (0 if model else 7)
    hydrogen_indices = [i for i, row in enumerate(meta) if row['local'] != 0]
    for mask in range(64):
        selected = [index for k, index in enumerate(hydrogen_indices) if mask & (1 << k)]
        # Independent allocation: construct the complete available name list per group.
        reserved = {g: {row['name'].strip() for i, row in enumerate(meta)
                        if row['group'] == g and i not in selected} for g in (0, 1)}
        available = {g: iter([f'H{n}' for n in range(1, 1000) if f'H{n}' not in reserved[g]]) for g in (0, 1)}
        expected = {}
        for index in selected:
            name = next(available[meta[index]['group']])
            expected[index] = (' ' + name).ljust(4) if len(name) < 4 else name
        result = ch.name_pdb_hydrogens_checked(text, selected[::-1])
        assert result.source_text == text and result.atom_count == 8
        assert [item.source_index for item in result.assignments] == selected
        assert {item.source_index: item.new_name for item in result.assignments} == expected
        assert len(result.pdb) == len(text)
        names_by_serial = {meta[i]['serial']: value for i, value in expected.items()}
        expected_updated = []
        for number, (a, b) in enumerate(zip(text.splitlines(keepends=True), result.pdb.splitlines(keepends=True)), 1):
            if a[:6].strip() in ('ATOM', 'HETATM', 'ANISOU') and int(a[6:11]) in names_by_serial:
                name = names_by_serial[int(a[6:11])]
                assert b == a[:12] + name + a[16:]
                if a[:6] == 'ANISOU' and a[12:16] != name:
                    expected_updated.append(number)
            else:
                assert a == b
            report['exact_lines'] += 1
        assert result.updated_anisou_lines == tuple(expected_updated)
        assert result.renamed_count == sum(meta[i]['name'] != name for i, name in expected.items())
        atoms, bonds = read(result.pdb)
        assert bonds == baseline_bonds
        for index, (a, b) in enumerate(zip(baseline_atoms, atoms)):
            assert b == a[:2] + (expected.get(index, a[2]),) + a[3:]
        assert ch.name_pdb_hydrogens_checked(result.pdb, selected).pdb == result.pdb
        report['idempotence_checks'] += 1
        if mask == 63:
            assert ch.name_pdb_hydrogens_checked(text).pdb == result.pdb
            report['default_selection_checks'] += 1
        report['cases'] += 1
        if model:
            report['explicit_model_coordinate_checks'] += 1
        else:
            report['implicit_model_bond_checks'] += 1
        report['atom_cross_reads'] += len(atoms)
        report['selected_hydrogens'] += len(selected)
        report['changed_anisou'] += len(expected_updated)

# Live Open Babel source is separately projected only with explicit metadata consent.
original = (ROOT/'openbabel-generated.pdb').read_text()
projection = ch.select_pdb_atoms(original, allow_metadata_loss=True)
result = ch.name_pdb_hydrogens_checked(projection['pdb'])
assert result.renamed_count == 4
assert [a.new_name.strip() for a in result.assignments] == ['H1', 'H2', 'H3', 'H4']
for before, after in zip(projection['pdb'].splitlines(keepends=True), result.pdb.splitlines(keepends=True)):
    assert before[:12] == after[:12] and before[16:] == after[16:]
(ROOT/'openbabel-coordinate-source.pdb').write_text(projection['pdb'])
(ROOT/'openbabel-named.pdb').write_text(result.pdb)
report['openbabel_omitted_metadata'] = projection['losses']
report['openbabel_renamed_hydrogens'] = result.renamed_count
report['rdkit_explicit_model_limitation'] = (
    'RDKit ignores global CONECT after ENDMDL in these controls (zero source/output bonds); '
    'bond cross-read applies to implicit models only. All original CONECT bytes '
    'and native checked references are verified for every case.'
)
report['oracle_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(ROOT/'oracle.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, sort_keys=True))
