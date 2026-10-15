#!/usr/bin/env python3
"""Original PDB subsets, pinned RDKit cross-read and Open Babel #2481 control."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import platform
import subprocess
import tempfile

import chematic
from rdkit import Chem, rdBase, RDLogger

ROOT = Path(__file__).resolve().parents[1]
FOLDER = ROOT/'validation/openbabel_issues/pdb-selection'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def atom_lines(text):
    return [s for s in text.splitlines() if s.startswith(('ATOM  ', 'HETATM'))]


def conect(text):
    return [tuple(int(s[i:i+5]) for i in range(6, 31, 5) if s[i:i+5].strip()) for s in text.splitlines() if s.startswith('CONECT')]


def recipe(original, order, reverse, explicit_model, crlf):
    rows = original.splitlines()
    groups = [rows[i:i+5] for i in (0, 5, 10)]
    mapping = {}
    coordinates = []
    for position, index in enumerate(order):
        pairs = [(groups[index][0], groups[index][1]), (groups[index][2], groups[index][3])]
        if reverse[index]:
            pairs.reverse()
        for j, (atom, anisou) in enumerate(pairs):
            serial = position*4+j+1
            mapping[int(atom[6:11])] = serial
            coordinates.extend([atom[:6]+f'{serial:5d}'+atom[11:], anisou[:6]+f'{serial:5d}'+anisou[11:]])
        coordinates.append(f'TER   {serial+1:5d}      '+atom[17:27])
    connections = ['CONECT'+''.join(f'{mapping[v]:5d}' for v in values) for values in conect(original)]
    if explicit_model:
        coordinates = ['MODEL        7', *coordinates, 'ENDMDL']
    text = '\n'.join([*coordinates, *connections, 'END'])+'\n'
    return text.replace('\n', '\r\n') if crlf else text


def rdkit_check(output, expected):
    mol = Chem.MolFromPDBBlock(output, sanitize=False, removeHs=False, proximityBonding=False)
    if not expected:
        assert mol is None or mol.GetNumAtoms() == 0
        return
    assert mol is not None and mol.GetNumAtoms() == len(expected)
    for atom, line in zip(mol.GetAtoms(), expected):
        info = atom.GetPDBResidueInfo()
        assert info is not None
        assert (info.GetSerialNumber(), info.GetName(), info.GetAltLoc(), info.GetResidueName(), info.GetChainId(), info.GetResidueNumber(), info.GetInsertionCode(), info.GetIsHeteroAtom()) == (
            int(line[6:11]), line[12:16], line[16], line[17:20], line[21], int(line[22:26]), line[26], line.startswith('HETATM'))
        point = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        assert [point.x, point.y, point.z] == [float(line[i:i+8]) for i in (30, 38, 46)]


def native_check(text, indices):
    original = atom_lines(text)
    selected = sorted(indices)
    expected = [original[i] for i in selected]
    r = chematic.select_pdb_atoms(text, list(reversed(selected)))
    assert r['source_text'] == text and r['source_atom_count'] == 6
    assert r['selected_source_indices'] == selected and r['losses'] == []
    assert r['old_to_new'] == [selected.index(i) if i in selected else None for i in range(6)]
    assert atom_lines(r['pdb']) == expected
    serials = {int(s[6:11]) for s in expected}
    expected_anisou = [s for s in text.splitlines() if s.startswith('ANISOU') and int(s[6:11]) in serials]
    assert [s for s in r['pdb'].splitlines() if s.startswith('ANISOU')] == expected_anisou
    links = []
    for values in conect(text):
        if values[0] in serials:
            targets = [v for v in values[1:] if v in serials]
            if targets:
                links.append((values[0], *targets))
    assert conect(r['pdb']) == links
    expected_ter = []
    last = None
    for line in text.splitlines():
        if line.startswith(('ATOM  ', 'HETATM')) and int(line[6:11]) in serials:
            last = line
        elif line.startswith('TER'):
            if last is not None:
                expected_ter.append((int(last[6:11])+1, last[17:27]))
            last = None
    assert [(int(s[6:11]), s[17:27]) for s in r['pdb'].splitlines() if s.startswith('TER')] == expected_ter
    assert chematic.select_pdb_atoms(r['pdb'])['pdb'] == r['pdb']
    rdkit_check(r['pdb'], expected)
    return r


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe', type=Path, default=Path('/private/tmp/chematic-ob2481-delete-probe'))
    parser.add_argument('--output', type=Path, default=ROOT/'validation/results/openbabel-pdb-selection-2026-10-15.json')
    args = parser.parse_args()
    assert rdBase.rdkitVersion == '2026.09.1', rdBase.rdkitVersion
    RDLogger.DisableLog('rdApp.*')
    manifest = json.loads((FOLDER/'fixtures.json').read_text())
    for record in manifest['records']:
        text = (FOLDER/record['file']).read_text()
        assert digest(text.encode()) == record['sha256']
        native_check(text, range(4))
    checksum = hashlib.sha256()
    cases = atoms = anisou = ter_updates = filtered = rdkit_nonempty = 0
    records = []
    original = (FOLDER/'rich.pdb').read_text()
    for order, reverse, model, crlf in itertools.product(itertools.permutations(range(3)), itertools.product((False, True), repeat=3), (False, True), (False, True)):
        text = recipe(original, order, reverse, model, crlf)
        for mask in range(64):
            indices = [i for i in range(6) if mask & (1 << i)]
            r = native_check(text, indices)
            checksum.update(r['pdb'].encode())
            cases += 1
            atoms += len(r['atoms'])
            anisou += len([s for s in r['pdb'].splitlines() if s.startswith('ANISOU')])
            ter_updates += len(r['retargeted_ter_lines'])
            filtered += len(r['filtered_conect_lines'])
            rdkit_nonempty += bool(indices)
        records.append(dict(id=f'g{"".join(map(str, order))}-r{"".join(str(int(v)) for v in reverse)}-m{int(model)}-crlf{int(crlf)}',
                            group_order=order, reversed_atom_pairs=reverse, explicit_model=model, crlf=crlf,
                            provenance='Original rich.pdb fixture transformed by scripts/validate_pdb_selection.py, authored 2026-10-15',
                            license='MIT OR Apache-2.0', redistribution='allowed; project-authored data',
                            minimization='Six atoms cover two residue fragments and one ligand, associated records and identity fields',
                            oracle='Exact raw-record/selection contract; RDKit 2026.09.1 PDB read with sanitize=False, removeHs=False, proximityBonding=False',
                            expected=dict(status='all 64 subsets preserve declared atom identity and associated-record contracts', units='Angstrom; source row indices', tolerance=0),
                            sha256=digest(text.replace('\r\n', '\n').encode()), exact_subsets=64, status='pass'))
    with tempfile.TemporaryDirectory(prefix='chematic-ob2481-') as directory:
        target = Path(directory)/'output.pdb'
        controls = []
        for deleted in ('', '4,5'):
            p = subprocess.run([str(args.probe), str(FOLDER/'basic.pdb'), deleted, str(target)], capture_output=True, text=True, check=True, timeout=10)
            data = json.loads(p.stdout)
            assert data['version'] == '3.2.1'
            before = {a['id']: a for a in data['before']}
            changes = [a['id'] for a in data['after'] if a != before[a['id']]]
            if not deleted:
                assert not changes and data['before'] == data['after_write']
            else:
                assert changes == [0, 1, 2, 3]
                assert [a['residue_name'] for a in data['after']] == ['UNL']*4
                assert [a['chain'] for a in data['after']] == [' ']*4
                assert data['after'] == data['after_write']
            controls.append(dict(deleted_source_indices=[] if not deleted else [4, 5], changed_retained_ids=changes, snapshots=data, output_sha256=digest(target.read_bytes())))
    report = dict(schema_version=1, date='2026-10-15', issue=dict(number=2481, url='https://github.com/openbabel/openbabel/issues/2481', state_verified='OPEN', scope='Original synthetic reproduction of metadata mutation after ligand deletion; original issue attachment not used'),
                  runtime=dict(python=platform.python_version(), rdkit=rdBase.rdkitVersion, openbabel='3.2.1'),
                  fixture_manifest_sha256=digest((FOLDER/'fixtures.json').read_bytes()),
                  oracle_settings=dict(native='All 64 atom subsets; 6 residue orders x 8 within-residue orders x LF/CRLF x implicit/explicit single model; exact atom/ANISOU contents and CONECT multiplicity; TER endpoint/serial; idempotence', rdkit='Chem.MolFromPDBBlock(sanitize=False, removeHs=False, proximityBonding=False), exact identities and coordinates', openbabel='Public OBMol.DeleteAtom in descending original index order; OBAtom.GetResidue before/after deletion and after PDB serialization; no-deletion control'),
                  provenance='Original synthetic MIT OR Apache-2.0 fixtures, derived only by atom/record order, serial remapping, model wrappers and line endings; no original issue attachment/upstream code copied',
                  results=dict(source_variants=len(records), exact_subsets=cases, selected_atom_lines=atoms, retained_anisou_lines=anisou, retargeted_ter_records=ter_updates, filtered_conect_records=filtered, rdkit_nonempty_cross_reads=rdkit_nonempty, rdkit_empty_checks=cases-rdkit_nonempty, unexpected_mismatches=0, output_sha256=checksum.hexdigest(), openbabel_no_deletion_changed_atoms=0, openbabel_ligand_deletion_changed_retained_atoms=4),
                  records=records, openbabel_controls=controls,
                  scope=['Single-model coordinate selection; no chemical or residue perception', 'No atom renumbering/reordering; TER serials/residue endpoints may update with collision checks', 'Unknown record loss needs explicit consent and original input is retained', 'Open Babel control reproduces one synthetic workflow, not all structures from #2481', 'Not a complete PDB entry validator or upstream fix/publication'])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['results'], indent=2))


if __name__ == '__main__':
    main()
