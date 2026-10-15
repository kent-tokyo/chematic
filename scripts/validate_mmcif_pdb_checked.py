#!/usr/bin/env python3
"""Pinned Open Babel #3027 diagnostic plus independent PDB field checks."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import chematic

ROOT = Path(__file__).resolve().parents[1]
GLYCAN = ROOT / 'validation/openbabel_issues/mmcif/branched-residues.cif'
LONG = ROOT / 'validation/openbabel_issues/mmcif/long-component-name.cif'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def records(text):
    return [s for s in text.splitlines() if s.startswith(('ATOM  ', 'HETATM'))]


def convert(text, **kwargs):
    return chematic.mmcif_to_pdb_checked(text, allow_metadata_loss=True, **kwargs)


def field_checks(line, source):
    assert len(line) == 80
    assert line[:6].strip() == source[0]
    assert int(line[6:11]) == int(source[1])
    assert line[12:16].strip() == source[18]
    assert line[16] == (' ' if source[4] in ('.', '?') else source[4])
    assert line[17:20].strip() == source[16]
    assert line[21] == source[17]
    assert int(line[22:26]) == int(source[15])
    assert line[26] == (' ' if source[9] in ('.', '?') else source[9])
    assert [float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54))] == [float(v) for v in source[10:13]]
    assert float(line[54:60]) == float(source[13])
    assert float(line[60:66]) == float(source[14])
    assert line[76:78].strip() == source[2].upper()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--obabel', default='/opt/homebrew/bin/obabel')
    parser.add_argument('--output', type=Path, default=ROOT / 'validation/results/openbabel-mmcif-pdb-2026-10-15.json')
    args = parser.parse_args()
    version = subprocess.run([args.obabel, '-V'], capture_output=True, text=True, check=True, timeout=10).stdout.strip()
    if not version.startswith('Open Babel 3.2.1 '):
        raise RuntimeError(f'expected Open Babel 3.2.1, got {version}')
    text = GLYCAN.read_text()
    header = text.split('HETATM')[0]
    rows = [s.split() for s in text.splitlines() if s.startswith('HETATM')]
    checked = []
    for order in itertools.permutations(range(6)):
        source = header + '\n'.join(' '.join(rows[i]) for i in order) + '\n#\n'
        result = convert(source)
        output = records(result['pdb'])
        assert len(output) == 6
        assert result['source_text'] == source
        assert result['source'] == chematic.parse_mmcif_residues(source)
        for line, old in zip(output, order):
            field_checks(line, rows[old])
        assert len(set((line[21], line[22:26], line[26]) for line in output)) == 3
        checked.append({'order': order, 'output_sha256': digest(result['pdb'].encode())})
    long_text = LONG.read_text()
    # Include a numeric label-sequence control so residue merging (#3028)
    # cannot explain the independent five-character name truncation (#3027).
    long_rows = [s.split() for s in long_text.splitlines() if s.startswith('HETATM')]
    for r in long_rows: r[8] = r[15]
    long_control = long_text.split('HETATM')[0] + '\n'.join(' '.join(r) for r in long_rows) + '\n#\n'
    rejection = []
    for namespace in ('auth', 'label'):
        for allow in (False, True):
            try:
                chematic.mmcif_to_pdb_checked(long_control, namespace=namespace, allow_metadata_loss=allow)
            except ValueError as error:
                assert 'residue name' in str(error)
                rejection.append({'namespace': namespace, 'allow_metadata_loss': allow, 'error': str(error)})
            else:
                raise AssertionError('five-character residue unexpectedly accepted')
    with tempfile.TemporaryDirectory(prefix='chematic-pdb-checked-') as tmp:
        tmp = Path(tmp)
        fixture = tmp / 'long.cif'; fixture.write_text(long_control)
        process = subprocess.run([args.obabel, '-immcif', str(fixture), '-opdb'], capture_output=True, text=True, check=True, timeout=10)
        output = records(process.stdout)
        assert len(output) == 2
        assert [float(s[30:38]) for s in output] == [20, 22]
        assert [s[17:20].strip() for s in output] == ['A1A', 'A1A'], 'upstream behavior changed; inspect before rebaseline'
        control = tmp / 'short.cif'; control.write_text(long_control.replace('A1AJ1', 'UNL'))
        short = subprocess.run([args.obabel, '-immcif', str(control), '-opdb'], capture_output=True, text=True, check=True, timeout=10)
        assert [s[17:20].strip() for s in records(short.stdout)] == ['UNL', 'UNL']
        native_short = convert(control.read_text())
        assert [s[17:20].strip() for s in records(native_short['pdb'])] == ['UNL', 'UNL']
        pdb = tmp / 'native.pdb'; pdb.write_text(convert(text)['pdb'])
        cross = subprocess.run([args.obabel, '-ipdb', str(pdb), '-opdb'], capture_output=True, text=True, check=True, timeout=10)
        cross_rows = records(cross.stdout)
        assert len(cross_rows) == 6
        assert [s[17:20].strip() for s in cross_rows] == ['NAG'] * 4 + ['BMA'] * 2
        assert [int(s[22:26]) for s in cross_rows] == [1, 1, 2, 2, 3, 3]
        assert [float(s[30:38]) for s in cross_rows] == [20, 22, 40, 42, 60, 62]
    paths = ['crates/chematic-mol/src/mmcif_pdb.rs', 'crates/chematic-py/src/formats/mmcif_pdb.rs', 'crates/chematic-py/src/formats/mmcif_residues.rs', 'scripts/validate_mmcif_pdb_checked.py']
    record = {
        'schema_version': 1, 'status': 'local-source-candidate',
        'issue': 'https://github.com/openbabel/openbabel/issues/3027',
        'base': 'cac480fb source snapshot plus previously preserved RDKit/Open Babel candidates; patch isolates this adapter',
        'scope': 'model-1 ATOM/HETATM coordinate block plus END; no bonds, TER, full deposition export or arbitrary mmCIF round trip',
        'runtime': {'python': sys.version.split()[0], 'openbabel': version, 'chematic_extension_sha256': digest(Path(chematic.__file__).resolve().read_bytes())},
        'source_sha256': {s: digest((ROOT / s).read_bytes()) for s in paths},
        'fixture': {
            'path': str(LONG.relative_to(ROOT)), 'sha256': digest(long_text.encode()),
            'provenance': {'source': 'original synthetic atom_site rows derived from the project branched-residues fixture', 'date': '2026-10-15', 'transformations': 'retain first two rows; replace both component namespaces with A1AJ1; numeric-label control assigns source auth_seq_id to label_seq_id', 'external_data': 'no downloaded structure or upstream attachment'},
            'license': 'MIT OR Apache-2.0', 'redistribution': {'status': 'allowed', 'rationale': 'original project fixture'},
            'minimization': 'two isolated C/O coordinate records, one residue, five-character component ID; not a chemical model of A1AJ1',
            'oracle': {'name': 'wwPDB PDB format 3.3 fixed columns / source-authored identifiers', 'version': '3.3', 'url': 'https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html', 'settings': 'no renumbering, no truncation, auth or label selected explicitly'},
            'expected': {'native': 'InvalidField residue name in both namespaces regardless of metadata consent', 'external_Open_Babel_3_2_1': 'A1AJ1 truncates to A1A; UNL control preserved'},
        },
        'results': {'permutations': len(checked), 'field_mismatches': 0, 'source_sidecar_mismatches': 0, 'cases_sha256': digest(json.dumps(checked, sort_keys=True, separators=(',', ':')).encode()), 'rejections': rejection, 'native_short_pdb': native_short['pdb'], 'cross_read_preserved_residues': 3},
        'permutation_fixture': {'path': str(GLYCAN.relative_to(ROOT)), 'sha256': digest(text.encode()), 'expected': 'six rows, three NAG/NAG/BMA residues; all 720 input orders preserve auth names, residue keys, serials, coordinates and reader sidecar'},
        'openbabel': {'command': ['obabel', '-immcif', '<numeric-label-long-name-control>', '-opdb'], 'input_sha256': digest(long_control.encode()), 'output_names': [s[17:20].strip() for s in output], 'output_sha256': digest(process.stdout.encode()), 'stderr': process.stderr.strip(), 'short_name_control': {'input_sha256': digest(long_control.replace('A1AJ1', 'UNL').encode()), 'output_sha256': digest(short.stdout.encode()), 'passed': True}, 'native_PDB_cross_read': {'output_sha256': digest(cross.stdout.encode()), 'stderr': cross.stderr.strip(), 'passed': True}},
        'limitations': 'Inherited integer-only auth_seq_id/tokenizer; missing selected aliases fail. Parsed source sidecar covers the existing atom_site reader only; source_text retains the entire original Python input. losses lists known parsed-view omissions, not an exhaustive mmCIF category diff. Rust callers retain their original input text separately.',
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({'permutations': len(checked), 'field_mismatches': 0, 'upstream_truncation_reproduced': True, 'native_rejections': len(rejection), 'cross_read_preserved_residues': 3}))


if __name__ == '__main__':
    main()
