#!/usr/bin/env python3
"""Bounded native MOL2 regeneration diagnostic, with RDKit as a cross-reader.

Authored tokens supply the expected binary64 values. Open Babel and chematic
are judged by the same raw observer and RDKit parser. No speed/coverage ranking.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import struct
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run, identity, require, strict_json

WRITERS = ('record', 'record_adjacent', 'graph', 'source_graph')
PRODUCTION = ('crates/chematic-mol/src/mol2_tripos.rs', 'crates/chematic-mol/src/mol2_source.rs')


def bits(value):
    return struct.pack('>d', float(value)).hex()


def raw_snapshot(text):
    """Inspect ATOM real columns without either engine's chemical inference."""
    active = False
    rows = []
    for line in text.splitlines():
        if line.startswith('@<TRIPOS>'):
            active = line == '@<TRIPOS>ATOM'
        elif active and line.strip():
            row = line.split()
            require(len(row) >= 9, 'missing MOL2 numeric columns')
            rows.append(row)
    require(rows, 'missing MOL2 atoms')
    ids = {int(row[0]): i for i, row in enumerate(rows)}
    require(len(ids) == len(rows), 'duplicate MOL2 atom id')
    symbols = {'C': 6}
    active = False
    bonds = []
    for line in text.splitlines():
        if line.startswith('@<TRIPOS>'):
            active = line == '@<TRIPOS>BOND'
        elif active and line.strip():
            row = line.split()
            bonds.append([*sorted((ids[int(row[1])], ids[int(row[2])])), int(row[3])])
    return dict(coordinates=[[bits(v) for v in row[2:5]] for row in rows],
                partial_charges=[bits(row[8]) for row in rows],
                atomic_numbers=[symbols[row[5].split('.')[0]] for row in rows], bonds=sorted(bonds))


def rdkit_snapshot(text):
    from rdkit import Chem
    mol = Chem.MolFromMol2Block(text, sanitize=False, removeHs=False, cleanupSubstructures=False)
    require(mol is not None and mol.GetNumConformers() == 1, 'RDKit cannot read MOL2')
    conf = mol.GetConformer()
    return dict(coordinates=[[bits(v) for v in conf.GetAtomPosition(i)] for i in range(mol.GetNumAtoms())],
                partial_charges=[bits(a.GetProp('_TriposPartialCharge')) for a in mol.GetAtoms()],
                atomic_numbers=[a.GetAtomicNum() for a in mol.GetAtoms()],
                bonds=sorted([*sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx())), b.GetBondTypeAsDouble()]
                             for b in mol.GetBonds()))


def observe(expected, text, reader):
    try:
        actual = reader(text)
    except (ValueError, RuntimeError, KeyError, IndexError, OverflowError) as error:
        return dict(readable=False, coordinates_exact=False, partial_charges_exact=False, indexed_graph_exact=False, error=str(error))
    return dict(readable=True, indexed_graph_exact=all(actual[k] == expected[k] for k in ('atomic_numbers', 'bonds')),
                coordinates_exact=actual['coordinates'] == expected['coordinates'],
                partial_charges_exact=actual['partial_charges'] == expected['partial_charges'], actual=actual)


def expected_snapshot(case):
    return dict(coordinates=[[bits(v) for v in row] for row in case['expected_coordinates']],
                partial_charges=[bits(v) for v in case['expected_partial_charges']],
                atomic_numbers=case['expected_atomic_numbers'], bonds=case['expected_bonds'])


def terminal(run):
    return {key: value.decode('utf8', errors='replace') if isinstance(value, bytes) else value
            for key, value in run.items()}


def checked_run(command, data=b''):
    run = terminal(bounded_run(command, data))
    run['command'] = command
    require(run['resource_status'] == 'completed' and run['exit_code'] == 0, f'failed command: {command}')
    return run


def source_identities(revision):
    result = {}
    for path in PRODUCTION:
        data = subprocess.check_output(['git', 'show', f'{revision}:{path}'], cwd=ROOT) if revision else (ROOT / path).read_bytes()
        result[path] = dict(sha256=hashlib.sha256(data).hexdigest(), size_bytes=len(data))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native-executable', required=True, type=Path)
    parser.add_argument('--probe-source', type=Path, default=ROOT / 'crates/chematic-mol/examples/mol2_numeric_dump.rs')
    parser.add_argument('--source-revision', help='Pinned production source revision for a saved baseline executable')
    parser.add_argument('--obabel', type=Path, default=Path('/opt/homebrew/bin/obabel'))
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--baseline', action='store_true')
    args = parser.parse_args()
    from rdkit import rdBase
    from rdkit.Chem import rdmolfiles
    packet_path = ROOT / 'validation/mol2_numeric_preservation/cases.json'
    packet = strict_json(packet_path.read_text())
    native = checked_run([str(args.native_executable.resolve()), str(packet_path)])
    outputs = strict_json(native['stdout'])['rows']
    require(len(outputs) == len(packet['cases']), 'native row count mismatch')
    rows = []
    for case, output in zip(packet['cases'], outputs):
        require(case['id'] == output['id'], 'native case order mismatch')
        expected = expected_snapshot(case)
        require(raw_snapshot(case['input']) == expected, 'authored source/oracle disagreement')
        row = dict(id=case['id'], scope=case['scope'], input=case['input'], expected=expected,
                   native_output=output, native_raw={k: observe(expected, output[k], raw_snapshot) for k in WRITERS})
        if case['scope'] == 'paired':
            row['source_rdkit'] = observe(expected, case['input'], rdkit_snapshot)
            require(all(row['source_rdkit'][key] for key in ('readable', 'indexed_graph_exact', 'coordinates_exact', 'partial_charges_exact')), 'source not qualified by external reader')
            row['native_rdkit'] = {k: observe(expected, output[k], rdkit_snapshot) for k in WRITERS}
            command = [str(args.obabel.resolve()), '-imol2', '-omol2']
            converted = terminal(bounded_run(command, case['input'].encode()))
            converted['command'] = command
            row['openbabel_run'] = converted
            successful = converted['resource_status'] == 'completed' and converted['exit_code'] == 0
            text = converted['stdout'] if successful else ''
            row['openbabel_raw'] = observe(expected, text, raw_snapshot)
            row['openbabel_rdkit'] = observe(expected, text, rdkit_snapshot)
        rows.append(row)
    paired = [row for row in rows if row['scope'] == 'paired']
    metrics = ('readable', 'indexed_graph_exact', 'coordinates_exact', 'partial_charges_exact')
    summary = {k: {m: sum(r['native_rdkit'][k][m] for r in paired) for m in metrics} for k in WRITERS}
    summary['openbabel'] = {m: sum(r['openbabel_rdkit'][m] for r in paired) for m in metrics}
    for k in (*WRITERS, 'openbabel'):
        observations = [r['openbabel_rdkit'] if k == 'openbabel' else r['native_rdkit'][k] for r in paired]
        summary[k]['indexed_numeric_exact'] = sum(all(o[m] for m in metrics) for o in observations)
    result = dict(schema='chematic.mol2-numeric-preservation.v1', baseline=args.baseline,
                  head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  source_revision=args.source_revision or 'working_tree', production_sources=source_identities(args.source_revision),
                  packet=identity(packet_path), harness=identity(__file__), probe_source=identity(args.probe_source), native_executable=identity(args.native_executable),
                  native_run=native, openbabel=identity(args.obabel), openbabel_version=checked_run([str(args.obabel.resolve()), '-V']),
                  rdkit_version=rdBase.rdkitVersion, rdkit_parser=identity(rdmolfiles.__file__),
                  rdkit_options=dict(sanitize=False, removeHs=False, cleanupSubstructures=False),
                  paired_case_count=len(paired), native_only_case_count=len(rows)-len(paired), summary=summary, cases=rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(gzip.compress((json.dumps(result, indent=2)+'\n').encode(), mtime=0))
    print(json.dumps(summary))
    if not args.baseline:
        for row in rows:
            for k in WRITERS:
                require(row['native_raw'][k]['indexed_graph_exact'], 'native indexed graph drift')
                require(row['native_raw'][k]['coordinates_exact'], 'native coordinate drift')
                if k.startswith('record'):
                    require(row['native_raw'][k]['partial_charges_exact'], 'native partial charge drift')
            if row['scope'] == 'paired':
                for k in WRITERS:
                    require(row['native_rdkit'][k]['indexed_graph_exact'], 'external indexed graph drift')
                    require(row['native_rdkit'][k]['coordinates_exact'], 'external coordinate drift')
                    if k.startswith('record'):
                        require(row['native_rdkit'][k]['partial_charges_exact'], 'external partial charge drift')


if __name__ == '__main__':
    main()
