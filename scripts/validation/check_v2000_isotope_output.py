#!/usr/bin/env python3
"""Bounded SMILES -> regenerated V2000 comparison; isotope labels are explicit."""
from __future__ import annotations

import argparse
import gzip
import importlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run, identity, require, strict_json

NODE = """
const fs=require('node:fs'),w=require(process.argv[1]);
const p=JSON.parse(fs.readFileSync(0,'utf8'));
process.stdout.write(JSON.stringify(p.cases.map(c=>{
 const m=w.parse_smiles(c.smiles);
 try{return {id:c.id,outputs:{direct:w.to_mol_block(m),
  common_mol:w.convert_common_format(c.smiles,'smiles','mol'),
  common_sdf:w.convert_common_format(c.smiles,'smiles','sdf'),
  array_sdf:w.smiles_array_to_sdf(JSON.stringify([c.smiles]))}};
 }finally{m.free();}
})));
"""


def terminal(run):
    return {k: v.decode('utf8', errors='replace') if isinstance(v, bytes) else v
            for k, v in run.items()}


def read_fields(text):
    from rdkit import Chem
    mol = Chem.MolFromMolBlock(text, sanitize=False, removeHs=False, strictParsing=True)
    require(mol is not None, 'unreadable V2000 output')
    return dict(numbers=[a.GetAtomicNum() for a in mol.GetAtoms()],
                charges=[a.GetFormalCharge() for a in mol.GetAtoms()],
                isotopes=[a.GetIsotope() for a in mol.GetAtoms()],
                atom_maps=[a.GetAtomMapNum() for a in mol.GetAtoms()],
                bonds=sorted([*sorted([b.GetBeginAtomIdx(), b.GetEndAtomIdx()]),
                              b.GetBondTypeAsDouble()] for b in mol.GetBonds()))


def isotope_layout(text):
    """Check only atom mass-difference bounds and fixed-width M ISO records."""
    try:
        lines = text.splitlines()
        require('V2000' in lines[3], 'not V2000')
        atoms = int(lines[3][:3])
        require(len(lines) >= 4 + atoms, 'missing atom lines')
        for line in lines[4:4 + atoms]:
            require(len(line) >= 36 and -3 <= int(line[34:36]) <= 4,
                    'legacy mass difference outside -3..4')
        seen = set()
        for line in lines[4 + atoms:]:
            if not line.startswith('M  ISO'):
                continue
            count = int(line[6:9])
            require(1 <= count <= 8 and len(line.rstrip()) == 9 + 8 * count,
                    'bad M ISO count/width')
            for i in range(count):
                start = 9 + 8 * i
                idx, mass = int(line[start:start + 4]), int(line[start + 4:start + 8])
                require(1 <= idx <= atoms and 1 <= mass <= 9999 and idx not in seen,
                        'invalid or duplicate M ISO assignment')
                seen.add(idx)
        return True
    except (ValueError, IndexError, TypeError):
        return False


def observe(expected, text, reader=None):
    result = dict(isotope_layout=isotope_layout(text), readable=False, fields_exact=False)
    try:
        actual = (reader or read_fields)(text)
        result.update(readable=True, actual=actual, fields_exact=actual == expected)
    except (ValueError, RuntimeError, KeyError, IndexError, OverflowError) as error:
        result['error'] = str(error)
    return result


def passed(observation):
    return observation['fields_exact'] and observation['isotope_layout']


def production_sources(revision):
    paths = ['crates/chematic-mol/src/mol2000.rs',
             'crates/chematic-mol/src/mol2000/reader.rs',
             'crates/chematic-mol/src/mol2000/writer.rs',
             'crates/chematic-mol/src/mol2000/writer/isotope_output.rs',
             'crates/chematic-py/src/mol_methods/three_d.rs']
    present = set(subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', revision,
                                          '--', *paths], cwd=ROOT, text=True).splitlines()) if revision else set(paths)
    return {p: identity(ROOT / p, subprocess.check_output(['git', 'show', f'{revision}:{p}'], cwd=ROOT)
                        if revision else None) if p in present else dict(absent=True) for p in paths}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wasm-module', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--baseline', action='store_true')
    parser.add_argument('--runtime-revision')
    args = parser.parse_args()
    import chematic
    from rdkit import rdBase
    from rdkit.Chem import rdmolfiles
    extension = importlib.import_module('chematic.chematic')
    path = ROOT / 'validation/v2000_isotope_output/cases.json'
    packet = strict_json(path.read_text())
    command = ['node', '-e', NODE, str(args.wasm_module.resolve())]
    run = terminal(bounded_run(command, json.dumps(packet).encode()))
    run['command'] = command
    require(run['resource_status'] == 'completed' and run['exit_code'] == 0, 'Node probe failed')
    wasm = strict_json(run['stdout'])
    require(len(wasm) == len(packet['cases']), 'Node row count mismatch')
    rows = []
    for case, js in zip(packet['cases'], wasm):
        require(case['id'] == js['id'], 'Node row order mismatch')
        smiles, expected = case['smiles'], case['expected']
        mol = chematic.from_smiles(smiles)
        coords = [[float(i), 0.0] for i in range(len(expected['numbers']))]
        py = dict(direct=mol.to_mol_block(), coords=mol.to_mol_block_2d(coords),
                  common_mol=chematic.convert_format(smiles, 'smiles', 'mol'),
                  common_sdf=chematic.convert_format(smiles, 'smiles', 'sdf'),
                  charges_sdf=mol.to_sdf_with_charges([xy + [0.5] for xy in coords], [0.0] * len(coords)))
        row = dict(id=case['id'], scope=case['scope'], smiles=smiles, expected=expected,
                   python=py, wasm=js['outputs'])
        for name in ['python', 'wasm']:
            row[name + '_observations'] = {k: observe(expected, text) for k, text in row[name].items()}
        row['openbabel'] = {}
        for profile, flags in [('default', []), ('atom_classes', ['-xa'])]:
            command = ['/opt/homebrew/bin/obabel', '-ismi', '-omol', *flags]
            converted = terminal(bounded_run(command, (smiles + '\n').encode()))
            converted['command'] = command
            good = converted['resource_status'] == 'completed' and converted['exit_code'] == 0 and bool(converted['stdout'].strip())
            row['openbabel'][profile] = dict(run=converted, observation=observe(expected, converted['stdout']) if good else
                                            dict(readable=False, fields_exact=False, isotope_layout=False))
        rows.append(row)
    summary = {}
    for scope in ['ordinary', 'storage_boundary']:
        group = [r for r in rows if r['scope'] == scope]
        summary[scope] = dict(cases=len(group), **{
            name: sum(all(passed(o) for o in r[name + '_observations'].values()) for r in group)
            for name in ['python', 'wasm']}, **{
            'openbabel_' + name: sum(passed(r['openbabel'][name]['observation']) for r in group)
            for name in ['default', 'atom_classes']})
    result = dict(schema='chematic.v2000-isotope-observation.v1', baseline=args.baseline,
                  head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  runtime_revision=args.runtime_revision or 'working_tree', production_sources=production_sources(args.runtime_revision),
                  packet=identity(path), harness=identity(__file__), python_package=identity(chematic.__file__),
                  python_extension=identity(extension.__file__), wasm_js=identity(args.wasm_module),
                  wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),
                  rdkit_version=rdBase.rdkitVersion, rdkit_parser=identity(rdmolfiles.__file__),
                  rdkit_options=dict(sanitize=False, removeHs=False, strictParsing=True),
                  openbabel=identity('/opt/homebrew/bin/obabel'),
                  openbabel_version=terminal(bounded_run(['/opt/homebrew/bin/obabel', '-V'])),
                  node_run=run, summary=summary, cases=rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(gzip.compress((json.dumps(result, indent=2) + '\n').encode(), mtime=0))
    print(json.dumps(summary))
    if not args.baseline:
        require(all(summary[s][name] == summary[s]['cases'] for s in summary for name in ['python', 'wasm']),
                'candidate isotope or neighboring field preservation failed')


if __name__ == '__main__':
    main()
