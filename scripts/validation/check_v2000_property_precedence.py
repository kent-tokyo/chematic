#!/usr/bin/env python3
"""Spec-authored V2000 property precedence, with output and input-reader lanes."""
from __future__ import annotations
import argparse
import gzip
import importlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run, identity, require, strict_json

NODE = """
const fs=require('node:fs'),w=require(process.argv[1]);
const p=JSON.parse(fs.readFileSync(0,'utf8'));
process.stdout.write(JSON.stringify(p.cases.map(c=>{
 const m=w.mol_from_sdf_block(c.input);
 try{return {id:c.id,outputs:{direct:m.canonical_smiles(),
  common:w.convert_common_format(c.input,'mol','smiles'),
  sdf:JSON.parse(w.sdf_to_smiles_json(c.input+'$$$$\\n'))[0]}};
 }finally{m.free();}
})));
"""


def terminal(run):
    return {k: v.decode('utf8', errors='replace') if isinstance(v, bytes) else v
            for k, v in run.items()}


def normalize(fields):
    """This packet has exactly one atom per atomic number, not general isomorphism."""
    numbers = fields['numbers']
    size = len(numbers)
    require(size > 0 and len(set(numbers)) == size, 'packet requires unique atomic numbers')
    for key in ['numbers', 'charges', 'isotopes', 'atom_maps']:
        require(len(fields[key]) == size and all(type(v) is int for v in fields[key]), 'invalid atom fields')
    atoms = sorted(zip(numbers, fields['charges'], fields['isotopes'], fields['atom_maps']))
    bonds = []
    for a, b, order in fields['bonds']:
        require(type(a) is int and type(b) is int and 0 <= a < size and 0 <= b < size and a != b,
                'invalid bond endpoints')
        require(type(order) in [int, float] and order in [1, 2, 3], 'unsupported observed bond order')
        bonds.append([*sorted([numbers[a], numbers[b]]), order])
    return dict(atoms=atoms, bonds=sorted(bonds))


def molecule_fields(mol):
    require(mol is not None, 'unreadable molecule')
    return dict(numbers=[a.GetAtomicNum() for a in mol.GetAtoms()],
                charges=[a.GetFormalCharge() for a in mol.GetAtoms()],
                isotopes=[a.GetIsotope() for a in mol.GetAtoms()],
                atom_maps=[a.GetAtomMapNum() for a in mol.GetAtoms()],
                bonds=[[b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble()] for b in mol.GetBonds()])


def read_smiles(text):
    from rdkit import Chem
    params = Chem.SmilesParserParams()
    params.sanitize = False
    params.removeHs = False
    return molecule_fields(Chem.MolFromSmiles(text, params))


def read_input(text):
    from rdkit import Chem
    return molecule_fields(Chem.MolFromMolBlock(text, sanitize=False, removeHs=False, strictParsing=True))


def observe(expected, text, reader=None):
    try:
        actual = (reader or read_smiles)(text)
        return dict(readable=True, fields_exact=normalize(expected) == normalize(actual), actual=actual)
    except (ValueError, RuntimeError, KeyError, IndexError, TypeError, OverflowError) as error:
        return dict(readable=False, fields_exact=False, error=str(error))


def production_source(revision):
    path = 'crates/chematic-mol/src/mol2000/reader.rs'
    data = subprocess.check_output(['git', 'show', f'{revision}:{path}'], cwd=ROOT) if revision else None
    return {path: identity(ROOT / path, data)}


def python_outputs(chematic, text, path):
    sdf = text + '$$$$\n'
    path.write_text(sdf)
    records = list(chematic.iter_sdf_str(sdf))
    fast = list(chematic.SDMolSupplier(path))
    forward = list(chematic.ForwardSDMolSupplier(io.BytesIO(sdf.encode())))
    require(len(records) == len(fast) == len(forward) == 1, 'SDF record accounting mismatch')
    return dict(direct=chematic.from_mol_block(text).smiles,
                coords=chematic.from_mol_block_with_coords(text)[0].smiles,
                common=chematic.convert_format(text, 'mol', 'smiles'),
                sdf=records[0].mol.smiles, fast=fast[0].smiles, forward=forward[0].smiles)


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
    path = ROOT / 'validation/v2000_property_precedence/cases.json'
    packet = strict_json(path.read_text())
    command = ['node', '-e', NODE, str(args.wasm_module.resolve())]
    node = terminal(bounded_run(command, json.dumps(packet).encode()))
    node['command'] = command
    require(node['resource_status'] == 'completed' and node['exit_code'] == 0, 'Node probe failed')
    wasm = strict_json(node['stdout'])
    require(len(wasm) == len(packet['cases']), 'Node row count mismatch')
    rows = []
    with tempfile.TemporaryDirectory(prefix='property-precedence-', dir=ROOT / 'target') as temp:
        for case, js in zip(packet['cases'], wasm):
            require(case['id'] == js['id'], 'Node order mismatch')
            expected = case['expected']
            normalize(expected)  # Enforce the authored packet's bounded correspondence rule.
            row = dict(id=case['id'], scope=case['scope'], input=case['input'], expected=expected,
                       python=python_outputs(chematic, case['input'], Path(temp) / 'input.sdf'), wasm=js['outputs'])
            for name in ['python', 'wasm']:
                row[name + '_observations'] = {key: observe(expected, value) for key, value in row[name].items()}
            command = ['/opt/homebrew/bin/obabel', '-imol', '-osmi', '-xn', '-xa']
            run = terminal(bounded_run(command, case['input'].encode()))
            run['command'] = command
            good = run['resource_status'] == 'completed' and run['exit_code'] == 0 and bool(run['stdout'].strip())
            row['openbabel'] = dict(run=run, observation=observe(expected, run['stdout'].strip()) if good else
                                   dict(readable=False, fields_exact=False))
            # Diagnostic only: source readers can disagree with the authored CTfile precedence oracle.
            row['rdkit_input'] = observe(expected, case['input'], read_input)
            rows.append(row)
    summary = {}
    for scope in ['control', 'precedence']:
        group = [r for r in rows if r['scope'] == scope]
        summary[scope] = dict(cases=len(group), **{
            name: sum(all(o['fields_exact'] for o in r[name + '_observations'].values()) for r in group)
            for name in ['python', 'wasm']},
            openbabel=sum(r['openbabel']['observation']['fields_exact'] for r in group),
            rdkit_input_diagnostic=sum(r['rdkit_input']['fields_exact'] for r in group))
    result = dict(schema='chematic.v2000-property-precedence-observation.v1', baseline=args.baseline,
                  head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  runtime_revision=args.runtime_revision or 'working_tree', production_sources=production_source(args.runtime_revision),
                  packet=identity(path), harness=identity(__file__), python_package=identity(chematic.__file__),
                  python_extension=identity(extension.__file__), wasm_js=identity(args.wasm_module),
                  wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),
                  rdkit_version=rdBase.rdkitVersion, rdkit_parser=identity(rdmolfiles.__file__),
                  rdkit_options=dict(sanitize=False, removeHs=False, mol_strictParsing=True),
                  openbabel=identity('/opt/homebrew/bin/obabel'),
                  openbabel_version=terminal(bounded_run(['/opt/homebrew/bin/obabel', '-V'])),
                  node_run=node, summary=summary, cases=rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(gzip.compress((json.dumps(result, indent=2) + '\n').encode(), mtime=0))
    print(json.dumps(summary))
    if not args.baseline:
        require(all(summary[s][name] == summary[s]['cases'] for s in summary for name in ['python', 'wasm']),
                'candidate property precedence failed')


if __name__ == '__main__':
    main()
