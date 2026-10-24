#!/usr/bin/env python3
"""Local MolJSON H-authority diagnostic; RDKit observes the resulting SMILES.

The packet supplies expected non-node H counts and achiral reference SMILES.
Open Babel's installed format inventory is recorded, not assigned a score for
this unsupported dialect. This is not an independent MolJSON cross-reader.
"""
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
const fs=require('node:fs'), w=require(process.argv[1]);
const packet=JSON.parse(fs.readFileSync(0,'utf8'));
const rows=packet.cases.map(c=>{
 const m=w.mol_from_moljson(c.input);
 try{return {id:c.id,smiles:m.canonical_smiles(),output:w.to_moljson(m),
 common_smiles:w.convert_common_format(c.input,'moljson','smiles')};}finally{m.free();}
});
const refusals=packet.invalid_hydrogens.map(h=>{
 try{const m=w.mol_from_moljson(JSON.stringify({atoms:[{id:'a1',element:'C',hydrogens:h}],bonds:[]}));m.free();return false;}
 catch(e){return String(e).includes('hydrogens');}
});
process.stdout.write(JSON.stringify({rows,refusals}));
"""


def external_snapshot(smiles):
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    params = Chem.SmilesParserParams()
    params.removeHs = False
    mol = Chem.MolFromSmiles(smiles, params)
    require(mol is not None, 'RDKit could not interpret emitted SMILES')
    return dict(smiles=Chem.MolToSmiles(mol, canonical=True), formula=rdMolDescriptors.CalcMolFormula(mol),
                atoms=mol.GetNumAtoms(), radical_electrons=sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms()))


def observe(case, row, reader=external_snapshot):
    encoded = strict_json(row['output'])
    counts = [atom.get('hydrogens') for atom in encoded['atoms']]
    result = dict(non_node_h_exact=all(type(count) is int for count in counts) and counts == case['expected_hydrogens'],
                  smiles_graph_exact=False, common_smiles_graph_exact=False)
    expected = reader(case['expected_smiles'])
    for field, metric in [('smiles', 'smiles_graph_exact'), ('common_smiles', 'common_smiles_graph_exact')]:
        try:
            actual = reader(row[field])
            result[metric] = actual == expected
            result[field + '_observation'] = actual
        except (ValueError, RuntimeError) as error:
            result[field + '_error'] = str(error)
    result['preserved'] = all(result[k] for k in ('non_node_h_exact', 'smiles_graph_exact', 'common_smiles_graph_exact'))
    return result


def terminal(run):
    return {k: v.decode('utf8', errors='replace') if isinstance(v, bytes) else v for k, v in run.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wasm-module', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--baseline-runtime-revision')
    parser.add_argument('--baseline', action='store_true')
    args = parser.parse_args()
    import chematic
    from rdkit import rdBase
    from rdkit.Chem import rdmolfiles
    extension = importlib.import_module('chematic.chematic')
    packet_path = ROOT / 'validation/moljson_hydrogens/cases.json'
    packet = strict_json(packet_path.read_text())
    command = ['node', '-e', NODE, str(args.wasm_module.resolve())]
    run = terminal(bounded_run(command, json.dumps(packet).encode()))
    run['command'] = command
    require(run['resource_status'] == 'completed' and run['exit_code'] == 0, 'Node runtime failed')
    wasm = strict_json(run['stdout'])
    require(len(wasm['rows']) == len(packet['cases']), 'Node case count mismatch')
    rows = []
    refusals = []
    for count in packet['invalid_hydrogens']:
        source = json.dumps(dict(atoms=[dict(id='a1', element='C', hydrogens=count)], bonds=[]))
        try:
            chematic.from_moljson(source)
            refusals.append(False)
        except ValueError as error:
            refusals.append('hydrogens' in str(error))
    for case, js in zip(packet['cases'], wasm['rows']):
        require(case['id'] == js['id'], 'Node case order mismatch')
        mol = chematic.from_moljson(case['input'])
        py = dict(id=case['id'], smiles=mol.smiles, output=mol.to_moljson(),
                  common_smiles=chematic.convert_format(case['input'], 'moljson', 'smiles'))
        rows.append(dict(id=case['id'], input=case['input'], expected=external_snapshot(case['expected_smiles']),
                         python=py, wasm=js, python_observation=observe(case, py), wasm_observation=observe(case, js)))
    summary = {engine: dict(preserved=sum(row[engine+'_observation']['preserved'] for row in rows),
                            invalid_h_refused=sum(refusals if engine == 'python' else wasm['refusals']))
               for engine in ('python', 'wasm')}
    path = 'crates/chematic-mol/src/moljson.rs'
    source = subprocess.check_output(['git', 'show', f'{args.baseline_runtime_revision}:{path}'], cwd=ROOT) if args.baseline_runtime_revision else None
    ob_command = ['/opt/homebrew/bin/obabel', '-L', 'formats']
    result = dict(schema='chematic.moljson-h-observation.v1', baseline=args.baseline,
                  head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  runtime_revision=args.baseline_runtime_revision or 'working_tree',
                  packet=identity(packet_path), harness=identity(__file__), production_source=identity(ROOT/path, source),
                  python_package=identity(chematic.__file__), python_extension=identity(extension.__file__),
                  wasm_module=identity(args.wasm_module), wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),
                  rdkit_version=rdBase.rdkitVersion, rdkit_parser=identity(rdmolfiles.__file__),
                  external_settings='RDKit SMILES parse: sanitize=True, removeHs=False; canonical SMILES, formula, atom count, radical electrons. Not a MolJSON cross-reader.',
                  openbabel_binary=identity(ob_command[0]),
                  openbabel_version=terminal(bounded_run([ob_command[0], '-V'])),
                  openbabel_format_inventory=dict(command=ob_command, **terminal(bounded_run(ob_command))),
                  node_run=run, positive_cases=len(rows), invalid_h_cases=len(refusals),
                  summary=summary, python_refusals=refusals, wasm_refusals=wasm['refusals'], cases=rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(gzip.compress((json.dumps(result, indent=2)+'\n').encode(), mtime=0))
    print(json.dumps(summary))
    if not args.baseline:
        require(all(row[e+'_observation']['preserved'] for row in rows for e in ('python', 'wasm')), 'H or graph contract failed')
        require(all(refusals) and all(wasm['refusals']), 'invalid H count accepted')


if __name__ == '__main__':
    main()
