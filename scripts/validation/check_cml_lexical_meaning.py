#!/usr/bin/env python3
"""Bounded CML lexical-layout comparison against an authored graph/H oracle."""
from __future__ import annotations
import argparse
import importlib
import json
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run, identity, require, strict_json
from scripts.validation.check_v3000_coordinate_projection import terminal

GRAPH_KEYS = ('numbers', 'charges', 'isotopes', 'hydrogens', 'atom_maps', 'bonds')
NODE = r"""
const fs=require('node:fs'),w=require(process.argv[1]);
const attempt=f=>{try{return {output:f()}}catch(e){return {error:String(e)}}};
process.stdout.write(JSON.stringify(JSON.parse(fs.readFileSync(0,'utf8')).map(c=>{
 const row={id:c.id};
 for(const [key,read] of [['wasm_default',w.mol_from_cml],['wasm_strict',w.mol_from_cml_strict]])
  row[key]=attempt(()=>{const m=read(c.input);try{return m.canonical_smiles()}finally{m.free()}});
 row.wasm_common=attempt(()=>w.convert_common_format(c.input,'cml','smiles'));
 row.wasm_source=attempt(()=>JSON.parse(w.convert_with_report_json(c.input,'cml','cml','preserve_source')).output);
 return row;
})));
"""


def compare_fields(expected, actual):
    """Unique-element correspondence is sufficient only for this tiny packet."""
    result = dict.fromkeys(GRAPH_KEYS, False)
    try:
        numbers = actual['numbers']
        require(len(set(expected['numbers'])) == len(expected['numbers']), 'oracle requires unique elements')
        require(len(numbers) == len(expected['numbers']) and set(numbers) == set(expected['numbers']), 'atom mismatch')
        order = [numbers.index(n) for n in expected['numbers']]
        inverse = {old: new for new, old in enumerate(order)}
        result['numbers'] = True
        for key in GRAPH_KEYS[1:-1]:
            require(len(actual[key]) == len(numbers), 'atom field length mismatch')
            result[key] = [actual[key][i] for i in order] == expected[key]
        bonds = sorted([min(inverse[a], inverse[b]), max(inverse[a], inverse[b]), o]
                       for a, b, o in actual['bonds'])
        result['bonds'] = bonds == sorted(expected['bonds'])
    except (ValueError, TypeError, KeyError, IndexError):
        pass
    return result


def observe(expected, output):
    """RDKit observes emitted SMILES, never serves as the CML reader/oracle."""
    result = dict(outcome='no_output', preserved=False)
    if not isinstance(output, str) or not output.strip():
        return result
    result['outcome'] = 'invalid_output'
    lines = [s.strip() for s in output.splitlines() if s.strip()]
    if len(lines) != 1:
        return result
    from rdkit import Chem
    params = Chem.SmilesParserParams()
    params.sanitize = False
    params.removeHs = False
    mol = Chem.MolFromSmiles(lines[0].split()[0], params)
    if mol is None:
        return result
    try:
        mol.UpdatePropertyCache(strict=False)
        actual = dict(numbers=[a.GetAtomicNum() for a in mol.GetAtoms()],
                      charges=[a.GetFormalCharge() for a in mol.GetAtoms()],
                      isotopes=[a.GetIsotope() for a in mol.GetAtoms()],
                      hydrogens=[a.GetTotalNumHs(includeNeighbors=True) for a in mol.GetAtoms()],
                      atom_maps=[a.GetAtomMapNum() for a in mol.GetAtoms()],
                      bonds=[[b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble()]
                             for b in mol.GetBonds()])
        fields = compare_fields(expected, actual)
        result.update(outcome='smiles_output', actual=actual, fields=fields, preserved=all(fields.values()))
    except (ValueError, RuntimeError):
        pass
    return result


def source_fields(text):
    root = ET.fromstring(text)
    atoms = root.findall('.//{*}atom')
    ids = {a.attrib['id']: i for i, a in enumerate(atoms)}
    return dict(numbers=[{'C': 6, 'N': 7, 'O': 8}[a.attrib['elementType']] for a in atoms],
                charges=[int(a.attrib.get('formalCharge', 0)) for a in atoms],
                isotopes=[int(a.attrib.get('isotopeNumber', 0)) for a in atoms],
                hydrogens=[int(a.attrib['hydrogenCount']) for a in atoms],
                atom_maps=[0] * len(atoms),
                bonds=[[*(ids[r] for r in b.attrib['atomRefs2'].split()), int(b.attrib['order'])]
                       for b in root.findall('.//{*}bond')])


def attempt(call):
    try:
        return dict(output=call())
    except ValueError as error:
        return dict(error=str(error))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--wasm-module', type=Path, required=True)
    p.add_argument('--obabel', type=Path, default=Path('/opt/homebrew/bin/obabel'))
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--baseline', action='store_true')
    args = p.parse_args()
    import chematic
    from rdkit import rdBase
    extension = importlib.import_module('chematic.chematic')
    fixture = ROOT / 'validation/cml_lexical_meaning/cases.json'
    cases = strict_json(fixture.read_text())['cases']
    for case in cases:
        require(source_fields(case['input']) == {k: case['expected'][k] for k in GRAPH_KEYS}, 'source oracle mismatch')
    raw = terminal(bounded_run(['node', '-e', NODE, str(args.wasm_module.resolve())], json.dumps(cases).encode()))
    require(raw['resource_status'] == 'completed' and raw['exit_code'] == 0, 'Node runner failed')
    node = strict_json(raw['stdout'])
    require(len(node) == len(cases), 'Node count mismatch')
    results = []
    routes = ('python_default', 'python_strict', 'python_default_bytes', 'python_strict_bytes',
              'python_common', 'python_3d', 'wasm_default', 'wasm_strict', 'wasm_common', 'openbabel')
    command = [str(args.obabel.resolve()), '-icml', '-osmi', '-xn', '-xa']
    for case, w in zip(cases, node, strict=True):
        require(w.pop('id') == case['id'], 'Node row order mismatch')
        row = dict(id=case['id'], scope=case['scope'], input=case['input'], expected=case['expected'], **w)
        for label, reader in [('default', chematic.from_cml), ('strict', chematic.from_cml_strict)]:
            for suffix, text in [('', case['input']), ('_bytes', case['input'].encode())]:
                row['python_' + label + suffix] = attempt(lambda: reader(text).smiles)
        row['python_common'] = attempt(lambda: chematic.convert_format(case['input'], 'cml', 'smiles'))
        row['python_3d'] = attempt(lambda: chematic.from_cml_with_3d(case['input'])[0].smiles)
        row['python_source'] = attempt(lambda: chematic.convert_with_report(
            case['input'], 'cml', 'cml', policy='preserve_source')['output'])
        value = terminal(bounded_run(command, case['input'].encode()))
        if value['resource_status'] == 'completed' and value['exit_code'] == 0:
            value['output'] = value['stdout']
        row['openbabel'] = value
        for route in routes:
            row[route]['observation'] = observe(case['expected'], row[route].get('output'))
        for route in ('python_source', 'wasm_source'):
            row[route]['exact_source_bytes'] = row[route].get('output') == case['input']
        results.append(row)
    summary = {route: {group: dict(count=sum(r['scope'] == group for r in results),
                                  preserved=sum(r['scope'] == group and r[route]['observation']['preserved'] for r in results),
                                  no_output=sum(r['scope'] == group and r[route]['observation']['outcome'] == 'no_output' for r in results))
                       for group in ('control', 'reproduction')} for route in routes}
    source_counts = {route: sum(r[route]['exact_source_bytes'] for r in results)
                     for route in ('python_source', 'wasm_source')}
    source_paths = ['crates/chematic-mol/src/cml.rs', 'crates/chematic-mol/src/cml/attributes.rs',
                    'crates/chematic-mol/src/cml/lexical.rs',
                    'crates/chematic-mol/src/cml/geometry.rs', 'crates/chematic-mol/src/cml/geometry_xml.rs']
    baseline_commit = '342d1d5f9ad09aa6ade4e3079a902cba856209e2'
    comparison_base = 'fcbf6d0a8c901dee5220cca965cf119c0a3b7b46'
    baseline_sources = [path for path in source_paths if not path.endswith('/lexical.rs')]
    if args.baseline:
        subprocess.run(['git', 'diff', '--exit-code', baseline_commit, comparison_base, '--', *baseline_sources], cwd=ROOT, check=True)
    record = dict(schema='chematic.cml-lexical-meaning-report.v1', baseline=args.baseline,
                  baseline_runtime_commit=baseline_commit if args.baseline else None,
                  reader_comparison_base=comparison_base,
                  baseline_reader_sources=[identity(ROOT / path, subprocess.check_output(
                      ['git', 'show', baseline_commit + ':' + path], cwd=ROOT)) for path in baseline_sources] if args.baseline else [],
                  scope='synthetic development packet; graph/H observation only; native tests cover 2D coordinates separately; no schema/stereo/aromaticity/query/speed claim',
                  head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  sources=[] if args.baseline else [identity(ROOT / path) for path in source_paths],
                  runner=identity(__file__), observer=identity(ROOT / 'scripts/checks/check_interchange_preservation.py'),
                  fixture=identity(fixture), python_package=identity(chematic.__file__), python_extension=identity(extension.__file__),
                  wasm_module=identity(args.wasm_module), wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),
                  openbabel_binary=identity(args.obabel), openbabel_command=command,
                  openbabel_version=terminal(bounded_run([str(args.obabel.resolve()), '-V'])), rdkit_observer_version=rdBase.rdkitVersion,
                  summary=summary, exact_source_bytes=source_counts, cases=results)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps(dict(summary=summary, exact_source_bytes=source_counts)))
    if not args.baseline:
        for route in routes[:-1]:
            require(all(r[route]['observation']['preserved'] for r in results), route + ' failed')
        require(all(v == len(cases) for v in source_counts.values()), 'source retention failed')


if __name__ == '__main__':
    main()
