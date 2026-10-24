#!/usr/bin/env python3
"""Paired CDXML class/direction observations against an authored graph oracle."""
from __future__ import annotations
import argparse
import importlib
import json
from pathlib import Path
import subprocess
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run, identity, require, strict_json
from scripts.validation.cdxml_bond_observer import cdxml_fields, compare, observe

NODE = r"""
const fs=require('node:fs'),w=require(process.argv[1]);
const attempt=f=>{try{return {output:f()}}catch(e){return {error:String(e)}}};
process.stdout.write(JSON.stringify(JSON.parse(fs.readFileSync(0,'utf8')).map(c=>{
 const r={id:c.id};
 for(const infer of [false,true])r['wasm_direct_'+infer]=attempt(()=>{
  const m=w.mol_from_cdxml_with_options(c.input,infer);try{return m.canonical_smiles()}finally{m.free()}});
 r.wasm_common=attempt(()=>w.convert_common_format(c.input,'cdxml','smiles'));
 r.wasm_roundtrip=attempt(()=>w.convert_common_format(c.input,'cdxml','cdxml'));
 if(c.expected)r.wasm_writer=attempt(()=>w.convert_common_format(c.expected.smiles,'smiles','cdxml'));
 return r;
})));
"""


def terminal(raw):
    return {key: value.decode('utf-8', errors='replace') if isinstance(value, bytes) else value for key, value in raw.items()}


def attempt(call):
    try:
        return dict(output=call())
    except ValueError as error:
        return dict(error=str(error))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--wasm-module', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--obabel', type=Path, default=Path('/opt/homebrew/bin/obabel'))
    p.add_argument('--baseline', action='store_true')
    args = p.parse_args()
    import chematic
    from rdkit import rdBase
    extension = importlib.import_module('chematic.chematic')
    fixture = ROOT / 'validation/cdxml_bond_meaning/cases.json'
    cases = strict_json(fixture.read_text())['cases']
    for case in cases:
        if case['group'] != 'refusal':
            atoms, bond = cdxml_fields(case['input'])
            require(all(value for value in compare(case['expected'], atoms, bond).values() if value is not None), 'source oracle mismatch')
            require(observe(case, dict(output=case['expected']['smiles']), 'smiles')['preserved'], 'SMILES observer calibration failed')
    raw = terminal(bounded_run(['node', '-e', NODE, str(args.wasm_module.resolve())], json.dumps(cases).encode()))
    require(raw['resource_status'] == 'completed' and raw['exit_code'] == 0, 'Node runner failed')
    node = strict_json(raw['stdout'])
    require(len(node) == len(cases), 'Node row count mismatch')
    routes = {f'{lang}_direct_{infer}': 'smiles' for lang in ('python', 'wasm') for infer in ('false', 'true')}
    for lang in ('python', 'wasm'):
        routes.update({lang + '_common': 'smiles', lang + '_roundtrip': 'cdxml', lang + '_writer': 'cdxml'})
    routes['openbabel'] = 'smiles'
    command = [str(args.obabel.resolve()), '-icdxml', '-osmi', '-xn', '-xa', '-l1']
    rows = []
    for case, w in zip(cases, node, strict=True):
        require(w.pop('id') == case['id'], 'Node order mismatch')
        row = dict(id=case['id'], group=case['group'], input=case['input'], expected=case.get('expected'), **w)
        for infer in (False, True):
            row['python_direct_' + str(infer).lower()] = attempt(lambda: chematic.from_cdxml(case['input'], infer_nondirectional_stereo=infer).smiles)
        row['python_common'] = attempt(lambda: chematic.convert_format(case['input'], 'cdxml', 'smiles'))
        row['python_roundtrip'] = attempt(lambda: chematic.convert_format(case['input'], 'cdxml', 'cdxml'))
        if case.get('expected'):
            row['python_writer'] = attempt(lambda: chematic.convert_format(case['expected']['smiles'], 'smiles', 'cdxml'))
        value = terminal(bounded_run(command, case['input'].encode()))
        if value['resource_status'] == 'completed' and value['exit_code'] == 0:
            value['output'] = value['stdout']
        row['openbabel'] = value
        for route, format in routes.items():
            if route in row:
                row[route]['observation'] = observe(case, row[route], format)
        rows.append(row)
    summary = {route: {group: dict(count=sum(r['group'] == group and route in r for r in rows),
        preserved=sum(r['group'] == group and r.get(route, {}).get('observation', {}).get('preserved', False) for r in rows),
        profile_refused=sum(r['group'] == group and r.get(route, {}).get('observation', {}).get('profile_refused', False) for r in rows),
        output_outside_profile=sum(r['group'] == group and r.get(route, {}).get('observation', {}).get('outcome') == 'output_outside_profile' for r in rows))
        for group in ('control', 'quadruple', 'dative', 'refusal')} for route in routes}
    baseline = '342d1d5f9ad09aa6ade4e3079a902cba856209e2'
    base = '08c2da7100b9eb9ef0f50403e1900ec8a298d722'
    paths = ['crates/chematic-mol/src/cdxml.rs', 'crates/chematic-mol/src/cml/attributes.rs']
    if args.baseline:
        subprocess.run(['git', 'diff', '--exit-code', baseline, base, '--', *paths], cwd=ROOT, check=True)
    record = dict(schema='chematic.cdxml-bond-meaning-report.v1', baseline=args.baseline,
        baseline_runtime_commit=baseline if args.baseline else None, reader_comparison_base=base,
        head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        scope='two-atom authored development packet; bond class/dative direction and graph scalars; independent XML output observer and RDKit SMILES observer; no source/layout/coordinate comparison, schema/physical/stereo/aromaticity/general ranking/speed claim',
        sources=[] if args.baseline else [identity(ROOT / path) for path in paths],
        baseline_sources=[identity(ROOT / path, subprocess.check_output(['git', 'show', baseline + ':' + path], cwd=ROOT)) for path in paths] if args.baseline else [],
        runner=identity(__file__), observer=identity(ROOT / 'scripts/validation/cdxml_bond_observer.py'),
        process_bounds=identity(ROOT / 'scripts/checks/check_interchange_preservation.py'), fixture=identity(fixture),
        python_package=identity(chematic.__file__), python_extension=identity(extension.__file__),
        wasm_module=identity(args.wasm_module), wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),
        openbabel_binary=identity(args.obabel), openbabel_command=command,
        openbabel_version=terminal(bounded_run([str(args.obabel.resolve()), '-V'])), rdkit_observer_version=rdBase.rdkitVersion,
        process_limits=dict(seconds_per_process=10, combined_output_bytes=1048576),
        summary=summary, cases=rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps(summary))
    for route in routes:
        require(summary[route]['control']['preserved'] == 4, route + ' control failed')
        if not args.baseline and route != 'openbabel':
            require(summary[route]['quadruple']['preserved'] == 3 and summary[route]['dative']['preserved'] == 5, route + ' graph failed')
            if not route.endswith('_writer'):
                require(summary[route]['refusal']['profile_refused'] == 12, route + ' refusal failed')


if __name__ == '__main__':
    main()
