#!/usr/bin/env python3
"""Bounded CJSON common-bridge graph validation development comparison."""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import math
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run, identity, require, strict_json

NODE = """
const fs=require('node:fs'),w=require(process.argv[1]);
const cases=JSON.parse(fs.readFileSync(0,'utf8'));
process.stdout.write(JSON.stringify(cases.map(c=>{try{return {id:c.id,output:w.convert_common_format(c.input,'cjson','cml')}}catch(e){return {id:c.id,error:String(e)}}})));
"""


def compact_json(text):
    """Strip only JSON whitespace outside strings; keep duplicate keys intact."""
    result, quoted, escaped = [], False, False
    for c in text:
        if quoted:
            result.append(c)
            if escaped:
                escaped = False
            elif c == '\\':
                escaped = True
            elif c == '"':
                quoted = False
        elif c == '"':
            quoted = True
            result.append(c)
        elif c not in ' \t\r\n':
            result.append(c)
    return ''.join(result)


def integer_array(values, count, field):
    require(isinstance(values, list) and len(values) == count, field+' length/type')
    require(all(isinstance(v, int) and not isinstance(v, bool) for v in values), field+' integer type')
    return values


def read_graph(data):
    """Independent narrow observer: typed indexed arrays, no chemistry perception."""
    atoms = data['atoms']
    numbers = atoms['elements']['number']
    require(isinstance(numbers, list), 'atomic numbers array')
    n = len(numbers)
    integer_array(numbers, n, 'numbers')
    require(all(1 <= v <= 118 for v in numbers), 'atomic number range')
    charges = integer_array(atoms.get('formalCharges', [0]*n), n, 'charges')
    isotopes = integer_array(atoms.get('isotopes', [0]*n), n, 'isotopes')
    bonds = data.get('bonds', {})
    indices = bonds.get('connections', {}).get('index', [])
    require(isinstance(indices, list) and len(indices)%2 == 0, 'bond indices pairs')
    integer_array(indices, len(indices), 'bond indices')
    require(all(0 <= v < n for v in indices), 'bond index range')
    orders = bonds.get('order', [1]*(len(indices)//2))
    require(isinstance(orders, list) and len(orders)*2 == len(indices), 'bond order count')
    require(all(isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) for v in orders), 'bond order type')
    return dict(numbers=numbers, charges=charges, isotopes=isotopes, indices=indices, orders=orders)


def read_cml_graph(text):
    """Narrow standard-library XML observer for packet elements/indexed bonds."""
    root = ET.fromstring(text)
    local = lambda node: node.tag.rsplit('}',1)[-1]
    molecules = [node for node in root.iter() if local(node)=='molecule']
    require(len(molecules)==1, 'one output molecule required')
    atoms = [node for node in molecules[0].iter() if local(node)=='atom']
    symbols = {'H':1,'C':6,'N':7,'O':8,'Na':11,'Cl':17}
    ids, numbers, charges, isotopes = {}, [], [], []
    for atom in atoms:
        id = atom.attrib['id']
        require(id not in ids, 'duplicate output atom ID')
        ids[id] = len(numbers)
        numbers.append(symbols[atom.attrib['elementType']])
        charges.append(int(atom.get('formalCharge','0')))
        isotopes.append(int(atom.get('isotopeNumber',atom.get('isotope','0'))))
    indices, orders = [], []
    tokens = {'1':1,'S':1,'2':2,'D':2,'3':3,'T':3,'A':1.5,'4':4}
    for bond in molecules[0].iter():
        if local(bond)!='bond':
            continue
        refs = bond.attrib['atomRefs2'].split()
        require(len(refs)==2, 'bond reference pair required')
        indices.extend(ids[id] for id in refs)
        token = bond.get('order','1')
        order = tokens[token] if token in tokens else float(token)
        require(math.isfinite(order), 'finite output bond order required')
        orders.append(order)
    return dict(numbers=numbers,charges=charges,isotopes=isotopes,indices=indices,orders=orders)


def canonical_graph(graph):
    edges = sorted((min(a,b), max(a,b), order) for a,b,order in
                   zip(graph['indices'][::2], graph['indices'][1::2], graph['orders']))
    return {key:graph[key] for key in ['numbers','charges','isotopes']}, edges


def observe(case, result):
    if result.get('resource_status', 'completed') != 'completed':
        return {'outcome':'resource_failure', 'selected_graph':False}
    output = result.get('output', '')
    if not output.strip():
        error_status = 'error' in result or result.get('exit_code', 0) != 0
        return {'outcome':'refused_no_output' if error_status else 'empty_output', 'selected_graph':False}
    if result.get('exit_code', 0) != 0:
        return {'outcome':'error_with_output', 'selected_graph':False}
    try:
        if result.get('format')=='cml':
            # Emitting parseable XML from invalid source is separate from
            # whether the output fits this narrow typed graph observer.
            ET.fromstring(output)
            graph = read_cml_graph(output) if case['group']!='invalid' else None
        else:
            data = strict_json(output)
            require(isinstance(data, dict), 'CJSON object required')
            graph = read_graph(data) if case['group']!='invalid' else None
        if case['group'] == 'invalid':
            return {'outcome':'accepted_outside_profile', 'selected_graph':False}
        matches = canonical_graph(graph) == canonical_graph(case['expected'])
        return {'outcome':'graph_match' if matches else 'graph_mismatch', 'selected_graph':matches}
    except (ValueError, KeyError, TypeError, AttributeError, ET.ParseError) as error:
        return {'outcome':'invalid_output', 'selected_graph':False, 'reason':str(error)}


def terminal(raw):
    return {key:value.decode('utf-8', errors='replace') if isinstance(value,bytes) else value for key,value in raw.items()}


def cli_run(command, text):
    raw = terminal(bounded_run(command, text.encode(), timeout=3))
    return {**raw, 'command':command, 'output':raw['stdout']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wasm-module', type=Path, required=True)
    parser.add_argument('--cli', type=Path)
    parser.add_argument('--obabel', type=Path, default=Path('/opt/homebrew/bin/obabel'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--baseline', action='store_true')
    args = parser.parse_args()
    require(args.baseline or args.cli is not None, 'candidate CLI required')
    import chematic
    extension = importlib.import_module('chematic.chematic')
    packet_path = ROOT/'validation/cjson_checked/cases.json'
    packet = strict_json(packet_path.read_text())
    cases = []
    for group in ['valid','compatibility','invalid']:
        for c in packet[group]:
            require(hashlib.sha256(c['input'].encode()).hexdigest()==c['sha256'], 'fixture SHA mismatch')
            layouts = [('lf',c['input']),('crlf',c['input'].replace('\n','\r\n')),('compact',compact_json(c['input']))]
            for layout,text in layouts:
                cases.append({**c,'group':group,'id':c['id']+'/'+layout,'input':text})
    command = ['node','-e',NODE,str(args.wasm_module.resolve())]
    node_input = json.dumps(cases).encode()
    node_run = bounded_run(command,node_input)
    require(node_run['resource_status']=='completed' and node_run['exit_code']==0, 'Node failed')
    nodes = strict_json(node_run['stdout'])
    require(len(nodes)==len(cases), 'Node count mismatch')
    labels = ['python','wasm','openbabel']+(['cli'] if args.cli else [])
    results = []
    for case,node in zip(cases,nodes):
        require(case['id']==node['id'], 'Node order mismatch')
        row = dict(id=case['id'],group=case['group'],input_sha256=hashlib.sha256(case['input'].encode()).hexdigest(),wasm={**node,'format':'cml'})
        try:
            row['python'] = {'output':chematic.convert_format(case['input'],'cjson','cml'),'format':'cml'}
        except ValueError as error:
            row['python'] = {'error':str(error)}
        row['openbabel'] = cli_run([str(args.obabel.resolve()),'-icjson','-ocml'],case['input'])
        if args.cli:
            row['cli'] = cli_run([str(args.cli.resolve()),'convert','--input-format','cjson','--output-format','cml'],case['input'])
        for label in labels:
            row[label]['format'] = 'cml'
            row[label]['observation'] = observe(case,row[label])
            if not args.baseline and label != 'openbabel':
                expected = 'refused_no_output' if case['group']=='invalid' else 'graph_match'
                require(row[label]['observation']['outcome']==expected, f'{label} {case["id"]}: {row[label]["observation"]}')
        results.append(row)
    summary = {}
    for label in labels:
        summary[label] = {}
        for group in ['valid','compatibility','invalid']:
            rows = [r[label]['observation'] for r in results if r['group']==group]
            summary[label][group] = {'cases':len(rows),'selected_graph':sum(r['selected_graph'] for r in rows),
                'outcomes':{outcome:sum(r['outcome']==outcome for r in rows) for outcome in sorted({r['outcome'] for r in rows})}}
    record = dict(schema='chematic.cjson-checked-results.v1',baseline=args.baseline,route='cjson-to-cml',packet=identity(packet_path),
                  python_package=identity(chematic.__file__),python_extension=identity(extension.__file__),
                  wasm_module=identity(args.wasm_module),wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),
                  openbabel=identity(args.obabel),openbabel_version=terminal(bounded_run([str(args.obabel.resolve()),'-V'])),
                  node_run={**terminal(node_run),'command':['node','-e','<NODE source in hashed runner>',str(args.wasm_module.resolve())],'input_sha256':hashlib.sha256(node_input).hexdigest()},
                  summary=summary,cases=results)
    if args.cli:
        record['cli'] = identity(args.cli)
    args.output.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(summary))


if __name__=='__main__':
    main()
