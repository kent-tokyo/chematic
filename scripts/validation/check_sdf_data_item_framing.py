#!/usr/bin/env python3
"""Paired SDF field read/write development cases; select native build via PYTHONPATH."""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import platform
import re
import sys

from scripts.check_interchange_preservation import bounded_run, identity, require

ROOT = Path(__file__).resolve().parents[2]
PYTHON_CALL = '''import io,json,pathlib,sys,tempfile
import chematic
data=sys.stdin.buffer.read()
records=list(chematic.ForwardSDMolSupplier(io.BytesIO(data)))
if len(records)!=1 or records[0] is None: raise ValueError('expected one record')
mol=records[0]
props={k:mol.GetProp(k) for k in mol.GetPropNames() if not k.startswith('_')}
with tempfile.TemporaryDirectory(prefix='chematic-sdf-fields-') as directory:
 path=pathlib.Path(directory)/'output.sdf'
 writer=chematic.SDWriter(path,compute2d=False)
 writer.SetForceV3000(sys.argv[1]=='v3000')
 writer.write(mol);writer.close()
 output=path.read_text()
print(json.dumps({'parsed_properties':props,'output':output}))'''


def properties(text):
    # Independent block observer: split data items on empty physical lines,
    # then parse only the first line of each item as its header.
    normalized = text.replace('\r\n', '\n')
    require(normalized.count('$$$$') == 1, 'expected one record delimiter')
    body = normalized.split('M  END\n', 1)[1].rsplit('$$$$', 1)[0]
    result = {}
    for item in re.split(r'\n{2,}', body.strip('\n')):
        if not item:
            continue
        header, _, value = item.partition('\n')
        match = re.fullmatch(r'>[^<]*<([^>]*)>[^\n]*', header)
        require(match is not None, 'unresolved SD header')
        key = match.group(1)
        require(key not in result, 'duplicate key outside this packet')
        result[key] = value
    return result


def graph(text):
    lines = text.replace('\r\n', '\n').split('M  END\n',1)[0].splitlines()
    if lines[3].endswith('V3000'):
        counts = [line.split()[3:5] for line in lines if line.startswith('M  V30 COUNTS ')]
        require(counts == [['2','1']], 'V3000 counts changed')
        atoms = [line.split()[3] for line in lines if re.match(r'^M  V30 [12] [CO] ', line)]
        bonds = [line.split()[3:6] for line in lines if line.startswith('M  V30 1 1 1 2')]
        require(atoms == ['C', 'O'] and bonds == [['1', '1', '2']], 'V3000 graph changed')
        require(not any('CHG=' in line or 'MASS=' in line for line in lines), 'unexpected atom annotations')
    else:
        require((int(lines[3][:3]), int(lines[3][3:6])) == (2, 1), 'counts changed')
        require([line[31:34].strip() for line in lines[4:6]] == ['C', 'O'], 'elements changed')
        require([int(lines[6][a:b]) for a,b in [(0,3),(3,6),(6,9)]] == [1,2,1], 'bond changed')
        require(not any(line.startswith(('M  CHG', 'M  ISO')) for line in lines), 'charge/isotope changed')
        require(all(int(line[34:36]) == int(line[36:39]) == 0 for line in lines[4:6]), 'atom charge/isotope changed')
    return 'C-O single bond, neutral unspecified isotopes'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--obabel', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), 'output already exists')
    manifest = ROOT/'validation/sdf_data_item_framing/cases.json'
    packet = json.loads(manifest.read_text())
    version = bounded_run([str(args.obabel), '-V'])
    require(version['exit_code'] == 0 and b'Open Babel 3.2.1' in version['stdout'], 'expected Open Babel 3.2.1')
    import importlib
    extension = importlib.import_module('chematic.chematic')
    rows, cases = [], []
    for fmt, ctab in packet['ctabs'].items():
        for ending in ['lf', 'crlf']:
            for case in packet['cases']:
                source = ctab+case['data']+'$$$$\n'
                if ending == 'crlf': source = source.replace('\n', '\r\n')
                case_id = '/'.join([fmt, ending, case['id']])
                expected = case['expected_properties']
                require(properties(source) == expected, 'bad fixture oracle')
                graph(source)
                cases.append({'id':case_id, 'source':source, 'expected_properties':expected})
                for engine, command in [
                    ('chematic', [sys.executable, '-c', PYTHON_CALL, fmt]),
                    ('openbabel', [str(args.obabel), '-isdf', '-osdf', *(['-x3'] if fmt == 'v3000' else [])]),
                ]:
                    run = bounded_run(command, source.encode(), timeout=15, max_output=262144)
                    row = {'case':case_id, 'engine':engine, 'command':command,
                           **{k:v.decode('utf8', errors='replace') if isinstance(v,bytes) else v for k,v in run.items()},
                           'terminal':'process_failure'}
                    if run['resource_status'] != 'completed': row['terminal'] = run['resource_status']
                    elif run['exit_code'] == 0:
                        try:
                            if engine == 'chematic':
                                result = json.loads(row['stdout'])
                                row['parsed_properties'] = result['parsed_properties']
                                output = result['output']
                            else: output = row['stdout']
                            row['output'] = output
                            row['graph'] = graph(output)
                            row['output_properties'] = properties(output)
                            row['terminal'] = 'fields_preserved' if row['output_properties'] == expected and row.get('parsed_properties',expected) == expected else 'field_loss'
                        except (ValueError, KeyError, IndexError) as error:
                            row['terminal'],row['error'] = 'unresolved_output',str(error)
                    rows.append(row)
    result = {'schema':'chematic.sdf-data-item-framing-development.v1',
              'scope':'48 original V2000/V3000 LF/CRLF read/write variants, unique field names and exact values plus explicit C/O graph only. No header metadata/order/repeated-key, full CTAB chemistry, independent holdout, timing or toolkit-wide claim. Failures are not coverage.',
              'manifest':identity(manifest), 'fixture_metadata':{k:v for k,v in packet.items() if k not in ('cases','ctabs')},
              'oracle':{'name':'Authored strings plus independent blank-separated block observer','version':platform.python_version(),'script':identity(__file__),'tolerance':'exact field strings, no whitespace trimming'},
              'engines':{'chematic_extension':identity(extension.__file__),'openbabel':identity(args.obabel)},
              'openbabel_version':version['stdout'].decode().strip(),'cases':cases,'rows':rows,
              'counts':{e:dict(Counter(r['terminal'] for r in rows if r['engine']==e)) for e in ['chematic','openbabel']}}
    with args.output.open('x') as output: output.write(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result['counts'],sort_keys=True))


if __name__ == '__main__': main()
