#!/usr/bin/env python3
"""Indexed CJSON→MOL2 formal/partial charge diagnostic with an external reader.

The Open Babel C++ probe reads candidate MOL2, keeping stored partial charges
without lazy charge calculation. No engine supplies the expected source values.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import math
from pathlib import Path
import struct
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run, identity, require, strict_json

NODE="""
const fs=require('node:fs'), w=require(process.argv[1]);
const cases=JSON.parse(fs.readFileSync(0,'utf8'));
process.stdout.write(JSON.stringify(cases.map(c=>({id:c.id,
 report:JSON.parse(w.convert_with_report_json(c.input,'cjson','mol2','allow_loss'))}))));
"""

ELEMENTS = [
    "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
    "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
    "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",
    "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
    "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv",
    "Ts", "Og",
]
ATOMIC_NUMBER = {symbol: number for number, symbol in enumerate(ELEMENTS) if symbol}


def _mol2_section(lines, name):
    marker = f"@<TRIPOS>{name}"
    try:
        start = lines.index(marker) + 1
    except ValueError:
        return []
    end = next((i for i in range(start, len(lines)) if lines[i].startswith("@<TRIPOS>")), len(lines))
    return [line for line in lines[start:end] if line.strip()]


def mol2_snapshot(text):
    """Read the explicit graph fields emitted by the bounded CJSON adapter."""
    lines = text.replace("\r\n", "\n").splitlines()
    atoms = [line.split() for line in _mol2_section(lines, "ATOM")]
    require(all(len(row) >= 6 for row in atoms), "invalid MOL2 atom row")
    ids = [int(row[0]) for row in atoms]
    require(len(ids) == len(set(ids)), "duplicate MOL2 atom id")
    positions = {source_id: index for index, source_id in enumerate(ids)}
    symbols = [row[5].split(".", 1)[0] for row in atoms]
    require(all(symbol in ATOMIC_NUMBER for symbol in symbols), "unknown MOL2 atom type")
    charges = [1 if row[5] == "N.4" else 0 for row in atoms]
    unity = _mol2_section(lines, "UNITY_ATOM_ATTR")
    cursor = 0
    while cursor < len(unity):
        header = unity[cursor].split()
        require(len(header) == 2 and int(header[1]) >= 0, "invalid UNITY header")
        source_id, count = int(header[0]), int(header[1])
        require(source_id in positions and cursor + count < len(unity), "invalid UNITY atom id")
        for offset in range(1, count + 1):
            key, value = unity[cursor + offset].split(maxsplit=1)
            if key == "charge":
                charges[positions[source_id]] = int(value)
        cursor += count + 1
    bonds = []
    for row in (line.split() for line in _mol2_section(lines, "BOND")):
        require(len(row) >= 4, "invalid MOL2 bond row")
        bonds.extend([positions[int(row[1])], positions[int(row[2])]])
    return {
        "elements": [ATOMIC_NUMBER[symbol] for symbol in symbols],
        "formal_charges": charges,
        "coordinates": [float(value) for row in atoms for value in row[2:5]],
        "bonds": bonds,
    }


def read_probe(text):
    lines=text.splitlines()
    require(bool(lines),'missing probe count')
    n,b=map(int,lines[0].split())
    require(0<=n<=64 and 0<=b<=128 and len(lines)==1+n+b,'probe count mismatch')
    atoms=[s.split() for s in lines[1:1+n]]
    bonds=[list(map(int,s.split())) for s in lines[1+n:]]
    require(all(len(a)==6 for a in atoms) and all(len(edge)==3 for edge in bonds),'probe row shape')
    value=dict(atomic_numbers=[int(a[0]) for a in atoms],formal_charges=[int(a[1]) for a in atoms],
               partial_charges=[float(a[2]) for a in atoms],coords=[[float(v) for v in a[3:]] for a in atoms],bonds=bonds)
    require(all(math.isfinite(v) for v in value['partial_charges']+[v for p in value['coords'] for v in p]),'nonfinite probe data')
    require(all(0<=a<n and 0<=z<n and a!=z for a,z,_ in bonds),'invalid probe bond')
    return value


def observe(expected, actual):
    if actual is None:
        return dict(indexed_graph=False,formal_charges=False,partial_charges_exact=False,coordinates_exact=False)
    edges=lambda rows:sorted((min(a,b),max(a,b),order) for a,b,order in rows)
    graph=actual['atomic_numbers']==expected['atomic_numbers'] and edges(actual['bonds'])==edges(expected['bonds'])
    bits=lambda xs:[struct.pack('>d',v) for v in xs]
    return dict(indexed_graph=graph,formal_charges=graph and actual['formal_charges']==expected['formal_charges'],
                partial_charges_exact=graph and bits(actual['partial_charges'])==bits(expected['partial_charges']),
                coordinates_exact=graph and bits([v for p in actual['coords'] for v in p])==bits([v for p in expected['coords'] for v in p]))


def terminal(raw):
    return {k:v.decode('utf-8',errors='replace') if isinstance(v,bytes) else v for k,v in raw.items()}


def external_read(probe,text):
    raw=terminal(bounded_run([str(probe.resolve())],text.encode()))
    if raw['resource_status']=='completed' and raw['exit_code']==0:
        raw['data']=read_probe(raw['stdout'])
    return raw


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--openbabel-probe',required=True,type=Path)
    parser.add_argument('--obabel',default=Path('/opt/homebrew/bin/obabel'),type=Path)
    parser.add_argument('--wasm-module',type=Path)
    parser.add_argument('--baseline',action='store_true')
    args=parser.parse_args()
    import chematic
    extension=importlib.import_module('chematic.chematic')
    packet_path=ROOT/'validation/cjson_mol2_formal_charges/cases.json'
    packet=strict_json(packet_path.read_text())
    cases=[]
    for c in packet['cases']:
        require(hashlib.sha256(c['input'].encode()).hexdigest()==c['sha256'],'fixture hash mismatch')
        for name,text in [('lf',c['input']),('crlf',c['input'].replace('\n','\r\n'))]:
            cases.append(dict(id=c['id']+'/'+name,input=text,expected=c['expected']))
    wasm=None
    if args.wasm_module:
        run=bounded_run(['node','-e',NODE,str(args.wasm_module.resolve())],json.dumps(cases).encode())
        require(run['resource_status']=='completed' and run['exit_code']==0,'Node/WASM failed')
        wasm=strict_json(run['stdout'])
        require(len(wasm)==len(cases),'WASM count mismatch')
    results=[]
    for i,c in enumerate(cases):
        report=chematic.convert_with_report(c['input'],'cjson','mol2',policy='allow_loss')
        require(report['output'] is not None,'native candidate missing')
        row=dict(id=c['id'],input_sha256=hashlib.sha256(c['input'].encode()).hexdigest(),expected=c['expected'],python_report=report)
        row['candidate_cross_read']=external_read(args.openbabel_probe,report['output'])
        row['candidate_observation']=observe(c['expected'],row['candidate_cross_read'].get('data'))
        if wasm is not None:
            require(wasm[i]['id']==c['id'],'WASM order mismatch')
            row['wasm_report']=wasm[i]['report']
            row['wasm_python_equal']=wasm[i]['report']==report
        command=[str(args.obabel.resolve()),'-icjson','-omol2']
        converted=terminal(bounded_run(command,c['input'].encode()))
        converted['command']=command
        row['openbabel_conversion']=converted
        after=external_read(args.openbabel_probe,converted['stdout']) if converted['resource_status']=='completed' and converted['exit_code']==0 else {}
        row['openbabel_cross_read']=after
        row['openbabel_observation']=observe(c['expected'],after.get('data'))
        results.append(row)
    metrics=['indexed_graph','formal_charges','partial_charges_exact','coordinates_exact']
    summary={e:{m:sum(r[e+'_observation'][m] for r in results) for m in metrics} for e in ['candidate','openbabel']}
    result=dict(schema='chematic.cjson-mol2-formal-charge-validation.v1',baseline=args.baseline,packet=identity(packet_path),
        python_package=identity(chematic.__file__),python_extension=identity(extension.__file__),
        openbabel_probe=identity(args.openbabel_probe),openbabel=identity(args.obabel),
        openbabel_version=terminal(bounded_run([str(args.obabel.resolve()),'-V'])),case_count=len(results),summary=summary,cases=results)
    if args.wasm_module:
        result['wasm_module']=identity(args.wasm_module);result['wasm_binary']=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm'))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(summary))
    if not args.baseline:
        require(all(all(r['candidate_observation'].values()) for r in results),'candidate indexed data contract failed')
        require(all(r.get('wasm_python_equal',True) for r in results),'WASM/Python report mismatch')


if __name__=='__main__':
    main()
