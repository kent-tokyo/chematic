#!/usr/bin/env python3
"""Paired CJSON -> CML development comparison with independent XML observation."""
from __future__ import annotations
import argparse
import gzip
import hashlib
import importlib
import json
import math
from pathlib import Path
import struct
import sys
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run, identity, require, strict_json
NS='{http://www.xml-cml.org/schema}'
SYMBOLS={'H':1,'C':6,'N':7,'O':8}
NODE="""
const fs=require('node:fs'),w=require(process.argv[1]);
process.stdout.write(JSON.stringify(JSON.parse(fs.readFileSync(0,'utf8')).map(c=>{
try{return {id:c.id,output:w.convert_common_format(c.input,'cjson','cml')}}
catch(e){return {id:c.id,error:String(e)}}})));
"""


def bits(values):
    return [struct.pack('>d',float(v)).hex() for v in values]


def observe(expected, output):
    result=dict(xml=False,indexed_graph=False,isotopes=False,coordinate_presence=False,coordinate_bits=False,name=False)
    if not output:return result
    try:
        root=ET.fromstring(output)
        if root.tag==NS+'cml':
            molecules=root.findall(NS+'molecule');require(len(molecules)==1,'molecule count');root=molecules[0]
        require(root.tag==NS+'molecule','wrong molecule root')
        result['xml']=True
        atoms=root.findall('./'+NS+'atomArray/'+NS+'atom')
        numbers=[SYMBOLS[a.attrib['elementType']] for a in atoms]
        require(len(set(a.attrib['id'] for a in atoms))==len(atoms),'duplicate atom ID')
        ids={a.attrib['id']:i for i,a in enumerate(atoms)}
        charges=[int(a.get('formalCharge','0')) for a in atoms]
        isotopes=[int(a.get('isotopeNumber',a.get('isotope','0'))) for a in atoms]
        bonds=[]
        for b in root.findall('./'+NS+'bondArray/'+NS+'bond'):
            pair=[ids[r] for r in b.attrib['atomRefs2'].split()];require(len(pair)==2,'bond endpoints')
            token=b.attrib['order'];order={'S':1,'D':2,'T':3,'A':1.5}.get(token)
            if order is None:order=float(token)
            bonds.append([*sorted(pair),order])
        correspondence=numbers==expected['numbers']
        result['indexed_graph']=correspondence and charges==expected['charges'] and sorted(bonds)==sorted([min(a,b),max(a,b),o] for a,b,o in expected['bonds'])
        result['isotopes']=correspondence and isotopes==expected['isotopes']
        presence=[[k in a.attrib for k in ['x3','y3','z3']] for a in atoms]
        if expected['coordinates'] is None:
            result['coordinate_presence']=not any(any(row) for row in presence)
            result['coordinate_bits']=correspondence and result['coordinate_presence']
        else:
            result['coordinate_presence']=len(atoms)==len(expected['coordinates']) and all(all(row) for row in presence)
            if result['coordinate_presence']:
                values=[float(a.attrib[k]) for a in atoms for k in ['x3','y3','z3']]
                result['coordinate_bits']=correspondence and all(math.isfinite(v) for v in values) and bits(values)==bits(v for row in expected['coordinates'] for v in row)
        names=[''.join(n.itertext()) for n in root.findall(NS+'name')]
        if 'title' in root.attrib:names.append(root.attrib['title'])
        result['name']=expected['name'] is None or bool(names) and all(name==expected['name'] for name in names)
    except (ValueError,KeyError,TypeError,IndexError,OverflowError,ET.ParseError):pass
    return result


def terminal(result):
    return {k:v.decode('utf8',errors='replace') if isinstance(v,bytes) else v for k,v in result.items()}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cli',type=Path,required=True)
    parser.add_argument('--wasm-module',type=Path,required=True)
    parser.add_argument('--obabel',type=Path,default=Path('/opt/homebrew/bin/obabel'))
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--baseline',action='store_true')
    args=parser.parse_args()
    import chematic
    extension=importlib.import_module('chematic.chematic')
    packet_path=ROOT/'validation/cjson_cml_projection/cases.json';cases=[]
    for case in strict_json(packet_path.read_text())['valid']:
        require(hashlib.sha256(case['input'].encode()).hexdigest()==case['sha256'],'fixture hash')
        for label,source in [('lf',case['input']),('crlf',case['input'].replace('\n','\r\n')),('compact',json.dumps(strict_json(case['input']),ensure_ascii=False,separators=(',',':')))]:
            cases.append(dict(id=case['id']+'/'+label,input=source,expected=case['expected'],group=case['group']))
    raw=terminal(bounded_run(['node','-e',NODE,str(args.wasm_module.resolve())],json.dumps(cases).encode()))
    require(raw['resource_status']=='completed' and raw['exit_code']==0,'Node failed')
    node=strict_json(raw['stdout']);require(len(node)==len(cases),'Node count')
    rows=[]
    for case,wasm in zip(cases,node):
        require(case['id']==wasm['id'],'Node order');row=dict(case=case,wasm=wasm)
        try:row['python']=dict(output=chematic.convert_format(case['input'],'cjson','cml'))
        except ValueError as error:row['python']=dict(error=str(error))
        commands=[('cli',[str(args.cli.resolve()),'convert','--input-format','cjson','--output-format','cml']),
                  ('openbabel',[str(args.obabel.resolve()),'-icjson','-ocml']),
                  ('openbabel_metadata',[str(args.obabel.resolve()),'-icjson','-ocml','-xp','-xm'])]
        for label,command in commands:
            result=terminal(bounded_run(command,case['input'].encode()));result['command']=command
            if result['resource_status']=='completed' and result['exit_code']==0:result['output']=result['stdout']
            row[label]=result
        for label in ['python','wasm','cli','openbabel','openbabel_metadata']:
            row[label]['observation']=observe(case['expected'],row[label].get('output'))
        rows.append(row)
    labels=['python','wasm','cli','openbabel','openbabel_metadata']
    summary={group:{label:{metric:sum(row[label]['observation'][metric] for row in rows if row['case']['group']==group) for metric in rows[0][label]['observation']} for label in labels} for group in ['common','compatibility']}
    for group in summary:
        for label in labels:
            selected=[row[label] for row in rows if row['case']['group']==group]
            summary[group][label].update(complete_projection=sum(all(r['observation'].values()) for r in selected),
                no_output=sum(not r.get('output') for r in selected),
                malformed_output=sum(bool(r.get('output')) and not r['observation']['xml'] for r in selected))
    sources=['crates/chematic-mol/src/cml.rs','crates/chematic-mol/src/conversion/cjson_cml.rs','crates/chematic-py/python/chematic/__init__.py','crates/chematic-py/src/formats/interchange.rs','crates/chematic-wasm/src/mol_io/common_formats.rs','crates/chematic-cli/src/conversion_format.rs','crates/chematic-cli/src/main.rs']
    record=dict(schema='chematic.cjson-cml-projection-results.v1',baseline=args.baseline,
                baseline_provenance='local c693dfb7 snapshots from prior CML retention work; not a released baseline' if args.baseline else None,
                packet=identity(packet_path),runner=identity(__file__),implementation_sources=[] if args.baseline else [identity(ROOT/p) for p in sources],
                python_package=identity(chematic.__file__),python_extension=identity(extension.__file__),cli=identity(args.cli),
                wasm_module=identity(args.wasm_module),wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),
                openbabel=identity(args.obabel),openbabel_version=terminal(bounded_run([str(args.obabel.resolve()),'-V'])),
                node_terminal=raw,summary=summary,cases=rows)
    args.output.parent.mkdir(exist_ok=True,parents=True)
    args.output.write_bytes(gzip.compress((json.dumps(record,ensure_ascii=False,indent=2)+'\n').encode(),mtime=0));print(json.dumps(summary))
    if not args.baseline:require(all(all(row[label]['observation'].values()) for row in rows for label in labels[:3]),'candidate mismatch')


if __name__=='__main__':main()
