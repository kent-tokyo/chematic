#!/usr/bin/env python3
"""Compare namespace-aware CML projection; retain every output and refusal."""
from __future__ import annotations
import argparse
import importlib
import json
from pathlib import Path
import struct
import subprocess
import sys
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run, identity, require, strict_json
from scripts.validation.check_v3000_coordinate_projection import observe, terminal

NODE="""
const fs=require('node:fs'),w=require(process.argv[1]);
process.stdout.write(JSON.stringify(JSON.parse(fs.readFileSync(0,'utf8')).map(c=>{
 try{return {output:w.convert_common_format(c.input,'cml','cjson')}}
 catch(e){return {error:String(e)}}})));
"""


def source_fields(text):
    """Independent ElementTree expanded-name observer, positive fixtures only."""
    root=ET.fromstring(text)
    mol=root if root.tag.endswith('molecule') else next(c for c in root if c.tag.endswith('molecule'))
    ns=mol.tag[:-len('molecule')]
    atoms=[];bonds=[]
    for child in mol:
        if child.tag==ns+'atom':atoms.append(child)
        if child.tag==ns+'bond':bonds.append(child)
        if child.tag==ns+'atomArray':atoms.extend(c for c in child if c.tag==ns+'atom')
        if child.tag==ns+'bondArray':bonds.extend(c for c in child if c.tag==ns+'bond')
    ids={a.attrib['id']:i for i,a in enumerate(atoms)}
    numbers={'C':6,'O':8,'N':7,'Cl':17}
    return dict(numbers=[numbers[a.attrib['elementType']] for a in atoms],
        coordinates=[[float(a.attrib[k]) for k in ('x3','y3','z3')] for a in atoms],
        formal_charges=[int(a.attrib.get('formalCharge',0)) for a in atoms],
        isotopes=[int(a.attrib.get('isotopeNumber',0)) for a in atoms],
        bonds=[[*(ids[r] for r in b.attrib['atomRefs2'].split()),int(b.attrib['order'])] for b in bonds])


def cases(packet):
    result=[]
    for group in ('valid','invalid'):
        for c in packet[group]:
            if group=='valid':
                observed=source_fields(c['input'])
                require(observed==c['expected'],'independent source fields disagree')
                bits=lambda xyz:[struct.pack('>d',x) for row in xyz for x in row]
                require(bits(observed['coordinates'])==bits(c['expected']['coordinates']),'source coordinate bits disagree')
            lf=c['input'].replace('><','>\n<')
            for variant,text in [('compact',c['input']),('lf',lf),('crlf',lf.replace('\n','\r\n'))]:
                result.append(dict(id=c['id']+'/'+variant,input=text,group=group,expected=c.get('expected')))
    return result


def observation(case,result):
    output=result.get('output','')
    if not output.strip():return dict(outcome='no_output',exact=False,graph_coordinates=False)
    if case['group']=='invalid':return dict(outcome='accepted_outside_profile',exact=False,graph_coordinates=False)
    fields=observe(case['expected'],output)
    return dict(outcome='json_output' if fields['json_object'] else 'invalid_output',
                exact=all(fields.values()),graph_coordinates=fields['indexed_graph'] and fields['coordinate_bits'],fields=fields)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--wasm-module',type=Path,required=True)
    p.add_argument('--obabel',type=Path,default=Path('/opt/homebrew/bin/obabel'))
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--baseline',action='store_true');args=p.parse_args()
    import chematic
    extension=importlib.import_module('chematic.chematic')
    fixture=ROOT/'validation/cml_namespace_projection/cases.json'
    inputs=cases(strict_json(fixture.read_text()))
    raw=bounded_run(['node','-e',NODE,str(args.wasm_module.resolve())],json.dumps(inputs).encode())
    require(raw['resource_status']=='completed' and raw['exit_code']==0,'Node failed')
    node=strict_json(raw['stdout']);require(len(node)==len(inputs),'Node row count mismatch')
    results=[]
    for c,w in zip(inputs,node,strict=True):
        row=dict(id=c['id'],group=c['group'],input=c['input'],expected=c['expected'],wasm=w)
        try:row['python']=dict(output=chematic.convert_format(c['input'],'cml','cjson'))
        except ValueError as error:row['python']=dict(error=str(error))
        for label,options in [('openbabel',[]),('openbabel_single',['-xs'])]:
            command=[str(args.obabel.resolve()),'-icml','-ocjson']+options
            value=terminal(bounded_run(command,c['input'].encode()))
            value['command']=command
            if value['resource_status']=='completed' and value['exit_code']==0:value['output']=value['stdout']
            row[label]=value
        for label in ('python','wasm','openbabel','openbabel_single'):
            row[label]['observation']=observation(c,row[label])
        results.append(row)
    summary={}
    for label in ('python','wasm','openbabel','openbabel_single'):
        positive=[r[label]['observation'] for r in results if r['group']=='valid']
        negative=[r[label]['observation'] for r in results if r['group']=='invalid']
        summary[label]=dict(valid_count=len(positive),valid_exact=sum(r['exact'] for r in positive),
            valid_graph_coordinates=sum(r['graph_coordinates'] for r in positive),
            valid_no_output=sum(r['outcome']=='no_output' for r in positive),
            invalid_count=len(negative),invalid_no_output=sum(r['outcome']=='no_output' for r in negative))
    paths=['crates/chematic-mol/src/cml.rs','crates/chematic-mol/src/cml/geometry.rs','crates/chematic-mol/src/cml/geometry_xml.rs',
           'crates/chematic-py/src/formats/cml.rs','crates/chematic-py/python/chematic/__init__.py',
           'crates/chematic-wasm/src/mol_io/common_formats.rs',
           'scripts/validation/check_cml_namespace_projection.py','scripts/validation/check_v3000_coordinate_projection.py']
    record=dict(schema='chematic.cml-namespace-projection.v1',baseline=args.baseline,
        scope='bounded molecular CML -> local CJSON; no XML-schema/stereo/metadata or speed claim',
        head=subprocess.check_output(['git','rev-parse','HEAD'],text=True,cwd=ROOT).strip(),
        source_files=[] if args.baseline else [identity(ROOT/path) for path in paths],
        runner=identity(__file__),fixture=identity(fixture),python_package=identity(chematic.__file__),python_extension=identity(extension.__file__),
        wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')),
        openbabel=identity(args.obabel),openbabel_version=terminal(bounded_run([str(args.obabel.resolve()),'-V'])),
        case_count=len(results),summary=summary,cases=results)
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(summary))
    if not args.baseline:
        for label in ('python','wasm'):
            require(summary[label]['valid_exact']==33 and summary[label]['valid_no_output']==27
                    and summary[label]['invalid_no_output']==84,'candidate failed')
        # 27 reader-profile annotations are now refused by checked conversion;
        # refusals are not counted as successful preservation.


if __name__=='__main__':main()
