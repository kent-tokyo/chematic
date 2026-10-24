#!/usr/bin/env python3
"""Native/OB 3D V3000 regeneration with exact coordinate and RDKit observers.

The authored numeric tokens supply expected bits. RDKit independently reads
source and output indexed graphs/stereo; no coordinate tolerance or speed score.
"""
from __future__ import annotations
import argparse
import gzip
import json
from pathlib import Path
import struct
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run, identity, require, strict_json


def bits(value):
    return struct.pack('>d',float(value)).hex()


def atom_rows(text):
    rows=[]
    pending=''
    active=False
    for line in text.splitlines():
        if not line.startswith('M  V30 '):
            continue
        payload=line[7:]
        pending+=payload[:-1] if payload.endswith('-') else payload
        if payload.endswith('-'):
            continue
        if pending=='END ATOM': active=False
        elif active: rows.append(pending.split())
        if pending=='BEGIN ATOM': active=True
        pending=''
    require(not pending and rows,'missing or incomplete V3000 atoms')
    require(all(len(row)>=6 for row in rows),'short V3000 atom row')
    return rows


def raw_coordinates(text):
    return [[bits(v) for v in row[2:5]] for row in atom_rows(text)]


def external_snapshot(text):
    from rdkit import Chem
    mol=Chem.MolFromMolBlock(text,sanitize=False,removeHs=False,strictParsing=True)
    require(mol is not None and mol.GetNumConformers()==1,'RDKit cannot read V3000')
    Chem.AssignStereochemistry(mol,cleanIt=True,force=True)
    conf=mol.GetConformer()
    return dict(coordinates=[[bits(v) for v in conf.GetAtomPosition(i)] for i in range(mol.GetNumAtoms())],
                graph=dict(atoms=[[a.GetAtomicNum(),a.GetFormalCharge(),a.GetIsotope(),a.GetAtomMapNum(),a.GetIsAromatic()] for a in mol.GetAtoms()],
                           bonds=sorted([*sorted((b.GetBeginAtomIdx(),b.GetEndAtomIdx())),b.GetBondTypeAsDouble()] for b in mol.GetBonds())),
                stereo=[a.GetProp('_CIPCode') if a.HasProp('_CIPCode') else None for a in mol.GetAtoms()])


def observe(expected,reference,text,reader=external_snapshot):
    result=dict(readable=False,coordinates_exact=False,indexed_graph_exact=False,source_stereo_exact=False)
    try:
        actual=reader(text)
        result.update(readable=True,coordinates_exact=actual['coordinates']==expected,
                      indexed_graph_exact=actual['graph']==reference['graph'],source_stereo_exact=actual['stereo']==reference['stereo'],actual=actual)
    except (ValueError,RuntimeError,KeyError,IndexError,OverflowError) as error:
        result['error']=str(error)
    result['joint_preserved']=all(result[k] for k in ('readable','coordinates_exact','indexed_graph_exact','source_stereo_exact'))
    return result


def terminal(run):
    return {k:v.decode('utf8',errors='replace') if isinstance(v,bytes) else v for k,v in run.items()}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native-executable',type=Path,required=True)
    p.add_argument('--probe-source',type=Path,default=ROOT/'crates/chematic-mol/examples/v3000_3d_numeric_dump.rs')
    p.add_argument('--source-revision')
    p.add_argument('--baseline',action='store_true')
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    from rdkit import rdBase
    from rdkit.Chem import rdmolfiles
    packet_path=ROOT/'validation/v3000_3d_numeric/cases.json'
    packet=strict_json(packet_path.read_text())
    command=[str(args.native_executable.resolve()),str(packet_path)]
    run=terminal(bounded_run(command));run['command']=command
    require(run['resource_status']=='completed' and run['exit_code']==0,'native probe failed')
    native=strict_json(run['stdout'])['rows']
    require(len(native)==len(packet['cases']),'native case count mismatch')
    rows=[]
    for case,output in zip(packet['cases'],native):
        require(case['id']==output['id'],'native order mismatch')
        expected=[[bits(v) for v in row] for row in case['expected_coordinates']]
        require(raw_coordinates(case['input'])==expected,'authored tokens disagree')
        row=dict(id=case['id'],scope=case['scope'],input=case['input'],expected_coordinates=expected,native_output=output,
                 native_raw_exact={k:raw_coordinates(output[k])==expected for k in ('plain','checked')},
                 numeric_packet_lines_within_80={k:all(len(line.encode())<=80 for line in output[k].splitlines() if line.startswith('M  V30 ')) for k in ('plain','checked')})
        if case['scope']=='paired':
            reference=external_snapshot(case['input'])
            require(reference['coordinates']==expected,'source external reader does not retain authored coordinates')
            row['source_rdkit']=reference
            row['native_observation']={k:observe(expected,reference,output[k]) for k in ('plain','checked')}
            ob_command=['/opt/homebrew/bin/obabel','-imol','-omol','-x3']
            ob=terminal(bounded_run(ob_command,case['input'].encode()));ob['command']=ob_command
            row['openbabel_run']=ob
            text=ob['stdout'] if ob['resource_status']=='completed' and ob['exit_code']==0 else ''
            row['openbabel_observation']=observe(expected,reference,text)
            try: row['openbabel_raw_exact']=raw_coordinates(text)==expected
            except (ValueError,IndexError): row['openbabel_raw_exact']=False
        rows.append(row)
    paired=[r for r in rows if r['scope']=='paired']
    stereo_controls=[r['source_rdkit']['stereo'] for r in paired if r['id'].startswith('tetra_')]
    require(len(stereo_controls)==2 and all(sum(code is not None for code in codes)==1 for codes in stereo_controls)
            and sorted(codes[0] for codes in stereo_controls)==['R','S'],'mirrored source CIP controls failed qualification')
    metrics=('readable','coordinates_exact','indexed_graph_exact','source_stereo_exact','joint_preserved')
    summary={k:{m:sum((r['openbabel_observation'] if k=='openbabel' else r['native_observation'][k])[m] for r in paired) for m in metrics} for k in ('plain','checked','openbabel')}
    production='crates/chematic-mol/src/mol3000/writer.rs'
    data=subprocess.check_output(['git','show',f'{args.source_revision}:{production}'],cwd=ROOT) if args.source_revision else None
    result=dict(schema='chematic.v3000-3d-numeric-observation.v1',baseline=args.baseline,head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                source_revision=args.source_revision or 'working_tree',production_source=identity(ROOT/production,data),
                packet=identity(packet_path),harness=identity(__file__),native_executable=identity(args.native_executable),probe_source=identity(args.probe_source),native_run=run,
                openbabel=identity('/opt/homebrew/bin/obabel'),openbabel_version=terminal(bounded_run(['/opt/homebrew/bin/obabel','-V'])),
                rdkit_version=rdBase.rdkitVersion,rdkit_parser=identity(rdmolfiles.__file__),rdkit_options=dict(sanitize=False,removeHs=False,strictParsing=True),
                paired_case_count=len(paired),native_only_case_count=len(rows)-len(paired),summary=summary,cases=rows)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_bytes(gzip.compress((json.dumps(result,indent=2)+'\n').encode(),mtime=0))
    print(json.dumps(summary))
    if not args.baseline:
        require(all(all(r['native_raw_exact'].values()) and all(r['numeric_packet_lines_within_80'].values()) for r in rows),'native bits or packet line width failed')
        require(all(r['native_observation'][k]['joint_preserved'] for r in paired for k in ('plain','checked')),'external contract failed')


if __name__=='__main__':main()
