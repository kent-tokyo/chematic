#!/usr/bin/env python3
"""Paired fallback-off/on check on previously exposed B/Se development cases."""
import argparse
import hashlib
import json
from pathlib import Path
import chematic
from rdkit import Chem, rdBase
from cod_robust_scoring import annotate
from platinum_development_probe import graph_matches
from platinum_geometry import prepare_geometry_checks
from summarize_cod_candidate_v2 import recovered


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument('--sdf',required=True);p.add_argument('--selected',required=True);p.add_argument('--output',required=True);args=p.parse_args()
    source=json.loads(Path(args.selected).read_text())
    if sha(args.sdf)!=source['source_sha256']: raise ValueError('SDF mismatch')
    supplier=Chem.SDMolSupplier(args.sdf,removeHs=False)
    out={'source_sha256':sha(args.sdf),'selection_sha256':sha(args.selected),'script_sha256':sha(__file__),'api_sha256':sha(chematic._timed_conformers.__file__),'rdkit':rdBase.rdkitVersion,'count':20,'seconds':30,'rows':[]}
    for row in source['rows']:
        entry={'index':row['index'],'name':row['name'],'engines':{}}
        try:
            ref=supplier[row['index']]
            if ref is None or ref.GetProp('_Name')!=row['name']: raise ValueError('identity mismatch')
            ref=Chem.RemoveHs(ref); smiles=Chem.MolToSmiles(ref,isomericSmiles=True); mol=Chem.AddHs(Chem.MolFromSmiles(smiles)); other=chematic.from_smiles(smiles).add_hydrogens()
            if not graph_matches(other,mol): raise ValueError('indexed graph mismatch')
            checks=prepare_geometry_checks(mol)
            entry['geometry_reference']={'bonds_supported':len(checks[0]),'bonds_total':mol.GetNumBonds(),'angles_supported':len(checks[1]),'angles_total':sum(a.GetDegree()*(a.GetDegree()-1)//2 for a in mol.GetAtoms())}
            for name,enabled in [('off',False),('on',True)]:
                v=chematic.timed_conformer_batch([smiles],seconds=30,num_confs=20,repair_mmff_iters=400,allow_unoptimized=enabled)[0]
                item={'status':v['status'],'error':v['error'],'fallback_reason':v['fallback_reason'],'elapsed':v['elapsed'],'candidates':[(status,coords) for status,_,coords in v['conformers']]}
                entry['engines'][name]=annotate(item,mol,ref,checks)
        except Exception as exc:
            entry['input_error']=repr(exc)
        out['rows'].append(entry)
        Path(args.output).write_text(json.dumps(out,allow_nan=False)+'\n')
        print(len(out['rows']),'/',len(source['rows']),entry['name'],flush=True)
    transitions={'recovered':[],'lost':[],'unchanged_failure':[],'unchanged_success':[]}
    for row in out['rows']:
        a=recovered(row,'off');b=recovered(row,'on')
        key='recovered' if b and not a else 'lost' if a and not b else 'unchanged_success' if a else 'unchanged_failure'
        transitions[key].append(row['name'])
    out['transitions']=transitions
    Path(args.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print({k:len(v) for k,v in transitions.items()},flush=True)

if __name__=='__main__':main()
