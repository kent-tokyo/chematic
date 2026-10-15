#!/usr/bin/env python3
"""Independent, pinned Open Babel #3028 residue-partition diagnostic."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import chematic

ROOT=Path(__file__).resolve().parents[1]
FIXTURE=ROOT/'validation/openbabel_issues/mmcif/branched-residues.cif'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def pdb_records(text):
    return [s for s in text.splitlines() if s.startswith(('ATOM  ','HETATM'))]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--obabel',default='/opt/homebrew/bin/obabel')
    parser.add_argument('--output',type=Path,default=ROOT/'validation/results/openbabel-mmcif-residues-2026-10-14.json')
    args=parser.parse_args()
    version=subprocess.run([args.obabel,'-V'],capture_output=True,text=True,check=True,timeout=10).stdout.strip()
    if not version.startswith('Open Babel 3.2.1 '):raise RuntimeError(f'expected Open Babel 3.2.1, got {version}')
    text=FIXTURE.read_text();header=text.split('HETATM')[0]
    rows=[s.split() for s in text.splitlines() if s.startswith('HETATM')]
    base=chematic.parse_mmcif_residues(text)
    expected_names=['NAG','NAG','NAG','NAG','BMA','BMA']
    assert [a['res_name'] for a in base['atoms']]==expected_names
    checked=[]
    for order in itertools.permutations(range(6)):
        shuffled=header+'\n'.join(' '.join(rows[i]) for i in order)+'\n#\n'
        doc=chematic.parse_mmcif_residues(shuffled)
        assert len(doc['residues'])==3
        for i,old in enumerate(order):
            r=doc['residues'][doc['atom_to_residue'][i]]
            assert r['auth_seq_id']==str(old//2+1)
            assert r['label_seq_id'] is None
            assert r['label_asym_id']=='E' and r['auth_asym_id']=='A'
            assert r['label_comp_id']==r['auth_comp_id']==expected_names[old]
            assert r['atom_indices']==[j for j,v in enumerate(order) if v//2==old//2]
            assert doc['atoms'][i]['x']==float(rows[old][10])
        checked.append({'new_to_old':order,'atom_to_residue':doc['atom_to_residue'],'residue_identity':[{'auth_seq_id':r['auth_seq_id'],'label_comp_id':r['label_comp_id'],'atom_indices':r['atom_indices']} for r in doc['residues']]})
    proc=subprocess.run([args.obabel,'-immcif',str(FIXTURE),'-opdb'],capture_output=True,text=True,check=True,timeout=10)
    output=pdb_records(proc.stdout)
    assert len(output)==6
    names=[r[17:20].strip() for r in output]
    assert names==['NAG']*6, 'upstream behavior changed; investigate before rebaseline'
    assert [float(r[30:38]) for r in output]==[20,22,40,42,60,62]
    control_rows=[r.copy() for r in rows]
    for r in control_rows:r[8]=r[15]
    control=header+'\n'.join(' '.join(r) for r in control_rows)+'\n#\n'
    with tempfile.TemporaryDirectory(prefix='chematic-mmcif-control-') as tmp:
        path=Path(tmp)/'label-sequence.cif';path.write_text(control)
        control_run=subprocess.run([args.obabel,'-immcif',str(path),'-opdb'],capture_output=True,text=True,check=True,timeout=10)
    control_output=pdb_records(control_run.stdout)
    assert [r[17:20].strip() for r in control_output]==expected_names
    control_native=chematic.parse_mmcif_residues(control)
    assert len(control_native['residues'])==3
    source_paths=['crates/chematic-mol/src/mmcif/residues.rs','crates/chematic-py/src/formats/mmcif_residues.rs']
    record={
        'schema_version':1,'status':'local-source-candidate','issue':'https://github.com/openbabel/openbabel/issues/3028',
        'base':'cac480fb source snapshot plus preserved RDKit #9675/#9629 and Open Babel #2845 candidates; patch isolates this reader',
        'scope':'atom_site residue partition; no PDB writer, glycan connectivity, bond perception or real-protein frequency claim',
        'runtime':{'python':sys.version.split()[0],'openbabel':version,'chematic_extension_sha256':digest(Path(chematic.__file__).resolve().read_bytes())},
        'source_sha256':{s:digest((ROOT/s).read_bytes()) for s in source_paths},
        'fixture':{'id':'three-synthetic-nonpolymer-residues','path':str(FIXTURE.relative_to(ROOT)),'sha256':digest(text.encode()),
            'provenance':{'source':'original synthetic atom_site rows from wwPDB identity definitions','retrieved':'2026-10-14','issue':'https://github.com/openbabel/openbabel/issues/3028','transformations':'two atom names per residue; NAG/NAG/BMA with absent label_seq_id and distinct numeric auth_seq_id; no downloaded structure or attachment'},
            'license':'MIT OR Apache-2.0','redistribution':{'status':'allowed','rationale':'original project fixture'},
            'minimization':'six rows isolate repeated atom names, repeated component names, missing polymer sequence, two chain namespaces and one renamed component',
            'oracle':{'name':'source-authored residue identifiers / independent fixed-column PDB decoder','version':'wwPDB PDBx v5 identity definitions; Open Babel 3.2.1 external diagnostic','settings':'all 6! source-row permutations; external conversion defaults'},
            'expected':{'residue_count':3,'component_names':['NAG','NAG','BMA'],'auth_seq_id':['1','2','3'],'atom_counts':[2,2,2],'coordinate_unit':'angstrom','count_tolerance':0},
        },
        'results':{'permutations':len(checked),'partition_mismatches':0,'identity_mismatches':0,'cases_sha256':digest(json.dumps(checked,sort_keys=True,separators=(',',':')).encode()),'native_residues':base['residues']},
        'openbabel':{'command':['obabel','-immcif','<fixture>','-opdb'],'output_names':names,'output_residue_keys':sorted(set((r[21],r[22:26].strip(),r[17:20].strip()) for r in output)),'renamed_BMA_atoms':2,'stderr':proc.stderr.strip(),'output_sha256':digest(proc.stdout.encode()),'numeric_label_control':{'source_sha256':digest(control.encode()),'output_names':[r[17:20].strip() for r in control_output],'stderr':control_run.stderr.strip(),'output_sha256':digest(control_run.stdout.encode())}},
        'limits':'new strict grouping rejects absent chain/sequence/component/atom identity, duplicate atoms and incomplete rows; inherited atom reader supports integer auth_seq_id only',
    }
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'permutations':len(checked),'partition_mismatches':0,'native_residue_count':3,'upstream_residue_count':1,'control_passed':True}))


if __name__=='__main__':main()
