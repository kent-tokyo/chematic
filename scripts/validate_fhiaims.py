#!/usr/bin/env python3
"""Original FHI-aims geometry recipes against Open Babel and ASE oracles."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import tempfile

import ase
from ase.io import read
from ase.constraints import FixAtoms, FixCartesian
import chematic
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
FIELDS=('coordinate_kinds','lattice_angstrom','constraints','initial_charges','initial_moments','lattice_constraints')


def digest(data):return hashlib.sha256(data).hexdigest()


def emit(g):return chematic.write_fhiaims(g['elements'],g['coordinates'],**{k:g[k] for k in FIELDS})


def recipe(seed,dimension,kind,metadata,crlf):
    n=seed%4+1
    symbols=['Si','Cl','Fe','H'][:n]
    coordinates=[[seed/7+i/3,-i/5+0.125,(-1)**i*(seed+1)/11] for i in range(n)]
    lattice=[[4+seed/8,1,0],[-1,5,2],[1,0,(-1)**seed*6]][:dimension]
    kinds=[kind if kind!='mixed' else ('fractional' if i%2==0 else 'cartesian') for i in range(n)]
    mask=[bool(seed&1),bool(seed&2),bool(seed&4)]
    g={'elements':symbols,'coordinates':coordinates,'coordinate_kinds':kinds,'lattice_angstrom':lattice,
       'constraints':[mask[:] if metadata else [False]*3 for _ in symbols],
       'lattice_constraints':[[False,True,False] if metadata else [False]*3 for _ in lattice],
       'initial_charges':[0.25 if metadata and i%2==0 else None for i in range(n)],
       'initial_moments':[(-1)**i*1.5 if metadata else None for i in range(n)]}
    text=emit(g)
    if crlf:text=text.replace('\n','\r\n')
    return g,text


def native_check(expected,text):
    g=chematic.parse_fhiaims(text)
    for key,value in expected.items():assert g[key]==value,key
    assert g['source_text']==text
    back=chematic.parse_fhiaims(emit(g))
    for key in g.keys()-{'source_text'}:assert back[key]==g[key],key
    points=[]
    for kind,row in zip(g['coordinate_kinds'],g['coordinates']):
        points.append((np.array(row)@np.array(g['lattice_angstrom'])).tolist() if kind=='fractional' else row)
    np.testing.assert_allclose(g['positions_angstrom'],points,rtol=0,atol=2e-12)
    return g


def ase_check(g,path):
    projected=dict(g)
    projected['lattice_constraints']=[[False]*3 for _ in g['lattice_angstrom']]
    if len(set(g['coordinate_kinds']))>1:
        projected['coordinates']=g['positions_angstrom']
        projected['coordinate_kinds']=['cartesian']*len(g['elements'])
    path.write_text(emit(projected))
    a=read(str(path),format='aims')
    assert a.get_chemical_symbols()==g['elements']
    np.testing.assert_allclose(a.positions,g['positions_angstrom'],rtol=0,atol=2e-12)
    np.testing.assert_allclose(a.get_initial_charges(),[0 if x is None else x for x in g['initial_charges']],rtol=0,atol=0)
    np.testing.assert_allclose(a.get_initial_magnetic_moments(),[0 if x is None else x for x in g['initial_moments']],rtol=0,atol=0)
    if g['periodicity']==3:np.testing.assert_allclose(a.cell.array,g['lattice_angstrom'],rtol=0,atol=2e-12)
    mask=np.zeros((len(g['elements']),3),dtype=bool)
    for constraint in a.constraints:
        if isinstance(constraint,FixAtoms):mask[constraint.get_indices()]=True
        elif isinstance(constraint,FixCartesian):mask[constraint.get_indices()]=constraint.mask
        else:raise AssertionError(type(constraint))
    assert mask.tolist()==g['constraints']


def obabel(executable,path,fmt):
    r=subprocess.run([executable,'-ifhiaims',str(path),'-o'+fmt,'-ab'],capture_output=True,text=True,check=True,timeout=10)
    assert '1 molecule converted' in r.stderr,r.stderr
    return r.stdout


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--obabel',default='/opt/homebrew/bin/obabel')
    p.add_argument('--output',type=Path,default=ROOT/'validation/results/openbabel-fhiaims-2026-10-15.json')
    args=p.parse_args()
    version=subprocess.check_output([args.obabel,'-V'],text=True).strip()
    assert version.split()[2]=='3.2.1' and ase.__version__=='3.27.0'
    assert np.__version__=='2.5.4'
    manifest_path=ROOT/'validation/openbabel_issues/fhiaims/contract-fixtures.json'
    manifest=json.loads(manifest_path.read_text())
    for fixture in manifest['records']:
        assert digest(fixture['payload'].encode())==fixture['sha256']
        actual=chematic.parse_fhiaims(fixture['payload'])
        for key,value in fixture['expected'].items():
            if key not in ('units','tolerance'):assert actual[key]==value
    records=[];ase_cases=ob_cases=fractional_cases=0;max_ob_error=max_cell_error=0.;metadata_dropped=0
    layouts=[(0,'cartesian'),(1,'cartesian'),(2,'cartesian'),(3,'cartesian'),(3,'fractional'),(3,'mixed')]
    with tempfile.TemporaryDirectory(prefix='chematic-fhiaims-oracle-') as directory:
        path=Path(directory)/'geometry.in'
        for seed,(dimension,kind),metadata,crlf in itertools.product(range(8),layouts,[False,True],[False,True]):
            expected,text=recipe(seed,dimension,kind,metadata,crlf)
            g=native_check(expected,text)
            fractional_cases+=any(k=='fractional' for k in g['coordinate_kinds'])
            if dimension in (0,3):ase_check(g,path);ase_cases+=1
            if dimension in (0,3) and kind=='cartesian':
                path.write_bytes(text.encode());doc=json.loads(obabel(args.obabel,path,'cjson'))
                assert doc['atoms']['elements']['number']==[{'Si':14,'Cl':17,'Fe':26,'H':1}[s] for s in g['elements']]
                np.testing.assert_allclose(np.array(doc['atoms']['coords']['3d']).reshape(-1,3),g['positions_angstrom'],rtol=0,atol=2e-12)
                if dimension:
                    error=float(np.max(np.abs(np.array(doc['unitCell']['cellVectors']).reshape(3,3)-g['lattice_angstrom'])))
                    max_cell_error=max(max_cell_error,error);assert error<2e-12
                back=chematic.parse_fhiaims(obabel(args.obabel,path,'fhiaims'))
                error=float(np.max(np.abs(np.array(back['positions_angstrom'])-g['positions_angstrom'])))
                if dimension:error=max(error,float(np.max(np.abs(np.array(back['lattice_angstrom'])-g['lattice_angstrom']))))
                max_ob_error=max(max_ob_error,error);assert error<=5.00001e-6
                assert all(v is None for v in back['initial_charges']+back['initial_moments'])
                assert not any(itertools.chain.from_iterable(back['constraints']+back['lattice_constraints']))
                metadata_dropped+=sum(x is not None for x in g['initial_charges']+g['initial_moments'])
                ob_cases+=1
            records.append({'id':f's{seed}-d{dimension}-{kind}-metadata{int(metadata)}-crlf{int(crlf)}',
                'provenance':'Original synthetic recipe() in scripts/validate_fhiaims.py; authored 2026-10-15',
                'license':'MIT OR Apache-2.0','redistribution':'allowed; original project data',
                'minimization':'1-4 canonical elements; coordinate frame, dimensionality, metadata and CRLF varied independently',
                'oracle':'ASE 3.27.0, Open Babel 3.2.1, NumPy 2.5.4 as scoped by this gate',
                'expected':{'structured_roundtrip':'exact','Cartesian_tolerance_angstrom':2e-12,'OpenBabel_output_tolerance_angstrom':5.00001e-6},
                'sha256':digest(text.replace('\r\n','\n').encode()),'status':'passed'})
        diagnostic='lattice_vector 4 1 0\nlattice_vector -1 5 2\nlattice_vector 1 0 -6\natom_frac 0.5 0.25 -0.5 Si\n'
        path.write_text(diagnostic);fractional_ob=json.loads(obabel(args.obabel,path,'cjson'))['atoms']['coords']['3d']
        assert fractional_ob==[0.5,0.25,-0.5]
        assert chematic.parse_fhiaims(diagnostic)['positions_angstrom']==[[1.25,1.75,3.5]]
        path.write_text('lattice_vector 4 1 0\natom 1 2 3 C\n')
        assert 'unitCell' not in json.loads(obabel(args.obabel,path,'cjson'))
    report={'schema_version':1,'status':'local-source-candidate; uncommitted and unpushed',
        'contract_fixture_manifest':str(manifest_path.relative_to(ROOT)),
        'oracles':{'openbabel':version,'ase':ase.__version__,'numpy':np.__version__,'FHIaims_manual':'261002'},
        'results':{'structured_roundtrips':len(records),'ase_cases':ase_cases,'openbabel_cartesian_cases':ob_cases,
            'fractional_or_mixed_cases':fractional_cases,'max_openbabel_cell_error_angstrom':max_cell_error,
            'max_openbabel_written_error_angstrom':max_ob_error,'initial_density_fields_dropped_by_openbabel':metadata_dropped,'unexpected_mismatches':0},
        'records':records,'normalized_records_sha256':digest(json.dumps(records,sort_keys=True,separators=(',',':')).encode()),
        'scope':['Open Babel comparison only for Cartesian 0D/3D geometry','ASE checks homogeneous input or Cartesian projection for mixed frames; lattice constraints omitted from ASE input',
            'Partial periodicity and lattice constraints checked by exact internal roundtrip','Fractional Cartesian projection independently checked with NumPy matrix multiplication',
            'No FHI-aims calculation executable or control.in validation'],
        'observed_openbabel_deltas':{'fractional_input_reported_as_cartesian':fractional_ob,'native_and_ase_expected':[1.25,1.75,3.5],
            'partial_lattice':'Tested 1D cell is omitted from CJSON','metadata':'Writer omits atom/lattice constraints and initial charge/moment'}}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['results'],indent=2))


if __name__=='__main__':main()
