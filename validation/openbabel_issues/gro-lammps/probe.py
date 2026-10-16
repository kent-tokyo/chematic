"""Original Open Babel box-loss probe and independent geometry checks."""
import argparse
import hashlib
import itertools
import json
import re
import subprocess
from pathlib import Path


def openbabel():
    fixture=Path(__file__).with_name('orthogonal.gro')
    result=subprocess.run(['obabel','-igro',str(fixture),'-olmpdat'],capture_output=True,text=True,check=True)
    bounds={axis:[float(a),float(b)] for a,b,axis in re.findall(r'([-+\d.eE]+)\s+([-+\d.eE]+)\s+([xyz])lo\s+[xyz]hi',result.stdout)}
    assert bounds=={'x':[-.5,9.5],'y':[-.5,9.5],'z':[-.5,9.5]}
    return {'version':subprocess.check_output(['obabel','-V'],text=True).strip(),
            'actual_bounds_angstrom':bounds,'expected_bounds_angstrom':{'x':[0,20],'y':[0,30],'z':[0,40]},
            'input_sha256':hashlib.sha256(fixture.read_bytes()).hexdigest()}


def compare():
    import numpy as np
    import chematic
    cells=[np.diag([2.,3.,4.]),np.diag([1.5,2.25,3.75]),
           np.array([[2,0,0],[.5,3,0],[-.25,.75,4]]),
           np.array([[2.5,0,0],[-.75,3.5,0],[.125,-.5,4.25]])]
    fractions=np.array([[.125,.25,.5],[-.25,1.25,.5],[.75,.5,.25]])
    velocity=np.array([[.125,-.25,.5],[.25,.5,-.75],[-.5,.25,.125]])
    cases=[]
    for cell_index,cell in enumerate(cells):
        for units,vel_factor in [('real',.01),('metal',10.)]:
            for has_velocity in [False,True]:
                for order in itertools.permutations(range(3)):
                    positions=fractions@cell
                    rows=[]
                    for original in order:
                        row=f'{7:5d}{"ION":<5}{["CA","NA","XX"][original]:>5}{[99999,0,0][original]:5d}'
                        row+=''.join(f'{v:11.6f}' for v in positions[original])
                        if has_velocity:row+=''.join(f'{v:11.7f}' for v in velocity[original])
                        rows.append(row)
                    a,b,c=cell
                    box=[a[0],b[1],c[2],a[1],a[2],b[0],b[2],c[0],c[1]]
                    text='original analytic geometry\n3\n'+'\n'.join(rows)+'\n'+' '.join(map(str,box))+'\n'
                    result=chematic.gro_to_lammps_atomic(text,[[1,2,1][i] for i in order],[12.,14.],units=units)
                    data=chematic.parse_lammps_data(result['data_text'],'atomic')
                    out=data['box'];tilt=out['tilt'] or [0.,0.,0.]
                    transformed=np.array([[out['hi'][0],0,0],[tilt[0],out['hi'][1],0],[tilt[1],tilt[2],out['hi'][2]]])
                    assert np.allclose(transformed,cell*10,rtol=0,atol=1e-12)
                    actual=np.array([[atom[k] for k in ['x','y','z']] for atom in data['atoms']])
                    recovered=np.linalg.solve(transformed.T,actual.T).T
                    assert np.allclose(recovered,fractions[list(order)],rtol=0,atol=1e-12)
                    assert np.isclose(np.linalg.det(transformed)/np.linalg.det(cell),1000.,rtol=1e-12)
                    if has_velocity:
                        actual_velocity=np.array([[v[k] for k in ['vx','vy','vz']] for v in data['velocities']])
                        assert np.allclose(actual_velocity,velocity[list(order)]*vel_factor,rtol=0,atol=1e-12)
                    else:assert data['velocities']==[]
                    assert result['source_text']==text
                    assert [i['atom_number'] for i in result['atom_identity']]==[[99999,0,0][i] for i in order]
                    assert [a['atom_type'] for a in data['atoms']]==[[1,2,1][i] for i in order]
                    cases.append({'cell':cell_index,'units':units,'velocities':has_velocity,'order':order,
                                  'input_sha256':hashlib.sha256(text.encode()).hexdigest(),
                                  'output_sha256':hashlib.sha256(result['data_text'].encode()).hexdigest()})
    module=Path(__import__('chematic.chematic',fromlist=['__file__']).__file__)
    return {'passed':len(cases),'cases':cases,'module_sha256':hashlib.sha256(module.read_bytes()).hexdigest()}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--compare',action='store_true');args=parser.parse_args()
    result={'openbabel':openbabel()}
    if args.compare:result['analytic_geometry']=compare()
    print(json.dumps(result,indent=2))
