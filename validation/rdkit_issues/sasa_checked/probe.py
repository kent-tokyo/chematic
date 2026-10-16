"""Analytic sphere controls, explicit RDKit FreeSASA diagnostics, no backend parity claim."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import resource
import subprocess
import sys

p=argparse.ArgumentParser()
p.add_argument('--output',type=Path)
p.add_argument('--short-radii',action='store_true')
a=p.parse_args()
if a.short_radii:
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    from rdkit import Chem,rdBase
    from rdkit.Chem import rdFreeSASA as F
    mol=Chem.MolFromSmiles('CO')
    conf=Chem.Conformer(2);conf.SetAtomPosition(0,(0,0,0));conf.SetAtomPosition(1,(0,1,0));mol.AddConformer(conf)
    print(json.dumps({'rdkit_version':rdBase.rdkitVersion}),flush=True)
    try:
        print(json.dumps({'returned':F.CalcSASA(mol,[1.0])}),flush=True)
    except Exception as exc:
        print(json.dumps({'exception':type(exc).__name__,'message':str(exc)}),flush=True)
    sys.exit(0)

import chematic as c
from rdkit import Chem,rdBase
from rdkit.Chem import rdFreeSASA as F


def analytic(r,s,d):
    full=[4*math.pi*r*r,4*math.pi*s*s]
    if d>=r+s:return full
    if d<=abs(r-s):return [full[0] if r>s else 0.0,full[1] if s>r else 0.0]
    x=(d*d+r*r-s*s)/(2*d)
    y=d-x
    return [2*math.pi*r*(r+x),2*math.pi*s*(s+y)]

mol=c.from_smiles('CO')
rows=[]
for radii,probe,direction,d in itertools.product(
    [(1.,1.),(1.,1.5),(1.2,1.8)], [0.,1.4],
    [(1.,0.,0.),(0.,1.,0.),(0.,0.,1.),(.6,.8,0.)], [.4,1.,2.,4.]):
    xyz=[[0.,0.,0.],[d*v for v in direction]]
    inflated=[v+probe for v in radii]
    expected=analytic(*inflated,d)
    full=[4*math.pi*r*r for r in inflated]
    native=mol.sasa_with_radii_checked(xyz,radii,probe_radius=probe,sphere_points=8192)
    err=[abs(x-y)/f for x,y,f in zip(native['per_atom'],expected,full)]
    assert max(err)<=.005,(radii,probe,direction,d,err)
    assert native['total']==math.fsum(native['per_atom']) or abs(native['total']-math.fsum(native['per_atom']))<1e-12
    permuted=mol.renumber_atoms([1,0]).sasa_with_radii_checked(xyz[::-1],radii[::-1],probe_radius=probe,sphere_points=8192)
    assert permuted['per_atom']==native['per_atom'][::-1]
    translated=mol.sasa_with_radii_checked([[v+128 for v in row] for row in xyz],radii,probe_radius=probe,sphere_points=8192)
    # Translation rounds nonbinary input decimals, so compare sampled area
    # within one sample per sphere rather than assuming exact float equality.
    assert all(abs(x-y)<=f/8192+1e-12 for x,y,f in zip(translated['per_atom'],native['per_atom'],full))
    rm=Chem.MolFromSmiles('CO');conf=Chem.Conformer(2)
    for i,row in enumerate(xyz):conf.SetAtomPosition(i,row)
    rm.AddConformer(conf)
    opts=F.SASAOpts(F.LeeRichards,F.Protor,probe)
    oracle=F.CalcSASA(rm,list(radii),opts=opts)
    oracle_atom=[atom.GetDoubleProp('SASA') for atom in rm.GetAtoms()]
    oracle_finite=all(math.isfinite(x) for x in [oracle,*oracle_atom])
    oracle_err=[abs(x-y)/f if math.isfinite(x) else None for x,y,f in zip(oracle_atom,expected,full)]
    # FreeSASA is a diagnostic, not the correctness oracle: some valid
    # geometries return NaN at this default Lee-Richards resolution.
    oracle_status=('nonfinite' if not oracle_finite else
                   'outside_diagnostic_tolerance' if max(oracle_err)>.025 else 'within_diagnostic_tolerance')
    rows.append({'radii':radii,'probe_radius':probe,'direction':direction,'distance':d,
                 'analytic':expected,'native':native,'native_relative_full_sphere_error':err,
                 'rdkit_total':oracle if math.isfinite(oracle) else None,
                 'rdkit_per_atom':[x if math.isfinite(x) else None for x in oracle_atom],
                 'rdkit_relative_full_sphere_error':oracle_err,'rdkit_status':oracle_status})

extra=[]
for scale in [1e-100,1e-50,1.,1e50,1e100]:
    native=mol.sasa_with_radii_checked([[0,0,0],[0,scale,0]],[scale,scale],probe_radius=0,sphere_points=4096)
    ratio=native['total']/(scale*scale)
    assert abs(ratio-6*math.pi)<1e-12
    extra.append({'scale':scale,'total':native['total'],'normalized_total':ratio})

short=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--short-radii'],capture_output=True,text=True,timeout=15)
import chematic.chematic as ext
result={'rdkit_version':rdBase.rdkitVersion,'rdkit_algorithm':'LeeRichards (FreeSASA default resolution)',
        'native_algorithm':'midpoint Fibonacci Shrake-Rupley, 8192 points',
        'native_path':ext.__file__,'native_sha256':hashlib.sha256(Path(ext.__file__).read_bytes()).hexdigest(),
        'analytic_geometry_cases':len(rows),'native_geometry_runs':len(rows)*3,
        'max_native_error_relative_to_full_sphere':max(max(r['native_relative_full_sphere_error']) for r in rows),
        'max_finite_rdkit_error_relative_to_full_sphere':max(e for r in rows for e in r['rdkit_relative_full_sphere_error'] if e is not None),
        'rdkit_diagnostic_counts':{status:sum(r['rdkit_status']==status for r in rows) for status in ['within_diagnostic_tolerance','outside_diagnostic_tolerance','nonfinite']},
        'native_tolerance_relative_to_full_sphere':.005,'rdkit_diagnostic_tolerance_relative_to_full_sphere':.025,
        'extreme_scale_cases':extra,'short_radii_subprocess':{'returncode':short.returncode,'stdout':short.stdout,'stderr':short.stderr},
        'native_analytic_failures':0,'rows':rows}
assert a.output is not None
a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
