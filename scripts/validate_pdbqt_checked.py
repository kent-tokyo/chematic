#!/usr/bin/env python3
"""Reproduce #2101; check fixed columns, real WASM parity and OB cross-read.

Use a Python package built from the candidate via PYTHONPATH. The probe is the
original C++ program next to the XYZ fixture, linked to installed Open Babel.
This checks serialization only, not atom typing quality or docking validity.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import platform
import subprocess
from pathlib import Path
import chematic

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "validation/openbabel_issues/pdbqt-checked/coordinate-overflow.xyz"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fixtures():
    cases = []
    # Element order is authored independently of the writer/reader. Explicit H
    # remains a graph atom; implicit H is intentionally not materialized.
    profiles = [("C", [6]), ("CCO", [6, 6, 8]), ("c1ccccc1", [6]*6),
                ("[H]O[H]", [1, 8, 1]), ("Cl.Br", [17, 35]),
                ("[Mg+2].[Zn+2]", [12, 30]), ("NPSF", [7, 15, 16, 9]),
                ("[Ca+2].[Mn+2].[Fe+2].I", [20, 25, 26, 53])]
    for profile, (smiles, elements) in enumerate(profiles):
        for translation in [-125.125, 0.0, 9990.125]:
            n = len(elements)
            cases.append(dict(id=f"profile-{profile}-{translation}", smiles=smiles,
                elements=elements, coords=[[translation+i*.71234, -50.23123+i*.52345, i*.31415] for i in range(n)],
                charges=[-.1234+i*.03123 for i in range(n)], name="source order", accepted=True))
    for i, (coords, q) in enumerate([([9999.9994, -999.9994, -.0001], 99.9994),
                                     ([0, 0, 0], -9.9994), ([0, 0, 0], -.0001)]):
        cases.append(dict(id=f"edge-{i}", smiles="C", elements=[6], coords=[coords], charges=[q], name="", accepted=True))
    for axis in range(3):
        for bad in [2261608.729, 10000, -1000, 9999.9999, -999.9999]:
            p=[0,0,0];p[axis]=bad
            cases.append(dict(id=f"coordinate-{axis}-{bad}", smiles="C", coords=[p], charges=[0], name="LIG", accepted=False))
    for q in [100,-10,99.9999,-9.9999]:
        cases.append(dict(id=f"charge-{q}", smiles="C", coords=[[0,0,0]], charges=[q], name="LIG", accepted=False))
    for field in ["coords","charges"]:
        for count in [0,2]:
            c=dict(id=f"count-{field}-{count}", smiles="C", coords=[[0,0,0]], charges=[0], name="LIG", accepted=False)
            c[field]=([[0,0,0]] if field=="coords" else [0])*count
            cases.append(c)
    for i,name in enumerate(["\nATOM","a\rb","\0","\t","\x7f","分子","x"*1025]):
        cases.append(dict(id=f"name-{i}", smiles="C", coords=[[0,0,0]], charges=[0], name=name, accepted=False))
    for smiles in ["*","[B]","[Si]"]:
        cases.append(dict(id=f"type-{smiles}", smiles=smiles, coords=[[0,0,0]], charges=[0], name="LIG", accepted=False))
    return cases


def fixed_columns(text):
    rows=[]
    for line in text.splitlines():
        if not line.startswith("ATOM"):
            continue
        assert len(line)==79, repr(line)
        rows.append(dict(serial=int(line[6:11]), coords=[float(line[a:b]) for a,b in [(30,38),(38,46),(46,54)]],
                         charge=float(line[70:76]), atom_type=line[77:79].strip()))
    return rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wasm-module', required=True, type=Path)
    parser.add_argument('--openbabel-probe', required=True, type=Path)
    parser.add_argument('--obabel', default='/opt/homebrew/bin/obabel')
    parser.add_argument('--output', type=Path, default=ROOT/'validation/results/pdbqt-checked-2026-10-16')
    args=parser.parse_args()
    cases=fixtures()
    node=subprocess.run(['node',str(ROOT/'scripts/validate_pdbqt_checked_node.mjs'),str(args.wasm_module.resolve())],
                        input=json.dumps(cases),text=True,capture_output=True,check=True)
    wasm=json.loads(node.stdout)
    results=[]
    for c,w in zip(cases,wasm,strict=True):
        assert c['id']==w['id']
        m=chematic.from_smiles(c['smiles'])
        try:
            text=chematic.write_pdbqt_checked(m,c['coords'],c['charges'],c['name'])
            py={'output':text}
        except ValueError as e:
            py={'error':str(e)}
        result={'input':c,'python':py,'wasm':w}
        if c['accepted']:
            assert 'output' in py and w['output']==text, c['id']
            assert text==m.to_pdbqt([tuple(p) for p in c['coords']],c['charges'],c['name'])
            decoded=fixed_columns(text)
            assert len(decoded)==len(c['elements'])
            for i,row in enumerate(decoded):
                assert row['serial']==i+1
                assert row['coords']==[float(f'{v:.3f}') for v in c['coords'][i]]
                assert row['charge']==float(f"{c['charges'][i]:.3f}")
            cross=subprocess.run([str(args.openbabel_probe)],input=text,text=True,capture_output=True,check=True)
            data=cross.stdout.splitlines()
            assert int(data[0])==len(decoded)
            for i,line in enumerate(data[1:]):
                values=line.split()
                assert int(values[0])==c['elements'][i]
                assert list(map(float,values[1:4]))==decoded[i]['coords']
                assert float(values[4])==decoded[i]['charge']
            assert len(data)==1+len(decoded)
            result.update(fixed_columns=decoded,openbabel_stdout=cross.stdout,openbabel_stderr=cross.stderr)
        else:
            assert 'error' in py and 'error' in w, c['id']
        results.append(result)
    # The authored XYZ has the reported magnitude pattern, without depending on
    # the issue's attachment or gen3D randomness. Exit success alone is unsafe.
    command=[args.obabel,'-ixyz',str(FIXTURE),'-opdbqt','-xr','-xp','-xh']
    upstream=subprocess.run(command,text=True,capture_output=True,check=True)
    line=next(s for s in upstream.stdout.splitlines() if s.startswith('ATOM'))
    assert '2261608.7291334600.251' in line
    legacy=chematic.from_smiles('N').to_pdbqt([(2261608.729,1334600.251,0)],[0],'probe')
    assert '2261608.7291334600.251' in legacy
    assert len(line)>79
    packet={'schema':'chematic.pdbqt-checked.v1','cases':results,
            'overflow_reproduction':{'command':command,'stdout':upstream.stdout,'stderr':upstream.stderr,'returncode':upstream.returncode,'legacy_chematic':legacy}}
    output=args.output;output.mkdir(parents=True,exist_ok=True)
    raw=gzip.compress(json.dumps(packet,sort_keys=True,ensure_ascii=True).encode(),mtime=0)
    (output/'comparisons.json.gz').write_bytes(raw)
    extension=Path(chematic.chematic.__file__)
    paths=['crates/chematic-mol/src/pdbqt/checked.rs','crates/chematic-mol/src/pdbqt.rs',
           'crates/chematic-py/src/formats/pdbqt_checked.rs','crates/chematic-wasm/src/pdbqt_checked.rs',
           'scripts/validate_pdbqt_checked.py','scripts/validate_pdbqt_checked_node.mjs',
           'validation/openbabel_issues/pdbqt-checked/openbabel_read_probe.cpp']
    summary={'schema':'chematic.pdbqt-checked-summary.v1',
             'source_base':subprocess.check_output(['git','merge-base','HEAD','origin/main'],cwd=ROOT,text=True).strip(),
             'source_sha256':{p:sha(ROOT/p) for p in paths},
             'runtime':{'python':platform.python_version(),'platform':platform.platform(),
                        'openbabel':subprocess.check_output([args.obabel,'-V'],text=True).strip(),
                        'python_extension_sha256':sha(extension),'wasm_sha256':sha(args.wasm_module.with_name('chematic_wasm_bg.wasm')),
                        'openbabel_probe_sha256':sha(args.openbabel_probe)},
             'results':{'accepted_python_wasm_byte_parity':sum(c['accepted'] for c in cases),
                        'accepted_openbabel_element_coordinate_charge_crossreads':sum(c['accepted'] for c in cases),
                        'refused_python_and_wasm':sum(not c['accepted'] for c in cases),'failures':0},
             'issue':'https://github.com/openbabel/openbabel/issues/2101','fixture_sha256':sha(FIXTURE),
             'packet_sha256':hashlib.sha256(raw).hexdigest(),
             'scope':'Fixed-column serialization only; no atom typing, docking, protonation, torsion tree or release claim.'}
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary['results'],indent=2))


if __name__=='__main__':
    main()
