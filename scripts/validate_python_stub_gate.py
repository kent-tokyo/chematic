#!/usr/bin/env python3
"""Exercise the stub CLI on an actual wheel and original archive mutations.

The baseline wheel is read only. Mutations live in a temporary directory and
are deleted on completion. Vendor stub code is never copied into the report.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import warnings
import zipfile
from pathlib import Path
from check_python_stubs import check_syntax

ROOT=Path(__file__).resolve().parents[1]


def check(wheel, report, expected):
    command=[sys.executable,str(ROOT/'scripts/check_python_stubs.py'),'--wheel',str(wheel),'--report',str(report)]
    result=subprocess.run(command,text=True,capture_output=True)
    assert result.returncode==expected,(command,result.returncode,result.stderr)
    return dict(command=command,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr,
                report=json.loads(report.read_text()))


def mutate(source, destination, mode):
    target='chematic/__init__.pyi';marker='chematic/py.typed'
    endings={'annotation':b'\nNone: int\n','default-order':b'\ndef broken(default=0, required): ...\n',
             'duplicate-argument':b'\ndef broken(value, value): ...\n','stale-valid':b'\nnot_in_source: int\n'}
    with zipfile.ZipFile(source) as old,zipfile.ZipFile(destination,'w') as new:
        for info in old.infolist():
            name=info.filename
            if (mode=='missing-stub' and name==target) or (mode=='missing-marker' and name==marker):
                continue
            if name==target and mode in endings:
                new.writestr(copy.copy(info),old.read(info)+endings[mode])
            elif name==marker and mode=='stale-marker':
                new.writestr(copy.copy(info),b'partial\n')
            else:
                with old.open(info) as src,new.open(copy.copy(info),'w') as dst:
                    shutil.copyfileobj(src,dst,1024*1024)
        if mode=='duplicate-member':
            with warnings.catch_warnings():
                warnings.simplefilter('ignore',UserWarning)
                new.writestr(target,b'None: int\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--wheel',required=True,type=Path)
    p.add_argument('--rdkit-stubs',type=Path)
    p.add_argument('--rdkit-version',default='unspecified')
    p.add_argument('--output',type=Path,default=ROOT/'validation/results/python-stub-gate-2026-10-16')
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    runs=[check(args.wheel,args.output/'real-wheel.json',0)]
    with tempfile.TemporaryDirectory(prefix='chematic-stub-negative-') as temp:
        wheel=Path(temp)/args.wheel.name
        for mode in ['annotation','default-order','duplicate-argument','missing-stub','missing-marker','stale-valid','stale-marker','duplicate-member']:
            mutate(args.wheel,wheel,mode)
            result=check(wheel,Path(temp)/'report.json',1)
            result['mutation']=mode
            runs.append(result)
    upstream=[];count=0
    if args.rdkit_stubs:
        files=sorted(args.rdkit_stubs.rglob('*.pyi'));count=len(files)
        assert count>0,'no RDKit stubs found'
        for file in files:
            data=file.read_bytes();relative=file.relative_to(args.rdkit_stubs).as_posix()
            issue=check_syntax(data,relative)
            if issue:upstream.append(dict(path=relative,sha256=hashlib.sha256(data).hexdigest(),diagnostic=issue))
    paths=['scripts/check_python_stubs.py','scripts/tests/test_check_python_stubs.py','scripts/validate_python_stub_gate.py','.github/workflows/ci.yml','.github/workflows/publish-pypi.yml']
    result=dict(schema='chematic.python-stub-gate-validation.v1',host_python=sys.version.split()[0],
                source_base=subprocess.check_output(['git','merge-base','HEAD','origin/main'],cwd=ROOT,text=True).strip(),
                source_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in paths},
                runs=runs,rdkit=dict(version=args.rdkit_version,files=count,failures=upstream),
                scope='Existing retained macOS wheel, not a new build or publication. Grammar/compiler policy only; native modules and stubs are not imported or executed; no type/API/ABI parity claim.')
    (args.output/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Actual wheel accepted; 8 mutated archives refused with CLI exit 1')
    print('RDKit parse failures:',sum(e['diagnostic']['phase']=='parse' for e in upstream),'additional compiler-policy failures:',sum(e['diagnostic']['phase']=='compile' for e in upstream))


if __name__=='__main__':
    main()
