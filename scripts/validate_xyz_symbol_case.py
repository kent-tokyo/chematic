#!/usr/bin/env python3
"""Check XYZ species against RDKit 2026.09.1's periodic table; diagnose Open Babel #2333.

Provide a fresh source extension and pinned RDKit on PYTHONPATH. Optional
--baseline-pythonpath checks a retained pre-fix source runtime in a subprocess.
Optional --wasm-module repeats the identity checks in real Node/WASM.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import gzip
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ["crates/chematic-3d/src/xyz.rs", "crates/chematic-mol/src/xyz.rs",
          "crates/chematic-mol/src/xyz/extxyz.rs", "crates/chematic-py/src/formats/xyz.rs"]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--obabel",default="/opt/homebrew/bin/obabel")
    p.add_argument("--wasm-module",type=Path)
    p.add_argument("--baseline-pythonpath",type=Path)
    p.add_argument("--output",type=Path,default=ROOT/"validation/results/openbabel-xyz-elements-2026-10-16.json")
    args = p.parse_args()
    import chematic
    from rdkit import Chem,rdBase
    assert rdBase.rdkitVersion == "2026.09.1"
    table = Chem.GetPeriodicTable()
    cases = []
    for number in range(1,119):
        canonical = table.GetElementSymbol(number)
        for mode,token in [("canonical",canonical),("uppercase",canonical.upper())]:
            cases.append(dict(atomic_number=number,canonical=canonical,mode=mode,token=token,
                              text=f"1\noriginal species probe\n{token} 0.125 -2.5 1.25\n"))
    counts = Counter()
    for case in cases:
        mol,coords = chematic.from_xyz(case["text"])
        ext = chematic.from_extxyz(case["text"])
        for name,m,c in [("plain",mol,coords),("extxyz",ext["mol"],ext["coords"])]:
            assert m.atom_count == 1 and c == [[0.125,-2.5,1.25]]
            symbol = m.to_xyz(c).splitlines()[2].split()[0]
            assert symbol == case["canonical"]
            assert table.GetAtomicNumber(symbol) == case["atomic_number"]
            counts[f"python_{name}_matches"] += 1
    baseline = None
    if args.baseline_pythonpath:
        code = """import sys,json,chematic
cases=json.load(sys.stdin);out=[]
for c in cases:
 try:
  m,p=chematic.from_xyz(c['text']);out.append({'token':c['token'],'accepted':True,'symbol':m.to_xyz(p).splitlines()[2].split()[0]})
 except ValueError as e:out.append({'token':c['token'],'accepted':False,'error':str(e)})
print(json.dumps(out))
"""
        env = dict(os.environ,PYTHONPATH=str(args.baseline_pythonpath))
        rows=json.loads(subprocess.check_output([sys.executable,"-c",code],input=json.dumps(cases),text=True,env=env,timeout=30))
        assert len(rows)==len(cases)
        baseline=dict(scope="Retained pre-fix source runtime; not a released-package comparison",
                      extension_sha256=sha(args.baseline_pythonpath/"chematic/chematic.so"),rows=rows,
                      rejected_uppercase=sum(not r["accepted"] for c,r in zip(cases,rows) if c["mode"]=="uppercase"))
    wasm = None
    if args.wasm_module:
        js = """const fs=require('fs'),w=require(process.argv[1]);const cases=JSON.parse(fs.readFileSync(0,'utf8'));let n=0;for(const c of cases){for(const read of [w.mol_from_xyz,w.mol_from_extxyz]){const m=read(c.text);if(m.atom_count()!==1||w.to_xyz(m).split('\\n')[2].trim().split(/\\s+/)[0]!==c.canonical)throw Error('element mismatch');m.free();n++;}const f=JSON.parse(w.extxyz_frame_json(c.text));if(JSON.stringify(f.coords)!==JSON.stringify([[0.125,-2.5,1.25]]))throw Error('coordinate mismatch');}console.log(n);"""
        n=int(subprocess.check_output(["node","-e",js,str(args.wasm_module.resolve())],input=json.dumps(cases),text=True,timeout=60))
        assert n == 472
        counts["wasm_plain_and_extxyz_matches"]=n
        wasm={f.name:sha(f) for f in args.wasm_module.parent.iterdir() if f.is_file()}
    ob_rows=[]
    for mode in ["canonical","uppercase"]:
        selected=[c for c in cases if c["mode"]==mode]
        text="118\noriginal periodic-table probe\n"+"".join(f"{c['token']} {i*100} 0 0\n" for i,c in enumerate(selected))
        r=subprocess.run([args.obabel,"-ixyz","-oxyz"],input=text,text=True,capture_output=True,timeout=10)
        assert r.returncode==0,r.stderr
        actual=[line.split()[0] for line in r.stdout.splitlines()[2:] if line.strip()]
        assert len(actual)==118,(mode,r.stdout,r.stderr)
        mismatches=[dict(atomic_number=c["atomic_number"],input=c["token"],expected=c["canonical"],actual=a) for c,a in zip(selected,actual) if a!=c["canonical"]]
        ob_rows.append(dict(mode=mode,input=text,input_sha256=hashlib.sha256(text.encode()).hexdigest(),
                            stdout=r.stdout,symbols=actual,mismatches=mismatches,stderr=r.stderr))
    binary=Path(importlib.import_module("chematic.chematic").__file__)
    fixture=ROOT/"validation/openbabel_issues/xyz-elements/chloromethane-uppercase.xyz"
    args.output.parent.mkdir(parents=True,exist_ok=True)
    packet=args.output.with_suffix(".packet.json.gz")
    payload=dict(cases=cases,baseline=baseline,openbabel=ob_rows)
    with packet.open("wb") as f:
        with gzip.GzipFile(fileobj=f,mode="wb",mtime=0) as z:
            z.write((json.dumps(payload,separators=(",",":"))+"\n").encode())
    baseline_summary = None if baseline is None else {k:v for k,v in baseline.items() if k!="rows"}
    result=dict(schema="chematic.xyz-species-validation.v1",utc=datetime.now(timezone.utc).isoformat(),
                issue="https://github.com/openbabel/openbabel/issues/2333",
                source_base=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
                validation_runner_sha256=sha(__file__),
                source_sha256={f:sha(ROOT/f) for f in SOURCE},python_extension_sha256=sha(binary),
                rdkit_version=rdBase.rdkitVersion,rdkit_runtime_sha256=sha(rdBase.__file__),wasm_runtime_sha256=wasm,
                fixture_sha256=sha(fixture),casing_cases=len(cases),distinct_species_tokens=len({c["token"] for c in cases}),
                packet=packet.name,packet_sha256=sha(packet),counts=dict(counts),candidate_failures=0,baseline=baseline_summary,
                openbabel_version=subprocess.check_output([args.obabel,"-V"],text=True).strip(),
                openbabel_binary_sha256=sha(args.obabel),
                openbabel=[dict(mode=r["mode"],matches=118-len(r["mismatches"]),mismatch_count=len(r["mismatches"]),samples=r["mismatches"][:8]) for r in ob_rows],
                boundary="Element/coordinate identity in native XYZ/extxyz only; no bond inference, chemical perception or exact source preservation. Core symbol lookup, PDB names and SMILES grammar remain case-sensitive. Source candidate only.")
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(dict(casing_cases=len(cases),counts=dict(counts),baseline_rejected_uppercase=baseline and baseline["rejected_uppercase"],
                         openbabel_mismatches={r["mode"]:len(r["mismatches"]) for r in ob_rows},candidate_failures=0)))


if __name__ == "__main__": main()
