#!/usr/bin/env python3
"""Observe exact source Molden geometry/vibration fields across round trips.

Derived Atoms geometry added to a vibration-only file is allowed. Atom-unit
changes are source-encoding changes, not by themselves physical errors. No
bond, frequency, intensity or displacement normalization is performed.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[2]
AXES=("atom_units","atom_rows","fr_coords_bohr","frequencies","normal_modes_bohr","ir_values","raman_values")


def numbers(words):
    result=[float(v.replace("D","e").replace("d","e")) for v in words]
    if not all(map(math.isfinite,result)):raise ValueError("nonfinite output")
    return result


def observe(text):
    lines=[line.strip() for line in text.splitlines() if line.strip()]
    if not lines or lines.pop(0).lower()!="[molden format]":raise ValueError("no Molden output")
    sections={};current=None
    for line in lines:
        if line.startswith("["):
            key,tail=line[1:].split("]",1);key=key.upper()
            if key in sections:raise ValueError("duplicate output section")
            sections[key]=(tail.strip(),[]);current=key
        elif current is None:raise ValueError("unscoped data")
        else:sections[current][1].append(line.split())
    r=dict(atoms=None,fr_coords=None,frequencies=None,normal_modes=None,intensities=None)
    if "ATOMS" in sections:
        unit,rows=sections["ATOMS"];unit={"angs":"angstrom","au":"bohr"}[unit.lower()]
        atoms=[]
        for row in rows:
            if len(row)!=6:raise ValueError("invalid atom row")
            atoms.append(dict(symbol=row[0],index=int(row[1]),atomic_number=int(row[2]),position=numbers(row[3:])))
        if not atoms:raise ValueError("empty atom geometry")
        r["atoms"]=dict(unit=unit,atoms=atoms)
    if "FR-COORD" in sections:
        coords=[]
        for row in sections["FR-COORD"][1]:
            if len(row)!=4:raise ValueError("invalid FR coordinate")
            coords.append(dict(symbol=row[0],position_bohr=numbers(row[1:])))
        r["fr_coords"]=coords
    if "FREQ" in sections:
        freq=[]
        for row in sections["FREQ"][1]:
            if len(row)!=1:raise ValueError("invalid frequency row")
            freq.append(numbers(row)[0])
        r["frequencies"]=freq
    if "FR-NORM-COORD" in sections:
        if not r["fr_coords"]:raise ValueError("modes without coordinates")
        rows=iter(sections["FR-NORM-COORD"][1]);modes=[]
        try:
            for row in rows:
                if len(row)!=2 or row[0].lower()!="vibration":raise ValueError("invalid mode header")
                vectors=[]
                for _ in r["fr_coords"]:
                    vector=numbers(next(rows))
                    if len(vector)!=3:raise ValueError("invalid mode vector")
                    vectors.append(vector)
                modes.append(dict(index=int(row[1]),displacements_bohr=vectors))
        except StopIteration as error:raise ValueError("truncated mode") from error
        r["normal_modes"]=modes
    if "INT" in sections:
        r["intensities"]=[]
        for row in sections["INT"][1]:
            v=numbers(row)
            if len(v) not in (1,2):raise ValueError("invalid intensity row")
            r["intensities"].append(dict(ir=v[0],raman=v[1] if len(v)==2 else None))
    if r["atoms"] is None and not r["fr_coords"]:raise ValueError("missing geometry")
    return r


def axes(r):
    a=r["atoms"];i=r["intensities"]
    raman=[[index,v["raman"]] for index,v in enumerate(i or []) if v["raman"] is not None]
    return dict(atom_units=a["unit"] if a else None,atom_rows=a["atoms"] if a else None,
                fr_coords_bohr=r["fr_coords"],frequencies=r["frequencies"],normal_modes_bohr=r["normal_modes"],
                ir_values=[v["ir"] for v in i] if i is not None else None,raman_values=raman or None)


def run(command,text):
    try:r=subprocess.run(command,input=text,text=True,encoding="utf-8",capture_output=True,timeout=30,check=False)
    except (OSError,subprocess.TimeoutExpired) as e:return dict(status="execution_failed",error=str(e))
    receipt=dict(returncode=r.returncode,stdout=r.stdout,stderr=r.stderr)
    if r.returncode:return dict(status="rejected",**receipt)
    try:record=observe(r.stdout)
    except (ValueError,KeyError,IndexError) as e:return dict(status="unobservable",error=str(e),**receipt)
    return dict(status="observed",record=record,**receipt)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--cli",type=Path,required=True)
    parser.add_argument("--obabel",type=Path,default=Path("/opt/homebrew/bin/obabel"));parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args();packet_path=ROOT/"validation/molden_vibrations/cases.json"
    packet=json.loads(packet_path.read_text(encoding="utf-8"))
    commands=dict(chematic=[str(args.cli.resolve()),"molden"],openbabel=[str(args.obabel.resolve()),"-imolden","-omolden","-ab","-l1"])
    totals={k:{a:0 for a in AXES} for k in commands};eligible={a:0 for a in AXES};rows=[]
    for case in packet["cases"]:
        expected=axes(case["expected"]);row=dict(id=case["id"],input_sha256=hashlib.sha256(case["input"].encode()).hexdigest())
        for a in AXES:eligible[a]+=expected[a] is not None
        for engine,command in commands.items():
            result=run(command,case["input"]);got=axes(result["record"]) if result["status"]=="observed" else {}
            result["preserved"]={a:None if expected[a] is None else a in got and got[a]==expected[a] for a in AXES}
            for a,v in result["preserved"].items():totals[engine][a]+=v is True
            result["all_source_fields_preserved"]=all(v is not False for v in result["preserved"].values()) and result["status"]=="observed"
            row[engine]=result
        rows.append(row)
    sources=["crates/chematic-mol/src/molden.rs","crates/chematic-mol/src/molden/reader.rs","crates/chematic-mol/src/molden/writer.rs","crates/chematic-cli/src/main.rs","crates/chematic-cli/src/dispatch.rs"]
    report=dict(schema_version=1,scope="original Molden geometry/vibration profile; exact declared source fields; extra derived geometry allowed",
                case_count=len(rows),timeout_seconds=30,commands=commands,eligible_counts=eligible,preserved_counts=totals,
                observed_counts={k:sum(c[k]["status"]=="observed" for c in rows) for k in commands},
                complete_source_preservation={k:sum(c[k]["all_source_fields_preserved"] for c in rows) for k in commands},
                versions={"chematic":subprocess.check_output([str(args.cli.resolve()),"--version"],text=True).strip(),"openbabel":subprocess.check_output([str(args.obabel.resolve()),"-V"],text=True).strip()},
                packet_sha256=hashlib.sha256(packet_path.read_bytes()).hexdigest(),source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
                binary_sha256={k:hashlib.sha256(Path(v[0]).read_bytes()).hexdigest() for k,v in commands.items()},cases=rows)
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:report[k] for k in ["case_count","observed_counts","eligible_counts","preserved_counts","complete_source_preservation"]},indent=2))
    return int(report["complete_source_preservation"]["chematic"]!=len(rows))


if __name__=="__main__":raise SystemExit(main())
