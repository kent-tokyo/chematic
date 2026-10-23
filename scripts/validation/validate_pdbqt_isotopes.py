#!/usr/bin/env python3
"""Check isotope sidecars, fixed columns and real Python/WASM parity.

Run with candidate PYTHONPATH. RDKit is an independent source-label observer;
this does not test docking quality or flexible-ligand torsion preparation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
from pathlib import Path

import chematic
import chematic.chematic as native

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "validation/openbabel_issues/pdbqt-isotopes/cases.json"


def run(argv, **kwargs):
    return subprocess.run(argv, capture_output=True, text=True, check=True, **kwargs).stdout


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def observe_source(rdkit_python, cases):
    program = """
import json, sys
from rdkit import Chem, rdBase
result=[]
for c in json.load(sys.stdin):
    p=Chem.SmilesParserParams();p.removeHs=False
    m=Chem.MolFromSmiles(c['smiles'],p)
    assert m is not None
    result.append(dict(atomic_numbers=[a.GetAtomicNum() for a in m.GetAtoms()],
        isotopes_omitted=[dict(source_atom_index=a.GetIdx(),pdbqt_serial=a.GetIdx()+1,
            atomic_number=a.GetAtomicNum(),mass_number=a.GetIsotope())
            for a in m.GetAtoms() if a.GetIsotope()]))
print(json.dumps(dict(version=rdBase.rdkitVersion,cases=result)))
"""
    result = json.loads(run([rdkit_python, "-c", program], input=json.dumps(cases)))
    for case, observed in zip(cases, result["cases"], strict=True):
        assert observed == {k: case[k] for k in ("atomic_numbers", "isotopes_omitted")}
    return result["version"]


def wasm_reports(module, cases):
    program = """
const fs=require('node:fs'),w=require(process.argv[1]);
const result=JSON.parse(fs.readFileSync(0,'utf8')).map(c=>{
  const xyz=c.atomic_numbers.map((_,i)=>[i,-2.345,6.789]);
  const q=xyz.map(()=>-.123);
  let refused=false;
  try {w.smiles_to_pdbqt_checked(c.smiles,JSON.stringify(xyz),JSON.stringify(q),'LIG');}
  catch(e) {if(!String(e).includes('isotope'))throw e;refused=true;}
  return {refused,report:JSON.parse(w.smiles_to_pdbqt_isotope_projection(
      c.smiles,JSON.stringify(xyz),JSON.stringify(q),'LIG'))};
});
process.stdout.write(JSON.stringify(result));
"""
    return json.loads(run(["node", "-e", program, str(Path(module).resolve())], input=json.dumps(cases)))


def decode(text, case):
    """Independent fixed-column reader; no chematic readback is used here."""
    type_elements = {"H": 1, "HD": 1, "C": 6, "A": 6, "N": 7, "NA": 7,
                     "O": 8, "OA": 8, "S": 16, "SA": 16, "P": 15,
                     "F": 9, "Cl": 17, "Br": 35, "I": 53, "Fe": 26}
    lines = [line for line in text.splitlines() if line.startswith("ATOM")]
    assert len(lines) == len(case["atomic_numbers"])
    for i, line in enumerate(lines):
        assert len(line) == 79
        assert int(line[6:11]) == i+1
        assert type_elements[line[77:79].strip()] == case["atomic_numbers"][i]
        assert [float(line[a:b]) for a, b in [(30,38),(38,46),(46,54)]] == [i,-2.345,6.789]
        assert float(line[70:76]) == -.123


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wasm-module", required=True)
    parser.add_argument("--rdkit-python", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cases = json.loads(FIXTURE.read_text())
    rdkit = observe_source(args.rdkit_python, cases)
    wasm = wasm_reports(args.wasm_module, cases)
    results = []
    for case, js in zip(cases, wasm, strict=True):
        mol = chematic.from_smiles(case["smiles"])
        coords = [[i,-2.345,6.789] for i in range(mol.atom_count)]
        q = [-.123] * mol.atom_count
        refused = False
        try:
            chematic.write_pdbqt_checked(mol, coords, q)
        except ValueError as exc:
            assert "isotope" in str(exc)
            refused = True
        assert refused == js["refused"] == bool(case["isotopes_omitted"])
        report = chematic.project_pdbqt_isotopes(mol, coords, q)
        assert report == js["report"]
        assert report["isotopes_omitted"] == case["isotopes_omitted"]
        decode(report["pdbqt"], case)
        results.append(dict(name=case["name"], checked_refused=refused,
                            atoms=mol.atom_count, **report))
    sources = ["crates/chematic-mol/src/pdbqt.rs", "crates/chematic-mol/src/pdbqt/checked.rs",
               "crates/chematic-mol/src/pdbqt/isotopes.rs", "crates/chematic-py/src/formats/pdbqt.rs",
               "crates/chematic-wasm/src/pdbqt_checked.rs", str(FIXTURE.relative_to(ROOT)),
               str(Path(__file__).resolve().relative_to(ROOT))]
    evidence = dict(scope="rigid PDBQT isotope omission; no docking/torsion quality claim",
        head=run(["git","rev-parse","HEAD"], cwd=ROOT).strip(),
        worktree_status=run(["git","status","--short","--untracked-files=no"], cwd=ROOT),
        host=platform.platform(), python=platform.python_version(), rdkit=rdkit,
        python_extension_sha256=sha(native.__file__),
        node=run(["node","--version"]).strip(),
        wasm_sha256=sha(Path(args.wasm_module).with_suffix(".wasm").with_name("chematic_wasm_bg.wasm")),
        source_sha256={path:sha(ROOT/path) for path in sources},
        accepted_projections=len(results), checked_refusals=sum(r["checked_refused"] for r in results),
        isotope_labels=sum(len(r["isotopes_omitted"]) for r in results), results=results)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(evidence, indent=2)+"\n")
    print(f"PDBQT isotope evidence passed: {len(results)} projections, "
          f"{evidence['checked_refusals']} refusals, {evidence['isotope_labels']} labels")


if __name__ == "__main__":
    main()
