#!/usr/bin/env python3
"""Compare constructed rooted SMILES with RDKit 2026.09.1 indexed CIP labels.

Generate --cases with cargo run -p chematic-smiles --example rooted_stereo_probe.
Provide a fresh chematic Python extension and the pinned RDKit on PYTHONPATH.
Optional --wasm-module runs the same parsed inputs in a fresh real Node/WASM build.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import importlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ["crates/chematic-smiles/src/rooted.rs", "crates/chematic-smiles/src/writer.rs",
           "crates/chematic-core/src/molecule/renumber.rs", "crates/chematic-smiles/tests/support/rooted_bridgehead.rs",
           "crates/chematic-smiles/examples/rooted_stereo_probe.rs", "crates/chematic-py/src/mol_methods/rooted.rs",
           "crates/chematic-wasm/src/rooted.rs"]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build(Chem, case):
    rw = Chem.RWMol()
    for i, symbol in enumerate(["C","C","C","O","C","C","C"]):
        atom = Chem.Atom(symbol)
        atom.SetAtomMapNum(i+1)
        if i in (2,4):
            ccw = (i == 2) ^ bool(case["flips"] & (1 if i == 2 else 2))
            atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW if ccw else Chem.ChiralType.CHI_TETRAHEDRAL_CW)
            if case["explicit_h"]:
                atom.SetNumExplicitHs(1)
                atom.SetNoImplicit(True)
        rw.AddAtom(atom)
    edges = [(2,3,1),(3,4,1),(2,1,1),(4,0,1),(1,0,2),(5,6,1),(2,6,1),(4,5,1)]
    k = case["rotation"]
    edges = edges[k:]+edges[:k]
    if case["reverse"]: edges.reverse()
    for a,b,order in edges:
        rw.AddBond(a,b,Chem.BondType.DOUBLE if order == 2 else Chem.BondType.SINGLE)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    return mol


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cases", type=Path, required=True)
    p.add_argument("--wasm-module", type=Path)
    p.add_argument("--output", type=Path, default=ROOT/"validation/results/rdkit-rooted-stereo-2026-10-16.json")
    args = p.parse_args()
    import chematic
    from rdkit import Chem, rdBase
    from rdkit.Chem import rdCIPLabeler
    assert rdBase.rdkitVersion == "2026.09.1"

    def labels(mol):
        assert mol is not None
        rdCIPLabeler.AssignCIPLabels(mol)
        return {str(a.GetAtomMapNum()):a.GetProp("_CIPCode") for a in mol.GetAtoms() if a.HasProp("_CIPCode")}

    records = json.loads(args.cases.read_text())
    assert len(records) == 896
    assert len({(c["explicit_h"],c["flips"],c["rotation"],c["reverse"],c["root"]) for c in records}) == 896
    counts, failures = Counter(), []
    for case in records:
        original = build(Chem, case)
        expected = labels(original)
        parsed = Chem.MolFromSmiles(case["smiles"])
        actual = labels(parsed)
        assert actual == expected, (case, expected, actual)
        assert [a.GetAtomMapNum()-1 for a in parsed.GetAtoms()] == case["atom_order"]
        assert case["atom_order"][0] == case["root"]
        legacy_text = Chem.MolToSmiles(original, rootedAtAtom=case["root"], canonical=False)
        legacy = labels(Chem.MolFromSmiles(legacy_text))
        key = "explicit_h" if case["explicit_h"] else "implicit_h"
        counts[f"rdkit_{key}_matches"] += legacy == expected
        counts[f"rdkit_{key}_mismatches"] += legacy != expected
        counts["native_matches"] += 1
        py_text, py_order = chematic.from_smiles(case["smiles"]).smiles_rooted_checked(case["root"])
        py_mol = Chem.MolFromSmiles(py_text)
        assert labels(py_mol) == expected
        input_maps = [a.GetAtomMapNum() for a in parsed.GetAtoms()]
        assert [a.GetAtomMapNum() for a in py_mol.GetAtoms()] == [input_maps[i] for i in py_order]
        counts["python_matches"] += 1
        case.update(expected_cip=expected,native_cip=actual,rdkit_smiles=legacy_text,rdkit_cip=legacy,
                    python_smiles=py_text,python_atom_order=py_order)
        if legacy != expected and len(failures)<4: failures.append(case.copy())
    wasm_runtime = None
    if args.wasm_module:
        js = """const fs=require('fs');const w=require(process.argv[1]);const cases=JSON.parse(fs.readFileSync(0,'utf8'));let n=0;for(const c of cases){const m=w.parse_smiles(c.smiles);const o=JSON.parse(m.smiles_rooted_checked_json(c.root));if(o.smiles!==c.python_smiles||JSON.stringify(o.atom_order)!==JSON.stringify(c.python_atom_order))throw Error('binding mismatch');m.free();n++;}console.log(n);"""
        n = int(subprocess.check_output(["node","-e",js,str(args.wasm_module.resolve())],input=json.dumps(records),text=True,timeout=60))
        assert n == 896
        counts["wasm_exact_python_matches"] = n
        wasm_runtime = {f.name:sha(f) for f in args.wasm_module.parent.glob("*") if f.is_file()}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    packet = args.output.with_suffix(".packet.json.gz")
    with packet.open("wb") as f:
        with gzip.GzipFile(fileobj=f,mode="wb",mtime=0) as z:
            z.write((json.dumps(records,separators=(",",":"))+"\n").encode())
    binary = Path(importlib.import_module("chematic.chematic").__file__)
    result = dict(schema="chematic.rooted-stereo-validation.v1",utc=datetime.now(timezone.utc).isoformat(),
                  issue="https://github.com/rdkit/rdkit/issues/9646",source_base=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
                  source_sha256={f:sha(ROOT/f) for f in SOURCES},rdkit_version=rdBase.rdkitVersion,
                  rdkit_runtime_sha256={"rdBase":sha(rdBase.__file__),"rdCIPLabeler":sha(rdCIPLabeler.__file__)},
                  python_extension_sha256=sha(binary),wasm_runtime_sha256=wasm_runtime,
                  generated_input_sha256=sha(args.cases),packet=packet.name,packet_sha256=sha(packet),
                  cases=len(records),counts=dict(counts),candidate_failures=0,rdkit_mismatch_samples=failures,
                  boundary="Formal tetrahedral configuration and indexed CIP preservation only; not proof that every assignment embeds in 3D. Native checked opt-in API, not a byte-for-byte RDKit writer or a change to its compatibility lane; maximum 256 atoms.")
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(dict(cases=len(records),counts=dict(counts),candidate_failures=0)))


if __name__ == "__main__": main()
