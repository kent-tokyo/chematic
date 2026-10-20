#!/usr/bin/env python3
"""Synthetic CDJSON projection/oracle packet for Open Babel #2763."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import subprocess
from concurrent.futures import ThreadPoolExecutor

import chematic
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "validation/openbabel_issues/cdjson/fixtures.json"
OUTPUT = ROOT / "validation/results/openbabel-cdjson-2026-10-15.json"


def normalized(v):
    return json.dumps(v, sort_keys=True, separators=(",", ":")).encode()


def obabel(input_format, output_format, text, *options):
    result = subprocess.run(["obabel", f"-i{input_format}", f"-o{output_format}", *options],
                            input=text, text=True, capture_output=True, timeout=60, check=True)
    return result.stdout, result.stderr


def oracle(molecule):
    model = Chem.RWMol()
    for atom in molecule["a"]:
        a = Chem.Atom(atom.get("l", "C"))
        a.SetFormalCharge(atom.get("c", 0))
        if atom.get("m", -1) != -1:
            a.SetIsotope(atom["m"])
        if atom.get("h", -1) >= 0:
            a.SetNoImplicit(True)
            a.SetNumExplicitHs(atom["h"])
        model.AddAtom(a)
    for bond in molecule.get("b", []):
        model.AddBond(bond["b"], bond["e"], {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}[bond.get("o", 1)])
    Chem.SanitizeMol(model)
    return model


def variants(record):
    molecule = record["payload"]
    count = len(molecule["a"])
    for variant in range(32):
        order = list(range(count))
        if variant & 1:
            order.reverse()
        shift = (variant // 2) % count
        order = order[shift:] + order[:shift]
        mapping = {old: new for new, old in enumerate(order)}
        atoms = []
        for index, old in enumerate(order):
            atom = dict(molecule["a"][old])
            atom.update(x=variant*.12345678901234567+index*20., y=-index*3.125,
                        z=(variant % 3)*.9876543210123456, i=f"{record['id']}-{variant}-a{index}")
            if variant & 2 and atom.get("l", "C") == "C":
                atom["l"] = "C"
            atoms.append(atom)
        bonds = [dict(b, b=mapping[b["b"]], e=mapping[b["e"]], i=f"{record['id']}-{variant}-b{i}")
                 for i, b in enumerate(molecule.get("b", []))]
        if variant & 4:
            bonds.reverse()
        yield f"{record['id']}-{variant:02d}",dict(a=atoms,b=bonds,i=f"{record['id']}-{variant}"),record


def validate_projection(source, projection):
    mol = projection["mol"]
    model = oracle(source)
    assert mol.formula == rdMolDescriptors.CalcMolFormula(model)
    assert mol.atom_count == model.GetNumAtoms()
    assert projection["coordinates"] == [[a["x"],a["y"],a.get("z",0)] for a in source["a"]]
    assert projection["id"] == source.get("i")
    assert projection["atom_ids"] == [a.get("i") for a in source["a"]]
    assert projection["bond_ids"] == [b.get("i") for b in source.get("b",[])]
    written = chematic.write_cdjson_graph(mol,projection["coordinates"])
    recovered = chematic.parse_cdjson(written)["molecules"][0]
    assert recovered["mol"].formula == mol.formula
    assert recovered["coordinates"] == projection["coordinates"]
    output = json.loads(written)["m"][0]
    for a,b,reference in zip(source["a"],output["a"],model.GetAtoms()):
        assert b["l"] == a.get("l","C")
        assert b["c"] == a.get("c",0)
        assert b.get("m",-1) == a.get("m",-1)
        assert b["h"] == reference.GetNumExplicitHs()+reference.GetNumImplicitHs()
    assert [(b["b"],b["e"],b["o"]) for b in output["b"]] == [(b["b"],b["e"],b.get("o",1)) for b in source.get("b",[])]
    return output, mol.formula


def main():
    assert rdBase.rdkitVersion == "2026.09.1", rdBase.rdkitVersion
    version = subprocess.check_output(["obabel","-V"],text=True).strip()
    assert version.startswith("Open Babel 3.2.1 "),version
    manifest=json.loads(FIXTURES.read_text())
    cases=[]
    for record in manifest["records"]:
        assert hashlib.sha256(normalized(record["payload"])).hexdigest()==record["sha256"]
        cases.extend(variants(record))
    source=json.dumps(dict(m=[m for _,m,_ in cases]))
    document=chematic.parse_cdjson(source)
    assert document["source_text"]==source
    assert len(document["molecules"])==len(cases)
    packet=hashlib.sha256();written=[];formulas=[]
    for (identifier,molecule,record),projection in zip(cases,document["molecules"]):
        output,formula=validate_projection(molecule,projection)
        assert formula == record["expected"]["formula"], identifier
        written.append(output);formulas.append(formula)
        packet.update(normalized([identifier,molecule,output,formula])+b"\n")
    batch_output,batch_stderr=obabel("cdjson","smi",source)
    assert len(batch_output.splitlines())==1,(len(batch_output.splitlines()),batch_stderr)
    def single_read(molecule):
        return obabel("cdjson","smi",json.dumps(dict(m=[molecule])))
    with ThreadPoolExecutor(max_workers=4) as pool:
        source_reads=list(pool.map(single_read,[m for _,m,_ in cases]))
        written_reads=list(pool.map(single_read,written))
    diagnostics=[]
    for lane,reads in [("source_read",source_reads),("native_write_cross_read",written_reads)]:
        exact=known=0;differences=[]
        for (identifier,_,record),(text,stderr),expected in zip(cases,reads,formulas):
            lines=text.splitlines();assert len(lines)==1,(identifier,text,stderr)
            line=lines[0]
            smiles=line.split()[0]
            model=Chem.MolFromSmiles(smiles);assert model is not None,(identifier,smiles)
            actual=rdMolDescriptors.CalcMolFormula(model)
            if actual==expected:exact+=1
            else:
                assert record["known_openbabel_zero_override_gap"],(identifier,actual,expected)
                known+=1;differences.append(dict(id=identifier,expected=expected,openbabel_formula=actual,openbabel_smiles=smiles))
        diagnostics.append(dict(lane=lane,exact_formula_cases=exact,known_zero_override_differences=known,
                                differences=differences))
    simple_smiles=[r["smiles"] for r in manifest["writer_controls"]]
    output,stderr=obabel("smi","cdjson","\n".join(simple_smiles)+"\n","--gen2D")
    imported=chematic.parse_cdjson(output)
    assert len(imported["molecules"])==len(simple_smiles),(len(imported["molecules"]),stderr)
    for smiles,item in zip(simple_smiles,imported["molecules"]):
        assert item["mol"].formula==rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(smiles))
    result=dict(schema_version=1,issue=dict(number=2763,state="OPEN",verified_date="2026-10-15",updated_at="2025-01-10T08:55:05Z",url="https://github.com/openbabel/openbabel/issues/2763"),
        scope="Strict molecular ChemDoodle JSON subset, explicit zero H overrides, source/ID preservation and graph writer",
        oracle=dict(openbabel=version,rdkit=rdBase.rdkitVersion,profile="Explicit source H overrides, RDKit ordinary atom defaults for -1/omitted H, graph reconstruction and formula comparison"),
        fixtures=dict(path=str(FIXTURES.relative_to(ROOT)),sha256=hashlib.sha256(FIXTURES.read_bytes()).hexdigest()),
        results=dict(native_rdkit_cases=len(cases),native_roundtrips=len(cases),source_id_coordinate_cases=len(cases),
                     normalized_packet_sha256=packet.hexdigest(),unexpected_mismatches=0,
                     openbabel_diagnostics=diagnostics,openbabel_writer_import_controls=len(simple_smiles),
                     multi_molecule_read_control=dict(input_molecules=len(cases),chematic_molecules=len(document["molecules"]),
                         openbabel_molecules=len(batch_output.splitlines()),openbabel_stderr=batch_stderr),
                     reported_default_carbon_bug_reproduced_on_3_2_1=False,
                     current_zero_h_override_gap_reproduced=True))
    OUTPUT.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({k:v for k,v in result["results"].items() if k!="openbabel_diagnostics"},indent=2))
    for d in diagnostics: print(d["lane"],d["exact_formula_cases"],d["known_zero_override_differences"])


if __name__=="__main__":main()
