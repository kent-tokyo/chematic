#!/usr/bin/env python3
"""Independent all-record, raw-distance and RDKit rigid-alignment controls."""
import argparse
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def coordinates(text):
    rows = text.split("@<TRIPOS>ATOM\n", 1)[1].split("@<TRIPOS>BOND", 1)[0].splitlines()
    return [[float(v) for v in row.split()[2:5]] for row in rows if row.strip()]


def replace_coordinates(text, points, name):
    prefix, atoms = text.split("@<TRIPOS>ATOM\n", 1)
    rows, tail = atoms.split("@<TRIPOS>BOND", 1)
    replaced = []
    for row, point in zip(rows.splitlines(), points):
        fields = row.split()
        fields[2:5] = [f"{v:.9f}" for v in point]
        replaced.append(" ".join(fields))
    header = prefix.splitlines()
    header[1] = name
    return "\n".join(header) + "\n@<TRIPOS>ATOM\n" + "\n".join(replaced) + "\n@<TRIPOS>BOND" + tail


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--baseline-package", required=True)
    parser.add_argument("--obabel", default="/opt/homebrew/bin/obabel")
    parser.add_argument("--obrms", default="/opt/homebrew/bin/obrms")
    args = parser.parse_args()
    import chematic
    import chematic.chematic as native
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    import numpy as np

    reference = (HERE / "references.mol2").read_text()
    probe = (HERE / "probe.mol2").read_text()
    baseline_script = '''import hashlib,json,chematic,chematic.chematic as native,sys
text=open(sys.argv[1]).read()
results={}
for name in ["from_mol2","parse_mol2_record_json"]:
 try:getattr(chematic,name)(text); results[name]="unexpected success"
 except ValueError as error:results[name]={"error":type(error).__name__,"message":str(error)}
results["new_apis_present"]=[hasattr(chematic,name) for name in ["parse_mol2_records_json","mol2_reference_rmsds_json"]]
results["module_sha256"]=hashlib.sha256(open(native.__file__,"rb").read()).hexdigest()
print(json.dumps(results))'''
    baseline = json.loads(subprocess.check_output(
        [sys.executable, "-c", baseline_script, str(HERE / "references.mol2")],
        env={**os.environ, "PYTHONPATH": args.baseline_package}, text=True))
    assert baseline["new_apis_present"] == [False, False]
    assert all(isinstance(baseline[name], dict) for name in ("from_mol2", "parse_mol2_record_json"))
    parsed = json.loads(chematic.parse_mol2_records_json(reference))["records"]
    assert len(parsed) == 5 and "".join(r["source"] for r in parsed) == reference
    raw = json.loads(chematic.mol2_reference_rmsds_json(reference, probe))["results"]
    aligned = json.loads(chematic.mol2_reference_rmsds_json(reference, probe, align=True))["results"]
    assert [r["rmsd"] for r in raw] == [0, 1, 2, 3, 5]
    assert all(r["rmsd"] < 1e-12 for r in aligned)

    template = probe.replace("40 H4", "40 N4").replace("1.000000 H 1", "1.000000 N.3 1")
    base = coordinates(template)
    rotations = [lambda x,y,z:(x,y,z), lambda x,y,z:(-y,x,z), lambda x,y,z:(z,y,-x),
                 lambda x,y,z:(x,-z,y), lambda x,y,z:(-x,y,-z),
                 lambda x,y,z:(math.cos(.37)*x-math.sin(.37)*y, math.sin(.37)*x+math.cos(.37)*y,z)]
    texts = []
    for r, rotation in enumerate(rotations):
        for t, translation in enumerate([(0,0,0),(1,0,0),(0,2,0),(0,0,3),(3,4,0)]):
            points = [[a+b for a,b in zip(rotation(*p),translation)] for p in base]
            texts.append(replace_coordinates(template, points, f"rigid_{r}_{t}"))
    for i in range(1,9):
        points = [p[:] for p in base]
        points[2][2] += i*.125
        texts.append(replace_coordinates(template, points, f"distortion_{i}"))
    texts.append(replace_coordinates(template, [[-x,y,z] for x,y,z in base], "mirror"))
    stream = "".join(texts)
    output_raw = json.loads(chematic.mol2_reference_rmsds_json(stream, template))["results"]
    output_aligned = json.loads(chematic.mol2_reference_rmsds_json(stream, template, align=True))["results"]
    assert len(output_raw) == len(output_aligned) == len(texts) == 39

    def rdkit_mol(points):
        mol = Chem.RWMol()
        for element in ["C","C","O","N"]: mol.AddAtom(Chem.Atom(element))
        for a,b in [(0,1),(1,2),(0,3)]: mol.AddBond(a,b,Chem.BondType.SINGLE)
        conf = Chem.Conformer(4)
        for i, point in enumerate(points): conf.SetAtomPosition(i, point)
        mol.AddConformer(conf)
        return mol.GetMol()

    rd_probe = rdkit_mol(base)
    comparisons = []
    for text, actual_raw, actual_aligned in zip(texts, output_raw, output_aligned):
        points = coordinates(text)
        expected_raw = math.sqrt(math.fsum((a-b)**2 for p,q in zip(points,base) for a,b in zip(p,q))/4)
        expected_aligned, _ = rdMolAlign.GetAlignmentTransform(
            rd_probe, rdkit_mol(points), atomMap=[(i,i) for i in range(4)], reflect=False, maxIters=1000)
        # Independent SVD solution with determinant +1; unlike the native
        # Jacobi solver and RDKit's iterative alignment, this uses NumPy/LAPACK.
        p = np.array(base) - np.mean(base, axis=0)
        q = np.array(points) - np.mean(points, axis=0)
        u, _, vt = np.linalg.svd(p.T @ q)
        rotation = u @ np.diag([1, 1, np.linalg.det(u @ vt)]) @ vt
        expected_svd = float(np.sqrt(np.sum((p @ rotation-q)**2)/4))
        assert abs(actual_raw["rmsd"]-expected_raw) < 1e-10, actual_raw
        assert abs(actual_aligned["rmsd"]-expected_svd) < 1e-9, actual_aligned
        # RDKit's iterative result retains ~1e-6 residual on the .37-radian
        # rigid control even with maxIters=1000. Keep this separate from the
        # tighter SVD gate, and record both raw oracle values rather than
        # treating RDKit's near-zero residual as exact mathematical zero.
        assert abs(actual_aligned["rmsd"]-expected_aligned) < 1e-5, actual_aligned
        comparisons.append({"name":actual_raw["reference_name"], "native_raw":actual_raw["rmsd"],
                            "math_raw":expected_raw, "native_aligned":actual_aligned["rmsd"],
                            "svd_aligned":expected_svd, "rdkit_aligned":expected_aligned})
    assert comparisons[-1]["native_aligned"] > .1, "mirror must not be accepted by a reflection"

    prefix, atom_part = template.split("@<TRIPOS>ATOM\n",1)
    atoms, tail = atom_part.split("@<TRIPOS>BOND",1)
    rows = atoms.splitlines()
    rejected = 0
    for permutation in itertools.permutations(range(4)):
        if permutation == (0,1,2,3): continue
        changed = prefix + "@<TRIPOS>ATOM\n" + "\n".join(rows[i] for i in permutation) + "\n@<TRIPOS>BOND" + tail
        try: chematic.mol2_reference_rmsds_json(changed, template)
        except ValueError as error:
            assert "identity/order" in str(error)
            rejected += 1
        else: raise AssertionError(permutation)
    assert rejected == 23

    commands = {
        "reported_order":[args.obrms,str(HERE/"references.mol2"),str(HERE/"probe.mol2")],
        "documented_all_vs_one":[args.obrms,"-f",str(HERE/"probe.mol2"),str(HERE/"references.mol2")],
        "documented_aligned":[args.obrms,"-m","-f",str(HERE/"probe.mol2"),str(HERE/"references.mol2")],
    }
    upstream = {}
    for name, command in commands.items():
        result = subprocess.run(command,capture_output=True,text=True,check=True,timeout=30)
        upstream[name] = {"command":command,"stdout":result.stdout,"stderr":result.stderr}
    assert len(upstream["reported_order"]["stdout"].splitlines()) == 1
    assert [float(line.split()[-1]) for line in upstream["documented_all_vs_one"]["stdout"].splitlines()] == [0,1,2,3,5]
    assert [float(line.split()[-1]) for line in upstream["documented_aligned"]["stdout"].splitlines()] == [0]*5
    report = {"schema":"chematic.mol2-rmsd-validation.v1", "scope":"Original bounded all-vs-one controls; fixed correspondence, no speed/release/general chemistry claim.",
              "issue":"https://github.com/openbabel/openbabel/issues/2671", "baseline":baseline,
              "candidate_module_sha256":digest(native.__file__),"rdkit_version":rdkit.__version__,
              "numpy_version":np.__version__,
              "tolerances":{"raw_math":1e-10,"aligned_svd":1e-9,"aligned_rdkit":1e-5,"rdkit_max_iters":1000},
              "max_alignment_delta_svd":max(abs(r["native_aligned"]-r["svd_aligned"]) for r in comparisons),
              "max_alignment_delta_rdkit":max(abs(r["native_aligned"]-r["rdkit_aligned"]) for r in comparisons),
              "python":sys.version,"openbabel_version":subprocess.check_output([args.obabel,"-V"],text=True).strip(),
              "fixture_sha256":{p.name:digest(p) for p in [HERE/"probe.mol2",HERE/"references.mol2",HERE/"expected.json"]},
              "analytical_raw":raw,"analytical_aligned":aligned,"generated_comparisons":comparisons,
              "generated_stream_sha256":hashlib.sha256(stream.encode()).hexdigest(),
              "identity_permutations_rejected":rejected,"upstream_diagnostic":upstream,
              "interpretation":"The default upstream command pairs nth records, per its manual. Reversing inputs with -f works; one-line output is not proof of an upstream parser defect."}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"analytical_records":5,"raw_and_rdkit_alignment_comparisons":39,"identity_permutations_rejected":23}))


if __name__ == "__main__":
    main()
