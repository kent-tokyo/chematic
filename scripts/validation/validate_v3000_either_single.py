#!/usr/bin/env python3
"""Fixed Open Babel #2429 input/output diagnostic; no literal SMILES oracle."""
import argparse
import datetime
import gzip
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", choices=["chematic", "rdkit", "openbabel"], required=True)
    parser.add_argument("--python-root")
    parser.add_argument("--label", required=True)
    parser.add_argument("--cases", type=Path, default=Path(__file__).resolve().parents[2] / "validation/openbabel-2429-2026-10-24/cases.json")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.python_root:
        sys.path.insert(0, args.python_root)
    if args.engine == "chematic":
        import chematic
        from chematic import chematic as native
        runtime = {"path": native.__file__, "sha256": sha(native.__file__)}

        def convert(text):
            mol, name, xy = chematic.from_mol_v3000_with_coords(text)
            before = mol.smiles
            output = mol.to_mol_v3000(xy, name=name, strict=True)
            after = chematic.from_mol_v3000(output).smiles
            return before, output, after, ""
    elif args.engine == "rdkit":
        from rdkit import Chem, rdBase
        runtime = {"version": rdBase.rdkitVersion, "path": rdBase.__file__, "sha256": sha(rdBase.__file__)}

        def convert(text):
            mol = Chem.MolFromMolBlock(text)
            assert mol is not None
            output = Chem.MolToV3KMolBlock(mol)
            reread = Chem.MolFromMolBlock(output)
            assert reread is not None
            return Chem.MolToSmiles(mol), output, Chem.MolToSmiles(reread), ""
    else:
        import shutil
        executable = shutil.which("obabel")
        assert executable
        runtime = {"version": subprocess.check_output([executable, "-V"], text=True).strip(), "path": executable, "sha256": sha(executable)}

        def run(text, output):
            return subprocess.run([executable, "-imol", "-o" + output, *(["-x3"] if output == "mol" else [])], input=text, capture_output=True, text=True, check=True, timeout=10)

        def convert(text):
            before = run(text, "can")
            written = run(text, "mol")
            after = run(written.stdout, "can")
            return before.stdout.split()[0], written.stdout, after.stdout.split()[0], before.stderr + written.stderr + after.stderr

    rows = []
    for case in json.loads(args.cases.read_text())["cases"]:
        before, output, after, diagnostics = convert(case["input"])
        defined = lambda smiles: "/" in smiles or "\\" in smiles
        matched = defined(before) == case["expected_defined_ez"] and defined(after) == case["expected_defined_ez"]
        rows.append({"id": case["id"], "expected_defined_ez": case["expected_defined_ez"], "input_smiles": before, "output": output, "reread_smiles": after, "matches_expected_stereo": matched, "diagnostics": diagnostics})
    result = {"schema": "chematic.v3000-either-validation.v1", "engine": args.engine, "label": args.label, "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(), "platform": platform.platform(), "runtime": runtime, "observer_sha256": sha(__file__), "cases_sha256": sha(args.cases), "matched": sum(r["matches_expected_stereo"] for r in rows), "total": len(rows), "cases": rows}
    args.output.write_bytes(gzip.compress((json.dumps(result, indent=2, allow_nan=False) + "\n").encode(), mtime=0))
    print(f"{args.label}: {result['matched']}/{result['total']} input/output stereo controls match")


if __name__ == "__main__":
    main()
