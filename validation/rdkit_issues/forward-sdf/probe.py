"""RDKit #9601 stream reproduction and optional native-reader comparisons.

Run with --compare after building chematic. The missing-read reproduction is
isolated in a child with a two-second timeout and disabled RDKit logging.
"""
import argparse
import gzip
import hashlib
import io
import json
import subprocess
import sys
from pathlib import Path

from rdkit import Chem, rdBase


def missing_read():
    code = (
        "from rdkit import Chem, rdBase; rdBase.DisableLog('rdApp.*'); "
        "s = Chem.ForwardSDMolSupplier(42); print('constructed', flush=True); "
        "print(len(list(s)), flush=True)"
    )
    try:
        result = subprocess.run([sys.executable, "-c", code], capture_output=True,
                                text=True, timeout=2, check=False)
        return {"timed_out": False, "stdout": result.stdout,
                "stderr": result.stderr, "returncode": result.returncode}
    except subprocess.TimeoutExpired as error:
        return {"timed_out": True, "timeout_seconds": 2,
                "stdout": (error.stdout or b"").decode("utf-8")}


def compare():
    import chematic
    from chematic import rdkit_compat

    cases = []
    smiles = ["CCO", "CN", "c1ccccc1", "F[C@H](Cl)Br", "[13CH3][NH3+]",
              "F/C=C/Cl", "F/C=C\\Cl", "C1CCCCC1", "CC(=O)O", "O=[N+]([O-])c1ccccc1"]
    for v3000 in (False, True):
        stream = io.StringIO()
        writer = Chem.SDWriter(stream)
        writer.SetForceV3000(v3000)
        for index, smi in enumerate(smiles):
            mol = Chem.MolFromSmiles(smi)
            mol.SetProp("_Name", f"case-{index}")
            mol.SetProp("Activity", str(index))
            writer.write(mol)
        writer.flush()
        data = stream.getvalue().encode("utf-8")
        expected = list(Chem.ForwardSDMolSupplier(io.BytesIO(data)))
        for api in (chematic, rdkit_compat):
            for compressed in (False, True):
                source = (gzip.GzipFile(fileobj=io.BytesIO(gzip.compress(data)))
                          if compressed else io.BytesIO(data))
                actual = list(api.ForwardSDMolSupplier(source))
                assert len(actual) == len(expected) == len(smiles)
                for left, right in zip(actual, expected):
                    smi = left.smiles if api is chematic else api.MolToSmiles(left)
                    actual_key = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
                    expected_key = Chem.MolToSmiles(right)
                    assert actual_key == expected_key, (smi, expected_key)
                    assert left.GetProp("_Name") == right.GetProp("_Name")
                    assert left.GetProp("Activity") == right.GetProp("Activity")
                    cases.append({"v3000": v3000, "gzip": compressed,
                                  "api": api.__name__, "smiles": actual_key})
                source.close()
    module = __import__("chematic.chematic", fromlist=["__file__"]).__file__
    return {"passed": len(cases), "cases": cases,
            "module_sha256": hashlib.sha256(Path(module).read_bytes()).hexdigest()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--compare", action="store_true")
    args = parser.parse_args()
    result = {"rdkit_version": rdBase.rdkitVersion, "missing_read": missing_read()}
    if args.compare:
        result["comparison"] = compare()
    print(json.dumps(result, ensure_ascii=False, indent=2))
