"""Compare checked SDF graph/metadata output across compression policies."""
import gzip
import hashlib
import json
import tempfile
from pathlib import Path

import chematic


def snapshot(path, mode):
    output = []
    for index, record, error in chematic.iter_sdf_checked(path, compression=mode):
        assert error is None, error
        output.append({"index": index, "smiles": record.smiles, "name": record.name,
                       "properties": record.properties(),
                       "stereo_diagnostics": record.stereo_diagnostics()})
    return output


def main():
    inputs = ["CCO", "CN", "c1ccccc1", "F[C@H](Cl)Br", "F[C@@H](Cl)Br",
              "[13CH3][NH3+]", "F/C=C/Cl", "F/C=C\\Cl", "C1CCCCC1", "CC(=O)O"]
    results = []
    with tempfile.TemporaryDirectory(prefix="chematic-gzip-compare-") as folder:
        root = Path(folder)
        for v3000 in [False, True]:
            source = root / "original.sdf"
            with chematic.SDWriter(source) as writer:
                writer.SetForceV3000(v3000)
                for index, smi in enumerate(inputs):
                    mol = chematic.from_smiles(smi)
                    mol.SetProp("_Name", f"compound-日本語-{index}")
                    mol.SetProp("Activity", str(index))
                    writer.write(mol)
            expected = snapshot(source, "plain")
            assert len(expected) == len(inputs)
            data = source.read_bytes()
            split = data.index("日".encode()) + 1
            members = gzip.compress(data[:split], mtime=0) + gzip.compress(data[split:], mtime=0)
            variants = [("plain_auto", data, "auto", "fake.gz"),
                        ("gzip_auto", gzip.compress(data, mtime=0), "auto", "opaque.dat"),
                        ("gzip_explicit", gzip.compress(data, mtime=0), "gzip", "opaque.dat"),
                        ("members_auto", members, "auto", "opaque.dat")]
            for name, payload, mode, filename in variants:
                path = root / filename
                path.write_bytes(payload)
                actual = snapshot(path, mode)
                assert actual == expected, (v3000, name)
                for record in actual:
                    results.append({"v3000": v3000, "variant": name, **record})
    module = Path(__import__("chematic.chematic", fromlist=["__file__"]).__file__)
    print(json.dumps({"comparisons_passed": len(results), "cases": results,
                      "module_sha256": hashlib.sha256(module.read_bytes()).hexdigest()},
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
