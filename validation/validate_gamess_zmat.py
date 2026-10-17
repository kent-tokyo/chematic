"""Original synthetic GAMESS geometry/IZMAT corpus and independent cross-reads."""

import hashlib
import json
import re
import subprocess
from pathlib import Path

import numpy as np
import chematic
import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolTransforms


def definitions(text):
    if "$ZMAT" not in text:
        return []
    body = text.split("$ZMAT IZMAT(1)=", 1)[1].split("$END", 1)[0]
    integers = [int(s.strip()) for s in body.replace("\n", "").split(",")]
    result = []
    while integers:
        code = integers.pop(0)
        assert code in (1, 2, 3)
        width = code + 1
        assert len(integers) >= width
        result.append((code, integers[:width]))
        del integers[:width]
    return result


def main():
    rng = np.random.default_rng(2523)
    corpus_hash = hashlib.sha256()
    records = checks = cross_reads = 0
    max_reference_difference = 0.0
    symbols = ["C", "N", "O", "F", "H", "Cl", "Br", "S"]
    for n in (1, 2, 3, 4, 8):
        for sample in range(64):
            coords = rng.normal(size=(n, 3)) + rng.uniform(-10.0, 10.0, 3)
            refs = [rng.choice(i, min(i, 3), replace=False).tolist() for i in range(n)]
            elements = symbols[:n]
            for digits in (8, 17):
                text = chematic.write_gamess_zmat(
                    elements,
                    coords.tolist(),
                    refs,
                    charge=-1,
                    multiplicity=2,
                    significant_digits=digits,
                )
                assert max(map(len, text.splitlines())) <= 80
                data = [
                    line
                    for line in text.split("C1\n", 1)[1]
                    .split("$END", 1)[0]
                    .splitlines()
                    if line.strip()
                ]
                assert [row.split()[0] for row in data] == elements
                rounded = np.array(
                    [[float(s) for s in row.split()[2:]] for row in data]
                )
                if digits == 17:
                    np.testing.assert_array_equal(rounded, coords)
                else:
                    np.testing.assert_allclose(rounded, coords, rtol=1e-7, atol=1e-12)
                parsed = definitions(text)
                expected = [
                    (kind, [i + 1] + [r + 1 for r in refs[i][:kind]])
                    for i in range(n)
                    for kind in range(1, min(i, 3) + 1)
                ]
                assert parsed == expected
                assert int(re.search(r"NZVAR=(\d+)", text)[1]) == len(expected)
                rows = chematic.cartesian_to_zmatrix(elements, rounded.tolist(), refs)
                mol = Chem.RWMol()
                for symbol in elements:
                    mol.AddAtom(Chem.Atom(symbol))
                conf = Chem.Conformer(n)
                for i, position in enumerate(rounded):
                    conf.SetAtomPosition(i, position)
                mol.AddConformer(conf)
                for code, ids in parsed:
                    indices = [i - 1 for i in ids]
                    func = (
                        rdMolTransforms.GetBondLength,
                        rdMolTransforms.GetAngleDeg,
                        rdMolTransforms.GetDihedralDeg,
                    )[code - 1]
                    actual = func(conf, *indices)
                    native = rows[indices[0]]["values"][code - 1]
                    difference = (
                        abs((actual - native + 180.0) % 360.0 - 180.0)
                        if code == 3
                        else abs(actual - native)
                    )
                    assert difference < 1e-7, (code, ids, actual, native)
                    max_reference_difference = max(max_reference_difference, difference)
                    checks += 1
                xyz = f"{n}\nOriginal synthetic geometry\n" + "".join(
                    f"{s} {p[0]:.17e} {p[1]:.17e} {p[2]:.17e}\n"
                    for s, p in zip(elements, coords)
                )
                direct = chematic.xyz_to_gamess_zmat(
                    xyz, refs, charge=-1, multiplicity=2, significant_digits=digits
                )
                assert direct == text
                if sample % 8 == 0:
                    proc = subprocess.run(
                        ["obabel", "-igamin", "-oxyz"],
                        input=text,
                        text=True,
                        capture_output=True,
                        check=True,
                    )
                    xyz_lines = proc.stdout.splitlines()
                    assert int(xyz_lines[0]) == n, proc.stderr
                    assert [s.split()[0] for s in xyz_lines[2 : 2 + n]] == elements
                    got = np.array(
                        [
                            [float(v) for v in s.split()[1:4]]
                            for s in xyz_lines[2 : 2 + n]
                        ]
                    )
                    np.testing.assert_allclose(got, rounded, rtol=0, atol=1e-4)
                    cross_reads += 1
                corpus_hash.update(text.encode())
                records += 1
    print(
        json.dumps(
            {
                "scope": "Cartesian frame and explicit IZMAT definitions; no GAMESS execution",
                "records": records,
                "direct_xyz_equivalence_checks": records,
                "rdkit_internal_value_checks": checks,
                "maximum_rdkit_internal_difference": max_reference_difference,
                "openbabel_cartesian_cross_reads": cross_reads,
                "corpus_sha256": corpus_hash.hexdigest(),
                "rdkit_version": rdkit.__version__,
                "numpy_version": np.__version__,
                "openbabel_version": subprocess.check_output(
                    ["obabel", "-V"], text=True
                ).strip(),
                "native_extension_sha256": hashlib.sha256(
                    Path(chematic.chematic.__file__).read_bytes()
                ).hexdigest(),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
