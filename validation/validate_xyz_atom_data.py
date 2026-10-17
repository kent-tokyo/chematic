#!/usr/bin/env python3
"""Seeded xyz+data -> extxyz cross-read with ASE 3.26.0 and native reader."""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import platform
import random
import subprocess

import ase
from ase.io import read
import numpy as np
import chematic


def validate() -> dict:
    if ase.__version__ != "3.26.0":
        raise RuntimeError("requires ASE 3.26.0")
    rng = random.Random(1703)
    digest = hashlib.sha256()
    records = atoms = scalar_values = 0
    comments = ["", "1 2", "T", "_JSON {}", '日本語 "quoted" ' + chr(92) + "path"]
    for width in (1, 3, 9):
        name = {1: "data_scalar", 3: "data_vector", 9: "data_tensor"}[width]
        for size in (1, 2, 5, 17, 64):
            symbols = [("C", "Cl", "O", "Na", "H")[i % 5] for i in range(size)]
            xyz = [[rng.uniform(-5, 5) for _ in range(3)] for _ in range(size)]
            xyz[0][0] = -0.0
            data = [[rng.uniform(-1e4, 1e4) for _ in range(width)] for _ in range(size)]
            data[0][0] = float.fromhex("0x0.0000000000001p-1022")
            if size > 1:
                data[1][0] = -0.0
            if size > 2:
                data[2][-1] = float.fromhex("0x1.fffffffffffffp+1023")
            rows = [
                symbol + " " + " ".join(format(v, ".17g") for v in point + values)
                for symbol, point, values in zip(symbols, xyz, data)
            ]
            for comment in comments:
                for newline in ("\n", "\r\n"):
                    source = (
                        newline.join([str(size), comment] + rows)
                        + newline
                        + "\t"
                        + newline
                    )
                    raw = source.encode("utf-8")
                    f = chematic.parse_xyz_atom_data(raw, components=width)
                    assert f.source_text.encode("utf-8") == raw
                    normalized = chematic.parse_xyz_atom_data(
                        f.to_xyz(), components=width
                    )
                    assert (
                        normalized.values == f.values
                        and normalized.positions == f.positions
                    )
                    ext = f.to_extxyz(name)
                    a = read(io.StringIO(ext), format="extxyz")
                    native = chematic.from_extxyz(ext)
                    assert a.get_chemical_symbols() == symbols
                    assert (
                        np.asarray(a.positions, dtype=np.float64).tobytes()
                        == np.asarray(xyz, dtype=np.float64).tobytes()
                    )
                    assert (
                        np.asarray(native["coords"], dtype=np.float64).tobytes()
                        == np.asarray(xyz, dtype=np.float64).tobytes()
                    )
                    expected = np.asarray(data, dtype=np.float64).reshape(size, width)
                    assert (
                        np.asarray(a.arrays[name], dtype=np.float64)
                        .reshape(size, width)
                        .tobytes()
                        == expected.tobytes()
                    )
                    assert (
                        np.asarray(
                            native["properties"][name], dtype=np.float64
                        ).tobytes()
                        == expected.tobytes()
                    )
                    assert a.info["comment"] == comment
                    assert json.loads(native["info"]["comment"][6:]) == comment
                    assert not a.pbc.any()
                    digest.update(raw)
                    digest.update(ext.encode("utf-8"))
                    records += 1
                    atoms += size
                    scalar_values += size * width
    # Independently verify the complete canonical element table with ASE.
    from ase.data import chemical_symbols

    elements = chemical_symbols[1:119]
    source = f"{len(elements)}\nall canonical elements\n" + "".join(
        f"{s} {i*10} 0 0 .25\n" for i, s in enumerate(elements)
    )
    f = chematic.parse_xyz_atom_data(source, components=1)
    assert list(f.elements) == elements
    assert (
        read(
            io.StringIO(f.to_extxyz("payload")), format="extxyz"
        ).get_chemical_symbols()
        == elements
    )
    probe = "3\nOriginal asymmetric vector fixture\nC 0 0 0 .125 -.25 .5\nCl 1.5 .1 0 -1 .75 .2\nO -.3 2 .2 .001 .002 .003\n"
    result = subprocess.run(
        ["/opt/homebrew/bin/obabel", "-ixyz", "-oxyz"],
        input=probe,
        capture_output=True,
        text=True,
        check=True,
    )
    return {
        "schema_version": 1,
        "scope": "single-frame data interchange, not topology or physical units",
        "seed": 1703,
        "python": platform.python_version(),
        "ase": ase.__version__,
        "numpy": np.__version__,
        "source_records": records,
        "ase_cross_reads": records,
        "native_cross_reads": records,
        "atom_rows": atoms,
        "payload_scalar_values": scalar_values,
        "canonical_element_check": len(elements),
        "value_check": "binary64 byte equality for coordinates and all payload values; negative zero, subnormal and largest finite payload included",
        "comment_check": "ASE decodes _JSON metadata; native raw metadata recovered with json.loads(value[6:])",
        "source_and_extxyz_sha256": digest.hexdigest(),
        "openbabel_probe": {
            "version": subprocess.check_output(
                ["/opt/homebrew/bin/obabel", "-V"], text=True
            ).strip(),
            "input": probe,
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_atom_column_counts": [
                len(row.split()) for row in result.stdout.splitlines()[2:5]
            ],
        },
        "exclusions": [
            "no topology parity",
            "no physical meaning or tensor reshaping",
            "no lattice/trajectory parsing",
            "no fresh native build",
            "no performance claim",
        ],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = json.dumps(validate(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(result)
    print(result, end="")
