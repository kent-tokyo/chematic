#!/usr/bin/env python3
"""Cross-read checked XSF geometry with ASE 3.26.0 and Open Babel 3.2.1."""
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
from ase.data import chemical_symbols
from ase.io import read
import numpy as np
from chematic import Lattice, PeriodicStructure, Site, write_xsf_checked

FIXTURE = Path(__file__).parent / "openbabel_issues/xsf/ordered-triclinic.cif"


def validate() -> dict:
    if ase.__version__ != "3.26.0":
        raise RuntimeError("requires ASE 3.26.0")
    rng = random.Random(2252)
    matrices = [
        [[3.4, 0.0, 0.0], [0.0, 3.4, 0.0], [0.0, 0.0, 3.4]],
        [[3.2, 0.12, -0.3], [0.4, 4.5, 0.28], [-0.15, 0.7, 5.6]],
        [[2.0, 0.1, 0.2], [0.0, 4.0, 0.3], [0.0, 0.0, 6.0]],
        [[2.0, 0.1, 0.0], [0.0, 3.0, 0.2], [0.3, 0.0, -4.0]],
    ]
    cases = atoms = ob_cases = 0
    max_ob_error = 0.0
    digest = hashlib.sha256()
    for matrix in matrices:
        for size in (1, 2, 5, 13, 64):
            sites = [
                Site(
                    [(("Na", "Cl", "C", "O", "H")[i % 5], 1.0)],
                    tuple(rng.uniform(-1, 2) for _ in range(3)),
                    f"original_site_{i}",
                )
                for i in range(size)
            ]
            original = PeriodicStructure(Lattice.from_matrix(matrix), sites)
            for structure in (
                original,
                original.wrap_into_cell(),
                original.make_supercell((2, 1, 1)),
            ):
                text = write_xsf_checked(structure, allow_metadata_loss=True)
                digest.update(text.encode())
                a = read(io.StringIO(text), format="xsf")
                assert a.cell.array.tobytes() == structure.lattice.matrix.tobytes()
                assert (
                    a.positions.tobytes() == structure.cartesian_positions().tobytes()
                )
                symbols = [s.species[0][0] for s in structure.sites]
                assert a.get_chemical_symbols() == symbols
                assert a.pbc.all()
                if size <= 13:
                    r = subprocess.run(
                        ["/opt/homebrew/bin/obabel", "-ixsf", "-oxyz"],
                        input=text,
                        capture_output=True,
                        text=True,
                        check=True,
                    )
                    rows = r.stdout.splitlines()
                    assert int(rows[0]) == structure.site_count()
                    assert [row.split()[0] for row in rows[2:]] == symbols
                    xyz = np.asarray(
                        [[float(x) for x in row.split()[1:4]] for row in rows[2:]]
                    )
                    delta = float(np.max(np.abs(xyz - structure.cartesian_positions())))
                    assert delta <= 5.1e-6, delta
                    max_ob_error = max(max_ob_error, delta)
                    ob_cases += 1
                cases += 1
                atoms += structure.site_count()
    cif = FIXTURE.read_text()
    structure = PeriodicStructure.from_cif(cif)
    assert structure.symmetry_status.is_complete_cell
    text = write_xsf_checked(structure, allow_metadata_loss=True)
    a = read(io.StringIO(text), format="xsf")
    assert a.positions.tobytes() == structure.cartesian_positions().tobytes()
    assert a.cell.array.tobytes() == structure.lattice.matrix.tobytes()
    numbers = PeriodicStructure(
        Lattice.cubic(2000.0),
        [
            Site([(symbol, 1.0)], (i / 118, 0.0, 0.0))
            for i, symbol in enumerate(chemical_symbols[1:119])
        ],
    )
    assert (
        read(
            io.StringIO(write_xsf_checked(numbers, allow_metadata_loss=True)),
            format="xsf",
        ).get_chemical_symbols()
        == chemical_symbols[1:119]
    )
    attempt = subprocess.run(
        ["/opt/homebrew/bin/obabel", "-icif", "-oxsf"],
        input=cif,
        capture_output=True,
        text=True,
    )
    formats = subprocess.check_output(
        ["/opt/homebrew/bin/obabel", "-L", "formats"], text=True
    )
    xsf_listing = [line for line in formats.splitlines() if line.startswith("xsf --")]
    return {
        "schema_version": 1,
        "seed": 2252,
        "scope": "ordered 3D-periodic geometry only; no topology or metadata parity",
        "python": platform.python_version(),
        "ase": ase.__version__,
        "numpy": np.__version__,
        "ase_corpus_crossreads": cases,
        "corpus_atom_rows": atoms,
        "openbabel_xyz_crossreads": ob_cases,
        "max_openbabel_xyz_error_angstrom": max_ob_error,
        "openbabel_xyz_tolerance_angstrom": 5.1e-6,
        "additional_cif_pipeline_crossreads": 1,
        "additional_all_element_check": 118,
        "cell_coordinate_check": "ASE binary64 byte equality to native geometry getters; no basis reduction or wrapping by writer",
        "corpus_xsf_sha256": digest.hexdigest(),
        "fixture_sha256": hashlib.sha256(FIXTURE.read_bytes()).hexdigest(),
        "openbabel": {
            "version": subprocess.check_output(
                ["/opt/homebrew/bin/obabel", "-V"], text=True
            ).strip(),
            "xsf_listing": xsf_listing,
            "cif_to_xsf_returncode": attempt.returncode,
            "cif_to_xsf_stdout": attempt.stdout[:2000],
            "cif_to_xsf_stderr": attempt.stderr[:2000],
        },
        "limits": [
            "no XCrySDen GUI run",
            "no metadata roundtrip",
            "no occupancy collapse",
            "no space-group operation inference",
            "no primitive/conventional reduction",
            "no fresh native build or speed claim",
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
