"""ASE cross-read of periodic conversion; separate Open Babel loss diagnostic."""

from io import StringIO
import hashlib
import json
from pathlib import Path
import subprocess

import ase
from ase.io import read
import chematic as ch
import chematic.chematic as native
import numpy as np


def main():
    rng = np.random.default_rng(282219)
    cells = [
        np.diag([4.0, 5.0, 6.0]),
        np.array([[4.0, 0.5, 0], [0, 5, 0.75], [0.25, 0, 6]]),
        np.array([[0, 4.0, 0.5], [0.75, 0, 5], [6, 0.25, 0]]),
        np.diag([-4.0, 5.0, 6.0]),
    ]
    comparisons = 0
    max_position_error = 0.0
    max_cell_error = 0.0
    for cell in cells:
        for _ in range(8):
            for count in (1, 3, 8):
                fractions = rng.uniform(-1.2, 1.8, (count, 3))
                positions = fractions @ cell
                elements = [["O", "Si", "O", "N"][i % 4] for i in range(count)]
                for mode in ("none", "all", "alternate"):
                    frozen = [
                        mode == "all" or (mode == "alternate" and i % 2 == 0)
                        for i in range(count)
                    ]
                    rows = ["$periodic 3", "$lattice angs"]
                    rows += [" ".join(map(str, row)) for row in cell]
                    rows += ["$coord frac"]
                    rows += [
                        " ".join(map(str, row))
                        + " "
                        + element.lower()
                        + (" f" if fixed else "")
                        for row, element, fixed in zip(fractions, elements, frozen)
                    ]
                    source = "\n".join(rows) + "\n$end\n"
                    converted = ch.turbomole_to_poscar_checked(source)
                    order = list(converted.source_indices)
                    assert sorted(order) == list(range(count))
                    atoms = read(StringIO(converted.text), format="vasp")
                    np.testing.assert_allclose(
                        atoms.cell.array, cell, atol=1e-12, rtol=0
                    )
                    np.testing.assert_allclose(
                        atoms.positions, positions[order], atol=1e-11, rtol=0
                    )
                    assert atoms.get_chemical_symbols() == [elements[i] for i in order]
                    fixed_indices = {
                        int(i)
                        for constraint in atoms.constraints
                        for i in constraint.get_indices()
                    }
                    assert fixed_indices == {
                        j for j, i in enumerate(order) if frozen[i]
                    }
                    back = ch.poscar_to_turbomole_checked(converted.text)
                    geometry = ch.parse_turbomole(back.text)
                    np.testing.assert_allclose(
                        geometry["lattice_angstrom"],
                        atoms.cell.array,
                        atol=1e-12,
                        rtol=0,
                    )
                    np.testing.assert_allclose(
                        geometry["positions_angstrom"],
                        atoms.positions,
                        atol=1e-11,
                        rtol=0,
                    )
                    assert geometry["elements"] == atoms.get_chemical_symbols()
                    assert geometry["frozen"] == [
                        i in fixed_indices for i in range(count)
                    ]
                    max_position_error = max(
                        max_position_error,
                        float(
                            np.max(
                                np.abs(
                                    np.array(geometry["positions_angstrom"])
                                    - positions[order]
                                )
                            )
                        ),
                    )
                    max_cell_error = max(
                        max_cell_error,
                        float(
                            np.max(
                                np.abs(np.array(geometry["lattice_angstrom"]) - cell)
                            )
                        ),
                    )
                    comparisons += 1
    # Original scalar Si geometry inspired by #2822, not a repository fixture copy.
    diagnostic = (
        "$periodic 3\n$lattice\n10 0 0\n0 10 0\n0 0 10\n$coord\n0 0 0 si\n$end\n"
    )
    probe = subprocess.run(
        ["obabel", "-itmol", "-ovasp"],
        input=diagnostic,
        text=True,
        capture_output=True,
        check=False,
    )
    lines = probe.stdout.splitlines()
    report = {
        "ase": ase.__version__,
        "numpy": np.__version__,
        "seed": 282219,
        "forward_cross_reads": comparisons,
        "reverse_geometry_checks": comparisons,
        "mismatches": 0,
        "max_position_error_angstrom": max_position_error,
        "max_lattice_component_error_angstrom": max_cell_error,
        "native_sha256": hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
        "openbabel_diagnostic": {
            "version": subprocess.run(
                ["obabel", "-V"], text=True, capture_output=True
            ).stdout.strip(),
            "returncode": probe.returncode,
            "stderr": probe.stderr,
            "vasp_output": probe.stdout,
            "zero_lattice": len(lines) >= 5
            and all(float(v) == 0 for row in lines[2:5] for v in row.split()),
        },
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
