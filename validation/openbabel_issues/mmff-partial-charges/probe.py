"""Independent graph BFS and analytic buffered-Coulomb controls, not GAFF parity."""
import argparse
from collections import deque
import hashlib
import itertools
import json
import math
from pathlib import Path

import chematic
import chematic.chematic as extension


def distances(adjacency):
    result = []
    for start in range(len(adjacency)):
        row = [math.inf] * len(adjacency)
        row[start] = 0
        queue = deque([start])
        while queue:
            i = queue.popleft()
            for j, bonded in enumerate(adjacency[i]):
                if bonded and row[j] == math.inf:
                    row[j] = row[i] + 1
                    queue.append(j)
        result.append(row)
    return result


def analytic(xyz, charges, topology, cutoff, ignore_fragments, epsilon):
    energy = 0.0
    gradient = [[0.0] * 3 for _ in xyz]
    for i in range(len(xyz)):
        for j in range(i + 1, len(xyz)):
            relation = topology[i][j]
            if relation <= 2 or (ignore_fragments and relation == math.inf):
                continue
            delta = [a - b for a, b in zip(xyz[i], xyz[j])]
            radius = math.hypot(*delta)
            if radius > cutoff:
                continue
            coefficient = 332.0716 * charges[i] * charges[j] / epsilon
            coefficient *= 0.75 if relation == 3 else 1.0
            energy += coefficient / (radius + 0.05)
            for k in range(3):
                value = -coefficient * delta[k] / (radius * (radius + 0.05)**2)
                gradient[i][k] += value
                gradient[j][k] -= value
    return energy, gradient


def flatten(rows):
    return [value for row in rows for value in row]


def main(output):
    fixture = Path(__file__).parents[1] / "mmff-dielectric/cases.json"
    fixtures = json.loads(fixture.read_text())["cases"]
    rows = []
    for case in fixtures:
        mol = chematic.from_smiles(case["smiles"]).add_hydrogens()
        xyz = case["coords"]
        topology = distances(mol.rdkit_adjacency_matrix())
        before = mol.smiles, mol.formula, mol.GetPropsAsDict()
        for variant, cutoff, ignore, model in itertools.product(
            ["MMFF94", "MMFF94s"], [3.0, 100.0], [False, True], range(3)
        ):
            charges = (
                [0.0] * len(xyz) if model == 0 else
                [0.2] * len(xyz) if model == 1 else
                [(-1)**i * (i + 1) * 0.03 for i in range(len(xyz))]
            )
            options = dict(variant=variant, non_bonded_thresh=cutoff,
                           ignore_interfrag_interactions=ignore)
            custom = mol.mmff_force_field(xyz, **options, partial_charges=charges)
            zero = mol.mmff_force_field(xyz, **options, partial_charges=[0.0]*len(xyz))
            zero_energy, zero_gradient = zero.energy_gradient()
            for epsilon in [1.0, 4.0, 78.0]:
                custom.set_dielectric_constant_(epsilon)
                expected_energy, expected_gradient = analytic(
                    xyz, charges, topology, cutoff, ignore, epsilon
                )
                energy, gradient = custom.energy_gradient()
                measured = energy - zero_energy
                measured_gradient = [a-b for a, b in zip(flatten(gradient), flatten(zero_gradient))]
                energy_error = abs(measured - expected_energy)
                gradient_error = max(abs(a-b) for a, b in zip(measured_gradient, flatten(expected_gradient)))
                assert energy_error <= 1e-8, (case["smiles"], energy_error)
                assert gradient_error <= 1e-8, (case["smiles"], gradient_error)
                rows.append(dict(smiles=case["smiles"], variant=variant,
                                 cutoff=cutoff, ignore_interfragment=ignore,
                                 charges=charges, epsilon=epsilon,
                                 measured_electrostatics=measured,
                                 analytic_electrostatics=expected_energy,
                                 energy_absolute_error=energy_error,
                                 gradient_max_absolute_error=gradient_error))
            assert (mol.smiles, mol.formula, mol.GetPropsAsDict()) == before

    # Differentiate the full hybrid model, including all unchanged MMFF terms.
    case = fixtures[2]
    mol = chematic.from_smiles(case["smiles"]).add_hydrogens()
    charges = [(-1)**i * (i + 1) * 0.03 for i in range(mol.atom_count)]
    field = mol.mmff_force_field(case["coords"], partial_charges=charges).set_dielectric_constant_(4)
    _, gradient = field.energy_gradient()
    finite_difference_errors = []
    step = 1e-5
    for atom, axis in itertools.product(range(mol.atom_count), range(3)):
        plus = [list(row) for row in case["coords"]]
        minus = [list(row) for row in case["coords"]]
        plus[atom][axis] += step
        minus[atom][axis] -= step
        field.set_positions_(plus)
        upper = field.energy()
        field.set_positions_(minus)
        lower = field.energy()
        error = abs((upper-lower)/(2*step) - gradient[atom][axis])
        assert error <= 1e-5, (atom, axis, error)
        finite_difference_errors.append(error)

    result = dict(
        native_path=extension.__file__,
        native_sha256=hashlib.sha256(Path(extension.__file__).read_bytes()).hexdigest(),
        fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest(),
        analytic_cases=len(rows), failures=0,
        max_energy_absolute_error=max(row["energy_absolute_error"] for row in rows),
        max_gradient_absolute_error=max(row["gradient_max_absolute_error"] for row in rows),
        analytic_tolerance=1e-8, finite_difference_axes=len(finite_difference_errors),
        max_finite_difference_error=max(finite_difference_errors),
        finite_difference_tolerance=1e-5, rows=rows,
    )
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "rows"}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    main(parser.parse_args().output)
