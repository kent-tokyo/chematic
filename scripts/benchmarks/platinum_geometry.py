"""Independent geometric screens for explicit-hydrogen conformers in Angstrom."""

import itertools
import math

from rdkit import Chem
from rdkit.Chem import AllChem


def prepare_geometry_checks(mol):
    """Prepare supported UFF reference bonds, angles, and heavy-atom clash pairs."""
    bonds = []
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        params = AllChem.GetUFFBondStretchParams(mol, i, j)
        if params is not None:
            bonds.append((i, j, params[1]))
    angles = []
    small_rings = [
        set(ring) for ring in mol.GetRingInfo().AtomRings() if len(ring) in (3, 4)
    ]
    for atom in mol.GetAtoms():
        center = atom.GetIdx()
        neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
        for i, k in itertools.combinations(neighbors, 2):
            params = AllChem.GetUFFAngleBendParams(mol, i, center, k)
            if params is not None:
                ring_sizes = [
                    len(ring) for ring in small_rings if {i, center, k}.issubset(ring)
                ]
                reference = (
                    {3: 60.0, 4: 90.0}[min(ring_sizes)] if ring_sizes else params[1]
                )
                angles.append((i, center, k, reference))
    distances = Chem.GetDistanceMatrix(mol)
    table = Chem.GetPeriodicTable()
    heavy = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    clashes = []
    for first, second in itertools.combinations(heavy, 2):
        i, j = first.GetIdx(), second.GetIdx()
        if distances[i, j] > 2:
            threshold = 0.65 * (
                table.GetRvdw(first.GetAtomicNum())
                + table.GetRvdw(second.GetAtomicNum())
            )
            clashes.append((i, j, threshold))
    return bonds, angles, clashes


def check_geometry(coords, checks):
    bonds, angles, clashes = checks
    if not all(math.isfinite(value) for point in coords for value in point):
        return {"finite": False}

    def distance(i, j):
        return math.dist(coords[i], coords[j])

    bad_bonds = sum(
        abs(distance(i, j) / reference - 1) > 0.5 for i, j, reference in bonds
    )
    bad_angles = 0
    for i, j, k, reference in angles:
        a = [coords[i][axis] - coords[j][axis] for axis in range(3)]
        b = [coords[k][axis] - coords[j][axis] for axis in range(3)]
        norm = math.sqrt(sum(v * v for v in a) * sum(v * v for v in b))
        if norm == 0:
            bad_angles += 1
            continue
        cosine = max(-1.0, min(1.0, sum(x * y for x, y in zip(a, b)) / norm))
        bad_angles += abs(math.degrees(math.acos(cosine)) - reference) > 35.0
    bad_clashes = sum(distance(i, j) < threshold for i, j, threshold in clashes)
    return {
        "finite": True,
        "bonds_checked": len(bonds),
        "bond_violations": bad_bonds,
        "angles_checked": len(angles),
        "angle_violations": bad_angles,
        "heavy_pairs_checked": len(clashes),
        "heavy_atom_clashes": bad_clashes,
    }
