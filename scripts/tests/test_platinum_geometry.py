import sys
import unittest
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
from platinum_geometry import check_geometry, prepare_geometry_checks


class PlatinumGeometryTest(unittest.TestCase):
    def test_bond_angle_and_heavy_clash_screens(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CCCC"))
        self.assertEqual(AllChem.EmbedMolecule(mol, randomSeed=42), 0)
        AllChem.UFFOptimizeMolecule(mol)
        checks = prepare_geometry_checks(mol)
        coords = [
            list(mol.GetConformer().GetAtomPosition(i))
            for i in range(mol.GetNumAtoms())
        ]
        normal = check_geometry(coords, checks)
        self.assertEqual(normal["bond_violations"], 0)
        self.assertEqual(normal["angle_violations"], 0)
        self.assertEqual(normal["heavy_atom_clashes"], 0)
        collapsed = [point[:] for point in coords]
        collapsed[3] = collapsed[0][:]
        broken = check_geometry(collapsed, checks)
        self.assertGreater(broken["heavy_atom_clashes"], 0)
        stretched = [point[:] for point in coords]
        stretched[0][0] += 10
        self.assertGreater(check_geometry(stretched, checks)["bond_violations"], 0)
        collinear = [point[:] for point in coords]
        collinear[0] = [2 * coords[1][axis] - coords[2][axis] for axis in range(3)]
        self.assertGreater(check_geometry(collinear, checks)["angle_violations"], 0)

    def test_nonfinite_coordinates_are_rejected(self):
        self.assertFalse(check_geometry([[float("nan"), 0, 0]], ([], [], []))["finite"])

    def test_small_ring_angle_reference(self):
        for smiles in ("C1CC1", "C1CCC1"):
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            self.assertEqual(AllChem.EmbedMolecule(mol, randomSeed=42), 0)
            AllChem.UFFOptimizeMolecule(mol)
            coords = [
                list(mol.GetConformer().GetAtomPosition(i))
                for i in range(mol.GetNumAtoms())
            ]
            self.assertEqual(
                check_geometry(coords, prepare_geometry_checks(mol))[
                    "angle_violations"
                ],
                0,
            )


if __name__ == "__main__":
    unittest.main()
