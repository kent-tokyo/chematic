"""External RDKit geometry checks for strict benchmark stereochemistry."""
import sys
import unittest
from unittest.mock import patch
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks'))
from platinum_stereo_strict import stereo_violations_strict  # noqa: E402
from summarize_platinum_stereo_strict import summarize  # noqa: E402


class StrictStereoTest(unittest.TestCase):
    def test_tetrahedral_and_mirror(self):
        mol = Chem.AddHs(Chem.MolFromSmiles('F[C@](Cl)(Br)I'))
        self.assertEqual(0, AllChem.EmbedMolecule(mol, randomSeed=42))
        conf = mol.GetConformer()
        coords = [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
        mirrored = [[-x, y, z] for x, y, z in coords]
        good, bad = stereo_violations_strict(mol, [coords, mirrored])
        self.assertEqual([], good['atoms'])
        self.assertEqual([1], bad['violated_atoms'])
        self.assertEqual([], bad['unevaluable_atoms'])

    def test_square_planar_permutation_and_ambiguous_ligands(self):
        mol = Chem.MolFromSmiles('[Pt@SP1](F)(Cl)(Br)I')
        good = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, -1, 0)]
        bad = [(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]
        bad_third = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, -1, 0), (-1, 0, 0)]
        a, b, c = stereo_violations_strict(mol, [good, bad, bad_third])
        self.assertEqual([], a['atoms'])
        self.assertEqual([0], b['violated_atoms'])
        self.assertEqual([0], c['violated_atoms'])
        ambiguous = Chem.MolFromSmiles('[Pt@SP1](F)(F)(Br)I')
        result = stereo_violations_strict(ambiguous, [good])[0]
        self.assertEqual([0], result['unevaluable_atoms'])
        self.assertEqual([0], result['atoms'])

    def test_declared_tetrahedral_centers_without_cip_are_checked(self):
        smiles = 'C[C@]1(C[NH2+]Cc2ccc[nH]2)CN(Cc2ccc[nH]2)[C@@H](c2ccc[nH]2)N(Cc2ccc[nH]2)C1'
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        center = mol.GetAtomWithIdx(1)
        self.assertFalse(center.HasProp('_CIPCode'))
        self.assertEqual(0, AllChem.EmbedMolecule(mol, randomSeed=42))
        conf = mol.GetConformer()
        coords = [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
        inverted = [p[:] for p in coords]
        neighbors = [atom.GetIdx() for atom in center.GetNeighbors()]
        inverted[neighbors[0]], inverted[neighbors[1]] = inverted[neighbors[1]], inverted[neighbors[0]]
        good, bad = stereo_violations_strict(mol, [coords, inverted])
        self.assertNotIn(1, good['atoms'])
        self.assertIn(1, bad['violated_atoms'])

    def test_unrecognized_raw_tetrahedral_tags_fail_closed_or_check_parity(self):
        mol = Chem.AddHs(Chem.MolFromSmiles('F[C@](Cl)(Br)I'))
        self.assertEqual(0, AllChem.EmbedMolecule(mol, randomSeed=42))
        conf = mol.GetConformer()
        coords = [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
        mirrored = [[-x, y, z] for x, y, z in coords]
        flattened = [[x, y, 0.0] for x, y, _ in coords]
        # Emulate a declared raw tag omitted by RDKit's potential-stereo list.
        with patch('platinum_stereo_strict.Chem.FindPotentialStereo', return_value=[]):
            good, inverted, flat = stereo_violations_strict(mol, [coords, mirrored, flattened])
        self.assertEqual([], good['atoms'])
        self.assertEqual([1], inverted['violated_atoms'])
        self.assertEqual([1], flat['unevaluable_atoms'])
        equivalent = Chem.MolFromSmiles('FC(F)(Cl)Br')
        equivalent.GetAtomWithIdx(1).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        with patch('platinum_stereo_strict.Chem.FindPotentialStereo', return_value=[]):
            result = stereo_violations_strict(equivalent, [coords])[0]
        self.assertEqual([1], result['unevaluable_atoms'])

    def test_tetrahedral_phosphorus_without_cip(self):
        smiles = 'Cc1cc(C)c(N2[P@H](=S)N(c3c(C)cc(C)cc3C)[P@H]2=S)c(C)c1'
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        centers = [a.GetIdx() for a in mol.GetAtoms()
                   if a.GetSymbol() == 'P' and a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED]
        self.assertEqual(2, len(centers))
        self.assertTrue(all(not mol.GetAtomWithIdx(i).HasProp('_CIPCode') for i in centers))
        self.assertEqual(0, AllChem.EmbedMolecule(mol, randomSeed=42))
        conf = mol.GetConformer()
        coords = [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
        mirrored = [[-x, y, z] for x, y, z in coords]
        good, bad = stereo_violations_strict(mol, [coords, mirrored])
        self.assertTrue(all(i not in good['atoms'] for i in centers))
        self.assertTrue(all(i in bad['violated_atoms'] for i in centers))

    def test_invalid_coordinates_fail_closed(self):
        mol = Chem.MolFromSmiles('F[C@](Cl)(Br)I')
        with self.assertRaises(ValueError):
            stereo_violations_strict(mol, [[(0.0, 0.0, 0.0)]])
        invalid = [(0.0, 0.0, 0.0)] * mol.GetNumAtoms()
        invalid[0] = (float('nan'), 0.0, 0.0)
        with self.assertRaises(ValueError):
            stereo_violations_strict(mol, [invalid])

    def test_summary_separates_unknown_from_contradiction(self):
        violated = {'atoms': [1], 'bonds': [], 'violated_atoms': [1],
                    'violated_bonds': [], 'unevaluable_atoms': [], 'unevaluable_bonds': []}
        unknown = {'atoms': [2], 'bonds': [], 'violated_atoms': [],
                   'violated_bonds': [], 'unevaluable_atoms': [2], 'unevaluable_bonds': []}
        data = {'stereo_policy': 'rdkit-declared-stereo-strict-v1', 'results': [
            {'chematic_stereo': [violated, unknown], 'rdkit_stereo': []}]}
        counts = summarize(data)['engines']['chematic']
        self.assertEqual(2, counts['conformers'])
        self.assertEqual(1, counts['contradicted'])
        self.assertEqual(1, counts['unevaluable'])
        self.assertEqual(2, counts['invalid_either'])
        unknown['atoms'] = []
        with self.assertRaises(ValueError):
            summarize(data)

    def test_double_bond_inversion(self):
        mol = Chem.AddHs(Chem.MolFromSmiles('F/C=C/F'))
        self.assertEqual(0, AllChem.EmbedMolecule(mol, randomSeed=42))
        conf = mol.GetConformer()
        coords = [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
        fluorine = 0
        carbon = mol.GetAtomWithIdx(1)
        hydrogen = next(a.GetIdx() for a in carbon.GetNeighbors() if a.GetSymbol() == 'H')
        inverted = [p[:] for p in coords]
        inverted[fluorine], inverted[hydrogen] = inverted[hydrogen], inverted[fluorine]
        good, bad = stereo_violations_strict(mol, [coords, inverted])
        self.assertEqual([], good['bonds'])
        self.assertTrue(bad['bonds'])


if __name__ == '__main__':
    unittest.main()
