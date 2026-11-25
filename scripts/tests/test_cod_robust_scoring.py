import sys
from pathlib import Path
from unittest.mock import patch

from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks'))
from cod_robust_scoring import annotate
from platinum_geometry import prepare_geometry_checks


def test_reverse_rmsd_match_does_not_skip_stereo_or_geometry():
    mol = Chem.AddHs(Chem.MolFromSmiles('CCO'))
    assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0
    reference = Chem.RemoveHs(Chem.Mol(mol))
    points = [list(mol.GetConformer().GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    with patch('cod_robust_scoring.rms_values', side_effect=RuntimeError('no match')):
        result = annotate({'candidates': [(0, points)]}, mol, reference, prepare_geometry_checks(mol))
    assert result['rmsd'][0] < 1e-5
    assert result['rmsd_errors'] == []
    assert len(result['stereo']) == len(result['geometry']) == 1
