"""RDKit 2026.03.1 reference values for the count-style RDKit fingerprints."""

import chematic


def _bits(raw: bytes) -> list[int]:
    return [i for i in range(len(raw) * 8) if raw[i // 8] >> (i % 8) & 1]


def test_morgan_count_simulation_and_chirality():
    mol = chematic.from_smiles("F/C=C/C[C@H](N)O")
    assert _bits(mol.rdkit_ecfp_config(2, 2048, include_chirality=True)) == [
        1, 53, 80, 227, 694, 724, 783, 786, 807, 894, 1002, 1157, 1171, 1589, 1649, 1928, 1938]
    aspirin = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
    assert _bits(aspirin.rdkit_ecfp_config(2, 2048, count_simulation=True))[:6] == [44, 92, 132, 256, 257, 320]


def test_hashed_count_fingerprints():
    aspirin = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
    assert aspirin.rdkit_atom_pair_counts(2048)[:3] == [(4, 1), (71, 2), (74, 1)]
    assert aspirin.rdkit_torsion_counts(2048)[:3] == [(7, 1), (94, 3), (102, 1)]
