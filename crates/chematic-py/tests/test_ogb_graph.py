"""OGB graph export contract at the Python boundary."""

import numpy as np
import pytest

import chematic


@pytest.mark.parametrize("smiles", ["CCO", "c1ccccc1", "[NH4+]", "C/C=C/C"])
def test_ogb_graph_shapes_and_edge_alignment(smiles):
    mol = chematic.from_smiles(smiles)
    graph = mol.ogb_graph_features()
    atom_count = len(mol.atom_table)
    bond_count = len(mol.bond_table)
    assert set(graph) == {"node_feat", "edge_index", "edge_feat", "num_nodes"}
    assert graph["num_nodes"] == atom_count
    assert graph["node_feat"].shape == (atom_count, 9)
    assert graph["edge_index"].shape == (2, 2 * bond_count)
    assert graph["edge_feat"].shape == (2 * bond_count, 3)
    for name in ("node_feat", "edge_index", "edge_feat"):
        assert graph[name].dtype == np.int64
    for pair in range(0, graph["edge_index"].shape[1], 2):
        np.testing.assert_array_equal(
            graph["edge_index"][:, pair], graph["edge_index"][::-1, pair + 1]
        )
        np.testing.assert_array_equal(graph["edge_feat"][pair], graph["edge_feat"][pair + 1])


def test_ogb_graph_from_sdf_record():
    source = chematic.from_smiles("CCO")
    sdf = source.to_mol_block() + "\n$$$$\n"
    record = next(iter(chematic.iter_sdf_str(sdf)))
    graph = record.mol.ogb_graph_features()
    assert graph["node_feat"].shape == (3, 9)
    assert graph["edge_index"].shape == (2, 4)


def test_sdf_record_retains_explicitly_unknown_double_bond_stereo():
    block = chematic.from_smiles("CC=CC").to_mol_block()
    assert "  2  3  2  3" in block  # V2000 crossed/either double bond
    record = next(iter(chematic.iter_sdf_str(block + "\n$$$$\n")))
    graph = record.ogb_graph_features()
    np.testing.assert_array_equal(graph["edge_feat"][2:4, 1], [5, 5])
    # A bare Mol has no source-record diagnostic for this CTAB marker.
    np.testing.assert_array_equal(record.mol.ogb_graph_features()["edge_feat"][2:4, 1], [0, 0])


def test_v3000_sdf_record_retains_unknown_double_bond_stereo(tmp_path):
    block = chematic.from_smiles("CC=CC").to_mol_v3000([])
    assert "M  V30 2 2 2 3" in block
    block = block.replace("M  V30 2 2 2 3", "M  V30 2 2 2 3 CFG=2")
    path = tmp_path / "either_v3000.sdf"
    path.write_text(block + "\n> <Activity>\n7.1\n\n$$$$\n")
    index, record, error = next(iter(chematic.iter_sdf_checked(str(path))))
    assert index == 0 and error is None and record is not None
    assert record.get("Activity") == "7.1"
    np.testing.assert_array_equal(record.ogb_graph_features()["edge_feat"][2:4, 1], [5, 5])
    indexed = chematic.SdfDataset.open(str(path))[0]
    assert indexed.get("Activity") == "7.1"
    assert indexed.ogb_graph_features()["edge_feat"][2, 1] == 5
    assert next(iter(chematic.iter_sdf_str(path.read_text(encoding="utf-8")))).get("Activity") == "7.1"
    assert next(iter(chematic.SDMolSupplier(str(path)))) is not None


def test_sdf_checked_stream_retains_rejected_record_index(tmp_path):
    good = chematic.from_smiles("CCO").to_mol_block()
    path = tmp_path / "mixed.sdf"
    path.write_text("invalid MOL block\n$$$$\n" + good + "\n$$$$\n")
    rejected, accepted = list(chematic.iter_sdf_checked(str(path)))
    assert rejected[0] == 0
    assert rejected[1] is None
    assert isinstance(rejected[2], str) and rejected[2]
    assert accepted[0] == 1
    assert accepted[2] is None
    assert accepted[1].mol.ogb_graph_features()["node_feat"].shape == (3, 9)


def test_sdf_checked_stream_recovers_after_non_utf8_record(tmp_path):
    good = chematic.from_smiles("CCO").to_mol_block().encode()
    path = tmp_path / "mixed_bytes.sdf"
    path.write_bytes(b"\xff\n$$$$\n" + good + b"\n$$$$\n")
    rejected, accepted = list(chematic.iter_sdf_checked(str(path)))
    assert rejected[0] == 0 and rejected[1] is None and rejected[2]
    assert accepted[0] == 1 and accepted[1] is not None and accepted[2] is None


def test_sdf_checked_stream_limit_is_reported(tmp_path):
    good = chematic.from_smiles("CCO").to_mol_block()
    path = tmp_path / "two.sdf"
    path.write_text(good + "\n$$$$\n" + good + "\n$$$$\n")
    first, limit = list(chematic.iter_sdf_checked(str(path), max_records=1))
    assert first[0] == 0 and first[1] is not None
    assert limit[0] == 1 and limit[1] is None and "records exceeds limit" in limit[2]
    assert all(record is not None for _, record, _ in chematic.iter_sdf_checked(str(path), max_records=2))
    with pytest.raises(ValueError, match="max_records"):
        chematic.iter_sdf_checked(str(path), max_records=0)


def test_ogb_graph_from_mol2():
    source = chematic.from_smiles("CCO")
    mol = chematic.from_mol2(source.to_mol2())
    np.testing.assert_array_equal(
        mol.ogb_graph_features()["node_feat"], source.ogb_graph_features()["node_feat"]
    )


def test_ogb_graph_hydrogen_policy():
    mol = chematic.from_smiles("[H]C([H])([H])[H]")
    default = mol.ogb_graph_features()
    assert default["node_feat"].shape == (1, 9)
    assert default["edge_index"].shape == (2, 0)
    explicit = mol.ogb_graph_features(remove_hs=False)
    assert explicit["node_feat"].shape == (5, 9)
    assert explicit["edge_index"].shape == (2, 8)


def test_ogb_graph_from_sdf_stereo():
    Chem = pytest.importorskip("rdkit.Chem")
    chiral = Chem.MolToMolBlock(Chem.MolFromSmiles("F[C@](Cl)(Br)I"))
    alkene = Chem.MolToMolBlock(Chem.MolFromSmiles("C/C=C/C"))
    assert chematic.from_mol_block(chiral).ogb_graph_features()["node_feat"][1, 1] == 2
    assert chematic.from_mol_block(alkene).ogb_graph_features()["edge_feat"][2, 1] == 2


def test_ogb_graph_matches_official_ogb_when_available():
    ogb = pytest.importorskip("ogb.utils.mol")
    for smiles in ("CCO", "c1ccccc1", "[NH4+]", "[CH3]", "F[C@](Cl)(Br)I", "C/C=C/C", "[H]C([H])([H])[H]"):
        actual = chematic.from_smiles(smiles).ogb_graph_features()
        expected = ogb.smiles2graph(smiles)
        assert actual["num_nodes"] == expected["num_nodes"]
        for name in ("node_feat", "edge_index", "edge_feat"):
            np.testing.assert_array_equal(actual[name], expected[name])


def test_ogb_graph_matches_rdkit_oracle():
    """Compare the OGB field definitions with live RDKit chemistry values."""
    Chem = pytest.importorskip("rdkit.Chem")
    # Categorical order from ogb/ogb/utils/features.py.
    chirality = ["CHI_UNSPECIFIED", "CHI_TETRAHEDRAL_CW", "CHI_TETRAHEDRAL_CCW", "CHI_OTHER", "misc"]
    hybridization = ["SP", "SP2", "SP3", "SP3D", "SP3D2", "misc"]
    bond_type = ["SINGLE", "DOUBLE", "TRIPLE", "AROMATIC", "misc"]
    bond_stereo = ["STEREONONE", "STEREOZ", "STEREOE", "STEREOCIS", "STEREOTRANS", "STEREOANY"]

    def index(allowed, value):
        return allowed.index(value) if value in allowed else len(allowed) - 1

    samples = (
        "CCO", "c1ccccc1", "[NH4+]", "[CH3]", "F[C@](Cl)(Br)I",
        "C/C=C/C", r"C/C=C\C", "CC(=O)N", "CC(=O)[O-]",
        "[N+](=O)[O-]", "c1ccncc1", "c1ccoc1", "C#N", "[Na+].[Cl-]",
        "C1CCCCC1", "CC(=O)Oc1ccccc1C(=O)O", "[13CH4]",
        "O=P(O)(O)O", "CS(=O)(=O)C", "[H]C([H])([H])[H]",
    )
    for smiles in samples:
        reference = Chem.MolFromSmiles(smiles)
        assert reference is not None
        actual = chematic.from_smiles(smiles).ogb_graph_features()
        nodes = []
        for atom in reference.GetAtoms():
            nodes.append([
                index([*range(1, 119), "misc"], atom.GetAtomicNum()),
                index(chirality, str(atom.GetChiralTag())),
                index([*range(11), "misc"], atom.GetTotalDegree()),
                index([*range(-5, 6), "misc"], atom.GetFormalCharge()),
                index([*range(9), "misc"], atom.GetTotalNumHs()),
                index([*range(5), "misc"], atom.GetNumRadicalElectrons()),
                index(hybridization, str(atom.GetHybridization())),
                int(atom.GetIsAromatic()),
                int(atom.IsInRing()),
            ])
        edges, features = [], []
        for bond in reference.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            feature = [
                index(bond_type, str(bond.GetBondType())),
                bond_stereo.index(str(bond.GetStereo())),
                int(bond.GetIsConjugated()),
            ]
            edges.extend(((i, j), (j, i)))
            features.extend((feature, feature))
        expected = {
            "node_feat": np.asarray(nodes, dtype=np.int64).reshape(-1, 9),
            "edge_index": np.asarray(edges, dtype=np.int64).reshape(-1, 2).T,
            "edge_feat": np.asarray(features, dtype=np.int64).reshape(-1, 3),
        }
        for field, values in expected.items():
            np.testing.assert_array_equal(actual[field], values, err_msg=f"{smiles}: {field}")
