"""Depictions show the declared stereo; SDWriter writes Å coordinates (batch 15)."""

import chematic


def _same_side(p, a, b, q):
    def side(r):
        return (b[0] - a[0]) * (r[1] - a[1]) - (b[1] - a[1]) * (r[0] - a[0])

    return side(p) * side(q) > 0


def _points(data):
    return [(a["x"], a["y"]) for a in data["atoms"]]


def test_cis_and_trans_are_drawn_as_declared():
    cis = _points(chematic.from_smiles("F/C=C\\F").depict_data())
    trans = _points(chematic.from_smiles("F/C=C/F").depict_data())
    assert _same_side(cis[0], cis[1], cis[2], cis[3])
    assert not _same_side(trans[0], trans[1], trans[2], trans[3])


def test_smiles_bond_directions_are_not_wedges():
    kinds = [b["kind"] for b in chematic.from_smiles("C/C=C/C").depict_data()["bonds"]]
    assert kinds == ["Single", "Double", "Single"]


def test_stereocentre_gets_one_wedge_from_the_centre():
    data = chematic.from_smiles("N[C@@H](C)C(=O)O").depict_data()
    wedges = [b for b in data["bonds"] if b["kind"] in ("Up", "Down")]
    assert len(wedges) == 1 and wedges[0]["atom1"] == 1
    assert "<polygon" in chematic.from_smiles("N[C@@H](C)C(=O)O").svg()


def test_sdwriter_writes_angstrom_coordinates_and_ez(tmp_path):
    path = tmp_path / "out.sdf"
    with chematic.SDWriter(str(path)) as w:
        w.write(chematic.from_smiles("F/C=C\\F"))
    block = path.read_text()
    atom_lines = block.splitlines()[4:8]
    xs = [float(line.split()[0]) for line in atom_lines]
    ys = [float(line.split()[1]) for line in atom_lines]
    bond = ((xs[1] - xs[2]) ** 2 + (ys[1] - ys[2]) ** 2) ** 0.5
    assert abs(bond - 1.5) < 1e-3
    back = chematic.from_mol_block(block.split("$$$$")[0])
    assert back.smiles == chematic.from_smiles("F/C=C\\F").smiles


def test_mol_block_without_stereo_is_laid_out():
    # Like RDKit's MolToMolBlock (and the WASM binding): a molecule without
    # stereo used to be written with every atom at the origin.
    block = chematic.from_smiles("CCO").to_mol_block()
    lines = block.splitlines()[4:7]
    xs = [float(line.split()[0]) for line in lines]
    ys = [float(line.split()[1]) for line in lines]
    bond = ((xs[0] - xs[1]) ** 2 + (ys[0] - ys[1]) ** 2) ** 0.5
    assert abs(bond - 1.5) < 1e-3
    report_block, _ = chematic.from_smiles("CCO").to_mol_block_with_report()
    assert report_block == block


def test_reaction_and_similarity_map_svgs_draw_wedges():
    # A stereocentre is drawn with its wedge (a filled polygon) in the
    # reaction and similarity-map SVGs too.
    rxn = chematic.reaction_svg("N[C@@H](C)C(=O)O>>N[C@@H](C)C(=O)OC")
    assert "<polygon" in rxn
    mol = chematic.from_smiles("N[C@@H](C)C(=O)O")
    assert "<polygon" in chematic.similarity_map_svg(mol, [0.1] * mol.heavy_atoms)
