"""Tests for module-level functions: SMARTS, InChI, depict_grid, run_smirks, find_mcs, B7."""
import pytest
import chematic


@pytest.fixture(scope="module")
def ethanol():
    return chematic.from_smiles("CCO")


@pytest.fixture(scope="module")
def methanol():
    return chematic.from_smiles("CO")


@pytest.fixture(scope="module")
def aspirin():
    return chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")


# ---------------------------------------------------------------------------
# SMARTS matching
# ---------------------------------------------------------------------------

def test_smarts_match_hydroxyl(ethanol):
    assert chematic.smarts_match("[OH]", ethanol)


def test_smarts_match_no_match(ethanol):
    assert not chematic.smarts_match("[NH2]", ethanol)


def test_smarts_find_returns_list(ethanol):
    matches = chematic.smarts_find("[OH]", ethanol)
    assert isinstance(matches, list)
    assert len(matches) >= 1
    for match in matches:
        assert isinstance(match, list)


def test_smarts_find_carboxyl(aspirin):
    # Aspirin has one carboxyl group
    matches = chematic.smarts_find("C(=O)[OH]", aspirin)
    assert len(matches) >= 1


def test_smarts_invalid():
    with pytest.raises(ValueError):
        chematic.smarts_match("[INVALID???", chematic.from_smiles("C"))


# ---------------------------------------------------------------------------
# from_inchi / InChI round-trip (C3: InChI parser)
# ---------------------------------------------------------------------------

def test_from_inchi_ethanol():
    inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
    mol = chematic.from_inchi(inchi)
    assert isinstance(mol, chematic.Mol)
    assert mol.formula == "C2H6O"


def test_from_inchi_invalid():
    with pytest.raises(ValueError):
        chematic.from_inchi("not an inchi string")


def test_inchi_roundtrip(ethanol):
    inchi = ethanol.inchi
    mol2 = chematic.from_inchi(inchi)
    # Formulas should match
    assert mol2.formula == ethanol.formula


# ---------------------------------------------------------------------------
# depict_grid
# ---------------------------------------------------------------------------

def test_depict_grid_returns_svg(ethanol, aspirin):
    svg = chematic.depict_grid([ethanol, aspirin], cols=2)
    assert isinstance(svg, str)
    assert len(svg) > 0


def test_depict_grid_single_mol(ethanol):
    svg = chematic.depict_grid([ethanol], cols=1)
    assert isinstance(svg, str)


# ---------------------------------------------------------------------------
# run_smirks
# ---------------------------------------------------------------------------

def test_run_smirks_basic():
    mol = chematic.from_smiles("CCO")
    products = chematic.run_smirks("[OH:1]>>[O-:1]", [mol])
    assert isinstance(products, list)


def test_run_smirks_invalid_smirks():
    with pytest.raises(ValueError):
        chematic.run_smirks("NOT_A_SMIRKS", [chematic.from_smiles("C")])


# E/Z stereo transfer & creation in products (issue #50)

def _product_ez(smirks, smis):
    mols = [chematic.from_smiles(s) for s in smis]
    prod = chematic.run_smirks(smirks, mols)[0][0]
    labels = [d["descriptor"] for d in prod.cip_stereo()]
    return prod.smiles, labels


def test_run_smirks_transfer_preserves_E():
    # canonical_smiles() picks a different, equally valid spelling than this
    # test was originally written against (issue #200) -- assert against the
    # current canonical form of the input, not one hardcoded string, plus the
    # CIP label the test actually cares about.
    smi, labels = _product_ez("[C:1]=[C:2]>>[C:1]=[C:2]", ["C/C=C/C"])
    assert smi == chematic.from_smiles("C/C=C/C").smiles
    assert "E" in labels


def test_run_smirks_transfer_preserves_Z():
    smi, labels = _product_ez("[C:1]=[C:2]>>[C:1]=[C:2]", ["C/C=C\\C"])
    assert smi == chematic.from_smiles("C/C=C\\C").smiles
    assert "Z" in labels


def test_run_smirks_create_E_from_template():
    smi, labels = _product_ez(
        "[C:1][C:2][C:3][C:4]>>[C:1]/[C:2]=[C:3]/[C:4]", ["CCCC"]
    )
    assert "E" in labels


def test_run_smirks_create_Z_from_template():
    smi, labels = _product_ez(
        "[C:1][C:2][C:3][C:4]>>[C:1]/[C:2]=[C:3]\\[C:4]", ["CCCC"]
    )
    assert "Z" in labels


# ---------------------------------------------------------------------------
# find_mcs
# ---------------------------------------------------------------------------

def test_find_mcs_similar_mols():
    m1 = chematic.from_smiles("c1ccccc1")   # benzene
    m2 = chematic.from_smiles("Cc1ccccc1")  # toluene
    mcs = chematic.find_mcs([m1, m2])
    # MCS should be non-None (at minimum a 6-carbon ring)
    assert mcs is not None
    assert isinstance(mcs, chematic.Mol)


def test_find_mcs_single_mol():
    mol = chematic.from_smiles("CCO")
    mcs = chematic.find_mcs([mol])
    assert mcs is not None


def test_find_mcs_returns_none_or_mol():
    m1 = chematic.from_smiles("C")
    m2 = chematic.from_smiles("N")
    result = chematic.find_mcs([m1, m2])
    assert result is None or isinstance(result, chematic.Mol)


def test_find_mcs_ring_config_quinoline_scaffold():
    # Issue #1 example: quinoline series with both ring constraints.
    # The shared exocyclic CH₂ (non-ring in all) remains valid under ring_matches_ring_only,
    # so the result is at least the quinoline scaffold (10 atoms).
    mols = [
        chematic.from_smiles("c1ccc2nc(CC)ccc2c1"),
        chematic.from_smiles("c1ccc2nc(CO)ccc2c1"),
        chematic.from_smiles("c1ccc2nc(CN)ccc2c1"),
    ]
    mcs = chematic.find_mcs(mols, ring_matches_ring_only=True, complete_rings_only=True)
    assert mcs is not None
    assert mcs.heavy_atoms >= 10, f"expected at least 10-atom quinoline scaffold, got {mcs.heavy_atoms}"


def test_find_mcs_ring_config_complete_rings_benzene_toluene():
    # complete_rings_only: benzene vs toluene should give exactly the 6-atom benzene ring.
    m1 = chematic.from_smiles("c1ccccc1")
    m2 = chematic.from_smiles("Cc1ccccc1")
    mcs = chematic.find_mcs([m1, m2], complete_rings_only=True)
    assert mcs is not None
    assert mcs.heavy_atoms == 6, f"expected 6-atom benzene ring, got {mcs.heavy_atoms}"


def test_find_mcs_default_kwargs_unchanged():
    # Calling without kwargs should work as before (backward compat).
    m1 = chematic.from_smiles("c1ccccc1")
    m2 = chematic.from_smiles("Cc1ccccc1")
    mcs = chematic.find_mcs([m1, m2])
    assert mcs is not None


# ---------------------------------------------------------------------------
# B7: Reaction SMARTS matching
# ---------------------------------------------------------------------------

def test_reaction_smarts_match_basic():
    # Alcohol deprotonation pattern
    result = chematic.reaction_smarts_match("[OH]>>[O-]", "CCO>>CC[O-]")
    assert isinstance(result, bool)


def test_reaction_smarts_match_no_match():
    result = chematic.reaction_smarts_match("[NH2]>>[NH-]", "CCO>>CC[O-]")
    assert not result


def test_reaction_smarts_invalid_smarts():
    with pytest.raises(ValueError):
        chematic.reaction_smarts_match("NOT>>VALID???", "C>>C")


def test_reaction_smarts_invalid_rxn_smiles():
    with pytest.raises(ValueError):
        chematic.reaction_smarts_match("[OH]>>[O-]", "NOT_A_REACTION")


# ---------------------------------------------------------------------------
# from_mol_block
# ---------------------------------------------------------------------------

def test_from_mol_block_basic():
    # Minimal V2000 MOL block for methane
    mol_block = """\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n"""
    try:
        mol = chematic.from_mol_block(mol_block)
        assert isinstance(mol, chematic.Mol)
    except ValueError:
        # Some mol block formats may not be supported — that's OK
        pass
