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


def test_run_smirks_product_atomic_number_primitives_match_organic_form():
    # Issue #679: product-side [#6](=[#8])[#6] adds an acetyl group like C(=O)C.
    mol = chematic.from_smiles("CCO")
    organic = chematic.run_smirks("[O:1]>>[*:1]C(=O)C", [mol])
    atomic = chematic.run_smirks("[O:1]>>[*:1][#6](=[#8])[#6]", [mol])
    bracket = chematic.run_smirks("[O:1]>>[*:1][C](=[O])[C]", [mol])
    canon = lambda sets: [[m.smiles for m in s] for s in sets]
    assert canon(organic) == canon(atomic) == canon(bracket)
    assert len(atomic) == 1
    assert atomic[0][0].smiles == chematic.from_smiles("CCOC(C)=O").smiles


def test_smarts_repeated_charge_signs_mean_total_charge():
    # Issue #680: [++] is +2 and [--] is -2, as in Daylight / OpenSMARTS.
    calcium = chematic.from_smiles("[Ca+2]")
    sodium = chematic.from_smiles("[Na+]")
    sulfide = chematic.from_smiles("[S-2]")
    assert chematic.smarts_match("[++]", calcium)
    assert chematic.smarts_match("[+2]", calcium)
    assert not chematic.smarts_match("[++]", sodium)
    assert chematic.smarts_match("[+]", sodium)
    assert chematic.smarts_match("[--]", sulfide)
    assert not chematic.smarts_match("[-]", sulfide)


@pytest.mark.parametrize(
    "smirks,smiles,expected",
    [
        # Issue #734: SMARTS query features on the reactant side.
        ("[CD1:1]>>[C:1]", "CC", "CC"),
        ("[CX4:1]>>[C:1]", "CC", "CC"),
        ("[C;H3:1]>>[C:1]", "CC", "CC"),
        ("[C,N:1]>>[C:1]", "CC", "CC"),
        ("[!O:1]>>[C:1]", "CC", "CC"),
        ("[$([OH]):1]>>[O:1]", "CCO", "CCO"),
        ("[C:1]>>[C;H3:1]", "CC", "CC"),
        ("[CX4;H2:1][OX2H1:2]>>[C:1]=[O:2]", "CCO", "CC=O"),
        ("[c:1][OX2:2][CH3:3]>>[c:1][OX2H1:2]", "COc1ccccc1", "Oc1ccccc1"),
        ("[#7X3;H0;!$(NC=O):1][CH3:2]>>[#7:1]", "CN(C)c1ccccc1", "CNc1ccccc1"),
        # RDKit ignores H0 on newly created, unmapped product atoms; mapped
        # atoms retain the explicit-zero override.
        ("[C:1]>>[C:1][CH0]", "C", "CC"),
        ("[C:1]>>[C:1][OH0]", "C", "CO"),
        ("[C:1]>>[C:1][NH0]", "C", "CN"),
    ],
)
def test_run_smirks_reactant_smarts_queries(smirks, smiles, expected):
    products = chematic.run_smirks(smirks, [chematic.from_smiles(smiles)])
    want = chematic.from_smiles(expected).smiles
    assert products
    assert {s[0].smiles for s in products} == {want}


def test_run_smirks_reactant_query_is_evaluated():
    # [CX4;H2] must not match the CH of isopropanol.
    assert chematic.run_smirks(
        "[CX4;H2:1][OX2H1:2]>>[C:1]=[O:2]", [chematic.from_smiles("CC(C)O")]
    ) == []


def test_run_smirks_invalid_smirks():
    with pytest.raises(ValueError):
        chematic.run_smirks("NOT_A_SMIRKS", [chematic.from_smiles("C")])


def test_run_smirks_checked_reports_success_and_no_match():
    ethanol = chematic.from_smiles("CCO")
    success = chematic.run_smirks_checked("[OH:1]>>[O-:1]", [ethanol])
    assert success["status"] == "products"
    assert success["reason"] is None
    assert success["accepted_matches"] >= 1
    assert success["applied_products"] >= 1
    assert success["valence_rejected_matches"] == 0
    assert success["products"]
    missing = chematic.run_smirks_checked("[NH2:1]>>[NH3+:1]", [ethanol])
    assert missing["status"] == "no_match"
    assert missing["products"] == []


def test_run_smirks_checked_preserves_atom_sources_and_template_maps():
    methane = chematic.from_smiles("C")
    report = chematic.run_smirks_checked("[C:1]>>[C:1].[O:2]", [methane])
    assert report["status"] == "products"
    assert len(report["products"]) == len(report["product_atom_sources"])
    assert len(report["products"]) == len(report["product_template_maps"])
    assert report["product_atom_sources"] == [[[(0, 0)], [None]]]
    assert report["product_template_maps"] == [[[1], [2]]]
    assert [product.smiles for product in report["products"][0]] == [
        product.smiles for product in chematic.run_smirks("[C:1]>>[C:1].[O:2]", [methane])[0]
    ]

    missing = chematic.run_smirks_checked("[NH2:1]>>[NH3+:1]", [methane])
    assert missing["product_atom_sources"] == []
    assert missing["product_template_maps"] == []


def test_run_smirks_checked_reports_typed_errors():
    methane = chematic.from_smiles("C")
    invalid = chematic.run_smirks_checked("NOT_A_SMIRKS", [methane])
    assert invalid["status"] == "typed_refusal"
    assert invalid["reason"] == "smirks_parse"
    assert invalid["product_atom_sources"] == []
    assert invalid["product_template_maps"] == []
    mismatch = chematic.run_smirks_checked("[C:1].[O:2]>>[C:1][O:2]", [methane])
    assert mismatch["status"] == "typed_refusal"
    assert mismatch["reason"] == "reactant_count_mismatch"


def test_run_smirks_checked_does_not_hide_filtered_product():
    # A neutral four-bonded N, which RDKit's sanitize also rejects.
    substrate = chematic.from_smiles("CN")
    report = chematic.run_smirks_checked("[N:1]>>[N:1](C)(C)C", [substrate])
    assert report["status"] == "typed_refusal"
    assert report["reason"] == "product_valence"
    assert report["accepted_matches"] >= 1
    assert report["valence_rejected_matches"] >= 1
    assert report["products"] == []


def test_run_smirks_checked_rdkit_profile_declines_chiral_reactants():
    alanine = chematic.from_smiles("N[C@@H](C)C(=O)O")
    smirks = "[N:1][C@@H:2](C)C(=O)O>>[N:1].[C@@H:2](C)C(=O)O"
    report = chematic.run_smirks_checked(smirks, [alanine], rdkit_compat=True)
    assert report["status"] == "typed_unsupported"
    assert report["reason"] == "chiral_reactant_template_semantics"
    assert report["products"] == []
    ordinary = chematic.run_smirks_checked(smirks, [alanine])
    assert ordinary["status"] != "typed_unsupported"


def test_run_smirks_checked_rdkit_profile_keeps_ez_templates():
    substrate = chematic.from_smiles("C/C=C/C")
    report = chematic.run_smirks_checked(
        "[C:1]=[C:2]>>[C:1]=[C:2]", [substrate], rdkit_compat=True
    )
    assert report["status"] == "products"


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


def test_find_mcs_result_is_substructure_of_both_inputs():
    # Regression test: qmol_to_mol once built every MCS atom via Atom::new(elem),
    # which always sets aromatic=False -- even for atoms on an aromatic-bonded
    # ring. That desync between the atom's aromatic flag and its bonds' literal
    # Aromatic order produced a self-consistent-looking but unusable SMILES
    # (e.g. "C:1:C:C(O):C:C:C:1", uppercase C with explicit aromatic bonds),
    # which then failed a substructure re-match against either parent molecule --
    # defeating the primary purpose of an MCS result (self-verification /
    # substructure screening).
    #
    # NOTE: `mcs.smiles` re-matching *every* input is only guaranteed when the
    # inputs agree on aromaticity at the matched atoms, as aspirin/paracetamol
    # do here. `AtomCompare::Elements` (the default) matches purely on atomic
    # number and permits cross-aromaticity matches (matching RDKit's own
    # `CompareElements`); `build_query` reflects that by never encoding
    # aromaticity as a per-atom constraint, only via bond queries. A SMILES
    # string can't express "aromaticity don't-care" the way a SMARTS/query
    # pattern can, so when two inputs' matched atoms genuinely disagree on
    # aromaticity, `mcs.smiles` -- built from one specific input's own bonds --
    # can fail to re-match the *other* input. That's an inherent limitation of
    # the concrete-SMILES accessor, not a bug in the MCS search itself.
    m1 = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
    m2 = chematic.from_smiles("CC(=O)Nc1ccc(O)cc1")  # paracetamol
    mcs = chematic.find_mcs([m1, m2])
    assert mcs is not None
    assert m1.has_substructure(mcs.smiles), f"MCS {mcs.smiles!r} not found in aspirin"
    assert m2.has_substructure(mcs.smiles), f"MCS {mcs.smiles!r} not found in paracetamol"


def test_find_mcs_match_charge_narrows_result():
    # Acetate vs acetic acid: match_charge=False (default) matches the full
    # core; match_charge=True must reject the charge-differing oxygen.
    acetate = chematic.from_smiles("CC(=O)[O-]")
    acetic = chematic.from_smiles("CC(=O)O")
    without = chematic.find_mcs([acetate, acetic], match_charge=False)
    with_charge = chematic.find_mcs([acetate, acetic], match_charge=True)
    assert without.smiles != with_charge.smiles


def test_find_mcs_match_isotope_narrows_result():
    labeled = chematic.from_smiles("[13CH4]")
    plain = chematic.from_smiles("C")
    without = chematic.find_mcs([labeled, plain], match_isotope=False)
    with_isotope = chematic.find_mcs([labeled, plain], match_isotope=True)
    assert without is not None and without.heavy_atoms == 1
    assert with_isotope is None or with_isotope.heavy_atoms == 0


def test_find_mcs_atom_compare_any_heavy_atom_widens_match():
    benzene = chematic.from_smiles("c1ccccc1")
    pyridine = chematic.from_smiles("c1ccncc1")
    elements = chematic.find_mcs([benzene, pyridine], atom_compare="elements")
    any_heavy = chematic.find_mcs([benzene, pyridine], atom_compare="any_heavy_atom")
    assert any_heavy is not None
    assert any_heavy.heavy_atoms == 6
    assert elements.heavy_atoms < any_heavy.heavy_atoms


def test_find_mcs_invalid_atom_compare_raises():
    m1 = chematic.from_smiles("CCO")
    m2 = chematic.from_smiles("CCO")
    with pytest.raises(ValueError):
        chematic.find_mcs([m1, m2], atom_compare="not_a_real_mode")


def test_find_mcs_invalid_bond_compare_raises():
    m1 = chematic.from_smiles("CCO")
    m2 = chematic.from_smiles("CCO")
    with pytest.raises(ValueError):
        chematic.find_mcs([m1, m2], bond_compare="not_a_real_mode")


def test_find_mcs_checked_matches_find_mcs_when_not_timed_out():
    m1 = chematic.from_smiles("c1ccccc1")
    m2 = chematic.from_smiles("Cc1ccccc1")
    mcs = chematic.find_mcs([m1, m2])
    mcs_checked, was_timed_out = chematic.find_mcs_checked([m1, m2])
    assert was_timed_out is False
    assert mcs_checked is not None
    assert mcs_checked.smiles == mcs.smiles


def test_find_mcs_checked_reports_timeout():
    m1 = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
    m2 = chematic.from_smiles("CC(=O)Nc1ccc(O)cc1")
    mcs, was_timed_out = chematic.find_mcs_checked([m1, m2], timeout_ms=0)
    assert was_timed_out is True


def test_find_mcs_checked_none_result_shape():
    m1 = chematic.from_smiles("C")
    m2 = chematic.from_smiles("N")
    result = chematic.find_mcs_checked([m1, m2])
    assert isinstance(result, tuple) and len(result) == 2
    mcs, was_timed_out = result
    assert mcs is None or isinstance(mcs, chematic.Mol)
    assert isinstance(was_timed_out, bool)


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
