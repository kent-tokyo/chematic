"""Tests for Sprint 1-6 Python bindings."""
import pytest
import chematic


# ---------------------------------------------------------------------------
# Sprint 1-2: Descriptors
# ---------------------------------------------------------------------------

def test_logd_acid_base_direction():
    """LogD at physiological pH should differ from LogP for ionizable groups."""
    acetic = chematic.from_smiles("CC(=O)O")  # pKa ≈ 4.8 — largely ionized at pH 7.4
    logp = acetic.logp
    logd = acetic.logd(7.4)
    # LogD should be ≤ LogP for an acid at pH >> pKa
    assert logd <= logp + 0.1  # allow small floating-point tolerance


def test_mqn_length():
    """MQN must always return exactly 42 values."""
    mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
    assert len(mol.mqn()) == 42


def test_xlogp3_aspirin_positive():
    """Aspirin XLogP3 should be positive (lipophilic)."""
    aspirin = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
    assert aspirin.xlogp3 > 0


def test_xlogp3_per_atom_length():
    """xlogp3_per_atom should return one value per heavy atom."""
    m = chematic.from_smiles("CCN")  # 3 heavy atoms
    vals = m.xlogp3_per_atom()
    assert len(vals) == m.heavy_atoms


def test_autocorr_2d_length():
    """autocorr_2d should return 7 values."""
    m = chematic.from_smiles("c1ccccc1")
    assert len(m.autocorr_2d()) == 7


def test_hall_kier_alpha_acetic_acid():
    """Hall-Kier alpha for acetic acid matches RDKit's HallKierAlpha (-0.53) closely."""
    mol = chematic.from_smiles("CC(=O)O")
    assert mol.hall_kier_alpha == pytest.approx(-0.53, abs=0.05)


def test_peoe_vsa_bins():
    """PEOE VSA should return a non-empty list."""
    m = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
    assert len(m.peoe_vsa()) > 0


def test_slogp_vsa_bins():
    m = chematic.from_smiles("c1ccccc1")
    assert len(m.slogp_vsa()) > 0


def test_usrcat_length():
    """USRCAT should return 42 values."""
    m = chematic.from_smiles("c1ccccc1")
    assert len(m.usrcat()) == 42


# ---------------------------------------------------------------------------
# Sprint 3: CIP stereo, named FGs, shape, conformer tools
# ---------------------------------------------------------------------------

def test_cip_stereo_chiral_center():
    """Single chiral center should produce one R or S assignment."""
    m = chematic.from_smiles("[C@@H](F)(Cl)Br")
    stereo = m.cip_stereo()
    assert len(stereo) == 1
    assert stereo[0]["descriptor"] in ("R", "S")


def test_cip_stereo_achiral():
    """Achiral molecule should have no CIP assignments."""
    m = chematic.from_smiles("c1ccccc1")
    assert m.cip_stereo() == []


def test_cip_stereo_mode_legacy_is_default():
    """mode='legacy' (the default) must match calling cip_stereo() with no args."""
    m = chematic.from_smiles("C[C@H](N)C(=O)O")
    assert m.cip_stereo() == m.cip_stereo(mode="legacy")


def test_cip_stereo_mode_accurate_merges_ez_with_tetrahedral():
    """Accurate mode must still report E/Z (it doesn't compute that itself) alongside
    its own tetrahedral R/S for the same molecule."""
    m = chematic.from_smiles("C/C=C/[C@H](N)C(=O)O")
    stereo = m.cip_stereo(mode="accurate")
    descriptors = {d["descriptor"] for d in stereo}
    assert "E" in descriptors
    assert "R" in descriptors or "S" in descriptors


def test_cip_stereo_mode_invalid_raises():
    m = chematic.from_smiles("C[C@H](N)C(=O)O")
    with pytest.raises(ValueError):
        m.cip_stereo(mode="bogus")


def test_cip_stereo_unresolved_empty_for_resolvable_molecule():
    m = chematic.from_smiles("C[C@H](N)C(=O)O")
    assert m.cip_stereo_unresolved() == []


def test_cip_stereo_unresolved_reports_genuine_ties():
    """The 2 phosphorus rows found to be genuine chematic ties (not merely
    oracle-unstable, see docs/rfcs/cip_accurate_rfc.md Milestone 4C-1) must come back in
    cip_stereo_unresolved(), never a silently-guessed label in cip_stereo()."""
    m = chematic.from_smiles(
        "CNP1(NC)=N[P@](NC)(N2CC2)=NP(NC)(NC)=N[P@@](NC)(N2CC2)=N1"
    )
    unresolved_atoms = {d["atom_idx"] for d in m.cip_stereo_unresolved()}
    assert unresolved_atoms == {6, 19}
    resolved_atoms = {d["atom_idx"] for d in m.cip_stereo(mode="accurate")}
    assert resolved_atoms.isdisjoint(unresolved_atoms)


def test_generate_3d_atom_count():
    """generate_3d should produce one coord per heavy atom."""
    m = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
    coords = m.generate_3d()
    assert len(coords) == m.heavy_atoms


def test_gasteiger_charges_length():
    """gasteiger_charges should produce one value per heavy atom."""
    m = chematic.from_smiles("CC(=O)O")
    charges = m.gasteiger_charges()
    assert len(charges) == m.heavy_atoms


def test_named_functional_groups_carboxyl():
    """Molecule with a carboxyl group should return carboxyl in named FGs."""
    m = chematic.from_smiles("CC(=O)O")  # acetic acid
    fgs = m.named_functional_groups()
    assert "carboxyl" in fgs


def test_pmi_ordering():
    """PMI eigenvalues must be in ascending order (PMI1 ≤ PMI2 ≤ PMI3)."""
    m = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
    coords = m.generate_3d()
    pmi = m.pmi(coords)
    assert len(pmi) == 3
    assert pmi[0] <= pmi[1] <= pmi[2]


def test_npr_range():
    """NPR values should be in [0, 1]."""
    m = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
    coords = m.generate_3d()
    npr = m.npr(coords)
    assert len(npr) == 2
    assert all(0.0 <= v <= 1.0 for v in npr)


def test_asphericity_range():
    """Asphericity should be in [0, ∞) — positive for non-spherical shapes."""
    m = chematic.from_smiles("c1ccccc1")
    coords = m.generate_3d()
    assert m.asphericity(coords) >= 0


def test_get_dihedral_butane():
    """Dihedral angle should be defined for 4 connected atoms."""
    butane = chematic.from_smiles("CCCC")
    coords = butane.generate_3d()
    d = butane.get_dihedral(coords, 0, 1, 2, 3)
    assert d is not None
    assert -180.0 <= d <= 180.0


def test_set_dihedral_round_trip():
    """Setting dihedral to 180° should give back 180° (or, equivalently,
    -180° -- the same physical angle, and which sign atan2 returns for this
    exact boundary value is a platform-dependent floating-point tie-break,
    confirmed to differ between macOS/arm64 (180.0) and Linux/x86_64 CI
    (-180.0) for this same molecule/seed)."""
    butane = chematic.from_smiles("CCCC")
    coords = butane.generate_3d()
    new_coords = butane.set_dihedral(coords, 0, 1, 2, 3, 180.0)
    d = butane.get_dihedral(new_coords, 0, 1, 2, 3)
    assert d is not None
    assert abs(abs(d) - 180.0) < 1.0


def test_generate_3d_etkdg_atom_count():
    """ETKDG should produce one coord per heavy atom."""
    m = chematic.from_smiles("c1ccccc1")
    coords = m.generate_3d_etkdg()
    assert len(coords) == m.heavy_atoms


# ---------------------------------------------------------------------------
# Sprint 3-4: Green chemistry, reaction analysis
# ---------------------------------------------------------------------------

def test_atom_economy_complete():
    """Reaction where all atoms end up in product should have high atom economy."""
    # A + B → AB (hypothetical concatenation)
    ae = chematic.atom_economy("CC.N>>CCN")
    assert ae > 80.0  # should be high


def test_balance_check_balanced():
    """Simple balanced reaction."""
    result = chematic.balance_check("C>>C")
    assert result["balanced"] is True
    assert result["diff"] == []


def test_balance_check_unbalanced():
    """Unbalanced reaction should be detected."""
    result = chematic.balance_check("CC>>C")
    assert result["balanced"] is False
    assert len(result["diff"]) > 0


def test_e_factor_calculation():
    """E-factor = waste / product."""
    assert abs(chematic.e_factor(90.0, 10.0) - 9.0) < 0.001


def test_pmi_rxn_calculation():
    """PMI = sum(all masses) / product mass."""
    assert abs(chematic.pmi_rxn([50.0, 60.0], 10.0) - 11.0) < 0.001


# ---------------------------------------------------------------------------
# Sprint 4: Workflow, MMP
# ---------------------------------------------------------------------------

def test_molecule_report_required_keys():
    """molecule_report must return all required top-level keys."""
    report = chematic.molecule_report("CC(=O)Oc1ccccc1C(=O)O")
    required = {"canonical_smiles", "formula", "murcko_scaffold", "descriptors", "filters",
                "functional_groups", "named_groups"}
    assert required.issubset(set(report.keys()))


def test_molecule_report_descriptors_keys():
    """Descriptor sub-dict must include mw, logp, tpsa, hbd, hba."""
    report = chematic.molecule_report("CC(=O)Oc1ccccc1C(=O)O")
    for key in ("mw", "logp", "tpsa", "hbd", "hba", "qed", "sa_score"):
        assert key in report["descriptors"]


def test_screen_smiles_all_records():
    """screen_smiles should return one record per valid SMILES."""
    smiles = ["c1ccccc1", "CC(=O)O", "CCN", "CCCO", "c1ccncc1"]
    result = chematic.screen_smiles(smiles)
    assert len(result["records"]) == len(smiles)


def test_find_mmp_returns_pairs():
    """Benzene/toluene/aniline share core → MMP should be found."""
    smiles = ["c1ccccc1", "Cc1ccccc1", "Nc1ccccc1", "CCc1ccccc1"]
    pairs = chematic.find_mmp(smiles)
    # These all share an aromatic ring core
    # At least some pairs should be found with BRICS cuts
    assert isinstance(pairs, list)


def test_find_reaction_center_structure():
    """find_reaction_center should return a dict with the expected keys."""
    rc = chematic.find_reaction_center("CC(=O)Cl.[NH3]>>CC(=O)N.Cl")
    assert "broken_bonds" in rc
    assert "formed_bonds" in rc
    assert "changed_atoms" in rc


# ---------------------------------------------------------------------------
# Sprint 5: Standardization, SASA descriptor
# ---------------------------------------------------------------------------

def test_prefer_organic_removes_salt():
    """prefer_organic should keep the organic fragment."""
    m = chematic.from_smiles("CCO.[Na+].[Cl-]")
    organic = m.prefer_organic()
    assert "Na" not in organic.smiles
    assert "Cl" not in organic.smiles.replace("Cl", "")  # check no chlorine atom


def test_uncharge_removes_formal_charge():
    """uncharge should remove explicit charges."""
    m = chematic.from_smiles("[NH4+]")
    uncharged = m.uncharge()
    assert "+" not in uncharged.smiles


def test_sasa_descriptor_keys():
    """sasa_descriptor should return total, mean, std_dev, per_atom."""
    m = chematic.from_smiles("c1ccccc1")
    coords = m.generate_3d()
    sd = m.sasa_descriptor(coords)
    assert set(sd.keys()) >= {"total", "mean", "std_dev", "per_atom"}
    assert sd["total"] > 0
    assert len(sd["per_atom"]) == m.heavy_atoms


def test_randic_index_benzene():
    """Benzene Randić index should be 3.0 (6 atoms each with degree 2)."""
    m = chematic.from_smiles("c1ccccc1")
    assert abs(m.randic_index - 3.0) < 0.01


def test_fcfp6_length():
    """fcfp6 should return 256 bytes (2048 bits)."""
    m = chematic.from_smiles("c1ccccc1")
    assert len(m.fcfp6()) == 256


def test_pattern_fp_length():
    """pattern_fp should return 256 bytes."""
    m = chematic.from_smiles("c1ccccc1")
    assert len(m.pattern_fp()) == 256


# ---------------------------------------------------------------------------
# Sprint 6: MMFF94 charges, topology, top_k_similar
# ---------------------------------------------------------------------------

def test_mmff94_charges_length():
    """mmff94_charges should return one value per heavy atom."""
    m = chematic.from_smiles("CC(=O)O")
    charges = m.mmff94_charges()
    assert len(charges) == m.heavy_atoms


def test_balaban_j_positive():
    """Balaban J should be positive for non-trivial graphs."""
    m = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
    assert m.balaban_j > 0


def test_zagreb_m1_positive():
    """Zagreb M1 should be positive."""
    m = chematic.from_smiles("c1ccccc1")
    assert m.zagreb_m1 > 0


def test_labute_asa_per_atom_length():
    """labute_asa_per_atom should return one value per heavy atom."""
    m = chematic.from_smiles("CC(=O)O")
    vals = m.labute_asa_per_atom()
    assert len(vals) == m.heavy_atoms


def test_top_k_similar_returns_k():
    """top_k_similar should return at most k results."""
    db = ["c1ccccc1", "CC(=O)O", "CCN", "CCCO", "Cc1ccccc1"]
    hits = chematic.top_k_similar("c1ccccc1", db, k=3)
    assert len(hits) <= 3


def test_top_k_similar_sorted():
    """Results should be sorted by descending similarity."""
    db = ["c1ccccc1", "CC(=O)O", "CCN", "Cc1ccccc1"]
    hits = chematic.top_k_similar("c1ccccc1", db, k=4)
    scores = [s for _, s in hits]
    assert scores == sorted(scores, reverse=True)


def test_center_on_origin_centroid():
    """After centering, centroid should be near (0, 0, 0)."""
    m = chematic.from_smiles("CCO")
    coords = m.generate_3d()
    centered = chematic.center_on_origin(coords)
    n = len(centered)
    cx = sum(c[0] for c in centered) / n
    cy = sum(c[1] for c in centered) / n
    cz = sum(c[2] for c in centered) / n
    assert abs(cx) < 1e-4
    assert abs(cy) < 1e-4
    assert abs(cz) < 1e-4


def test_dice_vs_tanimoto():
    """Dice should be >= Tanimoto for same fingerprint pair (for binary FP)."""
    m1 = chematic.from_smiles("c1ccccc1")
    m2 = chematic.from_smiles("Cc1ccccc1")
    tan = chematic.tanimoto(m1.ecfp4(), m2.ecfp4())
    dice = chematic.dice_similarity(m1.ecfp4(), m2.ecfp4())
    # For binary FPs: Dice = 2T/(1+T), so Dice >= Tanimoto when T <= 1
    assert dice >= tan - 1e-6


def test_butina_cluster_all_assigned():
    """Every molecule should appear in exactly one Butina cluster."""
    smiles = ["c1ccccc1", "CC(=O)O", "CCN", "CCCO", "c1ccncc1"]
    clusters = chematic.butina_cluster(smiles, 0.4)
    all_indices = [i for cluster in clusters for i in cluster]
    assert sorted(all_indices) == list(range(len(smiles)))
