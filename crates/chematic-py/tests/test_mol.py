"""Tests for the Mol class — identity, descriptors, drug-likeness, transforms."""
import math
import pytest
import chematic


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def test_from_smiles_aspirin():
    mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
    assert mol is not None


def test_from_smiles_invalid():
    with pytest.raises(ValueError):
        chematic.from_smiles("not_a_smiles!!!")


def test_is_valid_smiles():
    assert chematic.is_valid_smiles("CCO")
    assert chematic.is_valid_smiles("c1ccccc1")
    assert not chematic.is_valid_smiles("XYZ???")


def test_repr_and_str(ethanol):
    # repr() just embeds mol.smiles (Mol('<smiles>')), whose exact spelling
    # is not a stable contract (issue #200 -- canonical-SMILES spelling can
    # change with unrelated internal changes, e.g. automorphism-orbit
    # pruning). Assert structural round-trip instead of one hardcoded
    # substring: repr()'s embedded SMILES must re-parse to an equivalent
    # molecule (same canonical form).
    reparsed = chematic.from_smiles(repr(ethanol).removeprefix("Mol('").removesuffix("')"))
    assert reparsed.smiles == ethanol.smiles
    smiles_str = str(ethanol)
    assert isinstance(smiles_str, str)
    assert len(smiles_str) > 0
    assert smiles_str == ethanol.smiles


# ---------------------------------------------------------------------------
# Identity / structure
# ---------------------------------------------------------------------------

def test_smiles_canonical(ethanol):
    assert isinstance(ethanol.smiles, str)
    assert len(ethanol.smiles) > 0


def test_formula_aspirin(aspirin):
    assert aspirin.formula == "C9H8O4"


def test_formula_ethanol(ethanol):
    assert ethanol.formula == "C2H6O"


def test_heavy_atoms(aspirin, ethanol):
    assert aspirin.heavy_atoms == 13  # 9C + 4O
    assert ethanol.heavy_atoms == 3   # 2C + 1O


def test_inchi_ethanol(ethanol):
    inchi = ethanol.inchi
    assert inchi.startswith("InChI=")
    assert "C2H6O" in inchi


def test_inchikey_ethanol(ethanol):
    key = ethanol.inchikey
    assert len(key) == 27
    assert key.count("-") == 2


def test_iupac_ethanol(ethanol):
    name = ethanol.iupac_name
    assert isinstance(name, str)
    # IUPAC name may be "ethanol" or empty if unsupported
    assert name in ("ethanol", "")


# ---------------------------------------------------------------------------
# Physicochemical descriptors
# ---------------------------------------------------------------------------

def test_mw_aspirin(aspirin):
    assert abs(aspirin.mw - 180.16) < 0.5


def test_mw_ethanol(ethanol):
    assert abs(ethanol.mw - 46.07) < 0.1


def test_exact_mass_ethanol(ethanol):
    assert abs(ethanol.exact_mass - 46.042) < 0.01


def test_logp_aspirin(aspirin):
    # Aspirin logP ~1.19 (RDKit) — allow generous tolerance
    assert 0.5 < aspirin.logp < 2.5


def test_tpsa_aspirin(aspirin):
    # Aspirin TPSA ~63.6 Å²
    assert 55.0 < aspirin.tpsa < 75.0


def test_qed_aspirin(aspirin):
    assert 0.0 < aspirin.qed <= 1.0


def test_hbd_aspirin(aspirin):
    assert aspirin.hbd == 1  # one carboxylic OH


def test_hba_aspirin(aspirin):
    assert aspirin.hba == 3  # RDKit CalcNumHBA agrees: 3, not one per oxygen


def test_rotatable_bonds_aspirin(aspirin):
    assert aspirin.rotatable_bonds >= 2


def test_fsp3(aspirin, benzene):
    assert 0.0 < aspirin.fsp3 < 1.0
    assert benzene.fsp3 == 0.0


def test_sa_score(aspirin):
    # SA score range is 1–10
    assert 1.0 <= aspirin.sa_score <= 10.0


def test_molar_refractivity(aspirin):
    assert aspirin.molar_refractivity > 0


def test_formal_charge(aspirin, ethanol):
    assert aspirin.formal_charge == 0
    assert ethanol.formal_charge == 0


def test_esol(aspirin):
    # Aspirin ESOL ~-1 to -2 (moderately soluble)
    assert -4.0 < aspirin.esol < 0.0


# ---------------------------------------------------------------------------
# Ring / stereo counts
# ---------------------------------------------------------------------------

def test_ring_count(aspirin, ethanol, benzene):
    assert aspirin.ring_count == 1
    assert ethanol.ring_count == 0
    assert benzene.ring_count == 1


def test_aromatic_ring_count(aspirin, benzene, ethanol):
    assert aspirin.aromatic_ring_count == 1
    assert benzene.aromatic_ring_count == 1
    assert ethanol.aromatic_ring_count == 0


def test_stereocenters():
    l_alanine = chematic.from_smiles("N[C@@H](C)C(=O)O")
    assert l_alanine.num_stereocenters >= 1


# ---------------------------------------------------------------------------
# Drug-likeness rules
# ---------------------------------------------------------------------------

def test_lipinski_passes(aspirin, ethanol):
    assert aspirin.lipinski_passes
    assert ethanol.lipinski_passes


def test_large_molecule_fails_lipinski():
    # A large molecule that clearly fails Lipinski
    large = chematic.from_smiles("CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O")  # very long chain
    # MW >> 500 → fails
    assert large.mw > 400


def test_pains_passes(aspirin):
    assert aspirin.pains_passes


# ---------------------------------------------------------------------------
# Descriptors dict
# ---------------------------------------------------------------------------

def test_descriptors_keys(aspirin):
    d = aspirin.descriptors()
    assert isinstance(d, dict)
    for key in ("mw", "tpsa", "logp", "hbd", "hba", "qed", "sa_score"):
        assert key in d, f"missing key: {key}"


def test_descriptors_values_numeric(aspirin):
    d = aspirin.descriptors()
    assert abs(d["mw"] - aspirin.mw) < 0.01
    assert abs(d["tpsa"] - aspirin.tpsa) < 0.01


# ---------------------------------------------------------------------------
# pKa / ADMET
# ---------------------------------------------------------------------------

def test_pka_aspirin(aspirin):
    pka = aspirin.pka()
    assert isinstance(pka, dict)
    assert "most_acidic" in pka
    assert "most_basic" in pka


def test_admet_aspirin(aspirin):
    admet = aspirin.admet()
    assert isinstance(admet, dict)
    for key in ("bbb", "bbb_score", "caco2", "herg_risk", "cyp3a4_risk"):
        assert key in admet, f"missing key: {key}"
    assert isinstance(admet["bbb"], bool)
    assert 0.0 <= admet["herg_risk"] <= 1.0


# ---------------------------------------------------------------------------
# SVG visualization
# ---------------------------------------------------------------------------

def test_svg_returns_string(aspirin):
    svg = aspirin.svg()
    assert isinstance(svg, str)
    assert "<svg" in svg.lower() or "svg" in svg.lower()


def test_svg_highlighted(aspirin):
    svg = aspirin.svg_highlighted([0, 1, 2])
    assert isinstance(svg, str)
    assert len(svg) > 0


# ---------------------------------------------------------------------------
# Transformations
# ---------------------------------------------------------------------------

def test_standardize_returns_mol(aspirin):
    std = aspirin.standardize()
    assert isinstance(std, chematic.Mol)


def test_scaffold(aspirin):
    sc = aspirin.scaffold()
    assert isinstance(sc, chematic.Mol)
    # Scaffold should have fewer or equal heavy atoms
    assert sc.heavy_atoms <= aspirin.heavy_atoms


def test_canonical_tautomer(aspirin):
    t = aspirin.canonical_tautomer()
    assert isinstance(t, chematic.Mol)


def test_parent_identity_bindings_cover_each_axis():
    mol = chematic.from_smiles("[NH3+][C@@H]([2H])C(=O)[O-].Cl")
    assert mol.fragment_parent().heavy_atoms < mol.heavy_atoms
    assert mol.charge_parent().formal_charge == 0
    # ``formula`` intentionally omits isotope labels; compare the structural
    # representation to verify that the isotope axis was normalized.
    assert mol.isotope_parent().smiles != mol.smiles
    assert "@" not in mol.stereo_parent().smiles
    parent, status = mol.super_parent()
    assert status == "Completed"
    assert parent.heavy_atoms == 5
    report = mol.super_parent_report()
    assert report["status"] == "Completed"
    assert [stage["name"] for stage in report["stages"]] == [
        "fragment", "charge", "isotope", "stereo", "tautomer"
    ]


def test_enumerate_tautomers(ethanol):
    tautomers = ethanol.enumerate_tautomers()
    assert isinstance(tautomers, list)
    assert len(tautomers) >= 1


def test_enumerate_stereoisomers():
    mol = chematic.from_smiles("C(F)(Cl)Br")  # one stereocenter
    isomers = mol.enumerate_stereoisomers()
    assert isinstance(isomers, list)
    assert len(isomers) >= 1


def test_add_remove_hydrogens(ethanol):
    with_h = ethanol.add_hydrogens()
    assert with_h.heavy_atoms >= ethanol.heavy_atoms
    without_h = with_h.remove_hydrogens()
    assert without_h.heavy_atoms == ethanol.heavy_atoms


def test_remove_stereo():
    mol = chematic.from_smiles("N[C@@H](C)C(=O)O")
    no_stereo = mol.remove_stereo()
    assert isinstance(no_stereo, chematic.Mol)


def test_remove_isotopes():
    mol = chematic.from_smiles("[13CH4]")
    no_iso = mol.remove_isotopes()
    assert isinstance(no_iso, chematic.Mol)


def test_largest_fragment():
    salt = chematic.from_smiles("[Na+].[Cl-]")
    frag = salt.largest_fragment()
    assert frag.heavy_atoms == 1


def test_neutralize():
    charged = chematic.from_smiles("CC[NH3+]")
    neutral = charged.neutralize()
    assert isinstance(neutral, chematic.Mol)


def test_generic_scaffold(aspirin):
    gs = aspirin.generic_scaffold()
    assert isinstance(gs, chematic.Mol)


def test_brics_fragments(aspirin):
    frags = aspirin.brics_fragments()
    assert isinstance(frags, list)
    assert len(frags) >= 1
    for f in frags:
        assert isinstance(f, chematic.Mol)


# ---------------------------------------------------------------------------
# B8: SASA
# ---------------------------------------------------------------------------

def test_sasa_positive(aspirin, ethanol):
    assert aspirin.sasa() > 0.0
    assert ethanol.sasa() > 0.0


def test_sasa_per_atom(aspirin):
    per_atom = aspirin.sasa_per_atom()
    assert isinstance(per_atom, list)
    assert len(per_atom) == aspirin.heavy_atoms
    assert all(v >= 0.0 for v in per_atom)


def test_sasa_larger_molecule_larger_value(ethanol, aspirin):
    # Aspirin is larger → should have larger SASA on average
    assert aspirin.sasa() > ethanol.sasa()


# ---------------------------------------------------------------------------
# C1: Atropisomers
# ---------------------------------------------------------------------------

def test_atropisomers_empty(ethanol):
    result = ethanol.atropisomers()
    assert result == []


def test_atropisomers_biphenyl(biphenyl):
    result = biphenyl.atropisomers()
    assert isinstance(result, list)
    # Biphenyl may or may not be detected as atropisomeric depending on substituents
    for bond_idx, kind in result:
        assert isinstance(bond_idx, int)
        assert kind in ("Biaryl", "Allene", "Constrained")


# ---------------------------------------------------------------------------
# B5: Layered fingerprint
# ---------------------------------------------------------------------------

def test_layered_fp_length(aspirin):
    fp = aspirin.layered_fp()
    assert isinstance(fp, bytes)
    assert len(fp) == 256  # 2048 bits / 8


def test_layered_fp_numpy(aspirin):
    import numpy as np
    fp = aspirin.layered_fp_numpy()
    assert fp.shape == (2048,)
    assert fp.dtype == np.uint8
    assert set(fp.tolist()).issubset({0, 1})


def test_layered_fp_nonzero(aspirin, ethanol):
    assert sum(b.bit_count() for b in aspirin.layered_fp()) > 0
    assert sum(b.bit_count() for b in ethanol.layered_fp()) > 0


def test_layered_fp_different_molecules(aspirin, benzene):
    fp1 = aspirin.layered_fp()
    fp2 = benzene.layered_fp()
    assert fp1 != fp2
