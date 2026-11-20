"""MCS query output retains the alternatives omitted by concrete Mol projection."""

import pytest
import chematic as ch


@pytest.mark.parametrize("left,right,options", [
    ("NO", "CO", {"atom_compare": "any"}),
    ("C=N", "C#N", {"bond_compare": "any"}),
    ("C1CC1", "N1CC1", {"atom_compare": "any"}),
])
def test_mcs_smarts_keeps_alternatives(left, right, options):
    smarts, timed_out = ch.find_mcs_smarts_checked(
        [ch.from_smiles(left), ch.from_smiles(right)], **options
    )
    assert timed_out is False
    assert smarts is not None and "," in smarts


def test_mcs_smarts_reports_resource_limit():
    with pytest.raises(ch.ChematicInputError) as error:
        ch.find_mcs_smarts_checked(
            [ch.from_smiles("CCC"), ch.from_smiles("CCC")], max_search_nodes=0
        )
    assert error.value.code == "mcs_resource_limit"



@pytest.mark.parametrize("inputs", [
    ["[13CH3][NH3+]"],
    ["[13CH3][NH3+]", "[NH3+][13CH3]"],
])
def test_mcs_smarts_preserves_charge_and_isotope_predicates(inputs):
    molecules = [ch.from_smiles(smiles) for smiles in inputs]
    smarts, timed_out = ch.find_mcs_smarts_checked(
        molecules, match_charge=True, match_isotope=True
    )
    query = ch.compile_smarts(smarts)
    assert not timed_out and query.atom_count == 2
    assert all(query.matches(mol) for mol in molecules)
    for other in ("[12CH3][NH3+]", "[13CH3]N", "[13CH3][15NH3+]"):
        assert not query.matches(ch.from_smiles(other)), (smarts, other)


def test_mcs_smarts_applies_charge_to_every_element_alternative():
    molecules = [ch.from_smiles(smiles) for smiles in ("[NH4+]", "[PH4+]")]
    smarts, timed_out = ch.find_mcs_smarts_checked(
        molecules, atom_compare="any", match_charge=True
    )
    query = ch.compile_smarts(smarts)
    assert not timed_out and all(query.matches(mol) for mol in molecules)
    assert not query.matches(ch.from_smiles("N"))
    assert not query.matches(ch.from_smiles("P"))



def test_mcs_bond_first_can_choose_the_smaller_dense_component():
    molecule = ch.from_smiles("CCCCC.C12C3C1C23")
    for bond_count_first, atoms in ((False, 5), (True, 4)):
        smarts, timed_out = ch.find_mcs_smarts_checked(
            [molecule, molecule], bond_count_first=bond_count_first
        )
        query = ch.compile_smarts(smarts)
        assert not timed_out and query.atom_count == atoms
        assert query.matches(molecule)


def test_mcs_round_trips_dummy_atoms_in_element_comparator():
    dummy = ch.from_smiles("[*]O")
    carbon = ch.from_smiles("CO")
    smarts, timed_out = ch.find_mcs_smarts_checked([dummy, dummy])
    assert not timed_out and "[#0]" in smarts
    query = ch.compile_smarts(smarts)
    assert query.matches(dummy) and not query.matches(carbon)
    projected = ch.find_mcs([dummy, dummy])
    assert query.matches(projected)
    mixed, timed_out = ch.find_mcs_smarts_checked([dummy, carbon])
    assert not timed_out and ch.compile_smarts(mixed).atom_count == 1


@pytest.mark.parametrize("left,right,inverted", [
    ("[C@H](C)(O)Cl", "C[C@@H](O)Cl", "C[C@H](O)Cl"),
    ("[C@H]1(N)CCOC1", "N[C@H]1CCOC1", "N[C@@H]1CCOC1"),
    ("F[C@](Cl)(Br)I", "I[C@@](Cl)(Br)F", "F[C@@](Cl)(Br)I"),
])
def test_mcs_smarts_round_trips_a_complete_tetrahedral_frame(left, right, inverted):
    molecules = [ch.from_smiles(s) for s in (left, right)]
    smarts, timed_out = ch.find_mcs_smarts_checked(
        molecules, match_chiral_tag=True, bond_count_first=True
    )
    query = ch.compile_smarts(smarts)
    assert not timed_out and "@" in smarts
    assert all(query.matches(mol, use_chirality=True) for mol in molecules)
    assert not query.matches(ch.from_smiles(inverted), use_chirality=True)


def test_mcs_smarts_missing_fourth_ligand_rejects_inverted_orientation():
    left, right, inverted = [
        ch.from_smiles(s) for s in ("F[C@](Cl)(Br)I", "F[C@H](Cl)Br", "F[C@@H](Cl)Br")
    ]
    smarts, timed_out = ch.find_mcs_smarts_checked([left, right], match_chiral_tag=True)
    query = ch.compile_smarts(smarts)
    assert not timed_out and query.atom_count == 4
    assert query.matches(left, use_chirality=True) and query.matches(right, use_chirality=True)
    assert not query.matches(inverted, use_chirality=True)


def test_mcs_smarts_two_ligands_require_chirality_presence():
    molecules = [ch.from_smiles(s) for s in ("C[C@H](O)Cl", "C[C@H](I)Cl", "C[C@@H](I)Cl")]
    smarts, timed_out = ch.find_mcs_smarts_checked(molecules, match_chiral_tag=True)
    query = ch.compile_smarts(smarts)
    assert not timed_out and query.atom_count == 3
    assert all(query.matches(mol, use_chirality=True) for mol in molecules)
    assert not query.matches(ch.from_smiles("CC(O)Cl"), use_chirality=True)


def test_mcs_smarts_reports_unsupported_square_planar_stereo():
    molecule = ch.from_smiles("[Pt@SP1](F)(Cl)(Br)I")
    with pytest.raises(ValueError, match="representable absolute tetrahedral stereo"):
        ch.find_mcs_smarts_checked([molecule, molecule], match_chiral_tag=True)
