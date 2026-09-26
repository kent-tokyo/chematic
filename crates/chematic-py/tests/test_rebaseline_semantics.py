"""Issues #634/#635: CIP E/Z bond identity and SMARTS aromaticity semantics."""

import chematic


def test_cip_ez_entries_name_their_double_bond():
    mol = chematic.from_smiles("COc1cc2nc(N3CCN(C(=O)/C(F)=C/c4ccccc4)CC3)nc(N)c2cc1OC")
    ez = [item for item in mol.cip_stereo(mode="accurate") if item["descriptor"] in ("E", "Z")]
    assert len(ez) == 1
    item = ez[0]
    a1, a2 = item["bond_atoms"]
    assert item["atom_idx"] == a1
    assert mol.bond_table[item["bond_idx"]][:2] == (a1, a2)
    assert sorted((a1, a2)) == [13, 15]  # the C=C, not the C-F bond 13-14
    assert item["descriptor"] == "Z"


def test_cip_accurate_ez_ranks_aryl_above_tert_butyl():
    mol = chematic.from_smiles(
        "Cc1oc(-c2ccc(C(F)(F)F)cc2)nc1COc1cccc(/C(=C/Cn2oc(=O)[nH]c2=O)C(C)(C)C)c1"
    )
    ez = {tuple(sorted(i["bond_atoms"])): i["descriptor"]
          for i in mol.cip_stereo(mode="accurate") if "bond_atoms" in i}
    assert ez == {(23, 24): "E"}


def test_smarts_sees_perceived_aromaticity_for_kekule_input():
    kekule = chematic.from_smiles("C1=CC=CC=C1O")
    aromatic = chematic.from_smiles("c1ccccc1O")
    for query in ["c", "c1ccccc1", "[#6]=[#6]", "C", "[N,O]", "cO", "[OH]"]:
        assert chematic.smarts_find(query, kekule) == chematic.smarts_find(query, aromatic), query
        assert kekule.find_matches(query) == aromatic.find_matches(query), query
        assert kekule.has_substructure(query) == aromatic.has_substructure(query), query
        assert chematic.smarts_match(query, kekule) == chematic.smarts_match(query, aromatic), query
    assert len(chematic.smarts_find("c", kekule)) == 6
    assert chematic.smarts_find("[#6]=[#6]", kekule) == []
    # A non-aromatic ring keeps aliphatic semantics.
    assert len(chematic.smarts_find("[#6]=[#6]", chematic.from_smiles("C1=CCCCC1"))) == 1
    assert chematic.bulk.substructure_match("c1ccccc1", [kekule, aromatic]) == [0, 1]


def test_ring_count_primitive_uses_sssr_not_symmetrized_rings():
    # Documented boundary (#635): [Rn] counts SSSR rings (cycle rank = 2 for
    # bicyclo[2.2.2]octane); RDKit's symmetrized ring set has 3, so RDKit
    # reports the bridgeheads as [R3].
    mol = chematic.from_smiles("C1CC2CCC1CC2")
    assert chematic.smarts_find("[R3]", mol) == []
    r2 = {m[0] for m in chematic.smarts_find("[R2]", mol)}
    assert {2, 5} <= r2 and len(r2) == 4  # bridgeheads + the shared SSSR bridge
