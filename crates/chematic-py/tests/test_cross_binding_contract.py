"""Consume the shared Rust/Python/Node/WASM contract fixture."""

import json
from pathlib import Path

import pytest

import chematic


_FIXTURE_PATH = Path(__file__).parents[3] / "validation" / "cross_binding_contract.json"
_DOCUMENT = json.loads(_FIXTURE_PATH.read_text())


def test_shared_fixture_schema_is_stable():
    assert _DOCUMENT["schema_version"] == 1
    assert len(_DOCUMENT["fixtures"]) == 4
    operations = _DOCUMENT["operation_manifest"]["operations"]
    assert len(operations) == 56
    assert len({operation["id"] for operation in operations}) == len(operations)
    assert all(len(operation["bindings"]) == 4 for operation in operations)
    assert all(operation["test_anchors"] for operation in operations)
    assert _DOCUMENT["descriptor_contract"]["schema_version"] == 1
    assert _DOCUMENT["standardization_contract"]["schema_version"] == 1
    assert _DOCUMENT["descriptor_contract"]["fields"]["tpsa"]["unit"] == "A2"
    assert _DOCUMENT["descriptor_contract"]["fields"]["exact_mass"]["unit"] == "Da"
    assert _DOCUMENT["descriptor_contract"]["fields"]["formula"]["unit"] == "Hill notation"
    assert _DOCUMENT["fingerprint_contract"]["schema_version"] == 1
    assert _DOCUMENT["fingerprint_contract"]["operations"]["ecfp4"]["bytes"] == 256
    assert _DOCUMENT["fingerprint_contract"]["operations"]["maccs"]["bytes"] == 21
    assert _DOCUMENT["fingerprint_contract"]["operations"]["rdkit_rdk"]["bytes"] == 256
    assert _DOCUMENT["fingerprint_contract"]["operations"]["rdkit_path"]["bytes"] == 256
    assert _DOCUMENT["fingerprint_detail_contract"]["schema_version"] == 1
    assert _DOCUMENT["fingerprint_detail_contract"]["operations"]["rdkit_ecfp4_detail"]["configuration"]["radius"] == 2
    assert _DOCUMENT["batch_canonicalization_contract"]["schema_version"] == 1
    assert _DOCUMENT["smiles_validity_contract"]["schema_version"] == 1
    assert _DOCUMENT["inchi_contract"]["schema_version"] == 1
    assert _DOCUMENT["fragment_parent_contract"]["schema_version"] == 1
    assert _DOCUMENT["charge_parent_contract"]["schema_version"] == 1
    assert _DOCUMENT["isotope_parent_contract"]["schema_version"] == 1
    assert _DOCUMENT["stereo_parent_contract"]["schema_version"] == 1
    assert _DOCUMENT["canonical_tautomer_contract"]["schema_version"] == 1
    assert _DOCUMENT["neutralize_charges_contract"]["schema_version"] == 1
    assert _DOCUMENT["largest_fragment_contract"]["schema_version"] == 1
    assert _DOCUMENT["tautomer_parent_contract"]["schema_version"] == 1
    assert _DOCUMENT["super_parent_contract"]["schema_version"] == 1
    assert _DOCUMENT["hydrogen_contract"]["schema_version"] == 1
    assert _DOCUMENT["super_parent_report_contract"]["schema_version"] == 1
    assert _DOCUMENT["smarts_validity_contract"]["schema_version"] == 1
    assert _DOCUMENT["mol_v2000_contract"]["schema_version"] == 1
    assert _DOCUMENT["mol_v3000_contract"]["schema_version"] == 1
    assert _DOCUMENT["xyz_contract"]["schema_version"] == 1
    assert _DOCUMENT["extxyz_contract"]["schema_version"] == 1
    assert _DOCUMENT["pdb_contract"]["schema_version"] == 1
    assert _DOCUMENT["qcschema_contract"]["schema_version"] == 1
    assert _DOCUMENT["qcschema_atomic_input_contract"]["schema_version"] == 1
    assert _DOCUMENT["qcschema_atomic_result_contract"]["schema_version"] == 1
    assert _DOCUMENT["orca_input_contract"]["schema_version"] == 1
    assert _DOCUMENT["mol2_contract"]["schema_version"] == 1
    assert _DOCUMENT["cml_contract"]["schema_version"] == 1
    assert _DOCUMENT["cdxml_contract"]["schema_version"] == 1
    assert _DOCUMENT["mmcif_contract"]["schema_version"] == 1
    assert _DOCUMENT["moljson_contract"]["schema_version"] == 1
    assert _DOCUMENT["cjson_contract"]["schema_version"] == 1
    assert _DOCUMENT["xyz_batch_contract"]["schema_version"] == 1
    assert _DOCUMENT["rxn_document_contract"]["schema_version"] == 1
    assert _DOCUMENT["reaction_smarts_contract"]["schema_version"] == 1
    assert _DOCUMENT["semantic_expansion_contract"]["schema_version"] == 1


def test_python_binding_matches_shared_extxyz_contract():
    contract = _DOCUMENT["extxyz_contract"]
    actual = chematic.from_extxyz(contract["input"])
    expected = contract["expected"]
    assert all(
        observed == pytest.approx(reference)
        for observed, reference in zip(actual["coords"], expected["coords"])
    )
    assert actual["lattice"] == expected["lattice"]
    assert actual["properties"] == expected["properties"]
    assert actual["info"] == expected["info"]


def test_python_extxyz_writer_matches_shared_contract():
    contract = _DOCUMENT["extxyz_contract"]["writer"]
    mol = chematic.from_smiles(contract["smiles"])
    actual_text = chematic.to_extxyz(
        mol,
        contract["coords"],
        contract["lattice"],
        contract["properties"],
        contract["info"],
    )
    actual = chematic.from_extxyz(actual_text)
    expected = contract["expected"]
    assert all(
        observed == pytest.approx(reference)
        for observed, reference in zip(actual["coords"], expected["coords"])
    )
    assert actual["lattice"] == expected["lattice"]
    assert actual["properties"] == expected["properties"]
    assert actual["info"] == expected["info"]


@pytest.mark.parametrize("smiles", _DOCUMENT["smiles_validity_contract"]["accepted"])
def test_python_smiles_validity_contract_accepts(smiles):
    assert chematic.is_valid_smiles(smiles)


@pytest.mark.parametrize("smiles", _DOCUMENT["smiles_validity_contract"]["rejected"])
def test_python_smiles_validity_contract_rejects(smiles):
    assert not chematic.is_valid_smiles(smiles)


@pytest.mark.parametrize("case", _DOCUMENT["inchi_contract"]["cases"])
def test_python_inchi_contract(case):
    mol = chematic.from_smiles(case["smiles"])
    assert mol.inchi == case["inchi"]
    assert mol.inchikey == case["inchikey"]


@pytest.mark.parametrize("case", _DOCUMENT["fragment_parent_contract"]["cases"])
def test_python_fragment_parent_contract(case):
    parent = chematic.from_smiles(case["smiles"]).fragment_parent()
    assert parent.smiles == case["canonical_smiles"]


@pytest.mark.parametrize(
    "contract_name, method_name",
    [
        ("charge_parent_contract", "charge_parent"),
        ("isotope_parent_contract", "isotope_parent"),
        ("stereo_parent_contract", "stereo_parent"),
    ],
)
def test_python_parent_transform_contracts(contract_name, method_name):
    for case in _DOCUMENT[contract_name]["cases"]:
        parent = getattr(chematic.from_smiles(case["smiles"]), method_name)()
        assert parent.smiles == case["canonical_smiles"]


def test_python_canonical_tautomer_contract():
    for case in _DOCUMENT["canonical_tautomer_contract"]["cases"]:
        tautomer = chematic.from_smiles(case["smiles"]).canonical_tautomer()
        assert tautomer.smiles == case["canonical_smiles"]


def test_python_neutralize_charges_contract():
    for case in _DOCUMENT["neutralize_charges_contract"]["cases"]:
        neutral = chematic.from_smiles(case["smiles"]).neutralize()
        assert neutral.smiles == case["canonical_smiles"]


def test_python_largest_fragment_contract():
    for case in _DOCUMENT["largest_fragment_contract"]["cases"]:
        fragment = chematic.from_smiles(case["smiles"]).largest_fragment()
        assert fragment.smiles == case["canonical_smiles"]


@pytest.mark.parametrize("contract_name, method_name", [
    ("tautomer_parent_contract", "tautomer_parent"),
    ("super_parent_contract", "super_parent"),
])
def test_python_budgeted_parent_contracts(contract_name, method_name):
    limits = _DOCUMENT[contract_name]["limits"]
    for case in _DOCUMENT[contract_name]["cases"]:
        parent, status = getattr(chematic.from_smiles(case["smiles"]), method_name)(
            limits["max_transforms"], limits["max_tautomers"]
        )
        assert parent.smiles == case["canonical_smiles"]
        assert status == case["status"]


def test_python_hydrogen_contract():
    for case in _DOCUMENT["hydrogen_contract"]["cases"]:
        molecule = chematic.from_smiles(case["smiles"])
        assert molecule.add_hydrogens().smiles == case["added_canonical_smiles"]
        assert molecule.add_hydrogens().remove_hydrogens().smiles == case["removed_canonical_smiles"]


def test_python_super_parent_report_contract():
    contract = _DOCUMENT["super_parent_report_contract"]
    limits = contract["limits"]
    for case in contract["cases"]:
        report = chematic.from_smiles(case["smiles"]).super_parent_report(
            limits["max_transforms"], limits["max_tautomers"]
        )
        assert report["smiles"] == case["canonical_smiles"]
        assert report["status"] == case["status"]
        assert [stage["name"] for stage in report["stages"]] == contract["stage_names"]
        assert all(stage["smiles"] for stage in report["stages"])


@pytest.mark.parametrize("smarts", _DOCUMENT["smarts_validity_contract"]["accepted"])
def test_python_smarts_validity_source_only_contract_accepts(smarts):
    assert chematic.is_valid_smarts(smarts)


@pytest.mark.parametrize("smarts", _DOCUMENT["smarts_validity_contract"]["rejected"])
def test_python_smarts_validity_source_only_contract_rejects(smarts):
    assert not chematic.is_valid_smarts(smarts)


def test_python_binding_matches_shared_mol_v2000_contract():
    contract = _DOCUMENT["mol_v2000_contract"]
    mol = chematic.from_mol_block(contract["input"])
    assert mol.heavy_atoms == contract["expected"]["atom_count"]
    assert mol.smiles == contract["expected"]["canonical_smiles"]


def test_python_binding_matches_shared_mol_v2000_roundtrip_contract():
    contract = _DOCUMENT["mol_v2000_contract"]
    mol = chematic.from_mol_block(contract["input"])
    roundtripped = chematic.from_mol_block(mol.to_mol_block())
    assert roundtripped.heavy_atoms == contract["expected"]["atom_count"]
    assert roundtripped.smiles == contract["expected"]["canonical_smiles"]


def test_python_binding_matches_shared_mol_v3000_contract():
    contract = _DOCUMENT["mol_v3000_contract"]
    mol = chematic.from_mol_v3000(contract["input"])
    assert mol.heavy_atoms == contract["expected"]["atom_count"]
    assert mol.smiles == contract["expected"]["canonical_smiles"]


def test_python_binding_matches_shared_mol_v3000_roundtrip_contract():
    contract = _DOCUMENT["mol_v3000_contract"]
    mol = chematic.from_mol_v3000(contract["input"])
    roundtripped = chematic.from_mol_v3000(
        mol.to_mol_v3000([[0.0, 0.0], [1.5, 0.0]], name="shared-contract-roundtrip")
    )
    assert roundtripped.heavy_atoms == contract["expected"]["atom_count"]
    assert roundtripped.smiles == contract["expected"]["canonical_smiles"]


def test_python_binding_matches_shared_xyz_contract():
    contract = _DOCUMENT["xyz_contract"]
    mol, coords = chematic.from_xyz(contract["input"])
    expected = contract["expected"]
    assert mol.heavy_atoms == expected["heavy_atoms"]
    assert mol.smiles == expected["canonical_smiles"]
    assert all(
        observed == pytest.approx(reference)
        for observed, reference in zip(coords, expected["coords"])
    )


def test_python_binding_matches_shared_xyz_roundtrip_contract():
    contract = _DOCUMENT["xyz_contract"]
    mol, coords = chematic.from_xyz(contract["input"])
    serialized = mol.to_xyz(coords, "shared-contract-roundtrip")
    roundtripped, roundtripped_coords = chematic.from_xyz(serialized)
    assert roundtripped.heavy_atoms == contract["expected"]["heavy_atoms"]
    assert roundtripped.smiles == contract["expected"]["canonical_smiles"]
    assert all(
        observed == pytest.approx(reference)
        for observed, reference in zip(roundtripped_coords, coords)
    )


def test_python_binding_matches_shared_pdb_contract():
    contract = _DOCUMENT["pdb_contract"]
    mol, coords = chematic.from_pdb(contract["input"])
    assert mol.heavy_atoms == contract["expected"]["atom_count"]
    assert all(
        observed == pytest.approx(reference)
        for observed, reference in zip(coords, contract["expected"]["coords"])
    )


def test_python_binding_matches_shared_mol2_contract():
    contract = _DOCUMENT["mol2_contract"]
    mol = chematic.from_mol2(contract["input"])
    assert mol.heavy_atoms == contract["expected"]["atom_count"]


def test_python_binding_matches_shared_mol2_roundtrip_contract():
    contract = _DOCUMENT["mol2_contract"]
    mol = chematic.from_mol2(contract["input"])
    roundtripped = chematic.from_mol2(mol.to_mol2())
    assert roundtripped.heavy_atoms == contract["expected"]["atom_count"]
    assert len(roundtripped.bond_table) == contract["expected"]["bond_count"]


def test_python_binding_matches_shared_cml_contract():
    contract = _DOCUMENT["cml_contract"]
    mol = chematic.from_cml(contract["input"])
    assert mol.heavy_atoms == contract["expected"]["atom_count"]


def test_python_binding_matches_shared_cml_roundtrip_contract():
    contract = _DOCUMENT["cml_contract"]
    mol = chematic.from_cml(contract["input"])
    roundtripped = chematic.from_cml(mol.to_cml())
    assert roundtripped.heavy_atoms == contract["expected"]["atom_count"]


def test_python_cjson_topology_contract():
    contract = _DOCUMENT["cjson_contract"]
    mol, coords = chematic.from_cjson(contract["input"])
    assert mol.heavy_atoms == contract["expected"]["atom_count"]
    assert len(mol.bond_table) == contract["expected"]["bond_count"]
    assert coords == contract["expected"]["coords"]
    with pytest.raises(ValueError):
        chematic.from_cjson(contract["malformed_input"])


def test_python_cjson_roundtrip_contract():
    contract = _DOCUMENT["cjson_contract"]
    mol, coords = chematic.from_cjson(contract["input"])
    reparsed, reparsed_coords = chematic.from_cjson(mol.to_cjson(coords))
    assert reparsed.heavy_atoms == contract["expected"]["atom_count"]
    assert len(reparsed.bond_table) == contract["expected"]["bond_count"]
    assert reparsed_coords == coords


def test_python_binding_matches_source_only_pdb_strict_contract():
    contract = _DOCUMENT["pdb_strict_contract"]
    mol, coords = chematic.from_pdb_strict(contract["input"])
    assert mol.heavy_atoms == contract["expected"]["atom_count"]
    assert coords == contract["expected"]["coords"]
    with pytest.raises(ValueError, match=contract["malformed_error_field"]):
        chematic.from_pdb_strict(contract["malformed_input"])


def test_python_binding_matches_source_only_pdbqt_contract():
    contract = _DOCUMENT["pdbqt_contract"]
    mol = chematic.from_pdbqt(contract["input"])
    assert mol.heavy_atoms == contract["expected"]["atom_count"]
    with pytest.raises(ValueError):
        chematic.from_pdbqt(contract["malformed_input"])


def test_python_binding_matches_shared_qcschema_contract():
    contract = _DOCUMENT["qcschema_contract"]
    result = chematic.parse_qcschema_molecule(contract["input"])
    assert result["symbols"] == contract["expected"]["symbols"]
    assert result["mol"].formula == "H2O"
    assert len(result["symbols"]) == contract["expected"]["atom_count"]
    for observed, expected in zip(result["coords"], contract["expected"]["coords_angstrom"]):
        assert observed == pytest.approx(expected)
    assert result["molecular_charge"] == contract["expected"]["molecular_charge"]
    assert result["molecular_multiplicity"] == contract["expected"]["molecular_multiplicity"]


def test_python_binding_matches_shared_qcschema_atomic_input_contract():
    contract = _DOCUMENT["qcschema_atomic_input_contract"]
    result = chematic.parse_atomic_input(contract["input"])
    assert result["driver"] == contract["expected"]["driver"]
    assert result["model"]["method"] == contract["expected"]["model_method"]
    assert len(result["molecule"]["symbols"]) == contract["expected"]["atom_count"]
    assert result["extras"]["contract_marker"] == contract["expected"]["extra_marker"]
    assert result["vendor_extension"] == contract["expected"]["vendor_extension"]


def test_python_binding_matches_shared_qcschema_atomic_result_contract():
    contract = _DOCUMENT["qcschema_atomic_result_contract"]
    result = chematic.parse_atomic_result(contract["input"])
    assert result["driver"] == contract["expected"]["driver"]
    assert result["return_result"] == pytest.approx(contract["expected"]["return_result"])
    assert result["success"] is contract["expected"]["success"]
    assert len(result["molecule"]["symbols"]) == contract["expected"]["atom_count"]
    assert result["vendor_result_marker"] == contract["expected"]["vendor_result_marker"]


def test_python_binding_matches_shared_orca_input_contract():
    contract = _DOCUMENT["orca_input_contract"]
    result = chematic.parse_orca_input(contract["input"])
    assert result["keywords"] == contract["expected"]["keywords"]
    assert result["coords"]["kind"] == contract["expected"]["kind"]
    assert result["coords"]["charge"] == contract["expected"]["charge"]
    assert result["coords"]["multiplicity"] == contract["expected"]["multiplicity"]
    assert len(result["coords"]["atoms"]) == contract["expected"]["atom_count"]
    assert result["coords"]["atoms"][0]["element"] == contract["expected"]["element"]
    assert result["coords"]["atoms"][0]["x"] == pytest.approx(0.0)


def test_python_binding_matches_shared_orca_output_contract():
    contract = _DOCUMENT["orca_output_contract"]
    result = chematic.parse_orca_output(contract["input"])
    assert result["charge"] == contract["expected"]["charge"]
    assert result["multiplicity"] == contract["expected"]["multiplicity"]
    assert result["final_energy_hartree"] == pytest.approx(contract["expected"]["final_energy_hartree"])
    assert result["termination"]["kind"] == contract["expected"]["termination"]
    assert result["optimization_convergence"] == contract["expected"]["optimization_convergence"]
    assert len(result["trajectory"]) == contract["expected"]["trajectory_count"]
    with pytest.raises(ValueError, match=contract["malformed_error"]):
        chematic.parse_orca_output(contract["malformed_input"])


def test_python_binding_matches_shared_cdxml_contract():
    contract = _DOCUMENT["cdxml_contract"]
    mol = chematic.from_cdxml(contract["input"])
    assert mol.heavy_atoms == contract["expected"]["atom_count"]


def test_python_binding_matches_shared_cdxml_roundtrip_contract():
    contract = _DOCUMENT["cdxml_contract"]
    mol = chematic.from_cdxml(contract["input"])
    roundtripped = chematic.from_cdxml(mol.to_cdxml())
    assert roundtripped.heavy_atoms == contract["expected"]["atom_count"]


def test_python_binding_matches_shared_mmcif_contract():
    contract = _DOCUMENT["mmcif_contract"]
    result = chematic.parse_mmcif(contract["input"])
    assert len(result["atoms"]) == contract["expected"]["atom_count"]
    assert all(
        observed == pytest.approx(reference)
        for observed, reference in zip(result["coords"], contract["expected"]["coords"])
    )


def test_python_binding_matches_shared_moljson_contract():
    contract = _DOCUMENT["moljson_contract"]
    mol = chematic.from_moljson(contract["input"])
    assert mol.heavy_atoms == contract["expected"]["atom_count"]
    assert mol.smiles == contract["expected"]["canonical_smiles"]
    roundtripped = chematic.from_moljson(mol.to_moljson())
    assert roundtripped.heavy_atoms == contract["expected"]["atom_count"]
    assert roundtripped.smiles == contract["expected"]["canonical_smiles"]


def test_python_binding_matches_shared_rxn_document_contract():
    contract = _DOCUMENT["rxn_document_contract"]
    rxn = chematic.to_rxn_document_json(json.dumps(contract["document"]))
    decoded = json.loads(chematic.from_rxn_document_json(rxn))
    observed = [
        {"role": component["role"], "smiles": component["smiles"]}
        for component in decoded["steps"][0]["components"]
    ]
    assert observed == contract["expected_components"]


def test_python_binding_matches_shared_reaction_application_contract():
    contract = _DOCUMENT["reaction_application_contract"]
    reactants = [chematic.from_smiles(smiles) for smiles in contract["reactants"]]
    products = chematic.run_smirks(contract["smirks"], reactants)
    assert len(products) >= contract["expected"]["minimum_product_sets"]
    assert products[0][0].heavy_atoms == contract["expected"]["product_atom_count"]
    for case in contract["negative_cases"]:
        with pytest.raises(ValueError):
            chematic.run_smirks(
                case["smirks"],
                [chematic.from_smiles(smiles) for smiles in case["reactants"]],
            )
    for case in contract["additional_cases"]:
        reactants = [chematic.from_smiles(smiles) for smiles in case["reactants"]]
        products = chematic.run_smirks(case["smirks"], reactants)
        assert len(products) >= case["minimum_product_sets"]
        assert products[0][0].heavy_atoms == case["product_atom_count"]


@pytest.mark.parametrize("case", _DOCUMENT["reaction_smarts_contract"]["cases"])
def test_python_binding_matches_shared_reaction_smarts_contract(case):
    assert chematic.reaction_smarts_match(case["smarts"], case["reaction"]) == case["matches"]


def test_python_binding_matches_shared_reaction_balance_contract():
    for case in _DOCUMENT["reaction_balance_contract"]["cases"]:
        result = chematic.balance_check(case["reaction"])
        assert result["balanced"] == case["balanced"]
        assert result["diff"] == case["diff"]


def test_python_binding_matches_shared_reaction_center_contract():
    for case in _DOCUMENT["reaction_center_contract"]["cases"]:
        result = chematic.find_reaction_center(case["reaction"])
        assert result["broken_bonds"] == case["broken_bonds"]
        assert result["formed_bonds"] == case["formed_bonds"]
        assert result["changed_atoms"] == case["changed_atoms"]


@pytest.mark.parametrize(
    "fixture",
    _DOCUMENT["standardization_contract"]["fixtures"],
    ids=lambda item: item["id"],
)
def test_python_binding_matches_shared_standardization_profile(fixture):
    mol = chematic.from_smiles(fixture["smiles"])
    assert mol.standardize(largest_fragment_only=True).smiles == fixture["output_smiles"]


def test_python_binding_matches_shared_batch_canonicalization_contract():
    contract = _DOCUMENT["batch_canonicalization_contract"]
    actual = json.loads(chematic.canonicalize_smiles_batch_json(contract["inputs"]))
    assert actual["schema_version"] == 1
    assert actual["operation"] == "canonicalize_smiles"
    assert actual["status"] == "complete"
    assert actual["record_count"] == len(contract["expected"])
    observed = []
    for record in actual["records"]:
        result = {
            "input_index": record["input_index"],
            "input": record["input"],
            "status": record["status"],
        }
        if record["status"] == "accepted":
            result["canonical_smiles"] = record["canonical_smiles"]
        observed.append(result)
    assert observed == contract["expected"]


@pytest.mark.parametrize("format_name", ["xyz", "extxyz"])
def test_python_binding_matches_shared_xyz_batch_recovery_contract(tmp_path, format_name):
    contract = _DOCUMENT["xyz_batch_contract"]
    fixture = contract[format_name]
    path = tmp_path / f"recovery.{format_name}"
    path.write_text(fixture["input"])
    iterator = (
        chematic.iter_xyz_batched(str(path), batch_size=contract["batch_size"])
        if format_name == "xyz"
        else chematic.iter_extxyz_batched(str(path), batch_size=contract["batch_size"])
    )
    batches = list(iterator)
    expected = fixture["expected"]
    assert len(batches) == 1
    assert len(batches[0]) == 1
    manifest = json.loads(iterator.manifest_json())
    assert manifest["status"] == expected["status"]
    assert manifest["frames_seen"] == expected["record_count"]
    assert manifest["frames_emitted"] == 1
    assert manifest["rejected_frames"] == expected["rejected_count"]


@pytest.mark.parametrize(
    "fixture",
    _DOCUMENT["descriptor_contract"]["fixtures"],
    ids=lambda item: item["id"],
)
def test_python_binding_matches_shared_descriptor_contract(fixture):
    mol = chematic.from_smiles(fixture["smiles"])
    assert mol.mw == pytest.approx(fixture["molecular_weight"], abs=1e-6)
    assert mol.exact_mass == pytest.approx(fixture["exact_mass"], abs=1e-6)
    assert mol.formula == fixture["formula"]
    assert mol.tpsa == pytest.approx(fixture["tpsa"], abs=1e-6)
    assert mol.hbd == fixture["hbd"]
    assert mol.hba == fixture["hba"]
    assert mol.aromatic_ring_count == fixture["aromatic_ring_count"]
    assert mol.rotatable_bonds == fixture["rotatable_bonds"]
    assert mol.heavy_atoms == fixture["heavy_atoms"]


@pytest.mark.parametrize("fixture", _DOCUMENT["fixtures"], ids=lambda item: item["id"])
def test_python_binding_matches_shared_fixture(fixture):
    mol = chematic.from_smiles(fixture["smiles"])
    assert mol.smiles == fixture["canonical_smiles"]
    assert mol.heavy_atoms == fixture["heavy_atoms"]


@pytest.mark.parametrize(
    "fixture",
    _DOCUMENT["fingerprint_contract"]["fixtures"],
    ids=lambda item: item["id"],
)
def test_python_binding_matches_shared_fingerprint_shape(fixture):
    mol = chematic.from_smiles(fixture["smiles"])
    ecfp4 = mol.ecfp4()
    assert len(ecfp4) == 256
    assert [i for i in range(2048) if ecfp4[i // 8] & (1 << (i % 8))] == fixture["ecfp4_bits"]
    topo_path = mol.topo_path_fp()
    assert len(topo_path) == 256
    assert [i for i in range(2048) if topo_path[i // 8] & (1 << (i % 8))] == fixture["topo_path_bits"]
    torsion = mol.torsion_fp()
    assert len(torsion) == 256
    assert [i for i in range(2048) if torsion[i // 8] & (1 << (i % 8))] == fixture["torsion_bits"]
    rdkit_torsion = mol.rdkit_torsion_fp()
    assert len(rdkit_torsion) == 256
    assert [i for i in range(2048) if rdkit_torsion[i // 8] & (1 << (i % 8))] == fixture["rdkit_torsion_bits"]
    maccs = mol.maccs()
    assert len(maccs) == 21
    assert maccs.hex() == fixture["maccs_hex"]
    assert any(mol.ecfp4())
    assert any(maccs)


@pytest.mark.parametrize(
    "fixture",
    _DOCUMENT["fingerprint_contract"]["rdkit_rdk_fixtures"],
    ids=lambda item: item["id"],
)
def test_python_binding_matches_shared_rdkit_rdk_contract(fixture):
    mol = chematic.from_smiles(fixture["smiles"])
    actual = mol.rdkit_rdk_fp()
    assert len(actual) == 256
    assert [i for i in range(2048) if actual[i // 8] & (1 << (i % 8))] == fixture["rdkit_rdk_bits"]


@pytest.mark.parametrize(
    "fixture",
    _DOCUMENT["fingerprint_contract"]["rdkit_path_fixtures"],
    ids=lambda item: item["id"],
)
def test_python_binding_matches_shared_rdkit_path_contract(fixture):
    mol = chematic.from_smiles(fixture["smiles"])
    actual = mol.path_fp()
    assert len(actual) == 256
    assert [i for i in range(2048) if actual[i // 8] & (1 << (i % 8))] == fixture["rdkit_path_bits"]


@pytest.mark.parametrize(
    "fixture",
    _DOCUMENT["fingerprint_detail_contract"]["fixtures"],
    ids=lambda item: item["id"],
)
def test_python_binding_matches_shared_fingerprint_detail_contract(fixture):
    mol = chematic.from_smiles(fixture["smiles"])
    fp, sparse_counts, raw_bit_info, folded_bit_info = mol.rdkit_ecfp4_detail()
    assert len(fp) == 256
    assert sparse_counts
    assert raw_bit_info
    assert folded_bit_info
    for provenance in list(raw_bit_info.values()) + list(folded_bit_info.values()):
        for atom, radius in provenance:
            assert 0 <= atom < mol.heavy_atoms
            assert 0 <= radius <= 2


@pytest.mark.parametrize(
    "case",
    _DOCUMENT["adversarial"],
    ids=lambda item: item["id"],
)
def test_python_binding_rejects_shared_adversarial_fixture(case):
    parsers = {
        "smiles": chematic.from_smiles,
        "mol": chematic.from_mol_block,
        "mol_v3000": chematic.from_mol_v3000,
        "mol2": chematic.from_mol2,
        "cml": chematic.from_cml,
        "cjson": chematic.from_cjson,
        "moljson": chematic.from_moljson,
        "cdxml": chematic.from_cdxml,
    }
    with pytest.raises((ValueError, RuntimeError, TypeError)):
        parsers[case["format"]](case["input"])


def test_python_cdxml_document_edit_preserves_multi_page_presentation():
    cdxml = '<CDXML>\n<page id="p1">\n<arrow id="a1"/>\n</page>\n<page id="p2">\n<text id="t1"/>\n</page>\n</CDXML>'
    edited = chematic.edit_cdxml_document_json(
        cdxml,
        json.dumps({"kind": "set_page_attribute", "page_id": "p2", "key": "title", "value": "Page 2"}),
    )
    assert 'title="Page 2"' in edited
    assert '<arrow id="a1"/>' in edited


@pytest.mark.parametrize(
    "case",
    _DOCUMENT["semantic_expansion_contract"]["cases"],
    ids=lambda item: item["id"],
)
def test_python_semantic_expansion_matches_shared_contract(case):
    selected = chematic.semantic_apply_json_command(
        json.dumps(case["model"]), json.dumps(case["command"])
    )
    selected_model = json.loads(selected)
    if "expected_selected_alternative" in case:
        assert selected_model["r_groups"][0]["selected_alternative"] == case["expected_selected_alternative"]
    if "expected_repeat_count" in case:
        assert selected_model["polymer_units"][0]["repeat_count"] == case["expected_repeat_count"]
    expanded = json.loads(chematic.semantic_expand_json(case["base_smiles"], selected))
    assert expanded["schema"] == "chematic.semantic-expanded.v1"
    assert expanded["source_to_expanded"] == case["expected_source_to_expanded"]
    assert expanded["contracted_smiles"] == case["base_smiles"]


def test_python_rxn_document_contract_is_loss_aware():
    document = {
        "id": "rxn-contract",
        "steps": [{
            "id": "step-1",
            "components": [
                {"id": "reactant-1", "role": "reactant", "smiles": "CC", "coefficient": 1, "origin": "authored"},
                {"id": "product-1", "role": "product", "smiles": "CC", "coefficient": 1, "origin": "authored"},
            ],
            "conditions": [],
            "provenance": [],
            "origin": "authored",
        }],
        "provenance": [],
    }
    rxn = chematic.to_rxn_document_json(json.dumps(document))
    decoded = json.loads(chematic.from_rxn_document_json(rxn))
    assert [c["role"] for c in decoded["steps"][0]["components"]] == ["reactant", "product"]

    document["steps"][0]["conditions"] = [{"key": "temperature", "value": "25 C"}]
    with pytest.raises(ValueError):
        chematic.to_rxn_document_json(json.dumps(document))

    document["steps"][0]["conditions"] = []
    document["steps"][0]["components"][0]["smiles"] = "C>C"
    with pytest.raises(ValueError, match="invalid component"):
        chematic.to_rxn_document_json(json.dumps(document))
