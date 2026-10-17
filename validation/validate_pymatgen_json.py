"""Independent fixed-version pymatgen serialization/reader comparison."""

import hashlib
from importlib.metadata import version
import json

import numpy as np
from pymatgen.core import Element, Lattice, Structure
from chematic import parse_pymatgen_structures


def main():
    assert version("pymatgen") == "2025.6.14"
    rng = np.random.default_rng(2675)
    cases = sites_checked = species_checked = response_cases = 0
    max_coordinate_difference = 0.0
    checksum = hashlib.sha256()
    for n in [1, 2, 7, 31]:
        for lattice in [
            Lattice.cubic(5.25),
            Lattice.orthorhombic(4.1, 6.2, 8.3),
            Lattice.from_parameters(5.0, 6.0, 7.0, 75.0, 86.0, 110.0),
            Lattice([[-4.0, 0.5, 0.1], [0.2, 5.0, 0.7], [0.8, 0.1, 6.0]]),
        ]:
            for sample in range(6):
                coords = rng.uniform(-3.0, 4.0, (n, 3))
                labels = [f"synthetic α {i}" for i in range(n)]
                for profile in ["ordered", "disordered", "oxidation", "properties"]:
                    if profile == "disordered":
                        species = [{"Fe": 0.4, "Ni": 0.5} for _ in range(n)]
                    elif profile == "oxidation":
                        species = [
                            {"Fe2+": 0.4, "Fe3+": 0.6} if i % 2 == 0 else {"O2-": 1.0}
                            for i in range(n)
                        ]
                    else:
                        species = [["C", "O", "Cl", "Si"][i % 4] for i in range(n)]
                    props = (
                        {
                            "magmom": rng.uniform(-2.0, 2.0, n).tolist(),
                            "selective_dynamics": [
                                [True, False, True] for _ in range(n)
                            ],
                        }
                        if profile == "properties"
                        else None
                    )
                    model = Structure(
                        lattice,
                        species,
                        coords,
                        labels=labels,
                        site_properties=props,
                        properties=(
                            {"synthetic_profile": profile}
                            if profile == "properties"
                            else {}
                        ),
                    )
                    data = model.as_dict()
                    payload = json.dumps(data, ensure_ascii=False)
                    checksum.update(payload.encode())
                    record = parse_pymatgen_structures(payload.encode())[0]
                    assert record.source_text == payload
                    assert record.structure_dict() == json.loads(payload)
                    native = record.to_periodic_structure(allow_metadata_loss=True)
                    np.testing.assert_array_equal(
                        native.lattice.matrix, model.lattice.matrix
                    )
                    np.testing.assert_array_equal(
                        native.fractional_positions(), model.frac_coords
                    )
                    delta = float(
                        np.abs(native.cartesian_positions() - model.cart_coords).max()
                    )
                    max_coordinate_difference = max(max_coordinate_difference, delta)
                    np.testing.assert_allclose(
                        native.cartesian_positions(),
                        model.cart_coords,
                        rtol=0,
                        atol=1e-12,
                    )
                    assert [s.label for s in native.sites] == labels
                    for native_site, original_site in zip(native.sites, data["sites"]):
                        assert native_site.species == [
                            (s["element"], float(s["occu"]))
                            for s in original_site["species"]
                        ]
                        species_checked += len(native_site.species)
                    # Unmodified record metadata remains consumable by pymatgen,
                    # including separate oxidation states on the same element.
                    reread = Structure.from_dict(record.structure_dict())
                    assert reread.as_dict() == model.as_dict()
                    wrapped = json.dumps(
                        {
                            "data": [
                                {"material_id": "synthetic-a", "structure": data},
                                {"material_id": "synthetic-b", "structure": data},
                            ],
                            "meta": {"total_doc": 2},
                        }
                    )
                    records = parse_pymatgen_structures(wrapped)
                    assert len(records) == 2
                    for i, item in enumerate(records):
                        assert item.record_index == i and item.source_text == wrapped
                        assert item.structure_dict() == json.loads(payload)
                        assert f"data[{i}].material_id" in item.projection_losses
                        assert "response" in item.projection_losses
                    response_cases += 1
                    cases += 1
                    sites_checked += n
    elements = [Element.from_Z(i).symbol for i in range(1, 119)]
    all_elements = Structure(
        Lattice.cubic(100.0), elements, rng.uniform(-1.0, 2.0, (118, 3))
    )
    all_elements_record = parse_pymatgen_structures(json.dumps(all_elements.as_dict()))[
        0
    ]
    assert [
        s.species[0][0]
        for s in all_elements_record.to_periodic_structure(
            allow_metadata_loss=True
        ).sites
    ] == elements
    print(
        json.dumps(
            {
                "pymatgen_version": version("pymatgen"),
                "numpy_version": np.__version__,
                "seed": 2675,
                "additional_element_symbols_checked": 118,
                "structures": cases,
                "sites_checked": sites_checked,
                "species_entries_checked": species_checked,
                "two_record_responses": response_cases,
                "matrices_fractional_coordinates_labels_species": "exact equality",
                "metadata_dictionary_reread": "pymatgen Structure.from_dict(...).as_dict() equals original",
                "maximum_cartesian_absolute_difference_angstrom": max_coordinate_difference,
                "cartesian_atol_angstrom": 1e-12,
                "input_sha256": checksum.hexdigest(),
                "scope": "synthetic explicit periodic cells; no Materials Project request, symmetry inference, chemical bonds, or format-wide parity claim",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
