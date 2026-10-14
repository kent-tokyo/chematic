//! Open Babel #2655: explicit H neighbors and total CML hydrogenCount.
use chematic_core::{AtomIdx, Element, Molecule, implicit_hcount};
use chematic_mol::{parse_cml, parse_cml_strict, write_cml};
use chematic_smiles::{canonical_smiles, parse};

const CH2: &str = include_str!("../../../validation/openbabel_issues/cml/explicit-ch2.cml");

fn count_attribute(count: &str) -> String {
    CH2.replace(
        "id=\"carbon\" elementType=\"C\"",
        &format!("id=\"carbon\" elementType=\"C\" hydrogenCount=\"{count}\""),
    )
}

fn conserved(a: &Molecule, b: &Molecule) {
    assert_eq!(a.atom_count(), b.atom_count());
    assert_eq!(a.bond_count(), b.bond_count());
    assert_eq!(a.total_formula(), b.total_formula());
    for (idx, atom) in a.atoms() {
        let other = b.atom(idx);
        assert_eq!(atom.element, other.element);
        assert_eq!(atom.isotope, other.isotope);
        assert_eq!(atom.charge, other.charge);
        assert_eq!(implicit_hcount(a, idx), implicit_hcount(b, idx));
    }
}

#[test]
fn openbabel_2655_explicit_ch2_does_not_become_methane() {
    for parse_cml in [parse_cml, parse_cml_strict] {
        let (mol, _) = parse_cml(CH2).unwrap();
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.bond_count(), 2);
        assert_eq!(implicit_hcount(&mol, AtomIdx(0)), 0);
        assert_eq!(mol.total_formula(), "CH2");
        let smiles = canonical_smiles(&mol);
        assert_eq!(parse(&smiles).unwrap().total_formula(), "CH2");
    }
}

#[test]
fn total_hydrogen_count_subtracts_explicit_neighbors() {
    for (total, implicit, formula) in [(2, 0, "CH2"), (3, 1, "CH3"), (4, 2, "CH4")] {
        for reader in [parse_cml, parse_cml_strict] {
            let (mol, _) = reader(&count_attribute(&total.to_string())).unwrap();
            assert_eq!(implicit_hcount(&mol, AtomIdx(0)), implicit);
            assert_eq!(mol.total_formula(), formula);
        }
    }
}

#[test]
fn isotopic_hydrogens_count_as_explicit_neighbors_and_keep_labels() {
    for input in [CH2.to_owned(), count_attribute("3")] {
        let input = input.replace(
            "id=\"hydrogen7\" elementType=\"H\"",
            "id=\"hydrogen7\" elementType=\"H\" isotopeNumber=\"2\"",
        );
        let (mol, _) = parse_cml(&input).unwrap();
        assert_eq!(mol.atom(AtomIdx(1)).element, Element::H);
        assert_eq!(mol.atom(AtomIdx(1)).isotope, Some(2));
        let expected = if input.contains("hydrogenCount") {
            1
        } else {
            0
        };
        assert_eq!(implicit_hcount(&mol, AtomIdx(0)), expected);
        let (again, _) = parse_cml_strict(&write_cml(&mol, None)).unwrap();
        conserved(&mol, &again);
    }
}

#[test]
fn implicit_only_input_retains_normal_valence_inference() {
    for (element, charge, formula) in [("C", 0, "CH4"), ("O", 0, "H2O"), ("N", 1, "H4N")] {
        let input = format!(
            "<molecule>\n<atom id=\"a\" elementType=\"{element}\" formalCharge=\"{charge}\"/>\n</molecule>\n"
        );
        let (mol, _) = parse_cml(&input).unwrap();
        assert_eq!(mol.total_formula(), formula);
        let (again, _) = parse_cml(&write_cml(&mol, None)).unwrap();
        conserved(&mol, &again);
    }
}

#[test]
fn writer_keeps_explicit_and_implicit_h_mixtures_and_zero_h_counts() {
    for smiles in [
        "[H]C",
        "[H][C]",
        "[CH2]",
        "[13CH2][2H]",
        "[NH3+][H]",
        "[H]OC",
        "[H]c1ccccc1",
        "[H][H]",
        "[H]",
    ] {
        let mol = parse(smiles).unwrap();
        let written = write_cml(&mol, None);
        let (again, _) = parse_cml_strict(&written).unwrap();
        conserved(&mol, &again);
        let second = write_cml(&again, None);
        assert_eq!(written, second, "writer is not a fixed point: {smiles}");
    }
    let written = write_cml(&parse("[H]C").unwrap(), None);
    assert!(
        written
            .lines()
            .any(|line| line.contains("elementType=\"C\"") && line.contains("hydrogenCount=\"4\""))
    );
}

#[test]
fn malformed_or_contradictory_counts_are_rejected() {
    for count in [
        "-1",
        "",
        "NaN",
        "3.5",
        "bad",
        "0",
        "1",
        "258",
        "184467440737095516160",
    ] {
        for reader in [parse_cml, parse_cml_strict] {
            assert!(
                reader(&count_attribute(count)).is_err(),
                "accepted hydrogenCount={count}"
            );
        }
    }
}

#[test]
fn total_count_can_exceed_u8_when_explicit_neighbors_cover_the_difference() {
    let (mol, _) = parse_cml(&count_attribute("257")).unwrap();
    assert_eq!(implicit_hcount(&mol, AtomIdx(0)), 255);
    let (again, _) = parse_cml(&write_cml(&mol, None)).unwrap();
    conserved(&mol, &again);
}

#[test]
fn hydrogen_reconciliation_is_independent_of_atom_insertion_order() {
    let input = "<molecule>\n<atom id=\"h1\" elementType=\"H\"/>\n<atom id=\"c\" elementType=\"C\" hydrogenCount=\"3\"/>\n<atom id=\"h2\" elementType=\"H\"/>\n<bond atomRefs2=\"c h2\" order=\"1\"/>\n<bond atomRefs2=\"h1 c\" order=\"1\"/>\n</molecule>";
    let (mol, _) = parse_cml(input).unwrap();
    assert_eq!(mol.total_formula(), "CH3");
    assert_eq!(implicit_hcount(&mol, AtomIdx(1)), 1);
}
