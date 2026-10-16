use chematic_core::BondOrder;
use chematic_mol::{parse_mol2_record, parse_mol2_source_record, write_mol2_source_record};

const ACETATE: &str = "\
@<TRIPOS>MOLECULE
acetate
 4 3 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
1 C1 0.0 0.0 0.0 C.3 1 ACE 0.0
2 C2 1.5 0.0 0.0 C.2 1 ACE 0.0
3 O1 2.0 1.0 0.0 O.co2 1 ACE -0.5
4 O2 2.0 -1.0 0.0 O.co2 1 ACE -0.5
@<TRIPOS>UNITY_ATOM_ATTR
4 1
charge -1
@<TRIPOS>BOND
1 1 2 1
2 2 3 ar
3 2 4 ar
";

#[test]
fn typed_carboxylate_resonance_projects_to_formal_single_and_double_bonds() {
    let record = parse_mol2_record(ACETATE).unwrap();
    let source = parse_mol2_source_record(ACETATE).unwrap();
    for molecule in [&record.molecule, source.molecule()] {
        let orders: Vec<_> = molecule.bonds().map(|(_, bond)| bond.order).collect();
        assert_eq!(
            orders,
            [BondOrder::Single, BondOrder::Double, BondOrder::Single]
        );
    }
    assert_eq!(record.bonds[1].bond_type, "ar");
    assert_eq!(record.bonds[2].bond_type, "ar");
    assert_eq!(write_mol2_source_record(&source), ACETATE);
}

#[test]
fn ambiguous_or_untyped_oxygen_bonds_are_not_assigned_a_carboxylate_order() {
    for input in [
        ACETATE.replace("charge -1", "charge 0"),
        ACETATE.replacen("O.co2", "O.2", 1),
        ACETATE.replace("1 1 2 1", "1 1 2 ar"),
    ] {
        let record = parse_mol2_record(&input).unwrap();
        let source = parse_mol2_source_record(&input).unwrap();
        for molecule in [&record.molecule, source.molecule()] {
            assert_eq!(
                molecule.bonds().nth(1).unwrap().1.order,
                BondOrder::Aromatic
            );
            assert_eq!(
                molecule.bonds().nth(2).unwrap().1.order,
                BondOrder::Aromatic
            );
        }
    }
}
