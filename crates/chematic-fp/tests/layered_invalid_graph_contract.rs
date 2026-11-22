use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};
use chematic_fp::{RdkitLayeredConfig, rdkit_layered_config};

#[test]
fn layered_fingerprint_refuses_invalid_resource_caps_and_self_bonds() {
    let molecule = chematic_smiles::parse("CCO").unwrap();
    for config in [
        RdkitLayeredConfig {
            max_visits: 0,
            ..Default::default()
        },
        RdkitLayeredConfig {
            max_visits: 10_000_001,
            ..Default::default()
        },
        RdkitLayeredConfig {
            max_paths: 0,
            ..Default::default()
        },
        RdkitLayeredConfig {
            max_paths: 1_000_001,
            ..Default::default()
        },
    ] {
        let error = rdkit_layered_config(&molecule, &config).unwrap_err();
        assert_eq!(error.code, "resource_limit");
    }
    let mut builder = MoleculeBuilder::new();
    let atom = builder.add_atom(Atom::new(Element::C));
    builder.add_bond(atom, atom, BondOrder::Single).unwrap();
    let self_bond = builder.build();
    let error = rdkit_layered_config(&self_bond, &RdkitLayeredConfig::default()).unwrap_err();
    assert_eq!(error.code, "unsupported_graph");
    assert_eq!(self_bond.bond_count(), 1);
}
