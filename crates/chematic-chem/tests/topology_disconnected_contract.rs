use chematic_chem::{padmakar_ivan_index, wiener_index};
use chematic_smiles::parse;

#[test]
fn disconnected_counterion_and_explicit_hydrogen_do_not_change_heavy_graph_indices() {
    // Propane has heavy-atom distances 1, 1 and 2, and PI index 6.
    for smiles in ["CCC", "[H]CCC", "CCC.[Na+]", "[H]CCC.[Na+]"] {
        let molecule = parse(smiles).unwrap();
        assert_eq!(wiener_index(&molecule), 4.0, "{smiles}");
        assert_eq!(padmakar_ivan_index(&molecule), 6, "{smiles}");
    }
}
