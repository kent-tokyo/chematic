use chematic_depict::{RenderOptions, depict_svg_opts};
use chematic_smiles::parse;

fn svg(smiles: &str, kekulize: bool) -> String {
    depict_svg_opts(
        &parse(smiles).unwrap(),
        &RenderOptions {
            kekulize,
            ..Default::default()
        },
    )
}

#[test]
fn bonded_labels_typeset_charge_hydrogen_count_and_isotope() {
    for (smiles, expected) in [
        ("[O-]S(=O)(=O)[O-]", "O⁻"),
        ("C[NH3+]", "NH₃⁺"),
        ("C[N+](C)(C)C", "N⁺"),
        ("C[S+2]C", "S²⁺"),
        ("C[N-2]C", "N²⁻"),
        ("NCC", "NH₂"),
        ("C[13CH2]O", "¹³C"),
        ("C[15NH3+]", "¹⁵NH₃⁺"),
        ("[O-][n+]1ccccc1", "N⁺"),
    ] {
        for kekulize in [false, true] {
            let svg = svg(smiles, kekulize);
            assert!(
                svg.contains(&format!(">{expected}</text>")),
                "{smiles}: {svg}"
            );
        }
    }
}

#[test]
fn bonded_and_isolated_ions_share_one_style() {
    let svg = svg("[Na+].CC(=O)[O-]", false);
    assert!(svg.contains(">Na⁺</text>"), "{svg}");
    assert!(svg.contains(">O⁻</text>"), "{svg}");
    for ascii in [">O-</text>", ">Na+</text>"] {
        assert!(!svg.contains(ascii), "{ascii}: {svg}");
    }
}

#[test]
fn unchanged_labels_stay_as_they_were() {
    for (smiles, expected) in [("OCC", "OH"), ("c1ccncc1", "N"), ("CC(=O)O", "O")] {
        let svg = svg(smiles, false);
        assert!(
            svg.contains(&format!(">{expected}</text>")),
            "{smiles}: {svg}"
        );
    }
    // Plain interior carbons still get no label.
    assert!(!svg("CCC", false).contains("</text>"));
}
