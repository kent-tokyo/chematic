use chematic_cip::{CipBudget, assign_cip_accurate_experimental};
use chematic_core::{AtomIdx, CipCode};

#[test]
fn sulfoxide_auxiliaries_distinguish_a_pseudoasymmetric_carbon() {
    for (smiles, sulfur, carbon) in [
        (
            "O=[S@](CC)[C@H](O)[S@](CC)=O",
            [CipCode::R, CipCode::S],
            Some(CipCode::LowerR),
        ),
        (
            "O=[S@@](CC)[C@@H](O)[S@@](CC)=O",
            [CipCode::S, CipCode::R],
            Some(CipCode::LowerR),
        ),
        (
            "O=[S@](CC)[C@@H](O)[S@](CC)=O",
            [CipCode::R, CipCode::S],
            Some(CipCode::LowerS),
        ),
        (
            "O=[S@](CC)[C@H](O)[S@@](CC)=O",
            [CipCode::R, CipCode::R],
            None,
        ),
    ] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let result = assign_cip_accurate_experimental(&mol, CipBudget::default_budget()).unwrap();
        let code = |atom| {
            result
                .assignments
                .iter()
                .find(|(idx, _)| *idx == AtomIdx(atom))
                .map(|&(_, code)| code)
        };
        assert_eq!([code(1), code(6)], sulfur.map(Some), "{smiles}");
        assert_eq!(code(4), carbon, "{smiles}");
    }
}
