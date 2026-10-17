use super::{BondOrder, push_right_aligned_u32, push_v2000_bond_tail, v2000_bond_type};

#[test]
fn every_bond_kind_and_stereo_code_matches_independent_format_and_old_appends() {
    use BondOrder::*;
    for order in [
        Single,
        Up,
        Down,
        Dative,
        Double,
        Triple,
        Aromatic,
        Quadruple,
        QuerySingleOrDouble,
        QuerySingleOrAromatic,
        QueryDoubleOrAromatic,
        QueryAny,
        Zero,
    ] {
        for stereo in 0..=6 {
            let mut actual = String::from("prefix");
            push_v2000_bond_tail(&mut actual, order, stereo);
            let mut old = String::from("prefix");
            push_right_aligned_u32(&mut old, v2000_bond_type(order), 3);
            push_right_aligned_u32(&mut old, stereo, 3);
            old.push('\n');
            assert_eq!(actual, old);
            assert_eq!(
                actual,
                format!("prefix{:>3}{stereo:>3}\n", v2000_bond_type(order))
            );
        }
    }
}
