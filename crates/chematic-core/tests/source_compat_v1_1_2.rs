//! Compile-time source-compatibility fixture for the public v1.1.2
//! `Chirality` enum. This is an integration test on purpose: it compiles as a
//! downstream crate and therefore catches patch-release additions that break
//! an exhaustive consumer match.

use chematic_core::{Chirality, SquarePlanarPermutation};

fn legacy_exhaustive_match(value: Chirality) -> u8 {
    match value {
        Chirality::None => 0,
        Chirality::CounterClockwise => 1,
        Chirality::Clockwise => 2,
        Chirality::SquarePlanar(SquarePlanarPermutation::SP1) => 3,
        Chirality::SquarePlanar(SquarePlanarPermutation::SP2) => 4,
        Chirality::SquarePlanar(SquarePlanarPermutation::SP3) => 5,
        Chirality::SquarePlanarUnnumbered => 6,
        Chirality::TrigonalBipyramidal(permutation) => 10 + permutation,
        Chirality::Octahedral(permutation) => 40 + permutation,
    }
}

#[test]
fn v1_1_2_exhaustive_match_still_compiles() {
    assert_eq!(legacy_exhaustive_match(Chirality::Clockwise), 2);
    assert_eq!(Chirality::Allene1, Chirality::CounterClockwise);
    assert_eq!(Chirality::Allene2, Chirality::Clockwise);
}
