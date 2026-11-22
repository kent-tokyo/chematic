use chematic_cip::{AtomicNumberKey, RationalAtomicNumber};
use std::cmp::Ordering;

#[test]
fn mancude_means_sort_exactly_between_integer_atomic_numbers() {
    let six = RationalAtomicNumber::integer(6);
    let one_third = RationalAtomicNumber::mean(&[6, 6, 7]);
    let one_half = RationalAtomicNumber::mean(&[6, 7]);
    let two_thirds = RationalAtomicNumber::mean(&[6, 7, 7]);
    let seven = RationalAtomicNumber::integer(7);
    assert_eq!((one_third.numerator(), one_third.denominator()), (19, 3));
    let mut values = [seven, one_half, six, two_thirds, one_third];
    values.sort();
    assert_eq!(values, [six, one_third, one_half, two_thirds, seven]);
    assert_eq!(one_third.partial_cmp(&one_half), Some(Ordering::Less));
    assert_eq!(seven.cmp(&two_thirds), Ordering::Greater);

    let mut keys = [
        AtomicNumberKey::Integral(7),
        AtomicNumberKey::Rational(one_half),
        AtomicNumberKey::Integral(6),
    ];
    keys.sort();
    assert_eq!(
        keys,
        [
            AtomicNumberKey::Integral(6),
            AtomicNumberKey::Rational(one_half),
            AtomicNumberKey::Integral(7)
        ]
    );
    assert!(keys[0] < keys[1]);
    assert_eq!(keys[2].cmp(&keys[1]), Ordering::Greater);
}
