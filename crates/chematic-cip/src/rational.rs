//! Exact rational representation of a "mean atomic number".
//!
//! IUPAC's MANCUDE (maximum non-cumulated double bonds) treatment averages a ring atom's
//! atomic number across every valid Kekulé placement of its ring system, and that average
//! is not always an integer -- a heteroatom-adjacent position can land on a genuine
//! fraction like 6⅓ or 6½. Representing that as `f64` would reintroduce the exact kind of
//! determinism/order-dependence hazard this project has been repeatedly bitten by in other
//! areas (see `docs/rfcs/cip_accurate_rfc.md`'s earlier rounds) -- this type keeps the value as
//! an exact integer fraction instead, always stored in lowest terms so structural equality
//! (`#[derive(Eq, PartialEq)]`) matches mathematical equality: `mean(&[6, 6])` and
//! `integer(6)` both reduce to `6/1` and compare equal, not "equal in value but different
//! in representation."
//!
//! Milestone 3B-0 scope: this type is designed and unit-tested in isolation here. It is
//! not yet wired into [`crate::node::CipNodeKind`] or the comparator -- that's Milestone
//! 3B-1's job, once the corpus classification and digraph-diff evidence this milestone
//! collects has been reviewed.

use std::cmp::Ordering;

/// An exact, always-reduced rational atomic number: `numerator / denominator`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct RationalAtomicNumber {
    numerator: u32,
    denominator: u32,
}

impl RationalAtomicNumber {
    /// A plain integer atomic number (denominator 1) -- the degenerate case of `mean` over
    /// identical values, and the common case for atoms outside any MANCUDE ring system.
    pub fn integer(n: u32) -> Self {
        Self {
            numerator: n,
            denominator: 1,
        }
    }

    /// The IUPAC MANCUDE mean: `sum(values) / count`, reduced to lowest terms.
    ///
    /// # Panics
    /// Panics if `values` is empty. A mean over zero valid Kekulé placements is never a
    /// value CIP comparison legitimately needs -- every MANCUDE component has at least one
    /// valid placement by construction, or `enumerate_kekule_matchings` already errored
    /// before this would be called.
    pub fn mean(values: &[u32]) -> Self {
        assert!(
            !values.is_empty(),
            "RationalAtomicNumber::mean of an empty slice"
        );
        let sum: u64 = values.iter().map(|&v| u64::from(v)).sum();
        Self::reduced(sum, values.len() as u64)
    }

    fn reduced(numerator: u64, denominator: u64) -> Self {
        debug_assert!(denominator != 0, "RationalAtomicNumber: zero denominator");
        let g = gcd(numerator, denominator);
        Self {
            numerator: u32::try_from(numerator / g).expect("atomic number sum overflowed u32"),
            denominator: u32::try_from(denominator / g)
                .expect("atomic number count overflowed u32"),
        }
    }

    pub fn numerator(self) -> u32 {
        self.numerator
    }

    pub fn denominator(self) -> u32 {
        self.denominator
    }
}

fn gcd(a: u64, b: u64) -> u64 {
    if b == 0 { a } else { gcd(b, a % b) }
}

/// Compares two rational atomic numbers via cross-multiplication with a `u64`
/// intermediate: `a/b < c/d` iff `a*d < c*b` (both denominators are always positive, so
/// the inequality direction is preserved -- no float division needed).
pub fn cmp_atomic_number(a: RationalAtomicNumber, b: RationalAtomicNumber) -> Ordering {
    let lhs = u64::from(a.numerator) * u64::from(b.denominator);
    let rhs = u64::from(b.numerator) * u64::from(a.denominator);
    lhs.cmp(&rhs)
}

impl Ord for RationalAtomicNumber {
    fn cmp(&self, other: &Self) -> Ordering {
        cmp_atomic_number(*self, *other)
    }
}

impl PartialOrd for RationalAtomicNumber {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl core::fmt::Display for RationalAtomicNumber {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        if self.denominator == 1 {
            write!(f, "{}", self.numerator)
        } else {
            write!(f, "{}/{}", self.numerator, self.denominator)
        }
    }
}

/// The atomic number a [`crate::node::CipNode`] compares by: either a plain integer (every
/// node outside a MANCUDE resonance component, and always for real `Atom` nodes -- MANCUDE
/// only ever touches `MultipleBondDuplicate` nodes, never real atoms) or a
/// [`RationalAtomicNumber`] (a `MultipleBondDuplicate` whose *owner* atom sits in a MANCUDE
/// component -- see `mancude.rs`'s module docs for why the owner's, not the represented
/// atom's, fraction is the correct value).
///
/// `Integral` holds a `u8` to match [`chematic_core::Element::atomic_number`]'s existing
/// return type, not a wider int -- there is no atomic number this project needs to
/// represent that doesn't already fit in `u8`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AtomicNumberKey {
    Integral(u8),
    Rational(RationalAtomicNumber),
}

impl AtomicNumberKey {
    fn as_rational(self) -> RationalAtomicNumber {
        match self {
            AtomicNumberKey::Integral(n) => RationalAtomicNumber::integer(u32::from(n)),
            AtomicNumberKey::Rational(r) => r,
        }
    }
}

/// Compares two keys by promoting an `Integral` to its equivalent `RationalAtomicNumber`
/// (denominator 1) when either side is `Rational`, then delegating to
/// [`cmp_atomic_number`]. `Integral` vs `Integral` never leaves plain integer comparison
/// (the promotion is exact -- `n` and `n/1` compare identically either way).
pub fn cmp_atomic_number_key(a: AtomicNumberKey, b: AtomicNumberKey) -> Ordering {
    if let (AtomicNumberKey::Integral(x), AtomicNumberKey::Integral(y)) = (a, b) {
        return x.cmp(&y);
    }
    cmp_atomic_number(a.as_rational(), b.as_rational())
}

impl Ord for AtomicNumberKey {
    fn cmp(&self, other: &Self) -> Ordering {
        cmp_atomic_number_key(*self, *other)
    }
}

impl PartialOrd for AtomicNumberKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl core::fmt::Display for AtomicNumberKey {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            AtomicNumberKey::Integral(n) => write!(f, "{n}"),
            AtomicNumberKey::Rational(r) => write!(f, "{r}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mean_of_identical_values_reduces_to_integer() {
        assert_eq!(
            RationalAtomicNumber::mean(&[6, 6, 6]),
            RationalAtomicNumber::integer(6)
        );
    }

    #[test]
    fn mean_reduces_to_lowest_terms() {
        // 12/2 must store as 6/1, not 12/2 -- structural equality depends on this.
        let m = RationalAtomicNumber::mean(&[6, 6]);
        assert_eq!(m.numerator(), 6);
        assert_eq!(m.denominator(), 1);
    }

    #[test]
    fn fraction_examples_from_the_design_conversation() {
        // 6 + 7 + 7 over 3 placements = 20/3 = 6 2/3 -- NOT the 6 1/3 someone might guess
        // from "one heteroatom among three forms"; worked by hand: (6+7+7)/3 = 20/3.
        let six_and_two_thirds = RationalAtomicNumber::mean(&[6, 7, 7]);
        assert_eq!(six_and_two_thirds.numerator(), 20);
        assert_eq!(six_and_two_thirds.denominator(), 3);

        // 6 + 7 over 2 placements = 13/2 = 6 1/2.
        let six_and_a_half = RationalAtomicNumber::mean(&[6, 7]);
        assert_eq!(six_and_a_half.numerator(), 13);
        assert_eq!(six_and_a_half.denominator(), 2);
    }

    #[test]
    fn cross_multiply_comparison() {
        let half = RationalAtomicNumber {
            numerator: 1,
            denominator: 2,
        };
        let third = RationalAtomicNumber {
            numerator: 1,
            denominator: 3,
        };
        assert_eq!(cmp_atomic_number(third, half), Ordering::Less);
        assert_eq!(cmp_atomic_number(half, third), Ordering::Greater);

        let six = RationalAtomicNumber::integer(6);
        let six_and_a_half = RationalAtomicNumber::mean(&[6, 7]);
        assert_eq!(cmp_atomic_number(six, six_and_a_half), Ordering::Less);
        assert_eq!(cmp_atomic_number(six_and_a_half, six), Ordering::Greater);
    }

    #[test]
    fn reflexive_and_antisymmetric() {
        let values: &[RationalAtomicNumber] = &[
            RationalAtomicNumber::integer(6),
            RationalAtomicNumber::mean(&[6, 7]),
            RationalAtomicNumber::mean(&[6, 7, 7]),
            RationalAtomicNumber::mean(&[7, 7]),
        ];
        for &a in values {
            assert_eq!(
                cmp_atomic_number(a, a),
                Ordering::Equal,
                "reflexivity: {a:?}"
            );
            for &b in values {
                let ab = cmp_atomic_number(a, b);
                let ba = cmp_atomic_number(b, a);
                assert_eq!(ab, ba.reverse(), "antisymmetry: {a:?} vs {b:?}");
            }
        }
    }

    #[test]
    fn display_formats_integers_and_fractions() {
        assert_eq!(RationalAtomicNumber::integer(6).to_string(), "6");
        assert_eq!(RationalAtomicNumber::mean(&[6, 7]).to_string(), "13/2");
    }

    #[test]
    fn atomic_number_key_integral_vs_integral_matches_plain_u8_ordering() {
        assert_eq!(
            cmp_atomic_number_key(AtomicNumberKey::Integral(6), AtomicNumberKey::Integral(7)),
            Ordering::Less
        );
        assert_eq!(
            cmp_atomic_number_key(AtomicNumberKey::Integral(7), AtomicNumberKey::Integral(6)),
            Ordering::Greater
        );
        assert_eq!(
            cmp_atomic_number_key(AtomicNumberKey::Integral(6), AtomicNumberKey::Integral(6)),
            Ordering::Equal
        );
    }

    #[test]
    fn atomic_number_key_integral_vs_rational_promotion_at_the_divergence_table_boundaries() {
        // 6 vs 6.5 (pyridine's N-adjacent carbon fraction) -- integral loses.
        let six = AtomicNumberKey::Integral(6);
        let six_and_a_half = AtomicNumberKey::Rational(RationalAtomicNumber::mean(&[6, 7]));
        assert_eq!(cmp_atomic_number_key(six, six_and_a_half), Ordering::Less);
        assert_eq!(
            cmp_atomic_number_key(six_and_a_half, six),
            Ordering::Greater
        );

        // 6 vs 6.333 (quinoline's oracle-side value, 19/3) -- integral still loses, and 7
        // (nitrogen) still beats both.
        let six_and_a_third = AtomicNumberKey::Rational(RationalAtomicNumber::mean(&[6, 6, 7]));
        assert_eq!(cmp_atomic_number_key(six, six_and_a_third), Ordering::Less);
        let seven = AtomicNumberKey::Integral(7);
        assert_eq!(
            cmp_atomic_number_key(six_and_a_third, seven),
            Ordering::Less
        );

        // Exact-integer rational (denominator 1) must compare equal to the matching
        // Integral value -- the promotion must be lossless, not just order-preserving.
        let six_rational = AtomicNumberKey::Rational(RationalAtomicNumber::integer(6));
        assert_eq!(cmp_atomic_number_key(six, six_rational), Ordering::Equal);
    }

    #[test]
    fn atomic_number_key_reflexive_and_antisymmetric() {
        let values: &[AtomicNumberKey] = &[
            AtomicNumberKey::Integral(1),
            AtomicNumberKey::Integral(6),
            AtomicNumberKey::Integral(7),
            AtomicNumberKey::Rational(RationalAtomicNumber::integer(6)),
            AtomicNumberKey::Rational(RationalAtomicNumber::mean(&[6, 7])),
            AtomicNumberKey::Rational(RationalAtomicNumber::mean(&[6, 6, 7])),
        ];
        for &a in values {
            assert_eq!(
                cmp_atomic_number_key(a, a),
                Ordering::Equal,
                "reflexivity: {a:?}"
            );
            for &b in values {
                let ab = cmp_atomic_number_key(a, b);
                let ba = cmp_atomic_number_key(b, a);
                assert_eq!(ab, ba.reverse(), "antisymmetry: {a:?} vs {b:?}");
            }
        }
    }

    #[test]
    fn atomic_number_key_display() {
        assert_eq!(AtomicNumberKey::Integral(6).to_string(), "6");
        assert_eq!(
            AtomicNumberKey::Rational(RationalAtomicNumber::mean(&[6, 7])).to_string(),
            "13/2"
        );
    }
}
