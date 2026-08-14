//! Coordinate newtypes, occupancy, and the periodic site type.

use crate::error::CrystalError;
use crate::validation::require_finite3;
use chematic_core::Element;

// ---------------------------------------------------------------------------
// Coordinate newtypes
// ---------------------------------------------------------------------------

/// Fractional (lattice-relative, dimensionless) coordinates.
///
/// Not range-restricted by the type itself -- a raw fractional triple may
/// legally sit outside `[0, 1)` (e.g. mid-calculation, before wrapping).
/// Use [`FractionalCoord::wrapped`] to reduce into `[0, 1)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FractionalCoord(pub [f64; 3]);

impl FractionalCoord {
    /// Construct from raw components (no validation -- see [`Self::is_finite`]).
    #[inline]
    pub fn new(value: [f64; 3]) -> Self {
        Self(value)
    }

    /// `true` if every component is finite (no `NaN`, no `+-Infinity`).
    #[inline]
    pub fn is_finite(&self) -> bool {
        self.0.iter().all(|c| c.is_finite())
    }

    /// Reduce each component into `[0, 1)` via `rem_euclid(1.0)`.
    ///
    /// `1.0` itself maps to `0.0`; negative values wrap up (e.g. `-0.25` ->
    /// `0.75`), matching the usual crystallographic convention for "the
    /// representative image inside the unit cell". Non-finite input yields
    /// non-finite output (`rem_euclid` propagates `NaN`) -- callers that
    /// need a hard guarantee should validate first via [`Self::is_finite`].
    pub fn wrapped(&self) -> Self {
        Self([
            self.0[0].rem_euclid(1.0),
            self.0[1].rem_euclid(1.0),
            self.0[2].rem_euclid(1.0),
        ])
    }

    /// Translate by an integer lattice-image shift `[i, j, k]`.
    pub fn translated(&self, image: [i32; 3]) -> Self {
        Self([
            self.0[0] + f64::from(image[0]),
            self.0[1] + f64::from(image[1]),
            self.0[2] + f64::from(image[2]),
        ])
    }
}

/// Cartesian (orthogonal, Angstrom) coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CartesianCoord(pub [f64; 3]);

impl CartesianCoord {
    /// Construct from raw components (no validation -- see [`Self::is_finite`]).
    #[inline]
    pub fn new(value: [f64; 3]) -> Self {
        Self(value)
    }

    /// `true` if every component is finite (no `NaN`, no `+-Infinity`).
    #[inline]
    pub fn is_finite(&self) -> bool {
        self.0.iter().all(|c| c.is_finite())
    }

    /// Euclidean distance to another Cartesian point, in Angstrom.
    pub fn distance(&self, other: &CartesianCoord) -> f64 {
        let d = [
            self.0[0] - other.0[0],
            self.0[1] - other.0[1],
            self.0[2] - other.0[2],
        ];
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
    }
}

// ---------------------------------------------------------------------------
// Occupancy
// ---------------------------------------------------------------------------

/// A validated site occupancy fraction: finite and `>= 0`.
///
/// No fixed per-value upper bound -- a single species can legitimately have
/// a low partial occupancy (vacancy modeling); the upper bound that matters
/// is the *sum* over all species at one [`PeriodicSite`], enforced by
/// [`PeriodicSite::validate`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Occupancy(f64);

impl Occupancy {
    /// Tolerance (dimensionless -- occupancy is already a unitless 0..1-ish
    /// fraction) that a site's total occupancy sum may exceed `1.0` by
    /// before being rejected. Sized to absorb floating-point summation
    /// error (e.g. three species at exactly 1/3 each can sum to
    /// `0.9999999999999999` or `1.0000000000000002` depending on summation
    /// order), not to permit genuine over-occupancy.
    pub const SUM_TOLERANCE: f64 = 1e-6;

    /// Construct a validated occupancy. Rejects non-finite and negative
    /// values.
    pub fn new(value: f64) -> Result<Self, CrystalError> {
        if !value.is_finite() {
            return Err(CrystalError::NonFiniteOccupancy);
        }
        if value < 0.0 {
            return Err(CrystalError::NegativeOccupancy { value });
        }
        Ok(Self(value))
    }

    /// The underlying fraction.
    #[inline]
    pub fn value(&self) -> f64 {
        self.0
    }
}

// ---------------------------------------------------------------------------
// Species / site
// ---------------------------------------------------------------------------

/// One (element, occupancy) contribution at a [`PeriodicSite`].
///
/// Multiple `SiteSpecies` at one site model substitutional/positional
/// disorder (e.g. a site 60% Fe / 40% Ni). `v0.1` stores this data
/// faithfully but does not implement any disorder-aware chemistry on top of
/// it.
#[derive(Debug, Clone, PartialEq)]
pub struct SiteSpecies {
    pub element: Element,
    pub occupancy: Occupancy,
}

impl SiteSpecies {
    /// Convenience constructor for a fully-occupied single-element site.
    pub fn full(element: Element) -> Self {
        Self {
            element,
            // 1.0 is always finite and non-negative -- infallible.
            occupancy: Occupancy::new(1.0).expect("1.0 is a valid occupancy"),
        }
    }
}

/// A periodic site: one or more species (for disorder), a fractional
/// position, and an optional human-readable label (e.g. `"Na1"` from a CIF
/// file).
#[derive(Debug, Clone, PartialEq)]
pub struct PeriodicSite {
    pub species: Vec<SiteSpecies>,
    pub fractional: FractionalCoord,
    pub label: Option<String>,
}

impl PeriodicSite {
    /// Construct and validate a site: rejects an empty species list, a
    /// non-finite fractional position, and a species-occupancy sum above
    /// `1.0 + `[`Occupancy::SUM_TOLERANCE`].
    pub fn new(
        species: Vec<SiteSpecies>,
        fractional: FractionalCoord,
        label: Option<String>,
    ) -> Result<Self, CrystalError> {
        let site = Self {
            species,
            fractional,
            label,
        };
        site.validate()?;
        Ok(site)
    }

    /// Re-run this site's own invariants (species non-empty, fractional
    /// finite, occupancy sum within tolerance). Called automatically by
    /// [`PeriodicSite::new`]; exposed for structures assembled from raw
    /// struct-literal fields (all of `PeriodicSite`'s fields are `pub`) so
    /// [`crate::structure::PeriodicStructure::validate`] can check them
    /// too.
    pub fn validate(&self) -> Result<(), CrystalError> {
        require_finite3(self.fractional.0, "fractional")?;
        if self.species.is_empty() {
            return Err(CrystalError::EmptySpeciesList);
        }
        let sum: f64 = self.species.iter().map(|s| s.occupancy.value()).sum();
        if sum > 1.0 + Occupancy::SUM_TOLERANCE {
            return Err(CrystalError::OccupancySumExceeded {
                sum,
                tolerance: Occupancy::SUM_TOLERANCE,
            });
        }
        Ok(())
    }
}

/// Hand-written rather than derived throughout this module: every public
/// type here (`FractionalCoord`, `CartesianCoord`, `Occupancy`,
/// `SiteSpecies`, `PeriodicSite`) enforces an invariant at construction
/// (finite components, non-negative occupancy, non-empty species,
/// occupancy-sum tolerance), and a `#[derive(Deserialize)]` would silently
/// skip all of them by assigning fields directly. `SiteSpecies` additionally
/// can't derive at all: `chematic_core::Element` has no serde support (see
/// `docs/rfcs/chematic_crystal_foundation.md`), so it round-trips through
/// the element's existing `symbol()`/`from_symbol()` API as a string
/// instead.
#[cfg(feature = "serde")]
mod serde_impl {
    use super::{CartesianCoord, FractionalCoord, Occupancy, PeriodicSite, SiteSpecies};
    use crate::error::CrystalError;
    use chematic_core::Element;
    use serde::de::Error as _;
    use serde::ser::SerializeStruct;
    use serde::{Deserialize, Deserializer, Serialize, Serializer};

    impl Serialize for FractionalCoord {
        fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
            self.0.serialize(serializer)
        }
    }

    impl<'de> Deserialize<'de> for FractionalCoord {
        fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
            let value = <[f64; 3]>::deserialize(deserializer)?;
            if value.iter().all(|c| c.is_finite()) {
                Ok(FractionalCoord(value))
            } else {
                Err(D::Error::custom(
                    "fractional coordinate components must be finite",
                ))
            }
        }
    }

    impl Serialize for CartesianCoord {
        fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
            self.0.serialize(serializer)
        }
    }

    impl<'de> Deserialize<'de> for CartesianCoord {
        fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
            let value = <[f64; 3]>::deserialize(deserializer)?;
            if value.iter().all(|c| c.is_finite()) {
                Ok(CartesianCoord(value))
            } else {
                Err(D::Error::custom(
                    "Cartesian coordinate components must be finite",
                ))
            }
        }
    }

    impl Serialize for Occupancy {
        fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
            self.0.serialize(serializer)
        }
    }

    impl<'de> Deserialize<'de> for Occupancy {
        fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
            let value = f64::deserialize(deserializer)?;
            Occupancy::new(value).map_err(|e: CrystalError| D::Error::custom(e.to_string()))
        }
    }

    impl Serialize for SiteSpecies {
        fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
            let mut state = serializer.serialize_struct("SiteSpecies", 2)?;
            state.serialize_field("element", self.element.symbol())?;
            state.serialize_field("occupancy", &self.occupancy)?;
            state.end()
        }
    }

    #[derive(Deserialize)]
    struct SiteSpeciesData {
        element: String,
        occupancy: Occupancy,
    }

    impl<'de> Deserialize<'de> for SiteSpecies {
        fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
            let data = SiteSpeciesData::deserialize(deserializer)?;
            let element = Element::from_symbol(&data.element).ok_or_else(|| {
                D::Error::custom(format!("unknown element symbol {:?}", data.element))
            })?;
            Ok(SiteSpecies {
                element,
                occupancy: data.occupancy,
            })
        }
    }

    impl Serialize for PeriodicSite {
        fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
            let mut state = serializer.serialize_struct("PeriodicSite", 3)?;
            state.serialize_field("species", &self.species)?;
            state.serialize_field("fractional", &self.fractional)?;
            state.serialize_field("label", &self.label)?;
            state.end()
        }
    }

    #[derive(Deserialize)]
    struct PeriodicSiteData {
        species: Vec<SiteSpecies>,
        fractional: FractionalCoord,
        label: Option<String>,
    }

    impl<'de> Deserialize<'de> for PeriodicSite {
        fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
            let data = PeriodicSiteData::deserialize(deserializer)?;
            PeriodicSite::new(data.species, data.fractional, data.label)
                .map_err(|e: CrystalError| D::Error::custom(e.to_string()))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fractional_wrapped_reduces_into_unit_range() {
        let f = FractionalCoord::new([1.25, -0.25, 2.75]);
        let w = f.wrapped();
        assert!((w.0[0] - 0.25).abs() < 1e-12);
        assert!((w.0[1] - 0.75).abs() < 1e-12);
        assert!((w.0[2] - 0.75).abs() < 1e-12);
        for c in w.0 {
            assert!((0.0..1.0).contains(&c));
        }
    }

    #[test]
    fn fractional_wrapped_one_maps_to_zero() {
        let w = FractionalCoord::new([1.0, 1.0, 1.0]).wrapped();
        assert_eq!(w.0, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn fractional_translated_is_integer_shift() {
        let f = FractionalCoord::new([0.1, 0.2, 0.3]);
        let t = f.translated([1, -2, 0]);
        assert!((t.0[0] - 1.1).abs() < 1e-12);
        assert!((t.0[1] - (-1.8)).abs() < 1e-12);
        assert!((t.0[2] - 0.3).abs() < 1e-12);
    }

    #[test]
    fn occupancy_rejects_negative_and_non_finite() {
        assert!(Occupancy::new(0.5).is_ok());
        assert!(Occupancy::new(0.0).is_ok());
        assert_eq!(
            Occupancy::new(-0.1),
            Err(CrystalError::NegativeOccupancy { value: -0.1 })
        );
        assert_eq!(
            Occupancy::new(f64::NAN),
            Err(CrystalError::NonFiniteOccupancy)
        );
        assert_eq!(
            Occupancy::new(f64::INFINITY),
            Err(CrystalError::NonFiniteOccupancy)
        );
    }

    #[test]
    fn site_rejects_empty_species() {
        let err =
            PeriodicSite::new(vec![], FractionalCoord::new([0.0, 0.0, 0.0]), None).unwrap_err();
        assert_eq!(err, CrystalError::EmptySpeciesList);
    }

    #[test]
    fn site_rejects_non_finite_fractional() {
        let species = vec![SiteSpecies::full(Element::NA)];
        let err = PeriodicSite::new(species, FractionalCoord::new([f64::NAN, 0.0, 0.0]), None)
            .unwrap_err();
        assert_eq!(
            err,
            CrystalError::NonFinite {
                field: "fractional"
            }
        );
    }

    #[test]
    fn site_accepts_occupancy_sum_at_one() {
        let species = vec![
            SiteSpecies {
                element: Element::FE,
                occupancy: Occupancy::new(0.6).unwrap(),
            },
            SiteSpecies {
                element: Element::NI,
                occupancy: Occupancy::new(0.4).unwrap(),
            },
        ];
        assert!(PeriodicSite::new(species, FractionalCoord::new([0.0, 0.0, 0.0]), None).is_ok());
    }

    #[test]
    fn site_accepts_occupancy_sum_below_one_as_vacancy() {
        let species = vec![SiteSpecies {
            element: Element::FE,
            occupancy: Occupancy::new(0.5).unwrap(),
        }];
        assert!(PeriodicSite::new(species, FractionalCoord::new([0.0, 0.0, 0.0]), None).is_ok());
    }

    #[test]
    fn site_rejects_occupancy_sum_above_tolerance() {
        let species = vec![
            SiteSpecies {
                element: Element::FE,
                occupancy: Occupancy::new(0.7).unwrap(),
            },
            SiteSpecies {
                element: Element::NI,
                occupancy: Occupancy::new(0.5).unwrap(),
            },
        ];
        let err =
            PeriodicSite::new(species, FractionalCoord::new([0.0, 0.0, 0.0]), None).unwrap_err();
        assert!(matches!(err, CrystalError::OccupancySumExceeded { .. }));
    }

    #[test]
    fn site_accepts_occupancy_sum_within_float_tolerance() {
        // Three species at 1/3 each: floating-point summation lands just
        // above or below 1.0 depending on order, must not be rejected.
        let third = Occupancy::new(1.0 / 3.0).unwrap();
        let species = vec![
            SiteSpecies {
                element: Element::FE,
                occupancy: third,
            },
            SiteSpecies {
                element: Element::NI,
                occupancy: third,
            },
            SiteSpecies {
                element: Element::CU,
                occupancy: third,
            },
        ];
        assert!(PeriodicSite::new(species, FractionalCoord::new([0.0, 0.0, 0.0]), None).is_ok());
    }
}
