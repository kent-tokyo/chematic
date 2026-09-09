//! Experimental, deterministic geometry fingerprint for caller-supplied structures.
//! This is not symmetry detection, polymorph identification, or a retrieval score.

use crate::PeriodicStructure;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FingerprintLimits {
    pub max_sites: usize,
    pub max_species: usize,
}

impl Default for FingerprintLimits {
    fn default() -> Self {
        Self {
            max_sites: 100_000,
            max_species: 200_000,
        }
    }
}

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FingerprintError {
    MissingStructure,
    ResourceLimit,
    NonFiniteGeometry,
}

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeometryFingerprint {
    pub schema_version: u32,
    pub fingerprint: String,
    pub source: String,
    pub site_count: usize,
    pub species_count: usize,
    pub semantics: &'static str,
}

/// Hash the exact stored lattice/site order/species/occupancy/fractional coordinates.
/// Labels are intentionally excluded: they are source annotations, not geometry.
pub fn geometry_fingerprint(
    structure: Option<&PeriodicStructure>,
    source: impl Into<String>,
    limits: FingerprintLimits,
) -> Result<GeometryFingerprint, FingerprintError> {
    let Some(structure) = structure else {
        return Err(FingerprintError::MissingStructure);
    };
    if structure.site_count() > limits.max_sites {
        return Err(FingerprintError::ResourceLimit);
    }
    let species_count: usize = structure.sites().iter().map(|s| s.species.len()).sum();
    if species_count > limits.max_species {
        return Err(FingerprintError::ResourceLimit);
    }
    let mut hash = 0xcbf29ce484222325u64;
    let mut add = |bytes: &[u8]| {
        for b in bytes {
            hash ^= u64::from(*b);
            hash = hash.wrapping_mul(0x100000001b3);
        }
    };
    for row in structure.lattice().matrix() {
        for value in row {
            if !value.is_finite() {
                return Err(FingerprintError::NonFiniteGeometry);
            }
            add(&value.to_bits().to_le_bytes());
        }
    }
    for site in structure.sites() {
        for value in site.fractional.0 {
            if !value.is_finite() {
                return Err(FingerprintError::NonFiniteGeometry);
            }
            add(&value.to_bits().to_le_bytes());
        }
        for species in &site.species {
            add(&species.element.atomic_number().to_le_bytes());
            add(&species.occupancy.value().to_bits().to_le_bytes());
        }
    }
    Ok(GeometryFingerprint {
        schema_version: 1,
        fingerprint: format!("fnv1a64-{hash:016x}"),
        source: source.into(),
        site_count: structure.site_count(),
        species_count,
        semantics: "exact stored geometry descriptor; not symmetry or polymorph identity",
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{FractionalCoord, Lattice, PeriodicSite, SiteSpecies};
    use chematic_core::Element;
    fn structure(x: f64) -> PeriodicStructure {
        PeriodicStructure::new(
            Lattice::cubic(4.0).unwrap(),
            vec![
                PeriodicSite::new(
                    vec![SiteSpecies::full(Element::SI)],
                    FractionalCoord::new([x, 0.0, 0.0]),
                    None,
                )
                .unwrap(),
            ],
        )
        .unwrap()
    }
    #[test]
    fn deterministic_and_geometry_sensitive() {
        let a = geometry_fingerprint(
            Some(&structure(0.0)),
            "fixture",
            FingerprintLimits::default(),
        )
        .unwrap();
        let b = geometry_fingerprint(
            Some(&structure(0.0)),
            "fixture",
            FingerprintLimits::default(),
        )
        .unwrap();
        let c = geometry_fingerprint(
            Some(&structure(0.1)),
            "fixture",
            FingerprintLimits::default(),
        )
        .unwrap();
        assert_eq!(a, b);
        assert_ne!(a.fingerprint, c.fingerprint);
    }
    #[test]
    fn missing_and_limits_fail_closed() {
        assert_eq!(
            geometry_fingerprint(None, "x", FingerprintLimits::default()),
            Err(FingerprintError::MissingStructure)
        );
        let s = structure(0.0);
        assert_eq!(
            geometry_fingerprint(
                Some(&s),
                "x",
                FingerprintLimits {
                    max_sites: 0,
                    max_species: 1
                }
            ),
            Err(FingerprintError::ResourceLimit)
        );
    }
}
