//! [`PeriodicStructure`]: a [`Lattice`] plus its [`PeriodicSite`]s.

use crate::error::CrystalError;
use crate::lattice::Lattice;
use crate::neighbor::{self, PeriodicNeighbor};
use crate::site::{CartesianCoord, FractionalCoord, PeriodicSite};
use crate::supercell;

/// A periodic structure: a validated [`Lattice`] and an ordered list of
/// [`PeriodicSite`]s.
///
/// Immutable by design in `v0.1` -- operations that would change geometry
/// (wrapping, supercell expansion) return a new `PeriodicStructure` rather
/// than mutating in place, so an already-validated instance can never be
/// pushed into an invalid state after construction.
#[derive(Debug, Clone, PartialEq)]
pub struct PeriodicStructure {
    lattice: Lattice,
    sites: Vec<PeriodicSite>,
}

impl PeriodicStructure {
    /// Construct and validate a structure: every site must independently
    /// pass [`PeriodicSite::validate`].
    ///
    /// # Examples
    ///
    /// CsCl-type structure (one cation + one anion per cubic cell, anion
    /// at the body center -- e.g. CsCl, AlNi, beta-brass; illustrated here
    /// with Na/Cl for familiarity, not a claim about real NaCl, which is
    /// rock-salt-type and needs 8 sites in its conventional cubic cell):
    ///
    /// ```
    /// use chematic_core::Element;
    /// use chematic_crystal::{FractionalCoord, Lattice, PeriodicSite, PeriodicStructure, SiteSpecies};
    ///
    /// let lattice = Lattice::cubic(5.64)?;
    /// let sites = vec![
    ///     PeriodicSite::new(
    ///         vec![SiteSpecies::full(Element::NA)],
    ///         FractionalCoord::new([0.0, 0.0, 0.0]),
    ///         Some("Na1".to_string()),
    ///     )?,
    ///     PeriodicSite::new(
    ///         vec![SiteSpecies::full(Element::CL)],
    ///         FractionalCoord::new([0.5, 0.5, 0.5]),
    ///         Some("Cl1".to_string()),
    ///     )?,
    /// ];
    /// let structure = PeriodicStructure::new(lattice, sites)?;
    /// assert_eq!(structure.site_count(), 2);
    /// # Ok::<(), chematic_crystal::CrystalError>(())
    /// ```
    pub fn new(lattice: Lattice, sites: Vec<PeriodicSite>) -> Result<Self, CrystalError> {
        let structure = Self { lattice, sites };
        structure.validate()?;
        Ok(structure)
    }

    /// The structure's lattice.
    #[inline]
    pub fn lattice(&self) -> &Lattice {
        &self.lattice
    }

    /// The structure's sites, in construction order.
    #[inline]
    pub fn sites(&self) -> &[PeriodicSite] {
        &self.sites
    }

    /// Number of sites.
    #[inline]
    pub fn site_count(&self) -> usize {
        self.sites.len()
    }

    /// Fractional position of every site, in site order.
    pub fn fractional_positions(&self) -> Vec<FractionalCoord> {
        self.sites.iter().map(|s| s.fractional).collect()
    }

    /// Cartesian position of every site, in site order (`lattice.frac_to_cart`
    /// applied to each site's fractional position).
    pub fn cartesian_positions(&self) -> Vec<CartesianCoord> {
        self.sites
            .iter()
            .map(|s| self.lattice.frac_to_cart(s.fractional))
            .collect()
    }

    /// Re-check every site's own invariants. Called automatically by
    /// [`Self::new`]; exposed since `sites()` is a read-only view but
    /// structures can also be produced by other crate-internal
    /// constructors ([`Self::wrapped`], `make_supercell`) that reuse this
    /// check rather than re-deriving it.
    pub fn validate(&self) -> Result<(), CrystalError> {
        for (index, site) in self.sites.iter().enumerate() {
            site.validate()
                .map_err(|source| CrystalError::InvalidSite {
                    index,
                    source: Box::new(source),
                })?;
        }
        Ok(())
    }

    /// A new structure with every site's fractional position reduced into
    /// `[0, 1)` via [`FractionalCoord::wrapped`]. Species, occupancy, and
    /// labels are preserved; `self` is unmodified.
    pub fn wrapped(&self) -> Self {
        let sites = self
            .sites
            .iter()
            .map(|s| PeriodicSite {
                species: s.species.clone(),
                fractional: s.fractional.wrapped(),
                label: s.label.clone(),
            })
            .collect();
        // Wrapping a finite fractional coordinate (guaranteed by `self`
        // already having passed validation) can't un-finite it, and leaves
        // species/occupancy untouched -- constructing directly (rather than
        // re-running `Self::new`'s Result-returning validation) is safe.
        Self {
            lattice: self.lattice.clone(),
            sites,
        }
    }

    /// Every periodic neighbor pair within `cutoff` Angstrom (inclusive).
    /// See [`neighbor::neighbors_within`] for the full contract
    /// (self-image handling, ordering, error conditions).
    ///
    /// # Examples
    ///
    /// ```
    /// use chematic_core::Element;
    /// use chematic_crystal::{FractionalCoord, Lattice, PeriodicSite, PeriodicStructure, SiteSpecies};
    ///
    /// let lattice = Lattice::cubic(3.0)?;
    /// let structure = PeriodicStructure::new(
    ///     lattice,
    ///     vec![PeriodicSite::new(
    ///         vec![SiteSpecies::full(Element::AR)],
    ///         FractionalCoord::new([0.0, 0.0, 0.0]),
    ///         None,
    ///     )?],
    /// )?;
    /// // 6 face-adjacent periodic self-images at exactly 3.0 Angstrom.
    /// let neighbors = structure.neighbors_within(3.0)?;
    /// assert_eq!(neighbors.len(), 6);
    /// # Ok::<(), chematic_crystal::CrystalError>(())
    /// ```
    pub fn neighbors_within(&self, cutoff: f64) -> Result<Vec<PeriodicNeighbor>, CrystalError> {
        neighbor::neighbors_within(self, cutoff)
    }

    /// Build a diagonal `[nx, ny, nz]` supercell. See
    /// [`supercell::make_supercell`] for the full contract (site ordering,
    /// error conditions).
    ///
    /// # Examples
    ///
    /// ```
    /// use chematic_core::Element;
    /// use chematic_crystal::{FractionalCoord, Lattice, PeriodicSite, PeriodicStructure, SiteSpecies};
    ///
    /// let lattice = Lattice::cubic(4.0)?;
    /// let structure = PeriodicStructure::new(
    ///     lattice,
    ///     vec![PeriodicSite::new(
    ///         vec![SiteSpecies::full(Element::C)],
    ///         FractionalCoord::new([0.0, 0.0, 0.0]),
    ///         None,
    ///     )?],
    /// )?;
    /// let supercell = structure.make_supercell([2, 2, 2])?;
    /// assert_eq!(supercell.site_count(), 8);
    /// assert!((supercell.lattice().volume() - structure.lattice().volume() * 8.0).abs() < 1e-9);
    /// # Ok::<(), chematic_crystal::CrystalError>(())
    /// ```
    pub fn make_supercell(&self, mult: [u32; 3]) -> Result<Self, CrystalError> {
        supercell::make_supercell(self, mult)
    }
}

#[cfg(feature = "serde")]
mod serde_impl {
    use super::PeriodicStructure;
    use crate::error::CrystalError;
    use crate::lattice::Lattice;
    use crate::site::PeriodicSite;
    use serde::de::Error as _;
    use serde::ser::SerializeStruct;
    use serde::{Deserialize, Deserializer, Serialize, Serializer};

    impl Serialize for PeriodicStructure {
        fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
            let mut state = serializer.serialize_struct("PeriodicStructure", 2)?;
            state.serialize_field("lattice", &self.lattice)?;
            state.serialize_field("sites", &self.sites)?;
            state.end()
        }
    }

    /// Deserializes into raw `(lattice, sites)` fields, then re-validates
    /// through [`PeriodicStructure::new`] -- a deserialized structure can
    /// never skip the invariants a normally-constructed one goes through.
    #[derive(Deserialize)]
    struct PeriodicStructureData {
        lattice: Lattice,
        sites: Vec<PeriodicSite>,
    }

    impl<'de> Deserialize<'de> for PeriodicStructure {
        fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
            let data = PeriodicStructureData::deserialize(deserializer)?;
            PeriodicStructure::new(data.lattice, data.sites)
                .map_err(|e: CrystalError| D::Error::custom(e.to_string()))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::site::{Occupancy, SiteSpecies};
    use chematic_core::Element;

    fn cubic_two_site() -> PeriodicStructure {
        let lattice = Lattice::cubic(4.0).unwrap();
        let sites = vec![
            PeriodicSite::new(
                vec![SiteSpecies::full(Element::NA)],
                FractionalCoord::new([0.0, 0.0, 0.0]),
                Some("Na1".to_string()),
            )
            .unwrap(),
            PeriodicSite::new(
                vec![SiteSpecies::full(Element::CL)],
                FractionalCoord::new([0.5, 0.5, 0.5]),
                Some("Cl1".to_string()),
            )
            .unwrap(),
        ];
        PeriodicStructure::new(lattice, sites).unwrap()
    }

    #[test]
    fn new_accepts_valid_sites() {
        let s = cubic_two_site();
        assert_eq!(s.site_count(), 2);
        assert_eq!(s.sites().len(), 2);
    }

    #[test]
    fn new_rejects_invalid_site_with_index() {
        let lattice = Lattice::cubic(4.0).unwrap();
        let bad_site = PeriodicSite {
            species: vec![],
            fractional: FractionalCoord::new([0.0, 0.0, 0.0]),
            label: None,
        };
        let err = PeriodicStructure::new(lattice, vec![bad_site]).unwrap_err();
        match err {
            CrystalError::InvalidSite { index, source } => {
                assert_eq!(index, 0);
                assert_eq!(*source, CrystalError::EmptySpeciesList);
            }
            other => panic!("expected InvalidSite, got {other:?}"),
        }
    }

    #[test]
    fn fractional_and_cartesian_positions_match_site_order() {
        let s = cubic_two_site();
        let frac = s.fractional_positions();
        let cart = s.cartesian_positions();
        assert_eq!(frac.len(), 2);
        assert_eq!(cart.len(), 2);
        assert_eq!(frac[0].0, [0.0, 0.0, 0.0]);
        assert_eq!(cart[0].0, [0.0, 0.0, 0.0]);
        // second site at (0.5,0.5,0.5) in a cubic(4.0) cell -> (2,2,2) Cartesian.
        for c in cart[1].0 {
            assert!((c - 2.0).abs() < 1e-9);
        }
    }

    #[test]
    fn wrapped_reduces_out_of_cell_fractional_and_preserves_species() {
        let lattice = Lattice::cubic(4.0).unwrap();
        let site = PeriodicSite::new(
            vec![SiteSpecies {
                element: Element::FE,
                occupancy: Occupancy::new(0.6).unwrap(),
            }],
            FractionalCoord::new([1.25, -0.5, 3.0]),
            Some("Fe1".to_string()),
        )
        .unwrap();
        let s = PeriodicStructure::new(lattice, vec![site]).unwrap();
        let w = s.wrapped();
        let wf = w.sites()[0].fractional.0;
        assert!((wf[0] - 0.25).abs() < 1e-12);
        assert!((wf[1] - 0.5).abs() < 1e-12);
        assert!((wf[2] - 0.0).abs() < 1e-12);
        assert_eq!(w.sites()[0].species, s.sites()[0].species);
        assert_eq!(w.sites()[0].label, s.sites()[0].label);
        // original unchanged
        assert_eq!(s.sites()[0].fractional.0, [1.25, -0.5, 3.0]);
    }

    #[test]
    fn site_order_is_preserved_not_reordered() {
        let s = cubic_two_site();
        assert_eq!(s.sites()[0].label.as_deref(), Some("Na1"));
        assert_eq!(s.sites()[1].label.as_deref(), Some("Cl1"));
    }
}
