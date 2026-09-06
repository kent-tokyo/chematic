//! [`PeriodicStructure`]: a [`Lattice`] plus its [`PeriodicSite`]s.

use crate::error::CrystalError;
use crate::lattice::Lattice;
use crate::neighbor::{self, PeriodicNeighbor};
use crate::site::{CartesianCoord, FractionalCoord, PeriodicSite};
use crate::supercell;
use chematic_core::Element;

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

/// Deterministic, occupancy-weighted elemental composition of the stored
/// periodic cell.
///
/// Amounts are per *stored* cell, not normalized to one occupied site. A
/// species with zero occupancy remains present in the summary with amount
/// `0.0`, because it is explicit source data. Distinct sites are aggregated
/// only by element; sites at the same coordinate are still counted
/// independently before aggregation.
#[derive(Debug, Clone, PartialEq)]
pub struct CompositionSummary {
    amounts: std::collections::BTreeMap<Element, f64>,
    site_count: usize,
}

impl CompositionSummary {
    /// Occupancy-weighted amount for an element, or `0.0` if absent.
    pub fn amount(&self, element: Element) -> f64 {
        self.amounts.get(&element).copied().unwrap_or(0.0)
    }

    /// Deterministic iteration in atomic-number order.
    pub fn iter(&self) -> impl Iterator<Item = (Element, f64)> + '_ {
        self.amounts
            .iter()
            .map(|(&element, &amount)| (element, amount))
    }

    /// Number of sites in the stored cell, before elemental aggregation.
    pub fn site_count(&self) -> usize {
        self.site_count
    }

    /// Number of distinct element entries, including zero-occupancy species.
    pub fn element_count(&self) -> usize {
        self.amounts.len()
    }
}

impl PeriodicStructure {
    /// Version byte for [`Self::identity_bytes`].
    pub const IDENTITY_SERIALIZATION_VERSION: u8 = 1;

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

    /// Serialize the exact validated structure into a deterministic binary
    /// identity representation.
    ///
    /// The representation is versioned and intentionally preserves lattice
    /// matrix bits, site order, species order, occupancy bits, and labels.
    /// It is therefore suitable for exact cache keys and content-addressed
    /// storage, but is not a format-stability promise beyond the version byte.
    /// Use [`Self::IDENTITY_SERIALIZATION_VERSION`] when persisting it.
    pub fn identity_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(b"SCHEMATIC-PERIODIC\0");
        bytes.push(Self::IDENTITY_SERIALIZATION_VERSION);
        for row in self.lattice.matrix() {
            for value in row {
                bytes.extend_from_slice(&value.to_bits().to_le_bytes());
            }
        }
        push_u64(&mut bytes, self.sites.len());
        for site in &self.sites {
            push_u64(&mut bytes, site.species.len());
            for species in &site.species {
                bytes.push(species.element.atomic_number());
                bytes.extend_from_slice(&species.occupancy.value().to_bits().to_le_bytes());
            }
            for value in site.fractional.0 {
                bytes.extend_from_slice(&value.to_bits().to_le_bytes());
            }
            match &site.label {
                Some(label) => {
                    bytes.push(1);
                    push_bytes(&mut bytes, label.as_bytes());
                }
                None => bytes.push(0),
            }
        }
        bytes
    }

    /// Compute a SHA-256 digest of [`Self::identity_bytes`].
    ///
    /// The digest is an exact stored-representation key, not a symmetry-aware
    /// crystal fingerprint or a material-similarity score. The version tag in
    /// [`Self::identity_bytes`] therefore participates in the digest and
    /// prevents incompatible encodings from silently sharing keys.
    pub fn identity_digest(&self) -> [u8; 32] {
        sha256(&self.identity_bytes())
    }

    /// Compute the occupancy-weighted elemental composition of this stored
    /// cell. Calling this on a supercell reports the supercell amount (for
    /// example, a 2x2x2 cell reports eight times the amount of a one-cell
    /// structure); it never silently reduces back to the original cell.
    pub fn composition(&self) -> CompositionSummary {
        let mut amounts = std::collections::BTreeMap::new();
        for site in &self.sites {
            for species in &site.species {
                *amounts.entry(species.element).or_insert(0.0) += species.occupancy.value();
            }
        }
        CompositionSummary {
            amounts,
            site_count: self.sites.len(),
        }
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

fn push_u64(bytes: &mut Vec<u8>, value: usize) {
    bytes.extend_from_slice(&(value as u64).to_le_bytes());
}

fn push_bytes(bytes: &mut Vec<u8>, value: &[u8]) {
    push_u64(bytes, value.len());
    bytes.extend_from_slice(value);
}

fn sha256(input: &[u8]) -> [u8; 32] {
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
        0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
        0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
        0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
        0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
        0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
        0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
        0xc67178f2,
    ];
    let mut h: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
        0x5be0cd19,
    ];
    let bit_len = (input.len() as u64).wrapping_mul(8);
    let mut padded = input.to_vec();
    padded.push(0x80);
    while padded.len() % 64 != 56 {
        padded.push(0);
    }
    padded.extend_from_slice(&bit_len.to_be_bytes());

    for chunk in padded.chunks_exact(64) {
        let mut w = [0u32; 64];
        for (i, word) in chunk.chunks_exact(4).take(16).enumerate() {
            w[i] = u32::from_be_bytes([word[0], word[1], word[2], word[3]]);
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16]
                .wrapping_add(s0)
                .wrapping_add(w[i - 7])
                .wrapping_add(s1);
        }
        let (mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut hh): (
            u32,
            u32,
            u32,
            u32,
            u32,
            u32,
            u32,
            u32,
        ) = (h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7]);
        for i in 0..64 {
            let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let ch = (e & f) ^ ((!e) & g);
            let t1 = hh
                .wrapping_add(s1)
                .wrapping_add(ch)
                .wrapping_add(K[i])
                .wrapping_add(w[i]);
            let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let maj = (a & b) ^ (a & c) ^ (b & c);
            let t2 = s0.wrapping_add(maj);
            (hh, g, f, e, d, c, b, a) = (g, f, e, d.wrapping_add(t1), c, b, a, t1.wrapping_add(t2));
        }
        for (state, value) in h.iter_mut().zip([a, b, c, d, e, f, g, hh]) {
            *state = (*state).wrapping_add(value);
        }
    }

    let mut digest = [0u8; 32];
    for (i, word) in h.iter().enumerate() {
        digest[i * 4..i * 4 + 4].copy_from_slice(&word.to_be_bytes());
    }
    digest
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
    fn identity_serialization_is_deterministic_and_exact() {
        let structure = cubic_two_site();
        let first = structure.identity_bytes();
        let digest = structure.identity_digest();
        assert_eq!(first, structure.identity_bytes());
        assert_eq!(digest, structure.identity_digest());
        assert!(first.starts_with(b"SCHEMATIC-PERIODIC\0"));
        assert_eq!(
            first[b"SCHEMATIC-PERIODIC\0".len()],
            PeriodicStructure::IDENTITY_SERIALIZATION_VERSION
        );

        let mut changed = structure.sites().to_vec();
        changed[0].label = Some("Na2".to_string());
        let changed = PeriodicStructure::new(structure.lattice().clone(), changed).unwrap();
        assert_ne!(first, changed.identity_bytes());
        assert_ne!(digest, changed.identity_digest());
    }

    #[test]
    fn identity_digest_uses_standard_sha256() {
        assert_eq!(
            hex(&sha256(b"abc")),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
    }

    fn hex(bytes: &[u8]) -> String {
        bytes.iter().map(|b| format!("{b:02x}")).collect()
    }

    #[test]
    fn composition_is_occupancy_weighted_deterministic_and_keeps_zero_species() {
        let lattice = Lattice::cubic(4.0).unwrap();
        let sites = vec![
            PeriodicSite::new(
                vec![
                    SiteSpecies {
                        element: Element::FE,
                        occupancy: Occupancy::new(0.6).unwrap(),
                    },
                    SiteSpecies {
                        element: Element::NI,
                        occupancy: Occupancy::new(0.4).unwrap(),
                    },
                ],
                FractionalCoord::new([0.0, 0.0, 0.0]),
                Some("mixed".into()),
            )
            .unwrap(),
            PeriodicSite::new(
                vec![SiteSpecies {
                    element: Element::O,
                    occupancy: Occupancy::new(0.0).unwrap(),
                }],
                FractionalCoord::new([0.0, 0.0, 0.0]),
                Some("vacant".into()),
            )
            .unwrap(),
        ];
        let structure = PeriodicStructure::new(lattice, sites).unwrap();
        let composition = structure.composition();
        assert_eq!(composition.site_count(), 2);
        assert_eq!(composition.element_count(), 3);
        assert!((composition.amount(Element::FE) - 0.6).abs() < 1e-12);
        assert!((composition.amount(Element::NI) - 0.4).abs() < 1e-12);
        assert_eq!(composition.amount(Element::O), 0.0);
        assert_eq!(
            composition
                .iter()
                .map(|(element, _)| element)
                .collect::<Vec<_>>(),
            vec![Element::O, Element::FE, Element::NI]
        );
    }

    #[test]
    fn composition_scales_with_explicit_supercell() {
        let structure = cubic_two_site();
        let composition = structure.make_supercell([2, 2, 2]).unwrap().composition();
        assert_eq!(composition.site_count(), 16);
        assert_eq!(composition.amount(Element::NA), 8.0);
        assert_eq!(composition.amount(Element::CL), 8.0);
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
