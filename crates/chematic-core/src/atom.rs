//! Atom type: a single atom in a molecule.

use crate::element::Element;

/// A `@SP1`/`@SP2`/`@SP3` square-planar stereo tag (OpenSMILES's extended
/// chirality-class syntax, e.g. Pt(II)/Pd(II) complexes like cisplatin).
///
/// Each variant names which pair of the 4 explicit neighbor positions
/// (0-indexed, in the order recorded by [`crate::Molecule::stereo_neighbor_order`])
/// sit *trans* (~180°) to each other:
///
/// - `SP1`: positions (0,2) trans, (1,3) trans
/// - `SP2`: positions (0,1) trans, (2,3) trans
/// - `SP3`: positions (0,3) trans, (1,2) trans
///
/// Oracle-verified against RDKit 2026.03.3 (3D embedding bond angles, cross-checked
/// against RDKit's own documented cisplatin/transplatin example) — see
/// `docs/rfcs/square_planar_stereo_rfc.md`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SquarePlanarPermutation {
    SP1,
    SP2,
    SP3,
}

impl SquarePlanarPermutation {
    /// The two trans-pairs this permutation implies, as 0-indexed neighbor positions.
    pub fn trans_pairs(self) -> [(u8, u8); 2] {
        match self {
            Self::SP1 => [(0, 2), (1, 3)],
            Self::SP2 => [(0, 1), (2, 3)],
            Self::SP3 => [(0, 3), (1, 2)],
        }
    }
}

/// Chirality as specified in OpenSMILES: tetrahedral (`@`/`@@`) or, since this
/// crate also models coordination complexes, square-planar (`@SP1`/`@SP2`/`@SP3`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Chirality {
    /// No chirality specified.
    #[default]
    None,
    /// `@` — counterclockwise (looking from the first neighbor).
    CounterClockwise,
    /// `@@` — clockwise.
    Clockwise,
    /// `@SP1`/`@SP2`/`@SP3` — square-planar (4-coordinate) stereo.
    SquarePlanar(SquarePlanarPermutation),
}

impl Chirality {
    /// `true` only for [`Self::CounterClockwise`]/[`Self::Clockwise`] — the classic
    /// tetrahedral-parity forms every CIP/ECFP/dedup consumer written before
    /// square-planar existed assumes. Consumers that mean "is this a real
    /// tetrahedral stereocenter" must check this, not `!= Self::None`, now that a
    /// second non-tetrahedral kind of "not None" chirality exists.
    pub fn is_tetrahedral(&self) -> bool {
        matches!(self, Self::CounterClockwise | Self::Clockwise)
    }
}

/// Assigned CIP (Cahn–Ingold–Prelog) stereodescriptor.
///
/// Stored on [`Atom`] after running [`chematic_chem::assign_cip`] or
/// [`chematic_chem::cip::assign_cip`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CipCode {
    /// Tetrahedral center with *rectus* (right-handed) configuration.
    R,
    /// Tetrahedral center with *sinister* (left-handed) configuration.
    S,
    /// Double-bond *entgegen* (opposite, trans) geometry.
    E,
    /// Double-bond *zusammen* (together, cis) geometry.
    Z,
    /// Pseudoasymmetric center, *rectus*-like (Rule 5, lowercase `r`). Emitted only by
    /// `chematic_cip::assign_cip_accurate_experimental`'s Rule 5 pass; the default
    /// `chematic_chem::assign_cip` never produces this variant.
    LowerR,
    /// Pseudoasymmetric center, *sinister*-like (Rule 5, lowercase `s`). See [`Self::LowerR`].
    LowerS,
}

/// A single atom in a molecular graph.
///
/// - `isotope`: mass number (e.g. 13 for ¹³C). `None` = natural isotope abundance.
/// - `charge`: formal charge.
/// - `hydrogen_count`: explicit H count from a bracket atom `[...]`.
///   `None` for organic-subset atoms whose H count is inferred from valence.
/// - `aromatic`: set when the atom is written as a lowercase letter (c, n, …)
///   or connected via `:` bonds.
/// - `wildcard`: `true` for the SMILES `*` atom (any element, query context).
/// - `atom_map`: atom-mapping number used in reaction SMILES.
/// - `cip_code`: CIP stereodescriptor (R/S/E/Z). Populated by
///   `chematic_chem::assign_cip`; `None` until then.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Atom {
    pub element: Element,
    pub isotope: Option<u16>,
    pub charge: i8,
    /// Explicit H count (bracket atoms only). `None` for organic-subset atoms.
    pub hydrogen_count: Option<u8>,
    pub aromatic: bool,
    pub chirality: Chirality,
    /// True for the wildcard atom `*` or `[*]`.
    pub wildcard: bool,
    pub atom_map: Option<u16>,
    /// CIP stereodescriptor assigned by `chematic_chem::assign_cip`.
    /// `None` until explicitly computed.
    pub cip_code: Option<CipCode>,
}

impl Atom {
    /// Create a plain, neutral, non-aromatic atom.
    pub fn new(element: Element) -> Self {
        Self {
            element,
            isotope: None,
            charge: 0,
            hydrogen_count: None,
            aromatic: false,
            chirality: Chirality::None,
            wildcard: false,
            atom_map: None,
            cip_code: None,
        }
    }

    /// Organic-subset atom (charge=0, non-aromatic, implicit H from valence).
    pub fn organic(element: Element) -> Self {
        Self::new(element)
    }

    /// Aromatic organic atom (lowercase SMILES notation).
    pub fn aromatic(element: Element) -> Self {
        Self {
            aromatic: true,
            ..Self::new(element)
        }
    }

    /// Bracket atom with explicit properties.
    pub fn bracket(
        element: Element,
        isotope: Option<u16>,
        chirality: Chirality,
        hydrogen_count: u8,
        charge: i8,
        atom_map: Option<u16>,
    ) -> Self {
        Self {
            element,
            isotope,
            charge,
            hydrogen_count: Some(hydrogen_count),
            aromatic: false,
            chirality,
            wildcard: false,
            atom_map,
            cip_code: None,
        }
    }

    /// Wildcard atom `*` / `[*]` (matches any element in query contexts).
    pub fn wildcard() -> Self {
        Self {
            // Element is a placeholder; callers should check `wildcard` first.
            element: Element::C,
            wildcard: true,
            hydrogen_count: Some(0),
            ..Self::new(Element::C)
        }
    }

    /// Return the explicit H count for bracket atoms; `None` for organic-subset atoms.
    pub fn explicit_hcount(&self) -> Option<u8> {
        self.hydrogen_count
    }
}

impl core::fmt::Display for Atom {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        if self.wildcard {
            return write!(f, "*");
        }
        let symbol = if self.aromatic {
            self.element.symbol().to_lowercase()
        } else {
            self.element.symbol().to_string()
        };
        match self.isotope {
            Some(iso) => write!(f, "[{iso}{symbol}]"),
            None => write!(f, "{symbol}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_atom_new() {
        let a = Atom::new(Element::C);
        assert_eq!(a.element, Element::C);
        assert_eq!(a.charge, 0);
        assert!(!a.aromatic);
        assert!(!a.wildcard);
        assert_eq!(a.hydrogen_count, None);
    }

    #[test]
    fn test_aromatic_atom() {
        let a = Atom::aromatic(Element::C);
        assert!(a.aromatic);
    }

    #[test]
    fn test_wildcard_atom() {
        let a = Atom::wildcard();
        assert!(a.wildcard);
        assert_eq!(format!("{a}"), "*");
    }
}
