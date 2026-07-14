//! Atom type: a single atom in a molecule.

use crate::element::Element;

/// Tetrahedral chirality as specified in OpenSMILES.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Chirality {
    /// No chirality specified.
    #[default]
    None,
    /// `@` — counterclockwise (looking from the first neighbor).
    CounterClockwise,
    /// `@@` — clockwise.
    Clockwise,
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
