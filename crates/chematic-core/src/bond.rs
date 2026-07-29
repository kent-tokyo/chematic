//! Bond types and the bond-entry struct for molecular graphs.

/// Bond order, covering all SMILES bond types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondOrder {
    /// Single bond `-`
    Single,
    /// Double bond `=`
    Double,
    /// Triple bond `#`
    Triple,
    /// Quadruple bond `$` (used in some metal complexes)
    Quadruple,
    /// Aromatic bond `:` (pre-Kekulization representation)
    Aromatic,
    /// Stereo bond `/` (E/Z specification)
    Up,
    /// Stereo bond `\` (E/Z specification)
    Down,
    /// Zero-order bond used by CTAB/CXSMILES for non-valence connections.
    Zero,
    /// Dative/coordinate bond. The stored `atom1 → atom2` direction is the donor→acceptor direction.
    Dative,
    /// Query bond matching any bond (`~` / MDL type 8).
    QueryAny,
    /// Query bond matching single or double bonds (MDL type 5).
    QuerySingleOrDouble,
    /// Query bond matching single or aromatic bonds (MDL type 6).
    QuerySingleOrAromatic,
    /// Query bond matching double or aromatic bonds (MDL type 7).
    QueryDoubleOrAromatic,
}

impl BondOrder {
    /// Fractional bond order. Returns `None` for aromatic (non-integer).
    pub fn order_value(self) -> Option<f32> {
        match self {
            Self::Single => Some(1.0),
            Self::Double => Some(2.0),
            Self::Triple => Some(3.0),
            Self::Quadruple => Some(4.0),
            Self::Zero => Some(0.0),
            Self::Dative => Some(1.0),
            Self::Aromatic => None,
            Self::Up | Self::Down => Some(1.0),
            Self::QueryAny
            | Self::QuerySingleOrDouble
            | Self::QuerySingleOrAromatic
            | Self::QueryDoubleOrAromatic => None,
        }
    }

    /// Integer bond order used for valence calculations.
    /// Aromatic, stereo, dative, and query bonds count conservatively as 1.
    pub fn order_int(self) -> u8 {
        match self {
            Self::Zero => 0,
            Self::Single
            | Self::Up
            | Self::Down
            | Self::Aromatic
            | Self::Dative
            | Self::QueryAny
            | Self::QuerySingleOrDouble
            | Self::QuerySingleOrAromatic
            | Self::QueryDoubleOrAromatic => 1,
            Self::Double => 2,
            Self::Triple => 3,
            Self::Quadruple => 4,
        }
    }

    /// The SMILES/CXSMILES token that represents this bond order when available.
    pub fn smiles_token(self) -> &'static str {
        match self {
            Self::Single => "-",
            Self::Double => "=",
            Self::Triple => "#",
            Self::Quadruple => "$",
            Self::Aromatic => ":",
            Self::Up => "/",
            Self::Down => "\\",
            Self::Dative => "->",
            Self::Zero | Self::QueryAny => "~",
            Self::QuerySingleOrDouble
            | Self::QuerySingleOrAromatic
            | Self::QueryDoubleOrAromatic => "~",
        }
    }

    /// The **first** SMILES character of this bond order's token.
    ///
    /// Footgun: [`Self::smiles_token`] is not always one character.
    /// `Dative`'s token is `"->"`, so this returns a bare `'-'` for it —
    /// a plain single bond, which is both wrong and silently plausible.
    /// Use `smiles_token()` when emitting SMILES text (issue #194); a
    /// dative bond additionally needs its arrow oriented against the write
    /// order, see `chematic_smiles::writer::bond_token_from`.
    pub fn smiles_char(self) -> char {
        self.smiles_token().as_bytes()[0] as char
    }
}

impl core::fmt::Display for BondOrder {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "{}", self.smiles_token())
    }
}

/// A single bond edge in the molecular graph.
///
/// The graph is undirected; for most bond orders, `atom1`/`atom2` ordering has
/// no chemical meaning. For `BondOrder::Dative`, the ordering *is* meaningful:
/// `atom1` is the donor and `atom2` is the acceptor (see its own doc comment).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BondEntry {
    pub atom1: crate::molecule::AtomIdx,
    pub atom2: crate::molecule::AtomIdx,
    pub order: BondOrder,
}

impl BondEntry {
    /// Return the atom on the other side of this bond, or `None` if `idx` is not an endpoint.
    pub fn other(&self, idx: crate::molecule::AtomIdx) -> Option<crate::molecule::AtomIdx> {
        if self.atom1 == idx {
            Some(self.atom2)
        } else if self.atom2 == idx {
            Some(self.atom1)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::AtomIdx;

    #[test]
    fn test_bond_order_values() {
        assert_eq!(BondOrder::Single.order_value(), Some(1.0));
        assert_eq!(BondOrder::Double.order_value(), Some(2.0));
        assert_eq!(BondOrder::Triple.order_value(), Some(3.0));
        assert_eq!(BondOrder::Aromatic.order_value(), None);
    }

    #[test]
    fn test_bond_other() {
        let b = BondEntry {
            atom1: AtomIdx(0),
            atom2: AtomIdx(1),
            order: BondOrder::Single,
        };
        assert_eq!(b.other(AtomIdx(0)), Some(AtomIdx(1)));
        assert_eq!(b.other(AtomIdx(1)), Some(AtomIdx(0)));
        assert_eq!(b.other(AtomIdx(2)), None);
    }
}
