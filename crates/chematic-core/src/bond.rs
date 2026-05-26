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
}

impl BondOrder {
    /// Fractional bond order. Returns `None` for aromatic (non-integer).
    pub fn order_value(self) -> Option<f32> {
        match self {
            Self::Single    => Some(1.0),
            Self::Double    => Some(2.0),
            Self::Triple    => Some(3.0),
            Self::Quadruple => Some(4.0),
            Self::Aromatic  => None,
            Self::Up | Self::Down => Some(1.0),
        }
    }

    /// Integer bond order used for valence calculations.
    /// Aromatic and stereo bonds count as 1.
    pub fn order_int(self) -> u8 {
        match self {
            Self::Single | Self::Up | Self::Down | Self::Aromatic => 1,
            Self::Double    => 2,
            Self::Triple    => 3,
            Self::Quadruple => 4,
        }
    }

    /// The SMILES character that represents this bond order.
    pub fn smiles_char(self) -> char {
        match self {
            Self::Single    => '-',
            Self::Double    => '=',
            Self::Triple    => '#',
            Self::Quadruple => '$',
            Self::Aromatic  => ':',
            Self::Up        => '/',
            Self::Down      => '\\',
        }
    }
}

impl core::fmt::Display for BondOrder {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "{}", self.smiles_char())
    }
}

/// A single bond edge in the molecular graph.
///
/// The graph is undirected; `atom1`/`atom2` ordering is not guaranteed.
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
        let b = BondEntry { atom1: AtomIdx(0), atom2: AtomIdx(1), order: BondOrder::Single };
        assert_eq!(b.other(AtomIdx(0)), Some(AtomIdx(1)));
        assert_eq!(b.other(AtomIdx(1)), Some(AtomIdx(0)));
        assert_eq!(b.other(AtomIdx(2)), None);
    }
}
