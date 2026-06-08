//! Chemical abbreviation expansion.
//!
//! Provides a lookup table of common chemical shorthand symbols (Ph, Me, Et, …)
//! and functions to expand them into `Molecule` graphs.

use chematic_core::Molecule;

/// Expand a chemical abbreviation symbol into a `Molecule`.
///
/// Returns `Some(mol)` if the symbol is in the built-in table, or `None`
/// if it is not recognised.  The returned molecule is a fragment — it may
/// contain wildcard (`[*]`) attachment-point atoms depending on the
/// abbreviation definition.
///
/// # Examples
/// ```rust,ignore
/// use chematic_chem::expand_abbreviation;
/// let phenyl = expand_abbreviation("Ph").unwrap();
/// assert_eq!(phenyl.atom_count(), 6);
/// ```
pub fn expand_abbreviation(symbol: &str) -> Option<Molecule> {
    let smiles = ABBREVIATIONS
        .iter()
        .find(|(s, _)| *s == symbol)
        .map(|(_, sm)| *sm)?;
    chematic_smiles::parse(smiles).ok()
}

/// Return the built-in abbreviation table as a slice of `(symbol, SMILES)` pairs.
///
/// The SMILES strings describe the expanded fragment.  Attachment points are
/// encoded as `[*]` (wildcard) atoms.
pub fn abbreviations() -> &'static [(&'static str, &'static str)] {
    ABBREVIATIONS
}

// ---------------------------------------------------------------------------
// Static table
// ---------------------------------------------------------------------------

/// `(symbol, SMILES)` pairs — sorted alphabetically by symbol for readability.
static ABBREVIATIONS: &[(&str, &str)] = &[
    // Alkyl groups
    ("Me", "C"),
    ("Et", "CC"),
    ("Pr", "CCC"),
    ("iPr", "CC(C)"),
    ("Bu", "CCCC"),
    ("iBu", "CC(C)C"),
    ("tBu", "CC(C)(C)C"),
    ("nBu", "CCCC"),
    ("sBu", "CCC(C)"),
    ("Pent", "CCCCC"),
    ("Hex", "CCCCCC"),
    // Acyl / protecting groups
    ("Ac", "CC(=O)"),
    ("Boc", "CC(C)(C)OC(=O)"),
    ("Cbz", "O=C(OC[*])c1ccccc1"),
    ("Fmoc", "O=C(O[*])OCC1c2ccccc2-c2ccccc21"),
    // Aromatic / benzyl
    ("Ph", "c1ccccc1"),
    ("Bn", "Cc1ccccc1"),
    ("Np", "c1ccc2ccccc2c1"), // naphthyl
    ("Py", "c1ccncc1"),       // pyridyl
    // Functional groups
    ("OMe", "OC"),
    ("OEt", "OCC"),
    ("NHMe", "NC"),
    ("NMe2", "N(C)C"),
    ("CF3", "C(F)(F)F"),
    ("CCl3", "C(Cl)(Cl)Cl"),
    ("CN", "C#N"),
    ("NO2", "[N+](=O)[O-]"),
    // Silyl groups
    ("TMS", "[Si](C)(C)C"),
    ("TBS", "[Si](C)(C)CC(C)(C)C"),
    // Sulfonyl
    ("Ts", "Cc1ccc(cc1)S(=O)(=O)"),
    ("Ms", "CS(=O)(=O)"),
    ("Tf", "C(F)(F)(F)S(=O)(=O)"),
];

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_expand_me() {
        let mol = expand_abbreviation("Me").unwrap();
        assert_eq!(mol.atom_count(), 1);
    }

    #[test]
    fn test_expand_ph() {
        let mol = expand_abbreviation("Ph").unwrap();
        assert_eq!(mol.atom_count(), 6);
    }

    #[test]
    fn test_expand_cf3() {
        let mol = expand_abbreviation("CF3").unwrap();
        // C + 3 F = 4 heavy atoms
        assert_eq!(mol.atom_count(), 4);
    }

    #[test]
    fn test_unknown_returns_none() {
        assert!(expand_abbreviation("Xyz123").is_none());
    }

    #[test]
    fn test_abbreviations_slice_non_empty() {
        assert!(!abbreviations().is_empty());
    }
}
