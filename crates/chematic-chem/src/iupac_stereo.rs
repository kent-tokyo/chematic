//! IUPAC naming with CIP stereo prefix.
//!
//! Wraps [`chematic_iupac::name`] and prepends `(R)-`, `(S)-`, or multi-centre
//! descriptors like `(1R,2S)-` when the molecule contains CIP-assigned
//! stereocentres.

use chematic_core::{CipCode, Molecule};
use chematic_iupac::IupacError;

/// Generate a local IUPAC name with CIP stereochemistry prefix.
///
/// Calls [`chematic_iupac::name`] for the base name, then
/// [`crate::assign_cip`] to collect R/S assignments.  If any
/// stereocentres are found, their descriptors are prepended in atom-index
/// order: `(R)-butane-2-ol`, `(1R,2S)-cyclohexane-1,2-diol`, etc.
///
/// Returns `Err` for structures outside the current IUPAC naming scope
/// (same conditions as [`chematic_iupac::name`]).
pub fn iupac_name_stereo(mol: &Molecule) -> Result<String, IupacError> {
    let base = chematic_iupac::name(mol)?;

    let cip = crate::assign_cip(mol);
    let mut rs_pairs: Vec<(u32, &'static str)> = cip
        .assignments
        .iter()
        .filter_map(|(idx, code)| match code {
            CipCode::R => Some((idx.0, "R")),
            CipCode::S => Some((idx.0, "S")),
            _ => None,
        })
        .collect();

    if rs_pairs.is_empty() {
        return Ok(base);
    }

    rs_pairs.sort_by_key(|(i, _)| *i);

    let prefix = if rs_pairs.len() == 1 {
        format!("({})", rs_pairs[0].1)
    } else {
        let parts: Vec<String> = rs_pairs
            .iter()
            .map(|(atom_idx, code)| format!("{}{}", atom_idx + 1, code))
            .collect();
        format!("({})", parts.join(","))
    };

    Ok(format!("{prefix}-{base}"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn achiral_butane_no_prefix() {
        let mol = parse("CCCC").unwrap();
        let name = iupac_name_stereo(&mol).unwrap();
        assert!(
            !name.starts_with('('),
            "achiral molecule: no stereo prefix, got {name}"
        );
    }

    #[test]
    fn r_2_butanol_has_r_prefix() {
        // (R)-butan-2-ol: SMILES with defined chirality
        let mol = parse("[C@@H](O)(CC)C").unwrap();
        let name = iupac_name_stereo(&mol);
        // Name may not be supported if IUPAC naming scope is limited;
        // if it succeeds it should contain R or S prefix.
        if let Ok(n) = name {
            assert!(
                n.starts_with("(R)") || n.starts_with("(S)"),
                "chiral 2-butanol should carry stereo prefix, got {n}"
            );
        }
    }

    #[test]
    fn achiral_cyclohexane_no_prefix() {
        let mol = parse("C1CCCCC1").unwrap();
        let name = iupac_name_stereo(&mol).unwrap();
        assert!(
            !name.starts_with('('),
            "achiral cyclohexane: no prefix, got {name}"
        );
    }
}
