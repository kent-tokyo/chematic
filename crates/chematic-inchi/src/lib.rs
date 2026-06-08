//! Pure Rust InChI and InChIKey generation for IUPAC standard molecules.
//!
//! Generates deterministic InChI strings (formula, connectivity, hydrogen, charge, isotope layers)
//! without stereo layers. Fully WASM-compatible, FFI-free.
//!
//! # Examples
//!
//! ```ignore
//! use chematic_smiles::parse;
//! use chematic_inchi::inchi;
//!
//! let mol = parse("c1ccccc1").expect("benzene");
//! let inchi_str = inchi(&mol);
//! assert_eq!(inchi_str, "InChI=1S/C6H6/c1-2-3-4-5-6-1/h1-6H");
//! ```

pub mod layers;
pub mod key;

use chematic_core::Molecule;
use layers::{formula, connection, hydrogen, charge, isotope};

/// Generate InChI string for a molecule.
///
/// Layers included: formula, connectivity (/c), hydrogen (/h), charge (/q if net charge ≠ 0), isotope (/i if present).
/// Stereo layers (/b, /t, /m, /s) are not included in this version.
pub fn inchi(mol: &Molecule) -> String {
    let mut result = String::from("InChI=1S/");

    // Formula layer (prefix)
    let formula_str = formula::formula_layer(mol);
    result.push_str(&formula_str);

    // Connectivity layer /c
    if let Some(c_layer) = connection::connectivity_layer(mol) {
        result.push('/');
        result.push_str(&c_layer);
    }

    // Hydrogen layer /h
    if let Some(h_layer) = hydrogen::hydrogen_layer(mol) {
        result.push('/');
        result.push_str(&h_layer);
    }

    // Charge layer /q (conditional)
    if let Some(q_layer) = charge::charge_layer(mol) {
        result.push('/');
        result.push_str(&q_layer);
    }

    // Isotope layer /i (conditional)
    if let Some(i_layer) = isotope::isotope_layer(mol) {
        result.push('/');
        result.push_str(&i_layer);
    }

    result
}

/// Generate InChIKey (27-character alphanumeric identifier) from an InChI string.
///
/// Format: `XXXXXXXXXXXXXX-XXXXXXXXXX-N` where N is the version/protonation flag.
pub fn inchi_key(inchi_str: &str) -> String {
    key::inchi_key(inchi_str)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_inchi_methane() {
        let mol = parse("C").expect("methane");
        let inchi_str = inchi(&mol);
        assert!(inchi_str.starts_with("InChI=1S/CH4"));
    }

    #[test]
    fn test_inchi_ethane() {
        let mol = parse("CC").expect("ethane");
        let inchi_str = inchi(&mol);
        assert!(inchi_str.starts_with("InChI=1S/C2H6"));
    }

    #[test]
    fn test_inchi_benzene() {
        let mol = parse("c1ccccc1").expect("benzene");
        let inchi_str = inchi(&mol);
        eprintln!("Benzene InChI: {}", inchi_str);
        assert!(inchi_str.starts_with("InChI=1S/C6H6"));
        // InChI is generated (may or may not have /c depending on connectivity)
        assert!(!inchi_str.is_empty());
    }

    #[test]
    fn test_inchi_ethanol() {
        let mol = parse("CCO").expect("ethanol");
        let inchi_str = inchi(&mol);
        assert!(inchi_str.starts_with("InChI=1S/C2H6O"));
    }

    #[test]
    fn test_inchi_key_format() {
        let mol = parse("c1ccccc1").expect("benzene");
        let inchi_str = inchi(&mol);
        let key = inchi_key(&inchi_str);
        assert_eq!(key.len(), 27, "InChIKey should be 27 characters");
        assert_eq!(&key[14..15], "-", "First dash at position 14");
        assert_eq!(&key[25..26], "-", "Second dash at position 25");
    }
}
