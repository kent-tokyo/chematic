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

use chematic_core::{Molecule, AtomIdx};
use chematic_smiles::canonical::canonical_atom_order;
use layers::{formula, connection, hydrogen, charge, isotope, stereo};
use std::collections::HashMap;

/// Build a mapping from AtomIdx to InChI 1-indexed atom numbers (excluding H).
fn build_inchi_index(mol: &Molecule) -> HashMap<AtomIdx, usize> {
    let canonical_order = canonical_atom_order(mol);
    let mut inchi_index: HashMap<AtomIdx, usize> = HashMap::new();
    let mut inchi_num = 0;
    for &canon_idx in &canonical_order {
        let atom_idx = AtomIdx(canon_idx as u32);
        let atom = mol.atom(atom_idx);
        if atom.element.atomic_number() != 1 {
            inchi_num += 1;
            inchi_index.insert(atom_idx, inchi_num);
        }
    }
    inchi_index
}

/// Generate InChI string for a molecule.
///
/// Layers included: formula, connectivity (/c), hydrogen (/h), double-bond stereo (/b),
/// tetrahedral stereo (/t), charge (/q if net charge ≠ 0), isotope (/i if present).
pub fn inchi(mol: &Molecule) -> String {
    let mut result = String::from("InChI=1S/");
    let inchi_index = build_inchi_index(mol);

    // Formula layer (prefix)
    let formula_str = formula::formula_layer(mol);
    result.push_str(&formula_str);

    // Connectivity layer /c
    if let Some(c_layer) = connection::connectivity_layer(mol) {
        result.push_str("/c");
        result.push_str(&c_layer);
    }

    // Hydrogen layer /h
    if let Some(h_layer) = hydrogen::hydrogen_layer(mol) {
        result.push_str("/h");
        result.push_str(&h_layer);
    }

    // Double-bond stereo layer /b (E/Z)
    if let Some(b_layer) = stereo::ez_stereo_layer(mol, &inchi_index) {
        result.push_str("/b");
        result.push_str(&b_layer);
    }

    // Tetrahedral stereo layer /t (R/S)
    if let Some(t_layer) = stereo::tetrahedral_stereo_layer(mol, &inchi_index) {
        result.push_str("/t");
        result.push_str(&t_layer);
    }

    // Charge layer /q (conditional)
    if let Some(q_layer) = charge::charge_layer(mol) {
        result.push_str("/q");
        result.push_str(&q_layer);
    }

    // Isotope layer /i (conditional)
    if let Some(i_layer) = isotope::isotope_layer(mol) {
        result.push_str("/i");
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
        // Benzene should have ring closure: /c1-2-3-4-5-6-1/h1-6H
        assert!(inchi_str.contains("/c1-2-3-4-5-6-1"), "Benzene should have ring closure in connectivity");
        assert!(inchi_str.contains("/h1-6H"), "Benzene should have hydrogen layer");
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
