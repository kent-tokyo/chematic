//! InChI and InChIKey generation and parsing.
//!
//! # Default (pure Rust, WASM-compatible)
//!
//! [`inchi`] and [`inchi_key`] are a pure-Rust, FFI-free approximation.
//! They are **not bit-exact** with the IUPAC reference implementation.
//! Use them when WASM compatibility or zero C dependencies are required.
//!
//! # Standard InChI (`native-inchi` feature)
//!
//! Enable the `native-inchi` Cargo feature to build with the vendored
//! IUPAC InChI C library (v1.07.5, MIT license). This provides
//! [`standard_inchi`] and [`standard_inchi_key`] which produce output
//! bit-exact with the IUPAC reference. Requires a C compiler; not
//! available for WASM targets.
//!
//! ```toml
//! [dependencies]
//! chematic-inchi = { version = "0.4", features = ["native-inchi"] }
//! ```
//!
//! # Examples
//!
//! ```ignore
//! use chematic_smiles::parse;
//! use chematic_inchi::{inchi, parse_inchi};
//!
//! let mol = parse("c1ccccc1").expect("benzene");
//! let inchi_str = inchi(&mol);
//!
//! // Parse InChI back to Molecule
//! let mol2 = parse_inchi(&inchi_str).expect("parse");
//! assert_eq!(mol2.atom_count(), 6);
//! ```
//!
//! # Deduplication (`dedup` module)
//!
//! [`dedup`] provides a high-confidence molecule deduplication API: fast
//! [`chematic_smiles::canonical_smiles`] candidate grouping, verified by
//! native standard InChI (never the pure-Rust approximation above). See the
//! [`dedup`] module docs for the full identity-policy rationale.

pub mod dedup;
pub mod key;
pub mod layers;
pub mod parser;

#[cfg(feature = "native-inchi")]
pub mod native;
#[cfg(feature = "native-inchi")]
pub use native::{InchiError, standard_inchi, standard_inchi_key};

use chematic_core::{AtomIdx, Molecule};
use chematic_smiles::canonical::canonical_atom_order;
use layers::{charge, connection, formula, hydrogen, isotope, stereo};
use std::collections::HashMap;

/// Build a mapping from AtomIdx to InChI 1-indexed atom numbers (excluding H).
pub fn build_inchi_index(mol: &Molecule) -> HashMap<AtomIdx, usize> {
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

/// Generate an InChI-like string for a molecule using the pure-Rust implementation.
///
/// **Note**: this output is **not standard IUPAC InChI**. It uses a different canonical
/// ordering and does not match the IUPAC reference library. For standard-compliant
/// InChI, enable the `native-inchi` feature and use [`standard_inchi`] instead.
///
/// Layers: formula, /c, /h, /b, /t, /m, /s, /q, /i.
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

    // Relative stereo parity layer /m (for 2+ stereocenters)
    if let Some(m_layer) = stereo::relative_stereo_parity_layer(mol, &inchi_index) {
        result.push_str("/m");
        result.push_str(&m_layer);
    }

    // Stereo type layer /s (absolute=1, relative=2, racemic=3)
    if let Some(s_layer) = stereo::stereo_type_layer(mol) {
        result.push_str("/s");
        result.push_str(&s_layer);
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

/// Generate an InChIKey from an InChI string using the pure-Rust implementation.
///
/// **Note**: this uses a non-standard SHA-256-based algorithm and does **not** match
/// IUPAC standard InChIKeys. For standard-compliant keys, enable the `native-inchi`
/// feature and use [`standard_inchi_key`] instead.
pub fn inchi_key(inchi_str: &str) -> String {
    key::inchi_key(inchi_str)
}

/// Parse an InChI string back into a Molecule representation.
///
/// Parse an InChI string back into a Molecule.
///
/// Supports organic molecules with stereo layers (`/b`, `/t`, `/m`, `/s`),
/// isotope layers (`/i`), and charge layers (`/q`).
pub use parser::{InchiParseError, parse_inchi};

// ─── Reaction InChI ──────────────────────────────────────────────────────────

/// Simplified Reaction InChI (not IUPAC-standard, WASM-compatible).
///
/// Generates a unique reaction identifier by combining InChI strings for
/// reactants, products, and agents.  Components within each role are sorted
/// for canonicality.
///
/// **Format** (pseudo-IUPAC, not bit-identical to the official IUPAC RInChI library):
///
/// ```text
/// RInChI=1.00.1S/<R1>!<R2>/d+!<P1>!<P2>/d-[!<A1>/d=]
/// ```
///
/// where `<Rx>` is an InChI string without the `"InChI=1S/"` prefix,
/// `d+` marks the reactant block, `d-` the product block, and `d=` the agent block.
///
/// For standard-compliant RInChI, use the separate IUPAC RInChI C library.
pub fn reaction_inchi(rxn: &chematic_rxn::Reaction) -> String {
    fn strip(s: &str) -> String {
        s.strip_prefix("InChI=1S/")
            .or_else(|| s.strip_prefix("InChI=1/"))
            .unwrap_or(s)
            .to_owned()
    }

    let mut r_parts: Vec<String> = rxn.reactants.iter().map(|m| strip(&inchi(m))).collect();
    let mut p_parts: Vec<String> = rxn.products.iter().map(|m| strip(&inchi(m))).collect();
    let mut a_parts: Vec<String> = rxn.agents.iter().map(|m| strip(&inchi(m))).collect();

    r_parts.sort();
    p_parts.sort();
    a_parts.sort();

    let mut out = "RInChI=1.00.1S/".to_string();
    if r_parts.is_empty() {
        out.push_str("d+");
    } else {
        out.push_str(&r_parts.join("!"));
        out.push_str("/d+");
    }
    // Only emit the product block when products are non-empty.
    // An unconditional '!' before an empty product block would produce the
    // malformed token `!d-` that looks like a second reactant named "d-".
    if !p_parts.is_empty() {
        out.push('!');
        out.push_str(&p_parts.join("!"));
        out.push_str("/d-");
    }
    if !a_parts.is_empty() {
        out.push('!');
        out.push_str(&a_parts.join("!"));
        out.push_str("/d=");
    }
    out
}

/// Compute an RInChIKey from a `reaction_inchi()` string.
///
/// Uses the same SHA-256-based hashing as `inchi_key()`.  Not identical to the
/// IUPAC-standard RInChIKey.
pub fn reaction_inchi_key(rinchi: &str) -> String {
    key::inchi_key(rinchi)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_rxn::parse_reaction;
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
        // Benzene's 6 ring atoms are fully symmetric, so no single numbering is
        // uniquely "correct" (see test_connectivity_benzene in layers/connection.rs)
        // -- just check a /c layer with a ring closure is present.
        let c_layer = inchi_str
            .split('/')
            .find(|s| s.starts_with('c'))
            .unwrap_or_else(|| panic!("no /c layer in {inchi_str:?}"));
        assert!(
            c_layer[1..].split('-').count() == 7,
            "Benzene should have ring closure in connectivity: {c_layer:?}"
        );
        assert!(
            inchi_str.contains("/h1-6H"),
            "Benzene should have hydrogen layer"
        );
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

    #[test]
    fn test_inchi_l_alanine_with_stereo_layers() {
        let mol = parse("N[C@@H](C)C(=O)O").expect("L-alanine");
        let inchi_str = inchi(&mol);
        eprintln!("L-alanine InChI: {}", inchi_str);
        assert!(
            inchi_str.contains("/t"),
            "L-alanine should have /t layer (R/S)"
        );
        assert!(
            inchi_str.contains("/s1"),
            "L-alanine should have /s1 layer (absolute stereo)"
        );
        // Single stereocenter should NOT have /m layer
        assert!(
            !inchi_str.contains("/m"),
            "Single stereocenter should not have /m layer"
        );
    }

    #[test]
    fn test_inchi_tartaric_acid_with_relative_parity() {
        // Tartaric acid: 2R,3S configuration
        let mol = parse("C[C@H](O)[C@@H](O)C(=O)O").expect("tartaric acid");
        let inchi_str = inchi(&mol);
        eprintln!("Tartaric acid InChI: {}", inchi_str);
        assert!(
            inchi_str.contains("/t"),
            "Tartaric acid should have /t layer"
        );
        assert!(
            inchi_str.contains("/m"),
            "Two stereocenters should have /m layer"
        );
        assert!(
            inchi_str.contains("/s1"),
            "With chirality markers, should have /s1"
        );
    }

    // ── RInChI ───────────────────────────────────────────────────────────────

    #[test]
    fn rinchi_starts_with_prefix() {
        let rxn = parse_reaction("CC>>CCO").unwrap();
        let ri = reaction_inchi(&rxn);
        assert!(ri.starts_with("RInChI=1.00.1S/"), "got: {ri}");
    }

    #[test]
    fn rinchi_contains_d_plus_and_d_minus() {
        let rxn = parse_reaction("CC>>CCO").unwrap();
        let ri = reaction_inchi(&rxn);
        assert!(ri.contains("d+"), "missing d+ in: {ri}");
        assert!(ri.contains("d-"), "missing d- in: {ri}");
    }

    #[test]
    fn rinchi_two_reactions_differ() {
        let r1 = reaction_inchi(&parse_reaction("CC>>CCO").unwrap());
        let r2 = reaction_inchi(&parse_reaction("c1ccccc1>>c1ccccc1O").unwrap());
        assert_ne!(r1, r2, "different reactions must produce different RInChI");
    }

    #[test]
    fn rinchi_reactant_order_canonical() {
        // Reactant order should not matter.
        let r1 = reaction_inchi(&parse_reaction("CC.c1ccccc1>>CCc1ccccc1").unwrap());
        let r2 = reaction_inchi(&parse_reaction("c1ccccc1.CC>>CCc1ccccc1").unwrap());
        assert_eq!(r1, r2, "reactant order must not affect RInChI");
    }

    #[test]
    fn rinchi_key_is_non_empty() {
        let rxn = parse_reaction("CC>>CCO").unwrap();
        let ri = reaction_inchi(&rxn);
        let rk = reaction_inchi_key(&ri);
        assert!(!rk.is_empty());
    }

    #[test]
    fn test_inchi_ethane_no_stereo() {
        let mol = parse("CC").expect("ethane");
        let inchi_str = inchi(&mol);
        eprintln!("Ethane InChI: {}", inchi_str);
        assert!(!inchi_str.contains("/t"), "Ethane should not have /t layer");
        assert!(!inchi_str.contains("/m"), "Ethane should not have /m layer");
        assert!(
            inchi_str.contains("/s3"),
            "Achiral ethane should have /s3 (racemic)"
        );
    }
}
