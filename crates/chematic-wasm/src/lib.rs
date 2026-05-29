//! `chematic-wasm` — WebAssembly bindings for the chematic cheminformatics library.
//!
//! Exposes a small, ergonomic API for parsing SMILES and computing molecular
//! descriptors from JavaScript/TypeScript via `wasm-bindgen`.

use wasm_bindgen::prelude::*;

// ---------------------------------------------------------------------------
// MolHandle
// ---------------------------------------------------------------------------

/// A handle to a parsed molecule.  Owns the molecule behind an `Rc` so that
/// it can be cheaply cloned on the JS side without copying atom/bond data.
#[wasm_bindgen]
pub struct MolHandle {
    inner: std::rc::Rc<chematic_core::Molecule>,
}

#[wasm_bindgen]
impl MolHandle {
    /// Number of heavy atoms (explicit atoms in the graph; does not count implicit H).
    pub fn atom_count(&self) -> usize {
        self.inner.atom_count()
    }

    /// Number of bonds.
    pub fn bond_count(&self) -> usize {
        self.inner.bond_count()
    }

    /// Molecular formula string (Hill notation: C first, H second, then alphabetical).
    pub fn formula(&self) -> String {
        molecular_formula(&self.inner)
    }

    /// Canonical SMILES string.
    pub fn canonical_smiles(&self) -> String {
        chematic_smiles::canonical_smiles(&self.inner)
    }

    /// Average molecular weight (Da).
    pub fn molecular_weight(&self) -> f64 {
        chematic_chem::molecular_weight(&self.inner)
    }

    /// Topological polar surface area (Å²).
    pub fn tpsa(&self) -> f64 {
        chematic_chem::tpsa(&self.inner)
    }

    /// Returns `true` if the molecule satisfies Lipinski's Rule of Five.
    pub fn lipinski_passes(&self) -> bool {
        chematic_chem::lipinski_passes(&self.inner)
    }

    /// Number of non-hydrogen heavy atoms.
    pub fn heavy_atom_count(&self) -> usize {
        chematic_chem::heavy_atom_count(&self.inner)
    }

    /// Number of hydrogen bond donors (N-H or O-H groups).
    pub fn hbd_count(&self) -> usize {
        chematic_chem::hbd_count(&self.inner)
    }

    /// Number of hydrogen bond acceptors (Lipinski: all N and O atoms).
    pub fn hba_count(&self) -> usize {
        chematic_chem::hba_count(&self.inner)
    }

    /// Crippen–Wildman octanol/water partition coefficient (LogP).
    pub fn logp_crippen(&self) -> f64 {
        chematic_chem::logp_crippen(&self.inner)
    }

    /// Fraction of sp3 carbons (Fsp3).
    pub fn fsp3(&self) -> f64 {
        chematic_chem::fsp3(&self.inner)
    }

    /// Number of aromatic rings (all ring atoms aromatic).
    pub fn aromatic_ring_count(&self) -> usize {
        chematic_chem::aromatic_ring_count(&self.inner)
    }

    /// Quantitative Estimate of Drug-likeness (QED); range [0, 1].
    pub fn qed(&self) -> f64 {
        chematic_chem::qed(&self.inner)
    }

    /// Monoisotopic (exact) mass.
    pub fn exact_mass(&self) -> f64 {
        chematic_chem::exact_mass(&self.inner)
    }

    /// Number of rotatable bonds.
    pub fn rotatable_bond_count(&self) -> usize {
        chematic_chem::rotatable_bond_count(&self.inner)
    }

    /// Wildman–Crippen molar refractivity (MR).
    pub fn molar_refractivity(&self) -> f64 {
        chematic_chem::molar_refractivity(&self.inner)
    }

    /// Sum of formal charges.
    pub fn formal_charge_sum(&self) -> i32 {
        chematic_chem::formal_charge_sum(&self.inner)
    }

    /// Returns `true` if the molecule passes Veber's oral bioavailability criteria
    /// (TPSA ≤ 140 Å² and rotatable bonds ≤ 10).
    pub fn veber_passes(&self) -> bool {
        chematic_chem::veber_passes(&self.inner)
    }

    /// Returns `true` if the molecule passes Egan's absorption criteria
    /// (TPSA ≤ 131.6 Å² and LogP ≤ 5.88).
    pub fn egan_passes(&self) -> bool {
        chematic_chem::egan_passes(&self.inner)
    }

    /// Returns `true` if the molecule passes the REOS (Rapid Elimination Of Swill) filter.
    pub fn reos_passes(&self) -> bool {
        chematic_chem::reos_passes(&self.inner)
    }

    /// Returns `true` if the molecule passes Ghose's drug-likeness filter
    /// (MW 160–480, LogP −0.4–5.6, HeavyAtoms 20–70, MR 40–130).
    pub fn ghose_passes(&self) -> bool {
        chematic_chem::ghose_passes(&self.inner)
    }

    /// Number of heteroatoms (non-C, non-H heavy atoms).
    pub fn num_heteroatoms(&self) -> usize {
        chematic_chem::num_heteroatoms(&self.inner)
    }

    /// Total number of rings (SSSR count).
    pub fn ring_count(&self) -> usize {
        chematic_chem::ring_count(&self.inner)
    }

    /// Number of assigned stereocenters (R/S).
    pub fn num_stereocenters(&self) -> usize {
        chematic_chem::num_stereocenters(&self.inner)
    }

    /// Returns `true` if the molecule has no PAINS structural alerts.
    pub fn pains_passes(&self) -> bool {
        chematic_chem::pains_passes(&self.inner)
    }

    /// Number of aromatic rings containing at least one heteroatom (N, O, S, …).
    pub fn num_aromatic_heterocycles(&self) -> usize {
        chematic_chem::num_aromatic_heterocycles(&self.inner)
    }

    /// Number of non-aromatic rings containing at least one heteroatom.
    pub fn num_aliphatic_heterocycles(&self) -> usize {
        chematic_chem::num_aliphatic_heterocycles(&self.inner)
    }

    /// Number of fully saturated rings containing at least one heteroatom.
    pub fn num_saturated_heterocycles(&self) -> usize {
        chematic_chem::num_saturated_heterocycles(&self.inner)
    }

    /// Number of spiro atoms (sole shared atom between exactly 2 rings).
    pub fn num_spiro_atoms(&self) -> usize {
        chematic_chem::num_spiro_atoms(&self.inner)
    }

    /// Number of bridgehead atoms (shared by ≥2 rings with ≥3 ring bonds).
    pub fn num_bridgehead_atoms(&self) -> usize {
        chematic_chem::num_bridgehead_atoms(&self.inner)
    }

    /// 2D SVG depiction of the molecule (CPK coloring).
    pub fn depict_svg(&self) -> String {
        chematic_depict::depict_svg(&self.inner)
    }
}

// ---------------------------------------------------------------------------
// Free functions exported to JS
// ---------------------------------------------------------------------------

/// Parse a SMILES string into a `MolHandle`.
///
/// Returns a JS error string on parse failure.
#[wasm_bindgen]
pub fn parse_smiles(s: &str) -> Result<MolHandle, JsValue> {
    chematic_smiles::parse(s)
        .map(|mol| MolHandle { inner: std::rc::Rc::new(mol) })
        .map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Tanimoto similarity between two molecules using ECFP4 fingerprints.
#[wasm_bindgen]
pub fn tanimoto_ecfp4(a: &MolHandle, b: &MolHandle) -> f64 {
    chematic_fp::tanimoto_ecfp4(&a.inner, &b.inner)
}

/// Tanimoto similarity between two molecules using FCFP4 fingerprints (pharmacophore-based).
#[wasm_bindgen]
pub fn tanimoto_fcfp4(a: &MolHandle, b: &MolHandle) -> f64 {
    chematic_fp::tanimoto_fcfp4(&a.inner, &b.inner)
}

/// Tanimoto similarity between two molecules using AtomPair fingerprints.
#[wasm_bindgen]
pub fn tanimoto_atom_pair(a: &MolHandle, b: &MolHandle) -> f64 {
    let fa = chematic_fp::atom_pair_fp(&a.inner);
    let fb = chematic_fp::atom_pair_fp(&b.inner);
    fa.tanimoto(&fb)
}

/// Tanimoto similarity between two molecules using Topological Torsion fingerprints.
#[wasm_bindgen]
pub fn tanimoto_torsion(a: &MolHandle, b: &MolHandle) -> f64 {
    let fa = chematic_fp::torsion_fp(&a.inner);
    let fb = chematic_fp::torsion_fp(&b.inner);
    fa.tanimoto(&fb)
}

/// Number of BRICS fragments produced by fragmenting the molecule.
///
/// Returns 1 if no BRICS-breakable bonds exist (whole molecule is one fragment).
#[wasm_bindgen]
pub fn brics_fragment_count(mol: &MolHandle) -> usize {
    chematic_chem::brics_fragments(&mol.inner).len()
}

/// Compute the ECFP4 fingerprint as a bit-packed byte vector (256 bytes = 2048 bits).
#[wasm_bindgen]
pub fn ecfp4_bitvec(mol: &MolHandle) -> Vec<u8> {
    let fp = chematic_fp::ecfp4(&mol.inner);
    // BitVec2048 is 2048 bits; extract them byte-by-byte via the public `get` method.
    (0..256usize)
        .map(|byte_idx| {
            let mut byte = 0u8;
            for bit in 0..8usize {
                if fp.get(byte_idx * 8 + bit) {
                    byte |= 1 << bit;
                }
            }
            byte
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Private helper: molecular formula (Hill notation)
// ---------------------------------------------------------------------------

/// Build a molecular formula string in Hill notation.
///
/// Hill convention: carbon first, hydrogen second, remaining elements
/// in alphabetical order.  Implicit hydrogens (from valence model) are
/// included in the count.
fn molecular_formula(mol: &chematic_core::Molecule) -> String {
    use chematic_core::{Element, implicit_hcount};
    use std::collections::BTreeMap;

    let mut counts: BTreeMap<u8, u32> = BTreeMap::new();

    for (idx, atom) in mol.atoms() {
        let an = atom.element.atomic_number();
        if an != 1 {
            // Count the heavy atom.
            *counts.entry(an).or_insert(0) += 1;
            // Add its implicit hydrogens.
            let h = implicit_hcount(mol, idx) as u32;
            if h > 0 {
                *counts.entry(1).or_insert(0) += h;
            }
        } else {
            // Explicit hydrogen atom.
            *counts.entry(1).or_insert(0) += 1;
        }
    }

    // Collect into Hill order: C (6), H (1), then remaining by atomic number.
    let mut result = String::new();
    let append = |symbol: &str, count: u32, out: &mut String| {
        out.push_str(symbol);
        if count > 1 {
            out.push_str(&count.to_string());
        }
    };

    if let Some(&c_count) = counts.get(&6) {
        append("C", c_count, &mut result);
    }
    if let Some(&h_count) = counts.get(&1) {
        append("H", h_count, &mut result);
    }
    // Remaining elements in atomic-number order (BTreeMap is sorted by key).
    for (&an, &count) in &counts {
        if an == 1 || an == 6 {
            continue;
        }
        let elem = Element::from_atomic_number(an).unwrap();
        append(elem.symbol(), count, &mut result);
    }

    result
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(s: &str) -> MolHandle {
        MolHandle { inner: std::rc::Rc::new(chematic_smiles::parse(s).unwrap()) }
    }

    #[test]
    fn parse_benzene_atom_count() {
        assert_eq!(parse("c1ccccc1").atom_count(), 6);
    }

    #[test]
    fn canonical_smiles_benzene() {
        let mol = parse("c1ccccc1");
        let cs = mol.canonical_smiles();
        assert!(!cs.is_empty());
    }

    #[test]
    fn molecular_weight_aspirin() {
        let mw = parse("CC(=O)Oc1ccccc1C(=O)O").molecular_weight();
        assert!((mw - 180.16).abs() < 1.0);
    }

    #[test]
    fn lipinski_aspirin() {
        assert!(parse("CC(=O)Oc1ccccc1C(=O)O").lipinski_passes());
    }

    #[test]
    fn tanimoto_same_mol() {
        let a = parse("c1ccccc1");
        let b = parse("c1ccccc1");
        let sim = tanimoto_ecfp4(&a, &b);
        assert!((sim - 1.0).abs() < 1e-6);
    }

    #[test]
    fn tanimoto_different() {
        let a = parse("c1ccccc1");
        let b = parse("CC(=O)Oc1ccccc1C(=O)O");
        assert!(tanimoto_ecfp4(&a, &b) < 1.0);
    }

    #[test]
    fn heavy_atom_count_ethanol() {
        assert_eq!(parse("CCO").heavy_atom_count(), 3);
    }

    #[test]
    fn logp_crippen_aspirin_range() {
        let lp = parse("CC(=O)Oc1ccccc1C(=O)O").logp_crippen();
        assert!(lp > 0.5 && lp < 3.5, "aspirin LogP = {lp:.3}");
    }

    #[test]
    fn fsp3_benzene_zero() {
        assert_eq!(parse("c1ccccc1").fsp3(), 0.0, "benzene Fsp3 = 0");
    }

    #[test]
    fn fsp3_cyclohexane_one() {
        assert_eq!(parse("C1CCCCC1").fsp3(), 1.0, "cyclohexane Fsp3 = 1");
    }

    #[test]
    fn aromatic_ring_count_benzene() {
        assert_eq!(parse("c1ccccc1").aromatic_ring_count(), 1);
    }

    #[test]
    fn qed_aspirin_range() {
        let q = parse("CC(=O)Oc1ccccc1C(=O)O").qed();
        assert!(q > 0.0 && q <= 1.0, "aspirin QED = {q:.3}");
    }

    #[test]
    fn exact_mass_aspirin() {
        // Aspirin monoisotopic mass: 180.0423
        let em = parse("CC(=O)Oc1ccccc1C(=O)O").exact_mass();
        assert!((em - 180.042).abs() < 0.01, "aspirin exact mass = {em:.4}");
    }

    #[test]
    fn rotatable_bond_count_aspirin() {
        // Aspirin has 3 rotatable bonds (OC, C=O ester, and COOH)
        let rb = parse("CC(=O)Oc1ccccc1C(=O)O").rotatable_bond_count();
        assert!(rb >= 2 && rb <= 5, "aspirin rotbonds = {rb}");
    }

    #[test]
    fn tanimoto_atom_pair_same_mol() {
        let a = parse("c1ccccc1");
        let b = parse("c1ccccc1");
        assert!((tanimoto_atom_pair(&a, &b) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn tanimoto_torsion_same_mol() {
        let a = parse("CCCC");
        let b = parse("CCCC");
        assert!((tanimoto_torsion(&a, &b) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn brics_fragment_count_benzene() {
        assert_eq!(brics_fragment_count(&parse("c1ccccc1")), 1);
    }

    #[test]
    fn brics_fragment_count_aspirin() {
        assert!(brics_fragment_count(&parse("CC(=O)Oc1ccccc1C(=O)O")) >= 2);
    }
}
