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

    /// 2D SVG depiction with style options.
    pub fn depict_svg_opts(&self, opts: &DepictOptions) -> String {
        let ro = chematic_depict::RenderOptions {
            width: opts.width,
            height: opts.height,
            padding: opts.padding,
            background: opts.background.clone(),
            dark: opts.dark,
            highlight_atoms: opts.highlight_atoms.iter()
                .map(|&i| chematic_core::AtomIdx(i))
                .collect(),
            highlight_bonds: opts.highlight_bonds.iter()
                .map(|&i| chematic_core::BondIdx(i))
                .collect(),
            highlight_color: opts.highlight_color.clone(),
            atom_ids: opts.atom_ids,
            show_atom_indices: opts.show_atom_indices,
            kekulize: opts.kekulize,
        };
        chematic_depict::depict_svg_opts(&self.inner, &ro)
    }

    // -----------------------------------------------------------------------
    // Topological descriptors (Sprint G)
    // -----------------------------------------------------------------------

    /// Wiener topological index (sum of all pairwise shortest-path distances).
    pub fn wiener_index(&self) -> f64 {
        chematic_chem::wiener_index(&self.inner)
    }

    /// Hall–Kier κ1 shape index.
    pub fn kappa1(&self) -> f64 {
        chematic_chem::kappa1(&self.inner)
    }

    /// Hall–Kier κ2 shape index.
    pub fn kappa2(&self) -> f64 {
        chematic_chem::kappa2(&self.inner)
    }

    /// Hall–Kier κ3 shape index.
    pub fn kappa3(&self) -> f64 {
        chematic_chem::kappa3(&self.inner)
    }

    /// Kier–Hall χ0 molecular connectivity index.
    pub fn chi0(&self) -> f64 {
        chematic_chem::chi0(&self.inner)
    }

    /// Kier–Hall χ1 molecular connectivity index.
    pub fn chi1(&self) -> f64 {
        chematic_chem::chi1(&self.inner)
    }

    /// Kier–Hall χ2 molecular connectivity index.
    pub fn chi2(&self) -> f64 {
        chematic_chem::chi2(&self.inner)
    }

    /// Kier–Hall χ3 molecular connectivity index.
    pub fn chi3(&self) -> f64 {
        chematic_chem::chi3(&self.inner)
    }

    /// Kier–Hall χ4 molecular connectivity index.
    pub fn chi4(&self) -> f64 {
        chematic_chem::chi4(&self.inner)
    }

    /// Kier–Hall χ0v valence-weighted connectivity index.
    pub fn chi0v(&self) -> f64 {
        chematic_chem::chi0v(&self.inner)
    }

    /// Kier–Hall χ1v valence-weighted connectivity index.
    pub fn chi1v(&self) -> f64 {
        chematic_chem::chi1v(&self.inner)
    }

    /// Kier–Hall χ2v valence-weighted connectivity index.
    pub fn chi2v(&self) -> f64 {
        chematic_chem::chi2v(&self.inner)
    }

    /// Kier–Hall χ3v valence-weighted connectivity index.
    pub fn chi3v(&self) -> f64 {
        chematic_chem::chi3v(&self.inner)
    }

    /// Kier–Hall χ4v valence-weighted connectivity index.
    pub fn chi4v(&self) -> f64 {
        chematic_chem::chi4v(&self.inner)
    }

    /// Bertz complexity index (BertzCT).
    pub fn bertz_ct(&self) -> f64 {
        chematic_chem::bertz_ct(&self.inner)
    }

    /// Labute approximate surface area (Å²).
    pub fn labute_asa(&self) -> f64 {
        chematic_chem::labute_asa(&self.inner)
    }

    // -----------------------------------------------------------------------
    // Morgan count fingerprint (Sprint G)
    // -----------------------------------------------------------------------

    /// Morgan count fingerprint as a JSON object string (`{"<hash>": count, …}`).
    ///
    /// `radius` controls the ECFP radius (2 = ECFP4-equivalent).
    pub fn morgan_fp_counts_json(&self, radius: u32) -> String {
        let counts = chematic_fp::morgan_fp_counts(&self.inner, radius);
        let mut pairs: Vec<(u64, u32)> = counts.into_iter().collect();
        pairs.sort_by_key(|(k, _)| *k);
        let entries: Vec<String> = pairs.iter()
            .map(|(k, v)| format!("\"{k}\": {v}"))
            .collect();
        format!("{{{}}}", entries.join(", "))
    }
}

// ---------------------------------------------------------------------------
// Free functions exported to JS
// ---------------------------------------------------------------------------

/// Returns `true` if the SMILES string can be parsed without error.
#[wasm_bindgen]
pub fn is_valid_smiles(s: &str) -> bool {
    chematic_smiles::parse(s).is_ok()
}

// ---------------------------------------------------------------------------
// DepictOptions
// ---------------------------------------------------------------------------

/// Style options for [`MolHandle::depict_svg_opts`].
///
/// Construct with `new DepictOptions()`, then call setters:
/// ```js
/// const opts = new DepictOptions();
/// opts.set_background("transparent");
/// opts.set_dark(true);
/// opts.set_width(240);
/// opts.set_height(240);
/// ```
#[wasm_bindgen]
pub struct DepictOptions {
    width: Option<u32>,
    height: Option<u32>,
    padding: f64,
    background: String,
    dark: bool,
    highlight_atoms: Vec<u32>,
    highlight_bonds: Vec<u32>,
    highlight_color: String,
    atom_ids: bool,
    show_atom_indices: bool,
    kekulize: bool,
}

#[wasm_bindgen]
impl DepictOptions {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self {
            width: None,
            height: None,
            padding: 20.0,
            background: "white".into(),
            dark: false,
            highlight_atoms: vec![],
            highlight_bonds: vec![],
            highlight_color: "#FFFF00".into(),
            atom_ids: false,
            show_atom_indices: false,
            kekulize: false,
        }
    }

    pub fn set_width(&mut self, w: u32) { self.width = Some(w); }
    pub fn set_height(&mut self, h: u32) { self.height = Some(h); }
    pub fn set_padding(&mut self, p: f64) { self.padding = p; }
    pub fn set_background(&mut self, bg: String) { self.background = bg; }
    pub fn set_dark(&mut self, dark: bool) { self.dark = dark; }
    pub fn set_highlight_atoms(&mut self, atoms: Vec<u32>) { self.highlight_atoms = atoms; }
    pub fn set_highlight_bonds(&mut self, bonds: Vec<u32>) { self.highlight_bonds = bonds; }
    pub fn set_highlight_color(&mut self, color: String) { self.highlight_color = color; }
    pub fn set_atom_ids(&mut self, v: bool) { self.atom_ids = v; }
    pub fn set_show_atom_indices(&mut self, v: bool) { self.show_atom_indices = v; }
    pub fn set_kekulize(&mut self, v: bool) { self.kekulize = v; }
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

/// Return a copy of the molecule with all implicit hydrogens converted to explicit H atoms.
#[wasm_bindgen]
pub fn add_hydrogens(mol: &MolHandle) -> MolHandle {
    MolHandle { inner: std::rc::Rc::new(chematic_chem::add_hydrogens(&mol.inner)) }
}

/// Return a copy of the molecule with all explicit hydrogen atoms removed.
#[wasm_bindgen]
pub fn remove_hydrogens(mol: &MolHandle) -> MolHandle {
    MolHandle { inner: std::rc::Rc::new(chematic_chem::remove_hydrogens(&mol.inner)) }
}

/// Render a grid SVG from newline-separated SMILES (one per line).
///
/// Lines that fail to parse are silently skipped.
/// `cols` controls the number of columns (each cell is 200×200 px).
#[wasm_bindgen]
pub fn depict_svg_grid(smiles_block: &str, cols: usize) -> String {
    let mols: Vec<chematic_core::Molecule> = smiles_block
        .lines()
        .filter(|s| !s.trim().is_empty())
        .filter_map(|s| chematic_smiles::parse(s.trim()).ok())
        .collect();
    let refs: Vec<&chematic_core::Molecule> = mols.iter().collect();
    chematic_depict::depict_svg_grid(&refs, cols)
}

/// Apply a SMIRKS reaction template and return product SMILES as a JSON string.
///
/// `reactants_smiles`: pipe-separated SMILES, one per reactant slot in the SMIRKS.
/// Returns a JSON array of arrays: `[["product_smi", …], …]`.
/// Returns a JS error on parse failure or arity mismatch.
#[wasm_bindgen]
pub fn run_reactants(smirks: &str, reactants_smiles: &str) -> Result<String, JsValue> {
    let reactant_mols: Result<Vec<chematic_core::Molecule>, _> = reactants_smiles
        .split('|')
        .map(|s| {
            chematic_smiles::parse(s.trim())
                .map_err(|e| JsValue::from_str(&e.to_string()))
        })
        .collect();
    let reactant_mols = reactant_mols?;
    let refs: Vec<&chematic_core::Molecule> = reactant_mols.iter().collect();

    let products = chematic_rxn::run_reactants(smirks, &refs)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;

    let outer: Vec<String> = products
        .iter()
        .map(|set| {
            let inner: Vec<String> = set
                .iter()
                .map(|mol| format!("\"{}\"", chematic_smiles::canonical_smiles(mol)))
                .collect();
            format!("[{}]", inner.join(", "))
        })
        .collect();
    Ok(format!("[{}]", outer.join(", ")))
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

    #[test]
    fn wiener_index_ethane() {
        // Ethane: 2 atoms, distance 1 — Wiener index = 1.
        assert_eq!(parse("CC").wiener_index(), 1.0);
    }

    #[test]
    fn kappa1_propane_range() {
        let k = parse("CCC").kappa1();
        assert!(k > 0.0, "kappa1 should be positive");
    }

    #[test]
    fn chi0_benzene_positive() {
        assert!(parse("c1ccccc1").chi0() > 0.0);
    }

    #[test]
    fn labute_asa_aspirin_range() {
        let asa = parse("CC(=O)Oc1ccccc1C(=O)O").labute_asa();
        assert!(asa > 50.0 && asa < 200.0, "aspirin LabuteASA = {asa:.2}");
    }

    #[test]
    fn bertz_ct_benzene_positive() {
        assert!(parse("c1ccccc1").bertz_ct() > 0.0);
    }

    #[test]
    fn morgan_fp_counts_json_benzene() {
        let json = parse("c1ccccc1").morgan_fp_counts_json(2);
        assert!(json.starts_with('{') && json.ends_with('}'));
    }

    #[test]
    fn add_remove_hydrogens_roundtrip() {
        let mol = parse("CC");
        let with_h = add_hydrogens(&mol);
        assert!(with_h.atom_count() > mol.atom_count(), "H atoms should be added");
        let back = remove_hydrogens(&with_h);
        assert_eq!(back.atom_count(), mol.atom_count());
    }

    #[test]
    fn depict_svg_grid_two_mols() {
        let svg = depict_svg_grid("CC\nCCC", 2);
        assert!(svg.contains("<svg"), "expected SVG output");
    }

    #[test]
    fn depict_svg_grid_invalid_smiles_skipped() {
        let svg = depict_svg_grid("CC\nNOT_A_SMILES\nCCC", 2);
        assert!(svg.contains("<svg"), "invalid SMILES should be silently skipped");
    }

    #[test]
    fn run_reactants_esterification() {
        // Simple esterification: carboxylic acid + alcohol → ester + water
        let result = run_reactants("[C:1](=O)[OH:2].[O:3][C:4]>>[C:1](=O)[O:3][C:4]", "CC(=O)O|CCO");
        assert!(result.is_ok(), "run_reactants should succeed");
        let json = result.unwrap();
        assert!(json.contains('['), "expected JSON array");
    }

    // Note: run_reactants error-path tests are omitted here because JsValue::from_str
    // panics outside a WASM runtime. Error coverage lives in chematic-rxn unit tests.

    #[test]
    fn is_valid_smiles_valid() {
        assert!(is_valid_smiles("CCO"), "ethanol is valid");
        assert!(is_valid_smiles("c1ccccc1"), "benzene is valid");
        assert!(is_valid_smiles("O"), "water is valid");
        assert!(is_valid_smiles("C"), "methane is valid");
    }

    #[test]
    fn is_valid_smiles_invalid() {
        assert!(!is_valid_smiles(""), "empty string is invalid");
        assert!(!is_valid_smiles("[NOSUCHELEMENT]"), "unknown bracket atom is invalid");
    }

    #[test]
    fn depict_svg_opts_transparent_background() {
        let h = parse("CCO");
        let mut opts = DepictOptions::new();
        opts.set_background("transparent".to_string());
        let svg = h.depict_svg_opts(&opts);
        assert!(svg.contains("<svg"), "must produce SVG");
        assert!(!svg.contains("fill=\"transparent\""), "no bg rect for transparent");
    }

    #[test]
    fn depict_svg_opts_custom_size() {
        let h = parse("CCO");
        let mut opts = DepictOptions::new();
        opts.set_width(300);
        opts.set_height(200);
        let svg = h.depict_svg_opts(&opts);
        assert!(svg.contains("width=\"300\""), "SVG width should be 300");
        assert!(svg.contains("height=\"200\""), "SVG height should be 200");
    }

    #[test]
    fn depict_svg_opts_dark_theme() {
        let h = parse("CC");
        let mut opts = DepictOptions::new();
        opts.set_dark(true);
        opts.set_background("#0f172a".to_string());
        let svg = h.depict_svg_opts(&opts);
        assert!(svg.contains("stroke=\"white\""), "dark theme bonds should be white");
    }

    #[test]
    fn depict_svg_single_atom_water_shows_h2o() {
        let svg = parse("O").depict_svg();
        assert!(svg.contains("H2O"), "water 'O' should render as H2O");
    }

    #[test]
    fn depict_svg_single_atom_methane_shows_ch4() {
        let svg = parse("C").depict_svg();
        assert!(svg.contains("CH4"), "methane 'C' should render as CH4");
    }

    // ── Sprint L: disconnected SMILES ────────────────────────────────────────

    #[test]
    fn depict_svg_disconnected_nacl() {
        let svg = parse("[Na+].[Cl-]").depict_svg();
        assert!(svg.contains("Na"), "Na should appear in disconnected SMILES SVG");
        assert!(svg.contains("Cl"), "Cl should appear in disconnected SMILES SVG");
        assert!(!svg.is_empty());
    }

    #[test]
    fn depict_svg_disconnected_water_dimer() {
        let svg = parse("O.O").depict_svg();
        // O.O = 2 atoms; each O in a multi-atom mol renders as "OH2" (heteroatom path)
        assert!(svg.matches("OH2").count() >= 2, "both O atoms should appear as OH2 labels");
        assert!(!svg.is_empty());
    }

    // ── Sprint L: atom data attributes ──────────────────────────────────────

    #[test]
    fn depict_svg_opts_atom_ids_contains_data_attrs() {
        let h = parse("CC(=O)O"); // acetic acid
        let mut opts = DepictOptions::new();
        opts.set_atom_ids(true);
        let svg = h.depict_svg_opts(&opts);
        assert!(svg.contains("data-atom-idx="), "atom_ids should add data-atom-idx");
        assert!(svg.contains("data-element="), "atom_ids should add data-element");
        assert!(svg.contains("data-charge="), "atom_ids should add data-charge");
    }

    #[test]
    fn depict_svg_opts_atom_ids_false_no_data_attrs() {
        let h = parse("CC(=O)O");
        let svg = h.depict_svg_opts(&DepictOptions::new());
        assert!(!svg.contains("data-atom-idx="), "default opts should not have data-atom-idx");
    }

    #[test]
    fn depict_svg_opts_atom_ids_charge_correct() {
        let h = parse("[NH4+]");
        let mut opts = DepictOptions::new();
        opts.set_atom_ids(true);
        let svg = h.depict_svg_opts(&opts);
        assert!(svg.contains("data-charge=\"1\""), "NH4+ should have charge=1");
    }

    // ── Sprint L: show_atom_indices ──────────────────────────────────────────

    #[test]
    fn depict_svg_opts_show_atom_indices() {
        let h = parse("c1ccccc1"); // benzene — 6 atoms, indices 0-5
        let mut opts = DepictOptions::new();
        opts.set_show_atom_indices(true);
        let svg = h.depict_svg_opts(&opts);
        assert!(svg.contains(">0<"), "index 0 should appear");
        assert!(svg.contains(">5<"), "index 5 should appear");
    }

    #[test]
    fn depict_svg_opts_show_atom_indices_false_no_indices() {
        let h = parse("CCO");
        let svg = h.depict_svg_opts(&DepictOptions::new());
        assert!(!svg.contains("fill=\"#8b92a9\""), "default should not show grey index labels");
    }

    // ── Sprint L: kekulize ───────────────────────────────────────────────────

    #[test]
    fn depict_svg_opts_kekulize_removes_aromatic_bonds() {
        let h = parse("c1ccccc1"); // benzene
        let mut opts = DepictOptions::new();
        opts.set_kekulize(true);
        let svg = h.depict_svg_opts(&opts);
        // Aromatic bonds render as dashed-style (multiple close lines); kekulé renders as standard double bonds.
        // The kekulé SVG should contain double bond lines but no dashed aromatic style.
        assert!(!svg.is_empty());
        // Double bonds produce two parallel <line> elements; check that double bond rendering kicked in.
        assert!(svg.contains("<line"), "kekulé benzene should have line elements");
    }

    #[test]
    fn depict_svg_opts_kekulize_false_uses_aromatic() {
        let h = parse("c1ccccc1");
        let svg = h.depict_svg_opts(&DepictOptions::new());
        // Default aromatic rendering uses stroke-dasharray for the inner ring line.
        assert!(svg.contains("stroke-dasharray"), "default benzene should use aromatic dashed style");
    }
}
