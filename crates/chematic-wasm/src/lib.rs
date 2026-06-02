//! `chematic-wasm` — WebAssembly bindings for the chematic cheminformatics library.
//!
//! Exposes a small, ergonomic API for parsing SMILES and computing molecular
//! descriptors from JavaScript/TypeScript via `wasm-bindgen`.

use wasm_bindgen::prelude::*;

#[wasm_bindgen(start)]
pub fn start() {
    console_error_panic_hook::set_once();
}

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
        chematic_depict::depict_svg_opts(&self.inner, &opts.to_render_options())
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

    // -----------------------------------------------------------------------
    // EState descriptors (Sprint P)
    // -----------------------------------------------------------------------

    /// Sum of EState indices over all heavy atoms.
    pub fn sum_estate(&self) -> f64 {
        chematic_chem::sum_estate(&self.inner)
    }

    /// Maximum EState index across all heavy atoms.
    pub fn max_estate(&self) -> f64 {
        chematic_chem::max_estate(&self.inner)
    }

    /// Minimum EState index across all heavy atoms.
    pub fn min_estate(&self) -> f64 {
        chematic_chem::min_estate(&self.inner)
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
    // R5: store RenderOptions directly; only the JS-incompatible HashSet fields are separate.
    inner: chematic_depict::RenderOptions,
    highlight_atoms: Vec<u32>,
    highlight_bonds: Vec<u32>,
    atom_color_entries: Vec<(u32, String)>,
}

#[wasm_bindgen]
impl DepictOptions {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self {
            inner: chematic_depict::RenderOptions::default(),
            highlight_atoms: vec![],
            highlight_bonds: vec![],
            atom_color_entries: vec![],
        }
    }

    pub fn set_width(&mut self, w: u32)           { self.inner.width = Some(w); }
    pub fn set_height(&mut self, h: u32)          { self.inner.height = Some(h); }
    pub fn set_padding(&mut self, p: f64)         { self.inner.padding = p; }
    pub fn set_background(&mut self, bg: String)  { self.inner.background = bg; }
    pub fn set_dark(&mut self, dark: bool)        { self.inner.dark = dark; }
    pub fn set_highlight_atoms(&mut self, atoms: Vec<u32>) { self.highlight_atoms = atoms; }
    pub fn set_highlight_bonds(&mut self, bonds: Vec<u32>) { self.highlight_bonds = bonds; }
    pub fn set_highlight_color(&mut self, color: String)   { self.inner.highlight_color = color; }
    /// Set a per-atom color override (CSS color string).  Calling multiple times
    /// for the same `idx` uses the last value.  The atom is highlighted even if
    /// not in `set_highlight_atoms`.
    pub fn set_atom_color(&mut self, idx: u32, color: String) {
        self.atom_color_entries.push((idx, color));
    }
    pub fn set_atom_ids(&mut self, v: bool)        { self.inner.atom_ids = v; }
    pub fn set_show_atom_indices(&mut self, v: bool) { self.inner.show_atom_indices = v; }
    pub fn set_kekulize(&mut self, v: bool)        { self.inner.kekulize = v; }

    pub(crate) fn to_render_options(&self) -> chematic_depict::RenderOptions {
        let mut ro = self.inner.clone();
        ro.highlight_atoms = self.highlight_atoms.iter()
            .map(|&i| chematic_core::AtomIdx(i))
            .collect();
        ro.highlight_bonds = self.highlight_bonds.iter()
            .map(|&i| chematic_core::BondIdx(i))
            .collect();
        ro.atom_color_map = self.atom_color_entries.iter()
            .map(|(i, c)| (chematic_core::AtomIdx(*i), c.clone()))
            .collect(); // last write wins for duplicate indices
        ro
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

/// Find all substructure matches of a SMARTS pattern in `mol`.
///
/// Returns JSON array of arrays of atom indices (sorted, 0-based).
/// Example: `[[0,1,2],[3,4,5]]` — two matches.
/// Returns `"[]"` if no match. Returns a JS error on invalid SMARTS.
#[wasm_bindgen]
pub fn smarts_match_atoms(smarts: &str, mol: &MolHandle) -> Result<String, JsValue> {
    let query = chematic_smarts::parse_smarts(smarts)
        .map_err(|e| JsValue::from_str(&format!("{e:?}")))?;
    let matches = chematic_smarts::find_matches(&query, &mol.inner);
    let parts: Vec<String> = matches
        .into_iter()
        .map(|m| {
            let mut idxs: Vec<u32> = m.values().map(|a| a.0).collect();
            idxs.sort_unstable();
            format!("[{}]", idxs.iter().map(|i| i.to_string()).collect::<Vec<_>>().join(","))
        })
        .collect();
    Ok(format!("[{}]", parts.join(",")))
}

/// Generate 3D coordinates for the molecule and return a PDB string.
///
/// Coordinates are generated using distance-geometry placement with ring templates.
/// Returns heavy-atom PDB (HETATM records, no explicit H).
#[wasm_bindgen]
pub fn generate_3d_pdb(mol: &MolHandle) -> String {
    let coords = chematic_3d::generate_coords(&mol.inner);
    chematic_3d::write_pdb(&mol.inner, &coords)
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

/// Render a reaction SMILES string (e.g. `"CC(=O)O.CCO>>CC(=O)OCC.O"`) as a
/// single SVG showing reactants → products with `+` separators.
///
/// Returns a self-contained SVG string.  Returns a JS error on invalid input.
#[wasm_bindgen]
pub fn depict_reaction_svg(rxn_smiles: &str) -> Result<String, JsValue> {
    let rxn = chematic_rxn::parse_reaction(rxn_smiles)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;

    const MOL_W: u32 = 200;
    const MOL_H: u32 = 180;
    const SEP_PLUS: u32 = 40;
    const SEP_ARROW: u32 = 60;
    const TOP_PAD: u32 = 20;

    let opts = chematic_depict::RenderOptions {
        width: Some(MOL_W),
        height: Some(MOL_H),
        ..Default::default()
    };

    let mut frags: Vec<(u32, String)> = Vec::new();
    let mut seps: Vec<(u32, &'static str)> = Vec::new();
    let mut cursor: u32 = 0;

    for (i, mol) in rxn.reactants.iter().enumerate() {
        if i > 0 {
            seps.push((cursor + SEP_PLUS / 2, "+"));
            cursor += SEP_PLUS;
        }
        frags.push((cursor, chematic_depict::depict_svg_opts(mol, &opts)));
        cursor += MOL_W;
    }

    seps.push((cursor + SEP_ARROW / 2, "→"));
    cursor += SEP_ARROW;

    for (i, mol) in rxn.products.iter().enumerate() {
        if i > 0 {
            seps.push((cursor + SEP_PLUS / 2, "+"));
            cursor += SEP_PLUS;
        }
        frags.push((cursor, chematic_depict::depict_svg_opts(mol, &opts)));
        cursor += MOL_W;
    }

    let total_w = cursor;
    let total_h = MOL_H + TOP_PAD;
    let mid_y = MOL_H / 2 + TOP_PAD;

    let mut out = format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{total_w}" height="{total_h}" viewBox="0 0 {total_w} {total_h}">"#
    );
    for (x, svg) in &frags {
        out.push_str(&svg.replacen("<svg ", &format!(r#"<svg x="{x}" y="{TOP_PAD}" "#), 1));
    }
    for (cx, sym) in &seps {
        out.push_str(&format!(
            r##"<text x="{cx}" y="{mid_y}" text-anchor="middle" dominant-baseline="central" font-size="20" font-family="sans-serif" fill="#555">{sym}</text>"##
        ));
    }
    out.push_str("</svg>");
    Ok(out)
}

// ---------------------------------------------------------------------------
// SDF / MOL I/O (Sprint P)
// ---------------------------------------------------------------------------

/// Parse a MOL V2000 block and return a `MolHandle`.
///
/// Returns a JS error string on parse failure.
#[wasm_bindgen]
pub fn mol_from_sdf_block(block: &str) -> Result<MolHandle, JsValue> {
    let (mol, _meta) = chematic_mol::parse_mol(block)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(MolHandle { inner: std::rc::Rc::new(mol) })
}

/// Serialize a molecule to a MOL V2000 block.
///
/// All atom coordinates are written as 0.0 (the `Molecule` type has no 2D
/// coordinate storage; real coordinates would require a separate layout pass).
#[wasm_bindgen]
pub fn to_mol_block(mol: &MolHandle) -> String {
    let meta = chematic_mol::MolMetadata::default();
    chematic_mol::write_mol(&mol.inner, &meta)
}

/// Parse an SDF string and return a JSON array of canonical SMILES strings.
///
/// Invalid records are represented as `null` in the array.
#[wasm_bindgen]
pub fn sdf_to_smiles_json(sdf: &str) -> String {
    let entries: Vec<String> = chematic_mol::SdfReader::new(sdf)
        .map(|r| match r {
            Ok((mol, _)) => {
                let smi = chematic_smiles::canonical_smiles(&mol);
                format!("\"{}\"", smi.replace('"', "\\\""))
            }
            Err(_) => "null".to_string(),
        })
        .collect();
    format!("[{}]", entries.join(","))
}

// ---------------------------------------------------------------------------
// EState free functions (Sprint P)
// ---------------------------------------------------------------------------

/// Per-atom EState values as a JSON array of f64.
///
/// Indices match `mol.atoms()` order.  Hydrogen atoms get 0.0.
#[wasm_bindgen]
pub fn estate_indices_json(mol: &MolHandle) -> String {
    let vals = chematic_chem::estate_indices(&mol.inner);
    let parts: Vec<String> = vals.iter().map(|v| format!("{v:.4}")).collect();
    format!("[{}]", parts.join(","))
}

// ---------------------------------------------------------------------------
// Tanimoto with topological path FP (Sprint P)
// ---------------------------------------------------------------------------

/// Tanimoto similarity between two molecules using topological path fingerprints.
#[wasm_bindgen]
pub fn tanimoto_topo_path(a: &MolHandle, b: &MolHandle) -> f64 {
    chematic_fp::tanimoto_topo_path(&a.inner, &b.inner)
}

// ---------------------------------------------------------------------------
// Sprint Q: IFG, VSA descriptors, Gasteiger charges, SA Score, Diversity
// ---------------------------------------------------------------------------

/// Identify functional groups. Returns a JSON array of objects:
/// `[{"atoms":[0,2,3],"type":"C,N,O"}, …]`
#[wasm_bindgen]
pub fn identify_functional_groups(mol: &MolHandle) -> String {
    let groups = chematic_chem::identify_functional_groups(&mol.inner);
    let parts: Vec<String> = groups.iter().map(|g| {
        let atoms: Vec<String> = g.atom_indices.iter().map(|i| i.to_string()).collect();
        format!("{{\"atoms\":[{}],\"type\":\"{}\"}}", atoms.join(","), g.atom_types)
    }).collect();
    format!("[{}]", parts.join(","))
}

/// Gasteiger-Marsili PEOE partial charges as a JSON array of f64.
#[wasm_bindgen]
pub fn gasteiger_charges_json(mol: &MolHandle) -> String {
    let q = chematic_chem::gasteiger_charges(&mol.inner);
    let parts: Vec<String> = q.iter().map(|v| format!("{v:.4}")).collect();
    format!("[{}]", parts.join(","))
}

/// Synthetic Accessibility Score (1 = easy, 10 = hard).
#[wasm_bindgen]
pub fn sa_score(mol: &MolHandle) -> f64 {
    chematic_chem::sa_score(&mol.inner)
}

/// SlogP_VSA descriptors (12 bins) as a JSON array.
#[wasm_bindgen]
pub fn slogp_vsa_json(mol: &MolHandle) -> String {
    let v = chematic_chem::slogp_vsa(&mol.inner);
    let parts: Vec<String> = v.iter().map(|x| format!("{x:.4}")).collect();
    format!("[{}]", parts.join(","))
}

/// SMR_VSA descriptors (10 bins) as a JSON array.
#[wasm_bindgen]
pub fn smr_vsa_json(mol: &MolHandle) -> String {
    let v = chematic_chem::smr_vsa(&mol.inner);
    let parts: Vec<String> = v.iter().map(|x| format!("{x:.4}")).collect();
    format!("[{}]", parts.join(","))
}

/// PEOE_VSA descriptors (14 bins) as a JSON array.
#[wasm_bindgen]
pub fn peoe_vsa_json(mol: &MolHandle) -> String {
    let v = chematic_chem::peoe_vsa(&mol.inner);
    let parts: Vec<String> = v.iter().map(|x| format!("{x:.4}")).collect();
    format!("[{}]", parts.join(","))
}

// ---------------------------------------------------------------------------
// Sprint T: named functional groups, atom info, per-atom color
// ---------------------------------------------------------------------------

/// Detect named functional groups in `mol`.
///
/// Returns a JSON array of `{"name":"hydroxyl","atoms":[3]}` objects.
/// Multiple matches of the same group (e.g. two hydroxyl groups) each appear
/// as a separate entry.  Overlapping groups (carboxylic acid → "carboxyl" +
/// "hydroxyl" + "carbonyl") are all returned.
#[wasm_bindgen]
pub fn detect_functional_groups(mol: &MolHandle) -> String {
    let groups = chematic_chem::detect_named_functional_groups(&mol.inner);
    let parts: Vec<String> = groups.iter().map(|g| {
        let atoms: Vec<String> = g.atoms.iter().map(|a| a.0.to_string()).collect();
        format!("{{\"name\":\"{}\",\"atoms\":[{}]}}", g.name, atoms.join(","))
    }).collect();
    format!("[{}]", parts.join(","))
}

/// Return information about a single atom as a JSON object.
///
/// `idx` is the 0-based atom index (matching `atoms()` order).
/// Returns `"null"` if `idx` is out of range.
///
/// Fields: `element` (symbol), `hybridization` ("sp"/"sp2"/"sp3"),
/// `charge` (formal charge integer), `isAromatic` (bool),
/// `totalHydrogens` (explicit + implicit H count, integer).
/// sp3d/sp3d2 (hypervalent P/S) are not distinguished from sp3/sp2.
#[wasm_bindgen]
pub fn get_atom_info(mol: &MolHandle, idx: u32) -> String {
    use chematic_core::{AtomIdx, BondOrder, implicit_hcount};
    let mol = &*mol.inner;
    if idx as usize >= mol.atom_count() {
        return "null".to_string();
    }
    let atom_idx = AtomIdx(idx);
    let atom = mol.atom(atom_idx);
    let symbol = atom.element.symbol();
    let charge = atom.charge;
    let is_aromatic = atom.aromatic;
    let total_h = atom.hydrogen_count.unwrap_or(0) as u32
        + implicit_hcount(mol, atom_idx) as u32;

    let hybridization = if is_aromatic {
        Some("sp2")
    } else {
        let mut double_count = 0u32;
        let mut has_triple = false;
        for (_, bond_idx) in mol.neighbors(atom_idx) {
            match mol.bond(bond_idx).order {
                BondOrder::Double => double_count += 1,
                BondOrder::Triple => { has_triple = true; break; }
                _ => {}
            }
        }
        if has_triple { Some("sp") }
        else if double_count >= 2 { Some("sp") }
        else if double_count == 1 { Some("sp2") }
        else { Some("sp3") }
    };

    let hyb_json = match hybridization {
        Some(h) => format!("\"{}\"", h),
        None    => "null".to_string(),
    };
    format!(
        "{{\"element\":\"{}\",\"hybridization\":{},\"charge\":{},\"isAromatic\":{},\"totalHydrogens\":{}}}",
        symbol, hyb_json, charge, is_aromatic, total_h
    )
}

// ---------------------------------------------------------------------------
// Sprint U: convenience SMILES-string-in functions + bond info API
// ---------------------------------------------------------------------------

/// Render a highlighted SVG from a SMILES string in one call.
///
/// `atoms` — 0-based atom indices to highlight (Uint32Array in JS).
/// `bonds` — 0-based bond indices to highlight (Uint32Array in JS).
/// `color` — CSS color for highlights (e.g. `"#ef4444"`); empty string uses default yellow.
///
/// Returns a JS error on SMILES parse failure.
#[wasm_bindgen]
pub fn smiles_to_svg_highlighted(
    smiles: &str,
    atoms: Vec<u32>,
    bonds: Vec<u32>,
    color: &str,
) -> Result<String, JsValue> {
    let mol = chematic_smiles::parse(smiles)
        .map(|m| MolHandle { inner: std::rc::Rc::new(m) })
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let mut opts = DepictOptions::new();
    opts.set_highlight_atoms(atoms);
    opts.set_highlight_bonds(bonds);
    if !color.is_empty() {
        opts.set_highlight_color(color.to_string());
    }
    Ok(mol.depict_svg_opts(&opts))
}

/// Find all SMARTS matches in a molecule given only SMILES strings.
///
/// Convenience wrapper around `smarts_match_atoms` that accepts raw SMILES
/// instead of a `MolHandle`.  Returns the same JSON format: `[[0,1],[3,4]]`.
/// Returns a JS error on SMILES or SMARTS parse failure.
#[wasm_bindgen]
pub fn match_smarts_smiles(smiles: &str, smarts: &str) -> Result<String, JsValue> {
    let mol = chematic_smiles::parse(smiles)
        .map(|m| MolHandle { inner: std::rc::Rc::new(m) })
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    smarts_match_atoms(smarts, &mol)
}

/// Tanimoto similarity between two molecules given only SMILES strings (ECFP4).
///
/// Returns a JS error on parse failure.
#[wasm_bindgen]
pub fn tanimoto_smiles(smiles1: &str, smiles2: &str) -> Result<f64, JsValue> {
    let m1 = chematic_smiles::parse(smiles1)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let m2 = chematic_smiles::parse(smiles2)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(chematic_fp::tanimoto_ecfp4(&m1, &m2))
}

/// Serialize a SMILES string directly to a MOL V2000 block.
///
/// Convenience wrapper; all atom coordinates are 0.0.
/// Returns a JS error on SMILES parse failure.
#[wasm_bindgen]
pub fn mol_block_from_smiles(smiles: &str) -> Result<String, JsValue> {
    let mol = chematic_smiles::parse(smiles)
        .map(|m| MolHandle { inner: std::rc::Rc::new(m) })
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(to_mol_block(&mol))
}

/// Return bond information as a JSON object, looked up by bond index.
///
/// `idx` is the 0-based bond index (order matches `mol.bonds()` iteration).
/// Returns `"null"` if `idx` is out of range.
///
/// Fields: `bondOrder` (1.0/1.5/2.0/3.0), `isAromatic` (bool),
/// `isInRing` (bool), `atomFrom` (u32), `atomTo` (u32).
#[wasm_bindgen]
pub fn get_bond_info(mol: &MolHandle, idx: u32) -> String {
    use chematic_core::BondOrder;
    let mol_ref = &*mol.inner;
    if idx as usize >= mol_ref.bond_count() {
        return "null".to_string();
    }
    let bond_idx = chematic_core::BondIdx(idx);
    let bond = mol_ref.bond(bond_idx);
    let bond_order: f64 = match bond.order {
        BondOrder::Aromatic => 1.5,
        other => other.order_value().unwrap_or(1.0) as f64,
    };
    let is_aromatic = bond.order == BondOrder::Aromatic;
    let in_ring = bond_in_ring(mol_ref, bond.atom1, bond.atom2);
    format!(
        "{{\"bondOrder\":{:.1},\"isAromatic\":{},\"isInRing\":{},\"atomFrom\":{},\"atomTo\":{}}}",
        bond_order, is_aromatic, in_ring, bond.atom1.0, bond.atom2.0
    )
}

/// Return bond information as a JSON object, looked up by the two bonded atom indices.
///
/// Useful when you know the atom indices from SMARTS matching or `data-atom-idx` SVG
/// attributes but not the bond index.  Returns `"null"` if no bond exists between them.
///
/// Fields: same as `get_bond_info` plus `bondIdx` (u32).
#[wasm_bindgen]
pub fn get_bond_between(mol: &MolHandle, atom1: u32, atom2: u32) -> String {
    use chematic_core::{AtomIdx, BondOrder};
    let mol_ref = &*mol.inner;
    let a = AtomIdx(atom1);
    let b = AtomIdx(atom2);
    if atom1 as usize >= mol_ref.atom_count() || atom2 as usize >= mol_ref.atom_count() {
        return "null".to_string();
    }
    let Some((bond_idx, bond)) = mol_ref.bond_between(a, b) else {
        return "null".to_string();
    };
    let bond_order: f64 = match bond.order {
        BondOrder::Aromatic => 1.5,
        other => other.order_value().unwrap_or(1.0) as f64,
    };
    let is_aromatic = bond.order == BondOrder::Aromatic;
    let in_ring = bond_in_ring(mol_ref, a, b);
    format!(
        "{{\"bondIdx\":{},\"bondOrder\":{:.1},\"isAromatic\":{},\"isInRing\":{},\"atomFrom\":{},\"atomTo\":{}}}",
        bond_idx.0, bond_order, is_aromatic, in_ring, atom1, atom2
    )
}

fn bond_in_ring(mol: &chematic_core::Molecule, a: chematic_core::AtomIdx, b: chematic_core::AtomIdx) -> bool {
    let rings = chematic_perception::find_sssr(mol);
    for ring in rings.rings() {
        let n = ring.len();
        for i in 0..n {
            if (ring[i] == a && ring[(i + 1) % n] == b) ||
               (ring[i] == b && ring[(i + 1) % n] == a) {
                return true;
            }
        }
    }
    false
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
        // Degree-0 O atoms use isolated (Hill) notation: H2O
        assert!(svg.matches("H2O").count() >= 2, "both water O atoms should render as H2O");
        assert!(!svg.is_empty());
    }

    // ── Sprint L: atom data attributes ──────────────────────────────────────

    #[test]
    fn depict_svg_opts_atom_ids_contains_data_attrs() {
        let h = parse("CC(=O)O"); // acetic acid: 2 C (unlabelled) + 2 O (labelled)
        let mut opts = DepictOptions::new();
        opts.set_atom_ids(true);
        let svg = h.depict_svg_opts(&opts);
        assert!(svg.contains("data-atom-idx="), "atom_ids should add data-atom-idx");
        assert!(svg.contains("data-element="), "atom_ids should add data-element");
        assert!(svg.contains("data-charge="), "atom_ids should add data-charge");
        // B3: all 4 atoms should be addressable (unlabelled carbons get invisible anchor)
        assert_eq!(svg.matches("data-atom-idx=").count(), 4, "all atoms should have data-atom-idx");
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
        assert!(!svg.is_empty());
        assert!(svg.contains("<line"), "kekulé benzene should have line elements");
        assert!(!svg.contains("stroke-dasharray"), "kekulé benzene must not use aromatic dashed style"); // B1
    }

    #[test]
    fn depict_svg_opts_kekulize_false_uses_aromatic() {
        let h = parse("c1ccccc1");
        let svg = h.depict_svg_opts(&DepictOptions::new());
        // Default aromatic rendering uses stroke-dasharray for the inner ring line.
        assert!(svg.contains("stroke-dasharray"), "default benzene should use aromatic dashed style");
    }

    // ── Sprint M: smarts_match_atoms ─────────────────────────────────────────

    #[test]
    fn smarts_match_benzene_ring_returns_json() {
        let mol = parse("c1ccccc1");
        let result = smarts_match_atoms("c1ccccc1", &mol);
        assert!(result.is_ok(), "valid SMARTS should not error");
        let json = result.unwrap();
        assert!(!json.is_empty() && json != "[]", "benzene ring SMARTS should find a match");
        assert!(json.starts_with("[["), "result should be array of arrays");
    }

    #[test]
    fn smarts_match_no_match_returns_empty_array() {
        let mol = parse("CC"); // ethane has no aromatic ring
        let result = smarts_match_atoms("c1ccccc1", &mol);
        assert!(result.is_ok(), "valid SMARTS on non-matching mol should not error");
        assert_eq!(result.unwrap(), "[]", "no match should return empty JSON array");
    }

    // ── Sprint N: generate_3d_pdb ────────────────────────────────────────────

    #[test]
    fn generate_3d_pdb_benzene_has_6_hetatm_lines() {
        let mol = parse("c1ccccc1");
        let pdb = generate_3d_pdb(&mol);
        let hetatm_count = pdb.lines().filter(|l| l.starts_with("HETATM")).count();
        assert_eq!(hetatm_count, 6, "benzene should produce 6 HETATM lines");
    }

    #[test]
    fn generate_3d_pdb_aspirin_no_nan() {
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O");
        let pdb = generate_3d_pdb(&mol);
        assert!(!pdb.contains("nan") && !pdb.contains("inf"), "PDB must have no NaN/Inf coords");
        assert!(pdb.contains("HETATM"), "must produce HETATM records");
    }

    // Note: smarts_match_atoms error-path test (invalid SMARTS) is omitted here
    // because JsValue::from_str panics outside a WASM runtime.
    // The underlying chematic_smarts::parse_smarts error path is tested separately:
    #[test]
    fn smarts_parse_invalid_is_err() {
        assert!(chematic_smarts::parse_smarts("[invalid").is_err(),
            "invalid SMARTS should return Err from parse_smarts");
    }

    // ── Sprint P: SDF I/O, EState, topo path FP ─────────────────────────────

    const ETHANE_MOL_BLOCK: &str = "\
ethane
  chematic

  2  1  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
";

    #[test]
    fn mol_from_sdf_block_ethane_atom_count() {
        let h = mol_from_sdf_block(ETHANE_MOL_BLOCK).expect("ethane parse");
        assert_eq!(h.atom_count(), 2);
    }

    #[test]
    fn sdf_to_smiles_json_two_records() {
        let sdf = format!("{ETHANE_MOL_BLOCK}$$$$\n{ETHANE_MOL_BLOCK}$$$$\n");
        let json = sdf_to_smiles_json(&sdf);
        assert!(json.starts_with('[') && json.ends_with(']'));
        // Should contain 2 SMILES entries separated by comma.
        let count = json.matches("CC").count();
        assert_eq!(count, 2, "expected 2 ethane SMILES in JSON, got: {json}");
    }

    #[test]
    fn estate_indices_json_acetic_acid_nonempty() {
        let h = parse("CC(=O)O");
        let json = estate_indices_json(&h);
        assert!(json.starts_with('[') && json.ends_with(']'));
        assert!(!json.is_empty() && json != "[]");
    }

    #[test]
    fn tanimoto_topo_path_same_mol_is_one() {
        let a = parse("c1ccccc1");
        let b = parse("c1ccccc1");
        assert!((tanimoto_topo_path(&a, &b) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn tanimoto_topo_path_different_mols_lt_one() {
        let a = parse("c1ccccc1");
        let b = parse("CC(=O)Oc1ccccc1C(=O)O");
        assert!(tanimoto_topo_path(&a, &b) < 1.0);
    }

    #[test]
    fn sum_estate_aspirin_positive() {
        let h = parse("CC(=O)Oc1ccccc1C(=O)O");
        assert!(h.sum_estate() > 0.0);
    }

    #[test]
    fn max_min_estate_ordering() {
        let h = parse("CC(=O)O");
        assert!(h.max_estate() >= h.min_estate());
    }

    // ── Sprint O: depict_reaction_svg ────────────────────────────────────────

    #[test]
    fn depict_reaction_svg_esterification() {
        let svg = depict_reaction_svg("CC(=O)O.CCO>>CC(=O)OCC.O").unwrap();
        assert!(svg.contains("→"), "must contain arrow character");
        assert!(svg.contains("<svg"), "must be valid SVG");
    }

    #[test]
    fn depict_reaction_svg_single_step() {
        let svg = depict_reaction_svg("C>>CC").unwrap();
        assert!(svg.contains("→"));
        assert!(svg.contains("<svg"));
    }

    // Error-path tested via the underlying parse_reaction (JsValue::from_str panics outside WASM).
    #[test]
    fn rxn_parse_missing_arrow_is_err() {
        assert!(chematic_rxn::parse_reaction("not_a_reaction").is_err());
    }

    // ── Sprint Q: IFG, Gasteiger, SA Score, VSA ──────────────────────────────

    #[test]
    fn identify_functional_groups_pyridine_has_n() {
        let h = parse("c1ccncc1");
        let json = identify_functional_groups(&h);
        assert!(json.contains('N'), "pyridine FG JSON should contain N: {json}");
        assert!(json.starts_with('['), "should be JSON array");
    }

    #[test]
    fn identify_functional_groups_hexane_empty() {
        let h = parse("CCCCCC");
        let json = identify_functional_groups(&h);
        assert_eq!(json, "[]", "hexane should have no functional groups");
    }

    #[test]
    fn gasteiger_charges_json_oxygen_negative() {
        let h = parse("CC(=O)O"); // acetic acid
        let json = gasteiger_charges_json(&h);
        assert!(json.starts_with('[') && json.ends_with(']'));
        // Parse values and verify at least one is negative.
        let has_negative = json.trim_matches(|c| c == '[' || c == ']')
            .split(',')
            .any(|s| s.trim().parse::<f64>().map(|v| v < 0.0).unwrap_or(false));
        assert!(has_negative, "acetic acid should have at least one negative charge: {json}");
    }

    #[test]
    fn sa_score_range() {
        let h = parse("CC(=O)Oc1ccccc1C(=O)O"); // aspirin
        let score = sa_score(&h);
        assert!(score >= 1.0 && score <= 10.0, "SA score out of [1,10]: {score:.2}");
    }

    #[test]
    fn slogp_vsa_json_length_12() {
        let h = parse("c1ccccc1");
        let json = slogp_vsa_json(&h);
        let count = json.trim_matches(|c| c == '[' || c == ']').split(',').count();
        assert_eq!(count, 12, "SlogP_VSA should have 12 values");
    }

    #[test]
    fn smr_vsa_json_length_10() {
        let h = parse("c1ccccc1");
        let json = smr_vsa_json(&h);
        let count = json.trim_matches(|c| c == '[' || c == ']').split(',').count();
        assert_eq!(count, 10, "SMR_VSA should have 10 values");
    }

    #[test]
    fn peoe_vsa_json_length_14() {
        let h = parse("c1ccccc1");
        let json = peoe_vsa_json(&h);
        let count = json.trim_matches(|c| c == '[' || c == ']').split(',').count();
        assert_eq!(count, 14, "PEOE_VSA should have 14 values");
    }
}
