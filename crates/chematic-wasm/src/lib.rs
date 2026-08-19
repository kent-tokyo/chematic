//! `chematic-wasm` — WebAssembly bindings for the chematic cheminformatics library.
//!
//! Exposes a small, ergonomic API for parsing SMILES and computing molecular
//! descriptors from JavaScript/TypeScript via `wasm-bindgen`.

use wasm_bindgen::prelude::*;

/// Maximum atom count for WASM to prevent DoS via unbounded JSON output.
const WASM_MAX_ATOMS: usize = 10_000;
/// Maximum input string length accepted by WASM convenience APIs.
const WASM_MAX_INPUT_BYTES: usize = 1_000_000;
/// Maximum number of SMILES records accepted by WASM batch APIs.
const WASM_MAX_BATCH_ITEMS: usize = 1_024;
/// Maximum length of one SMILES/name/property JSON array string element.
const WASM_MAX_JSON_STRING_BYTES: usize = 100_000;
/// Maximum SMARTS matches returned by WASM APIs.
const WASM_MAX_SMARTS_MATCHES: usize = 10_000;

mod format_io;
mod mol_3d;
mod mol_depict;
mod mol_descriptors;
mod mol_edit;
mod mol_fingerprints;
mod mol_io;
mod mol_reactions;
mod pipeline_v2;
#[cfg(test)]
mod tests;

pub use format_io::*;
pub use mol_3d::*;
pub use mol_depict::*;
pub use mol_descriptors::*;
pub use mol_edit::*;
pub use mol_fingerprints::*;
pub use mol_io::*;
pub use mol_reactions::*;
pub use pipeline_v2::embed_pipeline_v2_json;

#[wasm_bindgen(start)]
pub fn start() {
    std::panic::set_hook(Box::new(|info| {
        web_sys::console::error_1(&info.to_string().into());
    }));
}

// High-level workflow APIs
pub mod workflow;

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

    /// 2D PNG depiction — not available in the WASM build (PNG stack disabled to reduce bundle size).
    /// Use `depict_svg()` in browser contexts; rasterize client-side if needed.
    pub fn depict_png(&self) -> Vec<u8> {
        // tiny-skia is excluded from the WASM build to save ~200 KB.
        // Browsers can display SVG directly; client-side rasterization via Canvas API is preferred.
        Vec::new()
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
        let entries: Vec<String> = pairs.iter().map(|(k, v)| format!("\"{k}\": {v}")).collect();
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

    // -----------------------------------------------------------------------
    // Ring / stereo descriptor completions (Sprint W)
    // -----------------------------------------------------------------------

    /// Count of aliphatic (non-aromatic) rings in the SSSR.
    pub fn num_aliphatic_rings(&self) -> usize {
        chematic_chem::num_aliphatic_rings(&self.inner)
    }

    /// Count of fully saturated rings in the SSSR.
    pub fn num_saturated_rings(&self) -> usize {
        chematic_chem::num_saturated_rings(&self.inner)
    }

    /// Count of tetrahedral stereocenters with unspecified configuration.
    pub fn num_unspecified_stereocenters(&self) -> usize {
        chematic_chem::num_unspecified_stereocenters(&self.inner)
    }

    /// InChI string representation of the molecule.
    pub fn to_inchi(&self) -> String {
        chematic_inchi::inchi(&self.inner)
    }

    /// InChIKey (27-character identifier) for the molecule.
    pub fn to_inchikey(&self) -> String {
        let inchi_str = chematic_inchi::inchi(&self.inner);
        chematic_inchi::inchi_key(&inchi_str)
    }

    /// LogD (distribution coefficient) at a specific pH.
    ///
    /// Accounts for ionization state: neutral molecules return LogP unchanged,
    /// ionizable molecules are adjusted by log(neutral_fraction).
    pub fn logd_at_ph(&self, ph: f64) -> f64 {
        chematic_chem::logd_simple(&self.inner, ph)
    }

    /// Most acidic pKa in the molecule, or NaN if no acidic site.
    pub fn pka_acid_value(&self) -> f64 {
        chematic_chem::pka_acid(&self.inner).unwrap_or(f64::NAN)
    }

    /// Most basic pKa in the molecule, or NaN if no basic site.
    pub fn pka_base_value(&self) -> f64 {
        chematic_chem::pka_base(&self.inner).unwrap_or(f64::NAN)
    }

    /// Clark (2000) blood-brain barrier logBB score.
    /// logBB > −1.0 = likely CNS penetrant.
    pub fn bbb_score(&self) -> f64 {
        chematic_chem::bbb_score(&self.inner)
    }

    /// Returns true when TPSA < 90 Å², MW < 400, HBD ≤ 3.
    pub fn bbb_passes(&self) -> bool {
        chematic_chem::bbb_passes(&self.inner)
    }

    /// Palm (1997) Caco-2 intestinal permeability (logPCaco2).
    /// > −5.5 = high permeability.
    pub fn caco2_permeability(&self) -> f64 {
        chematic_chem::caco2_permeability(&self.inner)
    }

    /// hERG cardiac toxicity risk score (0.0–1.0).
    pub fn herg_risk_score(&self) -> f64 {
        chematic_chem::herg_risk_score(&self.inner)
    }

    /// CYP3A4 metabolic inhibition risk score (0.0–1.0).
    pub fn cyp3a4_inhibition_risk(&self) -> f64 {
        chematic_chem::cyp3a4_inhibition_risk(&self.inner)
    }

    /// LogD profile across a pH range as JSON.
    ///
    /// Returns `[{"ph":0.0,"logd":2.5}, ...]` with `steps` evenly-spaced pH points.
    pub fn logd_profile_json(&self, ph_start: f64, ph_end: f64, steps: usize) -> String {
        let profile = chematic_chem::logd_profile(&self.inner, ph_start, ph_end, steps);
        let items: Vec<String> = profile
            .iter()
            .map(|(ph, ld)| format!(r#"{{"ph":{ph:.2},"logd":{ld:.4}}}"#))
            .collect();
        format!("[{}]", items.join(","))
    }

    /// Isotope distribution as JSON.
    ///
    /// Returns `[{"mass":100.0,"abundance":0.9},...]` sorted by mass.
    /// `resolution`: m/z bin width in Da (e.g. `0.1` for nominal, `0.01` for high-res).
    pub fn isotope_distribution_json(&self, resolution: f64) -> String {
        let dist = chematic_chem::isotope_distribution(&self.inner, resolution);
        let items: Vec<String> = dist
            .iter()
            .map(|(mass, abund)| format!(r#"{{"mass":{mass:.4},"abundance":{abund:.6}}}"#))
            .collect();
        format!("[{}]", items.join(","))
    }

    /// Randić connectivity index (χ₀).
    ///
    /// χ₀ = Σ 1/√(d_i × d_j) over all bonds, where d is heavy-atom degree.
    pub fn randic_index(&self) -> f64 {
        chematic_chem::randic_index(&self.inner)
    }

    /// Zagreb index M1: Σ d_i² over all heavy atoms.
    pub fn zagreb_index_m1(&self) -> u32 {
        chematic_chem::zagreb_index_m1(&self.inner)
    }

    /// Generate IUPAC systematic name for the molecule.
    ///
    /// Returns the name string on success, or an empty string when the
    /// structure is outside the supported naming scope (complex polycyclics,
    /// multi-functional groups, etc.).
    pub fn iupac_name(&self) -> String {
        chematic_iupac::name(&self.inner).unwrap_or_default()
    }

    /// Assign CIP (R/S/E/Z) stereocenters and return JSON.
    ///
    /// Format: `{"centers":[{"atom":0,"code":"R"},{"atom":3,"code":"E"}]}`
    pub fn assign_cip_json(&self) -> String {
        use chematic_chem::assign_cip;
        use chematic_core::CipCode;
        let assignment = assign_cip(&self.inner);
        let centers: Vec<String> = assignment
            .assignments
            .iter()
            .map(|(idx, code)| {
                let code_str = match code {
                    CipCode::R => "R",
                    CipCode::S => "S",
                    CipCode::E => "E",
                    CipCode::Z => "Z",
                    CipCode::LowerR => "r",
                    CipCode::LowerS => "s",
                };
                format!(r#"{{"atom":{},"code":"{}"}}"#, idx.0, code_str)
            })
            .collect();
        format!(r#"{{"centers":[{}]}}"#, centers.join(","))
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

impl Default for DepictOptions {
    fn default() -> Self {
        Self::new()
    }
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

    pub fn set_width(&mut self, w: u32) {
        self.inner.width = Some(w);
    }
    pub fn set_height(&mut self, h: u32) {
        self.inner.height = Some(h);
    }
    pub fn set_padding(&mut self, p: f64) {
        self.inner.padding = p;
    }
    pub fn set_background(&mut self, bg: String) {
        self.inner.background = bg;
    }
    pub fn set_dark(&mut self, dark: bool) {
        self.inner.dark = dark;
    }
    pub fn set_highlight_atoms(&mut self, atoms: Vec<u32>) {
        self.highlight_atoms = atoms;
    }
    pub fn set_highlight_bonds(&mut self, bonds: Vec<u32>) {
        self.highlight_bonds = bonds;
    }
    pub fn set_highlight_color(&mut self, color: String) {
        self.inner.highlight_color = color;
    }
    /// Set a per-atom color override (CSS color string).  Calling multiple times
    /// for the same `idx` uses the last value.  The atom is highlighted even if
    /// not in `set_highlight_atoms`.
    pub fn set_atom_color(&mut self, idx: u32, color: String) {
        self.atom_color_entries.push((idx, color));
    }
    pub fn set_atom_ids(&mut self, v: bool) {
        self.inner.atom_ids = v;
    }
    pub fn set_show_atom_indices(&mut self, v: bool) {
        self.inner.show_atom_indices = v;
    }
    pub fn set_kekulize(&mut self, v: bool) {
        self.inner.kekulize = v;
    }

    pub(crate) fn to_render_options(&self) -> chematic_depict::RenderOptions {
        let mut ro = self.inner.clone();
        ro.highlight_atoms = self
            .highlight_atoms
            .iter()
            .map(|&i| chematic_core::AtomIdx(i))
            .collect();
        ro.highlight_bonds = self
            .highlight_bonds
            .iter()
            .map(|&i| chematic_core::BondIdx(i))
            .collect();
        ro.atom_color_map = self
            .atom_color_entries
            .iter()
            .map(|(i, c)| (chematic_core::AtomIdx(*i), c.clone()))
            .collect(); // last write wins for duplicate indices
        ro
    }
}

// ---------------------------------------------------------------------------
// Free functions exported to JS
// ---------------------------------------------------------------------------

fn bond_in_ring(
    mol: &chematic_core::Molecule,
    a: chematic_core::AtomIdx,
    b: chematic_core::AtomIdx,
) -> bool {
    let rings = chematic_perception::find_sssr(mol);
    for ring in rings.rings() {
        let n = ring.len();
        for i in 0..n {
            if (ring[i] == a && ring[(i + 1) % n] == b) || (ring[i] == b && ring[(i + 1) % n] == a)
            {
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
// Sprint V: scaffold, tautomers, standardization, MACCS, bulk descriptors
// ---------------------------------------------------------------------------

/// Escape a string for use as a JSON string value.
fn escape_json_string(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out
}

fn json_option_string_array(values: &[Option<String>]) -> String {
    let items = values
        .iter()
        .map(|value| match value {
            Some(s) => format!(r#""{}""#, escape_json_string(s)),
            None => "null".to_string(),
        })
        .collect::<Vec<_>>()
        .join(",");
    format!("[{items}]")
}

fn json_option_u8_array(values: &[Option<u8>]) -> String {
    let items = values
        .iter()
        .map(|value| value.map_or_else(|| "null".to_string(), |v| v.to_string()))
        .collect::<Vec<_>>()
        .join(",");
    format!("[{items}]")
}

fn enforce_wasm_input_len(label: &str, input: &str) -> Result<(), JsValue> {
    if input.len() > WASM_MAX_INPUT_BYTES {
        return Err(JsValue::from_str(&format!(
            "{label} exceeds maximum input size ({} > {} bytes)",
            input.len(),
            WASM_MAX_INPUT_BYTES
        )));
    }
    Ok(())
}

fn enforce_wasm_molecule_size(mol: &chematic_core::Molecule) -> Result<(), JsValue> {
    if mol.atom_count() > WASM_MAX_ATOMS {
        return Err(JsValue::from_str(&format!(
            "Molecule exceeds maximum atom count ({} > {})",
            mol.atom_count(),
            WASM_MAX_ATOMS
        )));
    }
    Ok(())
}

fn parse_wasm_string_json_array(json: &str, label: &str) -> Result<Vec<String>, JsValue> {
    enforce_wasm_input_len(label, json)?;
    let values: Vec<String> = serde_json::from_str(json)
        .map_err(|e| JsValue::from_str(&format!("{label} must be a JSON array of strings: {e}")))?;
    if values.len() > WASM_MAX_BATCH_ITEMS {
        return Err(JsValue::from_str(&format!(
            "{label} exceeds maximum item count ({} > {})",
            values.len(),
            WASM_MAX_BATCH_ITEMS
        )));
    }
    if let Some((idx, value)) = values
        .iter()
        .enumerate()
        .find(|(_, value)| value.len() > WASM_MAX_JSON_STRING_BYTES)
    {
        return Err(JsValue::from_str(&format!(
            "{label}[{idx}] exceeds maximum string size ({} > {} bytes)",
            value.len(),
            WASM_MAX_JSON_STRING_BYTES
        )));
    }
    Ok(values)
}

/// Internal helper: extract SMILES strings from a JSON array like `["CC","c1ccccc1"]`.
fn parse_smiles_json_array(json: &str) -> Result<Vec<String>, JsValue> {
    parse_wasm_string_json_array(json, "smiles_json")
}

// ---------------------------------------------------------------------------
// Sprint BB — BB-3: ConformerEnsemble WASM wrapper
// ---------------------------------------------------------------------------

/// BFS from `start`, excluding `core_atoms`, build a sub-molecule, return canonical SMILES.
fn rgroup_fragment_smiles(
    mol: &chematic_core::Molecule,
    start: chematic_core::AtomIdx,
    core_atoms: &std::collections::HashSet<chematic_core::AtomIdx>,
) -> String {
    use chematic_core::{Atom, AtomIdx, MoleculeBuilder};
    use std::collections::{HashMap, HashSet, VecDeque};

    // BFS to collect non-core atoms reachable from `start`.
    let mut fragment: Vec<AtomIdx> = Vec::new();
    let mut visited: HashSet<AtomIdx> = HashSet::new();
    let mut queue: VecDeque<AtomIdx> = VecDeque::new();
    queue.push_back(start);
    while let Some(idx) = queue.pop_front() {
        if visited.contains(&idx) || core_atoms.contains(&idx) {
            continue;
        }
        visited.insert(idx);
        fragment.push(idx);
        for (neighbor, _) in mol.neighbors(idx) {
            if !visited.contains(&neighbor) {
                queue.push_back(neighbor);
            }
        }
    }

    if fragment.is_empty() {
        return String::new();
    }

    let fragment_set: HashSet<AtomIdx> = fragment.iter().copied().collect();

    // Build sub-molecule.
    let mut builder = MoleculeBuilder::new();
    let mut idx_map: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    for &orig in &fragment {
        let atom = mol.atom(orig);
        let mut a = Atom::new(atom.element);
        a.charge = atom.charge;
        a.isotope = atom.isotope;
        a.aromatic = atom.aromatic;
        a.chirality = atom.chirality;
        a.hydrogen_count = atom.hydrogen_count;
        a.atom_map = atom.atom_map;
        let new_idx = builder.add_atom(a);
        idx_map.insert(orig, new_idx);
    }

    for (_, bond) in mol.bonds() {
        if fragment_set.contains(&bond.atom1)
            && fragment_set.contains(&bond.atom2)
            && let (Some(&n1), Some(&n2)) = (idx_map.get(&bond.atom1), idx_map.get(&bond.atom2))
        {
            let _ = builder.add_bond(n1, n2, bond.order);
        }
    }

    chematic_smiles::canonical_smiles(&builder.build())
}

// ---------------------------------------------------------------------------
// Sprint AA: FCFP bitvecs, Dice ECFP6, write_smiles, reaction normalization
// ---------------------------------------------------------------------------
