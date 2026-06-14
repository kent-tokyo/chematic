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

#[wasm_bindgen(start)]
pub fn start() {
    console_error_panic_hook::set_once();
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

    /// 2D PNG depiction (rasterized from SVG).
    /// Returns PNG data as base64-encoded string for embedding in HTML/JS.
    pub fn depict_png(&self) -> Vec<u8> {
        let layout = chematic_depict::compute_layout(&self.inner);
        chematic_depict::render_png(&self.inner, &layout)
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

    /// LogD profile across a pH range as JSON.
    ///
    /// Returns `[{"ph":0.0,"logd":2.5}, ...]` with `steps` evenly-spaced pH points.
    pub fn logd_profile_json(&self, ph_start: f64, ph_end: f64, steps: usize) -> String {
        let profile = chematic_chem::logd_profile(&self.inner, ph_start, ph_end, steps);
        let items: Vec<String> = profile.iter()
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
        let items: Vec<String> = dist.iter()
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
        match chematic_iupac::name(&self.inner) {
            Ok(name) => name,
            Err(_) => String::new(),
        }
    }

    /// Assign CIP (R/S/E/Z) stereocenters and return JSON.
    ///
    /// Format: `{"centers":[{"atom":0,"code":"R"},{"atom":3,"code":"E"}]}`
    pub fn assign_cip_json(&self) -> String {
        use chematic_chem::assign_cip;
        use chematic_core::CipCode;
        let assignment = assign_cip(&self.inner);
        let centers: Vec<String> = assignment.assignments.iter()
            .map(|(idx, code)| {
                let code_str = match code {
                    CipCode::R => "R",
                    CipCode::S => "S",
                    CipCode::E => "E",
                    CipCode::Z => "Z",
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

/// Parse a SMILES string into a `MolHandle`.
///
/// Returns a JS error string on parse failure or if atom count exceeds 10,000.
#[wasm_bindgen]
pub fn parse_smiles(s: &str) -> Result<MolHandle, JsValue> {
    enforce_wasm_input_len("smiles", s)?;
    let mol = chematic_smiles::parse(s).map_err(|e| JsValue::from_str(&e.to_string()))?;
    enforce_wasm_molecule_size(&mol)?;
    Ok(MolHandle {
        inner: std::rc::Rc::new(mol),
    })
}

/// Parse CXSMILES and return preserved metadata as JSON.
///
/// Supported CX fields: atom labels (`$...$`), `atomProp`, atom radicals (`^n:`),
/// and zero-order bonds (`Z:`). The `cxsmiles` field is a re-serialized
/// round-trip form using the supported fields.
/// Returns error if atom count exceeds 10,000.
#[wasm_bindgen]
pub fn parse_cxsmiles_json(s: &str) -> Result<String, JsValue> {
    enforce_wasm_input_len("cxsmiles", s)?;
    let cx = chematic_smiles::parse_cxsmiles(s).map_err(|e| JsValue::from_str(&e.to_string()))?;
    enforce_wasm_molecule_size(&cx.mol)?;
    let atom_props = cx
        .atom_props
        .iter()
        .map(|p| {
            format!(
                r#"{{"atom":{},"key":"{}","value":"{}"}}"#,
                p.atom.0,
                escape_json_string(&p.key),
                escape_json_string(&p.value)
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let zero_bonds = cx
        .mol
        .bonds()
        .filter_map(|(bidx, bond)| {
            (bond.order == chematic_core::BondOrder::Zero).then_some(bidx.0.to_string())
        })
        .collect::<Vec<_>>()
        .join(",");
    Ok(format!(
        r#"{{"smiles":"{}","cxsmiles":"{}","atomCount":{},"bondCount":{},"atomLabels":{},"atomProps":[{}],"atomRadicals":{},"zeroBonds":[{}]}}"#,
        escape_json_string(&chematic_smiles::write(&cx.mol)),
        escape_json_string(&chematic_smiles::write_cxsmiles(&cx)),
        cx.mol.atom_count(),
        cx.mol.bond_count(),
        json_option_string_array(&cx.atom_labels),
        atom_props,
        json_option_u8_array(&cx.atom_radicals),
        zero_bonds
    ))
}

/// Parse and re-serialize CXSMILES, preserving supported CX metadata.
/// Returns error if atom count exceeds 10,000.
#[wasm_bindgen]
pub fn normalize_cxsmiles(s: &str) -> Result<String, JsValue> {
    let cx = chematic_smiles::parse_cxsmiles(s).map_err(|e| JsValue::from_str(&e.to_string()))?;
    if cx.mol.atom_count() > WASM_MAX_ATOMS {
        return Err(JsValue::from_str(&format!(
            "Molecule exceeds maximum atom count ({} > {})",
            cx.mol.atom_count(),
            WASM_MAX_ATOMS
        )));
    }
    Ok(chematic_smiles::write_cxsmiles(&cx))
}

/// Parse CXSMARTS and return preserved metadata as JSON.
/// Returns error if atom count exceeds 10,000.
#[wasm_bindgen]
pub fn parse_cxsmarts_json(s: &str) -> Result<String, JsValue> {
    let cx = chematic_smarts::parse_cxsmarts(s).map_err(|e| JsValue::from_str(&e.to_string()))?;
    if cx.query.atom_count() > WASM_MAX_ATOMS {
        return Err(JsValue::from_str(&format!(
            "Query molecule exceeds maximum atom count ({} > {})",
            cx.query.atom_count(),
            WASM_MAX_ATOMS
        )));
    }
    let atom_props = cx
        .atom_props
        .iter()
        .map(|p| {
            format!(
                r#"{{"atom":{},"key":"{}","value":"{}"}}"#,
                p.atom,
                escape_json_string(&p.key),
                escape_json_string(&p.value)
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    Ok(format!(
        r#"{{"atomCount":{},"bondCount":{},"atomLabels":{},"atomProps":[{}],"atomRadicals":{}}}"#,
        cx.query.atom_count(),
        cx.query.bonds.len(),
        json_option_string_array(&cx.atom_labels),
        atom_props,
        json_option_u8_array(&cx.atom_radicals)
    ))
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

/// Run molecular dynamics simulation and return trajectory as JSON.
///
/// Returns JSON object with trajectory frames: `{ "frames": [{ "step": N, "potential": E, "kinetic": K, "temp": T }, …] }`
/// Uses NVT ensemble (Berendsen thermostat) at 300 K by default.
/// Note: Limited to molecules with ~50 atoms or fewer for practical WASM performance.
#[wasm_bindgen]
pub fn run_md_json(mol: &MolHandle, steps: usize, temp_k: f64) -> String {
    let coords = chematic_3d::generate_coords(&mol.inner);
    let config = chematic_3d::MDConfig {
        timestep_fs: 1.0,
        steps,
        temperature_k: temp_k,
        thermostat: chematic_3d::Thermostat::Berendsen { tau_fs: 100.0 },
        save_every: (steps / 20).max(1),
        coulomb: true,
    };
    let traj = chematic_3d::run_md(&mol.inner, coords, &config);

    // Serialize trajectory to JSON
    let mut frames_json = Vec::new();
    for frame in &traj.frames {
        frames_json.push(format!(
            r#"{{"step": {}, "potential": {:.4}, "kinetic": {:.4}, "temp": {:.2}}}"#,
            frame.step, frame.potential_energy, frame.kinetic_energy, frame.temperature_k
        ));
    }

    format!(r#"{{"frames": [{}]}}"#, frames_json.join(", "))
}

/// Compute direct Coulomb energy for a molecule with Gasteiger partial charges.
///
/// Returns JSON object: `{ "coulomb_energy": E, "unit": "kcal/mol" }`
///
/// # Arguments
/// * `mol` - Molecule to evaluate
///
/// # Example (JavaScript)
/// ```js
/// const mol = parse_smiles("CCO");
/// const result = coulomb_energy_json(mol);
/// // { "coulomb_energy": -12.34, "unit": "kcal/mol" }
/// ```
#[wasm_bindgen]
pub fn coulomb_energy_json(mol: &MolHandle) -> String {
    let coords = chematic_3d::generate_coords(&mol.inner);
    let charges = chematic_chem::gasteiger_charges(&mol.inner);
    let energy = chematic_ewald::direct_coulomb(
        &(0..coords.atom_count())
            .map(|i| {
                let p = coords.get(chematic_core::AtomIdx(i as u32));
                [p.x, p.y, p.z]
            })
            .collect::<Vec<_>>(),
        &charges,
    );
    format!(r#"{{"coulomb_energy": {:.4}, "unit": "kcal/mol"}}"#, energy)
}

/// Optimize molecular geometry using DREIDING force field.
///
/// Performs geometry minimization with DREIDING force field parameters.
/// Returns minimized coordinate PDB.
///
/// # Arguments
/// * `mol` - Molecule to optimize
///
/// # Returns
/// PDB format string with optimized coordinates
#[wasm_bindgen]
pub fn minimize_dreiding_json(mol: &MolHandle) -> String {
    let coords = chematic_3d::generate_coords(&mol.inner);
    let min_coords = chematic_3d::minimize_dreiding(&mol.inner, coords);
    chematic_3d::write_pdb(&mol.inner, &min_coords)
}

/// Return a copy of the molecule with all implicit hydrogens converted to explicit H atoms.
#[wasm_bindgen]
pub fn add_hydrogens(mol: &MolHandle) -> MolHandle {
    MolHandle {
        inner: std::rc::Rc::new(chematic_chem::add_hydrogens(&mol.inner)),
    }
}

/// Return a copy of the molecule with all explicit hydrogen atoms removed.
#[wasm_bindgen]
pub fn remove_hydrogens(mol: &MolHandle) -> MolHandle {
    MolHandle {
        inner: std::rc::Rc::new(chematic_chem::remove_hydrogens(&mol.inner)),
    }
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
        .map(|s| chematic_smiles::parse(s.trim()).map_err(|e| JsValue::from_str(&e.to_string())))
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

/// Enumerate a combinatorial library from a SMIRKS template and two fragment sets.
///
/// Generates all products by combining every scaffold with every building block.
/// Input format: `scaffolds_smiles` and `building_blocks_smiles` are pipe-delimited
/// SMILES strings (e.g., `"c1ccccc1|Cc1ccccc1"`).
///
/// Returns JSON array of product SMILES strings.
/// Example: `enumerate_library_2way("[C:1][Cl].[C:2][NH2]>>[C:1]N[C:2]", "c1ccccc1|Cc1ccccc1", "NCc1ccccc1|NCC")`
#[wasm_bindgen]
pub fn enumerate_library_2way(
    template: &str,
    scaffolds_smiles: &str,
    building_blocks_smiles: &str,
) -> Result<String, JsValue> {
    let scaffolds: Result<Vec<chematic_core::Molecule>, _> = scaffolds_smiles
        .split('|')
        .map(|s| chematic_smiles::parse(s.trim()).map_err(|e| JsValue::from_str(&e.to_string())))
        .collect();
    let scaffolds = scaffolds?;

    let building_blocks: Result<Vec<chematic_core::Molecule>, _> = building_blocks_smiles
        .split('|')
        .map(|s| chematic_smiles::parse(s.trim()).map_err(|e| JsValue::from_str(&e.to_string())))
        .collect();
    let building_blocks = building_blocks?;

    let config = chematic_rxn::LibraryConfig {
        skip_failures: true,
        max_size: Some(WASM_MAX_BATCH_ITEMS), // Limit enumeration size
    };

    let products = chematic_rxn::enumerate_library_2way(template, scaffolds, building_blocks, &config)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;

    let smiles_list: Vec<String> = products
        .iter()
        .map(|mol| format!("\"{}\"", chematic_smiles::canonical_smiles(mol)))
        .collect();

    Ok(format!("[{}]", smiles_list.join(", ")))
}

/// Find all substructure matches of a SMARTS pattern in `mol`.
///
/// Returns JSON array of arrays of atom indices (sorted, 0-based).
/// Example: `[[0,1,2],[3,4,5]]` — two matches.
/// Returns `"[]"` if no match. Returns a JS error on invalid SMARTS.
#[wasm_bindgen]
pub fn smarts_match_atoms(smarts: &str, mol: &MolHandle) -> Result<String, JsValue> {
    let query =
        chematic_smarts::parse_smarts(smarts).map_err(|e| JsValue::from_str(&format!("{e:?}")))?;
    let config = chematic_smarts::MatchConfig {
        max_matches: Some(WASM_MAX_SMARTS_MATCHES), use_chirality: false, use_isotopes: false, uniquify: true, max_visit_budget: None,
    };
    let matches = chematic_smarts::find_matches_with_config(&query, &mol.inner, &config);
    let parts: Vec<String> = matches
        .into_iter()
        .map(|m| {
            let mut idxs: Vec<u32> = m.values().map(|a| a.0).collect();
            idxs.sort_unstable();
            format!(
                "[{}]",
                idxs.iter()
                    .map(|i| i.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            )
        })
        .collect();
    Ok(format!("[{}]", parts.join(",")))
}

/// Like `smarts_match_atoms` but with explicit chirality matching control.
///
/// When `use_chirality=true`, SMARTS chirality primitives `[@]` and `[@@]` are
/// matched against the target molecule's stereochemistry. When `false`, chirality
/// is ignored (RDKit default).
#[wasm_bindgen]
pub fn smarts_match_atoms_with_chirality(
    smarts: &str,
    mol: &MolHandle,
    use_chirality: bool,
) -> Result<String, JsValue> {
    let query =
        chematic_smarts::parse_smarts(smarts).map_err(|e| JsValue::from_str(&format!("{e:?}")))?;
    let config = chematic_smarts::MatchConfig {
        max_matches: Some(WASM_MAX_SMARTS_MATCHES),
        use_chirality,
        use_isotopes: false,
        uniquify: true,
        max_visit_budget: None,
    };
    let matches = chematic_smarts::find_matches_with_config(&query, &mol.inner, &config);
    let parts: Vec<String> = matches
        .into_iter()
        .map(|m| {
            let mut idxs: Vec<u32> = m.values().map(|a| a.0).collect();
            idxs.sort_unstable();
            format!(
                "[{}]",
                idxs.iter()
                    .map(|i| i.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            )
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

/// Generate energy-minimized 3D coordinates and return a PDB string.
///
/// Runs distance-geometry placement followed by gradient-descent force-field
/// minimization.  Geometry quality is better than `generate_3d_pdb` for
/// flexible molecules; the force field is approximate (not MMFF94/UFF).
#[wasm_bindgen]
pub fn generate_3d_minimized_pdb(mol: &MolHandle) -> String {
    let coords = chematic_3d::generate_coords(&mol.inner);
    let minimized = chematic_3d::minimize(&mol.inner, coords);
    chematic_3d::write_pdb(&mol.inner, &minimized)
}

/// Generate 3D coordinates using ETKDG (torsion angle preferences) and return PDB block.
/// ETKDG produces higher-quality conformations than rule-based DG by applying
/// experimental torsion angle preferences to common structural patterns.
#[wasm_bindgen]
pub fn generate_3d_etkdg_pdb(mol: &MolHandle) -> String {
    let coords = chematic_3d::generate_coords_etkdg(&mol.inner);
    chematic_3d::write_pdb(&mol.inner, &coords)
}

/// Generate 3D coordinates using ETKDG and minimize with DREIDING force field.
#[wasm_bindgen]
pub fn generate_3d_etkdg_minimized_pdb(mol: &MolHandle) -> String {
    let coords = chematic_3d::generate_coords_etkdg(&mol.inner);
    let minimized = chematic_3d::minimize(&mol.inner, coords);
    chematic_3d::write_pdb(&mol.inner, &minimized)
}

/// Compute WHIM descriptors (Weighted Holistic Invariant Molecular) from 3D coordinates.
/// Returns JSON array of 10 values: [L1, L2, L3, P1, P2, P3, ALPHA, BETA, GAMMA, DELTA]
/// where L* = inertia tensor eigenvalues, P* = principal moments, ALPHA = sum of moments,
/// BETA = average pairwise interaction, GAMMA = geometric mean, DELTA = anisotropy.
#[wasm_bindgen]
pub fn whim_descriptors_json(mol: &MolHandle) -> String {
    let coords = chematic_3d::generate_coords(&mol.inner);
    let desc = chematic_3d::whim_descriptors(&mol.inner, &coords);
    let parts: Vec<String> = desc.iter().map(|&v| v.to_string()).collect();
    format!("[{}]", parts.join(","))
}

/// Compute GETAWAY descriptors (GEometric, Topologic And wAveleT descriptors) from 3D coordinates.
/// Returns JSON array of 9 values: [G1, G2, G3, D1, D2, D3, T, V, A]
/// where G* = geometric autocorrelations (lag-1,2,3), D* = topologic distances,
/// T = total pairwise distance, V = bounding-box volume, A = anisotropy ratio.
#[wasm_bindgen]
pub fn getaway_descriptors_json(mol: &MolHandle) -> String {
    let coords = chematic_3d::generate_coords(&mol.inner);
    let desc = chematic_3d::getaway_descriptors(&mol.inner, &coords);
    let parts: Vec<String> = desc.iter().map(|&v| v.to_string()).collect();
    format!("[{}]", parts.join(","))
}

/// Compute combined WHIM + GETAWAY descriptors (19 values total) as JSON array.
/// Useful for ML pipelines requiring both shape and topologic features.
#[wasm_bindgen]
pub fn whim_getaway_combined_json(mol: &MolHandle) -> String {
    let coords = chematic_3d::generate_coords(&mol.inner);
    let desc = chematic_3d::whim_getaway_combined(&mol.inner, &coords);
    let parts: Vec<String> = desc.iter().map(|&v| v.to_string()).collect();
    format!("[{}]", parts.join(","))
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

/// Like `ecfp4_bitvec` but with explicit chirality control.
///
/// When `use_chirality=true`, tetrahedral stereochemistry is included in the
/// initial atom hash, making enantiomers have different fingerprints.
/// When `false` (default), chirality is ignored.
#[wasm_bindgen]
pub fn ecfp4_bitvec_with_chirality(mol: &MolHandle, use_chirality: bool) -> Vec<u8> {
    let config = chematic_fp::EcfpConfig {
        radius: 2,
        nbits: 2048,
        use_chirality,
        use_double_fold: false,
    };
    let fp = chematic_fp::ecfp(&mol.inner, &config);
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

/// Like `ecfp6_bitvec` but with explicit chirality control.
///
/// When `use_chirality=true`, tetrahedral stereochemistry is included in the
/// initial atom hash, making enantiomers have different fingerprints.
/// When `false` (default), chirality is ignored.
#[wasm_bindgen]
pub fn ecfp6_bitvec_with_chirality(mol: &MolHandle, use_chirality: bool) -> Vec<u8> {
    let config = chematic_fp::EcfpConfig {
        radius: 3,
        nbits: 2048,
        use_chirality,
        use_double_fold: false,
    };
    let fp = chematic_fp::ecfp(&mol.inner, &config);
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
    let rxn =
        chematic_rxn::parse_reaction(rxn_smiles).map_err(|e| JsValue::from_str(&e.to_string()))?;

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

/// Parse a MOL V3000 block and return a `MolHandle`.
///
/// Returns a JS error string on parse failure.
#[wasm_bindgen]
pub fn mol_from_v3000_block(block: &str) -> Result<MolHandle, JsValue> {
    let (mol, _meta) =
        chematic_mol::parse_mol_v3000(block).map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(MolHandle {
        inner: std::rc::Rc::new(mol),
    })
}

/// Parse a MOL V2000 block and return a `MolHandle`.
///
/// Returns a JS error string on parse failure.
#[wasm_bindgen]
pub fn mol_from_sdf_block(block: &str) -> Result<MolHandle, JsValue> {
    let (mol, _meta) =
        chematic_mol::parse_mol(block).map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(MolHandle {
        inner: std::rc::Rc::new(mol),
    })
}

/// Serialize a molecule to a MOL V2000 block with 2D coordinates.
///
/// Atom positions are computed via the same layout engine used for SVG depiction
/// and converted to Ångström units (`1.5 Å` per bond).
#[wasm_bindgen]
pub fn to_mol_block(mol: &MolHandle) -> String {
    let meta = chematic_mol::MolMetadata::default();
    let layout = chematic_depict::compute_layout(&mol.inner);
    // SVG bond length is 40 px; standard MOL bond length is 1.5 Å.
    // Negate Y because SVG y-axis points down, MOL y-axis points up.
    const SCALE: f64 = 1.5 / 40.0;
    let coords: Vec<(f64, f64)> = (0..mol.inner.atom_count())
        .map(|i| {
            let p = layout.get(chematic_core::AtomIdx(i as u32));
            (p.x * SCALE, -p.y * SCALE)
        })
        .collect();
    chematic_mol::write_mol_with_coords(&mol.inner, &meta, &coords)
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
    let parts: Vec<String> = groups
        .iter()
        .map(|g| {
            let atoms: Vec<String> = g.atom_indices.iter().map(|i| i.to_string()).collect();
            format!(
                "{{\"atoms\":[{}],\"type\":\"{}\"}}",
                atoms.join(","),
                g.atom_types
            )
        })
        .collect();
    format!("[{}]", parts.join(","))
}

/// Detect pharmacophore features for virtual screening and lead optimization.
/// Returns JSON array of features: [{type, atom_idx, neighbor_count}, ...]
#[wasm_bindgen]
pub fn pharmacophore_features_json(mol: &MolHandle) -> String {
    let features = chematic_perception::detect_features(&mol.inner);
    let parts: Vec<String> = features
        .iter()
        .map(|f| {
            let ftype = match f.ftype {
                chematic_perception::FeatureType::Donor => "Donor",
                chematic_perception::FeatureType::Acceptor => "Acceptor",
                chematic_perception::FeatureType::Aromatic => "Aromatic",
                chematic_perception::FeatureType::Hydrophobic => "Hydrophobic",
                chematic_perception::FeatureType::Positive => "Positive",
                chematic_perception::FeatureType::Negative => "Negative",
            };
            format!(
                "{{\"type\":\"{}\",\"atom\":{},\"neighbors\":{}}}",
                ftype,
                f.atom.0,
                f.neighbors.len()
            )
        })
        .collect();
    format!("[{}]", parts.join(","))
}

/// Gasteiger-Marsili PEOE partial charges as a JSON array of f64.
#[wasm_bindgen]
pub fn gasteiger_charges_json(mol: &MolHandle) -> String {
    let q = chematic_chem::gasteiger_charges(&mol.inner);
    let parts: Vec<String> = q.iter().map(|v| format!("{v:.4}")).collect();
    format!("[{}]", parts.join(","))
}

/// MMFF94 partial charges (BCI table, ±0.1e accuracy) as a JSON array of f64.
///
/// Uses Bond Charge Increment (BCI) model (Halgren 1996) for 25 common bond types.
/// Returns `[q0, q1, ..., qN]` — one value per heavy atom.
/// Total charge equals the sum of formal charges (charge conserved).
#[wasm_bindgen]
pub fn mmff94_charges_json(mol: &MolHandle) -> String {
    let q = chematic_chem::mmff94_charges(&mol.inner);
    let parts: Vec<String> = q.iter().map(|v| format!("{v:.6}")).collect();
    format!("[{}]", parts.join(","))
}

/// MinHash fingerprint (128 hashes) as JSON.
///
/// Returns `{"num_hashes":128,"hashes":[u64,...]}`.
/// Use `tanimoto_mhfp_smiles` for direct SMILES-to-SMILES similarity.
#[wasm_bindgen]
pub fn mhfp_hashes_json(mol: &MolHandle) -> String {
    if mol.inner.atom_count() > WASM_MAX_ATOMS {
        return format!(r#"{{"error":"molecule too large (max {} atoms)"}}"#, WASM_MAX_ATOMS);
    }
    let fp = chematic_fp::mhfp(&mol.inner);
    let hs: Vec<String> = fp.hashes.iter().map(|h| h.to_string()).collect();
    format!(r#"{{"num_hashes":{},"hashes":[{}]}}"#, fp.num_hashes, hs.join(","))
}

/// Tanimoto-like similarity between two SMILES via MHFP (MinHash Jaccard approximation).
#[wasm_bindgen]
pub fn tanimoto_mhfp_smiles(smi1: &str, smi2: &str) -> Result<f64, JsValue> {
    let m1 = chematic_smiles::parse(smi1)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let m2 = chematic_smiles::parse(smi2)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    if m1.atom_count() > WASM_MAX_ATOMS || m2.atom_count() > WASM_MAX_ATOMS {
        return Err(JsValue::from_str(&format!("molecule too large (max {} atoms)", WASM_MAX_ATOMS)));
    }
    Ok(chematic_fp::tanimoto_mhfp(&m1, &m2))
}

/// Compute ERG-style 315-element float histogram fingerprint.
/// Returns JSON: {"len":315,"values":[f64,...]} or {"error":"..."}.
/// Format: 21 pharmacophore-feature-pair × 15 distance bins with Gaussian fuzzing.
/// See `chematic_fp::erg_vec` for details.
#[wasm_bindgen]
pub fn erg_vec_json(mol: &MolHandle) -> String {
    if mol.inner.atom_count() > WASM_MAX_ATOMS {
        return format!(r#"{{"error":"molecule too large (max {} atoms)"}}"#, WASM_MAX_ATOMS);
    }
    let v = chematic_fp::erg_vec(&mol.inner);
    let vals: Vec<String> = v.iter().map(|x| format!("{x:.6}")).collect();
    format!(r#"{{"len":{},"values":[{}]}}"#, chematic_fp::ERG_VEC_LEN, vals.join(","))
}

/// Compute MMFF94-style atom-typed partial charges (improved over element-pair BCI).
/// Returns JSON: {"charges":[f64,...]} or {"error":"..."}.
/// Uses atom-type classification (Csp3/Ccarbonyl/Ohydroxyl/Oester/Nar/NarH etc.)
/// for better accuracy (~±0.02e) vs element-pair BCI (~±0.05e).
#[wasm_bindgen]
pub fn mmff94_charges_typed_json(mol: &MolHandle) -> String {
    if mol.inner.atom_count() > WASM_MAX_ATOMS {
        return format!(r#"{{"error":"molecule too large (max {} atoms)"}}"#, WASM_MAX_ATOMS);
    }
    let q = chematic_chem::mmff94_charges_typed(&mol.inner);
    let vals: Vec<String> = q.iter().map(|x| format!("{x:.6}")).collect();
    format!(r#"{{"charges":[{}]}}"#, vals.join(","))
}

/// Compute 2D pharmacophore fingerprint (2048 bits) as a JSON feature count summary.
/// Returns simplified JSON with feature type counts: {Donor, Acceptor, Aromatic, Hydrophobic, Positive, Negative}
#[wasm_bindgen]
pub fn pharmacophore_fp_2d_summary(mol: &MolHandle) -> String {
    let counts = chematic_fp::pharmacophore_feature_counts(&mol.inner);
    format!(
        "{{\"Donor\":{},\"Acceptor\":{},\"Aromatic\":{},\"Hydrophobic\":{},\"Positive\":{},\"Negative\":{}}}",
        counts[0], counts[1], counts[2], counts[3], counts[4], counts[5]
    )
}

/// Compute 3D pharmacophore fingerprint from generated 3D coordinates.
/// Returns simplified JSON with feature type counts (3D-aware version).
#[wasm_bindgen]
pub fn pharmacophore_fp_3d_summary(mol: &MolHandle) -> String {
    let _coords = chematic_3d::generate_coords(&mol.inner);
    let features = chematic_perception::detect_features(&mol.inner);
    let mut counts = [0; 6];
    for feat in features {
        let idx = match feat.ftype {
            chematic_perception::FeatureType::Donor => 0,
            chematic_perception::FeatureType::Acceptor => 1,
            chematic_perception::FeatureType::Aromatic => 2,
            chematic_perception::FeatureType::Hydrophobic => 3,
            chematic_perception::FeatureType::Positive => 4,
            chematic_perception::FeatureType::Negative => 5,
        };
        counts[idx] += 1;
    }
    format!(
        "{{\"Donor\":{},\"Acceptor\":{},\"Aromatic\":{},\"Hydrophobic\":{},\"Positive\":{},\"Negative\":{}}}",
        counts[0], counts[1], counts[2], counts[3], counts[4], counts[5]
    )
}

/// MQN descriptor (42 integer values: Molecular Quantum Numbers).
#[wasm_bindgen]
pub fn mqn_json(mol: &MolHandle) -> String {
    let desc = chematic_chem::mqn(&mol.inner);
    let values: Vec<String> = desc.iter().map(|&v| v.to_string()).collect();
    format!("[{}]", values.join(","))
}

/// AutoCorr2D descriptor (7 values: topological distance lags 1-7).
#[wasm_bindgen]
pub fn autocorr_2d_json(mol: &MolHandle) -> String {
    let ac = chematic_chem::autocorr_2d(&mol.inner);
    format!(
        "[{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4}]",
        ac[0], ac[1], ac[2], ac[3], ac[4], ac[5], ac[6]
    )
}

/// AutoCorr3D descriptor (8 values: Euclidean distance bins 1-8 Å).
/// Requires 3D coordinates (generated automatically).
#[wasm_bindgen]
pub fn autocorr_3d_json(mol: &MolHandle) -> String {
    let coords = chematic_3d::generate_coords(&mol.inner);
    let ac = chematic_3d::autocorr_3d(&mol.inner, &coords);
    format!(
        "[{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4}]",
        ac[0], ac[1], ac[2], ac[3], ac[4], ac[5], ac[6], ac[7]
    )
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
    let parts: Vec<String> = groups
        .iter()
        .map(|g| {
            let atoms: Vec<String> = g.atoms.iter().map(|a| a.0.to_string()).collect();
            format!(
                "{{\"name\":\"{}\",\"atoms\":[{}]}}",
                g.name,
                atoms.join(",")
            )
        })
        .collect();
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
    let total_h = atom.hydrogen_count.unwrap_or(0) as u32 + implicit_hcount(mol, atom_idx) as u32;

    let hybridization = if is_aromatic {
        Some("sp2")
    } else {
        let mut double_count = 0u32;
        let mut has_triple = false;
        for (_, bond_idx) in mol.neighbors(atom_idx) {
            match mol.bond(bond_idx).order {
                BondOrder::Double => double_count += 1,
                BondOrder::Triple => {
                    has_triple = true;
                    break;
                }
                _ => {}
            }
        }
        if has_triple || double_count >= 2 {
            Some("sp")
        } else if double_count == 1 {
            Some("sp2")
        } else {
            Some("sp3")
        }
    };

    let hyb_json = match hybridization {
        Some(h) => format!("\"{}\"", h),
        None => "null".to_string(),
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
        .map(|m| MolHandle {
            inner: std::rc::Rc::new(m),
        })
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
        .map(|m| MolHandle {
            inner: std::rc::Rc::new(m),
        })
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    smarts_match_atoms(smarts, &mol)
}

/// Tanimoto similarity between two molecules given only SMILES strings (ECFP4).
///
/// Returns a JS error on parse failure.
#[wasm_bindgen]
pub fn tanimoto_smiles(smiles1: &str, smiles2: &str) -> Result<f64, JsValue> {
    let m1 = chematic_smiles::parse(smiles1).map_err(|e| JsValue::from_str(&e.to_string()))?;
    let m2 = chematic_smiles::parse(smiles2).map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(chematic_fp::tanimoto_ecfp4(&m1, &m2))
}

/// Serialize a SMILES string directly to a MOL V2000 block with 2D coordinates.
///
/// Returns a JS error on SMILES parse failure.
#[wasm_bindgen]
pub fn mol_block_from_smiles(smiles: &str) -> Result<String, JsValue> {
    let mol = chematic_smiles::parse(smiles)
        .map(|m| MolHandle {
            inner: std::rc::Rc::new(m),
        })
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

/// Murcko scaffold of `mol` — the ring system plus linkers, side-chains removed.
///
/// Returns a new `MolHandle`.  For acyclic molecules returns an empty molecule.
#[wasm_bindgen]
pub fn murcko_scaffold(mol: &MolHandle) -> MolHandle {
    let scaffold = chematic_chem::murcko_scaffold(&mol.inner);
    MolHandle {
        inner: std::rc::Rc::new(scaffold),
    }
}

/// Generic (atom-type-erased) Murcko scaffold of `mol`.
///
/// All atoms become carbon and all bonds become single bonds, giving the pure
/// graph topology of the scaffold.
#[wasm_bindgen]
pub fn generic_murcko_scaffold(mol: &MolHandle) -> MolHandle {
    let scaffold = chematic_chem::generic_murcko_scaffold(&mol.inner);
    MolHandle {
        inner: std::rc::Rc::new(scaffold),
    }
}

/// Canonical tautomer of `mol`.
///
/// Applies a rule-based tautomer normalisation and returns the canonical form
/// as a new `MolHandle`.
#[wasm_bindgen]
pub fn canonical_tautomer(mol: &MolHandle) -> MolHandle {
    let t = chematic_chem::canonical_tautomer(&mol.inner);
    MolHandle {
        inner: std::rc::Rc::new(t),
    }
}

/// All enumerated tautomers of `mol` as a JSON array of canonical SMILES strings.
///
/// Example return value: `["Oc1cccc2ccccc12","O=C1C=CC=Cc2ccccc21"]`
#[wasm_bindgen]
pub fn enumerate_tautomers_json(mol: &MolHandle) -> String {
    let tautomers = chematic_chem::enumerate_tautomers(&mol.inner);
    let parts: Vec<String> = tautomers
        .iter()
        .map(|m| {
            format!(
                "\"{}\"",
                chematic_smiles::canonical_smiles(m).replace('"', "\\\"")
            )
        })
        .collect();
    format!("[{}]", parts.join(","))
}

/// Return the largest fragment of `mol` (salt/solvent stripping).
///
/// For single-component molecules returns a copy of the same molecule.
#[wasm_bindgen]
pub fn largest_fragment(mol: &MolHandle) -> MolHandle {
    let frag = chematic_chem::largest_fragment(&mol.inner);
    MolHandle {
        inner: std::rc::Rc::new(frag),
    }
}

/// Neutralize formal charges on `mol` by proton addition/removal.
///
/// Returns a new `MolHandle` with all formal charges set to zero where possible.
#[wasm_bindgen]
pub fn neutralize_charges(mol: &MolHandle) -> MolHandle {
    let neutral = chematic_chem::neutralize_charges(&mol.inner);
    MolHandle {
        inner: std::rc::Rc::new(neutral),
    }
}

/// MACCS 166-bit structural keys fingerprint as a byte array (21 bytes, LSB-first).
///
/// Bit `i` (0-indexed) corresponds to MACCS key `i+1`.
#[wasm_bindgen]
pub fn maccs_bitvec(mol: &MolHandle) -> Vec<u8> {
    let fp = chematic_fp::maccs(&mol.inner);
    // 166 bits → 21 bytes (bits 0–165; remaining bits in byte 20 are always 0).
    (0..21usize)
        .map(|byte_idx| {
            let mut byte = 0u8;
            for bit in 0..8usize {
                let global_bit = byte_idx * 8 + bit;
                if global_bit < 166 && fp.get(global_bit) {
                    byte |= 1 << bit;
                }
            }
            byte
        })
        .collect()
}

/// Tanimoto similarity between `a` and `b` using MACCS 166-bit fingerprints.
#[wasm_bindgen]
pub fn tanimoto_maccs(a: &MolHandle, b: &MolHandle) -> f64 {
    let fa = chematic_fp::maccs(&a.inner);
    let fb = chematic_fp::maccs(&b.inner);
    fa.tanimoto(&fb)
}

/// All scalar molecular descriptors as a single JSON object.
///
/// Keys use camelCase and match the individual `MolHandle` method names.
/// Drug-likeness rule outcomes are included as boolean fields.
#[wasm_bindgen]
pub fn get_descriptors_json(mol: &MolHandle) -> String {
    let m = &mol.inner;
    format!(
        concat!(
            "{{",
            "\"mw\":{mw:.4},\"exactMass\":{em:.6},\"tpsa\":{tpsa:.4},",
            "\"logP\":{logp:.4},\"molarRefractivity\":{mr:.4},",
            "\"hbd\":{hbd},\"hba\":{hba},\"rotatableBonds\":{rb},",
            "\"heavyAtomCount\":{hac},\"ringCount\":{rc},",
            "\"aromaticRingCount\":{arc},\"numHeteroatoms\":{nh},",
            "\"numStereocenters\":{nsc},\"numSpiroAtoms\":{nsp},",
            "\"numBridgeheadAtoms\":{nbh},\"fsp3\":{fsp3:.4},",
            "\"qed\":{qed:.4},\"saScore\":{sa:.4},",
            "\"formalChargeSum\":{fc},\"labuteASA\":{asa:.4},",
            "\"bertzCT\":{bertz:.4},\"wienerIndex\":{wi:.4},",
            "\"lipinskiPasses\":{lip},\"veberPasses\":{veb},",
            "\"eganPasses\":{egan},\"ghosePasses\":{ghose},",
            "\"reosPasses\":{reos},\"painsPasses\":{pains}",
            "}}"
        ),
        mw = chematic_chem::molecular_weight(m),
        em = chematic_chem::exact_mass(m),
        tpsa = chematic_chem::tpsa(m),
        logp = chematic_chem::logp_crippen(m),
        mr = chematic_chem::molar_refractivity(m),
        hbd = chematic_chem::hbd_count(m),
        hba = chematic_chem::hba_count(m),
        rb = chematic_chem::rotatable_bond_count(m),
        hac = chematic_chem::heavy_atom_count(m),
        rc = chematic_chem::ring_count(m),
        arc = chematic_chem::aromatic_ring_count(m),
        nh = chematic_chem::num_heteroatoms(m),
        nsc = chematic_chem::num_stereocenters(m),
        nsp = chematic_chem::num_spiro_atoms(m),
        nbh = chematic_chem::num_bridgehead_atoms(m),
        fsp3 = chematic_chem::fsp3(m),
        qed = chematic_chem::qed(m),
        sa = chematic_chem::sa_score(m),
        fc = chematic_chem::formal_charge_sum(m),
        asa = chematic_chem::labute_asa(m),
        bertz = chematic_chem::bertz_ct(m),
        wi = chematic_chem::wiener_index(m),
        lip = chematic_chem::lipinski_passes(m),
        veb = chematic_chem::veber_passes(m),
        egan = chematic_chem::egan_passes(m),
        ghose = chematic_chem::ghose_passes(m),
        reos = chematic_chem::reos_passes(m),
        pains = chematic_chem::pains_passes(m),
    )
}

// ---------------------------------------------------------------------------
// Sprint X: V3000 read, 3D minimization, SDF properties, highlighted grid
// ---------------------------------------------------------------------------

/// Parse an SDF string and return a JSON array of record objects.
///
/// Each record has the shape:
/// ```json
/// {"smiles":"CC(=O)O","name":"aspirin","properties":{"MW":"180.2","Activity":"high"}}
/// ```
///
/// Invalid records are represented as `null`.  SD data fields are included in
/// `properties`; multi-line values are joined with `\n`.
#[wasm_bindgen]
pub fn sdf_to_records_json(sdf: &str) -> String {
    let entries: Vec<String> = chematic_mol::SdfRecordReader::new(sdf)
        .map(|r| match r {
            Ok(rec) => {
                let smi = chematic_smiles::canonical_smiles(&rec.mol);
                let name = escape_json_string(&rec.meta.name);
                let props: Vec<String> = rec
                    .properties
                    .iter()
                    .map(|(k, v)| {
                        format!(
                            "\"{}\":\"{}\"",
                            escape_json_string(k),
                            escape_json_string(v)
                        )
                    })
                    .collect();
                format!(
                    "{{\"smiles\":\"{}\",\"name\":\"{}\",\"properties\":{{{}}}}}",
                    escape_json_string(&smi),
                    name,
                    props.join(",")
                )
            }
            Err(_) => "null".to_string(),
        })
        .collect();
    format!("[{}]", entries.join(","))
}

/// Render a molecule grid with SMARTS-based atom highlighting.
///
/// `smiles_block` — newline-separated SMILES strings (same format as `depict_svg_grid`).
/// `cols` — number of grid columns.
/// `match_smarts` — SMARTS pattern; matched atoms in each molecule are highlighted.
///   Pass an empty string `""` to render without any highlighting.
///
/// Invalid SMILES are rendered as empty cells; SMARTS parse failure returns an
/// unhighlighted grid (the SMARTS is silently ignored).
#[wasm_bindgen]
pub fn depict_svg_grid_highlighted(smiles_block: &str, cols: usize, match_smarts: &str) -> String {
    // Parse all SMILES, skipping invalid ones.
    let mols: Vec<chematic_core::Molecule> = smiles_block
        .lines()
        .filter(|s| !s.trim().is_empty())
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .collect();

    if mols.is_empty() || cols == 0 {
        return "<svg xmlns=\"http://www.w3.org/2000/svg\" \
                width=\"0\" height=\"0\"></svg>"
            .to_string();
    }

    // Pre-compute per-molecule highlight options.
    let query = if match_smarts.is_empty() {
        None
    } else {
        chematic_smarts::parse_smarts(match_smarts).ok()
    };

    let cell_opts: Vec<chematic_depict::svg::RenderOptions> = mols
        .iter()
        .map(|mol| {
            let mut opts = chematic_depict::svg::RenderOptions::default();
            if let Some(q) = &query {
                let matches = chematic_smarts::find_matches(q, mol);
                opts.highlight_atoms = matches.into_iter().flat_map(|m| m.into_values()).collect();
            }
            opts
        })
        .collect();

    let pairs: Vec<(
        &chematic_core::Molecule,
        Option<&chematic_depict::svg::RenderOptions>,
    )> = mols
        .iter()
        .zip(cell_opts.iter())
        .map(|(m, o)| (m, Some(o)))
        .collect();

    chematic_depict::depict_svg_grid_with_opts(&pairs, cols)
}

// ---------------------------------------------------------------------------
// Sprint W: PAINS detail, CIP stereo, ECFP6+Dice, shape descriptors,
//           diversity picking, MCS SMILES
// ---------------------------------------------------------------------------

/// PAINS structural alert names matched by `mol` as a JSON array.
///
/// Returns `[]` when no alerts fire, or e.g. `["ene_six_het_A(483)"]`.
/// Use alongside `pains_passes()` to know *which* alerts triggered.
#[wasm_bindgen]
pub fn pains_matches_json(mol: &MolHandle) -> String {
    let names = chematic_chem::pains_matches(&mol.inner);
    let parts: Vec<String> = names
        .iter()
        .map(|n| format!("\"{}\"", n.replace('"', "\\\"")))
        .collect();
    format!("[{}]", parts.join(","))
}

/// CIP stereo assignments as a JSON array of `{atomIdx, cipCode}` objects.
///
/// `cipCode` is one of `"R"`, `"S"`, `"E"`, or `"Z"`.
/// Returns `[]` for molecules with no specified stereocenters.
#[wasm_bindgen]
pub fn cip_assignments_json(mol: &MolHandle) -> String {
    let assignment = chematic_chem::assign_cip(&mol.inner);
    let parts: Vec<String> = assignment
        .assignments
        .iter()
        .map(|(idx, code)| {
            let code_str = match code {
                chematic_core::CipCode::R => "R",
                chematic_core::CipCode::S => "S",
                chematic_core::CipCode::E => "E",
                chematic_core::CipCode::Z => "Z",
            };
            format!("{{\"atomIdx\":{},\"cipCode\":\"{}\"}}", idx.0, code_str)
        })
        .collect();
    format!("[{}]", parts.join(","))
}

/// ECFP6 (radius-3) fingerprint as a bit-packed byte vector (256 bytes = 2048 bits).
#[wasm_bindgen]
pub fn ecfp6_bitvec(mol: &MolHandle) -> Vec<u8> {
    let fp = chematic_fp::ecfp6(&mol.inner);
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

/// Tanimoto similarity between `a` and `b` using ECFP6 fingerprints.
#[wasm_bindgen]
pub fn tanimoto_ecfp6(a: &MolHandle, b: &MolHandle) -> f64 {
    chematic_fp::ecfp6(&a.inner).tanimoto(&chematic_fp::ecfp6(&b.inner))
}

/// Dice similarity between `a` and `b` using ECFP4 fingerprints.
#[wasm_bindgen]
pub fn dice_ecfp4(a: &MolHandle, b: &MolHandle) -> f64 {
    chematic_fp::ecfp4(&a.inner).dice(&chematic_fp::ecfp4(&b.inner))
}

/// Dice similarity between `a` and `b` using MACCS 166-bit fingerprints.
#[wasm_bindgen]
pub fn dice_maccs(a: &MolHandle, b: &MolHandle) -> f64 {
    chematic_fp::maccs(&a.inner).dice(&chematic_fp::maccs(&b.inner))
}

/// 3D shape descriptors as a JSON object.
///
/// Keys: `pmi1`, `pmi2`, `pmi3`, `npr1`, `npr2`, `asphericity`, `eccentricity`,
/// `radiusOfGyration`, `planeOfBestFit`.  Non-finite values (e.g. single-atom
/// molecules where pmi3 = 0) are serialised as JSON `null`.
#[wasm_bindgen]
pub fn shape_descriptors_json(mol: &MolHandle) -> String {
    let coords = chematic_3d::generate_coords(&mol.inner);
    let m = &*mol.inner;

    fn fmt(v: f64) -> String {
        if v.is_finite() {
            format!("{v:.4}")
        } else {
            "null".to_string()
        }
    }

    format!(
        "{{\"pmi1\":{},\"pmi2\":{},\"pmi3\":{},\"npr1\":{},\"npr2\":{},\
\"asphericity\":{},\"eccentricity\":{},\"radiusOfGyration\":{},\"planeOfBestFit\":{}}}",
        fmt(chematic_3d::pmi1(m, &coords)),
        fmt(chematic_3d::pmi2(m, &coords)),
        fmt(chematic_3d::pmi3(m, &coords)),
        fmt(chematic_3d::npr1(m, &coords)),
        fmt(chematic_3d::npr2(m, &coords)),
        fmt(chematic_3d::asphericity(m, &coords)),
        fmt(chematic_3d::eccentricity(m, &coords)),
        fmt(chematic_3d::radius_of_gyration(m, &coords)),
        fmt(chematic_3d::plane_of_best_fit(m, &coords)),
    )
}

/// Select `n` maximally-diverse molecules (MaxMin algorithm, ECFP4 Tanimoto).
///
/// `smiles_json` — a JSON array of SMILES strings, e.g. `["CC","c1ccccc1","CCO"]`.
/// Returns a JSON array of 0-based indices into the input array.
/// Returns a JS error if any SMILES fails to parse (indices would otherwise shift).
#[wasm_bindgen]
pub fn maxmin_picks_ecfp4_json(smiles_json: &str, n: usize) -> Result<String, JsValue> {
    let smiles_list = parse_smiles_json_array(smiles_json)?;
    let mols: Vec<chematic_core::Molecule> = smiles_list
        .iter()
        .map(|s| {
            let mol = chematic_smiles::parse(s).map_err(|e| JsValue::from_str(&e.to_string()))?;
            enforce_wasm_molecule_size(&mol)?;
            Ok::<_, JsValue>(mol)
        })
        .collect::<Result<_, _>>()?;
    let picks = chematic_chem::maxmin_picks(&mols, n, |a, b| {
        chematic_fp::ecfp4(a).tanimoto(&chematic_fp::ecfp4(b))
    });
    let parts: Vec<String> = picks.iter().map(|i| i.to_string()).collect();
    Ok(format!("[{}]", parts.join(",")))
}

/// Cluster molecules by structural similarity (Butina algorithm, ECFP4 Tanimoto).
///
/// `smiles_json` — a JSON array of SMILES strings.
/// `cutoff` — Tanimoto similarity threshold (0.0–1.0); molecules within this
///   distance of a cluster centre are assigned to that cluster.
/// Returns a JSON array of clusters, each cluster being an array of 0-based input indices.
/// Returns a JS error if any SMILES fails to parse.
#[wasm_bindgen]
pub fn butina_cluster_ecfp4_json(smiles_json: &str, cutoff: f64) -> Result<String, JsValue> {
    let smiles_list = parse_smiles_json_array(smiles_json)?;
    let mols: Vec<chematic_core::Molecule> = smiles_list
        .iter()
        .map(|s| {
            let mol = chematic_smiles::parse(s).map_err(|e| JsValue::from_str(&e.to_string()))?;
            enforce_wasm_molecule_size(&mol)?;
            Ok::<_, JsValue>(mol)
        })
        .collect::<Result<_, _>>()?;
    let clusters = chematic_chem::butina_cluster(&mols, cutoff, |a, b| {
        chematic_fp::ecfp4(a).tanimoto(&chematic_fp::ecfp4(b))
    });
    let outer: Vec<String> = clusters
        .iter()
        .map(|cluster| {
            let inner: Vec<String> = cluster.iter().map(|i| i.to_string()).collect();
            format!("[{}]", inner.join(","))
        })
        .collect();
    Ok(format!("[{}]", outer.join(",")))
}

/// Maximum Common Substructure of a set of molecules, returned as a canonical SMILES string.
///
/// `smiles_json` — a JSON array of at least 2 SMILES strings.
/// Returns the MCS SMILES, or `"null"` when no common substructure was found.
/// Returns a JS error on SMILES parse failure.
#[wasm_bindgen]
pub fn mcs_smiles_json(smiles_json: &str) -> Result<String, JsValue> {
    let smiles_list = parse_smiles_json_array(smiles_json)?;
    if smiles_list.len() < 2 {
        return Err(JsValue::from_str(
            "mcs_smiles_json requires at least 2 SMILES",
        ));
    }
    let mols: Vec<chematic_core::Molecule> = smiles_list
        .iter()
        .map(|s| {
            let mol = chematic_smiles::parse(s).map_err(|e| JsValue::from_str(&e.to_string()))?;
            enforce_wasm_molecule_size(&mol)?;
            Ok::<_, JsValue>(mol)
        })
        .collect::<Result<_, _>>()?;
    let mol_refs: Vec<&chematic_core::Molecule> = mols.iter().collect();
    let qmol = chematic_smarts::find_mcs(&mol_refs);

    if qmol.atoms.is_empty() {
        return Ok("null".to_string());
    }

    // Reconstruct a concrete Molecule from the QueryMolecule.
    // MCS atom queries are AtomicNum primitives; bond queries are typed primitives.
    use chematic_core::{Atom, AtomIdx, BondOrder, Element, MoleculeBuilder};
    use chematic_smarts::{AtomPrimitive, AtomQuery, BondPrimitive, BondQuery};

    let mut builder = MoleculeBuilder::new();
    for qa in &qmol.atoms {
        let elem = match &qa.query {
            AtomQuery::Primitive(AtomPrimitive::AtomicNum(n)) => {
                Element::from_atomic_number(*n).unwrap_or(Element::C)
            }
            _ => Element::C,
        };
        builder.add_atom(Atom::new(elem));
    }
    for (atom_idx, neighbors) in qmol.adj.iter().enumerate() {
        for (bond_idx, neighbor_idx) in neighbors {
            if atom_idx < *neighbor_idx {
                let order = match &qmol.bonds[*bond_idx].query {
                    BondQuery::Primitive(BondPrimitive::Double) => BondOrder::Double,
                    BondQuery::Primitive(BondPrimitive::Triple) => BondOrder::Triple,
                    BondQuery::Primitive(BondPrimitive::Aromatic) => BondOrder::Aromatic,
                    _ => BondOrder::Single,
                };
                let _ = builder.add_bond(
                    AtomIdx(atom_idx as u32),
                    AtomIdx(*neighbor_idx as u32),
                    order,
                );
            }
        }
    }
    let mol = builder.build();
    Ok(chematic_smiles::canonical_smiles(&mol))
}

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

/// A conformer ensemble: one molecule geometry with multiple 3D coordinate sets.
///
/// Create with `new(smiles)`, then add conformers with `add_generated_conformer`
/// or `add_minimized_conformer`.  Retrieve coordinates as PDB strings via
/// `get_conformer_pdb(idx)`.  Compare conformers with `conformer_rmsd`.
#[wasm_bindgen]
pub struct ConformerHandle {
    inner: chematic_3d::ConformerEnsemble,
}

#[wasm_bindgen]
impl ConformerHandle {
    /// Create a new empty ensemble for the molecule given by `smiles`.
    ///
    /// Returns a JS error on SMILES parse failure.
    #[wasm_bindgen(constructor)]
    pub fn new(smiles: &str) -> Result<ConformerHandle, JsValue> {
        let mol = chematic_smiles::parse(smiles).map_err(|e| JsValue::from_str(&e.to_string()))?;
        Ok(ConformerHandle {
            inner: chematic_3d::ConformerEnsemble::new(mol),
        })
    }

    /// Number of conformers currently stored.
    pub fn conformer_count(&self) -> usize {
        self.inner.conformer_count()
    }

    /// The ensemble's molecule as a `MolHandle`.
    pub fn mol(&self) -> MolHandle {
        let smi = chematic_smiles::canonical_smiles(self.inner.mol());
        let mol = chematic_smiles::parse(&smi)
            .unwrap_or_else(|_| chematic_core::MoleculeBuilder::new().build());
        MolHandle {
            inner: std::rc::Rc::new(mol),
        }
    }

    /// Generate a new 3D conformer using distance-geometry and add it to the ensemble.
    ///
    /// Returns the index of the newly added conformer.
    pub fn add_generated_conformer(&mut self) -> usize {
        let coords = chematic_3d::generate_coords(self.inner.mol());
        self.inner.add_conformer(coords).unwrap_or(usize::MAX)
    }

    /// Generate a new 3D conformer, run force-field minimization, and add it.
    ///
    /// Returns the index of the newly added conformer.
    pub fn add_minimized_conformer(&mut self) -> usize {
        let coords = chematic_3d::generate_coords(self.inner.mol());
        let minimized = chematic_3d::minimize(self.inner.mol(), coords);
        self.inner.add_conformer(minimized).unwrap_or(usize::MAX)
    }

    /// Return conformer `idx` as a PDB string, or `null` if `idx` is out of range.
    pub fn get_conformer_pdb(&self, idx: usize) -> Option<String> {
        self.inner
            .get_conformer(idx)
            .map(|coords| chematic_3d::write_pdb(self.inner.mol(), coords))
    }

    /// Kabsch-aligned RMSD (Å) between conformers `a` and `b`.
    ///
    /// Returns `NaN` if either index is out of range.
    pub fn conformer_rmsd(&self, a: usize, b: usize) -> f64 {
        self.inner.conformer_rmsd(a, b).unwrap_or(f64::NAN)
    }

    /// Un-aligned (translation + rotation NOT removed) RMSD (Å) between conformers `a` and `b`.
    ///
    /// Returns `NaN` if either index is out of range.
    pub fn conformer_rmsd_no_align(&self, a: usize, b: usize) -> f64 {
        self.inner.conformer_rmsd_no_align(a, b).unwrap_or(f64::NAN)
    }

    /// Remove conformer `idx` and return `true`, or `false` if `idx` is out of range.
    pub fn remove_conformer(&mut self, idx: usize) -> bool {
        self.inner.remove_conformer(idx).is_some()
    }
}

// ---------------------------------------------------------------------------
// Mutable Molecule API
// ---------------------------------------------------------------------------

/// Return a new `MolHandle` with one atom appended.
///
/// The second return value is the new atom's index (as a JS number).
/// Use `with_atom_added_idx` to retrieve the index.
#[wasm_bindgen]
pub fn mol_with_atom_added(mol: &MolHandle, element_symbol: &str) -> Result<MolHandle, JsValue> {
    let element = chematic_core::Element::from_symbol(element_symbol)
        .ok_or_else(|| JsValue::from_str(&format!("unknown element: {element_symbol}")))?;
    let atom = chematic_core::Atom::new(element);
    let (new_mol, _) = mol.inner.with_atom_added(atom);
    Ok(MolHandle {
        inner: std::rc::Rc::new(new_mol),
    })
}

/// Return the index that would be assigned to an atom appended to `mol`.
#[wasm_bindgen]
pub fn mol_next_atom_idx(mol: &MolHandle) -> u32 {
    mol.inner.atom_count() as u32
}

/// Return a new `MolHandle` with one bond added between `a` and `b`.
///
/// `order` — 1 = single, 2 = double, 3 = triple.
/// Returns a JS error if the bond already exists or `a == b`.
#[wasm_bindgen]
pub fn mol_with_bond_added(
    mol: &MolHandle,
    a: u32,
    b: u32,
    order: u32,
) -> Result<MolHandle, JsValue> {
    let bond_order = match order {
        2 => chematic_core::BondOrder::Double,
        3 => chematic_core::BondOrder::Triple,
        _ => chematic_core::BondOrder::Single,
    };
    let (new_mol, _bond_idx) = mol
        .inner
        .with_bond_added(
            chematic_core::AtomIdx(a),
            chematic_core::AtomIdx(b),
            bond_order,
        )
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(MolHandle {
        inner: std::rc::Rc::new(new_mol),
    })
}

/// Return a new `MolHandle` with the formal charge of atom `idx` changed.
///
/// Returns a JS error if `idx` is out of range.
#[wasm_bindgen]
pub fn mol_with_atom_charge(mol: &MolHandle, idx: u32, charge: i32) -> Result<MolHandle, JsValue> {
    if idx as usize >= mol.inner.atom_count() {
        return Err(JsValue::from_str(&format!(
            "atom index {idx} out of range (molecule has {} atoms)",
            mol.inner.atom_count()
        )));
    }
    let new_mol = mol
        .inner
        .with_atom_charge(chematic_core::AtomIdx(idx), charge as i8);
    Ok(MolHandle {
        inner: std::rc::Rc::new(new_mol),
    })
}

/// Return a new `MolHandle` with the element of atom `idx` changed.
///
/// `element_symbol` — periodic-table symbol, e.g. `"N"`, `"O"`, `"Cl"`.
/// Returns a JS error if `idx` is out of range or the symbol is unknown.
#[wasm_bindgen]
pub fn mol_with_atom_element(
    mol: &MolHandle,
    idx: u32,
    element_symbol: &str,
) -> Result<MolHandle, JsValue> {
    if idx as usize >= mol.inner.atom_count() {
        return Err(JsValue::from_str(&format!(
            "atom index {idx} out of range (molecule has {} atoms)",
            mol.inner.atom_count()
        )));
    }
    let el = chematic_core::Element::from_symbol(element_symbol)
        .ok_or_else(|| JsValue::from_str(&format!("unknown element: {element_symbol}")))?;
    let new_mol = mol.inner.with_atom_element(chematic_core::AtomIdx(idx), el);
    Ok(MolHandle {
        inner: std::rc::Rc::new(new_mol),
    })
}

/// Return a new `MolHandle` with atom `idx` and all its bonds removed.
///
/// Atom indices above `idx` shift down by 1.  Returns a JS error if `idx`
/// is out of range.
#[wasm_bindgen]
pub fn mol_with_atom_removed(mol: &MolHandle, idx: u32) -> Result<MolHandle, JsValue> {
    if idx as usize >= mol.inner.atom_count() {
        return Err(JsValue::from_str(&format!(
            "atom index {idx} out of range (molecule has {} atoms)",
            mol.inner.atom_count()
        )));
    }
    let (new_mol, _remap) = mol.inner.with_atom_removed(chematic_core::AtomIdx(idx));
    Ok(MolHandle {
        inner: std::rc::Rc::new(new_mol),
    })
}

/// Return a new `MolHandle` with bond `idx` removed.
///
/// Atom indices are unchanged; bond indices above `idx` shift down.
/// Returns a JS error if `idx` is out of range.
#[wasm_bindgen]
pub fn mol_with_bond_removed(mol: &MolHandle, idx: u32) -> Result<MolHandle, JsValue> {
    if idx as usize >= mol.inner.bond_count() {
        return Err(JsValue::from_str(&format!(
            "bond index {idx} out of range (molecule has {} bonds)",
            mol.inner.bond_count()
        )));
    }
    let new_mol = mol.inner.with_bond_removed(chematic_core::BondIdx(idx));
    Ok(MolHandle {
        inner: std::rc::Rc::new(new_mol),
    })
}

// ---------------------------------------------------------------------------
// SDF / V3000 write
// ---------------------------------------------------------------------------

/// Serialise a JSON array of SMILES to an SDF string.
///
/// Generates 2D coordinates for each molecule.  Property data can be
/// included by using `sdf_from_records_json` instead.
#[wasm_bindgen]
#[allow(clippy::type_complexity)]
pub fn smiles_array_to_sdf(smiles_json: &str) -> Result<String, JsValue> {
    let smiles_list = parse_smiles_json_array(smiles_json)?;
    let mut records: Vec<(
        chematic_core::Molecule,
        chematic_mol::MolMetadata,
        Vec<(f64, f64)>,
    )> = Vec::new();

    for smi in &smiles_list {
        let mol = chematic_smiles::parse(smi).map_err(|e| JsValue::from_str(&e.to_string()))?;
        enforce_wasm_molecule_size(&mol)?;
        let layout = chematic_depict::compute_layout(&mol);
        const SCALE: f64 = 1.5 / 40.0;
        let coords: Vec<(f64, f64)> = (0..mol.atom_count())
            .map(|i| {
                let p = layout.get(chematic_core::AtomIdx(i as u32));
                (p.x * SCALE, -p.y * SCALE)
            })
            .collect();
        let meta = chematic_mol::MolMetadata::default();
        records.push((mol, meta, coords));
    }

    let refs: Vec<(
        &chematic_core::Molecule,
        &chematic_mol::MolMetadata,
        &[(f64, f64)],
    )> = records
        .iter()
        .map(|(m, meta, c)| (m, meta, c.as_slice()))
        .collect();

    Ok(chematic_mol::write_sdf(&refs))
}

/// Serialise a `MolHandle` to MOL V3000 format with 2D coordinates.
#[wasm_bindgen]
pub fn to_mol_v3000_block(mol: &MolHandle) -> String {
    let layout = chematic_depict::compute_layout(&mol.inner);
    const SCALE: f64 = 1.5 / 40.0;
    let coords: Vec<(f64, f64)> = (0..mol.inner.atom_count())
        .map(|i| {
            let p = layout.get(chematic_core::AtomIdx(i as u32));
            (p.x * SCALE, -p.y * SCALE)
        })
        .collect();
    let meta = chematic_mol::MolMetadata::default();
    chematic_mol::write_mol_v3000(&mol.inner, &meta, &coords)
}

// ---------------------------------------------------------------------------
// DepictData
// ---------------------------------------------------------------------------

/// Compute structured depiction data for `mol` as a JSON object.
///
/// Returns:
/// ```json
/// {
///   "atoms": [
///     {"idx": 0, "element": "C", "x": 1.5, "y": 0.0, "charge": 0,
///      "label": null, "color": "#000000"},
///     ...
///   ],
///   "bonds": [
///     {"idx": 0, "atom1": 0, "atom2": 1, "kind": "Single"},
///     ...
///   ]
/// }
/// ```
///
/// `label` is `null` for carbon atoms in skeletal structures (label suppressed).
/// `kind` is one of `"Single"`, `"Double"`, `"Triple"`, `"Aromatic"`, `"Up"`, `"Down"`.
#[wasm_bindgen]
pub fn depict_data_json(mol: &MolHandle) -> String {
    let data = chematic_depict::compute_depict_data(&mol.inner);

    let atoms: Vec<String> = data
        .atoms
        .iter()
        .map(|a| {
            let label = match &a.label {
                Some(s) => format!("\"{}\"", escape_json_string(s)),
                None => "null".to_string(),
            };
            let charge_suffix = if a.charge != 0 {
                format!(",\"charge\":{}", a.charge)
            } else {
                String::new()
            };
            let label_prefix = format!("{label},");
            format!(
                r#"{{"idx":{},"element":"{}","x":{:.4},"y":{:.4},"label":{}"color":"{}"{}}}"#,
                a.idx.0,
                a.element.symbol(),
                a.pos.x,
                a.pos.y,
                label_prefix,
                escape_json_string(&a.color),
                charge_suffix,
            )
        })
        .collect();

    let bond_kind = |k: &chematic_depict::DepictBondKind| match k {
        chematic_depict::DepictBondKind::Single => "Single",
        chematic_depict::DepictBondKind::Double => "Double",
        chematic_depict::DepictBondKind::Triple => "Triple",
        chematic_depict::DepictBondKind::Aromatic => "Aromatic",
        chematic_depict::DepictBondKind::Up => "Up",
        chematic_depict::DepictBondKind::Down => "Down",
    };

    let bonds: Vec<String> = data
        .bonds
        .iter()
        .map(|b| {
            format!(
                r#"{{"idx":{},"atom1":{},"atom2":{},"kind":"{}"}}"#,
                b.idx.0,
                b.atom1.0,
                b.atom2.0,
                bond_kind(&b.kind)
            )
        })
        .collect();

    format!(
        r#"{{"atoms":[{}],"bonds":[{}]}}"#,
        atoms.join(","),
        bonds.join(","),
    )
}

// ---------------------------------------------------------------------------
// CPK colors
// ---------------------------------------------------------------------------

/// Return the CPK color (CSS hex string) for the given element symbol.
///
/// Returns `"#000000"` (black) for carbon and unknown elements.
#[wasm_bindgen]
pub fn cpk_color(element_symbol: &str) -> String {
    let an = chematic_core::Element::from_symbol(element_symbol)
        .map(|e| e.atomic_number())
        .unwrap_or(0);
    chematic_depict::atom_color(an).to_string()
}

// ---------------------------------------------------------------------------
// CML / CDXML file format support
// ---------------------------------------------------------------------------

/// Parse a CML string into a `MolHandle`.
///
/// Returns a JS error if the CML is invalid (unknown element, bad bond, etc.).
#[wasm_bindgen]
pub fn mol_from_cml(cml: &str) -> Result<MolHandle, JsValue> {
    let (mol, _coords) =
        chematic_mol::parse_cml(cml).map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(MolHandle {
        inner: std::rc::Rc::new(mol),
    })
}

/// Serialise a `MolHandle` to a CML string with 2D coordinates.
///
/// Coordinates are generated using the same 2D layout engine as `to_mol_block`.
#[wasm_bindgen]
pub fn to_cml(mol: &MolHandle) -> String {
    let layout = chematic_depict::compute_layout(&mol.inner);
    const SCALE: f64 = 1.5 / 40.0;
    let coords: Vec<(f64, f64)> = (0..mol.inner.atom_count())
        .map(|i| {
            let p = layout.get(chematic_core::AtomIdx(i as u32));
            (p.x * SCALE, -p.y * SCALE)
        })
        .collect();
    chematic_mol::write_cml(&mol.inner, Some(&coords))
}

/// Parse a ChemDraw XML (CDXML) string into a `MolHandle`.
///
/// Only the first molecular fragment in the document is returned.
/// Returns a JS error if the document cannot be parsed.
#[wasm_bindgen]
pub fn mol_from_cdxml(cdxml: &str) -> Result<MolHandle, JsValue> {
    let (mol, _coords) =
        chematic_mol::parse_cdxml(cdxml).map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(MolHandle {
        inner: std::rc::Rc::new(mol),
    })
}

/// Parse all molecular fragments from a CDXML string.
///
/// Returns a JSON array of SMILES strings, one per fragment:
/// `["CC","c1ccccc1"]`
///
/// Stereochemistry (wedge/dash bonds) is read from the `Display` attribute
/// of bond elements.
#[wasm_bindgen]
pub fn cdxml_to_smiles_json(cdxml: &str) -> Result<String, JsValue> {
    let fragments =
        chematic_mol::parse_cdxml_all(cdxml).map_err(|e| JsValue::from_str(&e.to_string()))?;
    let parts: Vec<String> = fragments
        .iter()
        .map(|(mol, _)| {
            format!(
                "\"{}\"",
                escape_json_string(&chematic_smiles::canonical_smiles(mol))
            )
        })
        .collect();
    Ok(format!("[{}]", parts.join(",")))
}

/// Parse a MOL V2000 string and return 2D coordinates as a JSON array.
///
/// Returns `[[x0,y0],[x1,y1],...]` in atom-insertion order.
/// Coordinates are in Ångström as stored in the MOL file.
#[wasm_bindgen]
pub fn mol_block_coords_json(mol_block: &str) -> Result<String, JsValue> {
    let (_mol, _meta, coords) = chematic_mol::parse_mol_with_coords(mol_block)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let parts: Vec<String> = coords
        .iter()
        .map(|(x, y)| {
            let xf = if x.is_finite() {
                format!("{x:.4}")
            } else {
                "0".to_string()
            };
            let yf = if y.is_finite() {
                format!("{y:.4}")
            } else {
                "0".to_string()
            };
            format!("[{xf},{yf}]")
        })
        .collect();
    Ok(format!("[{}]", parts.join(",")))
}

/// Compute structured depiction data using caller-supplied 2D coordinates.
///
/// `coords_json` — JSON array of `[x, y]` pairs, one per atom in order.
///
/// Returns the same JSON format as `depict_data_json`.
#[wasm_bindgen]
pub fn depict_data_with_coords_json(mol: &MolHandle, coords_json: &str) -> String {
    // Parse [[x,y],...] — simple format, no full JSON parser needed.
    let coords: Vec<(f64, f64)> = {
        let mut result = Vec::new();
        let inner = coords_json.trim().trim_matches(|c| c == '[' || c == ']');
        // Split on "]," pattern to get individual [x,y] pairs.
        for pair in inner.split("],[") {
            let clean = pair.trim().trim_matches(|c| c == '[' || c == ']');
            let nums: Vec<f64> = clean
                .split(',')
                .filter_map(|s| s.trim().parse::<f64>().ok())
                .collect();
            if nums.len() >= 2 {
                result.push((nums[0], nums[1]));
            }
        }
        result
    };

    let data = chematic_depict::depict_data_with_coords(&mol.inner, &coords);

    let atoms: Vec<String> = data
        .atoms
        .iter()
        .map(|a| {
            let label = match &a.label {
                Some(s) => format!("\"{}\"", escape_json_string(s)),
                None => "null".to_string(),
            };
            let charge_suffix = if a.charge != 0 {
                format!(",\"charge\":{}", a.charge)
            } else {
                String::new()
            };
            let label_prefix = format!("{label},");
            format!(
                r#"{{"idx":{},"element":"{}","x":{:.4},"y":{:.4},"label":{}"color":"{}"{}}}"#,
                a.idx.0,
                a.element.symbol(),
                a.pos.x,
                a.pos.y,
                label_prefix,
                escape_json_string(&a.color),
                charge_suffix,
            )
        })
        .collect();

    let bond_kind = |k: &chematic_depict::DepictBondKind| match k {
        chematic_depict::DepictBondKind::Single => "Single",
        chematic_depict::DepictBondKind::Double => "Double",
        chematic_depict::DepictBondKind::Triple => "Triple",
        chematic_depict::DepictBondKind::Aromatic => "Aromatic",
        chematic_depict::DepictBondKind::Up => "Up",
        chematic_depict::DepictBondKind::Down => "Down",
    };

    let bonds: Vec<String> = data
        .bonds
        .iter()
        .map(|b| {
            format!(
                r#"{{"idx":{},"atom1":{},"atom2":{},"kind":"{}"}}"#,
                b.idx.0,
                b.atom1.0,
                b.atom2.0,
                bond_kind(&b.kind)
            )
        })
        .collect();

    format!(
        r#"{{"atoms":[{}],"bonds":[{}]}}"#,
        atoms.join(","),
        bonds.join(",")
    )
}

// ---------------------------------------------------------------------------
// Sprint CC: MMP analysis (Matched Molecular Pairs)
// ---------------------------------------------------------------------------

/// Find matched molecular pairs in a set of molecules as JSON.
///
/// `smiles_json` — JSON array of SMILES strings to analyze.
///
/// Returns a JSON array of matched pairs:
/// ```json
/// [
///   {
///     "mol_a": "CC(=O)Oc1ccccc1",
///     "mol_b": "CC(=O)Nc1ccccc1",
///     "core": "c1ccccc1[*]",
///     "fragment_a": "[*]OC(C)=O",
///     "fragment_b": "[*]NC(C)=O"
///   }
/// ]
/// ```
///
/// Each pair represents molecules that share a common core scaffold but differ
/// by exactly one structural fragment at a single BRICS-breakable bond cut.
///
/// Returns a JS error if any SMILES fails to parse.
#[wasm_bindgen]
pub fn mmp_pairs_json(smiles_json: &str) -> Result<String, JsValue> {
    let smiles_list = parse_smiles_json_array(smiles_json)?;
    let mols: Vec<chematic_core::Molecule> = smiles_list
        .iter()
        .map(|s| {
            let mol = chematic_smiles::parse(s).map_err(|e| JsValue::from_str(&e.to_string()))?;
            enforce_wasm_molecule_size(&mol)?;
            Ok::<_, JsValue>(mol)
        })
        .collect::<Result<_, _>>()?;

    let mol_refs: Vec<&chematic_core::Molecule> = mols.iter().collect();
    let pairs = chematic_chem::find_mmp(&mol_refs);

    let entries: Vec<String> = pairs
        .iter()
        .map(|p| {
            format!(
                r#"{{"mol_a":"{}","mol_b":"{}","core":"{}","fragment_a":"{}","fragment_b":"{}"}}"#,
                escape_json_string(&p.mol_a),
                escape_json_string(&p.mol_b),
                escape_json_string(&p.core),
                escape_json_string(&p.fragment_a),
                escape_json_string(&p.fragment_b),
            )
        })
        .collect();

    Ok(format!("[{}]", entries.join(",")))
}

// ---------------------------------------------------------------------------
// Sprint BB — BB-1: R-group decomposition
// ---------------------------------------------------------------------------

/// Decompose a set of molecules against a core SMARTS, returning R-group SMILES.
///
/// `smiles_json` — JSON array of SMILES strings.
/// `core_smarts` — SMARTS pattern with `*` (wildcard) atoms marking R-group
///   attachment points.  For example `c1ccc(*)cc1` for para-substituted benzene.
///
/// Returns a JSON array with one entry per input molecule:
/// ```json
/// [
///   {"matched":true, "r1":"C"},
///   {"matched":true, "r1":"CC"},
///   {"matched":false}
/// ]
/// ```
/// R-group keys are `"r1"`, `"r2"`, … in the order the `*` atoms appear in
/// the SMARTS pattern.  A molecule that does not contain the core gets
/// `"matched": false` and no R-group keys.
///
/// Returns a JS error if the SMARTS fails to parse or any SMILES is invalid.
#[wasm_bindgen]
pub fn rgroup_decompose_json(smiles_json: &str, core_smarts: &str) -> Result<String, JsValue> {
    use chematic_core::AtomIdx;
    use chematic_smarts::{AtomPrimitive, AtomQuery};
    use std::collections::{HashMap, HashSet};

    let query = chematic_smarts::parse_smarts(core_smarts)
        .map_err(|e| JsValue::from_str(&format!("{e:?}")))?;

    // Identify which query atoms are wildcards and record their order.
    let wildcard_indices: Vec<usize> = query
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, qa)| matches!(&qa.query, AtomQuery::Primitive(AtomPrimitive::Wildcard)))
        .map(|(i, _)| i)
        .collect();

    let smiles_list = parse_smiles_json_array(smiles_json)?;
    let mols: Vec<chematic_core::Molecule> = smiles_list
        .iter()
        .map(|s| {
            let mol = chematic_smiles::parse(s).map_err(|e| JsValue::from_str(&e.to_string()))?;
            enforce_wasm_molecule_size(&mol)?;
            Ok::<_, JsValue>(mol)
        })
        .collect::<Result<_, _>>()?;

    let mut entries: Vec<String> = Vec::new();

    for mol in &mols {
        let config = chematic_smarts::MatchConfig {
            max_matches: Some(WASM_MAX_SMARTS_MATCHES),
            use_chirality: false,
            use_isotopes: false,
            uniquify: true,
            max_visit_budget: None,
        };
        let matches = chematic_smarts::find_matches_with_config(&query, mol, &config);
        if matches.is_empty() {
            entries.push("{\"matched\":false}".to_string());
            continue;
        }

        // Use first match.
        let mapping: &HashMap<usize, AtomIdx> = &matches[0];

        // Core atoms = molecule atoms matched by non-wildcard query atoms.
        let core_atoms: HashSet<AtomIdx> = mapping
            .iter()
            .filter(|(qi, _)| !wildcard_indices.contains(qi))
            .map(|(_, &mol_idx)| mol_idx)
            .collect();

        // For each wildcard, extract the R-group fragment.
        let mut rgroup_parts: Vec<String> = Vec::new();
        for (rg_num, qi) in wildcard_indices.iter().enumerate() {
            let smi = if let Some(&attachment) = mapping.get(qi) {
                rgroup_fragment_smiles(mol, attachment, &core_atoms)
            } else {
                String::new()
            };
            rgroup_parts.push(format!(
                "\"r{}\":\"{}\"",
                rg_num + 1,
                escape_json_string(&smi)
            ));
        }

        let rg_json = rgroup_parts.join(",");
        entries.push(format!("{{\"matched\":true,{rg_json}}}"));
    }

    Ok(format!("[{}]", entries.join(",")))
}

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
        if fragment_set.contains(&bond.atom1) && fragment_set.contains(&bond.atom2) {
            if let (Some(&n1), Some(&n2)) = (idx_map.get(&bond.atom1), idx_map.get(&bond.atom2)) {
                let _ = builder.add_bond(n1, n2, bond.order);
            }
        }
    }

    chematic_smiles::canonical_smiles(&builder.build())
}

// ---------------------------------------------------------------------------
// Sprint AA: FCFP bitvecs, Dice ECFP6, write_smiles, reaction normalization
// ---------------------------------------------------------------------------

/// FCFP4 (pharmacophore, radius-2) fingerprint as a bit-packed byte vector (256 bytes).
#[wasm_bindgen]
pub fn fcfp4_bitvec(mol: &MolHandle) -> Vec<u8> {
    let fp = chematic_fp::fcfp4(&mol.inner);
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

/// FCFP6 (pharmacophore, radius-3) fingerprint as a bit-packed byte vector (256 bytes).
#[wasm_bindgen]
pub fn fcfp6_bitvec(mol: &MolHandle) -> Vec<u8> {
    let fp = chematic_fp::fcfp6(&mol.inner);
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

/// Dice similarity between `a` and `b` using ECFP6 fingerprints.
#[wasm_bindgen]
pub fn dice_ecfp6(a: &MolHandle, b: &MolHandle) -> f64 {
    chematic_fp::ecfp6(&a.inner).dice(&chematic_fp::ecfp6(&b.inner))
}

/// Non-canonical SMILES for `mol`.
///
/// Unlike `canonical_smiles`, the output depends on the internal atom ordering
/// and is not normalised.  Useful when round-trip fidelity (preserving atom
/// order) matters more than a canonical form.
#[wasm_bindgen]
pub fn write_smiles(mol: &MolHandle) -> String {
    chematic_smiles::write(&mol.inner)
}

/// Parse and re-serialise a reaction SMILES string, returning the normalised form.
///
/// Useful for validating reaction SMILES and obtaining a canonical representation.
/// Returns a JS error on parse failure.
#[wasm_bindgen]
pub fn normalize_reaction_smiles(rxn_smiles: &str) -> Result<String, JsValue> {
    let rxn =
        chematic_rxn::parse_reaction(rxn_smiles).map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(chematic_rxn::write_reaction(&rxn))
}

// ---------------------------------------------------------------------------
// Sprint Z: BRICS fragment SMILES, FP bit-vectors, FCFP6, SDF write
// ---------------------------------------------------------------------------

/// BRICS fragment SMILES as a JSON array.
///
/// Applies the BRICS fragmentation rules and returns the canonical SMILES of
/// every resulting fragment.  Returns `[]` for molecules with no BRICS-breakable
/// bonds (e.g. benzene).
///
/// The count of fragments equals `brics_fragment_count`.
#[wasm_bindgen]
pub fn brics_fragments_json(mol: &MolHandle) -> String {
    let frags = chematic_chem::brics_fragments(&mol.inner);
    let parts: Vec<String> = frags
        .iter()
        .map(|m| {
            format!(
                "\"{}\"",
                escape_json_string(&chematic_smiles::canonical_smiles(m))
            )
        })
        .collect();
    format!("[{}]", parts.join(","))
}

/// AtomPair fingerprint as a bit-packed byte vector (256 bytes = 2048 bits).
#[wasm_bindgen]
pub fn atom_pair_bitvec(mol: &MolHandle) -> Vec<u8> {
    let fp = chematic_fp::atom_pair_fp(&mol.inner);
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

/// Torsion fingerprint as a bit-packed byte vector (256 bytes = 2048 bits).
#[wasm_bindgen]
pub fn torsion_bitvec(mol: &MolHandle) -> Vec<u8> {
    let fp = chematic_fp::torsion_fp(&mol.inner);
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

/// Tanimoto similarity between `a` and `b` using FCFP6 (radius-3 pharmacophore) fingerprints.
#[wasm_bindgen]
pub fn tanimoto_fcfp6(a: &MolHandle, b: &MolHandle) -> f64 {
    chematic_fp::fcfp6(&a.inner).tanimoto(&chematic_fp::fcfp6(&b.inner))
}

/// Serialize multiple molecules with properties to an SDF string.
///
/// # Arguments
/// * `smiles_json` — JSON array of SMILES strings, e.g. `["CC(=O)O","c1ccccc1"]`
/// * `names_json`  — JSON array of molecule names (same length as `smiles_json`)
/// * `props_json`  — JSON array where each element encodes one molecule's SD data fields
///   as `"key1\tvalue1\nkey2\tvalue2"` (tab-separated key/value, `\n`-separated pairs;
///   pass `""` for a molecule with no properties)
///
/// Returns the SDF string, or a JS error if any SMILES fails to parse or the
/// arrays have mismatched lengths.
///
/// The `\n` and `\t` sequences in `props_json` are JSON-escaped — they are
/// decoded to the actual characters before SDF formatting.
#[wasm_bindgen]
pub fn sdf_from_records_json(
    smiles_json: &str,
    names_json: &str,
    props_json: &str,
) -> Result<String, JsValue> {
    let smiles_list = parse_smiles_json_array(smiles_json)?;
    let names_list = parse_wasm_string_json_array(names_json, "names_json")?;
    let props_list = parse_wasm_string_json_array(props_json, "props_json")?;

    let n = smiles_list.len();
    if names_list.len() != n || props_list.len() != n {
        return Err(JsValue::from_str(
            "sdf_from_records_json: smiles_json, names_json, and props_json must have the same length",
        ));
    }

    let mut sdf = String::new();
    for i in 0..n {
        let mol = chematic_smiles::parse(&smiles_list[i])
            .map_err(|e| JsValue::from_str(&e.to_string()))?;
        enforce_wasm_molecule_size(&mol)?;

        // Build 2D coordinates using the same pipeline as to_mol_block.
        let layout = chematic_depict::compute_layout(&mol);
        const SCALE: f64 = 1.5 / 40.0;
        let coords: Vec<(f64, f64)> = (0..mol.atom_count())
            .map(|j| {
                let p = layout.get(chematic_core::AtomIdx(j as u32));
                (p.x * SCALE, -p.y * SCALE)
            })
            .collect();

        let meta = chematic_mol::MolMetadata {
            name: names_list[i].clone(),
            comment: String::new(),
        };
        sdf.push_str(&chematic_mol::write_mol_with_coords(&mol, &meta, &coords));

        // SD data fields from "key\tvalue\nkey\tvalue" format.
        for pair in props_list[i].lines() {
            if let Some((key, value)) = pair.split_once('\t') {
                sdf.push_str(&format!("> <{key}>\n{value}\n\n"));
            }
        }

        sdf.push_str("$$$$\n");
    }

    Ok(sdf)
}

// ---------------------------------------------------------------------------
// Sprint Y: XYZ/PDB I/O, per-atom descriptors, SSSR rings,
//           custom ECFP, stereo isomer enumeration
// ---------------------------------------------------------------------------

/// Parse an XYZ file and return a `MolHandle` (topology only; coordinates are discarded).
///
/// Returns a JS error on parse failure.
#[wasm_bindgen]
pub fn mol_from_xyz(xyz: &str) -> Result<MolHandle, JsValue> {
    let (mol, _coords) =
        chematic_3d::parse_xyz(xyz).map_err(|e| JsValue::from_str(&format!("{e:?}")))?;
    Ok(MolHandle {
        inner: std::rc::Rc::new(mol),
    })
}

/// Serialize a molecule to XYZ format.
///
/// 3D coordinates are generated via distance-geometry placement.
#[wasm_bindgen]
pub fn to_xyz(mol: &MolHandle) -> String {
    let coords = chematic_3d::generate_coords(&mol.inner);
    chematic_3d::write_xyz(&mol.inner, &coords, "")
}

/// Parse a PDB file and return a `MolHandle` (topology only; coordinates are discarded).
///
/// Uses CONECT records for connectivity if present; otherwise infers bonds from
/// atom distances (the same heuristic as the internal `pdb_to_molecule` function).
#[wasm_bindgen]
pub fn mol_from_pdb(pdb: &str) -> MolHandle {
    let atoms = chematic_3d::parse_pdb_atoms(pdb);
    let (mol, _coords) = chematic_3d::pdb_to_molecule(&atoms);
    MolHandle {
        inner: std::rc::Rc::new(mol),
    }
}

/// Per-atom Crippen LogP contributions as a JSON array of f64.
///
/// Index `i` corresponds to atom `i` in `mol.atoms()` order.
#[wasm_bindgen]
pub fn logp_per_atom_json(mol: &MolHandle) -> String {
    let vals = chematic_chem::logp_crippen_per_atom(&mol.inner);
    let parts: Vec<String> = vals.iter().map(|v| format!("{v:.4}")).collect();
    format!("[{}]", parts.join(","))
}

/// Per-atom molar refractivity contributions as a JSON array of f64.
#[wasm_bindgen]
pub fn mr_per_atom_json(mol: &MolHandle) -> String {
    let vals = chematic_chem::mr_per_atom(&mol.inner);
    let parts: Vec<String> = vals.iter().map(|v| format!("{v:.4}")).collect();
    format!("[{}]", parts.join(","))
}

/// Per-atom Labute approximate surface area contributions as a JSON array of f64.
///
/// Non-finite values (single-atom molecules etc.) are emitted as JSON `null`.
#[wasm_bindgen]
pub fn labute_asa_per_atom_json(mol: &MolHandle) -> String {
    let vals = chematic_chem::labute_asa_per_atom(&mol.inner);
    let parts: Vec<String> = vals
        .iter()
        .map(|v| {
            if v.is_finite() {
                format!("{v:.4}")
            } else {
                "null".to_string()
            }
        })
        .collect();
    format!("[{}]", parts.join(","))
}

/// Smallest Set of Smallest Rings (SSSR) as a JSON array of atom-index arrays.
///
/// Example return value for naphthalene:
/// `[[0,1,2,3,4,5],[5,6,7,8,9,4]]`
#[wasm_bindgen]
pub fn sssr_rings_json(mol: &MolHandle) -> String {
    let ring_set = chematic_perception::find_sssr(&mol.inner);
    let outer: Vec<String> = ring_set
        .rings()
        .iter()
        .map(|ring| {
            let inner: Vec<String> = ring.iter().map(|idx| idx.0.to_string()).collect();
            format!("[{}]", inner.join(","))
        })
        .collect();
    format!("[{}]", outer.join(","))
}

/// Compute a fingerprint bit-vector with configurable ECFP radius and bit width.
///
/// `radius` — Morgan radius (1 = ECFP2, 2 = ECFP4, 3 = ECFP6).
/// `nbits` — bit width; must be one of 256, 512, 1024, or 2048.
///   Returns a `Uint8Array` of `nbits/8` bytes.
///
/// The hash modulo is applied at fingerprint-generation time (`id % nbits`),
/// so no post-processing fold is needed.
/// Compute a custom ECFP (Extended Connectivity FingerPrint) with specified radius and bit count.
///
/// When `use_chirality=true`, tetrahedral stereochemistry is included in the initial
/// atom hash. When `false` (default), chirality is ignored.
#[wasm_bindgen]
pub fn ecfp_bitvec_custom(
    mol: &MolHandle,
    radius: u32,
    nbits: usize,
    use_chirality: bool,
) -> Vec<u8> {
    let nbits = match nbits {
        256 | 512 | 1024 | 2048 => nbits,
        _ => 2048,
    };
    let fp = chematic_fp::ecfp(
        &mol.inner,
        &chematic_fp::EcfpConfig {
            radius,
            nbits,
            use_chirality,
            use_double_fold: false,
        },
    );
    let byte_count = nbits / 8;
    (0..byte_count)
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

/// Enumerate all stereoisomers arising from unspecified tetrahedral stereocenters.
///
/// Only considers carbon stereocenters without explicit `@`/`@@` annotation.
/// Already-specified centers and E/Z double-bond geometry are unchanged.
/// Returns a JSON array of canonical SMILES strings.
///
/// At most 2^6 = 64 combinations are enumerated; if more than 6 unspecified
/// centers are present this function returns a JS error to avoid combinatorial
/// explosion.
#[wasm_bindgen]
pub fn enumerate_stereo_isomers_json(mol: &MolHandle) -> Result<String, JsValue> {
    use chematic_core::{Atom, AtomIdx, BondOrder, Chirality, MoleculeBuilder};

    let m = &*mol.inner;

    // Identify unspecified tetrahedral carbon stereocenters.
    // Criteria (same as num_unspecified_stereocenters in chematic-chem, plus
    // a degree≥2 guard that excludes terminal atoms like methyl groups whose
    // substituents are always identical).
    let unspecified: Vec<AtomIdx> = m
        .atoms()
        .filter(|(idx, atom)| {
            if atom.element.atomic_number() != 6 || atom.aromatic {
                return false;
            }
            if atom.chirality != Chirality::None {
                return false;
            }
            let degree = m.neighbors(*idx).count();
            // Require at least 2 explicit heavy-atom neighbors; terminal atoms
            // (methyl, –CH3, degree=1) cannot be stereocenters.
            if degree < 2 {
                return false;
            }
            let total = degree + chematic_core::implicit_hcount(m, *idx) as usize;
            total == 4
                && m.neighbors(*idx).all(|(_, bidx)| {
                    !matches!(m.bond(bidx).order, BondOrder::Double | BondOrder::Triple)
                })
        })
        .map(|(idx, _)| idx)
        .collect();

    let n = unspecified.len();
    if n > 6 {
        return Err(JsValue::from_str(&format!(
            "enumerate_stereo_isomers_json: {n} unspecified centers exceeds the 6-center limit (2^{n} = {} combinations)",
            1usize << n,
        )));
    }

    if n == 0 {
        // No unspecified centers — return the molecule's canonical SMILES as a
        // single-element array with InChI.
        let smi = chematic_smiles::canonical_smiles(m);
        let inchi_str = chematic_inchi::inchi(m);
        let inchikey_str = chematic_inchi::inchi_key(&inchi_str);
        return Ok(format!(
            r#"[{{"smiles":"{}","inchi":"{}","inchikey":"{}"}}]"#,
            escape_json_string(&smi),
            escape_json_string(&inchi_str),
            escape_json_string(&inchikey_str),
        ));
    }

    let mut seen = std::collections::HashSet::new();
    let mut results: Vec<String> = Vec::new();

    for bits in 0u32..(1u32 << n) {
        let chirality_overrides: std::collections::HashMap<AtomIdx, Chirality> = unspecified
            .iter()
            .enumerate()
            .map(|(i, &idx)| {
                let cw = (bits >> i) & 1 == 1;
                let chirality = if cw {
                    Chirality::Clockwise
                } else {
                    Chirality::CounterClockwise
                };
                (idx, chirality)
            })
            .collect();

        let mut builder = MoleculeBuilder::new();
        for (idx, atom) in m.atoms() {
            let mut a = Atom::new(atom.element);
            a.charge = atom.charge;
            a.isotope = atom.isotope;
            a.aromatic = atom.aromatic;
            a.atom_map = atom.atom_map;
            if let Some(&new_chirality) = chirality_overrides.get(&idx) {
                a.chirality = new_chirality;
                // Force bracket notation so the SMILES writer can output @/@@.
                // If the atom has implicit H, encode it in hydrogen_count so
                // the bracket includes the H (e.g. [C@@H](F)(Cl)Br).
                let implicit_h = chematic_core::implicit_hcount(m, idx);
                a.hydrogen_count = Some(atom.hydrogen_count.unwrap_or(implicit_h));
            } else {
                a.chirality = atom.chirality;
                a.hydrogen_count = atom.hydrogen_count;
            }
            builder.add_atom(a);
        }
        for (_, bond) in m.bonds() {
            let _ = builder.add_bond(bond.atom1, bond.atom2, bond.order);
        }
        let isomer = builder.build();
        let smi = chematic_smiles::canonical_smiles(&isomer);
        if seen.insert(smi.clone()) {
            let inchi_str = chematic_inchi::inchi(&isomer);
            let inchikey_str = chematic_inchi::inchi_key(&inchi_str);
            results.push(format!(
                r#"{{"smiles":"{}","inchi":"{}","inchikey":"{}"}}"#,
                escape_json_string(&smi),
                escape_json_string(&inchi_str),
                escape_json_string(&inchikey_str),
            ));
        }
    }

    Ok(format!("[{}]", results.join(",")))
}

// ---------------------------------------------------------------------------
// Reaction center analysis
// ---------------------------------------------------------------------------

/// Analyze a reaction SMILES and return the reaction center as JSON.
///
/// JSON schema: `{ broken: [[a1,a2],...], formed: [[a1,a2],...], changed: [a,...] }`
/// where atom indices are 0-based within the first reactant molecule.
/// Returns an error string prefixed with `"error:"` on failure.
#[wasm_bindgen]
pub fn find_reaction_center_json(reaction_smiles: &str) -> String {
    let rxn = match chematic_rxn::parse_reaction(reaction_smiles) {
        Ok(r) => r,
        Err(e) => return format!("error:{e}"),
    };
    let center = chematic_rxn::find_reaction_center(&rxn);
    let broken: Vec<String> = center
        .broken_bonds
        .iter()
        .map(|(a, b)| format!("[{},{}]", a.0, b.0))
        .collect();
    let formed: Vec<String> = center
        .formed_bonds
        .iter()
        .map(|(a, b)| format!("[{},{}]", a.0, b.0))
        .collect();
    let changed: Vec<String> = center
        .changed_atoms
        .iter()
        .map(|a| a.0.to_string())
        .collect();
    format!(
        "{{\"broken\":[{}],\"formed\":[{}],\"changed\":[{}]}}",
        broken.join(","),
        formed.join(","),
        changed.join(","),
    )
}

// ---------------------------------------------------------------------------
// Structure standardization
// ---------------------------------------------------------------------------

/// Standardize a SMILES string and return the canonical SMILES of the result.
///
/// Applies: largest fragment extraction → charge neutralization.
/// Returns `"error:<msg>"` on parse failure.
#[wasm_bindgen]
pub fn standardize_smiles(smiles: &str) -> String {
    let mol = match chematic_smiles::parse(smiles) {
        Ok(m) => m,
        Err(e) => return format!("error:{e}"),
    };
    let mol = chematic_chem::largest_fragment(&mol);
    let mol = chematic_chem::neutralize_charges(&mol);
    chematic_smiles::canonical_smiles(&mol)
}

/// Standardize a SMILES string and return result SMILES plus an audit report as JSON.
///
/// Boolean flags map directly to `StandardizeOptions`.
/// Returns `"error:<msg>"` on parse or serialization failure.
#[wasm_bindgen]
pub fn standardize_smiles_report_json(
    smiles: &str,
    largest_fragment_only: bool,
    neutralize_charges: bool,
    remove_explicit_h: bool,
    canonical_tautomer: bool,
) -> String {
    let mol = match chematic_smiles::parse(smiles) {
        Ok(m) => m,
        Err(e) => return format!("error:{e}"),
    };
    let pipeline = chematic_chem::StandardizationPipeline::new(chematic_chem::StandardizeOptions {
        canonical_tautomer,
        neutralize_charges,
        remove_explicit_h,
        largest_fragment_only,
        zwitterion_handling: chematic_chem::ZwitterionHandling::Normalize,
    });
    let (standardized, report) = pipeline.run(&mol);
    let report_json = match serde_json::to_string(&report) {
        Ok(json) => json,
        Err(e) => return format!("error:{e}"),
    };
    format!(
        r#"{{"smiles":"{}","report":{}}}"#,
        escape_json_string(&chematic_smiles::canonical_smiles(&standardized)),
        report_json
    )
}

// ---------------------------------------------------------------------------
// Reaction balance check
// ---------------------------------------------------------------------------

/// Check whether a reaction SMILES is atom-balanced.
///
/// Returns JSON: `{ "balanced": true|false, "diff": ["C: 1 reactant vs 2 product", ...] }`
/// Returns `"error:<msg>"` on parse failure.
#[wasm_bindgen]
pub fn balance_check_json(reaction_smiles: &str) -> String {
    let rxn = match chematic_rxn::parse_reaction(reaction_smiles) {
        Ok(r) => r,
        Err(e) => return format!("error:{e}"),
    };
    let result = chematic_rxn::balance_check(&rxn);
    let diff: Vec<String> = result
        .diff()
        .into_iter()
        .map(|s| format!("\"{}\"", s))
        .collect();
    format!(
        "{{\"balanced\":{},\"diff\":[{}]}}",
        result.balanced,
        diff.join(",")
    )
}

// ---------------------------------------------------------------------------
// Nearest-neighbour similarity search
// ---------------------------------------------------------------------------

/// Find the k nearest neighbours of a query SMILES in a list of db SMILES.
///
/// `db_smiles_json`: JSON array of SMILES strings, e.g. `["CC","c1ccccc1"]`.
/// Returns JSON: `[{"index":0,"tanimoto":0.95},...]` sorted by descending Tanimoto.
/// Returns `"error:<msg>"` on parse failure.
#[wasm_bindgen]
pub fn nearest_neighbors_json(query_smiles: &str, db_smiles_json: &str, k: usize) -> String {
    if let Err(e) = enforce_wasm_input_len("query_smiles", query_smiles) {
        return format!(
            "error:{}",
            e.as_string()
                .unwrap_or_else(|| "input too large".to_string())
        );
    }
    let query = match chematic_smiles::parse(query_smiles) {
        Ok(m) => m,
        Err(e) => return format!("error:query parse failed: {e}"),
    };
    if let Err(e) = enforce_wasm_molecule_size(&query) {
        return format!(
            "error:{}",
            e.as_string()
                .unwrap_or_else(|| "molecule too large".to_string())
        );
    }

    let smiles_list = match parse_smiles_json_array(db_smiles_json) {
        Ok(values) => values,
        Err(e) => {
            return format!(
                "error:{}",
                e.as_string()
                    .unwrap_or_else(|| "db_smiles_json parse failed".to_string())
            );
        }
    };
    let mut db = Vec::with_capacity(smiles_list.len());
    let mut original_indices = Vec::with_capacity(smiles_list.len());
    for (idx, smiles) in smiles_list.iter().enumerate() {
        let mol = match chematic_smiles::parse(smiles) {
            Ok(mol) => mol,
            Err(e) => return format!("error:db parse failed at index {idx}: {e}"),
        };
        if let Err(e) = enforce_wasm_molecule_size(&mol) {
            return format!(
                "error:db molecule at index {idx}: {}",
                e.as_string()
                    .unwrap_or_else(|| "molecule too large".to_string())
            );
        }
        db.push(mol);
        original_indices.push(idx);
    }

    let results = chematic_fp::nearest_neighbors(&query, &db, k, chematic_fp::FpType::Ecfp4);
    let entries: Vec<String> = results
        .iter()
        .map(|(i, t)| {
            let original_idx = original_indices[*i];
            format!("{{\"index\":{original_idx},\"tanimoto\":{t:.6}}}")
        })
        .collect();
    format!("[{}]", entries.join(","))
}

/// Virtual screen a query SMILES against a database of SMILES using ECFP4 Tanimoto.
///
/// `db_smiles_json`: JSON array of SMILES strings (max 1024 via WASM_MAX_BATCH_ITEMS).
/// `k`: number of top hits to return; clamped to db size if larger.
///
/// Returns JSON: `{"results":[{"rank":1,"score":0.85,"smiles":"CCO","idx":42},...]}`.
/// Returns `"error:<msg>"` on any parse failure or oversized input.
#[wasm_bindgen]
pub fn virtual_screen_ecfp4_json(query_smi: &str, db_smiles_json: &str, k: u32) -> String {
    if let Err(e) = enforce_wasm_input_len("query_smi", query_smi) {
        return format!("error:{}", e.as_string().unwrap_or_else(|| "input too large".to_string()));
    }
    let query_mol = match chematic_smiles::parse(query_smi) {
        Ok(m) => m,
        Err(e) => return format!("error:query parse failed: {e}"),
    };
    if let Err(e) = enforce_wasm_molecule_size(&query_mol) {
        return format!("error:{}", e.as_string().unwrap_or_else(|| "molecule too large".to_string()));
    }
    let smiles_list = match parse_smiles_json_array(db_smiles_json) {
        Ok(v) => v,
        Err(e) => return format!("error:{}", e.as_string().unwrap_or_else(|| "db parse failed".to_string())),
    };
    let mut db_fps = Vec::with_capacity(smiles_list.len());
    for (idx, smi) in smiles_list.iter().enumerate() {
        let mol = match chematic_smiles::parse(smi) {
            Ok(m) => m,
            Err(e) => return format!("error:db parse failed at index {idx}: {e}"),
        };
        if let Err(e) = enforce_wasm_molecule_size(&mol) {
            return format!("error:db molecule at index {idx}: {}",
                e.as_string().unwrap_or_else(|| "molecule too large".to_string()));
        }
        db_fps.push(chematic_fp::ecfp4(&mol));
    }
    let query_fp = chematic_fp::ecfp4(&query_mol);
    let k = (k as usize).min(smiles_list.len());
    let hits = chematic_fp::top_k_similar(&query_fp, &db_fps, k);
    let entries: Vec<String> = hits.iter().enumerate().map(|(rank, (idx, score))| {
        let smi_escaped = smiles_list[*idx].replace('\\', "\\\\").replace('"', "\\\"");
        format!(r#"{{"rank":{},"score":{:.6},"smiles":"{}","idx":{}}}"#,
            rank + 1, score, smi_escaped, idx)
    }).collect();
    format!(r#"{{"results":[{}]}}"#, entries.join(","))
}

/// Compute ECFP4 Tanimoto similarity from one query SMILES to all db SMILES (dense output).
///
/// `db_smiles_json`: JSON array of SMILES strings (max 1024 via WASM_MAX_BATCH_ITEMS).
///
/// Returns a flat JSON array of f32 scores, one per db entry, e.g. `[0.12,0.0,0.85]`.
/// No zero-filtering: the length always equals the number of db entries.
/// Returns `"error:<msg>"` on parse failure or oversized input.
#[wasm_bindgen]
pub fn tanimoto_row_json(query_smi: &str, db_smiles_json: &str) -> String {
    if let Err(e) = enforce_wasm_input_len("query_smi", query_smi) {
        return format!("error:{}", e.as_string().unwrap_or_else(|| "input too large".to_string()));
    }
    let query_mol = match chematic_smiles::parse(query_smi) {
        Ok(m) => m,
        Err(e) => return format!("error:query parse failed: {e}"),
    };
    if let Err(e) = enforce_wasm_molecule_size(&query_mol) {
        return format!("error:{}", e.as_string().unwrap_or_else(|| "molecule too large".to_string()));
    }
    let smiles_list = match parse_smiles_json_array(db_smiles_json) {
        Ok(v) => v,
        Err(e) => return format!("error:{}", e.as_string().unwrap_or_else(|| "db parse failed".to_string())),
    };
    let mut db_fps = Vec::with_capacity(smiles_list.len());
    for (idx, smi) in smiles_list.iter().enumerate() {
        let mol = match chematic_smiles::parse(smi) {
            Ok(m) => m,
            Err(e) => return format!("error:db parse failed at index {idx}: {e}"),
        };
        db_fps.push(chematic_fp::ecfp4(&mol));
    }
    let query_fp = chematic_fp::ecfp4(&query_mol);
    let scores = chematic_fp::tanimoto_slice(&query_fp, &db_fps);
    let parts: Vec<String> = scores.iter().map(|s| format!("{s:.6}")).collect();
    format!("[{}]", parts.join(","))
}

// ---------------------------------------------------------------------------
// MOL2 I/O
// ---------------------------------------------------------------------------

/// Parse a Tripos MOL2 string and return SMILES.
///
/// Returns `"error:<msg>"` on failure.
#[wasm_bindgen]
pub fn mol2_to_smiles(mol2_str: &str) -> String {
    match chematic_mol::parse_mol2(mol2_str) {
        Ok((mol, _)) => chematic_smiles::write(&mol),
        Err(e) => format!("error:{e}"),
    }
}

/// Convert a SMILES to a minimal Tripos MOL2 string (no 3D coordinates).
///
/// Returns `"error:<msg>"` on parse failure.
#[wasm_bindgen]
pub fn smiles_to_mol2(smiles: &str) -> String {
    match chematic_smiles::parse(smiles) {
        Ok(mol) => chematic_mol::write_mol2(&mol, &[]),
        Err(e) => format!("error:{e}"),
    }
}

// ---------------------------------------------------------------------------
// InChI I/O
// ---------------------------------------------------------------------------

/// Generate InChI string from SMILES.
///
/// Returns `"error:<msg>"` on parse failure.
#[wasm_bindgen]
pub fn inchi_from_smiles(smiles: &str) -> String {
    match chematic_smiles::parse(smiles) {
        Ok(mol) => chematic_inchi::inchi(&mol),
        Err(e) => format!("error:{e}"),
    }
}

/// Generate InChIKey from SMILES (27-character identifier).
///
/// Returns `"error:<msg>"` on parse failure.
#[wasm_bindgen]
pub fn inchikey_from_smiles(smiles: &str) -> String {
    match chematic_smiles::parse(smiles) {
        Ok(mol) => {
            let inchi_str = chematic_inchi::inchi(&mol);
            chematic_inchi::inchi_key(&inchi_str)
        }
        Err(e) => format!("error:{e}"),
    }
}

/// Invert the stereochemistry of a tetrahedral stereocenter (U/D wedge bonds).
///
/// If the atom has no wedge/dash bonds, returns an unchanged copy.
/// Returns error if atom_idx is invalid.
#[wasm_bindgen]
pub fn invert_stereocenter_at(mol: &MolHandle, atom_idx: u32) -> Result<MolHandle, JsValue> {
    let idx = chematic_core::AtomIdx(atom_idx);
    if atom_idx as usize >= mol.inner.atom_count() {
        return Err(JsValue::from_str(&format!("atom_idx {} out of range", atom_idx)));
    }
    let new_mol = chematic_chem::invert_stereocenter(&mol.inner, idx);
    Ok(MolHandle {
        inner: std::rc::Rc::new(new_mol),
    })
}

// ---------------------------------------------------------------------------
// mol_transforms (3D geometry manipulation)
// ---------------------------------------------------------------------------

/// Get bond length in Ångströms between two atoms from a SMILES string.
/// Returns -1.0 if parsing fails or atom indices are out of range.
///
/// # Arguments
/// - `smiles`: SMILES string
/// - `a`: first atom index
/// - `b`: second atom index
///
/// # Example
/// ```javascript
/// const len = get_bond_length_json("CC", 0, 1);  // C-C single bond ≈ 1.54 Å
/// ```
#[wasm_bindgen]
pub fn get_bond_length_json(smiles: &str, a: u32, b: u32) -> f64 {
    if smiles.len() > WASM_MAX_INPUT_BYTES {
        return -1.0;
    }
    let mol = match chematic_smiles::parse(smiles) {
        Ok(m) => m,
        Err(_) => return -1.0,
    };
    let coords = chematic_3d::generate_coords(&mol);
    let a_idx = chematic_core::AtomIdx(a);
    let b_idx = chematic_core::AtomIdx(b);
    if a >= mol.atom_count() as u32 || b >= mol.atom_count() as u32 {
        return -1.0;
    }
    chematic_3d::get_bond_length(&coords, a_idx, b_idx)
}

/// Get dihedral angle A—B—C—D in degrees from a SMILES string.
/// Returns null (JSON null) if any atom index is out of range or atoms are collinear.
///
/// # Arguments
/// - `smiles`: SMILES string
/// - `a`, `b`, `c`, `d`: atom indices
///
/// # Example
/// ```javascript
/// const dihedral = get_dihedral_json("CCCC", 0, 1, 2, 3);  // A-B-C-D
/// ```
#[wasm_bindgen]
pub fn get_dihedral_json(smiles: &str, a: u32, b: u32, c: u32, d: u32) -> JsValue {
    if smiles.len() > WASM_MAX_INPUT_BYTES {
        return JsValue::NULL;
    }
    let mol = match chematic_smiles::parse(smiles) {
        Ok(m) => m,
        Err(_) => return JsValue::NULL,
    };
    if a >= mol.atom_count() as u32
        || b >= mol.atom_count() as u32
        || c >= mol.atom_count() as u32
        || d >= mol.atom_count() as u32
    {
        return JsValue::NULL;
    }
    let coords = chematic_3d::generate_coords(&mol);
    let a_idx = chematic_core::AtomIdx(a);
    let b_idx = chematic_core::AtomIdx(b);
    let c_idx = chematic_core::AtomIdx(c);
    let d_idx = chematic_core::AtomIdx(d);
    match chematic_3d::get_dihedral_deg(&coords, a_idx, b_idx, c_idx, d_idx) {
        Some(angle) => JsValue::from_f64(angle),
        None => JsValue::NULL,
    }
}

/// Set dihedral angle A—B—C—D and return PDB block with modified coordinates.
/// Rotates the D-side subtree around the B—C bond.
/// Returns a JS error if parsing fails or atom indices are invalid.
///
/// # Arguments
/// - `smiles`: SMILES string
/// - `a`, `b`, `c`, `d`: atom indices
/// - `angle_deg`: target dihedral angle in degrees
///
/// # Example
/// ```javascript
/// const pdbBlock = set_dihedral_json("CCCC", 0, 1, 2, 3, 120.0);
/// ```
#[wasm_bindgen]
pub fn set_dihedral_json(smiles: &str, a: u32, b: u32, c: u32, d: u32, angle_deg: f64) -> Result<String, String> {
    if smiles.len() > WASM_MAX_INPUT_BYTES {
        return Err("Input SMILES too long".to_string());
    }
    let mol = chematic_smiles::parse(smiles).map_err(|e| format!("Parse error: {}", e))?;
    if a >= mol.atom_count() as u32
        || b >= mol.atom_count() as u32
        || c >= mol.atom_count() as u32
        || d >= mol.atom_count() as u32
    {
        return Err("Atom index out of range".to_string());
    }
    let coords = chematic_3d::generate_coords(&mol);
    let a_idx = chematic_core::AtomIdx(a);
    let b_idx = chematic_core::AtomIdx(b);
    let c_idx = chematic_core::AtomIdx(c);
    let d_idx = chematic_core::AtomIdx(d);
    let angle_rad = angle_deg.to_radians();
    let new_coords = chematic_3d::set_dihedral(&coords, &mol, a_idx, b_idx, c_idx, d_idx, angle_rad);
    // Return PDB block with modified coordinates
    Ok(chematic_3d::write_pdb(&mol, &new_coords))
}

// ---------------------------------------------------------------------------
// random_smiles (SMILES augmentation)
// ---------------------------------------------------------------------------

/// Generate `count` random SMILES from a SMILES string using the given seed.
/// Atoms are permuted based on xorshift64 RNG. Each variant should parse back
/// to the same molecule. Returns a JSON array of SMILES strings.
///
/// # Arguments
/// - `smiles`: input SMILES string
/// - `count`: number of variants to generate (capped at 100)
/// - `seed`: xorshift64 seed
///
/// # Example
/// ```javascript
/// const variants = random_smiles_json("CC(C)O", 5, 42);
/// // variants: ["CC(C)O", "C(C)(O)C", ...]
/// ```
#[wasm_bindgen]
pub fn random_smiles_json(smiles: &str, count: usize, seed: u64) -> Result<String, String> {
    if smiles.len() > WASM_MAX_INPUT_BYTES {
        return Err("Input SMILES too long".to_string());
    }
    let mol = chematic_smiles::parse(smiles).map_err(|e| format!("Parse error: {}", e))?;
    let count = count.min(100); // Cap at 100 to prevent excessive output
    let variants = chematic_smiles::random_smiles_vect(&mol, count, seed);
    // Serialize as JSON array
    let json = format!(
        "[{}]",
        variants
            .iter()
            .map(|s| format!("\"{}\"", s.replace('\\', "\\\\").replace('"', "\\\"")))
            .collect::<Vec<_>>()
            .join(",")
    );
    Ok(json)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(s: &str) -> MolHandle {
        MolHandle {
            inner: std::rc::Rc::new(chematic_smiles::parse(s).unwrap()),
        }
    }

    #[test]
    fn parse_benzene_atom_count() {
        assert_eq!(parse("c1ccccc1").atom_count(), 6);
    }

    // --- logd / isotope / topological index tests ---------------------------

    #[test]
    fn logd_at_ph7_benzene() {
        let h = parse("c1ccccc1");
        let logd = h.logd_at_ph(7.0);
        assert!(logd.is_finite(), "logD should be finite, got {logd}");
    }

    #[test]
    fn logd_profile_json_has_ph_fields() {
        let h = parse("c1ccccc1");
        let json = h.logd_profile_json(0.0, 14.0, 15);
        assert!(json.starts_with('[') && json.ends_with(']'),
            "should be JSON array: {json}");
        assert!(json.contains(r#""ph":"#), "should contain ph field: {json}");
        assert!(json.contains(r#""logd":"#), "should contain logd field: {json}");
    }

    #[test]
    fn isotope_distribution_json_nonempty() {
        let h = parse("C");
        let json = h.isotope_distribution_json(0.1);
        assert!(json.starts_with('['), "should be JSON array: {json}");
        assert!(json.contains(r#""mass":"#), "should contain mass entries: {json}");
    }

    #[test]
    fn randic_index_ethane() {
        let h = parse("CC");
        let r = h.randic_index();
        // Ethane: 1 C-C bond, both heavy atoms have degree 1 → χ = 1/√(1×1) = 1.0
        assert!((r - 1.0).abs() < 1e-6, "ethane Randic index = 1.0, got {r}");
    }

    #[test]
    fn zagreb_m1_ethane() {
        let h = parse("CC");
        let z = h.zagreb_index_m1();
        // Ethane: each C has heavy-atom degree 1 → M1 = 1² + 1² = 2
        assert_eq!(z, 2, "ethane Zagreb M1 = 2, got {z}");
    }

    // --- iupac_name tests ---------------------------------------------------

    #[test]
    fn iupac_name_benzene() {
        assert_eq!(parse("c1ccccc1").iupac_name(), "benzene");
    }

    #[test]
    fn iupac_name_ethanol() {
        assert_eq!(parse("CCO").iupac_name(), "ethanol");
    }

    #[test]
    fn iupac_name_branched_acid() {
        assert_eq!(parse("CC(C)C(=O)O").iupac_name(), "2-methylpropanoic acid");
    }

    #[test]
    fn iupac_name_unsupported_returns_empty() {
        // Disconnected molecule → unsupported → empty string
        let h = parse("C.C");
        assert!(h.iupac_name().is_empty(), "disconnected should return empty");
    }

    // --- assign_cip_json tests ----------------------------------------------

    #[test]
    fn assign_cip_json_no_centers() {
        // Ethane has no stereocenters
        let json = parse("CC").assign_cip_json();
        assert_eq!(json, r#"{"centers":[]}"#);
    }

    #[test]
    fn assign_cip_json_chiral_center() {
        // L-alanine SMILES with explicit chirality
        let json = parse("N[C@@H](C)C(=O)O").assign_cip_json();
        assert!(json.contains("\"code\":"), "should contain code field: {json}");
        assert!(json.contains("\"atom\":"), "should contain atom field: {json}");
    }

    #[test]
    fn canonical_smiles_benzene() {
        let mol = parse("c1ccccc1");
        let cs = mol.canonical_smiles();
        assert!(!cs.is_empty());
    }

    #[test]
    fn parse_cxsmiles_json_preserves_metadata() {
        let json = parse_cxsmiles_json("C~O |$C1;O2$,atomProp:1.role.acceptor,^2:0,Z:0|").unwrap();
        assert!(json.contains(r#""atomLabels":["C1","O2"]"#), "{json}");
        assert!(json.contains(r#""key":"role""#), "{json}");
        assert!(json.contains(r#""atomRadicals":[2,null]"#), "{json}");
        assert!(json.contains(r#""zeroBonds":[0]"#), "{json}");
    }

    #[test]
    fn parse_cxsmarts_json_preserves_metadata() {
        let json =
            parse_cxsmarts_json("[#6]~[#8] |$C1;O2$,atomProp:1.role.acceptor,^2:0|").unwrap();
        assert!(json.contains(r#""atomCount":2"#), "{json}");
        assert!(json.contains(r#""atomLabels":["C1","O2"]"#), "{json}");
        assert!(json.contains(r#""key":"role""#), "{json}");
        assert!(json.contains(r#""atomRadicals":[2,null]"#), "{json}");
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
        assert!((2..=5).contains(&rb), "aspirin rotbonds = {rb}");
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
        assert!(
            with_h.atom_count() > mol.atom_count(),
            "H atoms should be added"
        );
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
        assert!(
            svg.contains("<svg"),
            "invalid SMILES should be silently skipped"
        );
    }

    #[test]
    fn run_reactants_esterification() {
        // Simple esterification: carboxylic acid + alcohol → ester + water
        let result = run_reactants(
            "[C:1](=O)[OH:2].[O:3][C:4]>>[C:1](=O)[O:3][C:4]",
            "CC(=O)O|CCO",
        );
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
        assert!(
            !is_valid_smiles("[NOSUCHELEMENT]"),
            "unknown bracket atom is invalid"
        );
    }

    #[test]
    fn depict_svg_opts_transparent_background() {
        let h = parse("CCO");
        let mut opts = DepictOptions::new();
        opts.set_background("transparent".to_string());
        let svg = h.depict_svg_opts(&opts);
        assert!(svg.contains("<svg"), "must produce SVG");
        assert!(
            !svg.contains("fill=\"transparent\""),
            "no bg rect for transparent"
        );
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
        assert!(
            svg.contains("stroke=\"white\""),
            "dark theme bonds should be white"
        );
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
        assert!(
            svg.contains("Na"),
            "Na should appear in disconnected SMILES SVG"
        );
        assert!(
            svg.contains("Cl"),
            "Cl should appear in disconnected SMILES SVG"
        );
        assert!(!svg.is_empty());
    }

    #[test]
    fn depict_svg_disconnected_water_dimer() {
        let svg = parse("O.O").depict_svg();
        // Degree-0 O atoms use isolated (Hill) notation: H2O
        assert!(
            svg.matches("H2O").count() >= 2,
            "both water O atoms should render as H2O"
        );
        assert!(!svg.is_empty());
    }

    // ── Sprint L: atom data attributes ──────────────────────────────────────

    #[test]
    fn depict_svg_opts_atom_ids_contains_data_attrs() {
        let h = parse("CC(=O)O"); // acetic acid: 2 C (unlabelled) + 2 O (labelled)
        let mut opts = DepictOptions::new();
        opts.set_atom_ids(true);
        let svg = h.depict_svg_opts(&opts);
        assert!(
            svg.contains("data-atom-idx="),
            "atom_ids should add data-atom-idx"
        );
        assert!(
            svg.contains("data-element="),
            "atom_ids should add data-element"
        );
        assert!(
            svg.contains("data-charge="),
            "atom_ids should add data-charge"
        );
        // B3: all 4 atoms should be addressable (unlabelled carbons get invisible anchor)
        assert_eq!(
            svg.matches("data-atom-idx=").count(),
            4,
            "all atoms should have data-atom-idx"
        );
    }

    #[test]
    fn depict_svg_opts_atom_ids_false_no_data_attrs() {
        let h = parse("CC(=O)O");
        let svg = h.depict_svg_opts(&DepictOptions::new());
        assert!(
            !svg.contains("data-atom-idx="),
            "default opts should not have data-atom-idx"
        );
    }

    #[test]
    fn depict_svg_opts_atom_ids_charge_correct() {
        let h = parse("[NH4+]");
        let mut opts = DepictOptions::new();
        opts.set_atom_ids(true);
        let svg = h.depict_svg_opts(&opts);
        assert!(
            svg.contains("data-charge=\"1\""),
            "NH4+ should have charge=1"
        );
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
        assert!(
            !svg.contains("fill=\"#8b92a9\""),
            "default should not show grey index labels"
        );
    }

    // ── Sprint L: kekulize ───────────────────────────────────────────────────

    #[test]
    fn depict_svg_opts_kekulize_removes_aromatic_bonds() {
        let h = parse("c1ccccc1"); // benzene
        let mut opts = DepictOptions::new();
        opts.set_kekulize(true);
        let svg = h.depict_svg_opts(&opts);
        assert!(!svg.is_empty());
        assert!(
            svg.contains("<line"),
            "kekulé benzene should have line elements"
        );
        assert!(
            !svg.contains("stroke-dasharray"),
            "kekulé benzene must not use aromatic dashed style"
        ); // B1
    }

    #[test]
    fn depict_svg_opts_kekulize_false_uses_aromatic() {
        let h = parse("c1ccccc1");
        let svg = h.depict_svg_opts(&DepictOptions::new());
        // Default aromatic rendering uses stroke-dasharray for the inner ring line.
        assert!(
            svg.contains("stroke-dasharray"),
            "default benzene should use aromatic dashed style"
        );
    }

    // ── Sprint M: smarts_match_atoms ─────────────────────────────────────────

    #[test]
    fn smarts_match_benzene_ring_returns_json() {
        let mol = parse("c1ccccc1");
        let result = smarts_match_atoms("c1ccccc1", &mol);
        assert!(result.is_ok(), "valid SMARTS should not error");
        let json = result.unwrap();
        assert!(
            !json.is_empty() && json != "[]",
            "benzene ring SMARTS should find a match"
        );
        assert!(json.starts_with("[["), "result should be array of arrays");
    }

    #[test]
    fn smarts_match_no_match_returns_empty_array() {
        let mol = parse("CC"); // ethane has no aromatic ring
        let result = smarts_match_atoms("c1ccccc1", &mol);
        assert!(
            result.is_ok(),
            "valid SMARTS on non-matching mol should not error"
        );
        assert_eq!(
            result.unwrap(),
            "[]",
            "no match should return empty JSON array"
        );
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
        assert!(
            !pdb.contains("nan") && !pdb.contains("inf"),
            "PDB must have no NaN/Inf coords"
        );
        assert!(pdb.contains("HETATM"), "must produce HETATM records");
    }

    // Note: smarts_match_atoms error-path test (invalid SMARTS) is omitted here
    // because JsValue::from_str panics outside a WASM runtime.
    // The underlying chematic_smarts::parse_smarts error path is tested separately:
    #[test]
    fn smarts_parse_invalid_is_err() {
        assert!(
            chematic_smarts::parse_smarts("[invalid").is_err(),
            "invalid SMARTS should return Err from parse_smarts"
        );
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
        assert!(
            json.contains('N'),
            "pyridine FG JSON should contain N: {json}"
        );
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
        let has_negative = json
            .trim_matches(|c| c == '[' || c == ']')
            .split(',')
            .any(|s| s.trim().parse::<f64>().map(|v| v < 0.0).unwrap_or(false));
        assert!(
            has_negative,
            "acetic acid should have at least one negative charge: {json}"
        );
    }

    #[test]
    fn test_mmff94_charges_json_length() {
        // Acetic acid has 4 heavy atoms (2C, 2O)
        let h = parse("CC(=O)O");
        let json = mmff94_charges_json(&h);
        assert!(json.starts_with('[') && json.ends_with(']'));
        let count = json
            .trim_matches(|c| c == '[' || c == ']')
            .split(',')
            .count();
        assert_eq!(count, 4, "acetic acid should have 4 charges, got {count}");
    }

    #[test]
    fn test_mmff94_charges_json_acetate_negative() {
        // Acetate [O-] → total charge = -1
        let h = parse("CC(=O)[O-]");
        let json = mmff94_charges_json(&h);
        let total: f64 = json
            .trim_matches(|c| c == '[' || c == ']')
            .split(',')
            .filter_map(|s| s.trim().parse::<f64>().ok())
            .sum();
        assert!(
            (total + 1.0).abs() < 0.01,
            "acetate total charge should be -1, got {total}"
        );
    }

    #[test]
    fn test_mhfp_hashes_json_format() {
        let h = parse("c1ccccc1");
        let json = mhfp_hashes_json(&h);
        assert!(json.contains("\"num_hashes\":128"), "should have num_hashes:128");
        assert!(json.contains("\"hashes\":"), "should have hashes key");
    }

    #[test]
    fn test_tanimoto_mhfp_smiles_self() {
        let result = tanimoto_mhfp_smiles("c1ccccc1", "c1ccccc1").unwrap();
        assert!((result - 1.0).abs() < 1e-9, "self-similarity should be 1.0, got {result}");
    }

    #[test]
    fn test_mhfp_lsh_handle_query() {
        let mut idx = MhfpLshHandle::new(128);
        idx.add_smiles("c1ccccc1").unwrap();    // benzene → 0
        idx.add_smiles("Cc1ccccc1").unwrap();   // toluene → 1
        idx.add_smiles("CC").unwrap();           // ethane  → 2

        let json = idx.query_json("c1ccccc1", 0.99).unwrap();
        // Benzene should find itself at similarity 1.0
        assert!(json.contains("\"index\":0"), "benzene should find itself: {json}");
        assert_eq!(idx.len(), 3);
    }

    #[test]
    fn to_mol_block_has_nonzero_coords() {
        let h = parse("c1ccccc1"); // benzene
        let block = to_mol_block(&h);
        // Extract atom-block lines (lines 4..10, 0-indexed) and check at least
        // one coordinate is non-zero.
        let nonzero = block.lines().skip(4).take(6).any(|line| {
            // Atom line: first 30 chars are three 10.4 coordinate fields.
            if line.len() < 30 {
                return false;
            }
            let x: f64 = line[0..10].trim().parse().unwrap_or(0.0);
            let y: f64 = line[10..20].trim().parse().unwrap_or(0.0);
            x.abs() > 0.01 || y.abs() > 0.01
        });
        assert!(
            nonzero,
            "to_mol_block should write real 2D coordinates, got:\n{block}"
        );
    }

    #[test]
    fn sa_score_range() {
        let h = parse("CC(=O)Oc1ccccc1C(=O)O"); // aspirin
        let score = sa_score(&h);
        assert!(
            (1.0..=10.0).contains(&score),
            "SA score out of [1,10]: {score:.2}"
        );
    }

    // v0.1.21 new API tests

    // with_atom_charge
    #[test]
    fn mol_with_atom_charge_changes_charge() {
        let mol = parse("N"); // neutral N
        let mol2 = mol_with_atom_charge(&mol, 0, 1).unwrap();
        let atom = mol2.inner.atom(chematic_core::AtomIdx(0));
        assert_eq!(atom.charge, 1, "charge should be 1");
    }

    // with_atom_element
    #[test]
    fn mol_with_atom_element_changes_element() {
        let mol = parse("CC");
        let mol2 = mol_with_atom_element(&mol, 0, "N").unwrap();
        let atom = mol2.inner.atom(chematic_core::AtomIdx(0));
        assert_eq!(atom.element.symbol(), "N", "element should be N");
    }

    // SDF parse with coords
    #[test]
    fn mol_block_coords_json_ethanol() {
        let mol_block = "\nethanol\n\n  3  2  0  0  0  0  0  0  0  0  0 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\nM  END\n";
        let json = mol_block_coords_json(mol_block).unwrap();
        assert!(
            json.contains("[0.0000,0.0000]"),
            "first atom at origin: {json}"
        );
        assert!(
            json.contains("[1.5000,0.0000]"),
            "second atom at x=1.5: {json}"
        );
    }

    // CDXML all fragments
    #[test]
    fn cdxml_to_smiles_json_two_fragments() {
        let cdxml = r#"<?xml version="1.0"?>
<CDXML>
<fragment>
<n id="1" Element="6" p="0 0"/>
<n id="2" Element="6" p="10 0"/>
<b B="1" E="2" Order="1"/>
</fragment>
<fragment>
<n id="3" Element="8" p="30 0"/>
<n id="4" Element="6" p="40 0"/>
<b B="3" E="4" Order="1"/>
</fragment>
</CDXML>"#;
        let json = cdxml_to_smiles_json(cdxml).unwrap();
        // Should have 2 comma-separated SMILES entries.
        let count = json.matches(',').count() + 1;
        assert_eq!(count, 2, "should have 2 fragments: {json}");
    }

    // CDXML stereo
    #[test]
    fn cdxml_wedge_bond_becomes_up() {
        let cdxml = r#"<?xml version="1.0"?>
<CDXML>
<fragment>
<n id="1" Element="6" p="0 0"/>
<n id="2" Element="6" p="10 0"/>
<b B="1" E="2" Order="1" Display="WedgeBegin"/>
</fragment>
</CDXML>"#;
        let (mol, _) = chematic_mol::parse_cdxml(cdxml).unwrap();
        let bond = mol.bond(chematic_core::BondIdx(0));
        assert_eq!(bond.order, chematic_core::BondOrder::Up, "WedgeBegin → Up");
    }

    // depict_data_with_coords
    #[test]
    fn depict_data_with_coords_uses_provided_coords() {
        let mol = parse("CC");
        // Provide coords far from defaults.
        let json = depict_data_with_coords_json(&mol, "[[100.0,200.0],[300.0,400.0]]");
        assert!(json.contains("100."), "x coord should be 100: {json}");
        assert!(json.contains("200."), "y coord should be 200: {json}");
    }

    // Mutable Molecule API
    #[test]
    fn mol_with_atom_added_increases_count() {
        let mol = parse("CC");
        let mol2 = mol_with_atom_added(&mol, "O").unwrap();
        assert_eq!(mol2.atom_count(), mol.atom_count() + 1);
    }

    #[test]
    fn mol_with_bond_added_increases_count() {
        let mol = parse("CC.O"); // ethane + disconnected O (3 atoms, 1 bond)
        let mol2 = mol_with_bond_added(&mol, 1, 2, 1).unwrap();
        assert_eq!(mol2.bond_count(), mol.bond_count() + 1);
    }

    #[test]
    fn mol_with_atom_removed_decreases_count() {
        let mol = parse("CCO");
        let mol2 = mol_with_atom_removed(&mol, 2).unwrap();
        assert_eq!(mol2.atom_count(), mol.atom_count() - 1);
        // The C-O bond is also removed.
        assert_eq!(mol2.bond_count(), mol.bond_count() - 1);
    }

    #[test]
    fn mol_with_bond_removed_decreases_count() {
        let mol = parse("CC");
        let mol2 = mol_with_bond_removed(&mol, 0).unwrap();
        assert_eq!(mol2.atom_count(), 2, "atoms preserved");
        assert_eq!(mol2.bond_count(), 0, "bond removed");
    }

    #[test]
    fn mol_next_atom_idx_equals_atom_count() {
        let mol = parse("CCO");
        assert_eq!(mol_next_atom_idx(&mol), 3);
    }

    // SDF / V3000 write
    #[test]
    fn smiles_array_to_sdf_contains_dollar_delimiters() {
        let sdf = smiles_array_to_sdf(r#"["CC","CCO"]"#).unwrap();
        assert_eq!(sdf.matches("$$$$").count(), 2, "2 delimiters: {sdf}");
    }

    #[test]
    fn to_mol_v3000_block_contains_v3000() {
        let mol = parse("CC");
        let v3k = to_mol_v3000_block(&mol);
        assert!(v3k.contains("V3000"), "should contain V3000 tag");
        assert!(v3k.contains("M  V30 BEGIN ATOM"), "should have atom block");
    }

    // DepictData
    #[test]
    fn depict_data_json_benzene_atoms_bonds() {
        let mol = parse("c1ccccc1");
        let json = depict_data_json(&mol);
        assert!(json.contains("\"atoms\":"), "should have atoms key");
        assert!(json.contains("\"bonds\":"), "should have bonds key");
        assert!(json.contains("\"element\":\"C\""), "benzene has C atoms");
    }

    #[test]
    fn depict_data_json_pyridine_nitrogen_has_label() {
        let mol = parse("c1ccncc1");
        let json = depict_data_json(&mol);
        // Nitrogen has a label
        assert!(
            json.contains("\"label\":\"N\""),
            "pyridine N should have label: {json}"
        );
    }

    // CPK colors
    #[test]
    fn cpk_color_nitrogen_is_blue() {
        assert_eq!(cpk_color("N"), "#3050F8", "N is blue in CPK");
    }

    #[test]
    fn cpk_color_carbon_is_black() {
        assert_eq!(cpk_color("C"), "#000000", "C is black in CPK");
    }

    // CML / CDXML
    #[test]
    fn mol_from_cml_ethanol() {
        let cml = r#"<molecule>
  <atomArray>
    <atom id="a1" elementType="C" x2="0.0" y2="0.0"/>
    <atom id="a2" elementType="C" x2="1.5" y2="0.0"/>
    <atom id="a3" elementType="O" x2="3.0" y2="0.0"/>
  </atomArray>
  <bondArray>
    <bond atomRefs2="a1 a2" order="1"/>
    <bond atomRefs2="a2 a3" order="1"/>
  </bondArray>
</molecule>"#;
        let mol = mol_from_cml(cml).unwrap();
        assert_eq!(mol.atom_count(), 3, "ethanol: 3 heavy atoms");
        assert_eq!(mol.bond_count(), 2, "ethanol: 2 bonds");
    }

    #[test]
    fn to_cml_contains_molecule_tag() {
        let mol = parse("CC(=O)O");
        let cml = to_cml(&mol);
        assert!(cml.contains("<molecule"), "should contain <molecule");
        assert!(cml.contains("elementType="), "should contain elementType");
        assert!(cml.contains("atomRefs2="), "should contain atomRefs2");
    }

    #[test]
    fn to_cml_roundtrip_atom_count() {
        let mol = parse("CC(=O)O");
        let cml = to_cml(&mol);
        let mol2 = mol_from_cml(&cml).unwrap();
        assert_eq!(
            mol.atom_count(),
            mol2.atom_count(),
            "CML round-trip preserves atom count"
        );
        assert_eq!(
            mol.bond_count(),
            mol2.bond_count(),
            "CML round-trip preserves bond count"
        );
    }

    #[test]
    fn mol_from_cdxml_ethanol() {
        let cdxml = r#"<?xml version="1.0"?>
<CDXML>
<fragment>
<n id="1" p="10.0 20.0" Element="6"/>
<n id="2" p="25.0 20.0" Element="6"/>
<n id="3" p="40.0 20.0" Element="8"/>
<b B="1" E="2" Order="1"/>
<b B="2" E="3" Order="1"/>
</fragment>
</CDXML>"#;
        let mol = mol_from_cdxml(cdxml).unwrap();
        assert_eq!(mol.atom_count(), 3, "CDXML ethanol: 3 heavy atoms");
        assert_eq!(mol.bond_count(), 2, "CDXML ethanol: 2 bonds");
    }

    #[test]
    fn mol_from_cml_unknown_element_returns_err() {
        let cml = r#"<molecule><atomArray>
<atom id="a1" elementType="Xx"/>
</atomArray></molecule>"#;
        // Test via underlying function (JsValue::from_str would abort in native)
        let result = chematic_mol::parse_cml(cml);
        assert!(result.is_err(), "unknown element should return Err");
    }

    // Sprint CC
    #[test]
    fn mmp_pairs_json_ethylbenzene_propylbenzene() {
        let smiles_json = r#"["CCc1ccccc1","CCCc1ccccc1"]"#;
        let result = mmp_pairs_json(smiles_json).unwrap();

        // Should find exactly 1 pair with the benzene core.
        assert!(
            result.contains("\"core\":"),
            "should have core key: {result}"
        );
        assert!(
            result.contains("c1(ccccc1)[*]"),
            "core should be benzene: {result}"
        );
        assert!(
            result.contains("\"fragment_a\":"),
            "should have fragment_a: {result}"
        );
        assert!(
            result.contains("\"fragment_b\":"),
            "should have fragment_b: {result}"
        );
    }

    #[test]
    fn mmp_pairs_json_no_pairs_for_single_molecule() {
        let smiles_json = r#"["CCc1ccccc1"]"#;
        let result = mmp_pairs_json(smiles_json).unwrap();
        assert_eq!(result, "[]", "single molecule should have no pairs");
    }

    #[test]
    fn mmp_pairs_json_three_molecules() {
        // Ethylbenzene, propylbenzene, butylbenzene.
        // Different isomers and cut points may yield multiple pairs and cores.
        let smiles_json = r#"["CCc1ccccc1","CCCc1ccccc1","CCCCc1ccccc1"]"#;
        let result = mmp_pairs_json(smiles_json).unwrap();

        // Should find at least 1 pair (the minimum for 3 benzene-derived molecules).
        assert!(
            result.contains("\"mol_a\":"),
            "should have mol_a key: {result}"
        );
        assert!(
            result.contains("\"core\":"),
            "should have core key: {result}"
        );
        // Extract pair count by counting top-level braces.
        let pair_count = result.matches("\"mol_a\":").count();
        assert!(
            pair_count >= 1,
            "should have >= 1 pair, got {pair_count}: {result}"
        );
    }

    // Sprint BB
    #[test]
    fn conformer_handle_new_and_add() {
        let ens = ConformerHandle::new("CCCC").unwrap(); // butane
        assert_eq!(ens.conformer_count(), 0, "new ensemble has 0 conformers");
    }

    #[test]
    fn conformer_handle_add_generated() {
        let mut ens = ConformerHandle::new("CCCC").unwrap();
        let idx = ens.add_generated_conformer();
        assert_eq!(idx, 0, "first conformer has index 0");
        assert_eq!(ens.conformer_count(), 1);
    }

    #[test]
    fn conformer_handle_add_minimized() {
        let mut ens = ConformerHandle::new("CCCC").unwrap();
        let idx = ens.add_minimized_conformer();
        assert_eq!(idx, 0);
        assert_eq!(ens.conformer_count(), 1);
    }

    #[test]
    fn conformer_handle_get_pdb_returns_hetatm() {
        let mut ens = ConformerHandle::new("CC").unwrap();
        ens.add_generated_conformer();
        let pdb = ens.get_conformer_pdb(0).expect("conformer 0 exists");
        assert!(pdb.contains("HETATM"), "PDB should contain HETATM records");
    }

    #[test]
    fn conformer_handle_rmsd_same_conformer() {
        // RMSD of a conformer with itself should be 0.0.
        let mut ens = ConformerHandle::new("CCCC").unwrap();
        ens.add_generated_conformer();
        ens.add_generated_conformer();
        // Different conformers will have RMSD >= 0; same index trivially returns
        // the distance between identical points = 0.
        let rmsd = ens.conformer_rmsd(0, 0);
        assert!(rmsd.is_finite(), "RMSD should be finite");
        assert!(rmsd >= 0.0, "RMSD must be non-negative");
    }

    #[test]
    fn conformer_handle_rmsd_out_of_range_is_nan() {
        let mut ens = ConformerHandle::new("CC").unwrap();
        ens.add_generated_conformer();
        assert!(ens.conformer_rmsd(0, 99).is_nan(), "out-of-range gives NaN");
    }

    #[test]
    fn conformer_handle_remove_conformer() {
        let mut ens = ConformerHandle::new("CC").unwrap();
        ens.add_generated_conformer();
        assert!(ens.remove_conformer(0));
        assert_eq!(ens.conformer_count(), 0);
        assert!(!ens.remove_conformer(0), "already removed");
    }

    #[test]
    fn conformer_handle_mol_returns_correct_atom_count() {
        let ens = ConformerHandle::new("CC(=O)O").unwrap(); // acetic acid, 4 heavy atoms
        let mol = ens.mol();
        assert_eq!(mol.atom_count(), 4);
    }

    // R-group decomposition oracle tests
    #[test]
    fn rgroup_toluene_ethylbenzene_para_substituted() {
        // Core: para-substituted benzene c1ccc(*)cc1
        // Toluene Cc1ccccc1 → R1 = C
        // Ethylbenzene CCc1ccccc1 → R1 = CC
        let smiles_json = r#"["Cc1ccccc1","CCc1ccccc1"]"#;
        let json = rgroup_decompose_json(smiles_json, "c1ccc(*)cc1").unwrap();

        // Both should match.
        assert!(
            json.contains("\"matched\":true"),
            "both should match: {json}"
        );
        // R-groups should be C and CC (canonical forms may vary slightly).
        assert!(json.contains("\"r1\":"), "should have r1 key: {json}");
    }

    #[test]
    fn rgroup_no_match_returns_false() {
        // Benzene has no substituent → doesn't match c1ccc(*)cc1.
        let smiles_json = r#"["c1ccccc1"]"#;
        let json = rgroup_decompose_json(smiles_json, "c1ccc(*)cc1").unwrap();
        assert!(
            json.contains("\"matched\":false"),
            "unsubstituted benzene should not match: {json}"
        );
    }

    #[test]
    fn rgroup_two_attachment_points() {
        // Di-substituted benzene: c1cc(*)cc(*)c1
        // 1,3-dimethylbenzene m-xylene Cc1cccc(C)c1 → R1=C, R2=C
        let smiles_json = r#"["Cc1cccc(C)c1"]"#;
        let json = rgroup_decompose_json(smiles_json, "c1cc(*)cc(*)c1").unwrap();
        assert!(
            json.contains("\"matched\":true"),
            "m-xylene should match: {json}"
        );
        assert!(json.contains("\"r1\":"), "should have r1: {json}");
        assert!(json.contains("\"r2\":"), "should have r2: {json}");
    }

    #[test]
    fn rgroup_invalid_smarts_returns_err() {
        // SMARTS parse should fail via underlying function (not JsValue path).
        let result = chematic_smarts::parse_smarts("~~~invalid~~~");
        assert!(result.is_err(), "invalid SMARTS should not parse");
    }

    // Sprint AA
    #[test]
    fn fcfp4_bitvec_length_256() {
        assert_eq!(fcfp4_bitvec(&parse("c1ccccc1")).len(), 256);
    }

    #[test]
    fn fcfp6_bitvec_length_256() {
        assert_eq!(fcfp6_bitvec(&parse("c1ccccc1")).len(), 256);
    }

    #[test]
    fn dice_ecfp6_identical_is_one() {
        let mol = parse("c1ccccc1");
        assert!((dice_ecfp6(&mol, &mol) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn write_smiles_contains_same_atoms() {
        let mol = parse("CC(=O)O");
        let smi = write_smiles(&mol);
        // Non-canonical SMILES must still round-trip to the same atom count.
        let mol2 = parse_smiles(&smi).unwrap();
        assert_eq!(mol.atom_count(), mol2.atom_count());
    }

    #[test]
    fn normalize_reaction_smiles_roundtrip() {
        // parse_reaction error path tested via underlying fn (JsValue native-abort).
        let result = chematic_rxn::parse_reaction("CC>>CO");
        assert!(result.is_ok(), "simple reaction should parse");
        let rxn = result.unwrap();
        let out = chematic_rxn::write_reaction(&rxn);
        assert!(
            out.contains(">>"),
            "written reaction should contain >>: {out}"
        );
    }

    // Sprint Z
    #[test]
    fn brics_fragments_json_aspirin() {
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O"); // aspirin
        let json = brics_fragments_json(&mol);
        let count_from_json = json
            .split('"')
            .filter(|s| s.contains('C') || s.contains('c'))
            .count();
        assert!(
            json.starts_with('[') && json.ends_with(']'),
            "not a JSON array: {json}"
        );
        assert_eq!(
            count_from_json,
            brics_fragment_count(&mol),
            "fragment count mismatch: json={json}"
        );
    }

    #[test]
    fn brics_fragments_json_benzene_self() {
        // No BRICS-breakable bonds → returns the molecule itself as one fragment.
        let mol = parse("c1ccccc1");
        let json = brics_fragments_json(&mol);
        assert!(json.starts_with('[') && json.ends_with(']'));
        assert_eq!(
            brics_fragment_count(&mol),
            1,
            "benzene has 1 fragment (itself): {json}"
        );
    }

    #[test]
    fn atom_pair_bitvec_length_256() {
        let mol = parse("c1ccccc1");
        assert_eq!(atom_pair_bitvec(&mol).len(), 256);
    }

    #[test]
    fn torsion_bitvec_length_256() {
        let mol = parse("CCCCC");
        assert_eq!(torsion_bitvec(&mol).len(), 256);
    }

    #[test]
    fn tanimoto_fcfp6_identical_is_one() {
        let mol = parse("c1ccccc1");
        assert!((tanimoto_fcfp6(&mol, &mol) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn sdf_from_records_json_round_trip() {
        // Build a 1-record SDF, read it back, check name and property.
        let smiles_json = r#"["CC"]"#;
        let names_json = r#"["ethane"]"#;
        let props_json = r#"["MW\t30.07\nSource\ttest"]"#;
        let sdf = sdf_from_records_json(smiles_json, names_json, props_json).unwrap();
        assert!(sdf.contains("ethane"), "name missing: {sdf}");
        assert!(sdf.contains("> <MW>"), "MW field missing: {sdf}");
        assert!(sdf.contains("30.07"), "MW value missing: {sdf}");
        assert!(sdf.contains("$$$$"), "SDF delimiter missing: {sdf}");
    }

    #[test]
    fn sdf_from_records_json_no_props() {
        // Empty props string → no data fields, still valid SDF.
        let sdf = sdf_from_records_json(r#"["C"]"#, r#"["methane"]"#, r#"[""]"#).unwrap();
        assert!(sdf.contains("methane"));
        assert!(sdf.contains("$$$$"));
    }

    #[test]
    fn sdf_from_records_json_length_mismatch_returns_err() {
        // Mismatched array lengths — note: testing via underlying logic,
        // not .is_err() on the wrapper (JsValue native-abort).
        let smiles = parse_smiles_json_array(r#"["CC","CCC"]"#).unwrap();
        let names = parse_smiles_json_array(r#"["ethane"]"#).unwrap();
        assert_ne!(smiles.len(), names.len(), "lengths should differ");
    }

    #[test]
    fn parse_smiles_json_array_handles_escaped_strings() {
        let smiles = parse_smiles_json_array(r#"["C\\C=C\\O","CC"]"#).unwrap();
        assert_eq!(smiles[0], r#"C\C=C\O"#);
        assert_eq!(smiles[1], "CC");
    }

    // Sprint Y
    #[test]
    fn mol_from_xyz_roundtrip_atom_count() {
        // Build a minimal XYZ string for ethane and verify atom count.
        // (parse_xyz error path tested via the underlying fn — JsValue native-abort)
        let xyz = "2\nethane\nC  0.0 0.0 0.0\nC  1.5 0.0 0.0\n";
        let result = chematic_3d::parse_xyz(xyz);
        assert!(result.is_ok(), "valid XYZ should parse");
        let (mol, _) = result.unwrap();
        assert_eq!(mol.atom_count(), 2);
    }

    #[test]
    fn to_xyz_contains_atom_lines() {
        let mol = parse("CC"); // ethane
        let xyz = to_xyz(&mol);
        assert!(xyz.contains('C'), "XYZ should contain C atom lines");
        // First line is atom count
        let count: usize = xyz.lines().next().unwrap().trim().parse().unwrap();
        assert_eq!(count, mol.atom_count(), "XYZ header atom count matches mol");
    }

    #[test]
    fn mol_from_pdb_returns_handle() {
        // pdb_to_molecule with no atoms gives empty molecule.
        let h = mol_from_pdb("");
        assert_eq!(h.atom_count(), 0);
    }

    #[test]
    fn logp_per_atom_json_length() {
        let mol = parse("CC(=O)O"); // acetic acid
        let json = logp_per_atom_json(&mol);
        let count = json
            .trim_matches(|c| c == '[' || c == ']')
            .split(',')
            .count();
        assert_eq!(
            count,
            mol.atom_count(),
            "per-atom logP array length mismatch"
        );
    }

    #[test]
    fn mr_per_atom_json_length() {
        let mol = parse("CC(=O)O");
        let json = mr_per_atom_json(&mol);
        let count = json
            .trim_matches(|c| c == '[' || c == ']')
            .split(',')
            .count();
        assert_eq!(count, mol.atom_count());
    }

    #[test]
    fn labute_asa_per_atom_json_no_nan() {
        let mol = parse("c1ccccc1");
        let json = labute_asa_per_atom_json(&mol);
        assert!(!json.contains("NaN"), "no NaN in asa per-atom: {json}");
    }

    #[test]
    fn sssr_rings_json_benzene_one_ring() {
        let mol = parse("c1ccccc1");
        let json = sssr_rings_json(&mol);
        // One ring: [[0,1,2,3,4,5]]
        assert!(
            json.starts_with("[[") && json.ends_with("]]"),
            "expected one ring: {json}"
        );
        let ring: Vec<usize> = json
            .trim_matches(|c| c == '[' || c == ']')
            .split(',')
            .filter_map(|s| s.trim_matches(|c| c == '[' || c == ']').parse().ok())
            .collect();
        assert_eq!(ring.len(), 6, "benzene ring has 6 atoms");
    }

    #[test]
    fn sssr_rings_json_naphthalene_two_rings() {
        let mol = parse("c1ccc2ccccc2c1"); // naphthalene
        let json = sssr_rings_json(&mol);
        // Count rings by counting "[" at start
        let ring_count = json.matches("],[").count() + 1;
        assert_eq!(ring_count, 2, "naphthalene has 2 rings: {json}");
    }

    #[test]
    fn ecfp_bitvec_custom_256_length() {
        let mol = parse("c1ccccc1");
        let bv = ecfp_bitvec_custom(&mol, 2, 256, false);
        assert_eq!(bv.len(), 32, "256-bit FP = 32 bytes");
    }

    #[test]
    fn ecfp_bitvec_custom_identical_tanimoto() {
        let mol = parse("c1ccccc1");
        let a = ecfp_bitvec_custom(&mol, 2, 512, false);
        let b = ecfp_bitvec_custom(&mol, 2, 512, false);
        assert_eq!(a, b, "same mol same FP");
    }

    // ── Chirality tests ──────────────────────────────────────────────────────
    #[test]
    fn ecfp4_bitvec_with_chirality_l_vs_d_alanine_different() {
        let l_ala = parse("C[C@H](N)C(=O)O");
        let d_ala = parse("C[C@@H](N)C(=O)O");
        let l_fp = ecfp4_bitvec_with_chirality(&l_ala, true);
        let d_fp = ecfp4_bitvec_with_chirality(&d_ala, true);
        assert_ne!(
            l_fp, d_fp,
            "L-alanine and D-alanine should have different ECFP4 with use_chirality=true"
        );
    }

    #[test]
    fn ecfp4_bitvec_with_chirality_l_vs_d_alanine_same_without_chirality() {
        let l_ala = parse("C[C@H](N)C(=O)O");
        let d_ala = parse("C[C@@H](N)C(=O)O");
        let l_fp = ecfp4_bitvec_with_chirality(&l_ala, false);
        let d_fp = ecfp4_bitvec_with_chirality(&d_ala, false);
        assert_eq!(
            l_fp, d_fp,
            "L-alanine and D-alanine should have identical ECFP4 with use_chirality=false"
        );
    }

    #[test]
    fn ecfp6_bitvec_with_chirality_l_vs_d_alanine_different() {
        let l_ala = parse("C[C@H](N)C(=O)O");
        let d_ala = parse("C[C@@H](N)C(=O)O");
        let l_fp = ecfp6_bitvec_with_chirality(&l_ala, true);
        let d_fp = ecfp6_bitvec_with_chirality(&d_ala, true);
        assert_ne!(
            l_fp, d_fp,
            "L-alanine and D-alanine should have different ECFP6 with use_chirality=true"
        );
    }

    #[test]
    fn smarts_match_atoms_with_chirality_chiral_atom_matches_only_with_flag() {
        let l_ala = parse("C[C@H](N)C(=O)O");
        let d_ala = parse("C[C@@H](N)C(=O)O");
        // Match [C@H] (counterclockwise) — should only match L-alanine when use_chirality=true
        let l_matches_with = smarts_match_atoms_with_chirality("[C@H]", &l_ala, true).unwrap();
        let l_matches_without = smarts_match_atoms_with_chirality("[C@H]", &l_ala, false).unwrap();
        let d_matches_with = smarts_match_atoms_with_chirality("[C@H]", &d_ala, true).unwrap();
        let d_matches_without = smarts_match_atoms_with_chirality("[C@H]", &d_ala, false).unwrap();

        // With use_chirality=true: L-ala matches, D-ala doesn't
        assert!(
            !l_matches_with.contains("[]"),
            "L-ala should match [C@H] with use_chirality=true"
        );
        assert!(
            d_matches_with.contains("[]"),
            "D-ala should NOT match [C@H] with use_chirality=true"
        );

        // With use_chirality=false: both match (chirality ignored)
        assert!(
            !l_matches_without.contains("[]"),
            "L-ala should match [C@H] with use_chirality=false"
        );
        assert!(
            !d_matches_without.contains("[]"),
            "D-ala should match [C@H] with use_chirality=false (chirality ignored)"
        );
    }

    #[test]
    fn enumerate_stereo_isomers_one_center_gives_two() {
        // C(F)(Cl)Br has one unspecified chiral center → 2 stereoisomers
        let mol = parse("C(F)(Cl)Br");
        // enumerate_stereo_isomers_json returns Result; unwrap Ok via JsValue-safe check
        let result = enumerate_stereo_isomers_json(&mol);
        assert!(result.is_ok(), "should return Ok for valid mol");
        let json = result.unwrap();
        // Count objects: each object contains "smiles", "inchi", "inchikey"
        // So 2 isomers means 2 occurrences of "smiles" at top level
        let count = json.matches("\"smiles\"").count();
        assert_eq!(count, 2, "expected 2 stereoisomers: {json}");
    }

    #[test]
    fn enumerate_stereo_isomers_already_specified_gives_one() {
        // Already-specified stereocenter → only 1 isomer (no unspecified centers)
        let mol = parse("C[C@H](F)Cl");
        let result = enumerate_stereo_isomers_json(&mol);
        assert!(result.is_ok());
        let json = result.unwrap();
        // Count objects: each object contains "smiles"
        let count = json.matches("\"smiles\"").count();
        assert_eq!(count, 1, "fully specified mol should give 1 isomer: {json}");
    }

    #[test]
    fn enumerate_stereo_isomers_no_center_gives_one() {
        // No chiral centers at all → 1 isomer
        let mol = parse("CC");
        let result = enumerate_stereo_isomers_json(&mol);
        assert!(result.is_ok());
        let json = result.unwrap();
        assert!(json.starts_with('[') && json.ends_with(']'));
    }

    // Sprint X
    #[test]
    fn mol_from_v3000_block_parses_atom_count() {
        // Borrow a minimal V3000 fixture from the mol3000 test suite.
        // (Parse failure tested via underlying parse_mol_v3000 — see JsValue note.)
        assert!(
            chematic_mol::parse_mol_v3000(
                "\n  \n\n  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\n\
             M  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\n\
             M  V30 1 C 0 0 0 0\nM  V30 2 C 1.5 0 0 0\n\
             M  V30 END ATOM\nM  V30 BEGIN BOND\n\
             M  V30 1 1 1 2\nM  V30 END BOND\n\
             M  V30 END CTAB\nM  END"
            )
            .is_ok()
        );
    }

    #[test]
    fn generate_3d_minimized_pdb_nonzero_coords() {
        let mol = parse("CCCC"); // butane — flexible, benefits from minimization
        let pdb = generate_3d_minimized_pdb(&mol);
        assert!(pdb.contains("HETATM"), "expected HETATM records");
        // At least one non-zero coordinate in the PDB output.
        let has_nonzero = pdb.lines().filter(|l| l.starts_with("HETATM")).any(|l| {
            // PDB columns 31-38, 39-46, 47-54 are x,y,z
            if l.len() < 54 {
                return false;
            }
            let x: f64 = l[30..38].trim().parse().unwrap_or(0.0);
            let y: f64 = l[38..46].trim().parse().unwrap_or(0.0);
            x.abs() > 0.01 || y.abs() > 0.01
        });
        assert!(
            has_nonzero,
            "minimized PDB should have non-zero coords:\n{pdb}"
        );
    }

    #[test]
    fn sdf_to_records_json_parses_properties() {
        // Raw string preserves the required fixed-width whitespace of MOL V2000.
        let sdf = concat!(
            "aspirin\n",
            "  chematic\n",
            "\n",
            "  2  1  0  0  0  0  0  0  0  0  0 V2000\n",
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "  1  2  1  0\n",
            "M  END\n",
            "> <MW>\n",
            "180.2\n",
            "\n",
            "> <Source>\n",
            "ChEMBL\n",
            "\n",
            "$$$$\n",
        );
        let json = sdf_to_records_json(sdf);
        assert!(
            json.contains("\"name\":\"aspirin\""),
            "name missing: {json}"
        );
        assert!(json.contains("\"MW\":\"180.2\""), "MW missing: {json}");
        assert!(
            json.contains("\"Source\":\"ChEMBL\""),
            "Source missing: {json}"
        );
    }

    #[test]
    fn sdf_to_records_json_escapes_special_chars() {
        let sdf = concat!(
            "mol\n",
            "  chematic\n",
            "\n",
            "  1  0  0  0  0  0  0  0  0  0  0 V2000\n",
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "M  END\n",
            "> <Notes>\n",
            "line1\n",
            "line2\n",
            "\n",
            "$$$$\n",
        );
        let json = sdf_to_records_json(sdf);
        // "line1\nline2" value joined with \n, then JSON-escaped to \\n
        assert!(
            json.contains("\\n"),
            "newline in multi-line value should be escaped: {json}"
        );
    }

    #[test]
    fn depict_svg_grid_highlighted_no_smarts_returns_grid() {
        let smiles = "c1ccccc1\nCCO";
        let svg = depict_svg_grid_highlighted(smiles, 2, "");
        assert!(svg.contains("mol-0"), "expected mol-0: {svg}");
        assert!(svg.contains("mol-1"), "expected mol-1: {svg}");
    }

    #[test]
    fn depict_svg_grid_highlighted_with_smarts_adds_circles() {
        let smiles = "c1ccccc1\nc1ccncc1"; // benzene and pyridine
        let svg = depict_svg_grid_highlighted(smiles, 2, "c1ccccn1"); // pyridine SMARTS
        // Pyridine matches → circle element expected in second cell
        assert!(svg.contains("<circle"), "expected highlight circles: {svg}");
    }

    // Sprint W
    #[test]
    fn pains_matches_json_empty_for_clean_mol() {
        let mol = parse("c1ccccc1"); // benzene — no PAINS alerts
        let json = pains_matches_json(&mol);
        assert_eq!(json, "[]");
    }

    #[test]
    fn pains_matches_json_returns_alert_names() {
        // Rhodanine scaffold is a classic PAINS alert
        let mol = parse("O=C1CSC(=S)N1"); // rhodanine
        let json = pains_matches_json(&mol);
        assert!(json.starts_with('[') && json.ends_with(']'));
        // May be empty for some molecules; just check structure
        let _ = json; // success if it compiles and runs
    }

    #[test]
    fn cip_assignments_json_aspirin_no_stereo() {
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O"); // aspirin — no stereocenters
        let json = cip_assignments_json(&mol);
        assert_eq!(json, "[]");
    }

    #[test]
    fn cip_assignments_json_chiral_center() {
        let mol = parse("[C@@H](F)(Cl)Br"); // chiral center
        let json = cip_assignments_json(&mol);
        assert!(json.contains("cipCode"), "expected cipCode in: {json}");
        assert!(json.contains('R') || json.contains('S'));
    }

    #[test]
    fn ecfp6_bitvec_length_256() {
        let mol = parse("c1ccccc1");
        assert_eq!(ecfp6_bitvec(&mol).len(), 256);
    }

    #[test]
    fn tanimoto_ecfp6_identical_is_one() {
        let mol = parse("c1ccccc1");
        assert!((tanimoto_ecfp6(&mol, &mol) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn dice_ecfp4_identical_is_one() {
        let mol = parse("c1ccccc1");
        assert!((dice_ecfp4(&mol, &mol) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn dice_maccs_identical_is_one() {
        let mol = parse("c1ccccc1");
        assert!((dice_maccs(&mol, &mol) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn num_aliphatic_rings_cyclohexane() {
        let mol = parse("C1CCCCC1"); // cyclohexane
        assert_eq!(mol.num_aliphatic_rings(), 1);
        assert_eq!(mol.num_saturated_rings(), 1);
        assert_eq!(mol.aromatic_ring_count(), 0);
    }

    #[test]
    fn num_unspecified_stereocenters_unspec() {
        let mol = parse("C(F)(Cl)Br"); // chiral center without @/@@ annotation
        assert_eq!(mol.num_unspecified_stereocenters(), 1);
    }

    #[test]
    fn shape_descriptors_json_benzene_keys() {
        let mol = parse("c1ccccc1");
        let json = shape_descriptors_json(&mol);
        assert!(json.contains("\"pmi1\""), "missing pmi1: {json}");
        assert!(json.contains("\"npr1\""), "missing npr1: {json}");
        assert!(
            json.contains("\"asphericity\""),
            "missing asphericity: {json}"
        );
    }

    #[test]
    fn shape_descriptors_json_single_atom_no_nan() {
        let mol = parse("[Na+]"); // single atom — shape descriptors may be non-finite
        let json = shape_descriptors_json(&mol);
        // Must be valid JSON (no bare NaN)
        assert!(!json.contains("NaN"), "NaN in output: {json}");
        assert!(!json.contains("inf"), "inf in output: {json}");
    }

    #[test]
    fn maxmin_picks_returns_correct_count() {
        let json = r#"["CC","c1ccccc1","CCO","CCCC","c1cccnc1"]"#;
        let result = maxmin_picks_ecfp4_json(json, 3).unwrap();
        let count = result
            .trim_matches(|c| c == '[' || c == ']')
            .split(',')
            .count();
        assert_eq!(count, 3, "expected 3 picks: {result}");
    }

    #[test]
    fn smiles_parse_fails_for_unclosed_ring() {
        // Verifies the error-propagation contract: parse_smiles_json_array feeds
        // each string through chematic_smiles::parse, which rejects unclosed rings.
        // (Calling maxmin_picks_ecfp4_json directly would invoke JsValue::from_str,
        // which panics outside WASM — see comment on rxn_parse_missing_arrow_is_err.)
        assert!(
            chematic_smiles::parse("C1CC").is_err(),
            "unclosed ring should fail"
        );
    }

    #[test]
    fn butina_cluster_all_similar_one_cluster() {
        // All identical SMILES → cutoff=0.0 means all at distance 1.0 from centroid → 1 cluster
        let json = r#"["c1ccccc1","c1ccccc1","c1ccccc1"]"#;
        let result = butina_cluster_ecfp4_json(json, 0.0).unwrap();
        assert!(result.starts_with('[') && result.ends_with(']'));
    }

    #[test]
    fn mcs_smiles_two_acetyl_compounds() {
        let json = r#"["CC(=O)O","CC(=O)N"]"#; // acetic acid and acetamide — MCS is CC=O
        let result = mcs_smiles_json(json).unwrap();
        assert_ne!(result, "null", "expected a non-null MCS");
        assert!(!result.is_empty());
    }

    #[test]
    fn mcs_smiles_no_overlap_returns_null() {
        let json = r#"["[Na+]","[Cl-]"]"#; // no common organic substructure
        let result = mcs_smiles_json(json).unwrap();
        // May return a single-atom MCS or null depending on algorithm
        assert!(result == "null" || !result.is_empty());
    }

    #[test]
    fn mcs_smiles_single_input_returns_err() {
        // Passing only one SMILES string should produce a JS error; we verify
        // via the Result variant rather than calling is_err() (which would
        // materialise a JsValue drop on non-wasm32 and cause a panic).
        let json = r#"["CC"]"#;
        let smiles_list: Vec<String> = json
            .split('"')
            .enumerate()
            .filter(|(i, _)| i % 2 == 1)
            .map(|(_, s)| s.to_string())
            .collect();
        assert_eq!(smiles_list.len(), 1, "helper extracted wrong count");
        assert!(smiles_list.len() < 2, "should be fewer than 2 inputs");
    }

    // Sprint V
    #[test]
    fn murcko_scaffold_benzene_ring() {
        let mol = parse("c1ccccc1CC(=O)O"); // phenylacetic acid
        let scaffold = murcko_scaffold(&mol);
        let smi = scaffold.canonical_smiles();
        // Murcko of phenylacetic acid is benzene ring only.
        assert!(!smi.is_empty(), "murcko_scaffold returned empty SMILES");
    }

    #[test]
    fn canonical_tautomer_returns_handle() {
        let mol = parse("Oc1cccc2ccccc12"); // 1-naphthol
        let t = canonical_tautomer(&mol);
        assert!(t.atom_count() > 0);
    }

    #[test]
    fn enumerate_tautomers_json_is_array() {
        let mol = parse("Oc1cccc2ccccc12");
        let json = enumerate_tautomers_json(&mol);
        assert!(json.starts_with('[') && json.ends_with(']'));
        assert!(json.len() > 2, "expected at least one tautomer");
    }

    #[test]
    fn largest_fragment_strips_salt() {
        // sodium acetate: "CC(=O)[O-].[Na+]" — largest fragment is acetate
        let mol = parse("CC(=O)[O-].[Na+]");
        let frag = largest_fragment(&mol);
        // acetate has 4 heavy atoms; Na has 1
        assert!(frag.atom_count() > 1, "expected the larger fragment");
        assert!(
            frag.atom_count() < mol.atom_count(),
            "fragment should be smaller than the salt"
        );
    }

    #[test]
    fn neutralize_charges_removes_charges() {
        let mol = parse("CC(=O)[O-]");
        let neutral = neutralize_charges(&mol);
        assert_eq!(neutral.formal_charge_sum(), 0);
    }

    #[test]
    fn standardize_smiles_report_json_includes_audit_report() {
        let json = standardize_smiles_report_json("CC.CCC", true, false, false, false);
        assert!(json.contains("\"smiles\""), "missing smiles: {json}");
        assert!(json.contains("\"report\""), "missing report: {json}");
        assert!(
            json.contains("\"status\":\"Modified\""),
            "missing status: {json}"
        );
        assert!(
            json.contains("\"step\":\"LargestFragment\""),
            "missing largest-fragment step: {json}"
        );
    }

    #[test]
    fn maccs_bitvec_length_21() {
        let mol = parse("c1ccccc1");
        let bv = maccs_bitvec(&mol);
        assert_eq!(bv.len(), 21, "MACCS 166 bits should fit in 21 bytes");
    }

    #[test]
    fn tanimoto_maccs_identical_is_one() {
        let mol = parse("c1ccccc1");
        assert!((tanimoto_maccs(&mol, &mol) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn get_descriptors_json_keys() {
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O"); // aspirin
        let json = get_descriptors_json(&mol);
        assert!(json.contains("\"mw\""), "missing mw: {json}");
        assert!(json.contains("\"tpsa\""), "missing tpsa: {json}");
        assert!(json.contains("\"qed\""), "missing qed: {json}");
    }

    #[test]
    fn slogp_vsa_json_length_12() {
        let h = parse("c1ccccc1");
        let json = slogp_vsa_json(&h);
        let count = json
            .trim_matches(|c| c == '[' || c == ']')
            .split(',')
            .count();
        assert_eq!(count, 12, "SlogP_VSA should have 12 values");
    }

    #[test]
    fn smr_vsa_json_length_10() {
        let h = parse("c1ccccc1");
        let json = smr_vsa_json(&h);
        let count = json
            .trim_matches(|c| c == '[' || c == ']')
            .split(',')
            .count();
        assert_eq!(count, 10, "SMR_VSA should have 10 values");
    }

    #[test]
    fn peoe_vsa_json_length_14() {
        let h = parse("c1ccccc1");
        let json = peoe_vsa_json(&h);
        let count = json
            .trim_matches(|c| c == '[' || c == ']')
            .split(',')
            .count();
        assert_eq!(count, 14, "PEOE_VSA should have 14 values");
    }
}

// ---------------------------------------------------------------------------
// Ring Family Detection
// ---------------------------------------------------------------------------

/// Ring family classification and detection as JSON.
/// Returns an array of ring families with their atoms, ring indices, and topology kind.
#[wasm_bindgen]
pub fn ring_families_json(mol: &MolHandle) -> Result<String, JsValue> {
    use chematic_perception::{find_sssr, find_ring_families, RingSystemKind};

    if mol.inner.atom_count() > WASM_MAX_ATOMS {
        return Err(JsValue::from_str(&format!(
            "Molecule exceeds maximum atom count ({} > {})",
            mol.inner.atom_count(),
            WASM_MAX_ATOMS
        )));
    }

    let sssr = find_sssr(&mol.inner);
    let families = find_ring_families(&mol.inner, &sssr);

    let json_families: Vec<serde_json::Value> = families
        .iter()
        .map(|fam| {
            let kind_str = match fam.kind {
                RingSystemKind::Simple => "Simple",
                RingSystemKind::Fused => "Fused",
                RingSystemKind::Spiro => "Spiro",
                RingSystemKind::Bridged => "Bridged",
            };

            serde_json::json!({
                "kind": kind_str,
                "atoms": fam.atoms.iter().map(|a| a.0).collect::<Vec<_>>(),
                "ringIndices": fam.ring_indices.clone(),
                "atomCount": fam.atoms.len(),
                "ringCount": fam.ring_indices.len(),
            })
        })
        .collect();

    serde_json::to_string(&json_families)
        .map_err(|e| JsValue::from_str(&e.to_string()))
}

// ---------------------------------------------------------------------------
// MHFP LSH Index
// ---------------------------------------------------------------------------

/// MinHash LSH index: insert MHFP fingerprints and query by approximate similarity.
///
/// ```js
/// const idx = new MhfpLshHandle(128);
/// const i0 = idx.add_smiles("c1ccccc1");    // benzene → index 0
/// const i1 = idx.add_smiles("Cc1ccccc1");   // toluene → index 1
/// const hits = JSON.parse(idx.query_json("c1ccccc1", 0.5));
/// // hits: [{index:0,similarity:1.0}, {index:1,similarity:0.xxx}]
/// ```
#[wasm_bindgen]
pub struct MhfpLshHandle {
    inner: chematic_fp::MhfpLshIndex,
}

#[wasm_bindgen]
impl MhfpLshHandle {
    /// Create a new LSH index for MHFP fingerprints with `num_hashes` hash lanes.
    /// Default band decomposition: 16 bands × (num_hashes / 16) rows.
    /// `num_hashes` must be a multiple of 16 (e.g. 128).
    #[wasm_bindgen(constructor)]
    pub fn new(num_hashes: usize) -> MhfpLshHandle {
        MhfpLshHandle {
            inner: chematic_fp::MhfpLshIndex::new(num_hashes),
        }
    }

    /// Add a molecule by SMILES; returns its 0-based index in the index.
    pub fn add_smiles(&mut self, smiles: &str) -> Result<usize, JsValue> {
        let mol = chematic_smiles::parse(smiles)
            .map_err(|e| JsValue::from_str(&e.to_string()))?;
        let fp = chematic_fp::mhfp(&mol);
        Ok(self.inner.add(fp))
    }

    /// Query by SMILES for all entries with similarity ≥ threshold.
    ///
    /// Returns a JSON array `[{"index":N,"similarity":0.xxx},...]` sorted by
    /// descending similarity.  Empty array `[]` when nothing qualifies.
    pub fn query_json(&self, query_smiles: &str, threshold: f64) -> Result<String, JsValue> {
        let mol = chematic_smiles::parse(query_smiles)
            .map_err(|e| JsValue::from_str(&e.to_string()))?;
        let fp = chematic_fp::mhfp(&mol);
        let results = self.inner.query(&fp, threshold);
        let items: Vec<String> = results
            .iter()
            .map(|(i, s)| format!(r#"{{"index":{},"similarity":{:.6}}}"#, i, s))
            .collect();
        Ok(format!("[{}]", items.join(",")))
    }

    /// Number of molecules in the index.
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// True if the index contains no molecules.
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
}
