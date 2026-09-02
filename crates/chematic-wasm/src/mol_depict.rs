//! 2D depiction (SVG), SMARTS matching, and misc reporting bindings.

use crate::{
    DepictOptions, MolHandle, WASM_MAX_ATOMS, WASM_MAX_BATCH_ITEMS, WASM_MAX_INPUT_BYTES,
    WASM_MAX_SMARTS_MATCHES, bond_in_ring, enforce_wasm_input_len, enforce_wasm_molecule_size,
    escape_json_string,
};
use wasm_bindgen::prelude::*;

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

/// Render a grid SVG from newline-separated SMILES (one per line).
///
/// Lines that fail to parse are silently skipped.
/// `cols` controls the number of columns (each cell is 200×200 px).
/// The input is limited to 1 MiB, 1,024 non-empty records, and 10,000 atoms
/// per successfully parsed molecule. Limit violations return an empty SVG
/// containing an error title.
#[wasm_bindgen]
pub fn depict_svg_grid(smiles_block: &str, cols: usize) -> String {
    if smiles_block.len() > WASM_MAX_INPUT_BYTES {
        return wasm_error_svg("input exceeds maximum size");
    }
    if nonempty_line_count(smiles_block) > WASM_MAX_BATCH_ITEMS {
        return wasm_error_svg("molecule count exceeds maximum (1024)");
    }
    if has_oversized_molecule(smiles_block) {
        return wasm_error_svg("molecule exceeds maximum atom count (10000)");
    }
    let mols: Vec<chematic_core::Molecule> = smiles_block
        .lines()
        .filter(|s| !s.trim().is_empty())
        .filter_map(|s| chematic_smiles::parse(s.trim()).ok())
        .collect();
    let refs: Vec<&chematic_core::Molecule> = mols.iter().collect();
    chematic_depict::depict_svg_grid(&refs, cols)
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
    const WASM_MAX_SMARTS_VISITS: u64 = 1_000_000;
    let config = chematic_smarts::MatchConfig {
        max_matches: Some(WASM_MAX_SMARTS_MATCHES),
        use_chirality: false,
        use_isotopes: false,
        uniquify: true,
        max_visit_budget: Some(WASM_MAX_SMARTS_VISITS),
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
    const WASM_MAX_SMARTS_VISITS: u64 = 1_000_000;
    let config = chematic_smarts::MatchConfig {
        max_matches: Some(WASM_MAX_SMARTS_MATCHES),
        use_chirality,
        use_isotopes: false,
        uniquify: true,
        max_visit_budget: Some(WASM_MAX_SMARTS_VISITS),
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

/// Compute ERG-style 315-element float histogram fingerprint.
/// Returns JSON: {"len":315,"values":[f64,...]} or {"error":"..."}.
/// Format: 21 pharmacophore-feature-pair × 15 distance bins with Gaussian fuzzing.
/// See `chematic_fp::erg_vec` for details.
#[wasm_bindgen]
pub fn erg_vec_json(mol: &MolHandle) -> String {
    if mol.inner.atom_count() > WASM_MAX_ATOMS {
        return format!(
            r#"{{"error":"molecule too large (max {} atoms)"}}"#,
            WASM_MAX_ATOMS
        );
    }
    let v = chematic_fp::erg_vec(&mol.inner);
    let vals: Vec<String> = v.iter().map(|x| format!("{x:.6}")).collect();
    format!(
        r#"{{"len":{},"values":[{}]}}"#,
        chematic_fp::ERG_VEC_LEN,
        vals.join(",")
    )
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

/// Render a molecule grid with SMARTS-based atom highlighting.
///
/// `smiles_block` — newline-separated SMILES strings (same format as `depict_svg_grid`).
/// `cols` — number of grid columns.
/// `match_smarts` — SMARTS pattern; matched atoms in each molecule are highlighted.
///   Pass an empty string `""` to render without any highlighting.
///
/// Invalid SMILES are rendered as empty cells; SMARTS parse failure returns an
/// unhighlighted grid (the SMARTS is silently ignored).
/// The input is limited to 1 MiB, 1,024 non-empty records, and 10,000 atoms
/// per successfully parsed molecule. Limit violations return an empty SVG
/// containing an error title.
#[wasm_bindgen]
pub fn depict_svg_grid_highlighted(smiles_block: &str, cols: usize, match_smarts: &str) -> String {
    if smiles_block.len() > WASM_MAX_INPUT_BYTES || match_smarts.len() > WASM_MAX_INPUT_BYTES {
        return wasm_error_svg("input exceeds maximum size");
    }
    if nonempty_line_count(smiles_block) > WASM_MAX_BATCH_ITEMS {
        return wasm_error_svg("molecule count exceeds maximum (1024)");
    }
    if has_oversized_molecule(smiles_block) {
        return wasm_error_svg("molecule exceeds maximum atom count (10000)");
    }
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

/// Non-canonical SMILES for `mol`.
///
/// Unlike `canonical_smiles`, the output depends on the internal atom ordering
/// and is not normalised.  Useful when round-trip fidelity (preserving atom
/// order) matters more than a canonical form.
#[wasm_bindgen]
pub fn write_smiles(mol: &MolHandle) -> String {
    chematic_smiles::write(&mol.inner)
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
// Ring Family Detection
// ---------------------------------------------------------------------------

/// Ring family classification and detection as JSON.
/// Returns an array of ring families with their atoms, ring indices, and topology kind.
#[wasm_bindgen]
pub fn ring_families_json(mol: &MolHandle) -> Result<String, JsValue> {
    use chematic_perception::{RingSystemKind, find_ring_families, find_sssr};

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

    serde_json::to_string(&json_families).map_err(|e| JsValue::from_str(&e.to_string()))
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
/// Generate a self-contained HTML report for a newline-separated list of SMILES.
///
/// Empty lines and invalid SMILES are silently skipped.
/// Returns the same card-grid HTML as Python's `chematic.report()`.
/// The input is limited to 1 MiB, 1,024 non-empty records, and 10,000 atoms
/// per successfully parsed molecule. Limit violations return an HTML error.
///
/// ```js
/// const html = mod.batch_report_html("CCO\nc1ccccc1\nCC(=O)O");
/// const blob = new Blob([html], {type:'text/html'});
/// const url  = URL.createObjectURL(blob);
/// ```
#[wasm_bindgen]
pub fn batch_report_html(smiles_lines: &str) -> String {
    use chematic_chem::{
        brenk_passes, hba_count_lipinski, hbd_count, logp_and_mr, molecular_weight, pains_passes,
        qed_with_bundle, ring_bundle, tpsa,
    };

    if smiles_lines.len() > WASM_MAX_INPUT_BYTES {
        return wasm_error_html("input exceeds maximum size");
    }
    if nonempty_line_count(smiles_lines) > WASM_MAX_BATCH_ITEMS {
        return wasm_error_html("molecule count exceeds maximum (1024)");
    }
    if has_oversized_molecule(smiles_lines) {
        return wasm_error_html("molecule exceeds maximum atom count (10000)");
    }

    let mut cards: Vec<(f64, String)> = smiles_lines
        .lines()
        .filter_map(|line| {
            let s = line.trim();
            if s.is_empty() { return None; }
            chematic_smiles::parse(s).ok().map(|m| (s.to_string(), m))
        })
        .map(|(label, m)| {
            let mw = molecular_weight(&m);
            let (logp, _) = logp_and_mr(&m);
            let tpsa_val = tpsa(&m);
            let hbd = hbd_count(&m);
            let hba_lip = hba_count_lipinski(&m);
            let rb = ring_bundle(&m);
            let qed = qed_with_bundle(&m, &rb);
            let lip = mw <= 500.0 && hbd <= 5 && hba_lip <= 10 && logp <= 5.0;
            let pains_ok = pains_passes(&m);
            let brenk_ok = brenk_passes(&m);
            let svg = chematic_depict::depict_svg(&m);

            let lip_badge = if lip { r#"<span class="badge pass">Lipinski ✓</span>"# }
                            else   { r#"<span class="badge fail">Lipinski ✗</span>"# };
            let pains_badge = if pains_ok { r#"<span class="badge pass">PAINS ✓</span>"# }
                              else        { r#"<span class="badge fail">PAINS ✗</span>"# };
            let brenk_badge = if brenk_ok { r#"<span class="badge pass">Brenk ✓</span>"# }
                              else        { r#"<span class="badge warn">Brenk ⚠</span>"# };

            let card = format!(
                r#"<div class="card"><div class="svg">{svg}</div><div class="name">{label}</div><div class="desc">MW: {mw:.1} Da &nbsp;|&nbsp; LogP: {logp:.2} &nbsp;|&nbsp; TPSA: {tpsa_val:.1} Å²<br>HBD: {hbd} &nbsp;|&nbsp; HBA: {hba} &nbsp;|&nbsp; QED: {qed:.2}</div><div class="badges">{lip_badge}{pains_badge}{brenk_badge}</div></div>"#,
                hba = rb.hba_count,
            );
            (qed, card)
        })
        .collect();

    cards.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    let cards_html: String = cards.into_iter().map(|(_, c)| c).collect();
    let n = smiles_lines
        .lines()
        .filter(|l| !l.trim().is_empty())
        .count();
    let ver = env!("CARGO_PKG_VERSION");

    format!(
        r#"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>chematic report</title>
<style>
body{{font-family:system-ui,sans-serif;background:#f8f9fa;padding:24px;margin:0}}
h1{{font-size:1.4rem;color:#333;margin-bottom:4px}}
.meta{{font-size:.85rem;color:#666;margin-bottom:20px}}
.grid{{display:flex;flex-wrap:wrap;gap:16px}}
.card{{background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:12px;width:220px;box-shadow:0 1px 3px rgba(0,0,0,.06)}}
.svg{{width:100%;height:160px;overflow:hidden}}
.svg svg{{width:100%;height:100%}}
.name{{font-weight:600;font-size:.9rem;margin:6px 0 4px;color:#333;min-height:1.1em;word-break:break-all}}
.desc{{font-size:.78rem;color:#555;line-height:1.8;margin:4px 0}}
.badges{{display:flex;flex-wrap:wrap;gap:4px;margin-top:6px}}
.badge{{font-size:.7rem;padding:2px 7px;border-radius:10px;font-weight:500}}
.pass{{background:#d1e7dd;color:#0a3622}}
.fail{{background:#f8d7da;color:#58151c}}
.warn{{background:#fff3cd;color:#664d03}}
</style>
</head>
<body>
<h1>chematic report</h1>
<p class="meta">{n} molecule{plural} &middot; generated by chematic v{ver} (WASM)</p>
<div class="grid">{cards_html}</div>
</body>
</html>"#,
        plural = if n == 1 { "" } else { "s" },
    )
}

fn nonempty_line_count(input: &str) -> usize {
    input.lines().filter(|line| !line.trim().is_empty()).count()
}

fn has_oversized_molecule(input: &str) -> bool {
    input
        .lines()
        .filter_map(|line| chematic_smiles::parse(line.trim()).ok())
        .any(|mol| mol.atom_count() > WASM_MAX_ATOMS)
}

fn wasm_error_svg(message: &str) -> String {
    format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="0" height="0"><title>chematic error: {message}</title></svg>"#
    )
}

fn wasm_error_html(message: &str) -> String {
    format!(
        "<!doctype html><meta charset=\"utf-8\"><title>chematic error</title><p>chematic error: {message}</p>"
    )
}

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
        let mol = chematic_smiles::parse(smiles).map_err(|e| JsValue::from_str(&e.to_string()))?;
        let fp = chematic_fp::mhfp(&mol);
        Ok(self.inner.add(fp))
    }

    /// Query by SMILES for all entries with similarity ≥ threshold.
    ///
    /// Returns a JSON array `[{"index":N,"similarity":0.xxx},...]` sorted by
    /// descending similarity.  Empty array `[]` when nothing qualifies.
    pub fn query_json(&self, query_smiles: &str, threshold: f64) -> Result<String, JsValue> {
        let mol =
            chematic_smiles::parse(query_smiles).map_err(|e| JsValue::from_str(&e.to_string()))?;
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
