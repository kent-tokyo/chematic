//! Reaction (SMIRKS/MMP/BRICS/MCS), tautomer, and standardization bindings.

use crate::{
    MolHandle, WASM_MAX_ATOMS, WASM_MAX_BATCH_ITEMS, WASM_MAX_JSON_STRING_BYTES,
    WASM_MAX_SMARTS_MATCHES, enforce_wasm_input_len, enforce_wasm_molecule_size,
    escape_json_string, json_option_string_array, json_option_u8_array, parse_smiles_json_array,
    rgroup_fragment_smiles,
};
use serde::Deserialize;
use wasm_bindgen::prelude::*;

#[derive(Debug, Clone, Copy, Deserialize)]
enum AtomCompareJson {
    #[serde(rename = "elements")]
    Elements,
    #[serde(rename = "any_heavy_atom")]
    AnyHeavyAtom,
    #[serde(rename = "any")]
    Any,
}

impl From<AtomCompareJson> for chematic_smarts::AtomCompare {
    fn from(v: AtomCompareJson) -> Self {
        match v {
            AtomCompareJson::Elements => chematic_smarts::AtomCompare::Elements,
            AtomCompareJson::AnyHeavyAtom => chematic_smarts::AtomCompare::AnyHeavyAtom,
            AtomCompareJson::Any => chematic_smarts::AtomCompare::Any,
        }
    }
}

#[derive(Debug, Clone, Copy, Deserialize)]
enum BondCompareJson {
    #[serde(rename = "order_or_aromatic")]
    OrderOrAromatic,
    #[serde(rename = "any")]
    Any,
}

impl From<BondCompareJson> for chematic_smarts::BondCompare {
    fn from(v: BondCompareJson) -> Self {
        match v {
            BondCompareJson::OrderOrAromatic => chematic_smarts::BondCompare::OrderOrAromatic,
            BondCompareJson::Any => chematic_smarts::BondCompare::Any,
        }
    }
}

fn default_true() -> bool {
    true
}
fn default_min_atoms() -> usize {
    1
}
fn default_atom_compare() -> AtomCompareJson {
    AtomCompareJson::Elements
}
fn default_bond_compare() -> BondCompareJson {
    BondCompareJson::OrderOrAromatic
}

/// JSON config for [`mcs_smiles_json_with_config`] -- every field optional,
/// defaulting to exactly `McsConfig::default()`'s value.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields, rename_all = "camelCase")]
struct McsConfigJson {
    #[serde(default = "default_true")]
    match_bonds: bool,
    #[serde(default = "default_min_atoms")]
    min_atoms: usize,
    #[serde(default)]
    timeout_ms: Option<u64>,
    #[serde(default)]
    ring_matches_ring_only: bool,
    #[serde(default)]
    complete_rings_only: bool,
    #[serde(default = "default_atom_compare")]
    atom_compare: AtomCompareJson,
    #[serde(default = "default_bond_compare")]
    bond_compare: BondCompareJson,
    #[serde(default)]
    match_chiral_tag: bool,
    #[serde(default)]
    match_charge: bool,
    #[serde(default)]
    match_isotope: bool,
    #[serde(default = "default_true")]
    maximize_bonds: bool,
}

impl McsConfigJson {
    fn into_mcs_config(self) -> chematic_smarts::McsConfig {
        chematic_smarts::McsConfig {
            match_bonds: self.match_bonds,
            min_atoms: self.min_atoms,
            timeout_ms: self.timeout_ms,
            ring_matches_ring_only: self.ring_matches_ring_only,
            complete_rings_only: self.complete_rings_only,
            atom_compare: self.atom_compare.into(),
            bond_compare: self.bond_compare.into(),
            match_chiral_tag: self.match_chiral_tag,
            match_charge: self.match_charge,
            match_isotope: self.match_isotope,
            maximize_bonds: self.maximize_bonds,
        }
    }
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

/// Number of BRICS fragments produced by fragmenting the molecule.
///
/// Returns 1 if no BRICS-breakable bonds exist (whole molecule is one fragment).
#[wasm_bindgen]
pub fn brics_fragment_count(mol: &MolHandle) -> usize {
    chematic_chem::brics_fragments(&mol.inner).len()
}

// run_md_json and coulomb_energy_json removed to reduce WASM bundle size.
// MD simulation is impractical in browser; Coulomb energy was the only user of chematic-ewald.

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

    let products =
        chematic_rxn::enumerate_library_2way(template, scaffolds, building_blocks, &config)
            .map_err(|e| JsValue::from_str(&e.to_string()))?;

    let smiles_list: Vec<String> = products
        .iter()
        .map(|mol| format!("\"{}\"", chematic_smiles::canonical_smiles(mol)))
        .collect();

    Ok(format!("[{}]", smiles_list.join(", ")))
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

/// Compute the tautomer parent with explicit resource limits.
/// Returns `{"smiles":"...","status":"completed"}` (or a structured
/// error) so callers can distinguish a definite result from a budget-limited
/// one.
#[wasm_bindgen]
pub fn tautomer_parent_json(
    mol: &MolHandle,
    max_transforms: usize,
    max_tautomers: usize,
    timeout_ms: Option<u64>,
) -> String {
    if mol.inner.atom_count() > WASM_MAX_ATOMS {
        return format!(
            r#"{{"error":"molecule too large (max {} atoms)"}}"#,
            WASM_MAX_ATOMS
        );
    }
    let mut limits = chematic_chem::TautomerLimits::default();
    limits.max_transforms = max_transforms;
    limits.max_tautomers = max_tautomers;
    limits.timeout_ms = timeout_ms;
    let result = chematic_chem::tautomer_parent(&mol.inner, &limits);
    let status = match result.status {
        chematic_chem::ParentComputationStatus::Completed => "completed",
        chematic_chem::ParentComputationStatus::MaxTransformsReached => "max_transforms_reached",
        chematic_chem::ParentComputationStatus::MaxTautomersReached => "max_tautomers_reached",
        chematic_chem::ParentComputationStatus::TimedOut => "timed_out",
        chematic_chem::ParentComputationStatus::Abstained(_) => "abstained",
        chematic_chem::ParentComputationStatus::InvalidInput(_) => "invalid_input",
        _ => "unknown",
    };
    let smiles = escape_json_string(&chematic_smiles::canonical_smiles(&result.molecule));
    format!(r#"{{"smiles":"{smiles}","status":"{status}"}}"#)
}

/// Compute the canonical tautomer with specific atoms blocked from H-transfer.
///
/// `blocked_atom_indices_json`: JSON array of 0-based atom indices, e.g. `[0, 3]`.
/// Any tautomer move whose donor, bridge, or acceptor is in the blocked set is suppressed.
///
/// Returns canonical SMILES of the result, or `{"error":"..."}` on failure.
/// Out-of-range indices are silently ignored (no effect).
#[wasm_bindgen]
pub fn canonical_tautomer_with_blocked_atoms_json(
    mol: &MolHandle,
    blocked_atom_indices_json: &str,
) -> String {
    if blocked_atom_indices_json.len() > WASM_MAX_JSON_STRING_BYTES {
        return format!(
            r#"{{"error":"blocked_atom_indices_json too large ({} bytes)"}}"#,
            blocked_atom_indices_json.len()
        );
    }
    if mol.inner.atom_count() > WASM_MAX_ATOMS {
        return format!(
            r#"{{"error":"molecule too large (max {} atoms)"}}"#,
            WASM_MAX_ATOMS
        );
    }
    let indices: Vec<u32> = match serde_json::from_str(blocked_atom_indices_json) {
        Ok(v) => v,
        Err(e) => return format!(r#"{{"error":"invalid JSON: {e}"}}"#),
    };
    let blocked_atoms: std::collections::HashSet<chematic_core::AtomIdx> =
        indices.into_iter().map(chematic_core::AtomIdx).collect();
    let config = chematic_chem::TautomerConfig {
        blocked_atoms,
        ..chematic_chem::TautomerConfig::default()
    };
    let result = chematic_chem::canonical_tautomer_with_config(&mol.inner, &config);
    let smi = chematic_smiles::canonical_smiles(&result);
    let escaped = smi.replace('\\', "\\\\").replace('"', "\\\"");
    format!("\"{escaped}\"")
}

/// All enumerated tautomers of `mol` as a JSON array of canonical SMILES strings.
///
/// Example return value: `["Oc1cccc2ccccc12","O=C1C=CC=Cc2ccccc21"]`
#[wasm_bindgen]
pub fn enumerate_tautomers_json(mol: &MolHandle) -> String {
    if mol.inner.atom_count() > WASM_MAX_ATOMS {
        return format!(
            r#"["{{"error":"molecule too large (max {} atoms)"}}"]"#,
            WASM_MAX_ATOMS
        );
    }
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

/// Reconstruct a concrete `Molecule` from a `QueryMolecule` produced by MCS
/// search (atom queries are bare `AtomicNum` primitives; aromaticity is never a
/// per-atom constraint, only carried via aromatic bond queries -- see
/// `build_query`/`molecule_to_query` in `chematic-smarts`) -- shared by every
/// MCS binding below so they can't silently drift apart on how a query result
/// is turned back into a molecule.
fn qmol_to_molecule(qmol: &chematic_smarts::QueryMolecule) -> chematic_core::Molecule {
    use chematic_core::{Atom, AtomIdx, BondOrder, Element, MoleculeBuilder};
    use chematic_smarts::{AtomPrimitive, AtomQuery, BondPrimitive, BondQuery};

    fn extract_atomic_num(q: &AtomQuery) -> Option<u8> {
        match q {
            AtomQuery::Primitive(AtomPrimitive::AtomicNum(n)) => Some(*n),
            AtomQuery::And(lhs, rhs) => extract_atomic_num(lhs).or_else(|| extract_atomic_num(rhs)),
            _ => None,
        }
    }

    // `build_query`/`molecule_to_query` never encode aromaticity as a per-atom
    // constraint (matches RDKit's own `CompareElements` representation) -- it's
    // carried entirely by the aromatic bond queries, so an atom is aromatic here
    // iff at least one of its query bonds is `BondPrimitive::Aromatic`.
    let mut aromatic_atoms = vec![false; qmol.atoms.len()];
    for (atom_idx, neighbors) in qmol.adj.iter().enumerate() {
        for (bond_idx, neighbor_idx) in neighbors {
            if matches!(
                qmol.bonds[*bond_idx].query,
                BondQuery::Primitive(BondPrimitive::Aromatic)
            ) {
                aromatic_atoms[atom_idx] = true;
                aromatic_atoms[*neighbor_idx] = true;
            }
        }
    }

    let mut builder = MoleculeBuilder::new();
    for (idx, qa) in qmol.atoms.iter().enumerate() {
        let elem = extract_atomic_num(&qa.query)
            .and_then(Element::from_atomic_number)
            .unwrap_or(Element::C);
        let mut atom = Atom::new(elem);
        atom.aromatic = aromatic_atoms[idx];
        builder.add_atom(atom);
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
    builder.build()
}

fn parse_mcs_input(
    smiles_json: &str,
    fn_name: &str,
) -> Result<Vec<chematic_core::Molecule>, JsValue> {
    let smiles_list = parse_smiles_json_array(smiles_json)?;
    if smiles_list.len() < 2 {
        return Err(JsValue::from_str(&format!(
            "{fn_name} requires at least 2 SMILES"
        )));
    }
    smiles_list
        .iter()
        .map(|s| {
            let mol = chematic_smiles::parse(s).map_err(|e| JsValue::from_str(&e.to_string()))?;
            enforce_wasm_molecule_size(&mol)?;
            Ok::<_, JsValue>(mol)
        })
        .collect()
}

/// Maximum Common Substructure of a set of molecules, returned as a canonical SMILES string.
///
/// `smiles_json` — a JSON array of at least 2 SMILES strings.
/// Returns the MCS SMILES, or `"null"` when no common substructure was found.
/// Returns a JS error on SMILES parse failure.
#[wasm_bindgen]
pub fn mcs_smiles_json(smiles_json: &str) -> Result<String, JsValue> {
    let mols = parse_mcs_input(smiles_json, "mcs_smiles_json")?;
    let mol_refs: Vec<&chematic_core::Molecule> = mols.iter().collect();
    let qmol = chematic_smarts::find_mcs(&mol_refs);

    if qmol.atoms.is_empty() {
        return Ok("null".to_string());
    }
    Ok(chematic_smiles::canonical_smiles(&qmol_to_molecule(&qmol)))
}

/// MCS with ring-awareness constraints.
///
/// `smiles_json` — JSON array of at least 2 SMILES strings.
/// `ring_matches_ring_only` — ring atoms may only match ring atoms.
/// `complete_rings_only` — partial ring inclusion is removed from the result.
/// Returns the MCS SMILES, or `"null"` when no common substructure was found.
#[wasm_bindgen]
pub fn mcs_smiles_json_with_ring_config(
    smiles_json: &str,
    ring_matches_ring_only: bool,
    complete_rings_only: bool,
) -> Result<String, JsValue> {
    let mols = parse_mcs_input(smiles_json, "mcs_smiles_json_with_ring_config")?;
    let mol_refs: Vec<&chematic_core::Molecule> = mols.iter().collect();
    let config = chematic_smarts::McsConfig {
        ring_matches_ring_only,
        complete_rings_only,
        ..chematic_smarts::McsConfig::default()
    };
    let qmol = chematic_smarts::find_mcs_with_config(&mol_refs, &config);

    if qmol.atoms.is_empty() {
        return Ok("null".to_string());
    }
    Ok(chematic_smiles::canonical_smiles(&qmol_to_molecule(&qmol)))
}

/// Full `McsConfig` + `McsOutcome`-aware MCS search.
///
/// `smiles_json` — JSON array of at least 2 SMILES strings.
/// `config_json` — a JSON object with any subset of `McsConfig`'s fields
/// (camelCase keys, all optional, defaulting to `McsConfig::default()`):
/// `matchBonds`, `minAtoms`, `timeoutMs`, `ringMatchesRingOnly`,
/// `completeRingsOnly`, `atomCompare` (`"elements"` | `"any_heavy_atom"` |
/// `"any"`), `bondCompare` (`"order_or_aromatic"` | `"any"`),
/// `matchChiralTag`, `matchCharge`, `matchIsotope`, `maximizeBonds`.
///
/// Returns a JSON object `{"smiles": string|null, "wasTimedOut": bool}` --
/// `smiles` is `null` when there is no common substructure; `wasTimedOut` is
/// `true` if `timeoutMs` was reached before the search finished exhaustively
/// (the returned `smiles`, if any, is then the best result found so far, not
/// proven optimal).
#[wasm_bindgen]
pub fn mcs_smiles_json_with_config(
    smiles_json: &str,
    config_json: &str,
) -> Result<String, JsValue> {
    let mols = parse_mcs_input(smiles_json, "mcs_smiles_json_with_config")?;
    let config: McsConfigJson =
        serde_json::from_str(config_json).map_err(|e| JsValue::from_str(&e.to_string()))?;
    let config = config.into_mcs_config();

    let mol_refs: Vec<&chematic_core::Molecule> = mols.iter().collect();
    let outcome = chematic_smarts::find_mcs_with_config_checked(&mol_refs, &config);
    let was_timed_out = outcome.was_timed_out();
    let qmol = outcome.into_query();

    let smiles_json_value = if qmol.atoms.is_empty() {
        "null".to_string()
    } else {
        format!(
            "\"{}\"",
            escape_json_string(&chematic_smiles::canonical_smiles(&qmol_to_molecule(&qmol)))
        )
    };
    Ok(format!(
        r#"{{"smiles":{smiles_json_value},"wasTimedOut":{was_timed_out}}}"#
    ))
}

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
    use std::collections::HashSet;

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
        let mapping = &matches[0];

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

/// Single-step retrosynthetic disconnection (issue #91).
///
/// Thin wrapper around [`chematic_rxn::retro::retro_disconnect`] -- applies
/// the same built-in 60-template SMIRKS library and returns identical
/// disconnections (same templates, same precursor sets, same ordering:
/// fewest precursors first) as the Rust and Python (`Mol.retro_disconnect()`)
/// APIs. This function changes nothing about the underlying algorithm; it
/// only serializes the result to JSON.
///
/// `max_results` -- cap on returned disconnections (0 = unlimited).
///
/// `reaction_class` -- filter to a single reaction class, or `""` for all
/// classes. Valid values: `"AmideBond"`, `"Ester"`, `"Ether"`, `"CNBond"`,
/// `"CCBond"`, `"CSBond"`, `"Other"`. An unrecognized non-empty value is a
/// JS error (not silently ignored).
///
/// JSON schema: array of
/// `{"template":str,"reaction_class":str,"precursors":[str,...],"sa_scores":[number,...],"max_sa_score":number}`
/// -- same field names as the Python binding's dict output. Returns `[]`
/// when no template matches the molecule (e.g. it has no disconnectable
/// bond the template library recognizes) -- a valid, non-error result,
/// distinct from the `reaction_class` validation error above.
#[wasm_bindgen]
pub fn retro_disconnect_json(
    mol: &MolHandle,
    max_results: u32,
    reaction_class: &str,
) -> Result<String, JsValue> {
    use chematic_rxn::retro::{DEFAULT_TEMPLATES, RetroClass, RetroTemplate, retro_disconnect};

    enforce_wasm_molecule_size(&mol.inner)?;

    let filter_class: Option<RetroClass> = match reaction_class {
        "" => None,
        "AmideBond" => Some(RetroClass::AmideBond),
        "Ester" => Some(RetroClass::Ester),
        "Ether" => Some(RetroClass::Ether),
        "CNBond" => Some(RetroClass::CNBond),
        "CCBond" => Some(RetroClass::CCBond),
        "CSBond" => Some(RetroClass::CSBond),
        "Other" => Some(RetroClass::Other),
        other => {
            return Err(JsValue::from_str(&format!(
                "unknown reaction_class '{other}'; valid: AmideBond, Ester, Ether, CNBond, CCBond, CSBond, Other"
            )));
        }
    };

    let owned: Vec<RetroTemplate> = DEFAULT_TEMPLATES
        .iter()
        .filter(|t| filter_class.map(|c| c == t.reaction_class).unwrap_or(true))
        .map(|t| RetroTemplate {
            name: t.name,
            smirks: t.smirks,
            reaction_class: t.reaction_class,
        })
        .collect();

    let results = retro_disconnect(&mol.inner, &owned, max_results as usize);

    let parts: Vec<String> = results
        .iter()
        .map(|r| {
            let sa_scores: Vec<f64> = r.precursors.iter().map(chematic_chem::sa_score).collect();
            let max_sa = sa_scores
                .iter()
                .cloned()
                .fold(f64::NEG_INFINITY, f64::max);
            let max_sa = if max_sa.is_finite() { max_sa } else { 0.0 };
            let precursors_json: Vec<String> = r
                .precursor_smiles
                .iter()
                .map(|s| format!("\"{}\"", escape_json_string(s)))
                .collect();
            let sa_scores_json: Vec<String> =
                sa_scores.iter().map(|v| format!("{v}")).collect();
            format!(
                r#"{{"template":"{}","reaction_class":"{}","precursors":[{}],"sa_scores":[{}],"max_sa_score":{}}}"#,
                escape_json_string(&r.template_name),
                r.reaction_class.as_str(),
                precursors_json.join(","),
                sa_scores_json.join(","),
                max_sa,
            )
        })
        .collect();

    Ok(format!("[{}]", parts.join(",")))
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

/// Invert the stereochemistry of a tetrahedral stereocenter (U/D wedge bonds).
///
/// If the atom has no wedge/dash bonds, returns an unchanged copy.
/// Returns error if atom_idx is invalid.
#[wasm_bindgen]
pub fn invert_stereocenter_at(mol: &MolHandle, atom_idx: u32) -> Result<MolHandle, JsValue> {
    let idx = chematic_core::AtomIdx(atom_idx);
    if atom_idx as usize >= mol.inner.atom_count() {
        return Err(JsValue::from_str(&format!(
            "atom_idx {} out of range",
            atom_idx
        )));
    }
    let new_mol = chematic_chem::invert_stereocenter(&mol.inner, idx);
    Ok(MolHandle {
        inner: std::rc::Rc::new(new_mol),
    })
}

// ---------------------------------------------------------------------------
// mol_transforms (3D geometry manipulation)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod mcs_config_tests {
    use super::*;

    #[test]
    fn default_config_matches_bare_mcs_smiles_json() {
        let smiles = r#"["CC(=O)Oc1ccccc1C(=O)O","CC(=O)Nc1ccc(O)cc1"]"#;
        let full = mcs_smiles_json_with_config(smiles, "{}").unwrap();
        let bare = mcs_smiles_json(smiles).unwrap();
        let value: serde_json::Value = serde_json::from_str(&full).unwrap();
        assert_eq!(value["smiles"], serde_json::Value::String(bare));
        assert_eq!(value["wasTimedOut"], false);
    }

    #[test]
    fn no_common_substructure_returns_null_smiles() {
        // Two disconnected, chemically unrelated single-heavy-atom molecules.
        let smiles = r#"["[He]","[Ne]"]"#;
        let result = mcs_smiles_json_with_config(smiles, "{}").unwrap();
        let value: serde_json::Value = serde_json::from_str(&result).unwrap();
        assert_eq!(value["smiles"], serde_json::Value::Null);
        assert_eq!(value["wasTimedOut"], false);
    }

    #[test]
    fn match_charge_true_prevents_charged_neutral_match() {
        // Acetate vs acetic acid: MCS with match_charge=false matches the
        // full 4-heavy-atom core; match_charge=true must reject the
        // charge-differing carboxylate oxygen, shrinking the MCS.
        let smiles = r#"["CC(=O)[O-]","CC(=O)O"]"#;
        let without = mcs_smiles_json_with_config(smiles, r#"{"matchCharge":false}"#).unwrap();
        let with = mcs_smiles_json_with_config(smiles, r#"{"matchCharge":true}"#).unwrap();
        assert_ne!(
            without, with,
            "match_charge=true must change the MCS result"
        );
    }

    #[test]
    fn atom_compare_any_heavy_atom_widens_match() {
        // Benzene vs pyridine: default (element-exact) MCS excludes the N
        // position; any_heavy_atom compare should include it (6-atom ring).
        let smiles = r#"["c1ccccc1","c1ccncc1"]"#;
        let elements =
            mcs_smiles_json_with_config(smiles, r#"{"atomCompare":"elements"}"#).unwrap();
        let any_heavy =
            mcs_smiles_json_with_config(smiles, r#"{"atomCompare":"any_heavy_atom"}"#).unwrap();
        assert_ne!(
            elements, any_heavy,
            "any_heavy_atom compare must change the MCS result vs elements"
        );
    }

    #[test]
    fn timeout_zero_is_reported_as_timed_out() {
        let smiles = r#"["CC(=O)Oc1ccccc1C(=O)O","CC(=O)Nc1ccc(O)cc1"]"#;
        let result = mcs_smiles_json_with_config(smiles, r#"{"timeoutMs":0}"#).unwrap();
        let value: serde_json::Value = serde_json::from_str(&result).unwrap();
        assert_eq!(value["wasTimedOut"], true);
    }

    // Error paths that surface as a real thrown `JsValue` are only safely
    // constructible on the wasm32 target (constructing a `JsValue` natively
    // aborts the process, per wasm-bindgen's own design) -- covered instead
    // by `tests/mcs.test.mjs`, run under the actual wasm runtime via CI's
    // `Test (WASM)` job. What's safe and worth testing natively here is the
    // pure-serde config-parsing layer, which never touches `JsValue`.

    #[test]
    fn invalid_atom_compare_string_rejected_at_parse_time() {
        let err = serde_json::from_str::<McsConfigJson>(r#"{"atomCompare":"bogus"}"#).unwrap_err();
        assert!(err.to_string().contains("atomCompare") || err.to_string().contains("unknown"));
    }

    #[test]
    fn unknown_config_field_rejected_at_parse_time() {
        let err = serde_json::from_str::<McsConfigJson>(r#"{"notAField":true}"#).unwrap_err();
        assert!(err.to_string().contains("notAField") || err.to_string().contains("unknown"));
    }
}
