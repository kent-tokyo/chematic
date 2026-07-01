//! Reaction (SMIRKS/MMP/BRICS/MCS), tautomer, and standardization bindings.

use crate::{
    MolHandle, WASM_MAX_ATOMS, WASM_MAX_BATCH_ITEMS, WASM_MAX_JSON_STRING_BYTES,
    WASM_MAX_SMARTS_MATCHES, enforce_wasm_input_len, enforce_wasm_molecule_size,
    escape_json_string, json_option_string_array, json_option_u8_array, parse_smiles_json_array,
    rgroup_fragment_smiles,
};
use wasm_bindgen::prelude::*;

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
    let smiles_list = parse_smiles_json_array(smiles_json)?;
    if smiles_list.len() < 2 {
        return Err(JsValue::from_str(
            "mcs_smiles_json_with_ring_config requires at least 2 SMILES",
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
    let config = chematic_smarts::McsConfig {
        ring_matches_ring_only,
        complete_rings_only,
        ..chematic_smarts::McsConfig::default()
    };
    let qmol = chematic_smarts::find_mcs_with_config(&mol_refs, &config);

    if qmol.atoms.is_empty() {
        return Ok("null".to_string());
    }

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
