//! Reaction (SMIRKS/MMP/BRICS/MCS), tautomer, and standardization bindings.

use crate::{
    MolHandle, WASM_MAX_ATOMS, WASM_MAX_BATCH_ITEMS, WASM_MAX_INPUT_BYTES,
    WASM_MAX_JSON_STRING_BYTES, WASM_MAX_OUTPUT_BYTES, WASM_MAX_SMARTS_MATCHES,
    enforce_wasm_input_len, enforce_wasm_molecule_size, escape_json_string, json_error,
    json_option_string_array, json_option_u8_array, json_string, parse_smiles_json_array,
    rgroup_fragment_smiles,
};
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

mod molecule_workflows;
pub use molecule_workflows::*;

fn bounded_json_output(output: String) -> Result<String, JsValue> {
    if output.len() > WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str(&format!(
            "JSON output exceeds maximum size ({} > {WASM_MAX_OUTPUT_BYTES} bytes)",
            output.len()
        )));
    }
    Ok(output)
}

fn parse_wasm_reaction(
    reaction_smiles: &str,
) -> Result<chematic_rxn::Reaction, chematic_rxn::RxnError> {
    chematic_rxn::parse_reaction_with_limits(
        reaction_smiles,
        &chematic_rxn::ReactionParseLimits {
            max_input_bytes: WASM_MAX_INPUT_BYTES,
            max_components_per_side: WASM_MAX_BATCH_ITEMS,
            max_atoms_per_molecule: WASM_MAX_ATOMS,
            max_bonds_per_molecule: WASM_MAX_ATOMS.saturating_mul(2),
        },
    )
}

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
            ..chematic_smarts::McsConfig::default()
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
    if smirks.len() > WASM_MAX_INPUT_BYTES || reactants_smiles.len() > WASM_MAX_INPUT_BYTES {
        return Err(JsValue::from_str(&format!(
            "reaction input exceeds maximum size of {WASM_MAX_INPUT_BYTES} bytes"
        )));
    }
    if reactants_smiles.split('|').count() > WASM_MAX_BATCH_ITEMS {
        return Err(JsValue::from_str(&format!(
            "reactant count exceeds maximum ({WASM_MAX_BATCH_ITEMS})"
        )));
    }
    let reactant_mols: Result<Vec<chematic_core::Molecule>, _> = reactants_smiles
        .split('|')
        .map(|s| {
            let mol =
                chematic_smiles::parse(s.trim()).map_err(|e| JsValue::from_str(&e.to_string()))?;
            enforce_wasm_molecule_size(&mol)?;
            Ok::<_, JsValue>(mol)
        })
        .collect();
    let reactant_mols = reactant_mols?;
    let refs: Vec<&chematic_core::Molecule> = reactant_mols.iter().collect();

    let products = chematic_rxn::run_reactants(smirks, &refs)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    if products.len() > WASM_MAX_BATCH_ITEMS
        || products.iter().any(|set| set.len() > WASM_MAX_BATCH_ITEMS)
        || products
            .iter()
            .flatten()
            .any(|mol| mol.atom_count() > WASM_MAX_ATOMS)
    {
        return Err(JsValue::from_str(&format!(
            "reaction result exceeds maximum of {WASM_MAX_BATCH_ITEMS} products per set"
        )));
    }

    let outer: Vec<Vec<String>> = products
        .iter()
        .map(|set| set.iter().map(chematic_smiles::canonical_smiles).collect())
        .collect();
    bounded_json_output(
        serde_json::to_string(&outer).map_err(|e| JsValue::from_str(&e.to_string()))?,
    )
}

#[derive(Serialize)]
struct ReactantAtomJson {
    reactant: usize,
    atom: u32,
}

#[derive(Serialize)]
struct ProductRejectionAtomJson {
    atom_index: u32,
    atom_map: Option<u16>,
    element: &'static str,
    explicit_hydrogen_count: Option<u8>,
    observed_valence: Option<i16>,
    max_allowed_valence: Option<i16>,
}

#[derive(Serialize)]
struct ProductRejectionJson {
    product_index: usize,
    accepted: bool,
    reason: Option<&'static str>,
    atoms: Vec<ProductRejectionAtomJson>,
}

#[derive(Serialize)]
struct CheckedReactionResponse {
    profile: &'static str,
    status: &'static str,
    reason: Option<&'static str>,
    detail: Option<String>,
    products: Vec<Vec<String>>,
    product_atom_sources: Vec<Vec<Vec<Option<ReactantAtomJson>>>>,
    product_template_maps: Vec<Vec<Vec<Option<u16>>>>,
    rejection_diagnostics: Vec<Vec<ProductRejectionJson>>,
    accepted_matches: usize,
    applied_products: usize,
    valence_rejected_matches: usize,
    truncated_matches: bool,
}

impl CheckedReactionResponse {
    fn refusal(
        profile: &'static str,
        status: &'static str,
        reason: &'static str,
        detail: impl Into<String>,
    ) -> Self {
        Self {
            profile,
            status,
            reason: Some(reason),
            detail: Some(detail.into()),
            products: Vec::new(),
            product_atom_sources: Vec::new(),
            product_template_maps: Vec::new(),
            rejection_diagnostics: Vec::new(),
            accepted_matches: 0,
            applied_products: 0,
            valence_rejected_matches: 0,
            truncated_matches: false,
        }
    }

    fn json(&self) -> String {
        serde_json::to_string(self).expect("checked reaction response is JSON-serializable")
    }
}

fn prepare_checked_reaction(
    smirks: &str,
    reactants_smiles: &str,
    rdkit_compat: bool,
    profile: &'static str,
) -> Result<
    (
        Vec<chematic_core::Molecule>,
        Option<std::sync::Arc<chematic_rxn::PreparedReaction>>,
    ),
    String,
> {
    if smirks.len() > WASM_MAX_INPUT_BYTES || reactants_smiles.len() > WASM_MAX_INPUT_BYTES {
        return Err(CheckedReactionResponse::refusal(
            profile,
            "typed_refusal",
            "input_too_large",
            format!("reaction input exceeds {WASM_MAX_INPUT_BYTES} bytes"),
        )
        .json());
    }
    if reactants_smiles.split('|').count() > WASM_MAX_BATCH_ITEMS {
        return Err(CheckedReactionResponse::refusal(
            profile,
            "typed_refusal",
            "reactant_count_limit",
            format!("reactant count exceeds {WASM_MAX_BATCH_ITEMS}"),
        )
        .json());
    }
    let mut reactants = Vec::new();
    for smiles in reactants_smiles.split('|') {
        let mol = chematic_smiles::parse(smiles.trim()).map_err(|error| {
            CheckedReactionResponse::refusal(
                profile,
                "typed_refusal",
                "smiles_parse",
                error.to_string(),
            )
            .json()
        })?;
        if mol.atom_count() > WASM_MAX_ATOMS {
            return Err(CheckedReactionResponse::refusal(
                profile,
                "typed_refusal",
                "reactant_too_large",
                format!("reactant exceeds {WASM_MAX_ATOMS} atoms"),
            )
            .json());
        }
        reactants.push(mol);
    }
    // The RDKit profile prepares its own reading of the template.
    if rdkit_compat {
        return Ok((reactants, None));
    }
    let prepared = chematic_rxn::PreparedReaction::shared(smirks).map_err(|error| {
        CheckedReactionResponse::refusal(
            profile,
            "typed_refusal",
            "smirks_parse",
            error.to_string(),
        )
        .json()
    })?;
    Ok((reactants, Some(prepared)))
}

fn checked_reaction_response(
    profile: &'static str,
    report: chematic_rxn::TracedReactionTransformReport,
    rejection_diagnostics: Vec<Vec<chematic_rxn::ProductSanitizationDiagnostic>>,
) -> CheckedReactionResponse {
    let diagnostics = report.diagnostics;
    let mut result = CheckedReactionResponse {
        profile,
        status: if diagnostics.valence_rejected_matches > 0 || diagnostics.truncated_matches {
            if report.products.is_empty() {
                "typed_refusal"
            } else {
                "partial_products"
            }
        } else if report.products.is_empty() {
            "no_match"
        } else {
            "products"
        },
        reason: if diagnostics.truncated_matches {
            Some("truncated_matches")
        } else if diagnostics.valence_rejected_matches > 0 {
            Some("product_valence")
        } else {
            None
        },
        detail: None,
        products: Vec::new(),
        product_atom_sources: Vec::new(),
        product_template_maps: Vec::new(),
        rejection_diagnostics: rejection_diagnostics
            .into_iter()
            .map(|set| {
                set.into_iter()
                    .map(|product| ProductRejectionJson {
                        product_index: product.product_index,
                        accepted: product.accepted,
                        reason: product.reason.map(|reason| reason.reason_code()),
                        atoms: product
                            .atoms
                            .into_iter()
                            .map(|atom| ProductRejectionAtomJson {
                                atom_index: atom.atom.0,
                                atom_map: atom.template_map,
                                element: atom.element.symbol(),
                                explicit_hydrogen_count: atom.explicit_hydrogen_count,
                                observed_valence: atom.observed_valence,
                                max_allowed_valence: atom.max_allowed_valence,
                            })
                            .collect(),
                    })
                    .collect()
            })
            .collect(),
        accepted_matches: diagnostics.accepted_matches,
        applied_products: diagnostics.applied_products,
        valence_rejected_matches: diagnostics.valence_rejected_matches,
        truncated_matches: diagnostics.truncated_matches,
    };
    if report.products.len() > WASM_MAX_BATCH_ITEMS
        || report.products.iter().any(|set| {
            set.len() > WASM_MAX_BATCH_ITEMS
                || set
                    .iter()
                    .any(|product| product.molecule.atom_count() > WASM_MAX_ATOMS)
        })
    {
        result.status = "typed_refusal";
        result.reason = Some("result_too_large");
        result.detail = Some("reaction products exceed WASM limits".to_string());
    } else {
        for set in report.products {
            let mut smiles_set = Vec::with_capacity(set.len());
            let mut source_set = Vec::with_capacity(set.len());
            let mut map_set = Vec::with_capacity(set.len());
            for product in set {
                let (smiles, output_order) =
                    chematic_smiles::canonical_smiles_with_atom_order(&product.molecule);
                smiles_set.push(smiles);
                source_set.push(
                    output_order
                        .iter()
                        .map(|index| {
                            let source = product.atom_sources[index.0 as usize];
                            source.map(|source| ReactantAtomJson {
                                reactant: source.reactant,
                                atom: source.atom.0,
                            })
                        })
                        .collect(),
                );
                map_set.push(
                    output_order
                        .iter()
                        .map(|index| product.template_maps[index.0 as usize])
                        .collect(),
                );
            }
            result.products.push(smiles_set);
            result.product_atom_sources.push(source_set);
            result.product_template_maps.push(map_set);
        }
    }
    result
}

/// Apply a SMIRKS template with explicit outcome accounting as a JSON string.
///
/// `reactants_smiles` is pipe-separated, as for [`run_reactants`]. With
/// `rdkit_compat = true`, matching and product stereochemistry follow pinned
/// RDKit 2026.03.6; inputs it cannot reproduce return `typed_unsupported` and
/// a stable reason code (`ambiguous_stereo_bond_order`,
/// `ez_reactant_template_semantics`, `chiral_reactant_template_semantics`)
/// instead of an apparently compatible product. Native semantics and the
/// existing [`run_reactants`] API are unchanged. Product graph/origin/
/// template-map parity remains a separate oracle gate. Source indices and
/// product-template map labels are returned alongside the product SMILES;
/// their atom positions match the canonical SMILES parse order.
///
/// The JSON fields mirror Python's `run_smirks_checked`: `status`, `reason`,
/// `detail`, `products`, `accepted_matches`, `applied_products`,
/// `valence_rejected_matches`, and `truncated_matches`, plus `profile`,
/// `product_atom_sources`, and `product_template_maps`.
#[wasm_bindgen]
pub fn run_reactants_checked(smirks: &str, reactants_smiles: &str, rdkit_compat: bool) -> String {
    let profile = if rdkit_compat {
        "rdkit-2026.03.6"
    } else {
        "native"
    };
    let (reactants, prepared) =
        match prepare_checked_reaction(smirks, reactants_smiles, rdkit_compat, profile) {
            Ok(prepared) => prepared,
            Err(refusal_json) => return refusal_json,
        };
    let refs: Vec<_> = reactants.iter().collect();
    let limits = chematic_rxn::ReactionTransformLimits {
        max_matches: WASM_MAX_BATCH_ITEMS,
    };
    let (report, rejection_diagnostics) = match &prepared {
        Some(prepared) => match prepared.run_reactants_traced_with_diagnostics(&refs, &limits) {
            Ok(report) => (report, Vec::new()),
            Err(error) => {
                let reason = match &error {
                    chematic_rxn::TransformError::ResourceLimit { .. } => "resource_limit",
                    chematic_rxn::TransformError::ReactantCountMismatch { .. } => {
                        "reactant_count_mismatch"
                    }
                    chematic_rxn::TransformError::SmirksParse(_) => "smirks_parse",
                };
                return CheckedReactionResponse::refusal(
                    profile,
                    "typed_refusal",
                    reason,
                    error.to_string(),
                )
                .json();
            }
        },
        None => match chematic_rxn::run_reactants_traced_rdkit_2026_03_6_detailed(
            smirks, &refs, &limits,
        ) {
            Ok(chematic_rxn::RdkitDetailedProfileOutcome::Report(report)) => {
                (report.report, report.rejection_diagnostics)
            }
            Ok(chematic_rxn::RdkitDetailedProfileOutcome::Unsupported(unsupported)) => {
                return CheckedReactionResponse::refusal(
                    profile,
                    "typed_unsupported",
                    unsupported.reason_code(),
                    "RDKit 2026.03.6 reaction semantics cannot be reproduced for this input",
                )
                .json();
            }
            Err(error) => {
                let reason = match &error {
                    chematic_rxn::TransformError::ResourceLimit { .. } => "resource_limit",
                    chematic_rxn::TransformError::ReactantCountMismatch { .. } => {
                        "reactant_count_mismatch"
                    }
                    chematic_rxn::TransformError::SmirksParse(_) => "smirks_parse",
                };
                return CheckedReactionResponse::refusal(
                    profile,
                    "typed_refusal",
                    reason,
                    error.to_string(),
                )
                .json();
            }
        },
    };
    let mut result = checked_reaction_response(profile, report, rejection_diagnostics);
    let json = result.json();
    if json.len() > WASM_MAX_OUTPUT_BYTES {
        result.status = "typed_refusal";
        result.reason = Some("output_too_large");
        result.detail = Some("JSON output exceeds WASM limit".to_string());
        result.products.clear();
        result.product_atom_sources.clear();
        result.product_template_maps.clear();
        result.rejection_diagnostics.clear();
        result.json()
    } else {
        json
    }
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
    if template.len() > WASM_MAX_INPUT_BYTES
        || scaffolds_smiles.len() > WASM_MAX_INPUT_BYTES
        || building_blocks_smiles.len() > WASM_MAX_INPUT_BYTES
    {
        return Err(JsValue::from_str(&format!(
            "library input exceeds maximum size of {WASM_MAX_INPUT_BYTES} bytes"
        )));
    }
    if scaffolds_smiles.split('|').count() > WASM_MAX_BATCH_ITEMS {
        return Err(JsValue::from_str(&format!(
            "scaffold count exceeds maximum ({WASM_MAX_BATCH_ITEMS})"
        )));
    }
    if building_blocks_smiles.split('|').count() > WASM_MAX_BATCH_ITEMS {
        return Err(JsValue::from_str(&format!(
            "building-block count exceeds maximum ({WASM_MAX_BATCH_ITEMS})"
        )));
    }
    let scaffolds: Result<Vec<chematic_core::Molecule>, _> = scaffolds_smiles
        .split('|')
        .map(|s| {
            let mol =
                chematic_smiles::parse(s.trim()).map_err(|e| JsValue::from_str(&e.to_string()))?;
            enforce_wasm_molecule_size(&mol)?;
            Ok::<_, JsValue>(mol)
        })
        .collect();
    let scaffolds = scaffolds?;

    let building_blocks: Result<Vec<chematic_core::Molecule>, _> = building_blocks_smiles
        .split('|')
        .map(|s| {
            let mol =
                chematic_smiles::parse(s.trim()).map_err(|e| JsValue::from_str(&e.to_string()))?;
            enforce_wasm_molecule_size(&mol)?;
            Ok::<_, JsValue>(mol)
        })
        .collect();
    let building_blocks = building_blocks?;

    let config = chematic_rxn::LibraryConfig {
        skip_failures: true,
        max_size: Some(WASM_MAX_BATCH_ITEMS), // Limit enumeration size
    };

    let products =
        chematic_rxn::enumerate_library_2way(template, scaffolds, building_blocks, &config)
            .map_err(|e| JsValue::from_str(&e.to_string()))?;
    if products.iter().any(|mol| mol.atom_count() > WASM_MAX_ATOMS) {
        return Err(JsValue::from_str(&format!(
            "library product exceeds maximum atom count ({WASM_MAX_ATOMS})"
        )));
    }

    let smiles_list: Vec<String> = products
        .iter()
        .map(chematic_smiles::canonical_smiles)
        .collect();

    bounded_json_output(
        serde_json::to_string(&smiles_list).map_err(|e| JsValue::from_str(&e.to_string()))?,
    )
}

/// Render a reaction SMILES string (e.g. `"CC(=O)O.CCO>>CC(=O)OCC.O"`) as a
/// single SVG showing reactants → products with `+` separators.
///
/// Returns a self-contained SVG string.  Returns a JS error on invalid input.
#[wasm_bindgen]
pub fn depict_reaction_svg(rxn_smiles: &str) -> Result<String, JsValue> {
    let rxn = parse_wasm_reaction(rxn_smiles).map_err(|e| JsValue::from_str(&e.to_string()))?;

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
        frags.push((
            cursor,
            chematic_mol::stereo_depiction::with_stereo_depiction(mol, |m, layout| {
                chematic_depict::render_svg_opts(m, layout, &opts)
            }),
        ));
        cursor += MOL_W;
    }

    seps.push((cursor + SEP_ARROW / 2, "→"));
    cursor += SEP_ARROW;

    for (i, mol) in rxn.products.iter().enumerate() {
        if i > 0 {
            seps.push((cursor + SEP_PLUS / 2, "+"));
            cursor += SEP_PLUS;
        }
        frags.push((
            cursor,
            chematic_mol::stereo_depiction::with_stereo_depiction(mol, |m, layout| {
                chematic_depict::render_svg_opts(m, layout, &opts)
            }),
        ));
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

    enforce_wasm_input_len("core_smarts", core_smarts)?;
    let query = chematic_smarts::parse_smarts(core_smarts)
        .map_err(|e| JsValue::from_str(&format!("{e:?}")))?;
    if query.atoms.len() > WASM_MAX_ATOMS {
        return Err(JsValue::from_str(&format!(
            "core_smarts exceeds maximum query atom count ({WASM_MAX_ATOMS})"
        )));
    }

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
    let rxn = parse_wasm_reaction(rxn_smiles).map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(chematic_rxn::write_reaction(&rxn))
}

/// Check whether a reaction SMILES matches a reaction SMARTS query.
///
/// The middle section of the query supports agent alternatives separated by
/// `|`. This source-level API remains bounded and returns a typed JS error for
/// invalid input; the generated Node artifact is updated separately when the
/// wasm-bindgen toolchain is available.
#[wasm_bindgen]
pub fn reaction_smarts_match(smarts: &str, reaction_smiles: &str) -> Result<bool, JsValue> {
    enforce_wasm_input_len("reaction SMARTS", smarts)?;
    enforce_wasm_input_len("reaction SMILES", reaction_smiles)?;
    let query = chematic_rxn::parse_reaction_query(smarts)
        .map_err(|e| JsValue::from_str(&format!("invalid reaction SMARTS: {e}")))?;
    let reaction = parse_wasm_reaction(reaction_smiles)
        .map_err(|e| JsValue::from_str(&format!("invalid reaction SMILES: {e}")))?;
    Ok(chematic_rxn::has_reaction_substructure_match(
        &reaction, &query,
    ))
}

/// Validate and normalize a typed Markush/polymer semantic model JSON.
#[wasm_bindgen]
pub fn semantic_model_json(model_json: &str) -> Result<String, JsValue> {
    let value: serde_json::Value = serde_json::from_str(model_json)
        .map_err(|e| JsValue::from_str(&format!("invalid semantic JSON: {e}")))?;
    let model = chematic_mol::SemanticModel::from_json(&value)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    serde_json::to_string(&model.to_json()).map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Apply an explicit Markush selection command to a semantic model JSON.
#[wasm_bindgen]
pub fn semantic_apply_json_command(
    model_json: &str,
    command_json: &str,
) -> Result<String, JsValue> {
    let model_value: serde_json::Value = serde_json::from_str(model_json)
        .map_err(|e| JsValue::from_str(&format!("invalid semantic JSON: {e}")))?;
    let command: serde_json::Value = serde_json::from_str(command_json)
        .map_err(|e| JsValue::from_str(&format!("invalid semantic command JSON: {e}")))?;
    let model = chematic_mol::SemanticModel::from_json(&model_value)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let next = model
        .apply_json_command(&command)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    serde_json::to_string(&next.to_json()).map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Expand a validated semantic model against a base SMILES and return the
/// expanded graph plus source-to-expanded atom mapping.
#[wasm_bindgen]
pub fn semantic_expand_json(base_smiles: &str, model_json: &str) -> Result<String, JsValue> {
    semantic_expand_json_with_limits(
        base_smiles,
        model_json,
        chematic_mol::SemanticExpansionLimits::default().max_atoms as u32,
        chematic_mol::SemanticExpansionLimits::default().max_repeat_count,
    )
}

/// Expand a validated semantic model with explicit finite atom/repeat budgets.
#[wasm_bindgen]
pub fn semantic_expand_json_with_limits(
    base_smiles: &str,
    model_json: &str,
    max_atoms: u32,
    max_repeat_count: u32,
) -> Result<String, JsValue> {
    let value: serde_json::Value = serde_json::from_str(model_json)
        .map_err(|e| JsValue::from_str(&format!("invalid semantic JSON: {e}")))?;
    let model = chematic_mol::SemanticModel::from_json(&value)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let base =
        chematic_smiles::parse(base_smiles).map_err(|e| JsValue::from_str(&e.to_string()))?;
    let expanded = model
        .expand_with_limits(
            &base,
            &chematic_mol::SemanticExpansionLimits {
                max_atoms: max_atoms as usize,
                max_repeat_count,
            },
        )
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    serde_json::to_string(&expanded.to_json()).map_err(|e| JsValue::from_str(&e.to_string()))
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
/// `max_results` -- cap on returned disconnections (0 = the WASM safety cap).
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
    if max_results as usize > WASM_MAX_BATCH_ITEMS {
        return Err(JsValue::from_str(&format!(
            "retro result count exceeds maximum ({WASM_MAX_BATCH_ITEMS})"
        )));
    }
    let effective_max_results = if max_results == 0 {
        WASM_MAX_BATCH_ITEMS
    } else {
        max_results as usize
    };

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

    let results = retro_disconnect(&mol.inner, &owned, effective_max_results);

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
    let rxn = match parse_wasm_reaction(reaction_smiles) {
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
    if smiles.len() > WASM_MAX_INPUT_BYTES {
        return format!(
            "error:SMILES exceeds maximum input size ({} > {WASM_MAX_INPUT_BYTES} bytes)",
            smiles.len()
        );
    }
    let mol = match chematic_smiles::parse(smiles) {
        Ok(m) => m,
        Err(e) => return format!("error:{e}"),
    };
    if mol.atom_count() > WASM_MAX_ATOMS {
        return format!("error:SMILES exceeds maximum atom count ({WASM_MAX_ATOMS})");
    }
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
    if smiles.len() > WASM_MAX_INPUT_BYTES {
        return format!(
            "error:SMILES exceeds maximum input size ({} > {WASM_MAX_INPUT_BYTES} bytes)",
            smiles.len()
        );
    }
    let mol = match chematic_smiles::parse(smiles) {
        Ok(m) => m,
        Err(e) => return format!("error:{e}"),
    };
    if mol.atom_count() > WASM_MAX_ATOMS {
        return format!("error:SMILES exceeds maximum atom count ({WASM_MAX_ATOMS})");
    }
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
    let rxn = match parse_wasm_reaction(reaction_smiles) {
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

/// Analyze a rich reaction document's explicit atom/isotope inventory and
/// formal charges. The returned JSON includes evidence scope, per-step
/// diagnostics, and a status that does not imply chemical completeness.
/// Returns `error:<msg>` when the document or a component is invalid.
#[wasm_bindgen]
pub fn stoichiometry_report_json(document_json: &str) -> String {
    let document = match chematic_rxn::ReactionDocument::from_json_str(document_json) {
        Ok(document) => document,
        Err(error) => return format!("error:{error}"),
    };
    match chematic_rxn::analyze_reaction_document(&document).and_then(|report| {
        serde_json::to_string(&report).map_err(|error| chematic_rxn::StoichiometryError::Parse {
            path: "/report".into(),
            message: error.to_string(),
        })
    }) {
        Ok(json) => json,
        Err(error) => format!("error:{error}"),
    }
}

// ---------------------------------------------------------------------------
// Nearest-neighbour similarity search
// ---------------------------------------------------------------------------

/// Invert a tetrahedral stereocenter and synchronize incident U/D wedge bonds.
///
/// Recorded chirality is inverted even without a wedge; when a wedge/hash is
/// present its depiction marker is inverted too.
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

#[cfg(test)]
include!("mol_reactions_tests.rs");
