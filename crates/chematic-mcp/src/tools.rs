//! MCP tool implementations for chematic.
//!
//! Each `tool_*` function computes its chemistry exactly once and returns a
//! bare JSON payload (`Value`) — never a protocol envelope. The
//! presentation layer (`server.rs`) wraps that payload differently per
//! protocol era (legacy `content`-only vs. modern `content` +
//! `structuredContent`), so the underlying computation is never run twice
//! and never duplicated between eras.

#![forbid(unsafe_code)]

use std::collections::{HashMap, HashSet, VecDeque};

use chematic_3d::{generate_and_minimize_dreiding, write_xyz};
use chematic_core::{Atom, AtomIdx, BondOrder, Element, MoleculeBuilder};
use chematic_fp::{BitVec2048, ecfp4, tanimoto_ecfp4};
use chematic_inchi::inchi;
use chematic_mol::{parse_moljson, write_cml, write_moljson};
use chematic_smarts::{
    AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, McsConfig, find_matches,
    find_mcs_with_config, parse_smarts,
};
use serde_json::{Value, json};

use chematic_chem::{
    admet_profile, boiled_egg, brenk_matches, brenk_passes, brics_bonds, exact_mass, hba_count,
    hbd_count, heavy_atom_count, lipinski_passes, logp_crippen, molecular_weight, pains_matches,
    pains_passes, qed, rotatable_bond_count, sa_score, tpsa,
};

// ── tool-call error taxonomy ───────────────────────────────────────────────

/// Why a `tools/call` failed, split along the line section 9 of the RC
/// implementation draws between transport-level and domain-level failure:
///
/// - `InvalidArgs`: the arguments themselves are missing, wrong-typed, or
///   otherwise violate the tool's declared `inputSchema` (including an
///   unknown tool name). This is a request-shape problem the *server*
///   rejects before running any chemistry.
/// - `Domain`: the arguments were well-formed, but the requested chemistry
///   computation failed (unparseable SMILES/SMARTS, a molecule too large or
///   disconnected for a bounded algorithm, a failed PubChem lookup, ...).
///   The tool call itself succeeded as an RPC; the result reports failure.
///
/// Legacy-era wire behavior is unaffected by this split — both variants
/// carry the exact same message text `call_tool` returned before this
/// refactor (see `legacy_message`), preserving the current
/// `{"code":-32000,...}` shape byte-for-byte. Only the modern era uses the
/// split: `InvalidArgs` becomes `-32602`, `Domain` becomes a successful
/// `CallToolResult` with `isError: true`.
#[derive(Debug, Clone)]
pub enum ToolCallError {
    InvalidArgs(String),
    Domain {
        code: &'static str,
        message: String,
        details: Value,
    },
}

impl ToolCallError {
    fn invalid_args(msg: impl Into<String>) -> Self {
        ToolCallError::InvalidArgs(msg.into())
    }

    fn domain(code: &'static str, msg: impl Into<String>) -> Self {
        ToolCallError::Domain {
            code,
            message: msg.into(),
            details: json!({}),
        }
    }

    /// The message text this tool call would have produced before the
    /// `ToolCallError` split existed — used verbatim by the legacy
    /// presentation layer so wire output is unchanged.
    pub fn legacy_message(&self) -> &str {
        match self {
            ToolCallError::InvalidArgs(m) => m,
            ToolCallError::Domain { message, .. } => message,
        }
    }

    /// `true` for argument/schema-shape problems (maps to `-32602` in the
    /// modern era); `false` for domain/chemistry failures (maps to a
    /// successful result with `isError: true`).
    pub fn is_invalid_args(&self) -> bool {
        matches!(self, ToolCallError::InvalidArgs(_))
    }

    /// Machine-readable `{code, message, details}` object for
    /// `structuredContent.error` in the modern era.
    pub fn to_structured_error(&self) -> Value {
        match self {
            ToolCallError::InvalidArgs(m) => json!({
                "code": "INVALID_ARGUMENTS",
                "message": m,
                "details": {}
            }),
            ToolCallError::Domain {
                code,
                message,
                details,
            } => json!({ "code": code, "message": message, "details": details }),
        }
    }
}

// ── helpers ──────────────────────────────────────────────────────────────────

fn get_str<'a>(args: &'a Value, key: &str) -> Result<&'a str, ToolCallError> {
    args.get(key)
        .and_then(|v| v.as_str())
        .ok_or_else(|| ToolCallError::invalid_args(format!("Missing required argument: {key}")))
}

/// Parse a SMILES tool argument. A parse failure is a *domain* error
/// (`INVALID_SMILES`), not an argument-shape error — the argument was a
/// well-formed string, it just doesn't describe a valid molecule.
fn parse_mol_arg(smiles: &str) -> Result<chematic_core::Molecule, ToolCallError> {
    chematic_smiles::parse(smiles).map_err(|e| {
        ToolCallError::domain("INVALID_SMILES", format!("Invalid SMILES '{smiles}': {e}"))
    })
}

fn round3(x: f64) -> f64 {
    (x * 1000.0).round() / 1000.0
}

fn round2(x: f64) -> f64 {
    (x * 100.0).round() / 100.0
}

/// Convert a BitVec2048 to a lowercase hex string (256 bytes = 2048 bits).
fn bitvec_to_hex(fp: &BitVec2048) -> String {
    let mut bytes = [0u8; 256];
    for i in 0..2048_usize {
        if fp.get(i) {
            bytes[i / 8] |= 1u8 << (i % 8);
        }
    }
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

/// Reconstruct a concrete `Molecule` from an MCS `QueryMolecule`.
fn qmol_to_molecule(qmol: &chematic_smarts::QueryMolecule) -> chematic_core::Molecule {
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
    builder.build()
}

// ── retrosynthesis helpers ────────────────────────────────────────────────────

/// BFS from `start`, skipping the directed edge (`excl_a` → `excl_b`) in both
/// directions.  Returns the set of atoms reachable without crossing that bond.
fn atoms_reachable_excl_bond(
    mol: &chematic_core::Molecule,
    start: AtomIdx,
    excl_a: AtomIdx,
    excl_b: AtomIdx,
) -> HashSet<AtomIdx> {
    let mut visited: HashSet<AtomIdx> = HashSet::new();
    let mut queue: VecDeque<AtomIdx> = VecDeque::new();
    queue.push_back(start);
    visited.insert(start);
    while let Some(curr) = queue.pop_front() {
        for (nb, _) in mol.neighbors(curr) {
            if (curr == excl_a && nb == excl_b) || (curr == excl_b && nb == excl_a) {
                continue;
            }
            if visited.insert(nb) {
                queue.push_back(nb);
            }
        }
    }
    visited
}

/// Build a sub-molecule from a subset of atoms, preserving all internal bonds.
/// `hydrogen_count` is cleared so implicit Hs are re-derived from standard valence.
fn build_submol(
    mol: &chematic_core::Molecule,
    atom_set: &HashSet<AtomIdx>,
) -> chematic_core::Molecule {
    let mut builder = MoleculeBuilder::new();
    let mut old_to_new: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    let mut sorted: Vec<AtomIdx> = atom_set.iter().cloned().collect();
    sorted.sort();

    for &old_idx in &sorted {
        let mut atom = mol.atom(old_idx).clone();
        // Reset explicit H count so the SMILES writer infers Hs from valence;
        // otherwise the cut atom would retain a stale bracket-H count.
        atom.hydrogen_count = None;
        atom.cip_code = None; // stereo may be invalid after the bond is removed
        let new_idx = builder.add_atom(atom);
        old_to_new.insert(old_idx, new_idx);
    }

    for &old_a in &sorted {
        for (old_b, bidx) in mol.neighbors(old_a) {
            if old_a < old_b && atom_set.contains(&old_b) {
                let bond = mol.bond(bidx);
                let new_a = old_to_new[&old_a];
                let new_b = old_to_new[&old_b];
                let _ = builder.add_bond(new_a, new_b, bond.order);
            }
        }
    }

    builder.build()
}

/// Return the number of connected components in `mol`.
fn component_count(mol: &chematic_core::Molecule) -> usize {
    if mol.atom_count() == 0 {
        return 0;
    }
    let mut visited: HashSet<AtomIdx> = HashSet::new();
    let mut count = 0;
    for (start, _) in mol.atoms() {
        if visited.contains(&start) {
            continue;
        }
        count += 1;
        let mut queue: VecDeque<AtomIdx> = VecDeque::new();
        queue.push_back(start);
        visited.insert(start);
        while let Some(curr) = queue.pop_front() {
            for (nb, _) in mol.neighbors(curr) {
                if visited.insert(nb) {
                    queue.push_back(nb);
                }
            }
        }
    }
    count
}

// ── tool list schema ──────────────────────────────────────────────────────────

/// Number of tools this server exposes. A single source of truth so
/// registry-invariant tests (`exactly 20 unique names`, negative controls
/// for "duplicate a tool name" / "remove a tool") never hardcode `20` twice.
pub const TOOL_COUNT: usize = 20;

pub fn list_tools() -> Value {
    json!({ "tools": [
        {
            "name": "parse_smiles",
            "description": "Parse a SMILES string and return basic molecular information (atom count, bond count, molecular weight).",
            "inputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "smiles": { "type": "string", "description": "SMILES string to parse", "minLength": 1, "maxLength": 100000 }
                },
                "required": ["smiles"]
            },
            "outputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "valid": { "type": "boolean" },
                    "atoms": { "type": "integer", "minimum": 0 },
                    "bonds": { "type": "integer", "minimum": 0 },
                    "mol_weight": { "type": "number", "minimum": 0 },
                    "smiles": { "type": "string" }
                },
                "required": ["valid", "atoms", "bonds", "mol_weight", "smiles"]
            }
        },
        {
            "name": "calc_properties",
            "description": "Calculate molecular properties: MW, exact mass, LogP (Crippen), TPSA, HBD, HBA, rotatable bonds, heavy atom count, and QED drug-likeness.",
            "inputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "smiles": { "type": "string", "description": "SMILES string", "minLength": 1, "maxLength": 100000 }
                },
                "required": ["smiles"]
            },
            "outputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "mw": { "type": "number", "minimum": 0 },
                    "exact_mass": { "type": "number", "minimum": 0 },
                    "logp": { "type": "number" },
                    "tpsa": { "type": "number", "minimum": 0 },
                    "hbd": { "type": "integer", "minimum": 0 },
                    "hba": { "type": "integer", "minimum": 0 },
                    "rotatable_bonds": { "type": "integer", "minimum": 0 },
                    "heavy_atom_count": { "type": "integer", "minimum": 0 },
                    "qed": { "type": "number", "minimum": 0, "maximum": 1 }
                },
                "required": ["mw", "exact_mass", "logp", "tpsa", "hbd", "hba", "rotatable_bonds", "heavy_atom_count", "qed"]
            }
        },
        {
            "name": "ecfp4",
            "description": "Compute the ECFP4 (Morgan radius-2) circular fingerprint as a 2048-bit hex string, plus popcount.",
            "inputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "smiles": { "type": "string", "description": "SMILES string", "minLength": 1, "maxLength": 100000 }
                },
                "required": ["smiles"]
            },
            "outputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "fingerprint": { "type": "string", "minLength": 512, "maxLength": 512, "description": "512-char lowercase hex string (256 bytes = 2048 bits)" },
                    "popcount": { "type": "integer", "minimum": 0, "maximum": 2048 }
                },
                "required": ["fingerprint", "popcount"]
            }
        },
        {
            "name": "tanimoto",
            "description": "Compute the Tanimoto (Jaccard) similarity between two molecules using ECFP4 fingerprints.",
            "inputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "smiles1": { "type": "string", "description": "First molecule SMILES", "minLength": 1, "maxLength": 100000 },
                    "smiles2": { "type": "string", "description": "Second molecule SMILES", "minLength": 1, "maxLength": 100000 }
                },
                "required": ["smiles1", "smiles2"]
            },
            "outputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "similarity": { "type": "number", "minimum": 0, "maximum": 1 },
                    "similarity_percent": { "type": "number", "minimum": 0, "maximum": 100 }
                },
                "required": ["similarity", "similarity_percent"]
            }
        },
        {
            "name": "smarts_match",
            "description": "Perform SMARTS substructure search and return whether the pattern matches, match count, and atom index maps.",
            "inputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "smarts": { "type": "string", "description": "SMARTS pattern", "minLength": 1, "maxLength": 100000 },
                    "smiles": { "type": "string", "description": "Molecule SMILES to search in", "minLength": 1, "maxLength": 100000 }
                },
                "required": ["smarts", "smiles"]
            },
            "outputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "matches": { "type": "boolean" },
                    "match_count": { "type": "integer", "minimum": 0 },
                    "atom_maps": {
                        "type": "array",
                        "items": { "type": "array", "items": { "type": "integer", "minimum": 0 } }
                    }
                },
                "required": ["matches", "match_count", "atom_maps"]
            }
        },
        {
            "name": "canonical_smiles",
            "description": "Return the canonical SMILES representation of a molecule.",
            "inputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "smiles": { "type": "string", "description": "Input SMILES string", "minLength": 1, "maxLength": 100000 }
                },
                "required": ["smiles"]
            },
            "outputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "canonical": { "type": "string" }
                },
                "required": ["canonical"]
            }
        },
        {
            "name": "find_mcs",
            "description": "Find the maximum common substructure (MCS) across a list of molecules. Returns the MCS as a canonical SMILES string.",
            "inputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "smiles_list": {
                        "type": "array",
                        "items": { "type": "string", "minLength": 1, "maxLength": 100000 },
                        "description": "List of SMILES strings (2–20 molecules)",
                        "minItems": 2,
                        "maxItems": 20
                    }
                },
                "required": ["smiles_list"]
            },
            "outputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "mcs": { "type": ["string", "null"] },
                    "atom_count": { "type": "integer", "minimum": 0 },
                    "bond_count": { "type": "integer", "minimum": 0 }
                },
                "required": ["mcs", "atom_count", "bond_count"]
            }
        },
        {
            "name": "generate_3d",
            "description": "Generate 3D coordinates for a molecule using rule-based placement and DREIDING force-field minimization. Returns XYZ format.",
            "inputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "smiles": { "type": "string", "description": "SMILES string", "minLength": 1, "maxLength": 100000 }
                },
                "required": ["smiles"]
            },
            "outputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "xyz": { "type": "string" },
                    "atom_count": { "type": "integer", "minimum": 0 }
                },
                "required": ["xyz", "atom_count"]
            }
        },
        {
            "name": "pains_check",
            "description": "Check whether a molecule contains Pan-Assay Interference Compounds (PAINS) structural alerts. PAINS compounds often produce false positives in high-throughput screening.",
            "inputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "smiles": { "type": "string", "description": "SMILES string", "minLength": 1, "maxLength": 100000 }
                },
                "required": ["smiles"]
            },
            "outputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "passes": { "type": "boolean" },
                    "alert_count": { "type": "integer", "minimum": 0 },
                    "alerts": { "type": "array", "items": { "type": "string" } }
                },
                "required": ["passes", "alert_count", "alerts"]
            }
        },
        {
            "name": "brenk_check",
            "description": "Check whether a molecule contains Brenk structural alerts (unwanted functional groups associated with toxicity, metabolic instability, or undesirable reactivity).",
            "inputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "smiles": { "type": "string", "description": "SMILES string", "minLength": 1, "maxLength": 100000 }
                },
                "required": ["smiles"]
            },
            "outputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "passes": { "type": "boolean" },
                    "alert_count": { "type": "integer", "minimum": 0 },
                    "alerts": { "type": "array", "items": { "type": "string" } }
                },
                "required": ["passes", "alert_count", "alerts"]
            }
        },
        {
            "name": "sa_score",
            "description": "Estimate synthetic accessibility (SA Score, Ertl & Schuffenhauer 2009). Returns a score from 1 (easy to synthesize) to 10 (very difficult). Drug-like molecules typically score 2–4.",
            "inputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "smiles": { "type": "string", "description": "SMILES string", "minLength": 1, "maxLength": 100000 }
                },
                "required": ["smiles"]
            },
            "outputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "sa_score": { "type": "number", "minimum": 1, "maximum": 10 },
                    "easy_to_synthesize": { "type": "boolean" },
                    "note": { "type": "string" }
                },
                "required": ["sa_score", "easy_to_synthesize", "note"]
            }
        },
        {
            "name": "admet_profile",
            "description": "Compute a full ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity) profile including BBB penetration, Caco-2 permeability, hERG risk, CYP3A4 inhibition risk, AMES mutagenicity risk, plasma protein binding, and hepatic clearance class.",
            "inputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "smiles": { "type": "string", "description": "SMILES string", "minLength": 1, "maxLength": 100000 }
                },
                "required": ["smiles"]
            },
            "outputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "bbb_score": { "type": "number" },
                    "bbb_passes": { "type": "boolean" },
                    "caco2_logp": { "type": "number" },
                    "herg_risk": { "type": "number" },
                    "cyp3a4_risk": { "type": "number" },
                    "pka_acid": { "type": ["number", "null"] },
                    "pka_base": { "type": ["number", "null"] },
                    "esol_logs": { "type": "number" },
                    "logd_74": { "type": "number" },
                    "mw": { "type": "number", "minimum": 0 },
                    "logp": { "type": "number" },
                    "tpsa": { "type": "number", "minimum": 0 },
                    "hbd": { "type": "integer", "minimum": 0 },
                    "hba": { "type": "integer", "minimum": 0 },
                    "rotatable_bonds": { "type": "integer", "minimum": 0 },
                    "ames_risk": { "type": "number" },
                    "ppb_percent": { "type": "number" },
                    "clearance_class": { "type": "string" }
                },
                "required": [
                    "bbb_score", "bbb_passes", "caco2_logp", "herg_risk", "cyp3a4_risk",
                    "pka_acid", "pka_base", "esol_logs", "logd_74", "mw", "logp", "tpsa",
                    "hbd", "hba", "rotatable_bonds", "ames_risk", "ppb_percent", "clearance_class"
                ]
            }
        },
        {
            "name": "boiled_egg",
            "description": "Predict passive gastrointestinal (GI) absorption and blood-brain barrier (BBB) penetration using the BOILED-Egg method (Daina & Zoete 2016). Uses LogP and TPSA thresholds to classify molecules into the egg-white (GI absorbed) and egg-yolk (BBB penetrant) zones.",
            "inputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "smiles": { "type": "string", "description": "SMILES string", "minLength": 1, "maxLength": 100000 }
                },
                "required": ["smiles"]
            },
            "outputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "gi_absorbed": { "type": "boolean" },
                    "bbb_penetrant": { "type": "boolean" },
                    "logp": { "type": "number" },
                    "tpsa": { "type": "number", "minimum": 0 },
                    "method": { "type": "string" },
                    "thresholds": {
                        "type": "object",
                        "additionalProperties": false,
                        "properties": {
                            "gi_white": { "type": "string" },
                            "bbb_yolk": { "type": "string" }
                        },
                        "required": ["gi_white", "bbb_yolk"]
                    }
                },
                "required": ["gi_absorbed", "bbb_penetrant", "logp", "tpsa", "method", "thresholds"]
            }
        },
        {
            "name": "lipinski_check",
            "description": "Check Lipinski's Rule of Five for oral drug-likeness (MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10). Returns whether the molecule passes and individual property values.",
            "inputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "smiles": { "type": "string", "description": "SMILES string", "minLength": 1, "maxLength": 100000 }
                },
                "required": ["smiles"]
            },
            "outputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "passes": { "type": "boolean" },
                    "mw": { "type": "number", "minimum": 0 },
                    "logp": { "type": "number" },
                    "hbd": { "type": "integer", "minimum": 0 },
                    "hba": { "type": "integer", "minimum": 0 },
                    "rules": {
                        "type": "object",
                        "additionalProperties": false,
                        "properties": {
                            "mw_le_500": { "type": "boolean" },
                            "logp_le_5": { "type": "boolean" },
                            "hbd_le_5": { "type": "boolean" },
                            "hba_le_10": { "type": "boolean" }
                        },
                        "required": ["mw_le_500", "logp_le_5", "hbd_le_5", "hba_le_10"]
                    }
                },
                "required": ["passes", "mw", "logp", "hbd", "hba", "rules"]
            }
        },
        {
            "name": "name_to_smiles",
            "description": "Convert a chemical name (IUPAC, trivial, or trade name) to an isomeric SMILES string via the PubChem REST API. Requires internet access. Examples: 'aspirin', 'caffeine', 'ibuprofen', '2-propanol'.",
            "inputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "name": { "type": "string", "description": "Chemical name (IUPAC, common, or trade name)", "minLength": 1, "maxLength": 500 }
                },
                "required": ["name"]
            },
            "outputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "name": { "type": "string" },
                    "smiles": { "type": "string" },
                    "source": { "type": "string" }
                },
                "required": ["name", "smiles", "source"]
            }
        },
        {
            "name": "retrosynthesis",
            "description": "One-step retrosynthetic disconnection via BRICS (Breaking of Retrosynthetically Interesting Chemical Substructures, Dien 2008). Identifies all BRICS-breakable bonds, cuts each one individually, and returns the resulting fragment pairs ranked by their maximum SA Score (1=easy to synthesize, 10=hard). Lower max-SA means both building blocks are easier to make. Useful for identifying practical synthetic disconnections for drug-like molecules.",
            "inputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "smiles": { "type": "string", "description": "SMILES string of the target molecule", "minLength": 1, "maxLength": 100000 }
                },
                "required": ["smiles"]
            },
            "outputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "target": { "type": "string" },
                    "target_sa_score": { "type": "number" },
                    "disconnections": {
                        "type": "array",
                        "items": {
                            "type": "object",
                            "additionalProperties": false,
                            "properties": {
                                "fragments": { "type": "array", "items": { "type": "string" }, "minItems": 2, "maxItems": 2 },
                                "fragment_sa_scores": { "type": "array", "items": { "type": "number" }, "minItems": 2, "maxItems": 2 },
                                "max_fragment_sa": { "type": "number" }
                            },
                            "required": ["fragments", "fragment_sa_scores", "max_fragment_sa"]
                        }
                    },
                    "total_brics_bonds": { "type": "integer", "minimum": 0 },
                    "note": { "type": "string" }
                },
                "required": ["target", "target_sa_score", "disconnections", "total_brics_bonds", "note"]
            }
        },
        {
            "name": "smiles_to_moljson",
            "description": "Convert a SMILES string to MolJSON — a JSON-based molecular representation designed for LLM compatibility. MolJSON makes atoms, bonds, and connectivity explicit without domain-specific parsing rules.",
            "inputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "smiles": { "type": "string", "description": "SMILES string to convert", "minLength": 1, "maxLength": 100000 }
                },
                "required": ["smiles"]
            },
            "outputSchema": {
                "type": "string",
                "description": "MolJSON, serialized as a JSON string (the tool's own output IS a JSON document, so the MCP structuredContent value is a string containing it, not an already-parsed object)."
            }
        },
        {
            "name": "moljson_to_smiles",
            "description": "Convert a MolJSON string to canonical SMILES.",
            "inputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "json": { "type": "string", "description": "MolJSON string to convert", "minLength": 1, "maxLength": 1000000 }
                },
                "required": ["json"]
            },
            "outputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "canonical_smiles": { "type": "string" }
                },
                "required": ["canonical_smiles"]
            }
        },
        {
            "name": "representation_router",
            "description": "Convert SMILES to the best molecular text representation for an LLM task. Based on arXiv 2026: CML/MolJSON outperform SMILES on structural reasoning; InChI is best for identification.",
            "inputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "smiles": { "type": "string", "description": "SMILES string", "minLength": 1, "maxLength": 100000 },
                    "task":   { "type": "string",
                                "enum": ["structural_reasoning","shortest_path","identification","property_prediction","generation","editing","default"],
                                "description": "LLM task type — omit to use default (canonical_smiles)" }
                },
                "required": ["smiles"]
            },
            "outputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "task": { "type": "string" },
                    "format": { "type": "string", "enum": ["moljson", "inchi", "cml", "canonical_smiles"] },
                    "representation": { "type": "string" }
                },
                "required": ["task", "format", "representation"]
            }
        },
        {
            "name": "molecule_context_pack",
            "description": "Assemble a rich molecular context for LLM/RAG use. Returns identifiers, physicochemical properties, drug-likeness flags, ADMET profile, structural alerts, and MolJSON representation in a single JSON object.",
            "inputSchema": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "smiles": { "type": "string", "description": "SMILES string", "minLength": 1, "maxLength": 100000 },
                    "format": { "type": "string",
                                "enum": ["json", "markdown", "prompt"],
                                "description": "Output format: json (default), markdown (for LLM prompts), prompt (compact one-liner)" }
                },
                "required": ["smiles"]
            },
            "outputSchema": {
                "description": "Shape depends on the `format` argument: `markdown`/`prompt` return a short text summary; the default `json` format returns the full structured context object.",
                "oneOf": [
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "properties": {
                            "format": { "const": "markdown" },
                            "context": { "type": "string" }
                        },
                        "required": ["format", "context"]
                    },
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "properties": {
                            "format": { "const": "prompt" },
                            "context": { "type": "string" }
                        },
                        "required": ["format", "context"]
                    },
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "properties": {
                            "identifiers": {
                                "type": "object",
                                "additionalProperties": false,
                                "properties": {
                                    "smiles": { "type": "string" },
                                    "inchi": { "type": "string" },
                                    "mw": { "type": "number" },
                                    "exact_mass": { "type": "number" }
                                },
                                "required": ["smiles", "inchi", "mw", "exact_mass"]
                            },
                            "properties": {
                                "type": "object",
                                "additionalProperties": false,
                                "properties": {
                                    "logp": { "type": "number" },
                                    "tpsa": { "type": "number" },
                                    "hbd": { "type": "integer" },
                                    "hba": { "type": "integer" },
                                    "rotatable_bonds": { "type": "integer" },
                                    "heavy_atoms": { "type": "integer" },
                                    "qed": { "type": "number" },
                                    "sa_score": { "type": "number" }
                                },
                                "required": ["logp", "tpsa", "hbd", "hba", "rotatable_bonds", "heavy_atoms", "qed", "sa_score"]
                            },
                            "drug_likeness": {
                                "type": "object",
                                "additionalProperties": false,
                                "properties": {
                                    "lipinski_passes": { "type": "boolean" },
                                    "pains_passes": { "type": "boolean" },
                                    "brenk_passes": { "type": "boolean" },
                                    "pains_alerts": { "type": "array", "items": { "type": "string" } },
                                    "brenk_alerts": { "type": "array", "items": { "type": "string" } }
                                },
                                "required": ["lipinski_passes", "pains_passes", "brenk_passes", "pains_alerts", "brenk_alerts"]
                            },
                            "representations": {
                                "type": "object",
                                "additionalProperties": false,
                                "properties": {
                                    "canonical_smiles": { "type": "string" },
                                    "moljson": { "type": "string" }
                                },
                                "required": ["canonical_smiles", "moljson"]
                            }
                        },
                        "required": ["identifiers", "properties", "drug_likeness", "representations"]
                    }
                ]
            }
        }
    ]})
}

// ── tool dispatch ─────────────────────────────────────────────────────────────

/// Dispatch a `tools/call` by name. Returns the bare JSON payload on
/// success — callers (the presentation layer) decide how to wrap it per
/// protocol era.
pub fn call_tool(name: &str, args: &Value) -> Result<Value, ToolCallError> {
    match name {
        "parse_smiles" => tool_parse_smiles(args),
        "calc_properties" => tool_calc_properties(args),
        "ecfp4" => tool_ecfp4(args),
        "tanimoto" => tool_tanimoto(args),
        "smarts_match" => tool_smarts_match(args),
        "canonical_smiles" => tool_canonical_smiles(args),
        "find_mcs" => tool_find_mcs(args),
        "generate_3d" => tool_generate_3d(args),
        "pains_check" => tool_pains_check(args),
        "brenk_check" => tool_brenk_check(args),
        "sa_score" => tool_sa_score(args),
        "admet_profile" => tool_admet_profile(args),
        "boiled_egg" => tool_boiled_egg(args),
        "lipinski_check" => tool_lipinski_check(args),
        "name_to_smiles" => tool_name_to_smiles(args),
        "retrosynthesis" => tool_retrosynthesis(args),
        "smiles_to_moljson" => tool_smiles_to_moljson(args),
        "moljson_to_smiles" => tool_moljson_to_smiles(args),
        "representation_router" => tool_representation_router(args),
        "molecule_context_pack" => tool_molecule_context_pack(args),
        _ => Err(ToolCallError::invalid_args(format!("Unknown tool: {name}"))),
    }
}

// ── individual tools ──────────────────────────────────────────────────────────

fn tool_parse_smiles(args: &Value) -> Result<Value, ToolCallError> {
    let smiles = get_str(args, "smiles")?;
    let mol = parse_mol_arg(smiles)?;
    Ok(json!({
        "valid": true,
        "atoms": mol.atom_count(),
        "bonds": mol.bond_count(),
        "mol_weight": round3(molecular_weight(&mol)),
        "smiles": chematic_smiles::write(&mol)
    }))
}

fn tool_calc_properties(args: &Value) -> Result<Value, ToolCallError> {
    let smiles = get_str(args, "smiles")?;
    let mol = parse_mol_arg(smiles)?;
    Ok(json!({
        "mw":               round3(molecular_weight(&mol)),
        "exact_mass":       round2(exact_mass(&mol) * 100.0) / 100.0,
        "logp":             round3(logp_crippen(&mol)),
        "tpsa":             round2(tpsa(&mol)),
        "hbd":              hbd_count(&mol),
        "hba":              hba_count(&mol),
        "rotatable_bonds":  rotatable_bond_count(&mol),
        "heavy_atom_count": heavy_atom_count(&mol),
        "qed":              round3(qed(&mol))
    }))
}

fn tool_ecfp4(args: &Value) -> Result<Value, ToolCallError> {
    let smiles = get_str(args, "smiles")?;
    let mol = parse_mol_arg(smiles)?;
    let fp = ecfp4(&mol);
    Ok(json!({
        "fingerprint": bitvec_to_hex(&fp),
        "popcount": fp.popcount()
    }))
}

fn tool_tanimoto(args: &Value) -> Result<Value, ToolCallError> {
    let smiles1 = get_str(args, "smiles1")?;
    let smiles2 = get_str(args, "smiles2")?;
    let mol1 = parse_mol_arg(smiles1)?;
    let mol2 = parse_mol_arg(smiles2)?;
    let sim = tanimoto_ecfp4(&mol1, &mol2);
    Ok(json!({
        "similarity": round3(sim),
        "similarity_percent": round2(sim * 100.0)
    }))
}

fn tool_smarts_match(args: &Value) -> Result<Value, ToolCallError> {
    let smarts = get_str(args, "smarts")?;
    let smiles = get_str(args, "smiles")?;
    let query = parse_smarts(smarts).map_err(|e| {
        ToolCallError::domain("INVALID_SMARTS", format!("Invalid SMARTS '{smarts}': {e}"))
    })?;
    let mol = parse_mol_arg(smiles)?;
    let matches = find_matches(&query, &mol);
    let atom_maps: Vec<Vec<u32>> = matches
        .iter()
        .map(|m| {
            let mut atoms: Vec<(usize, u32)> = m.iter().map(|(&q, &a)| (q, a.0)).collect();
            atoms.sort_by_key(|(q, _)| *q);
            atoms.into_iter().map(|(_, a)| a).collect()
        })
        .collect();
    Ok(json!({
        "matches": !matches.is_empty(),
        "match_count": matches.len(),
        "atom_maps": atom_maps
    }))
}

fn tool_canonical_smiles(args: &Value) -> Result<Value, ToolCallError> {
    let smiles = get_str(args, "smiles")?;
    let mol = parse_mol_arg(smiles)?;
    Ok(json!({ "canonical": chematic_smiles::canonical_smiles(&mol) }))
}

fn tool_find_mcs(args: &Value) -> Result<Value, ToolCallError> {
    let smiles_list = args
        .get("smiles_list")
        .and_then(|v| v.as_array())
        .ok_or_else(|| ToolCallError::invalid_args("Missing or invalid smiles_list argument"))?;
    if smiles_list.len() < 2 {
        return Err(ToolCallError::invalid_args(
            "find_mcs requires at least 2 molecules",
        ));
    }
    if smiles_list.len() > 20 {
        return Err(ToolCallError::invalid_args(
            "find_mcs accepts at most 20 molecules",
        ));
    }
    let mols: Result<Vec<_>, ToolCallError> = smiles_list
        .iter()
        .map(|v| {
            v.as_str()
                .ok_or_else(|| ToolCallError::invalid_args("smiles_list must contain strings"))
                .and_then(parse_mol_arg)
        })
        .collect();
    let mols = mols?;
    for mol in &mols {
        if mol.atom_count() > 200 {
            return Err(ToolCallError::domain(
                "MOLECULE_TOO_LARGE",
                "find_mcs: molecule exceeds 200-atom limit",
            ));
        }
    }
    let mol_refs: Vec<&chematic_core::Molecule> = mols.iter().collect();
    let config = McsConfig {
        timeout_ms: Some(5_000),
        ..McsConfig::default()
    };
    let qmol = find_mcs_with_config(&mol_refs, &config);
    if qmol.atoms.is_empty() {
        return Ok(json!({ "mcs": null, "atom_count": 0, "bond_count": 0 }));
    }
    let mol = qmol_to_molecule(&qmol);
    Ok(json!({
        "mcs": chematic_smiles::canonical_smiles(&mol),
        "atom_count": qmol.atoms.len(),
        "bond_count": qmol.bonds.len()
    }))
}

fn tool_generate_3d(args: &Value) -> Result<Value, ToolCallError> {
    let smiles = get_str(args, "smiles")?;
    let mol = parse_mol_arg(smiles)?;
    let coords = generate_and_minimize_dreiding(&mol);
    let xyz = write_xyz(&mol, &coords, smiles);
    Ok(json!({
        "xyz": xyz,
        "atom_count": mol.atom_count()
    }))
}

fn tool_pains_check(args: &Value) -> Result<Value, ToolCallError> {
    let smiles = get_str(args, "smiles")?;
    let mol = parse_mol_arg(smiles)?;
    let passes = pains_passes(&mol);
    let alerts: Vec<&str> = pains_matches(&mol);
    Ok(json!({
        "passes": passes,
        "alert_count": alerts.len(),
        "alerts": alerts
    }))
}

fn tool_brenk_check(args: &Value) -> Result<Value, ToolCallError> {
    let smiles = get_str(args, "smiles")?;
    let mol = parse_mol_arg(smiles)?;
    let passes = brenk_passes(&mol);
    let alerts: Vec<&str> = brenk_matches(&mol);
    Ok(json!({
        "passes": passes,
        "alert_count": alerts.len(),
        "alerts": alerts
    }))
}

fn tool_sa_score(args: &Value) -> Result<Value, ToolCallError> {
    let smiles = get_str(args, "smiles")?;
    let mol = parse_mol_arg(smiles)?;
    let score = sa_score(&mol);
    Ok(json!({
        "sa_score": round3(score),
        "easy_to_synthesize": score < 6.0,
        "note": "1 = easiest, 10 = hardest; < 6 = synthesizable, > 6 = challenging"
    }))
}

fn tool_admet_profile(args: &Value) -> Result<Value, ToolCallError> {
    let smiles = get_str(args, "smiles")?;
    let mol = parse_mol_arg(smiles)?;
    let p = admet_profile(&mol);
    Ok(json!({
        "bbb_score": round3(p.bbb_score),
        "bbb_passes": p.bbb_passes,
        "caco2_logp": round3(p.caco2),
        "herg_risk": round3(p.herg_risk),
        "cyp3a4_risk": round3(p.cyp3a4_risk),
        "pka_acid": p.pka_acid.map(round2),
        "pka_base": p.pka_base.map(round2),
        "esol_logs": round3(p.esol),
        "logd_74": round3(p.logd74),
        "mw": round2(p.mw),
        "logp": round3(p.logp),
        "tpsa": round2(p.tpsa),
        "hbd": p.hbd,
        "hba": p.hba,
        "rotatable_bonds": p.rotatable_bonds,
        "ames_risk": round3(p.ames_risk),
        "ppb_percent": round2(p.ppb),
        "clearance_class": format!("{:?}", p.clearance)
    }))
}

fn tool_boiled_egg(args: &Value) -> Result<Value, ToolCallError> {
    let smiles = get_str(args, "smiles")?;
    let mol = parse_mol_arg(smiles)?;
    let e = boiled_egg(&mol);
    Ok(json!({
        "gi_absorbed": e.gi_absorbed,
        "bbb_penetrant": e.bbb_penetrant,
        "logp": round3(e.logp),
        "tpsa": round2(e.tpsa),
        "method": "BOILED-Egg (Daina & Zoete 2016)",
        "thresholds": {
            "gi_white": "logP ≤ 5.88 AND TPSA ≤ 131.6 Å²",
            "bbb_yolk": "logP ∈ [-0.3, 6.1] AND TPSA ≤ 71.1 Å²"
        }
    }))
}

fn tool_lipinski_check(args: &Value) -> Result<Value, ToolCallError> {
    let smiles = get_str(args, "smiles")?;
    let mol = parse_mol_arg(smiles)?;
    let passes = lipinski_passes(&mol);
    let mw = round2(molecular_weight(&mol));
    let logp = round3(logp_crippen(&mol));
    let hbd = hbd_count(&mol);
    let hba = hba_count(&mol);
    Ok(json!({
        "passes": passes,
        "mw": mw,
        "logp": logp,
        "hbd": hbd,
        "hba": hba,
        "rules": {
            "mw_le_500": mw <= 500.0,
            "logp_le_5": logp <= 5.0,
            "hbd_le_5": hbd <= 5,
            "hba_le_10": hba <= 10
        }
    }))
}

fn tool_name_to_smiles(args: &Value) -> Result<Value, ToolCallError> {
    let name = get_str(args, "name")?;
    if name.len() > 500 {
        return Err(ToolCallError::invalid_args(
            "compound name too long (max 500 characters)",
        ));
    }
    // Percent-encode the name for the URL path segment.
    // Iterate over UTF-8 bytes of each char so that multi-byte code points
    // (e.g. accented letters, CJK) are encoded correctly as %XX%XX sequences
    // rather than truncated to their low byte via `c as u8`.
    let mut encoded = String::with_capacity(name.len() * 3);
    let mut buf = [0u8; 4];
    for c in name.chars() {
        if c.is_ascii_alphanumeric() || "-_.~".contains(c) {
            encoded.push(c);
        } else if c == ' ' {
            encoded.push_str("%20");
        } else {
            for &byte in c.encode_utf8(&mut buf).as_bytes() {
                encoded.push('%');
                encoded.push(char::from_digit((byte >> 4) as u32, 16).unwrap_or('0'));
                encoded.push(char::from_digit((byte & 0xf) as u32, 16).unwrap_or('0'));
            }
        }
    }
    let encoded = encoded;

    let url = format!(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/IsomericSMILES/JSON",
        encoded
    );

    let agent = ureq::config::Config::builder()
        .timeout_global(Some(std::time::Duration::from_secs(10)))
        .build()
        .new_agent();
    let mut resp = agent.get(&url).call().map_err(|e| {
        ToolCallError::domain(
            "PUBCHEM_LOOKUP_FAILED",
            format!("PubChem request failed: {e}"),
        )
    })?;

    let raw = resp.body_mut().read_to_string().map_err(|e| {
        ToolCallError::domain(
            "PUBCHEM_LOOKUP_FAILED",
            format!("PubChem response read error: {e}"),
        )
    })?;
    let body: Value = serde_json::from_str(&raw).map_err(|e| {
        ToolCallError::domain(
            "PUBCHEM_LOOKUP_FAILED",
            format!("PubChem response parse error: {e}"),
        )
    })?;

    let smiles = body
        .pointer("/PropertyTable/Properties/0/IsomericSMILES")
        .and_then(|v| v.as_str())
        .ok_or_else(|| {
            ToolCallError::domain(
                "PUBCHEM_LOOKUP_FAILED",
                format!("Name not found in PubChem: {name}"),
            )
        })?;

    Ok(json!({
        "name": name,
        "smiles": smiles,
        "source": "PubChem"
    }))
}

fn tool_retrosynthesis(args: &Value) -> Result<Value, ToolCallError> {
    let smiles = get_str(args, "smiles")?;
    let mol = parse_mol_arg(smiles)?;

    // Guard against DoS from very large molecules: BRICS runs find_sssr (O(V³))
    // and per-bond BFS (O(V)) for each breakable bond.
    if mol.atom_count() > 500 {
        return Err(ToolCallError::domain(
            "MOLECULE_TOO_LARGE",
            format!(
                "molecule too large for retrosynthesis ({} atoms; maximum is 500)",
                mol.atom_count()
            ),
        ));
    }

    if component_count(&mol) > 1 {
        return Err(ToolCallError::domain(
            "DISCONNECTED_MOLECULE",
            "retrosynthesis requires a single connected molecule; \
             input appears to be a mixture or salt",
        ));
    }

    let target_sa = round3(sa_score(&mol));
    let target_canon = chematic_smiles::canonical_smiles(&mol);
    let bonds = brics_bonds(&mol);

    if bonds.is_empty() {
        return Ok(json!({
            "target": target_canon,
            "target_sa_score": target_sa,
            "disconnections": [],
            "total_brics_bonds": 0,
            "note": "No BRICS-breakable bonds found. Molecule may already be a simple building block."
        }));
    }

    let mut disconnections: Vec<Value> = Vec::new();

    for (a_idx, b_idx) in &bonds {
        let atoms_a = atoms_reachable_excl_bond(&mol, *a_idx, *a_idx, *b_idx);
        let atoms_b = atoms_reachable_excl_bond(&mol, *b_idx, *a_idx, *b_idx);

        if atoms_a.len() + atoms_b.len() != mol.atom_count() {
            continue; // defensive: skip if the split doesn't partition cleanly
        }

        let frag_a = build_submol(&mol, &atoms_a);
        let frag_b = build_submol(&mol, &atoms_b);

        let smiles_a = chematic_smiles::canonical_smiles(&frag_a);
        let smiles_b = chematic_smiles::canonical_smiles(&frag_b);
        let sa_a = round3(sa_score(&frag_a));
        let sa_b = round3(sa_score(&frag_b));
        let max_sa = if sa_a > sa_b { sa_a } else { sa_b };

        disconnections.push(json!({
            "fragments": [smiles_a, smiles_b],
            "fragment_sa_scores": [sa_a, sa_b],
            "max_fragment_sa": round3(max_sa)
        }));
    }

    // Rank by max SA score ascending (easiest disconnections first).
    disconnections.sort_by(|a, b| {
        a["max_fragment_sa"]
            .as_f64()
            .unwrap_or(10.0)
            .partial_cmp(&b["max_fragment_sa"].as_f64().unwrap_or(10.0))
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    Ok(json!({
        "target": target_canon,
        "target_sa_score": target_sa,
        "disconnections": disconnections,
        "total_brics_bonds": bonds.len(),
        "note": "Disconnections ranked by max SA score of fragments (1=easy, 10=hard). Lower = both building blocks easier to synthesize."
    }))
}

fn tool_smiles_to_moljson(args: &Value) -> Result<Value, ToolCallError> {
    let smiles = get_str(args, "smiles")?;
    let mol = parse_mol_arg(smiles)?;
    Ok(Value::String(write_moljson(&mol)))
}

fn tool_moljson_to_smiles(args: &Value) -> Result<Value, ToolCallError> {
    let json_str = get_str(args, "json")?;
    let mol = parse_moljson(json_str)
        .map_err(|e| ToolCallError::domain("INVALID_MOLJSON", e.to_string()))?;
    Ok(json!({ "canonical_smiles": chematic_smiles::canonical_smiles(&mol) }))
}

fn tool_representation_router(args: &Value) -> Result<Value, ToolCallError> {
    let smiles = get_str(args, "smiles")?;
    let task = args
        .get("task")
        .and_then(|v| v.as_str())
        .unwrap_or("default");
    let mol = parse_mol_arg(smiles)?;
    let (format, repr) = match task {
        "structural_reasoning" | "shortest_path" | "graph_reasoning" => {
            ("moljson", write_moljson(&mol))
        }
        "identification" | "exact_match" => ("inchi", inchi(&mol)),
        "editing" => ("cml", write_cml(&mol, None)),
        _ => ("canonical_smiles", chematic_smiles::canonical_smiles(&mol)),
    };
    Ok(json!({
        "task": task,
        "format": format,
        "representation": repr
    }))
}

fn tool_molecule_context_pack(args: &Value) -> Result<Value, ToolCallError> {
    use chematic_chem::{
        brenk_matches, brenk_passes, exact_mass, hba_count, hbd_count, heavy_atom_count,
        lipinski_passes, pains_matches, pains_passes, qed, rotatable_bond_count, sa_score,
    };

    let smiles = get_str(args, "smiles")?;
    let format = args
        .get("format")
        .and_then(|v| v.as_str())
        .unwrap_or("json");
    let mol = parse_mol_arg(smiles)?;
    let mw = round2(molecular_weight(&mol));
    let logp = round2(logp_crippen(&mol));
    let tp = round2(tpsa(&mol));
    let hbd = hbd_count(&mol);
    let hba = hba_count(&mol);
    let rb = rotatable_bond_count(&mol);
    let hac = heavy_atom_count(&mol);
    let q = round3(qed(&mol));
    let sa = round2(sa_score(&mol));
    let em = round3(exact_mass(&mol) as f64);
    let lip = lipinski_passes(&mol);
    let pains_ok = pains_passes(&mol);
    let brenk_ok = brenk_passes(&mol);
    let pains_list = pains_matches(&mol);
    let brenk_list = brenk_matches(&mol);
    let inchi_str = inchi(&mol);
    let moljson_str = write_moljson(&mol);
    let canonical = chematic_smiles::canonical_smiles(&mol);

    match format {
        "markdown" => {
            let alert_str = if pains_list.is_empty() && brenk_list.is_empty() {
                "none".to_owned()
            } else {
                let mut parts = Vec::new();
                if !pains_list.is_empty() {
                    parts.push(format!("PAINS: {}", pains_list.join(", ")));
                }
                if !brenk_list.is_empty() {
                    parts.push(format!("Brenk: {}", brenk_list.join(", ")));
                }
                parts.join("; ")
            };
            let admet = admet_profile(&mol);
            let md = format!(
                "## Molecule\n- **SMILES**: {canonical}\n- **MW**: {mw} Da  ExactMass: {em}\n\
                 ## Properties\n- LogP: {logp}, TPSA: {tp} Å², HBD: {hbd}, HBA: {hba}\n\
                 - RotBonds: {rb}, HAC: {hac}, QED: {q}, SA: {sa}\n\
                 ## Drug-Likeness\n- Lipinski: {lip}, PAINS: {pains_ok}, Brenk: {brenk_ok}\n\
                 - Alerts: {alert_str}\n\
                 ## ADMET\n- BBB passes: {bbb}, Caco-2: {caco2:.2}, CYP3A4 risk: {cyp:.2}\n\
                 ## MolJSON\n```json\n{moljson_str}\n```",
                bbb = admet.bbb_passes,
                caco2 = admet.caco2,
                cyp = admet.cyp3a4_risk,
            );
            Ok(json!({ "format": "markdown", "context": md }))
        }
        "prompt" => {
            let flags = if lip { "Lipinski✓" } else { "Lipinski✗" };
            let prompt = format!(
                "{canonical} | MW={mw} | LogP={logp} | TPSA={tp} | HBD={hbd} HBA={hba} | QED={q} | {flags}"
            );
            Ok(json!({ "format": "prompt", "context": prompt }))
        }
        _ => {
            // json (default)
            Ok(json!({
                "identifiers": {
                    "smiles": canonical,
                    "inchi": inchi_str,
                    "mw": mw,
                    "exact_mass": em,
                },
                "properties": {
                    "logp": logp, "tpsa": tp,
                    "hbd": hbd, "hba": hba,
                    "rotatable_bonds": rb, "heavy_atoms": hac,
                    "qed": q, "sa_score": sa,
                },
                "drug_likeness": {
                    "lipinski_passes": lip,
                    "pains_passes": pains_ok,
                    "brenk_passes": brenk_ok,
                    "pains_alerts": pains_list,
                    "brenk_alerts": brenk_list,
                },
                "representations": {
                    "canonical_smiles": canonical,
                    "moljson": moljson_str,
                }
            }))
        }
    }
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn args(pairs: &[(&str, &str)]) -> Value {
        let mut obj = serde_json::Map::new();
        for (k, v) in pairs {
            obj.insert(k.to_string(), json!(*v));
        }
        Value::Object(obj)
    }

    #[test]
    fn test_parse_benzene() {
        let v = tool_parse_smiles(&args(&[("smiles", "c1ccccc1")])).unwrap();
        assert_eq!(v["valid"], true);
        assert_eq!(v["atoms"], 6);
        assert_eq!(v["bonds"], 6);
    }

    #[test]
    fn test_calc_properties_benzene() {
        let v = tool_calc_properties(&args(&[("smiles", "c1ccccc1")])).unwrap();
        assert!(v["mw"].as_f64().unwrap() > 78.0);
        assert_eq!(v["hbd"], 0);
        assert_eq!(v["hba"], 0);
    }

    #[test]
    fn test_ecfp4_benzene() {
        let v = tool_ecfp4(&args(&[("smiles", "c1ccccc1")])).unwrap();
        let hex = v["fingerprint"].as_str().unwrap();
        assert_eq!(hex.len(), 512); // 256 bytes = 512 hex chars
        assert!(v["popcount"].as_u64().unwrap() > 0);
    }

    #[test]
    fn test_tanimoto_self_similarity() {
        let a = args(&[("smiles1", "c1ccccc1"), ("smiles2", "c1ccccc1")]);
        let v = tool_tanimoto(&a).unwrap();
        let sim = v["similarity"].as_f64().unwrap();
        assert!(
            (sim - 1.0).abs() < 1e-9,
            "self-similarity must be 1.0, got {sim}"
        );
    }

    #[test]
    fn test_tanimoto_different_molecules() {
        let a = args(&[("smiles1", "c1ccccc1"), ("smiles2", "CCO")]);
        let v = tool_tanimoto(&a).unwrap();
        let sim = v["similarity"].as_f64().unwrap();
        assert!(sim < 1.0);
    }

    #[test]
    fn test_smarts_match_hit() {
        let a = args(&[("smarts", "c1ccccc1"), ("smiles", "c1ccccc1")]);
        let v = tool_smarts_match(&a).unwrap();
        assert_eq!(v["matches"], true);
        assert!(v["match_count"].as_u64().unwrap() >= 1);
    }

    #[test]
    fn test_smarts_match_miss() {
        let a = args(&[("smarts", "N"), ("smiles", "c1ccccc1")]);
        let v = tool_smarts_match(&a).unwrap();
        assert_eq!(v["matches"], false);
    }

    #[test]
    fn test_canonical_smiles() {
        let v = tool_canonical_smiles(&args(&[("smiles", "C1=CC=CC=C1")])).unwrap();
        let canon = v["canonical"].as_str().unwrap();
        assert!(!canon.is_empty());
    }

    #[test]
    fn test_find_mcs_two_molecules() {
        let smiles_list = json!(["c1ccccc1", "c1ccccc1O"]);
        let mut args_obj = serde_json::Map::new();
        args_obj.insert("smiles_list".to_string(), smiles_list);
        let v = tool_find_mcs(&Value::Object(args_obj)).unwrap();
        assert!(v["atom_count"].as_u64().unwrap() >= 6);
    }

    #[test]
    fn test_find_mcs_requires_two_mols() {
        let smiles_list = json!(["c1ccccc1"]);
        let mut args_obj = serde_json::Map::new();
        args_obj.insert("smiles_list".to_string(), smiles_list);
        let result = tool_find_mcs(&Value::Object(args_obj));
        assert!(result.is_err());
        assert!(result.unwrap_err().is_invalid_args());
    }

    #[test]
    fn test_generate_3d_benzene() {
        let v = tool_generate_3d(&args(&[("smiles", "c1ccccc1")])).unwrap();
        assert_eq!(v["atom_count"], 6);
        let xyz = v["xyz"].as_str().unwrap();
        assert!(xyz.contains('C'));
    }

    #[test]
    fn test_parse_invalid_smiles_is_domain_error() {
        // Unbalanced ring closure — definitely invalid
        let result = tool_parse_smiles(&args(&[("smiles", "C1CC")]));
        let err = result.unwrap_err();
        assert!(
            !err.is_invalid_args(),
            "bad SMILES is a domain error, not an argument-shape error"
        );
    }

    #[test]
    fn test_missing_argument_is_invalid_args() {
        let result = tool_parse_smiles(&Value::Object(serde_json::Map::new()));
        let err = result.unwrap_err();
        assert!(err.is_invalid_args());
    }

    #[test]
    fn test_list_tools_count() {
        let tools = list_tools();
        let count = tools["tools"].as_array().unwrap().len();
        assert_eq!(count, TOOL_COUNT);
    }

    #[test]
    fn test_retrosynthesis_aspirin() {
        // Aspirin has 2 BRICS-breakable bonds (ester C-O, aryl C-O)
        let v = tool_retrosynthesis(&args(&[("smiles", "CC(=O)Oc1ccccc1C(=O)O")])).unwrap();
        assert!(
            v["total_brics_bonds"].as_u64().unwrap() >= 1,
            "aspirin should have ≥1 BRICS bond"
        );
        let discos = v["disconnections"].as_array().unwrap();
        assert!(!discos.is_empty(), "should have at least one disconnection");
        // First disconnection should have 2 fragments
        assert_eq!(
            discos[0]["fragments"].as_array().unwrap().len(),
            2,
            "each disconnection yields exactly 2 fragments"
        );
        // max_fragment_sa should be a reasonable number
        let max_sa = discos[0]["max_fragment_sa"].as_f64().unwrap();
        assert!(
            (1.0..=10.0).contains(&max_sa),
            "SA score out of range: {max_sa}"
        );
        // Disconnections should be sorted ascending by max_fragment_sa
        if discos.len() >= 2 {
            let sa0 = discos[0]["max_fragment_sa"].as_f64().unwrap();
            let sa1 = discos[1]["max_fragment_sa"].as_f64().unwrap();
            assert!(sa0 <= sa1, "disconnections not sorted: {sa0} > {sa1}");
        }
    }

    #[test]
    fn test_retrosynthesis_benzene_no_bonds() {
        // Benzene has no BRICS-breakable bonds
        let v = tool_retrosynthesis(&args(&[("smiles", "c1ccccc1")])).unwrap();
        assert_eq!(v["total_brics_bonds"], 0);
        assert_eq!(v["disconnections"].as_array().unwrap().len(), 0);
    }

    #[test]
    fn test_retrosynthesis_disconnected_mol_error() {
        // Mixture (salt): should return a domain error, not an argument-shape error
        let result = tool_retrosynthesis(&args(&[("smiles", "CC.OO")]));
        assert!(result.is_err(), "disconnected molecule should return error");
        assert!(!result.unwrap_err().is_invalid_args());
    }

    #[test]
    fn test_pains_check_clean() {
        let v = tool_pains_check(&args(&[("smiles", "CCO")])).unwrap();
        assert_eq!(v["passes"], true);
        assert_eq!(v["alert_count"], 0);
    }

    #[test]
    fn test_brenk_check_clean() {
        let v = tool_brenk_check(&args(&[("smiles", "CCO")])).unwrap();
        assert_eq!(v["passes"], true);
    }

    #[test]
    fn test_sa_score_ethanol() {
        let v = tool_sa_score(&args(&[("smiles", "CCO")])).unwrap();
        let score = v["sa_score"].as_f64().unwrap();
        assert!(
            (1.0..=10.0).contains(&score),
            "SA score out of range: {score}"
        );
    }

    #[test]
    fn test_sa_score_aspirin_easy() {
        let v = tool_sa_score(&args(&[("smiles", "CC(=O)Oc1ccccc1C(=O)O")])).unwrap();
        assert_eq!(v["easy_to_synthesize"], true);
    }

    #[test]
    fn test_admet_profile_benzene() {
        let v = tool_admet_profile(&args(&[("smiles", "c1ccccc1")])).unwrap();
        assert!(v.get("bbb_passes").is_some());
        assert!(v.get("clearance_class").is_some());
    }

    #[test]
    fn test_boiled_egg_aspirin() {
        let v = tool_boiled_egg(&args(&[("smiles", "CC(=O)Oc1ccccc1C(=O)O")])).unwrap();
        assert_eq!(v["gi_absorbed"], true);
    }

    #[test]
    fn test_lipinski_check_ethanol() {
        let v = tool_lipinski_check(&args(&[("smiles", "CCO")])).unwrap();
        assert_eq!(v["passes"], true);
        assert!(v["mw"].as_f64().unwrap() < 500.0);
    }

    #[test]
    fn test_smiles_to_moljson_is_bare_string() {
        let v = tool_smiles_to_moljson(&args(&[("smiles", "CCO")])).unwrap();
        assert!(
            v.is_string(),
            "smiles_to_moljson's payload must be a bare JSON string"
        );
    }

    #[test]
    fn test_unknown_tool_is_invalid_args() {
        let result = call_tool("does_not_exist", &json!({}));
        assert!(result.unwrap_err().is_invalid_args());
    }
}
