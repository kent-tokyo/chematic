//! MCP tool implementations for chematic.

#![forbid(unsafe_code)]

use chematic_core::{Atom, AtomIdx, BondOrder, Element, MoleculeBuilder};
use chematic_fp::{ecfp4, tanimoto_ecfp4, BitVec2048};
use chematic_smarts::{AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, find_matches, find_mcs, parse_smarts};
use chematic_3d::{generate_and_minimize_dreiding, write_xyz};
use serde_json::{Value, json};

use chematic_chem::{
    exact_mass, hba_count, hbd_count, heavy_atom_count, logp_crippen, molecular_weight, qed,
    rotatable_bond_count, tpsa,
};

// ── helpers ──────────────────────────────────────────────────────────────────

fn get_str<'a>(args: &'a Value, key: &str) -> Result<&'a str, String> {
    args.get(key)
        .and_then(|v| v.as_str())
        .ok_or_else(|| format!("Missing required argument: {key}"))
}

fn parse_mol(smiles: &str) -> Result<chematic_core::Molecule, String> {
    chematic_smiles::parse(smiles).map_err(|e| format!("Invalid SMILES '{smiles}': {e}"))
}

fn round3(x: f64) -> f64 {
    (x * 1000.0).round() / 1000.0
}

fn round2(x: f64) -> f64 {
    (x * 100.0).round() / 100.0
}

/// Wrap a JSON payload in the MCP content envelope.
fn content(value: &Value) -> Value {
    json!({ "content": [{ "type": "text", "text": value.to_string() }] })
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
                let _ = builder.add_bond(AtomIdx(atom_idx as u32), AtomIdx(*neighbor_idx as u32), order);
            }
        }
    }
    builder.build()
}

// ── tool list schema ──────────────────────────────────────────────────────────

pub fn list_tools() -> Value {
    json!({ "tools": [
        {
            "name": "parse_smiles",
            "description": "Parse a SMILES string and return basic molecular information (atom count, bond count, molecular weight).",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "smiles": { "type": "string", "description": "SMILES string to parse" }
                },
                "required": ["smiles"]
            }
        },
        {
            "name": "calc_properties",
            "description": "Calculate molecular properties: MW, exact mass, LogP (Crippen), TPSA, HBD, HBA, rotatable bonds, heavy atom count, and QED drug-likeness.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "smiles": { "type": "string", "description": "SMILES string" }
                },
                "required": ["smiles"]
            }
        },
        {
            "name": "ecfp4",
            "description": "Compute the ECFP4 (Morgan radius-2) circular fingerprint as a 2048-bit hex string, plus popcount.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "smiles": { "type": "string", "description": "SMILES string" }
                },
                "required": ["smiles"]
            }
        },
        {
            "name": "tanimoto",
            "description": "Compute the Tanimoto (Jaccard) similarity between two molecules using ECFP4 fingerprints.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "smiles1": { "type": "string", "description": "First molecule SMILES" },
                    "smiles2": { "type": "string", "description": "Second molecule SMILES" }
                },
                "required": ["smiles1", "smiles2"]
            }
        },
        {
            "name": "smarts_match",
            "description": "Perform SMARTS substructure search and return whether the pattern matches, match count, and atom index maps.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "smarts": { "type": "string", "description": "SMARTS pattern" },
                    "smiles": { "type": "string", "description": "Molecule SMILES to search in" }
                },
                "required": ["smarts", "smiles"]
            }
        },
        {
            "name": "canonical_smiles",
            "description": "Return the canonical SMILES representation of a molecule.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "smiles": { "type": "string", "description": "Input SMILES string" }
                },
                "required": ["smiles"]
            }
        },
        {
            "name": "find_mcs",
            "description": "Find the maximum common substructure (MCS) across a list of molecules. Returns the MCS as a canonical SMILES string.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "smiles_list": {
                        "type": "array",
                        "items": { "type": "string" },
                        "description": "List of SMILES strings (minimum 2 molecules required)",
                        "minItems": 2
                    }
                },
                "required": ["smiles_list"]
            }
        },
        {
            "name": "generate_3d",
            "description": "Generate 3D coordinates for a molecule using rule-based placement and DREIDING force-field minimization. Returns XYZ format.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "smiles": { "type": "string", "description": "SMILES string" }
                },
                "required": ["smiles"]
            }
        }
    ]})
}

// ── tool dispatch ─────────────────────────────────────────────────────────────

pub fn call_tool(name: &str, args: &Value) -> Result<Value, String> {
    match name {
        "parse_smiles" => tool_parse_smiles(args),
        "calc_properties" => tool_calc_properties(args),
        "ecfp4" => tool_ecfp4(args),
        "tanimoto" => tool_tanimoto(args),
        "smarts_match" => tool_smarts_match(args),
        "canonical_smiles" => tool_canonical_smiles(args),
        "find_mcs" => tool_find_mcs(args),
        "generate_3d" => tool_generate_3d(args),
        _ => Err(format!("Unknown tool: {name}")),
    }
}

// ── individual tools ──────────────────────────────────────────────────────────

fn tool_parse_smiles(args: &Value) -> Result<Value, String> {
    let smiles = get_str(args, "smiles")?;
    let mol = parse_mol(smiles)?;
    Ok(content(&json!({
        "valid": true,
        "atoms": mol.atom_count(),
        "bonds": mol.bond_count(),
        "mol_weight": round3(molecular_weight(&mol)),
        "smiles": chematic_smiles::write(&mol)
    })))
}

fn tool_calc_properties(args: &Value) -> Result<Value, String> {
    let smiles = get_str(args, "smiles")?;
    let mol = parse_mol(smiles)?;
    Ok(content(&json!({
        "mw":               round3(molecular_weight(&mol)),
        "exact_mass":       round2(exact_mass(&mol) * 100.0) / 100.0,
        "logp":             round3(logp_crippen(&mol)),
        "tpsa":             round2(tpsa(&mol)),
        "hbd":              hbd_count(&mol),
        "hba":              hba_count(&mol),
        "rotatable_bonds":  rotatable_bond_count(&mol),
        "heavy_atom_count": heavy_atom_count(&mol),
        "qed":              round3(qed(&mol))
    })))
}

fn tool_ecfp4(args: &Value) -> Result<Value, String> {
    let smiles = get_str(args, "smiles")?;
    let mol = parse_mol(smiles)?;
    let fp = ecfp4(&mol);
    Ok(content(&json!({
        "fingerprint": bitvec_to_hex(&fp),
        "popcount": fp.popcount()
    })))
}

fn tool_tanimoto(args: &Value) -> Result<Value, String> {
    let smiles1 = get_str(args, "smiles1")?;
    let smiles2 = get_str(args, "smiles2")?;
    let mol1 = parse_mol(smiles1)?;
    let mol2 = parse_mol(smiles2)?;
    let sim = tanimoto_ecfp4(&mol1, &mol2);
    Ok(content(&json!({
        "similarity": round3(sim),
        "similarity_percent": round2(sim * 100.0)
    })))
}

fn tool_smarts_match(args: &Value) -> Result<Value, String> {
    let smarts = get_str(args, "smarts")?;
    let smiles = get_str(args, "smiles")?;
    let query = parse_smarts(smarts).map_err(|e| format!("Invalid SMARTS '{smarts}': {e}"))?;
    let mol = parse_mol(smiles)?;
    let matches = find_matches(&query, &mol);
    let atom_maps: Vec<Vec<u32>> = matches
        .iter()
        .map(|m| {
            let mut atoms: Vec<(usize, u32)> =
                m.iter().map(|(&q, &a)| (q, a.0)).collect();
            atoms.sort_by_key(|(q, _)| *q);
            atoms.into_iter().map(|(_, a)| a).collect()
        })
        .collect();
    Ok(content(&json!({
        "matches": !matches.is_empty(),
        "match_count": matches.len(),
        "atom_maps": atom_maps
    })))
}

fn tool_canonical_smiles(args: &Value) -> Result<Value, String> {
    let smiles = get_str(args, "smiles")?;
    let mol = parse_mol(smiles)?;
    Ok(content(&json!({ "canonical": chematic_smiles::canonical_smiles(&mol) })))
}

fn tool_find_mcs(args: &Value) -> Result<Value, String> {
    let smiles_list = args
        .get("smiles_list")
        .and_then(|v| v.as_array())
        .ok_or_else(|| "Missing or invalid smiles_list argument".to_string())?;
    if smiles_list.len() < 2 {
        return Err("find_mcs requires at least 2 molecules".to_string());
    }
    let mols: Result<Vec<_>, _> = smiles_list
        .iter()
        .map(|v| {
            v.as_str()
                .ok_or_else(|| "smiles_list must contain strings".to_string())
                .and_then(parse_mol)
        })
        .collect();
    let mols = mols?;
    let mol_refs: Vec<&chematic_core::Molecule> = mols.iter().collect();
    let qmol = find_mcs(&mol_refs);
    if qmol.atoms.is_empty() {
        return Ok(content(&json!({ "mcs": null, "atom_count": 0, "bond_count": 0 })));
    }
    let mol = qmol_to_molecule(&qmol);
    Ok(content(&json!({
        "mcs": chematic_smiles::canonical_smiles(&mol),
        "atom_count": qmol.atoms.len(),
        "bond_count": qmol.bonds.len()
    })))
}

fn tool_generate_3d(args: &Value) -> Result<Value, String> {
    let smiles = get_str(args, "smiles")?;
    let mol = parse_mol(smiles)?;
    let coords = generate_and_minimize_dreiding(&mol);
    let xyz = write_xyz(&mol, &coords, smiles);
    Ok(content(&json!({
        "xyz": xyz,
        "atom_count": mol.atom_count()
    })))
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
        let result = tool_parse_smiles(&args(&[("smiles", "c1ccccc1")])).unwrap();
        let text = result["content"][0]["text"].as_str().unwrap();
        let v: Value = serde_json::from_str(text).unwrap();
        assert_eq!(v["valid"], true);
        assert_eq!(v["atoms"], 6);
        assert_eq!(v["bonds"], 6);
    }

    #[test]
    fn test_calc_properties_benzene() {
        let result = tool_calc_properties(&args(&[("smiles", "c1ccccc1")])).unwrap();
        let text = result["content"][0]["text"].as_str().unwrap();
        let v: Value = serde_json::from_str(text).unwrap();
        assert!(v["mw"].as_f64().unwrap() > 78.0);
        assert_eq!(v["hbd"], 0);
        assert_eq!(v["hba"], 0);
    }

    #[test]
    fn test_ecfp4_benzene() {
        let result = tool_ecfp4(&args(&[("smiles", "c1ccccc1")])).unwrap();
        let text = result["content"][0]["text"].as_str().unwrap();
        let v: Value = serde_json::from_str(text).unwrap();
        let hex = v["fingerprint"].as_str().unwrap();
        assert_eq!(hex.len(), 512); // 256 bytes = 512 hex chars
        assert!(v["popcount"].as_u64().unwrap() > 0);
    }

    #[test]
    fn test_tanimoto_self_similarity() {
        let a = args(&[("smiles1", "c1ccccc1"), ("smiles2", "c1ccccc1")]);
        let result = tool_tanimoto(&a).unwrap();
        let text = result["content"][0]["text"].as_str().unwrap();
        let v: Value = serde_json::from_str(text).unwrap();
        let sim = v["similarity"].as_f64().unwrap();
        assert!((sim - 1.0).abs() < 1e-9, "self-similarity must be 1.0, got {sim}");
    }

    #[test]
    fn test_tanimoto_different_molecules() {
        let a = args(&[("smiles1", "c1ccccc1"), ("smiles2", "CCO")]);
        let result = tool_tanimoto(&a).unwrap();
        let text = result["content"][0]["text"].as_str().unwrap();
        let v: Value = serde_json::from_str(text).unwrap();
        let sim = v["similarity"].as_f64().unwrap();
        assert!(sim < 1.0);
    }

    #[test]
    fn test_smarts_match_hit() {
        let a = args(&[("smarts", "c1ccccc1"), ("smiles", "c1ccccc1")]);
        let result = tool_smarts_match(&a).unwrap();
        let text = result["content"][0]["text"].as_str().unwrap();
        let v: Value = serde_json::from_str(text).unwrap();
        assert_eq!(v["matches"], true);
        assert!(v["match_count"].as_u64().unwrap() >= 1);
    }

    #[test]
    fn test_smarts_match_miss() {
        let a = args(&[("smarts", "N"), ("smiles", "c1ccccc1")]);
        let result = tool_smarts_match(&a).unwrap();
        let text = result["content"][0]["text"].as_str().unwrap();
        let v: Value = serde_json::from_str(text).unwrap();
        assert_eq!(v["matches"], false);
    }

    #[test]
    fn test_canonical_smiles() {
        let result = tool_canonical_smiles(&args(&[("smiles", "C1=CC=CC=C1")])).unwrap();
        let text = result["content"][0]["text"].as_str().unwrap();
        let v: Value = serde_json::from_str(text).unwrap();
        let canon = v["canonical"].as_str().unwrap();
        assert!(!canon.is_empty());
    }

    #[test]
    fn test_find_mcs_two_molecules() {
        let smiles_list = json!(["c1ccccc1", "c1ccccc1O"]);
        let mut args_obj = serde_json::Map::new();
        args_obj.insert("smiles_list".to_string(), smiles_list);
        let result = tool_find_mcs(&Value::Object(args_obj)).unwrap();
        let text = result["content"][0]["text"].as_str().unwrap();
        let v: Value = serde_json::from_str(text).unwrap();
        assert!(v["atom_count"].as_u64().unwrap() >= 6);
    }

    #[test]
    fn test_find_mcs_requires_two_mols() {
        let smiles_list = json!(["c1ccccc1"]);
        let mut args_obj = serde_json::Map::new();
        args_obj.insert("smiles_list".to_string(), smiles_list);
        let result = tool_find_mcs(&Value::Object(args_obj));
        assert!(result.is_err());
    }

    #[test]
    fn test_generate_3d_benzene() {
        let result = tool_generate_3d(&args(&[("smiles", "c1ccccc1")])).unwrap();
        let text = result["content"][0]["text"].as_str().unwrap();
        let v: Value = serde_json::from_str(text).unwrap();
        assert_eq!(v["atom_count"], 6);
        let xyz = v["xyz"].as_str().unwrap();
        assert!(xyz.contains('C'));
    }

    #[test]
    fn test_parse_invalid_smiles() {
        // Unbalanced ring closure — definitely invalid
        let result = tool_parse_smiles(&args(&[("smiles", "C1CC")]));
        assert!(result.is_err());
    }

    #[test]
    fn test_list_tools_count() {
        let tools = list_tools();
        let count = tools["tools"].as_array().unwrap().len();
        assert_eq!(count, 8);
    }
}
