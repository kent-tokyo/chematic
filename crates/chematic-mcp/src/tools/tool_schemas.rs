//! Static MCP tool registry and JSON Schema documents.

#![forbid(unsafe_code)]

use serde_json::{Value, json};

// ── tool list schema ──────────────────────────────────────────────────────────

/// Number of tools this server exposes. A single source of truth so
/// registry-invariant tests (`exactly 20 unique names`, negative controls
/// for "duplicate a tool name" / "remove a tool") never hardcode `20` twice.
pub const TOOL_COUNT: usize = 20;

pub fn list_tools() -> Value {
    let tools = vec![
        parse_smiles_schema(),
        calc_properties_schema(),
        ecfp4_schema(),
        tanimoto_schema(),
        smarts_match_schema(),
        canonical_smiles_schema(),
        find_mcs_schema(),
        generate_3d_schema(),
        pains_check_schema(),
        brenk_check_schema(),
        sa_score_schema(),
        admet_profile_schema(),
        boiled_egg_schema(),
        lipinski_check_schema(),
        name_to_smiles_schema(),
        retrosynthesis_schema(),
        smiles_to_moljson_schema(),
        moljson_to_smiles_schema(),
        representation_router_schema(),
        molecule_context_pack_schema(),
    ];
    debug_assert_eq!(tools.len(), TOOL_COUNT);
    json!({ "tools": tools })
}

fn parse_smiles_schema() -> Value {
    json!(
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
                }
    )
}

fn calc_properties_schema() -> Value {
    json!(
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
                }
    )
}

fn ecfp4_schema() -> Value {
    json!(
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
                }
    )
}

fn tanimoto_schema() -> Value {
    json!(
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
                }
    )
}

fn smarts_match_schema() -> Value {
    json!(
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
                }
    )
}

fn canonical_smiles_schema() -> Value {
    json!(
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
                }
    )
}

fn find_mcs_schema() -> Value {
    json!(
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
                }
    )
}

fn generate_3d_schema() -> Value {
    json!(
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
                }
    )
}

fn pains_check_schema() -> Value {
    json!(
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
                }
    )
}

fn brenk_check_schema() -> Value {
    json!(
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
                }
    )
}

fn sa_score_schema() -> Value {
    json!(
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
                }
    )
}

fn admet_profile_schema() -> Value {
    json!(
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
                }
    )
}

fn boiled_egg_schema() -> Value {
    json!(
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
                }
    )
}

fn lipinski_check_schema() -> Value {
    json!(
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
                }
    )
}

fn name_to_smiles_schema() -> Value {
    json!(
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
                }
    )
}

fn retrosynthesis_schema() -> Value {
    json!(
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
                }
    )
}

fn smiles_to_moljson_schema() -> Value {
    json!(
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
                }
    )
}

fn moljson_to_smiles_schema() -> Value {
    json!(
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
                }
    )
}

fn representation_router_schema() -> Value {
    json!(
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
                }
    )
}

fn molecule_context_pack_schema() -> Value {
    json!(
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
    )
}
