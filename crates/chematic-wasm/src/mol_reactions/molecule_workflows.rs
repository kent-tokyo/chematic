use super::*;

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

fn parent_transform_json(
    mol: &MolHandle,
    transform: fn(
        &chematic_core::Molecule,
    ) -> (chematic_core::Molecule, chematic_chem::TransformationRecord),
) -> String {
    if mol.inner.atom_count() > WASM_MAX_ATOMS {
        return format!(
            r#"{{"error":"molecule too large (max {} atoms)"}}"#,
            WASM_MAX_ATOMS
        );
    }
    let (parent, _) = transform(&mol.inner);
    let smiles = escape_json_string(&chematic_smiles::canonical_smiles(&parent));
    format!(r#"{{"smiles":"{smiles}","status":"completed"}}"#)
}

/// Compute the fragment Parent and return a status-shaped JSON result.
#[wasm_bindgen]
pub fn fragment_parent_json(mol: &MolHandle) -> String {
    parent_transform_json(mol, chematic_chem::fragment_parent)
}

/// Compute the charge Parent and return a status-shaped JSON result.
#[wasm_bindgen]
pub fn charge_parent_json(mol: &MolHandle) -> String {
    parent_transform_json(mol, chematic_chem::charge_parent)
}

/// Compute the isotope Parent and return a status-shaped JSON result.
#[wasm_bindgen]
pub fn isotope_parent_json(mol: &MolHandle) -> String {
    parent_transform_json(mol, chematic_chem::isotope_parent)
}

/// Compute the stereo Parent and return a status-shaped JSON result.
#[wasm_bindgen]
pub fn stereo_parent_json(mol: &MolHandle) -> String {
    parent_transform_json(mol, chematic_chem::stereo_parent)
}

/// Compute the composed Super Parent with explicit resource limits.
#[wasm_bindgen]
pub fn super_parent_json(
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
    let result = chematic_chem::super_parent(&mol.inner, &limits);
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

/// Compute the composed Super Parent and expose every ordered stage.
#[wasm_bindgen]
pub fn super_parent_report_json(
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
    let (fragment, _) = chematic_chem::fragment_parent(&mol.inner);
    let (charge, _) = chematic_chem::charge_parent(&fragment);
    let (isotope, _) = chematic_chem::isotope_parent(&charge);
    let (stereo, _) = chematic_chem::stereo_parent(&isotope);
    let result = chematic_chem::super_parent(&mol.inner, &limits);
    let status = match result.status {
        chematic_chem::ParentComputationStatus::Completed => "completed",
        chematic_chem::ParentComputationStatus::MaxTransformsReached => "max_transforms_reached",
        chematic_chem::ParentComputationStatus::MaxTautomersReached => "max_tautomers_reached",
        chematic_chem::ParentComputationStatus::TimedOut => "timed_out",
        chematic_chem::ParentComputationStatus::Abstained(_) => "abstained",
        chematic_chem::ParentComputationStatus::InvalidInput(_) => "invalid_input",
        _ => "unknown",
    };
    let stage = |name: &str, molecule: &chematic_core::Molecule| {
        format!(
            r#"{{"name":"{}","smiles":"{}"}}"#,
            name,
            escape_json_string(&chematic_smiles::canonical_smiles(molecule))
        )
    };
    let stages = [
        stage("fragment", &fragment),
        stage("charge", &charge),
        stage("isotope", &isotope),
        stage("stereo", &stereo),
        stage("tautomer", &result.molecule),
    ]
    .join(",");
    let smiles = escape_json_string(&chematic_smiles::canonical_smiles(&result.molecule));
    format!(r#"{{"smiles":"{smiles}","status":"{status}","stages":[{stages}]}}"#)
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
        Err(e) => return json_error(format!("invalid JSON: {e}")),
    };
    let blocked_atoms: std::collections::HashSet<chematic_core::AtomIdx> =
        indices.into_iter().map(chematic_core::AtomIdx).collect();
    let config = chematic_chem::TautomerConfig {
        blocked_atoms,
        ..chematic_chem::TautomerConfig::default()
    };
    let result = chematic_chem::canonical_tautomer_with_config(&mol.inner, &config);
    let smi = chematic_smiles::canonical_smiles(&result);
    json_string(&smi)
}

/// All enumerated tautomers of `mol` as a JSON array of canonical SMILES strings.
///
/// Example return value: `["Oc1cccc2ccccc12","O=C1C=CC=Cc2ccccc21"]`
#[wasm_bindgen]
pub fn enumerate_tautomers_json(mol: &MolHandle) -> String {
    if mol.inner.atom_count() > WASM_MAX_ATOMS {
        return json_error(format!("molecule too large (max {WASM_MAX_ATOMS} atoms)"));
    }
    let tautomers = chematic_chem::enumerate_tautomers(&mol.inner);
    if tautomers.len() > WASM_MAX_BATCH_ITEMS {
        return json_error(format!(
            "tautomer enumeration exceeds maximum result count ({} > {})",
            tautomers.len(),
            WASM_MAX_BATCH_ITEMS
        ));
    }
    let values: Vec<String> = tautomers
        .iter()
        .map(chematic_smiles::canonical_smiles)
        .collect();
    serde_json::to_string(&values).expect("serializing tautomer strings cannot fail")
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
/// `"rdkit_any_heavy_atom"` |
/// `"any"`), `bondCompare` (`"order_or_aromatic"` | `"any"`),
/// `matchChiralTag`, `matchCharge`, `matchIsotope`, `maximizeBonds`, and
/// `maxSearchNodes`.
///
/// Returns a JSON object with `smiles`, `wasTimedOut`, and
/// `wasResourceLimited`. `smiles` is `null` when there is no common
/// substructure; the two status fields distinguish deadline expiry from a
/// configured input/search limit. A returned partial result is not proven
/// optimal when either status is true.
#[wasm_bindgen]
pub fn mcs_smiles_json_with_config(
    smiles_json: &str,
    config_json: &str,
) -> Result<String, JsValue> {
    enforce_wasm_input_len("config_json", config_json)?;
    let mols = parse_mcs_input(smiles_json, "mcs_smiles_json_with_config")?;
    let config: McsConfigJson =
        serde_json::from_str(config_json).map_err(|e| JsValue::from_str(&e.to_string()))?;
    let config = config.into_mcs_config();

    let mol_refs: Vec<&chematic_core::Molecule> = mols.iter().collect();
    let outcome = chematic_smarts::find_mcs_with_config_checked(&mol_refs, &config);
    let was_timed_out = outcome.was_timed_out();
    let was_resource_limited = outcome.was_resource_limited();
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
        r#"{{"smiles":{smiles_json_value},"wasTimedOut":{was_timed_out},"wasResourceLimited":{was_resource_limited}}}"#
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
    if pairs.len() > WASM_MAX_BATCH_ITEMS {
        return Err(JsValue::from_str(&format!(
            "MMP result exceeds maximum pair count ({WASM_MAX_BATCH_ITEMS})"
        )));
    }

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

    bounded_json_output(format!("[{}]", entries.join(",")))
}
