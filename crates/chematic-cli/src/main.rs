//! Small, composable CLI for the common topology-bearing format bridge.

#![forbid(unsafe_code)]

use clap::{Args, Parser, Subcommand};
use std::fs;
use std::io::{self, Read, Write};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(
    name = "chematic",
    version,
    about = "Common chematic molecule utilities"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Args, Clone, Debug)]
struct BatchLimits {
    /// Maximum batch input size in bytes.
    #[arg(long, default_value_t = 64 * 1024 * 1024)]
    max_input_bytes: usize,
    /// Maximum number of non-empty, non-comment records.
    #[arg(long, default_value_t = 100_000)]
    max_records: usize,
    /// Maximum physical line size in bytes.
    #[arg(long, default_value_t = 1 * 1024 * 1024)]
    max_line_bytes: usize,
}

const MAX_OUTPUT_BYTES: usize = 64 * 1024 * 1024;
const MAX_INPUT_BYTES: usize = 64 * 1024 * 1024;
const MAX_SMILES_INPUT_BYTES: usize = 1 << 20;
const MAX_SMILES_ATOMS: usize = 10_000;

impl Default for BatchLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 * 1024 * 1024,
            max_records: 100_000,
            max_line_bytes: 1024 * 1024,
        }
    }
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Convert a topology-bearing molecule between common text formats.
    Convert {
        /// Input format: smiles, mol, mol_v3000, mol2, cml, cjson, moljson, or cdxml.
        #[arg(long)]
        input_format: String,
        /// Output format: smiles, mol, mol_v3000, mol2, cml, cjson, moljson, or cdxml.
        #[arg(long)]
        output_format: String,
        /// Read from this file instead of stdin.
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Write to this file instead of stdout.
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Maximum input size in bytes.
        #[arg(long, default_value_t = 64 * 1024 * 1024)]
        max_input_bytes: usize,
    },
    /// Parse and validate a SMILES molecule without descriptor calculation.
    Parse {
        /// SMILES to parse.
        smiles: String,
    },
    /// Calculate a compact JSON descriptor record from a SMILES string.
    Descriptors {
        /// SMILES to analyze.
        smiles: String,
    },
    /// Calculate a fingerprint and emit its set-bit indices as JSON.
    Fingerprint {
        /// SMILES to analyze.
        smiles: String,
        /// Algorithm: ecfp4, ecfp6, or maccs.
        #[arg(long, default_value = "ecfp4")]
        algorithm: String,
    },
    /// Compare two SMILES strings with fingerprint Tanimoto similarity.
    Similarity {
        /// First SMILES to analyze.
        smiles_a: String,
        /// Second SMILES to analyze.
        smiles_b: String,
        /// Algorithm: ecfp4, ecfp6, or maccs.
        #[arg(long, default_value = "ecfp4")]
        algorithm: String,
    },
    /// Find all SMARTS substructure matches in a SMILES molecule.
    Substructure {
        /// Molecule SMILES to search.
        smiles: String,
        /// SMARTS query.
        smarts: String,
    },
    /// Standardize a SMILES molecule and emit an auditable JSON report.
    Standardize {
        /// SMILES to standardize.
        smiles: String,
    },
    /// Generate a complete single-molecule analysis report as JSON.
    Report {
        /// SMILES to analyze.
        smiles: String,
    },
    /// Normalize a reaction SMILES and report its mapped reaction center.
    Reaction {
        /// Reaction SMILES in `reactants>agents>products` form.
        reaction_smiles: String,
    },
    /// Match a reaction against a reaction SMARTS query.
    ReactionMatch {
        /// Reaction SMILES in `reactants>agents>products` form.
        reaction_smiles: String,
        /// Reaction SMARTS in `reactant>>product` form.
        query: String,
    },
    /// Check atom balance for a reaction and report element counts.
    ReactionBalance {
        /// Reaction SMILES in `reactants>agents>products` form.
        reaction_smiles: String,
    },
    /// Generate an inspectable reaction fingerprint.
    ReactionFingerprint {
        /// Reaction SMILES in `reactants>agents>products` form.
        reaction_smiles: String,
        /// Combination mode: xor (transformation) or or (composition).
        #[arg(long, default_value = "xor")]
        mode: String,
    },
    /// Compare two reactions with reaction-fingerprint Tanimoto similarity.
    ReactionSimilarity {
        /// First reaction SMILES.
        reaction_a: String,
        /// Second reaction SMILES.
        reaction_b: String,
    },
    /// Process one SMILES per line and retain per-record errors in JSON.
    BatchReport {
        /// Read line-delimited SMILES from this file instead of stdin.
        #[arg(short, long)]
        input: Option<PathBuf>,
        #[command(flatten)]
        limits: BatchLimits,
    },
    /// Process one SMILES per line and return lightweight descriptor records.
    BatchDescriptors {
        /// Read line-delimited SMILES from this file instead of stdin.
        #[arg(short, long)]
        input: Option<PathBuf>,
        #[command(flatten)]
        limits: BatchLimits,
    },
    /// Process one SMILES per line and return fingerprint records.
    BatchFingerprints {
        /// Read line-delimited SMILES from this file instead of stdin.
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Algorithm: ecfp4, ecfp6, or maccs.
        #[arg(long, default_value = "ecfp4")]
        algorithm: String,
        #[command(flatten)]
        limits: BatchLimits,
    },
    /// Process one SMILES per line with the auditable standardization pipeline.
    BatchStandardize {
        /// Read line-delimited SMILES from this file instead of stdin.
        #[arg(short, long)]
        input: Option<PathBuf>,
        #[command(flatten)]
        limits: BatchLimits,
    },
    /// Compare one tab-separated pair of SMILES per line.
    BatchSimilarity {
        /// Read tab-separated pairs from this file instead of stdin.
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Algorithm: ecfp4, ecfp6, or maccs.
        #[arg(long, default_value = "ecfp4")]
        algorithm: String,
        #[command(flatten)]
        limits: BatchLimits,
    },
    /// Search one tab-separated SMILES/SMARTS pair per line.
    BatchSubstructure {
        /// Read tab-separated pairs from this file instead of stdin.
        #[arg(short, long)]
        input: Option<PathBuf>,
        #[command(flatten)]
        limits: BatchLimits,
    },
    /// Process one reaction SMILES per line and retain per-record errors.
    BatchReactions {
        /// Read line-delimited reaction SMILES from this file instead of stdin.
        #[arg(short, long)]
        input: Option<PathBuf>,
        #[command(flatten)]
        limits: BatchLimits,
    },
}

fn format_name(format: &str) -> Option<String> {
    match format
        .to_ascii_lowercase()
        .trim_start_matches('.')
        .replace('-', "_")
        .as_str()
    {
        "smi" | "smiles" => Some("smiles".to_string()),
        "mol" | "sdf" => Some("mol".to_string()),
        "v3000" | "mol_v3000" => Some("mol_v3000".to_string()),
        "mol2" => Some("mol2".to_string()),
        "cml" => Some("cml".to_string()),
        "cjson" => Some("cjson".to_string()),
        "moljson" => Some("moljson".to_string()),
        "cdxml" => Some("cdxml".to_string()),
        _ => None,
    }
}

fn convert_text(text: &str, input_format: &str, output_format: &str) -> Result<String, String> {
    let input = format_name(input_format)
        .ok_or_else(|| format!("unsupported input format: {input_format}"))?;
    let output = format_name(output_format)
        .ok_or_else(|| format!("unsupported output format: {output_format}"))?;
    let mol = match input.as_str() {
        "smiles" => chematic_smiles::parse(text).map_err(|e| e.to_string())?,
        "mol" => chematic_mol::parse_mol(text)
            .map(|(mol, _)| mol)
            .map_err(|e| e.to_string())?,
        "mol_v3000" => chematic_mol::parse_mol_v3000(text)
            .map(|(mol, _)| mol)
            .map_err(|e| e.to_string())?,
        "mol2" => chematic_mol::parse_mol2(text)
            .map(|(mol, _)| mol)
            .map_err(|e| e.to_string())?,
        "cml" => chematic_mol::parse_cml(text)
            .map(|(mol, _)| mol)
            .map_err(|e| e.to_string())?,
        "cjson" => chematic_mol::parse_cjson(text)
            .map(|(mol, _)| mol)
            .map_err(|e| e.to_string())?,
        "moljson" => chematic_mol::parse_moljson(text).map_err(|e| e.to_string())?,
        "cdxml" => chematic_mol::parse_cdxml(text)
            .map(|(mol, _)| mol)
            .map_err(|e| e.to_string())?,
        _ => unreachable!(),
    };
    if mol.atom_count() > MAX_SMILES_ATOMS {
        return Err(format!(
            "molecule exceeds maximum atom count ({MAX_SMILES_ATOMS})"
        ));
    }
    Ok(match output.as_str() {
        "smiles" => chematic_smiles::canonical_smiles(&mol),
        "mol" => chematic_mol::write_mol(&mol, &chematic_mol::MolMetadata::default()),
        "mol_v3000" => {
            chematic_mol::write_mol_v3000(&mol, &chematic_mol::MolMetadata::default(), &[])
        }
        "mol2" => chematic_mol::write_mol2(&mol, &[]),
        "cml" => chematic_mol::write_cml(&mol, None),
        "cjson" => chematic_mol::write_cjson(&mol, &[]),
        "moljson" => chematic_mol::write_moljson(&mol),
        "cdxml" => chematic_mol::write_cdxml(&mol, &[]),
        _ => unreachable!(),
    })
}

fn read_limited_input(path: Option<&PathBuf>, max_input_bytes: usize) -> Result<String, String> {
    if max_input_bytes > MAX_INPUT_BYTES {
        return Err(format!(
            "requested input limit exceeds CLI maximum ({MAX_INPUT_BYTES} bytes)"
        ));
    }
    let mut bytes = Vec::new();
    let read_result = match path {
        Some(path) => {
            let file = fs::File::open(path).map_err(|e| format!("read {}: {e}", path.display()))?;
            file.take(max_input_bytes.saturating_add(1) as u64)
                .read_to_end(&mut bytes)
        }
        None => io::stdin()
            .take(max_input_bytes.saturating_add(1) as u64)
            .read_to_end(&mut bytes),
    };
    read_result.map_err(|e| format!("read input: {e}"))?;
    if bytes.len() > max_input_bytes {
        return Err(format!(
            "input exceeds --max-input-bytes ({max_input_bytes})"
        ));
    }
    String::from_utf8(bytes).map_err(|e| format!("batch input is not UTF-8: {e}"))
}

fn parse_cli_smiles(smiles: &str) -> Result<chematic_core::Molecule, String> {
    if smiles.len() > MAX_SMILES_INPUT_BYTES {
        return Err(format!(
            "SMILES exceeds maximum input size ({} > {MAX_SMILES_INPUT_BYTES} bytes)",
            smiles.len()
        ));
    }
    let mol = chematic_smiles::parse(smiles).map_err(|e| e.to_string())?;
    if mol.atom_count() > MAX_SMILES_ATOMS {
        return Err(format!(
            "SMILES exceeds maximum atom count ({MAX_SMILES_ATOMS})"
        ));
    }
    Ok(mol)
}

fn batch_lines<'a>(text: &'a str, limits: &BatchLimits) -> Result<Vec<&'a str>, String> {
    if text.len() > limits.max_input_bytes {
        return Err(format!(
            "batch input exceeds --max-input-bytes ({})",
            limits.max_input_bytes
        ));
    }
    let mut lines = Vec::new();
    for line in text.lines() {
        if line.len() > limits.max_line_bytes {
            return Err(format!(
                "batch line exceeds --max-line-bytes ({})",
                limits.max_line_bytes
            ));
        }
        let line = line.trim();
        if !line.is_empty() && !line.starts_with('#') {
            if lines.len() == limits.max_records {
                return Err(format!(
                    "batch input exceeds --max-records ({})",
                    limits.max_records
                ));
            }
            lines.push(line);
        }
    }
    Ok(lines)
}

/// Add the stable envelope shared by every line-oriented batch command.
///
/// Records remain in input order and operation-specific fields remain at the
/// top level for backwards compatibility. The envelope makes a partial-result
/// manifest explicit: callers can persist the operation, applied limits, and
/// completion status without inferring them from record counts.
fn add_batch_manifest(value: &mut serde_json::Value, operation: &str, limits: &BatchLimits) {
    let Some(object) = value.as_object_mut() else {
        return;
    };
    object.insert("schema_version".to_string(), serde_json::json!(1));
    object.insert("operation".to_string(), serde_json::json!(operation));
    object.insert("status".to_string(), serde_json::json!("complete"));
    object.insert(
        "limits".to_string(),
        serde_json::json!({
            "max_input_bytes": limits.max_input_bytes,
            "max_records": limits.max_records,
            "max_line_bytes": limits.max_line_bytes,
        }),
    );
    if let Some(records) = object.get("records").and_then(serde_json::Value::as_array) {
        object.insert("record_count".to_string(), serde_json::json!(records.len()));
    }
}

fn write_output(path: Option<&PathBuf>, text: &str) -> Result<(), String> {
    if text.len() > MAX_OUTPUT_BYTES {
        return Err(format!(
            "output exceeds maximum size of {MAX_OUTPUT_BYTES} bytes"
        ));
    }
    match path {
        Some(path) => fs::write(path, text).map_err(|e| format!("write {}: {e}", path.display())),
        None => io::stdout()
            .write_all(text.as_bytes())
            .map_err(|e| format!("write stdout: {e}")),
    }
}

fn descriptors_json(smiles: &str) -> Result<String, String> {
    let mol = parse_cli_smiles(smiles)?;
    Ok(serde_json::json!({
        "smiles": chematic_smiles::canonical_smiles(&mol),
        "formula": chematic_chem::calc_mol_formula(&mol),
        "heavy_atoms": chematic_chem::heavy_atom_count(&mol),
        "molecular_weight": chematic_chem::molecular_weight(&mol),
        "exact_mass": chematic_chem::exact_mass(&mol),
        "logp": chematic_chem::logp_crippen(&mol),
        "tpsa": chematic_chem::tpsa(&mol),
        "hbd": chematic_chem::hbd_count(&mol),
        "hba": chematic_chem::hba_count(&mol),
        "rotatable_bonds": chematic_chem::rotatable_bond_count(&mol),
    })
    .to_string())
}

fn parse_json(smiles: &str) -> Result<String, String> {
    let mol = parse_cli_smiles(smiles)?;
    Ok(serde_json::json!({
        "input_smiles": smiles,
        "canonical_smiles": chematic_smiles::canonical_smiles(&mol),
        "formula": mol.total_formula(),
        "atoms": mol.atom_count(),
        "bonds": mol.bond_count(),
        "formal_charge": chematic_chem::formal_charge_sum(&mol),
    })
    .to_string())
}

fn fingerprint_json(smiles: &str, algorithm: &str) -> Result<String, String> {
    let mol = parse_cli_smiles(smiles)?;
    let algorithm = algorithm.to_ascii_lowercase();
    let fp = match algorithm.as_str() {
        "ecfp4" => chematic_fp::ecfp4(&mol),
        "ecfp6" => chematic_fp::ecfp6(&mol),
        "maccs" => chematic_fp::maccs(&mol),
        _ => return Err(format!("unsupported fingerprint algorithm: {algorithm}")),
    };
    let bits: Vec<usize> = (0..2048).filter(|&bit| fp.get(bit)).collect();
    Ok(serde_json::json!({
        "algorithm": algorithm,
        "n_bits": 2048,
        "set_bits": bits,
        "popcount": bits.len(),
        "smiles": chematic_smiles::canonical_smiles(&mol),
    })
    .to_string())
}

fn fingerprint_for(
    mol: &chematic_core::Molecule,
    algorithm: &str,
) -> Result<chematic_fp::BitVec2048, String> {
    match algorithm.to_ascii_lowercase().as_str() {
        "ecfp4" => Ok(chematic_fp::ecfp4(mol)),
        "ecfp6" => Ok(chematic_fp::ecfp6(mol)),
        "maccs" => Ok(chematic_fp::maccs(mol)),
        other => Err(format!("unsupported fingerprint algorithm: {other}")),
    }
}

fn similarity_json(smiles_a: &str, smiles_b: &str, algorithm: &str) -> Result<String, String> {
    let mol_a = parse_cli_smiles(smiles_a)?;
    let mol_b = parse_cli_smiles(smiles_b)?;
    let algorithm = algorithm.to_ascii_lowercase();
    let fp_a = fingerprint_for(&mol_a, &algorithm)?;
    let fp_b = fingerprint_for(&mol_b, &algorithm)?;
    Ok(serde_json::json!({
        "algorithm": algorithm,
        "similarity": fp_a.tanimoto(&fp_b),
        "smiles_a": chematic_smiles::canonical_smiles(&mol_a),
        "smiles_b": chematic_smiles::canonical_smiles(&mol_b),
    })
    .to_string())
}

fn substructure_json(smiles: &str, smarts: &str) -> Result<String, String> {
    let mol = parse_cli_smiles(smiles)?;
    let query = chematic_smarts::parse_smarts(smarts).map_err(|e| e.to_string())?;
    let matches = chematic_smarts::find_matches(&query, &mol);
    let matches: Vec<Vec<usize>> = matches
        .iter()
        .map(|mapping| {
            let mut pairs: Vec<(usize, usize)> = mapping
                .iter()
                .map(|(query, atom)| (*query, atom.0 as usize))
                .collect();
            pairs.sort_unstable_by_key(|(query, _)| *query);
            pairs.into_iter().map(|(_, atom)| atom).collect()
        })
        .collect();
    Ok(serde_json::json!({
        "smiles": chematic_smiles::canonical_smiles(&mol),
        "smarts": smarts,
        "match_count": matches.len(),
        "matches": matches,
    })
    .to_string())
}

fn standardize_json(smiles: &str) -> Result<String, String> {
    let mol = parse_cli_smiles(smiles)?;
    let input_smiles = chematic_smiles::canonical_smiles(&mol);
    let (standardized, report) = chematic_chem::StandardizationPipeline::default().run(&mol);
    let output_smiles = chematic_smiles::canonical_smiles(&standardized);
    let status = match report.status {
        chematic_chem::PipelineStatus::Unchanged => "unchanged",
        chematic_chem::PipelineStatus::Modified => "modified",
        chematic_chem::PipelineStatus::CompletedWithWarnings => "completed_with_warnings",
    };
    let steps: Vec<_> = report
        .steps
        .iter()
        .map(|step| {
            serde_json::json!({
                "step": step.step.as_str(),
                "enabled": step.enabled,
                "changed": step.changed,
                "before": {
                    "atoms": step.before.atoms,
                    "bonds": step.before.bonds,
                    "hash": step.before.hash,
                },
                "after": {
                    "atoms": step.after.atoms,
                    "bonds": step.after.bonds,
                    "hash": step.after.hash,
                },
            })
        })
        .collect();
    let warnings: Vec<_> = report
        .warnings
        .iter()
        .map(|warning| {
            serde_json::json!({
                "code": warning.code,
                "message": warning.message,
            })
        })
        .collect();
    Ok(serde_json::json!({
        "input_smiles": input_smiles,
        "output_smiles": output_smiles,
        "changed": report.changed(),
        "status": status,
        "input": {
            "atoms": report.input.atoms,
            "bonds": report.input.bonds,
            "hash": report.input.hash,
        },
        "output": {
            "atoms": report.output.atoms,
            "bonds": report.output.bonds,
            "hash": report.output.hash,
        },
        "steps": steps,
        "warnings": warnings,
    })
    .to_string())
}

fn report_json(smiles: &str) -> Result<String, String> {
    let _ = parse_cli_smiles(smiles)?;
    let report = chematic_chem::molecule_report(smiles).map_err(|e| e.to_string())?;
    serde_json::to_string(&report).map_err(|e| format!("serialize report: {e}"))
}

fn reaction_json(reaction_smiles: &str) -> Result<String, String> {
    let reaction = chematic_rxn::parse_reaction(reaction_smiles).map_err(|e| e.to_string())?;
    let canonical_part = |molecules: &[chematic_core::Molecule]| {
        molecules
            .iter()
            .map(chematic_smiles::canonical_smiles)
            .collect::<Vec<_>>()
            .join(".")
    };
    let reactants = canonical_part(&reaction.reactants);
    let agents = canonical_part(&reaction.agents);
    let products = canonical_part(&reaction.products);
    let center = chematic_rxn::find_reaction_center(&reaction);
    let pairs = |pairs: &[(chematic_core::AtomIdx, chematic_core::AtomIdx)]| {
        pairs
            .iter()
            .map(|(left, right)| vec![left.0 as usize, right.0 as usize])
            .collect::<Vec<_>>()
    };
    Ok(serde_json::json!({
        "reaction_smiles": format!("{reactants}>{agents}>{products}"),
        "reactants": reaction.reactants.len(),
        "agents": reaction.agents.len(),
        "products": reaction.products.len(),
        "mapped": !center.changed_atoms.is_empty()
            || !center.broken_bonds.is_empty()
            || !center.formed_bonds.is_empty(),
        "reaction_center": {
            "changed_atoms": center.changed_atoms.iter().map(|idx| idx.0 as usize).collect::<Vec<_>>(),
            "broken_bonds": pairs(&center.broken_bonds),
            "formed_bonds": pairs(&center.formed_bonds),
        },
    })
    .to_string())
}

fn reaction_match_json(reaction_smiles: &str, query: &str) -> Result<String, String> {
    let reaction = chematic_rxn::parse_reaction(reaction_smiles).map_err(|e| e.to_string())?;
    let parsed_query = chematic_rxn::parse_reaction_query(query)
        .map_err(|e| format!("invalid reaction SMARTS: {e}"))?;
    let matched = chematic_rxn::has_reaction_substructure_match(&reaction, &parsed_query);
    Ok(serde_json::json!({
        "reaction_smiles": chematic_rxn::write_reaction(&reaction),
        "query": query,
        "matched": matched,
    })
    .to_string())
}

fn reaction_balance_json(reaction_smiles: &str) -> Result<String, String> {
    let reaction = chematic_rxn::parse_reaction(reaction_smiles).map_err(|e| e.to_string())?;
    let balance = chematic_rxn::balance_check(&reaction);
    Ok(serde_json::json!({
        "reaction_smiles": chematic_rxn::write_reaction(&reaction),
        "balanced": balance.balanced,
        "reactant_formula": balance.reactant_formula,
        "product_formula": balance.product_formula,
        "differences": balance.diff(),
    })
    .to_string())
}

fn reaction_fingerprint_json(reaction_smiles: &str, mode: &str) -> Result<String, String> {
    let reaction = chematic_rxn::parse_reaction(reaction_smiles).map_err(|e| e.to_string())?;
    let mode = mode.to_ascii_lowercase();
    let use_xor = match mode.as_str() {
        "xor" => true,
        "or" => false,
        _ => return Err(format!("unsupported reaction fingerprint mode: {mode}")),
    };
    let fingerprint =
        chematic_fp::reaction_fp_with_config(&reaction, &chematic_fp::ReactionFpConfig { use_xor });
    let set_bits = |fp: &chematic_fp::BitVec2048| -> Vec<usize> {
        (0..2048).filter(|&bit| fp.get(bit)).collect()
    };
    let reactant_bits = set_bits(&fingerprint.reactant_fp);
    let product_bits = set_bits(&fingerprint.product_fp);
    let combined_bits = set_bits(&fingerprint.combined_fp);
    Ok(serde_json::json!({
        "reaction_smiles": chematic_rxn::write_reaction(&reaction),
        "mode": mode,
        "n_bits": 2048,
        "reactant_popcount": reactant_bits.len(),
        "product_popcount": product_bits.len(),
        "popcount": combined_bits.len(),
        "set_bits": combined_bits,
    })
    .to_string())
}

fn reaction_similarity_json(reaction_a: &str, reaction_b: &str) -> Result<String, String> {
    let first = chematic_rxn::parse_reaction(reaction_a).map_err(|e| e.to_string())?;
    let second = chematic_rxn::parse_reaction(reaction_b).map_err(|e| e.to_string())?;
    Ok(serde_json::json!({
        "similarity": chematic_fp::tanimoto_reaction_fp(&first, &second),
        "reaction_a": chematic_rxn::write_reaction(&first),
        "reaction_b": chematic_rxn::write_reaction(&second),
        "fingerprint": "reaction_ecfp4_xor",
    })
    .to_string())
}

fn batch_report_json(text: &str, limits: &BatchLimits) -> Result<String, String> {
    let smiles = batch_lines(text, limits)?;
    let refs: Vec<&str> = smiles.clone();
    let mut report = serde_json::to_value(chematic_chem::screen_smiles(&refs))
        .map_err(|e| format!("serialize batch report: {e}"))?;
    add_batch_manifest(&mut report, "report", limits);
    serde_json::to_string(&report).map_err(|e| format!("serialize batch report: {e}"))
}

fn batch_descriptors_json(text: &str, limits: &BatchLimits) -> Result<String, String> {
    let smiles = batch_lines(text, limits)?;
    let mut records = Vec::with_capacity(smiles.len());
    let mut valid_count = 0usize;
    for (input_index, input_smiles) in smiles.iter().enumerate() {
        match descriptors_json(input_smiles).and_then(|json| {
            serde_json::from_str::<serde_json::Value>(&json).map_err(|e| e.to_string())
        }) {
            Ok(descriptors) => {
                valid_count += 1;
                records.push(serde_json::json!({
                    "input_index": input_index,
                    "input_smiles": input_smiles,
                    "descriptors": descriptors,
                    "error": null,
                }));
            }
            Err(error) => records.push(serde_json::json!({
                "input_index": input_index,
                "input_smiles": input_smiles,
                "descriptors": null,
                "error": error,
            })),
        }
    }
    let mut output = serde_json::json!({
        "records": records,
        "valid_count": valid_count,
        "error_count": smiles.len() - valid_count,
    });
    add_batch_manifest(&mut output, "descriptors", limits);
    Ok(output.to_string())
}

fn batch_fingerprints_json(
    text: &str,
    algorithm: &str,
    limits: &BatchLimits,
) -> Result<String, String> {
    let smiles = batch_lines(text, limits)?;
    let mut records = Vec::with_capacity(smiles.len());
    let mut valid_count = 0usize;
    for (input_index, input_smiles) in smiles.iter().enumerate() {
        match fingerprint_json(input_smiles, algorithm).and_then(|json| {
            serde_json::from_str::<serde_json::Value>(&json).map_err(|e| e.to_string())
        }) {
            Ok(fingerprint) => {
                valid_count += 1;
                records.push(serde_json::json!({
                    "input_index": input_index,
                    "input_smiles": input_smiles,
                    "fingerprint": fingerprint,
                    "error": null,
                }));
            }
            Err(error) => records.push(serde_json::json!({
                "input_index": input_index,
                "input_smiles": input_smiles,
                "fingerprint": null,
                "error": error,
            })),
        }
    }
    let mut output = serde_json::json!({
        "algorithm": algorithm.to_ascii_lowercase(),
        "records": records,
        "valid_count": valid_count,
        "error_count": smiles.len() - valid_count,
    });
    add_batch_manifest(&mut output, "fingerprints", limits);
    Ok(output.to_string())
}

fn batch_standardize_json(text: &str, limits: &BatchLimits) -> Result<String, String> {
    let smiles = batch_lines(text, limits)?;
    let mut records = Vec::with_capacity(smiles.len());
    let mut valid_count = 0usize;
    for (input_index, input_smiles) in smiles.iter().enumerate() {
        match standardize_json(input_smiles).and_then(|json| {
            serde_json::from_str::<serde_json::Value>(&json).map_err(|e| e.to_string())
        }) {
            Ok(standardization) => {
                valid_count += 1;
                records.push(serde_json::json!({
                    "input_index": input_index,
                    "input_smiles": input_smiles,
                    "standardization": standardization,
                    "error": null,
                }));
            }
            Err(error) => records.push(serde_json::json!({
                "input_index": input_index,
                "input_smiles": input_smiles,
                "standardization": null,
                "error": error,
            })),
        }
    }
    let mut output = serde_json::json!({
        "records": records,
        "valid_count": valid_count,
        "error_count": smiles.len() - valid_count,
    });
    add_batch_manifest(&mut output, "standardize", limits);
    Ok(output.to_string())
}

fn batch_similarity_json(
    text: &str,
    algorithm: &str,
    limits: &BatchLimits,
) -> Result<String, String> {
    let lines = batch_lines(text, limits)?;
    let mut records = Vec::with_capacity(lines.len());
    let mut valid_count = 0usize;
    for (input_index, line) in lines.iter().enumerate() {
        let result = line
            .split_once('\t')
            .ok_or_else(|| "expected SMILES_A<TAB>SMILES_B".to_string())
            .and_then(|(smiles_a, smiles_b)| {
                similarity_json(smiles_a.trim(), smiles_b.trim(), algorithm)
            })
            .and_then(|json| {
                serde_json::from_str::<serde_json::Value>(&json).map_err(|e| e.to_string())
            });
        match result {
            Ok(similarity) => {
                valid_count += 1;
                records.push(serde_json::json!({
                    "input_index": input_index,
                    "input": line,
                    "similarity": similarity,
                    "error": null,
                }));
            }
            Err(error) => records.push(serde_json::json!({
                "input_index": input_index,
                "input": line,
                "similarity": null,
                "error": error,
            })),
        }
    }
    let mut output = serde_json::json!({
        "algorithm": algorithm.to_ascii_lowercase(),
        "records": records,
        "valid_count": valid_count,
        "error_count": lines.len() - valid_count,
    });
    add_batch_manifest(&mut output, "similarity", limits);
    Ok(output.to_string())
}

fn batch_substructure_json(text: &str, limits: &BatchLimits) -> Result<String, String> {
    let lines = batch_lines(text, limits)?;
    let mut records = Vec::with_capacity(lines.len());
    let mut valid_count = 0usize;
    for (input_index, line) in lines.iter().enumerate() {
        let result = line
            .split_once('\t')
            .ok_or_else(|| "expected SMILES<TAB>SMARTS".to_string())
            .and_then(|(smiles, smarts)| substructure_json(smiles.trim(), smarts.trim()))
            .and_then(|json| {
                serde_json::from_str::<serde_json::Value>(&json).map_err(|e| e.to_string())
            });
        match result {
            Ok(substructure) => {
                valid_count += 1;
                records.push(serde_json::json!({
                    "input_index": input_index,
                    "input": line,
                    "substructure": substructure,
                    "error": null,
                }));
            }
            Err(error) => records.push(serde_json::json!({
                "input_index": input_index,
                "input": line,
                "substructure": null,
                "error": error,
            })),
        }
    }
    let mut output = serde_json::json!({
        "records": records,
        "valid_count": valid_count,
        "error_count": lines.len() - valid_count,
    });
    add_batch_manifest(&mut output, "substructure", limits);
    Ok(output.to_string())
}

fn batch_reactions_json(text: &str, limits: &BatchLimits) -> Result<String, String> {
    let lines = batch_lines(text, limits)?;
    let mut records = Vec::with_capacity(lines.len());
    let mut valid_count = 0usize;
    for (input_index, line) in lines.iter().enumerate() {
        match reaction_json(line).and_then(|json| {
            serde_json::from_str::<serde_json::Value>(&json).map_err(|e| e.to_string())
        }) {
            Ok(reaction) => {
                valid_count += 1;
                records.push(serde_json::json!({
                    "input_index": input_index,
                    "input": line,
                    "reaction": reaction,
                    "error": null,
                }));
            }
            Err(error) => records.push(serde_json::json!({
                "input_index": input_index,
                "input": line,
                "reaction": null,
                "error": error,
            })),
        }
    }
    let mut output = serde_json::json!({
        "records": records,
        "valid_count": valid_count,
        "error_count": lines.len() - valid_count,
    });
    add_batch_manifest(&mut output, "reactions", limits);
    Ok(output.to_string())
}

fn run(cli: Cli) -> Result<(), String> {
    match cli.command {
        Command::Convert {
            input_format,
            output_format,
            input,
            output,
            max_input_bytes,
        } => {
            let text = read_limited_input(input.as_ref(), max_input_bytes)?;
            let converted = convert_text(&text, &input_format, &output_format)?;
            write_output(output.as_ref(), &converted)
        }
        Command::Descriptors { smiles } => {
            let json = descriptors_json(&smiles)?;
            write_output(None, &format!("{json}\n"))
        }
        Command::Parse { smiles } => {
            let json = parse_json(&smiles)?;
            write_output(None, &format!("{json}\n"))
        }
        Command::Fingerprint { smiles, algorithm } => {
            let json = fingerprint_json(&smiles, &algorithm)?;
            write_output(None, &format!("{json}\n"))
        }
        Command::Similarity {
            smiles_a,
            smiles_b,
            algorithm,
        } => {
            let json = similarity_json(&smiles_a, &smiles_b, &algorithm)?;
            write_output(None, &format!("{json}\n"))
        }
        Command::Substructure { smiles, smarts } => {
            let json = substructure_json(&smiles, &smarts)?;
            write_output(None, &format!("{json}\n"))
        }
        Command::Standardize { smiles } => {
            let json = standardize_json(&smiles)?;
            write_output(None, &format!("{json}\n"))
        }
        Command::Report { smiles } => {
            let json = report_json(&smiles)?;
            write_output(None, &format!("{json}\n"))
        }
        Command::Reaction { reaction_smiles } => {
            let json = reaction_json(&reaction_smiles)?;
            write_output(None, &format!("{json}\n"))
        }
        Command::ReactionMatch {
            reaction_smiles,
            query,
        } => {
            let json = reaction_match_json(&reaction_smiles, &query)?;
            write_output(None, &format!("{json}\n"))
        }
        Command::ReactionBalance { reaction_smiles } => {
            let json = reaction_balance_json(&reaction_smiles)?;
            write_output(None, &format!("{json}\n"))
        }
        Command::ReactionFingerprint {
            reaction_smiles,
            mode,
        } => {
            let json = reaction_fingerprint_json(&reaction_smiles, &mode)?;
            write_output(None, &format!("{json}\n"))
        }
        Command::ReactionSimilarity {
            reaction_a,
            reaction_b,
        } => {
            let json = reaction_similarity_json(&reaction_a, &reaction_b)?;
            write_output(None, &format!("{json}\n"))
        }
        Command::BatchReport { input, limits } => {
            let text = read_limited_input(input.as_ref(), limits.max_input_bytes)?;
            let json = batch_report_json(&text, &limits)?;
            write_output(None, &format!("{json}\n"))
        }
        Command::BatchDescriptors { input, limits } => {
            let text = read_limited_input(input.as_ref(), limits.max_input_bytes)?;
            let json = batch_descriptors_json(&text, &limits)?;
            write_output(None, &format!("{json}\n"))
        }
        Command::BatchFingerprints {
            input,
            algorithm,
            limits,
        } => {
            let text = read_limited_input(input.as_ref(), limits.max_input_bytes)?;
            let json = batch_fingerprints_json(&text, &algorithm, &limits)?;
            write_output(None, &format!("{json}\n"))
        }
        Command::BatchStandardize { input, limits } => {
            let text = read_limited_input(input.as_ref(), limits.max_input_bytes)?;
            let json = batch_standardize_json(&text, &limits)?;
            write_output(None, &format!("{json}\n"))
        }
        Command::BatchSimilarity {
            input,
            algorithm,
            limits,
        } => {
            let text = read_limited_input(input.as_ref(), limits.max_input_bytes)?;
            let json = batch_similarity_json(&text, &algorithm, &limits)?;
            write_output(None, &format!("{json}\n"))
        }
        Command::BatchSubstructure { input, limits } => {
            let text = read_limited_input(input.as_ref(), limits.max_input_bytes)?;
            let json = batch_substructure_json(&text, &limits)?;
            write_output(None, &format!("{json}\n"))
        }
        Command::BatchReactions { input, limits } => {
            let text = read_limited_input(input.as_ref(), limits.max_input_bytes)?;
            let json = batch_reactions_json(&text, &limits)?;
            write_output(None, &format!("{json}\n"))
        }
    }
}

fn main() {
    if let Err(error) = run(Cli::parse()) {
        eprintln!("chematic: {error}");
        std::process::exit(2);
    }
}

#[cfg(test)]
mod tests {
    use super::{
        BatchLimits, MAX_OUTPUT_BYTES, MAX_SMILES_ATOMS, MAX_SMILES_INPUT_BYTES,
        batch_descriptors_json, batch_fingerprints_json, batch_reactions_json, batch_report_json,
        batch_similarity_json, batch_standardize_json, batch_substructure_json, convert_text,
        descriptors_json, fingerprint_json, parse_cli_smiles, parse_json, reaction_balance_json,
        reaction_fingerprint_json, reaction_json, reaction_match_json, reaction_similarity_json,
        report_json, similarity_json, standardize_json, substructure_json,
    };

    fn default_batch_limits() -> BatchLimits {
        BatchLimits::default()
    }

    #[test]
    fn converts_smiles_to_mol2_and_back() {
        let mol2 = convert_text("CCO", "smiles", "mol2").unwrap();
        assert!(mol2.contains("@<TRIPOS>MOLECULE"));
        let smiles = convert_text(&mol2, ".mol2", "smi").unwrap();
        let reparsed = chematic_smiles::parse(&smiles).unwrap();
        let expected = chematic_smiles::parse("CCO").unwrap();
        assert_eq!(reparsed.atom_count(), expected.atom_count());
        assert_eq!(reparsed.bond_count(), expected.bond_count());
    }

    #[test]
    fn rejects_unknown_formats() {
        assert!(
            convert_text("CCO", "smiles", "nope")
                .unwrap_err()
                .contains("unsupported output format")
        );
    }

    #[test]
    fn descriptors_are_json_and_include_core_fields() {
        let json: serde_json::Value =
            serde_json::from_str(&descriptors_json("CCO").unwrap()).unwrap();
        assert_eq!(json["heavy_atoms"], 3);
        assert!(json["molecular_weight"].as_f64().unwrap() > 40.0);
        assert_eq!(json["formula"], "C2H6O");
    }

    #[test]
    fn parse_reports_structure_without_expensive_descriptors() {
        let json: serde_json::Value =
            serde_json::from_str(&parse_json("C[NH3+]").unwrap()).unwrap();
        assert_eq!(json["input_smiles"], "C[NH3+]");
        assert_eq!(json["canonical_smiles"], "C[NH3+]");
        assert_eq!(json["atoms"], 2);
        assert_eq!(json["bonds"], 1);
        assert_eq!(json["formal_charge"], 1);
    }

    #[test]
    fn parse_rejects_invalid_smiles() {
        assert!(parse_json("C1CC").is_err());
    }

    #[test]
    fn fingerprint_json_is_stable_and_reports_set_bits() {
        let json: serde_json::Value =
            serde_json::from_str(&fingerprint_json("CCO", "ecfp4").unwrap()).unwrap();
        assert_eq!(json["algorithm"], "ecfp4");
        assert_eq!(json["n_bits"], 2048);
        assert!(json["popcount"].as_u64().unwrap() > 0);
        assert_eq!(
            json["set_bits"].as_array().unwrap().len(),
            json["popcount"].as_u64().unwrap() as usize
        );
    }

    #[test]
    fn fingerprint_rejects_unknown_algorithms() {
        assert!(
            fingerprint_json("CCO", "unknown")
                .unwrap_err()
                .contains("unsupported fingerprint algorithm")
        );
    }

    #[test]
    fn similarity_is_one_for_identical_molecules() {
        let json: serde_json::Value =
            serde_json::from_str(&similarity_json("CCO", "CCO", "ecfp4").unwrap()).unwrap();
        assert_eq!(json["similarity"], 1.0);
        assert_eq!(json["algorithm"], "ecfp4");
    }

    #[test]
    fn similarity_rejects_unknown_algorithms() {
        assert!(
            similarity_json("CCO", "CCN", "unknown")
                .unwrap_err()
                .contains("unsupported fingerprint algorithm")
        );
    }

    #[test]
    fn substructure_reports_sorted_atom_mappings() {
        let json: serde_json::Value =
            serde_json::from_str(&substructure_json("c1ccccc1", "c").unwrap()).unwrap();
        assert_eq!(json["match_count"], 6);
        assert_eq!(
            json["matches"].as_array().unwrap()[0]
                .as_array()
                .unwrap()
                .len(),
            1
        );
    }

    #[test]
    fn substructure_rejects_invalid_smarts() {
        assert!(
            substructure_json("CCO", "[")
                .unwrap_err()
                .contains("SMARTS")
        );
    }

    #[test]
    fn standardize_reports_pipeline_and_canonical_output() {
        let json: serde_json::Value =
            serde_json::from_str(&standardize_json("C[NH3+]").unwrap()).unwrap();
        assert_eq!(json["input_smiles"], "C[NH3+]");
        assert_eq!(json["status"], "modified");
        assert!(json["changed"].as_bool().unwrap());
        assert_eq!(json["steps"].as_array().unwrap().len(), 6);
        assert!(
            json["steps"]
                .as_array()
                .unwrap()
                .iter()
                .any(|step| step["step"] == "neutralize_charges" && step["changed"] == true)
        );
        assert!(json["warnings"].is_array());
    }

    #[test]
    fn report_emits_complete_machine_readable_analysis() {
        let json: serde_json::Value =
            serde_json::from_str(&report_json("CC(=O)Oc1ccccc1C(=O)O").unwrap()).unwrap();
        assert_eq!(json["input_smiles"], "CC(=O)Oc1ccccc1C(=O)O");
        assert!(json["canonical_smiles"].as_str().is_some());
        assert!(json["descriptors"]["molecular_weight"].as_f64().unwrap() > 100.0);
        assert!(json["filters"]["pains_alerts"].is_array());
        assert!(json["functional_groups"].is_array());
    }

    #[test]
    fn report_rejects_invalid_smiles() {
        assert!(report_json("C1CC").is_err());
    }

    #[test]
    fn reaction_normalizes_components_and_reports_mapped_center() {
        let json: serde_json::Value =
            serde_json::from_str(&reaction_json("[CH3:1][OH:2]>>[CH3:1][O-:2]").unwrap()).unwrap();
        assert_eq!(json["reaction_smiles"], "[OH:2][CH3:1]>>[CH3:1][O-:2]");
        assert_eq!(json["reactants"], 1);
        assert_eq!(json["products"], 1);
        assert_eq!(json["mapped"], true);
        assert_eq!(
            json["reaction_center"]["changed_atoms"],
            serde_json::json!([1])
        );
    }

    #[test]
    fn reaction_rejects_missing_arrow() {
        assert!(
            reaction_json("CCO")
                .unwrap_err()
                .contains("reaction SMILES")
        );
    }

    #[test]
    fn reaction_match_reports_boolean_without_inference() {
        let json: serde_json::Value =
            serde_json::from_str(&reaction_match_json("CCO>>CCO", "[#6]>>[#6]").unwrap()).unwrap();
        assert_eq!(json["matched"], true);
        assert_eq!(json["query"], "[#6]>>[#6]");
    }

    #[test]
    fn reaction_match_rejects_invalid_query() {
        assert!(
            reaction_match_json("CCO>>CCO", "[#6]")
                .unwrap_err()
                .contains("invalid reaction SMARTS")
        );
    }

    #[test]
    fn reaction_balance_reports_counts_and_differences() {
        let balanced: serde_json::Value =
            serde_json::from_str(&reaction_balance_json("CO.CO>>COC.O").unwrap()).unwrap();
        assert_eq!(balanced["balanced"], true);
        assert!(balanced["differences"].as_array().unwrap().is_empty());

        let unbalanced: serde_json::Value =
            serde_json::from_str(&reaction_balance_json("C>>CC").unwrap()).unwrap();
        assert_eq!(unbalanced["balanced"], false);
        assert!(
            unbalanced["differences"]
                .as_array()
                .unwrap()
                .iter()
                .any(|difference| difference.as_str().unwrap().contains("C"))
        );
    }

    #[test]
    fn reaction_fingerprint_reports_transformation_bits() {
        let json: serde_json::Value =
            serde_json::from_str(&reaction_fingerprint_json("CCO>>CC=O", "xor").unwrap()).unwrap();
        assert_eq!(json["mode"], "xor");
        assert_eq!(json["n_bits"], 2048);
        assert!(json["popcount"].as_u64().unwrap() > 0);
        assert_eq!(
            json["set_bits"].as_array().unwrap().len(),
            json["popcount"].as_u64().unwrap() as usize
        );
    }

    #[test]
    fn reaction_fingerprint_rejects_unknown_mode() {
        assert!(
            reaction_fingerprint_json("CCO>>CCO", "bad")
                .unwrap_err()
                .contains("unsupported reaction fingerprint mode")
        );
    }

    #[test]
    fn reaction_similarity_is_one_for_identical_reactions() {
        let json: serde_json::Value =
            serde_json::from_str(&reaction_similarity_json("CCO>>CC=O", "CCO>>CC=O").unwrap())
                .unwrap();
        assert_eq!(json["similarity"], 1.0);
        assert_eq!(json["fingerprint"], "reaction_ecfp4_xor");
    }

    #[test]
    fn reaction_similarity_rejects_invalid_reaction() {
        assert!(reaction_similarity_json("CCO", "CCO>>CCO").is_err());
    }

    #[test]
    fn batch_report_retains_partial_errors_in_input_order() {
        let json: serde_json::Value = serde_json::from_str(
            &batch_report_json("CCO\nC1CC\n# comment\nCCN\n", &default_batch_limits()).unwrap(),
        )
        .unwrap();
        assert_eq!(json["schema_version"], 1);
        assert_eq!(json["operation"], "report");
        assert_eq!(json["status"], "complete");
        assert_eq!(json["record_count"], 3);
        let records = json["records"].as_array().unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records[0]["input_index"], 0);
        assert!(records[0]["report"].is_object());
        assert!(records[1]["error"].as_str().is_some());
        assert_eq!(records[2]["input_index"], 2);
        assert!(records[2]["report"].is_object());
    }

    #[test]
    fn batch_descriptors_returns_lightweight_partial_manifest() {
        let json: serde_json::Value = serde_json::from_str(
            &batch_descriptors_json("CCO\nC1CC\n# comment\nCCN\n", &default_batch_limits())
                .unwrap(),
        )
        .unwrap();
        assert_eq!(json["valid_count"], 2);
        assert_eq!(json["error_count"], 1);
        assert_eq!(json["schema_version"], 1);
        assert_eq!(json["operation"], "descriptors");
        assert_eq!(json["status"], "complete");
        assert_eq!(json["record_count"], 3);
        assert_eq!(
            json["limits"]["max_records"],
            default_batch_limits().max_records
        );
        let records = json["records"].as_array().unwrap();
        assert!(records[0]["descriptors"].is_object());
        assert!(records[1]["error"].as_str().is_some());
        assert!(records[2]["descriptors"].is_object());
    }

    #[test]
    fn batch_fingerprints_returns_algorithm_and_partial_manifest() {
        let json: serde_json::Value = serde_json::from_str(
            &batch_fingerprints_json("CCO\nC1CC\nCCN\n", "ecfp4", &default_batch_limits()).unwrap(),
        )
        .unwrap();
        assert_eq!(json["algorithm"], "ecfp4");
        assert_eq!(json["valid_count"], 2);
        assert_eq!(json["error_count"], 1);
        assert!(json["records"][0]["fingerprint"]["set_bits"].is_array());
        assert!(json["records"][1]["error"].as_str().is_some());
    }

    #[test]
    fn batch_fingerprints_rejects_unknown_algorithm_per_record() {
        let json: serde_json::Value = serde_json::from_str(
            &batch_fingerprints_json("CCO\nCCN\n", "bad", &default_batch_limits()).unwrap(),
        )
        .unwrap();
        assert_eq!(json["valid_count"], 0);
        assert_eq!(json["error_count"], 2);
    }

    #[test]
    fn batch_standardize_retains_audit_reports_and_errors() {
        let json: serde_json::Value = serde_json::from_str(
            &batch_standardize_json("C[NH3+]\nC1CC\nCCO\n", &default_batch_limits()).unwrap(),
        )
        .unwrap();
        assert_eq!(json["valid_count"], 2);
        assert_eq!(json["error_count"], 1);
        assert!(json["records"][0]["standardization"]["steps"].is_array());
        assert!(json["records"][1]["error"].as_str().is_some());
        assert_eq!(json["records"][2]["standardization"]["status"], "unchanged");
    }

    #[test]
    fn batch_similarity_processes_tsv_pairs_and_retains_errors() {
        let json: serde_json::Value = serde_json::from_str(
            &batch_similarity_json(
                "CCO\tCCO\nCCO\tC1CC\nmissing\n",
                "ecfp4",
                &default_batch_limits(),
            )
            .unwrap(),
        )
        .unwrap();
        assert_eq!(json["algorithm"], "ecfp4");
        assert_eq!(json["valid_count"], 1);
        assert_eq!(json["error_count"], 2);
        assert_eq!(json["records"][0]["similarity"]["similarity"], 1.0);
        assert!(json["records"][1]["error"].as_str().is_some());
        assert!(json["records"][2]["error"].as_str().is_some());
    }

    #[test]
    fn batch_substructure_processes_tsv_queries_and_retains_errors() {
        let json: serde_json::Value = serde_json::from_str(
            &batch_substructure_json("c1ccccc1\tc\nCCO\t[\n", &default_batch_limits()).unwrap(),
        )
        .unwrap();
        assert_eq!(json["valid_count"], 1);
        assert_eq!(json["error_count"], 1);
        assert_eq!(json["records"][0]["substructure"]["match_count"], 6);
        assert!(json["records"][1]["error"].as_str().is_some());
    }

    #[test]
    fn batch_reactions_returns_normalized_records_and_errors() {
        let json: serde_json::Value = serde_json::from_str(
            &batch_reactions_json("CCO>>CCO\nCCO\n", &default_batch_limits()).unwrap(),
        )
        .unwrap();
        assert_eq!(json["valid_count"], 1);
        assert_eq!(json["error_count"], 1);
        assert_eq!(json["records"][0]["reaction"]["reactants"], 1);
        assert!(json["records"][1]["error"].as_str().is_some());
    }

    #[test]
    fn batch_limits_reject_oversized_input_and_records() {
        let mut limits = default_batch_limits();
        limits.max_input_bytes = 3;
        assert!(
            batch_reactions_json("CCO\n", &limits)
                .unwrap_err()
                .contains("max-input-bytes")
        );

        limits = default_batch_limits();
        limits.max_records = 1;
        assert!(
            batch_reactions_json("CCO>>CCO\nCCN>>CCN\n", &limits)
                .unwrap_err()
                .contains("max-records")
        );

        limits = default_batch_limits();
        limits.max_line_bytes = 3;
        assert!(
            batch_reactions_json("CCO>>CCO\n", &limits)
                .unwrap_err()
                .contains("max-line-bytes")
        );
    }

    #[test]
    fn single_smiles_contract_rejects_oversized_bytes_before_parse() {
        let input = "C".repeat(MAX_SMILES_INPUT_BYTES + 1);
        let error = parse_cli_smiles(&input)
            .err()
            .expect("oversized input must fail");
        assert!(error.contains("maximum input size"));
    }

    #[test]
    fn single_smiles_contract_rejects_oversized_molecules() {
        let input = "C".repeat(MAX_SMILES_ATOMS + 1);
        let error = parse_cli_smiles(&input)
            .err()
            .expect("oversized molecule must fail");
        assert!(error.contains("maximum atom count"));
    }

    #[test]
    fn output_limit_is_explicit_and_bounded() {
        let oversized = "x".repeat(MAX_OUTPUT_BYTES + 1);
        assert!(
            super::write_output(Some(&std::path::PathBuf::from("/dev/null")), &oversized)
                .unwrap_err()
                .contains("maximum size")
        );
    }
}
