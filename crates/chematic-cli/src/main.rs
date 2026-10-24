//! Small, composable CLI for the common topology-bearing format bridge.

#![forbid(unsafe_code)]

use clap::{Args, Parser, Subcommand};
use std::fs;
use std::io::{self, Read, Write};
use std::path::PathBuf;

mod conversion_format;
mod dispatch;

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
    /// Preserve standard nonperiodic BIOGRF 200 atom and connectivity fields.
    Bgf {
        #[arg(short, long)]
        input: Option<PathBuf>,
        #[arg(short, long)]
        output: Option<PathBuf>,
        #[arg(long, conflicts_with = "from_json")]
        inspect: bool,
        #[arg(long, conflicts_with = "inspect")]
        from_json: bool,
        #[arg(long, default_value_t = 16 * 1024 * 1024)]
        max_input_bytes: usize,
    },
    /// Retain all explicit Cartesian GJF/COM Link1 jobs, or inspect them as JSON.
    GjfDocument {
        #[arg(short, long)]
        input: Option<PathBuf>,
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Emit ordered geometry/source ranges instead of exact source text.
        #[arg(long)]
        inspect: bool,
        #[arg(long, default_value_t = 64 * 1024 * 1024)]
        max_input_bytes: usize,
        /// Whole-document job limit, at most 128 (atom limit: 10000 total).
        #[arg(long, default_value_t = 128)]
        max_jobs: usize,
    },
    /// Emit a checked conversion report (MOL2, CJSON, or same-format CML input).
    ConvertWithReport {
        #[arg(long)]
        input_format: String,
        #[arg(long)]
        output_format: String,
        /// strict, allow_loss, or preserve_source.
        #[arg(long, default_value = "strict")]
        policy: String,
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Write the JSON report here; molecule output is only inside the report.
        #[arg(short, long)]
        output: Option<PathBuf>,
        #[arg(long, default_value_t = 16 * 1024 * 1024)]
        max_input_bytes: usize,
    },
    /// Convert common text formats; GJF/COM supports same-format source retention.
    Convert {
        /// Input format: smiles, mol, mol_v3000, mol2, cml, cjson, moljson, cdxml; gjf/com for same-format source retention.
        #[arg(long)]
        input_format: String,
        /// Output format: smiles, mol, mol_v3000, mol2, cml, cjson, moljson, cdxml; gjf/com for same-format source retention.
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
    /// Parse and re-serialize supported CXSMILES metadata as an auditable JSON record.
    Cxsmiles {
        /// CXSMILES to parse.
        cxsmiles: String,
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
        /// Algorithm: ecfp4, ecfp6, maccs, rdkit_morgan, rdkit_pattern, rdkit_rdk, or rdkit_torsion.
        #[arg(long, default_value = "ecfp4")]
        algorithm: String,
    },
    /// Compare two SMILES strings with fingerprint Tanimoto similarity.
    Similarity {
        /// First SMILES to analyze.
        smiles_a: String,
        /// Second SMILES to analyze.
        smiles_b: String,
        /// Algorithm: ecfp4, ecfp6, maccs, rdkit_morgan, rdkit_pattern, rdkit_rdk, or rdkit_torsion.
        #[arg(long, default_value = "ecfp4")]
        algorithm: String,
    },
    /// Find all SMARTS substructure matches in a SMILES molecule.
    Substructure {
        /// Molecule SMILES to search.
        smiles: String,
        /// SMARTS query.
        smarts: String,
        /// Enforce atom and bond stereochemistry in the SMARTS query.
        #[arg(long)]
        use_chirality: bool,
        /// Use the experimental RDKit-aligned aromaticity and ring-query model.
        #[arg(long)]
        rdkit_parity: bool,
    },
    /// Standardize a SMILES molecule and emit an auditable JSON report.
    Standardize {
        /// SMILES to standardize.
        smiles: String,
        /// Keep only the largest connected fragment (Phase 1 fragment-selection profile).
        #[arg(long)]
        largest_fragment_only: bool,
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
    /// Canonicalize one SMILES per line without calculating descriptors.
    BatchCanonicalize {
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
        /// Algorithm: ecfp4, ecfp6, maccs, rdkit_morgan, rdkit_pattern, rdkit_rdk, or rdkit_torsion.
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
        /// Keep only the largest connected fragment (Phase 1 fragment-selection profile).
        #[arg(long)]
        largest_fragment_only: bool,
        #[command(flatten)]
        limits: BatchLimits,
    },
    /// Compare one tab-separated pair of SMILES per line.
    BatchSimilarity {
        /// Read tab-separated pairs from this file instead of stdin.
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Algorithm: ecfp4, ecfp6, maccs, rdkit_morgan, rdkit_pattern, rdkit_rdk, or rdkit_torsion.
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
        /// Enforce atom and bond stereochemistry in each SMARTS query.
        #[arg(long)]
        use_chirality: bool,
        /// Use the experimental RDKit-aligned aromaticity and ring-query model.
        #[arg(long)]
        rdkit_parity: bool,
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

fn convert_text(text: &str, input_format: &str, output_format: &str) -> Result<String, String> {
    let input = conversion_format::format_name(input_format)
        .ok_or_else(|| format!("unsupported input format: {input_format}"))?;
    let output = conversion_format::format_name(output_format)
        .ok_or_else(|| format!("unsupported output format: {output_format}"))?;

    if let Some(result) = conversion_format::source_format_conversion(text, &input, &output) {
        return result;
    }
    if input == "mol2" && output == "mol2" {
        let record = chematic_mol::parse_mol2_source_record_with_limits(
            text,
            &chematic_mol::Mol2ParseLimits {
                max_atoms: MAX_SMILES_ATOMS,
                ..Default::default()
            },
        )
        .map_err(|error| error.to_string())?;
        return Ok(chematic_mol::write_mol2_source_record(&record));
    }

    // A V3000-to-V3000 conversion is a file-format round trip, not merely a
    // topology conversion.  The generic path below intentionally discards
    // layout metadata because most output formats cannot represent it, but
    // doing that here rewrites every atom at (0, 0, 0) and silently loses
    // wedge/direction-derived stereochemistry on re-read.  Preserve the
    // parsed coordinates and V3000 metadata when the target can carry them.
    if input == "mol_v3000" && output == "mol_v3000" {
        let (mol, metadata, coords) =
            chematic_mol::parse_mol_v3000_with_coords(text).map_err(|e| e.to_string())?;
        if mol.atom_count() > MAX_SMILES_ATOMS {
            return Err(format!(
                "molecule exceeds maximum atom count ({MAX_SMILES_ATOMS})"
            ));
        }
        return Ok(chematic_mol::write_mol_v3000(&mol, &metadata, &coords));
    }
    let mol = match input.as_str() {
        "smiles" => chematic_smiles::parse(text).map_err(|e| e.to_string())?,
        "mol" => chematic_mol::parse_mol(text)
            .map(|(mol, _)| mol)
            .map_err(|e| e.to_string())?,
        "mol_v3000" => {
            let (mol, metadata) = chematic_mol::parse_mol_v3000(text).map_err(|e| e.to_string())?;
            if v3000_has_unmodeled_query_semantics(&mol, &metadata) {
                return Err(format!(
                    "unsupported V3000 query semantics for output format '{output}'; use mol_v3000 for opaque round-trip preservation"
                ));
            }
            mol
        }
        "mol2" => chematic_mol::parse_mol2(text)
            .map(|(mol, _)| mol)
            .map_err(|e| e.to_string())?,
        "cml" => chematic_mol::parse_cml(text)
            .map(|(mol, _)| mol)
            .map_err(|e| e.to_string())?,
        "cjson" => chematic_mol::parse_cjson_checked_with_limits(
            text,
            &chematic_mol::CjsonParseLimits {
                max_atoms: MAX_SMILES_ATOMS,
                ..Default::default()
            },
        )
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
        // CTAB writers use Kekule bonds when aromatic bond type 4 would
        // discard an aromatic heteroatom's or charged carbon's hydrogen.
        // Recover the aromatic SMILES spelling after reading either CTAB;
        // leave the parsed graph unchanged for other output formats.
        "smiles" if matches!(input.as_str(), "mol" | "mol_v3000") => {
            chematic_smiles::canonical_smiles(&chematic_perception::apply_aromaticity(&mol))
        }
        "smiles" => chematic_smiles::canonical_smiles(&mol),
        "mol" => chematic_mol::write_mol(&mol, &chematic_mol::MolMetadata::default()),
        "mol_v3000" => {
            chematic_mol::write_mol_v3000(&mol, &chematic_mol::MolMetadata::default(), &[])
        }
        "mol2" => chematic_mol::write_mol2(&mol, &[]),
        "cml" => chematic_mol::write_cml(&mol, None),
        "cjson" => {
            return chematic_mol::write_cjson_checked(&mol, &[]).map_err(|error| error.to_string());
        }
        "moljson" => chematic_mol::write_moljson(&mol),
        "cdxml" => chematic_mol::write_cdxml(&mol, &[]),
        _ => unreachable!(),
    })
}

/// True when a V3000 record carries a query constraint outside chematic's
/// ordinary-molecule model.  The V3000 -> V3000 path retains those attributes
/// verbatim; conversion into an ordinary molecule format must instead refuse
/// so callers cannot mistake a weakened query for a concrete structure.
fn v3000_has_unmodeled_query_semantics(
    mol: &chematic_core::Molecule,
    metadata: &chematic_mol::MolMetadata,
) -> bool {
    use chematic_core::BondOrder;

    let query_bond = mol.bonds().any(|(_, bond)| {
        matches!(
            bond.order,
            BondOrder::QuerySingleOrDouble
                | BondOrder::QuerySingleOrAromatic
                | BondOrder::QueryDoubleOrAromatic
                | BondOrder::QueryAny
        )
    });
    let query_atom_attribute = metadata
        .v3000_atom_properties
        .iter()
        .any(|(_, properties)| {
            [
                "SUBST=",
                "RBCNT=",
                "UNSAT=",
                "RINGBONDCNT=",
                "ATTCHPT=",
                "RGROUPS=",
            ]
            .iter()
            .any(|key| properties.contains(key))
        });
    query_bond || query_atom_attribute
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
    if let Some(record_count) = object
        .get("records")
        .and_then(serde_json::Value::as_array)
        .map(Vec::len)
    {
        object.insert("record_count".to_string(), serde_json::json!(record_count));
        object.insert("input_count".to_string(), serde_json::json!(record_count));
    }
    // Keep the original valid_count/error_count fields for compatibility,
    // while exposing the shared BatchResult accounting names for callers that
    // need to prove every input has a terminal outcome.
    if let (Some(accepted), Some(rejected)) = (
        object
            .get("valid_count")
            .and_then(serde_json::Value::as_u64),
        object
            .get("error_count")
            .and_then(serde_json::Value::as_u64),
    ) {
        object.insert("accepted_count".to_string(), serde_json::json!(accepted));
        object.insert("rejected_count".to_string(), serde_json::json!(rejected));
        object.insert("refused_count".to_string(), serde_json::json!(0));
        object.insert("skipped_count".to_string(), serde_json::json!(0));
        object.insert(
            "all_succeeded".to_string(),
            serde_json::json!(rejected == 0),
        );
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
        "rdkit_molecular_weight": chematic_chem::rdkit_molecular_weight(&mol),
        "exact_mass": chematic_chem::exact_mass(&mol),
        "logp": chematic_chem::logp_crippen(&mol),
        "molar_refractivity": chematic_chem::molar_refractivity(&mol),
        "tpsa": chematic_chem::tpsa(&mol),
        "hbd": chematic_chem::hbd_count(&mol),
        "hba": chematic_chem::hba_count(&mol),
        "fsp3": chematic_chem::fsp3(&mol),
        "aromatic_ring_count": chematic_chem::aromatic_ring_count(&mol),
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

fn cxsmiles_json(cxsmiles: &str) -> Result<String, String> {
    if cxsmiles.len() > MAX_SMILES_INPUT_BYTES {
        return Err(format!(
            "CXSMILES input exceeds {MAX_SMILES_INPUT_BYTES} byte limit ({} bytes)",
            cxsmiles.len()
        ));
    }
    let cx = chematic_smiles::parse_cxsmiles(cxsmiles).map_err(|e| e.to_string())?;
    if cx.mol.atom_count() > MAX_SMILES_ATOMS {
        return Err(format!(
            "molecule exceeds maximum atom count ({MAX_SMILES_ATOMS})"
        ));
    }

    let atom_maps = cx
        .mol
        .atoms()
        .map(|(_, atom)| atom.atom_map)
        .collect::<Vec<_>>();
    let marked_attachment_points = (0..cx.mol.atom_count())
        .filter_map(|index| {
            let atom = chematic_core::AtomIdx(index as u32);
            cx.marked_attachment_point(atom).map(|identifier| {
                serde_json::json!({
                    "atom_index": index,
                    "atom_map": cx.mol.atom(atom).atom_map,
                    "identifier": identifier,
                })
            })
        })
        .collect::<Vec<_>>();
    Ok(serde_json::json!({
        "input_cxsmiles": cxsmiles,
        "cxsmiles": chematic_smiles::write_cxsmiles(&cx),
        "atom_count": cx.mol.atom_count(),
        "bond_count": cx.mol.bond_count(),
        "atom_maps": atom_maps,
        "atom_labels": cx.atom_labels,
        "marked_attachment_points": marked_attachment_points,
    })
    .to_string())
}

fn fingerprint_json(smiles: &str, algorithm: &str) -> Result<String, String> {
    let mol = parse_cli_smiles(smiles)?;
    let algorithm = algorithm.to_ascii_lowercase();
    let fp = fingerprint_for(&mol, &algorithm)?;
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
        "rdkit_morgan" => {
            chematic_fp::rdkit_morgan_ecfp4_bitvec(mol).map_err(|error| error.to_string())
        }
        "rdkit_pattern" => Ok(chematic_fp::rdkit_pattern_fp(mol)),
        "rdkit_rdk" => Ok(chematic_fp::rdkit_rdk_fp(mol)),
        "rdkit_torsion" => Ok(chematic_fp::rdkit_torsion_fp(mol)),
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

fn substructure_json(
    smiles: &str,
    smarts: &str,
    use_chirality: bool,
    rdkit_parity: bool,
) -> Result<String, String> {
    let parsed = parse_cli_smiles(smiles)?;
    // MolFromSmiles removes ordinary explicit H graph atoms before matching.
    // The parity matcher must see that same target graph: keeping the H atoms
    // changes degree queries and the indices in returned mappings.
    let mut mol = if rdkit_parity {
        chematic_smiles::rdkit_hydrogen_suppressed(&parsed).unwrap_or(parsed)
    } else {
        parsed
    };
    let has_stereo_markers = mol
        .atoms()
        .any(|(_, atom)| atom.chirality != chematic_core::Chirality::None)
        || mol.bonds().any(|(_, bond)| {
            matches!(
                bond.order,
                chematic_core::BondOrder::Up | chematic_core::BondOrder::Down
            )
        })
        || mol.has_bond_directions();
    let stereo = if use_chirality && has_stereo_markers {
        // MolFromSmiles discards a written @/@@ tag on a non-stereogenic atom.
        // Applying the same target-side perception avoids treating a symmetric
        // ring carbon as a chiral hit merely because its SMILES contained @.
        let stereo = chematic_smiles::rdkit_legacy_stereo(&mol).map_err(|e| e.to_string())?;
        for (index, tagged) in stereo.atom_tagged.iter().enumerate() {
            if !tagged {
                mol.set_chirality(
                    chematic_core::AtomIdx(index as u32),
                    chematic_core::Chirality::None,
                );
            }
        }
        Some(stereo)
    } else {
        None
    };
    let query = chematic_smarts::parse_smarts(smarts).map_err(|e| e.to_string())?;
    let config = chematic_smarts::MatchConfig {
        use_chirality,
        ..Default::default()
    };
    let mut matches = if rdkit_parity {
        let parity_config = chematic_smarts::RdkitParityConfig {
            base: config,
            use_rdkit_parity_aromaticity: true,
            ..Default::default()
        };
        let (matches, budget_exhausted) =
            chematic_smarts::find_matches_rdkit_parity(&query, &mol, &parity_config)
                .map_err(|e| e.to_string())?;
        if budget_exhausted {
            return Err("RDKit-parity matcher visit budget exhausted".to_string());
        }
        matches
    } else {
        chematic_smarts::find_matches_with_config(&query, &mol, &config)
    };
    if let Some(stereo) = &stereo {
        // A slash-written bond is not always stereogenic (e.g. identical
        // substituents at an alkene end). The matcher checks marker parity;
        // RDKit's target perception additionally rejects such double bonds.
        let stereo_doubles: Vec<_> = query
            .bonds
            .iter()
            .filter(|bond| {
                bond.query
                    == chematic_smarts::BondQuery::Primitive(chematic_smarts::BondPrimitive::Double)
                    && [bond.atom1, bond.atom2].into_iter().all(|end| {
                        query.adj[end].iter().any(|&(edge, other)| {
                            other != bond.atom1
                                && other != bond.atom2
                                && matches!(
                                    query.bonds[edge].query,
                                    chematic_smarts::BondQuery::Primitive(
                                        chematic_smarts::BondPrimitive::Up
                                            | chematic_smarts::BondPrimitive::Down
                                    )
                                )
                        })
                    })
            })
            .collect();
        matches.retain(|mapping| {
            stereo_doubles.iter().all(|bond| {
                let left = mapping[&bond.atom1];
                let right = mapping[&bond.atom2];
                mol.bond_between(left, right)
                    .is_some_and(|(idx, _)| stereo.bond_stereo[idx.0 as usize] > 1)
            })
        });
    }
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
        "use_chirality": use_chirality,
        "rdkit_parity": rdkit_parity,
        "match_count": matches.len(),
        "matches": matches,
    })
    .to_string())
}

fn standardize_json(smiles: &str, largest_fragment_only: bool) -> Result<String, String> {
    let mol = parse_cli_smiles(smiles)?;
    let input_smiles = chematic_smiles::canonical_smiles(&mol);
    let options = chematic_chem::StandardizeOptions {
        largest_fragment_only,
        ..Default::default()
    };
    let (standardized, report) = chematic_chem::StandardizationPipeline::new(options).run(&mol);
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
    let formed_bits = set_bits(&fingerprint.formed_fp);
    let broken_bits = set_bits(&fingerprint.broken_fp);
    let combined_bits = set_bits(&fingerprint.combined_fp);
    Ok(serde_json::json!({
        "reaction_smiles": chematic_rxn::write_reaction(&reaction),
        "mode": mode,
        "n_bits": 2048,
        "reactant_popcount": reactant_bits.len(),
        "product_popcount": product_bits.len(),
        "formed_popcount": formed_bits.len(),
        "broken_popcount": broken_bits.len(),
        "popcount": combined_bits.len(),
        "formed_bits": formed_bits,
        "broken_bits": broken_bits,
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

/// Canonicalize a bounded line-delimited SMILES batch without descriptor work.
///
/// Every accepted or rejected input remains in input order, so callers can use
/// it for identity and invariance checks without silently dropping an invalid
/// molecule or paying for unrelated descriptor calculations.
fn batch_canonicalize_json(text: &str, limits: &BatchLimits) -> Result<String, String> {
    let smiles = batch_lines(text, limits)?;
    let canonicalizer = chematic_smiles::SmilesBatchCanonicalizer::default();
    let result = canonicalizer.canonicalize_with_result(smiles.iter());
    let mut records = Vec::with_capacity(result.input_count);
    for record in result.records {
        match record.result {
            chematic_smiles::BatchCanonicalization::Accepted { canonical_smiles } => {
                records.push(serde_json::json!({
                    "input_index": record.input_index,
                    "input_smiles": record.input,
                    "status": "accepted",
                    "canonical_smiles": canonical_smiles,
                    "error": null,
                    "error_stage": null,
                }));
            }
            chematic_smiles::BatchCanonicalization::Rejected { error } => {
                records.push(serde_json::json!({
                    "input_index": record.input_index,
                    "input_smiles": record.input,
                    "status": "rejected",
                    "canonical_smiles": null,
                    "error": error,
                    "error_stage": "parse",
                }));
            }
        }
    }
    let mut output = serde_json::json!({
        "records": records,
        "valid_count": result.accepted_count,
        "error_count": result.rejected_count,
    });
    add_batch_manifest(&mut output, "canonicalize", limits);
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

fn batch_standardize_json(
    text: &str,
    largest_fragment_only: bool,
    limits: &BatchLimits,
) -> Result<String, String> {
    let smiles = batch_lines(text, limits)?;
    let mut records = Vec::with_capacity(smiles.len());
    let mut valid_count = 0usize;
    for (input_index, input_smiles) in smiles.iter().enumerate() {
        match standardize_json(input_smiles, largest_fragment_only).and_then(|json| {
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

fn batch_substructure_json(
    text: &str,
    limits: &BatchLimits,
    use_chirality: bool,
    rdkit_parity: bool,
) -> Result<String, String> {
    let lines = batch_lines(text, limits)?;
    let mut records = Vec::with_capacity(lines.len());
    let mut valid_count = 0usize;
    for (input_index, line) in lines.iter().enumerate() {
        let result = line
            .split_once('\t')
            .ok_or_else(|| "expected SMILES<TAB>SMARTS".to_string())
            .and_then(|(smiles, smarts)| {
                substructure_json(smiles.trim(), smarts.trim(), use_chirality, rdkit_parity)
            })
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

fn main() {
    if let Err(error) = dispatch::run(Cli::parse()) {
        eprintln!("chematic: {error}");
        std::process::exit(2);
    }
}

#[cfg(test)]
mod tests;
