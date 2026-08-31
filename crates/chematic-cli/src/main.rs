//! Small, composable CLI for the common topology-bearing format bridge.

#![forbid(unsafe_code)]

use clap::{Parser, Subcommand};
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

fn read_input(path: Option<&PathBuf>) -> Result<String, String> {
    match path {
        Some(path) => fs::read_to_string(path).map_err(|e| format!("read {}: {e}", path.display())),
        None => {
            let mut text = String::new();
            io::stdin()
                .read_to_string(&mut text)
                .map_err(|e| format!("read stdin: {e}"))?;
            Ok(text)
        }
    }
}

fn write_output(path: Option<&PathBuf>, text: &str) -> Result<(), String> {
    match path {
        Some(path) => fs::write(path, text).map_err(|e| format!("write {}: {e}", path.display())),
        None => io::stdout()
            .write_all(text.as_bytes())
            .map_err(|e| format!("write stdout: {e}")),
    }
}

fn descriptors_json(smiles: &str) -> Result<String, String> {
    let mol = chematic_smiles::parse(smiles).map_err(|e| e.to_string())?;
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

fn fingerprint_json(smiles: &str, algorithm: &str) -> Result<String, String> {
    let mol = chematic_smiles::parse(smiles).map_err(|e| e.to_string())?;
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

fn run(cli: Cli) -> Result<(), String> {
    match cli.command {
        Command::Convert {
            input_format,
            output_format,
            input,
            output,
        } => {
            let text = read_input(input.as_ref())?;
            let converted = convert_text(&text, &input_format, &output_format)?;
            write_output(output.as_ref(), &converted)
        }
        Command::Descriptors { smiles } => {
            let json = descriptors_json(&smiles)?;
            write_output(None, &format!("{json}\n"))
        }
        Command::Fingerprint { smiles, algorithm } => {
            let json = fingerprint_json(&smiles, &algorithm)?;
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
    use super::{convert_text, descriptors_json, fingerprint_json};

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
}
