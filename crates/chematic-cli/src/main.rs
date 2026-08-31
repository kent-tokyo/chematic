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
    use super::convert_text;

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
}
