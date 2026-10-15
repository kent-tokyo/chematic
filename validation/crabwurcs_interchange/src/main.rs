use std::env;
use std::error::Error;
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

use chematic_mol::{
    MolMetadata, read_mol_v3000_with_diagnostics, read_mol_with_diagnostics,
    write_mol_v3000_with_stereo_report, write_mol_with_stereo_report,
};
use chematic_smiles::{canonical_smiles, parse};

const TABLES: [&str; 2] = [
    "glycoshape_notations.tsv",
    "glycoshape_derived_notations.tsv",
];

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = env::args().skip(1);
    let data_dir = args
        .next()
        .ok_or("usage: probe <crabwurcs-data-dir> <output-dir>")?;
    let output_dir = args
        .next()
        .ok_or("usage: probe <crabwurcs-data-dir> <output-dir>")?;
    if args.next().is_some() {
        return Err("usage: probe <crabwurcs-data-dir> <output-dir>".into());
    }
    fs::create_dir_all(&output_dir)?;
    let mut outputs = [
        File::create(Path::new(&output_dir).join("v2000.sdf"))?,
        File::create(Path::new(&output_dir).join("v3000.sdf"))?,
    ];
    let mut rows = 0;
    let mut failures = Vec::new();

    for table in TABLES {
        let source = fs::read_to_string(Path::new(&data_dir).join(table))?;
        for (line_index, line) in source.lines().enumerate() {
            let label = format!("{table}:{}", line_index + 1);
            rows += 1;
            let Some(smiles) = line.split('\t').nth(4) else {
                failures.push(format!("{label}: missing fifth SMILES field"));
                continue;
            };
            let mol = match parse(smiles) {
                Ok(mol) => mol,
                Err(error) => {
                    failures.push(format!("{label}: SMILES parse: {error}"));
                    continue;
                }
            };
            let before = canonical_smiles(&mol);
            let metadata = MolMetadata::default();
            for (format, output) in outputs.iter_mut().enumerate() {
                let format_name = if format == 0 { "V2000" } else { "V3000" };
                let (block, loss) = if format == 0 {
                    write_mol_with_stereo_report(&mol, &metadata, &[])
                } else {
                    write_mol_v3000_with_stereo_report(&mol, &metadata, &[])
                };
                output.write_all(block.as_bytes())?;
                output.write_all(b"$$$$\n")?;
                if !loss.is_empty() {
                    failures.push(format!("{label} {format_name}: stereo loss: {loss}"));
                }
                let read = if format == 0 {
                    read_mol_with_diagnostics(&block)
                } else {
                    read_mol_v3000_with_diagnostics(&block)
                };
                match read {
                    Err(error) => {
                        failures.push(format!("{label} {format_name}: MOL read: {error}"))
                    }
                    Ok(report) => {
                        if !report.stereo_diagnostics.is_empty() {
                            failures.push(format!(
                                "{label} {format_name}: stereo diagnostics: {:?}",
                                report.stereo_diagnostics
                            ));
                        }
                        let after = canonical_smiles(&report.mol);
                        if before != after {
                            failures.push(format!(
                                "{label} {format_name}: canonical SMILES differ: {before} != {after}"
                            ));
                        }
                    }
                }
            }
        }
    }
    println!("source_rows={rows} failures={}", failures.len());
    for failure in failures.iter().take(20) {
        eprintln!("{failure}");
    }
    if failures.is_empty() {
        Ok(())
    } else {
        Err(format!("{} interchange checks failed", failures.len()).into())
    }
}
