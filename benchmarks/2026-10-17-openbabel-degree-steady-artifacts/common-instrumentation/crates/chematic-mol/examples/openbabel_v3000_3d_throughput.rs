//! Native 3D V3000 MOL APIs; fail closed if a parsed conformer is absent.
use std::{error::Error, fs, hint::black_box, io::Write};

use chematic_mol::{
    MolReadReport, read_mol_v3000_with_diagnostics, write_mol_v3000_with_conformer,
};
use serde_json::json;

#[path = "support/steady_timing.rs"]
mod steady_timing;

fn blocks(input: &str, expected: usize) -> Result<Vec<&str>, Box<dyn Error>> {
    let (mut start, mut position) = (0, 0);
    let mut result = Vec::new();
    for line in input.split_inclusive('\n') {
        if line.trim_end_matches(['\r', '\n']) == "$$$$" {
            result.push(&input[start..position]);
            start = position + line.len();
        }
        position += line.len();
    }
    if !input[start..].trim().is_empty() || result.len() != expected {
        return Err("V3000 framed record count differs".into());
    }
    Ok(result)
}

fn read(blocks: &[&str]) -> Result<Vec<MolReadReport>, Box<dyn Error>> {
    blocks
        .iter()
        .map(|block| read_mol_v3000_with_diagnostics(block).map_err(Into::into))
        .collect()
}

fn write(records: &[MolReadReport]) -> usize {
    records
        .iter()
        .map(|r| {
            black_box(write_mol_v3000_with_conformer(
                &r.mol,
                &r.metadata,
                r.conformer.as_ref().expect("3D corpus requires conformer"),
            ))
            .len()
        })
        .sum()
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 5 {
        return Err("usage: PATH MODE REPEATS EXPECTED_RECORDS".into());
    }
    let input = fs::read_to_string(&args[1])?;
    let repeats: usize = args[3].parse()?;
    let expected: usize = args[4].parse()?;
    if repeats == 0 || expected == 0 {
        return Err("positive repetitions and records required".into());
    }
    let blocks = blocks(&input, expected)?;
    let prepared = read(&blocks)?;
    if args[2] == "snapshot" {
        let rows: Vec<_> = prepared.iter().map(|r| json!({
            "canonical_smiles": chematic_smiles::rdkit_canonical_smiles(&r.mol).expect("supported corpus graph"),
            "atoms": r.mol.atoms().map(|(i,a)| format!("{i:?}:{a:?}:{:?}", r.mol.stereo_neighbor_order(i))).collect::<Vec<_>>(),
            "bonds": r.mol.bonds().map(|(i,b)| format!("{i:?}:{b:?}")).collect::<Vec<_>>(),
            "coords_bits": r.coords.iter().map(|(x,y)| [x.to_bits(),y.to_bits()]).collect::<Vec<_>>(),
            "metadata": format!("{:?}",r.metadata),
            "diagnostics": format!("{:?}:{:?}:{:?}:{:?}:{:?}:{:?}", r.stereo_diagnostics,r.ez_diagnostics,
                r.coordinate_dimension,r.geometry_rank,r.stereo3d_diagnostics,r.conformer),
        })).collect();
        println!("{}", serde_json::to_string(&rows)?);
        return Ok(());
    }
    write(&prepared);
    if args[2] == "dump" {
        if steady_timing::warmed_dump() {
            steady_timing::measure(steady_timing::dump_repeats()?, || Ok(write(&prepared)))?;
        }
        let mut stdout = std::io::stdout().lock();
        for r in &prepared {
            stdout.write_all(
                write_mol_v3000_with_conformer(
                    &r.mol,
                    &r.metadata,
                    r.conformer.as_ref().expect("3D corpus requires conformer"),
                )
                .as_bytes(),
            )?;
            stdout.write_all(b"$$$$\n")?;
        }
        return Ok(());
    }
    let timing = steady_timing::measure(repeats, || {
        Ok(match args[2].as_str() {
            "parse" => {
                let records = black_box(read(black_box(&blocks))?);
                records.iter().map(|r| r.mol.atom_count()).sum::<usize>()
            }
            "write" => write(&prepared),
            "roundtrip" => write(&read(black_box(&blocks))?),
            _ => return Err("unknown operation".into()),
        })
    })?;
    println!(
        "{}",
        json!({"engine":"chematic","mode":args[2],"records":expected*timing.repeats,
        "repeats":timing.repeats,"output_units":timing.units,"elapsed_ns":timing.elapsed_ns,
        "warmup_repeats":timing.warmup_repeats,"warmup_elapsed_ns":timing.warmup_elapsed_ns,
        "minimum_warmup_ns":timing.minimum_warmup_ns,"minimum_timed_ns":timing.minimum_timed_ns})
    );
    Ok(())
}
