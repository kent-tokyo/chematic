//! In-process SDF read/write throughput and exact source-regression snapshots.
use std::{collections::BTreeMap, error::Error, fs, hint::black_box, io::Write};

use chematic_mol::{SdfFileReader, SdfRecord, mol2000::write_sdf_record_into};
use serde_json::json;

#[path = "support/steady_timing.rs"]
mod steady_timing;

fn read(input: &str, expected: usize) -> Result<Vec<SdfRecord>, Box<dyn Error>> {
    let records = SdfFileReader::new(std::io::Cursor::new(input.as_bytes()))
        .collect::<Result<Vec<_>, _>>()?;
    if records.len() != expected {
        return Err(format!("expected {expected} records, got {}", records.len()).into());
    }
    Ok(records)
}

fn write(records: &[SdfRecord], output: &mut String) -> usize {
    let mut bytes = 0;
    for record in records {
        write_sdf_record_into(
            output,
            &record.mol,
            &record.meta,
            &record.coords,
            &record.properties,
        );
        bytes += black_box(output.len());
    }
    bytes
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 5 {
        return Err("usage: openbabel_sdf_throughput PATH MODE REPEATS EXPECTED_RECORDS".into());
    }
    let input = fs::read_to_string(&args[1])?;
    let mode = &args[2];
    let repeats: usize = args[3].parse()?;
    let expected: usize = args[4].parse()?;
    if repeats == 0 || expected == 0 {
        return Err("positive repetitions and record count required".into());
    }
    let prepared = read(&input, expected)?;
    let mut buffer = String::new();
    if mode == "snapshot" {
        let snapshot: Vec<_> = prepared.iter().map(|r| json!({
            "canonical_smiles": chematic_smiles::rdkit_canonical_smiles(&r.mol).expect("supported corpus molecule"),
            "atoms": r.mol.atoms().map(|(i,a)| format!("{i:?}:{a:?}:{:?}", r.mol.stereo_neighbor_order(i))).collect::<Vec<_>>(),
            "bonds": r.mol.bonds().map(|(i,b)| format!("{i:?}:{b:?}")).collect::<Vec<_>>(),
            "coords_bits": r.coords.iter().map(|(x,y)| [x.to_bits(), y.to_bits()]).collect::<Vec<_>>(),
            "metadata": format!("{:?}", r.meta),
            "properties": r.properties.iter().collect::<BTreeMap<_,_>>(),
            "stereo": format!("{:?}:{:?}:{:?}:{:?}:{:?}", r.stereo_diagnostics, r.ez_diagnostics,
                r.coordinate_dimension, r.geometry_rank, r.stereo3d_diagnostics),
        })).collect();
        println!("{}", serde_json::to_string(&snapshot)?);
        return Ok(());
    }
    if mode == "dump" {
        if steady_timing::warmed_dump() {
            steady_timing::measure(steady_timing::dump_repeats()?, || {
                Ok(write(&prepared, &mut buffer))
            })?;
        }
        let mut stdout = std::io::stdout().lock();
        for r in &prepared {
            write_sdf_record_into(&mut buffer, &r.mol, &r.meta, &r.coords, &r.properties);
            stdout.write_all(buffer.as_bytes())?;
        }
        return Ok(());
    }
    // Warm disjoint results; prepared write inputs stay unchanged.
    write(&prepared, &mut buffer);
    let timing = steady_timing::measure(repeats, || {
        Ok(match mode.as_str() {
            "parse" => {
                let records = black_box(read(black_box(&input), expected)?);
                records.iter().map(|r| r.mol.atom_count()).sum::<usize>()
            }
            "write" => write(&prepared, &mut buffer),
            "roundtrip" => write(&read(black_box(&input), expected)?, &mut buffer),
            _ => return Err("unknown operation".into()),
        })
    })?;
    println!(
        "{}",
        json!({"engine":"chematic", "mode":mode, "records":expected * timing.repeats,
        "repeats":timing.repeats, "output_units":timing.units, "elapsed_ns":timing.elapsed_ns,
        "warmup_repeats":timing.warmup_repeats,"warmup_elapsed_ns":timing.warmup_elapsed_ns,
        "minimum_warmup_ns":timing.minimum_warmup_ns,"minimum_timed_ns":timing.minimum_timed_ns})
    );
    Ok(())
}
