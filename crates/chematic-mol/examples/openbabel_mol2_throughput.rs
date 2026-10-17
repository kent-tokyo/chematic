//! Native MOL2 corpus timing and full indexed source snapshots.
use chematic_mol::{Mol2Record, parse_mol2_record, write_mol2_record};
use serde_json::json;
use std::{error::Error, fs, hint::black_box, io::Write};
// This probe uses untimed snapshots, leaving the shared dump warmup helpers unused.
#[allow(dead_code)]
#[path = "support/steady_timing.rs"]
mod steady_timing;

fn read(input: &str, expected: usize) -> Result<Vec<Mol2Record>, Box<dyn Error>> {
    let records = input
        .split("\n$$$$\n")
        .map(parse_mol2_record)
        .collect::<Result<Vec<_>, _>>()?;
    if records.len() != expected {
        return Err("record count differs".into());
    }
    Ok(records)
}
fn write(records: &[Mol2Record]) -> usize {
    records
        .iter()
        .map(|r| black_box(write_mol2_record(r)).len())
        .sum()
}
fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 5 {
        return Err("usage: PATH MODE REPEATS EXPECTED_RECORDS".into());
    }
    let input = fs::read_to_string(&args[1])?;
    let mode = &args[2];
    let repeats: usize = args[3].parse()?;
    let expected: usize = args[4].parse()?;
    if repeats == 0 || expected == 0 {
        return Err("positive counts required".into());
    }
    let prepared = read(&input, expected)?;
    if mode == "snapshot" {
        let result: Vec<_> = prepared.iter().map(|r| json!({
            "atoms":r.molecule.atoms().map(|(i,a)| format!("{i:?}:{a:?}:{:?}",r.molecule.stereo_neighbor_order(i))).collect::<Vec<_>>(),
            "bonds":r.molecule.bonds().map(|(i,b)| format!("{i:?}:{b:?}")).collect::<Vec<_>>(),
            "coords_bits":r.coords.iter().map(|(x,y,z)|[x.to_bits(),y.to_bits(),z.to_bits()]).collect::<Vec<_>>(),
            "name":r.name,"molecule_type":r.molecule_type,"charge_type":r.charge_type,
            "molecule_tail":r.molecule_tail,"source_atoms":format!("{:?}",r.atoms),
            "source_bonds":format!("{:?}",r.bonds),"duplicate_bonds":format!("{:?}",r.duplicate_bonds),
            "unity":format!("{:?}",r.unity_atom_attributes),"opaque":format!("{:?}",r.opaque_sections)
        })).collect();
        println!("{}", serde_json::to_string(&result)?);
        return Ok(());
    }
    if mode == "dump" {
        let mut out = std::io::stdout().lock();
        for (i, r) in prepared.iter().enumerate() {
            if i > 0 {
                out.write_all(b"\n$$$$\n")?;
            }
            out.write_all(write_mol2_record(r).as_bytes())?;
        }
        return Ok(());
    }
    black_box(write(&prepared));
    let timing = steady_timing::measure(repeats, || {
        Ok(match mode.as_str() {
            "parse" => black_box(read(black_box(&input), expected)?)
                .iter()
                .map(|r| r.molecule.atom_count())
                .sum(),
            "write" => write(&prepared),
            "roundtrip" => write(&read(black_box(&input), expected)?),
            _ => return Err("unknown operation".into()),
        })
    })?;
    println!(
        "{}",
        json!({"engine":"chematic","mode":mode,"records":expected*timing.repeats,
        "repeats":timing.repeats,"output_units":timing.units,"elapsed_ns":timing.elapsed_ns,
        "warmup_repeats":timing.warmup_repeats,"warmup_elapsed_ns":timing.warmup_elapsed_ns,
        "minimum_warmup_ns":timing.minimum_warmup_ns,"minimum_timed_ns":timing.minimum_timed_ns})
    );
    Ok(())
}
