//! Source-level speed gate. Build once before edits and retain that executable.
//! Run baseline/current processes alternately with identical arguments:
//! `hotpath_throughput corpus.smi corpus.sdf [repeats=30]`.
use std::{fs, hint::black_box, io::BufReader, time::Instant};

use chematic_mol::{SdfFileReader, mol2000::write_sdf_record_into};
use chematic_smiles::{canonical_smiles, parse};

fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert!(args.len() >= 3, "pass SMILES and SDF paths");
    let repeats: usize = args.get(3).map_or(30, |s| s.parse().unwrap());
    assert!(repeats > 0);
    let input = fs::read_to_string(&args[1]).unwrap();
    let smiles: Vec<_> = input
        .lines()
        .filter_map(|l| l.split_whitespace().next())
        .filter(|l| !l.starts_with('#'))
        .collect();
    let mols: Vec<_> = smiles.iter().map(|s| parse(s).unwrap()).collect();
    let read = || {
        SdfFileReader::fast(BufReader::new(fs::File::open(&args[2]).unwrap()))
            .map(Result::unwrap)
            .collect::<Vec<_>>()
    };
    let records = read();
    assert!(!smiles.is_empty() && !records.is_empty());
    // Output validation happens outside timing; FNV is an accidental-change
    // sentinel. Use --dump for an exact, length-delimited byte comparison.
    let mut digest = 0xcbf29ce484222325_u64;
    let mut output = Vec::new();
    let mut add = |s: &str| {
        for byte in (s.len() as u64).to_le_bytes().into_iter().chain(s.bytes()) {
            digest = (digest ^ u64::from(byte)).wrapping_mul(0x100000001b3);
            if args.iter().any(|s| s == "--dump") {
                output.push(byte);
            }
        }
    };
    for mol in &mols {
        add(&canonical_smiles(mol));
    }
    let mut buffer = String::new();
    for rec in &records {
        // HashMap property iteration is intentionally not a serialization
        // ordering contract: compare graph serialization and sorted properties.
        write_sdf_record_into(&mut buffer, &rec.mol, &rec.meta, &[], &Default::default());
        add(&buffer);
        let mut props: Vec<_> = rec.properties.iter().collect();
        props.sort_unstable();
        for (k, v) in props {
            add(k);
            add(v);
        }
    }
    if args.iter().any(|s| s == "--dump") {
        use std::io::Write;
        std::io::stdout().write_all(&output).unwrap();
        return;
    }
    let measure = |mut run: Box<dyn FnMut() + '_>, count: usize| {
        run();
        let mut times = Vec::new();
        for _ in 0..5 {
            let start = Instant::now();
            for _ in 0..repeats {
                run();
            }
            times.push(start.elapsed().as_secs_f64() * 1e6 / (count * repeats) as f64);
        }
        times.sort_by(f64::total_cmp);
        times[2]
    };
    let parse_us = measure(
        Box::new(|| {
            for s in &smiles {
                black_box(parse(black_box(s)).unwrap());
            }
        }),
        smiles.len(),
    );
    let canonical_us = measure(
        Box::new(|| {
            for mol in &mols {
                black_box(canonical_smiles(black_box(mol)));
            }
        }),
        mols.len(),
    );
    let read_us = measure(
        Box::new(|| {
            black_box(read());
        }),
        records.len(),
    );
    let write_us = measure(
        Box::new(|| {
            for rec in &records {
                write_sdf_record_into(&mut buffer, &rec.mol, &rec.meta, &[], &rec.properties);
                black_box(&buffer);
            }
        }),
        records.len(),
    );
    println!(
        "{}",
        serde_json::json!({"smiles": smiles.len(), "records": records.len(), "repeats": repeats,
        "output_fnv1a": format!("{digest:016x}"), "parse_us": parse_us,
        "canonical_us": canonical_us, "sdf_read_us": read_us, "sdf_write_us": write_us})
    );
}
