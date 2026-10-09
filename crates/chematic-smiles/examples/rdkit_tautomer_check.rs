//! Compares `rdkit_canonical_tautomer`, `rdkit_enumerate_tautomers` and
//! `rdkit_tautomer_score` with an RDKit reference dump (one JSON object per
//! line: `smi`, `canon`, `taut`, `score`, `status`, or `err`).
//!
//! `cargo run --release --example rdkit_tautomer_check -- REF.jsonl [LIMIT] [-v]`

use std::io::BufRead;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = &args[1];
    let limit: usize = args
        .get(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(usize::MAX);
    let verbose = args.iter().any(|a| a == "-v");
    let only: Option<usize> = args
        .iter()
        .position(|a| a == "--only")
        .and_then(|i| args.get(i + 1))
        .and_then(|s| s.parse().ok());
    let f = std::io::BufReader::new(std::fs::File::open(path).expect("open"));
    let (mut n, mut ok_c, mut ok_t, mut ok_s, mut ok_st, mut skipped) = (0, 0, 0, 0, 0, 0);
    let mut errs = 0;
    for line in f.lines().take(limit) {
        let line = line.expect("line");
        let v: serde_json::Value = serde_json::from_str(&line).expect("json");
        let i = v["i"].as_u64().unwrap() as usize;
        if only.is_some_and(|o| o != i) {
            continue;
        }
        if v.get("err").is_some() {
            skipped += 1;
            continue;
        }
        n += 1;
        let smi = v["smi"].as_str().unwrap();
        let mol = match chematic_smiles::parse(smi) {
            Ok(m) => m,
            Err(e) => {
                errs += 1;
                println!("{i}\tPARSE\t{smi}\t{e}");
                continue;
            }
        };
        let t0 = std::time::Instant::now();
        let canon = chematic_smiles::rdkit_canonical_tautomer(&mol);
        let en = chematic_smiles::rdkit_enumerate_tautomers(&mol);
        let score = chematic_smiles::rdkit_tautomer_score(&mol);
        let dt = t0.elapsed().as_secs_f64();
        if dt > 2.0 {
            eprintln!("{i}\tSLOW {dt:.1}s\t{smi}");
        }
        let rc = v["canon"].as_str().unwrap();
        match &canon {
            Ok(c) if c == rc => ok_c += 1,
            Ok(c) => println!("{i}\tCANON\t{smi}\tgot {c}\twant {rc}"),
            Err(e) => {
                errs += 1;
                println!("{i}\tCANON-ERR\t{smi}\t{e}")
            }
        }
        let rt: Vec<String> = v["taut"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_str().unwrap().to_string())
            .collect();
        let rstatus = v["status"].as_u64().unwrap();
        match &en {
            Ok(e) => {
                let mut got = e.smiles.clone();
                got.sort();
                if got == rt {
                    ok_t += 1;
                } else {
                    let extra: Vec<_> = got.iter().filter(|x| !rt.contains(x)).collect();
                    let missing: Vec<_> = rt.iter().filter(|x| !got.contains(x)).collect();
                    println!(
                        "{i}\tTAUT\t{smi}\tgot {} want {}\textra {:?}\tmissing {:?}",
                        got.len(),
                        rt.len(),
                        if verbose {
                            extra.clone()
                        } else {
                            extra.into_iter().take(3).collect()
                        },
                        if verbose {
                            missing.clone()
                        } else {
                            missing.into_iter().take(3).collect()
                        }
                    );
                }
                if e.status as u64 == rstatus {
                    ok_st += 1;
                } else {
                    println!("{i}\tSTATUS\t{smi}\tgot {:?} want {rstatus}", e.status);
                }
            }
            Err(e) => println!("{i}\tTAUT-ERR\t{smi}\t{e}"),
        }
        let rs = v["score"].as_i64().unwrap();
        match &score {
            Ok(s) if i64::from(*s) == rs => ok_s += 1,
            Ok(s) => println!("{i}\tSCORE\t{smi}\tgot {s} want {rs}"),
            Err(e) => println!("{i}\tSCORE-ERR\t{smi}\t{e}"),
        }
    }
    println!(
        "SUMMARY rows={n} skipped={skipped} canon={ok_c} taut={ok_t} status={ok_st} score={ok_s} errs={errs}"
    );
}
