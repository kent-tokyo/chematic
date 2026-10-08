//! Compare `rdkit_2d_coords` with RDKit's `rdDepictor.Compute2DCoords` for
//! a TSV of `smiles<TAB>n_atoms<TAB>x,y;x,y;...` lines (coordinates as
//! `float.hex()` strings, so they are read back exactly).
//!
//! `cargo run --release -p chematic-smiles --example rdkit_depict_check -- file.ref [max_shown] [tol]`

use std::io::BufRead;

fn parse_hex(s: &str) -> f64 {
    let s = s.trim();
    let (neg, s) = match s.strip_prefix('-') {
        Some(r) => (true, r),
        None => (false, s),
    };
    let s = s.strip_prefix("0x").expect("hex float");
    let (mant, exp) = s.split_once('p').expect("exponent");
    let exp: i32 = exp.parse().expect("exp");
    let (ip, fp) = mant.split_once('.').unwrap_or((mant, ""));
    let mut v: u64 = u64::from_str_radix(ip, 16).expect("int part");
    let mut e = exp;
    for c in fp.chars() {
        v = v * 16 + c.to_digit(16).expect("hex digit") as u64;
        e -= 4;
    }
    let r = (v as f64) * 2f64.powi(e);
    if neg { -r } else { r }
}

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args
        .next()
        .expect("usage: rdkit_depict_check FILE [max_shown] [tol]");
    let max_shown: usize = args.next().map_or(10, |s| s.parse().expect("number"));
    let tol: f64 = args.next().map_or(1e-6, |s| s.parse().expect("number"));
    let file = std::fs::File::open(&path).expect("open");
    let (mut exact, mut close, mut bad, mut err, mut skipped) =
        (0usize, 0usize, 0usize, 0usize, 0usize);
    let mut shown = 0;
    for (lineno, line) in std::io::BufReader::new(file).lines().enumerate() {
        let line = line.expect("read");
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            skipped += 1;
            continue;
        }
        let smi = parts[0];
        let n: usize = parts[1].parse().expect("n");
        let expected: Vec<[f64; 2]> = parts[2]
            .split(';')
            .filter(|s| !s.is_empty())
            .map(|p| {
                let (x, y) = p.split_once(',').expect("pair");
                [parse_hex(x), parse_hex(y)]
            })
            .collect();
        assert_eq!(expected.len(), n);
        let got = chematic_smiles::parse(smi)
            .map_err(|e| e.to_string())
            .and_then(|m| chematic_smiles::rdkit_2d_coords(&m).map_err(|e| e.to_string()));
        match got {
            Err(e) => {
                err += 1;
                if shown < max_shown {
                    shown += 1;
                    println!("ERROR #{lineno} {smi}: {e}");
                }
            }
            Ok(g) if g.len() != n => {
                err += 1;
                if shown < max_shown {
                    shown += 1;
                    println!("NATOMS #{lineno} {smi}: {} vs {n}", g.len());
                }
            }
            Ok(g) => {
                let mut maxd = 0.0f64;
                let mut worst = 0;
                let mut all_exact = true;
                for i in 0..n {
                    for k in 0..2 {
                        if g[i][k].to_bits() != expected[i][k].to_bits() {
                            all_exact = false;
                        }
                        let d = (g[i][k] - expected[i][k]).abs();
                        if d > maxd {
                            maxd = d;
                            worst = i;
                        }
                    }
                }
                if all_exact {
                    exact += 1;
                } else if maxd <= tol {
                    close += 1;
                } else {
                    bad += 1;
                    if shown < max_shown {
                        shown += 1;
                        println!(
                            "MISMATCH #{lineno} {smi}: max|d|={maxd:.3e} at atom {worst}: ours {:?} rdkit {:?}",
                            g[worst], expected[worst]
                        );
                    }
                }
            }
        }
    }
    println!(
        "exact {exact} close(<= {tol:e}) {close} mismatch {bad} error {err} skipped {skipped}"
    );
}
