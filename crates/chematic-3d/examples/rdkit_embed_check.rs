//! Compare `rdkit_embed_molecule` with RDKit reference coordinates.
//! Input JSONL lines: {"smiles", "mode": "dg"|"etkdg", "coords": [str] | null}
//! (from rdkit_embed_gen.py). Prints per-row status and a summary.
use chematic_3d::rdkit_embed::{RdkitEmbedError, RdkitEmbedOptions, rdkit_embed_molecule};
use std::io::BufRead;

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args.next().expect("input");
    let verbose = args.next().is_some();
    let mut counts = std::collections::BTreeMap::<String, usize>::new();
    for line in std::io::BufReader::new(std::fs::File::open(path).unwrap()).lines() {
        let line = line.unwrap();
        let v: serde_json::Value = serde_json::from_str(&line).unwrap();
        let smi = v["smiles"].as_str().unwrap();
        let mut o = RdkitEmbedOptions::default();
        if v["mode"] == "dg" {
            o.use_exp_torsion_angle_prefs = false;
            o.use_basic_knowledge = false;
        }
        let expect: Option<Vec<f64>> = v["coords"].as_array().map(|a| {
            a.iter()
                .map(|x| x.as_str().unwrap().parse::<f64>().unwrap())
                .collect()
        });
        let mol = chematic_smiles::parse(smi).unwrap();
        let mh = chematic_chem::add_hydrogens(&mol);
        let got = rdkit_embed_molecule(&mh, &o);
        let status = match (&got, &expect) {
            (Ok(g), Some(e)) => {
                let md = g
                    .iter()
                    .flatten()
                    .zip(e)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0, f64::max);
                if g.len() * 3 == e.len() && md <= 1e-8 {
                    "match".to_string()
                } else {
                    format!("mismatch {md:e}")
                }
            }
            (Err(RdkitEmbedError::Failed), None) => "match".to_string(),
            (Err(RdkitEmbedError::Unsupported(w)), _) => format!("unsupported: {w}"),
            (Err(e), _) => format!("error: {e}"),
            (Ok(_), None) => "rdkit failed, we succeeded".to_string(),
        };
        if verbose && status != "match" {
            eprintln!("{smi}\t{status}");
        }
        let key = if status.starts_with("mismatch") {
            "mismatch".to_string()
        } else {
            status
        };
        *counts.entry(key).or_default() += 1;
    }
    for (k, c) in counts {
        println!("{c}\t{k}");
    }
}
