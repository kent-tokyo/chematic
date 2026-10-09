//! Dumps RDKit-port SMILES surfaces for a SMILES file, one JSON object per
//! line: `MolToCXSmiles`, `MolToRandomSmilesVect(5, 42)` and
//! `MolFragmentToSmiles` of the first half of the atoms.
//!
//!     cargo run --release -p chematic-smiles --example rdkit_surfaces_dump -- FILE.smi

use chematic_smiles::{
    RdkitSmilesParams, parse, rdkit_cx_smiles, rdkit_fragment_smiles, rdkit_random_smiles,
};
use serde_json::{Value, json};

fn wrap<T: Into<Value>>(r: Result<T, chematic_smiles::RdkitSmilesError>) -> Value {
    match r {
        Ok(v) => json!({"ok": v.into()}),
        Err(e) => json!({"err": e.to_string()}),
    }
}

fn main() {
    let path = std::env::args().nth(1).expect("SMILES file");
    let text = std::fs::read_to_string(path).expect("read");
    let p = RdkitSmilesParams::default();
    let noncanon = RdkitSmilesParams {
        canonical: false,
        ..p
    };
    let kek = RdkitSmilesParams {
        kekule: true,
        all_hs_explicit: true,
        ..p
    };
    for line in text.lines() {
        let smi = line.split_whitespace().next().unwrap_or("");
        if smi.is_empty() {
            continue;
        }
        let Ok(mol) = parse(smi) else {
            println!("{}", json!({"smiles": smi, "parse": false}));
            continue;
        };
        let n = chematic_smiles::rdkit_num_atoms(&mol).unwrap_or(0);
        let atoms: Vec<usize> = (0..(n / 2).max(1)).collect();
        println!(
            "{}",
            json!({
                "smiles": smi,
                "cx": wrap(rdkit_cx_smiles(&mol, &p)),
                "rand": wrap(rdkit_random_smiles(&mol, 5, 42, &p)),
                "frag": wrap(rdkit_fragment_smiles(&mol, &atoms, None, &p)),
                "frag2": wrap(rdkit_fragment_smiles(&mol, &(0..n).step_by(2).collect::<Vec<_>>(), None, &p)),
                "frag3": wrap(rdkit_fragment_smiles(&mol, &(n / 3..n).collect::<Vec<_>>(), None, &noncanon)),
                "randk": wrap(rdkit_random_smiles(&mol, 3, 7, &kek)),
            })
        );
    }
}
