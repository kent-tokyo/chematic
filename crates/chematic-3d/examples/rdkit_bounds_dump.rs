//! Dump RDKit-port bounds matrices of `add_hydrogens(parse(smiles))`:
//! one JSON line per SMILES with the raw (`setTopolBounds`) matrix for
//! useMacrocycle14config false/true and the triangle-smoothed ones, floats
//! as shortest round-trip strings.
use chematic_3d::rdkit_embed::{
    BoundsMatrix, init_bounds_mat, set_topol_bounds, triangle_smooth_bounds,
};
use std::io::BufRead;

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args.next().expect("smiles file");
    let n: usize = args
        .next()
        .map(|s| s.parse().unwrap())
        .unwrap_or(usize::MAX);
    for line in std::io::BufReader::new(std::fs::File::open(path).unwrap())
        .lines()
        .take(n)
    {
        let line = line.unwrap();
        let smi = line.split_whitespace().next().unwrap_or("").to_string();
        let Ok(m) = chematic_smiles::parse(&smi) else {
            println!("{}", serde_json::json!({"smiles": smi, "error": "parse"}));
            continue;
        };
        let mh = chematic_chem::add_hydrogens(&m);
        let view = match chematic_smiles::rdkit_mol_view(&mh) {
            Ok(v) => v,
            Err(e) => {
                println!(
                    "{}",
                    serde_json::json!({"smiles": smi, "error": e.to_string()})
                );
                continue;
            }
        };
        let labels = chematic_ff::rdkit_uff::rdkit_uff_atom_labels(&mh);
        let mut out = serde_json::json!({"smiles": smi});
        for mac in [false, true] {
            let mut bm = BoundsMatrix::new(view.num_atoms());
            init_bounds_mat(&mut bm);
            let key = if mac { "mac" } else { "std" };
            match set_topol_bounds(&view, &labels, &mut bm, true, mac, true) {
                Ok(()) => {
                    out[format!("{key}_raw")] = bm
                        .raw()
                        .iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<_>>()
                        .into();
                    let ok = triangle_smooth_bounds(&mut bm, 0.0);
                    out[format!("{key}_smooth_ok")] = ok.into();
                    out[format!("{key}_smooth")] = bm
                        .raw()
                        .iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<_>>()
                        .into();
                }
                Err(e) => out[format!("{key}_error")] = format!("{e:?}").into(),
            }
        }
        println!("{out}");
    }
}
