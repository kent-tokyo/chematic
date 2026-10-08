//! Dump `chematic_smiles::rdkit_mol_view` of `add_hydrogens(parse(smiles))`
//! for every SMILES of a file (first N), one JSON line each, in the format
//! of `rdkit_view_gen.py`, to diff RDKit's molecule model.
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
        let out = match chematic_smiles::parse(&smi) {
            Ok(m) => {
                let mh = chematic_chem::add_hydrogens(&m);
                match chematic_smiles::rdkit_mol_view(&mh) {
                    Ok(v) => serde_json::json!({
                        "smiles": smi,
                        "atoms": v.atoms.iter().map(|a| serde_json::json!([a.atomic_num, a.hybridization, a.chiral_tag, a.total_num_hs_with_neighbors])).collect::<Vec<_>>(),
                        "bonds": v.bonds.iter().map(|b| serde_json::json!([b.begin, b.end, b.bond_type, b.is_conjugated, b.stereo, b.stereo_atoms])).collect::<Vec<_>>(),
                        "atom_bonds": v.atom_bonds,
                        "atom_rings": v.atom_rings,
                        "bond_rings": v.bond_rings,
                    }),
                    Err(e) => serde_json::json!({"smiles": smi, "error": e.to_string()}),
                }
            }
            Err(e) => serde_json::json!({"smiles": smi, "error": e.to_string()}),
        };
        println!("{out}");
    }
}
