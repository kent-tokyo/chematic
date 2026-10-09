//! Fast `rdkit_num_rings` against RDKit's full symmetrized SSSR on SMILES files.
fn main() {
    let (mut n, mut bad) = (0usize, 0usize);
    for path in std::env::args().skip(1) {
        for line in std::fs::read_to_string(&path).unwrap().lines() {
            let Some(s) = line.split_whitespace().next() else {
                continue;
            };
            let Ok(m) = chematic_smiles::parse(s) else {
                continue;
            };
            n += 1;
            let full = chematic_perception::rdkit_sssr_ring_order(&m)
                .map(|r| r.len())
                .unwrap_or_else(|| chematic_perception::find_symmetrized_sssr(&m).rings().len());
            if chematic_chem::rdkit_num_rings(&m) != full {
                bad += 1;
                if bad < 20 {
                    println!(
                        "MISMATCH {s} fast {} full {full}",
                        chematic_chem::rdkit_num_rings(&m)
                    );
                }
            }
        }
    }
    println!("{n} molecules, {bad} mismatches");
}
