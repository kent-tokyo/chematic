//! Profiling driver: the RDKit-parity aromaticity view over a SMILES file.
fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(300);
    let what = args.get(3).map(String::as_str).unwrap_or("view");
    let smis: Vec<String> = std::fs::read_to_string(&args[1])
        .unwrap()
        .lines()
        .filter_map(|l| l.split_whitespace().next().map(str::to_string))
        .take(n)
        .collect();
    let mols: Vec<_> = smis
        .iter()
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .collect();
    for m in &mols {
        match what {
            "view" => {
                std::hint::black_box(chematic_perception::with_rdkit_parity_view(m, |v| {
                    v.is_ok()
                }));
            }
            "labute" => {
                std::hint::black_box(chematic_chem::labute_asa(m));
            }
            "tpsa" => {
                std::hint::black_box(chematic_chem::rdkit_tpsa(m));
            }
            _ => {
                std::hint::black_box(chematic_chem::rdkit_num_rings(m));
            }
        }
    }
}
