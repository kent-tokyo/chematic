//! `rdkit_hydrogen_suppressed` must keep RDKit's molecule: its RDKit SMILES
//! equals the input's. `cargo run --example h_suppress_check -- corpus.smi`
fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let smis: Vec<String> = if args.len() == 1 && args[0].ends_with(".smi") {
        std::fs::read_to_string(&args[0])
            .unwrap()
            .lines()
            .filter_map(|l| l.split_whitespace().next().map(String::from))
            .collect()
    } else {
        args
    };
    let (mut n, mut bad) = (0, 0);
    for s in &smis {
        let Ok(m) = chematic_smiles::parse(s) else {
            continue;
        };
        let Some(h) = chematic_smiles::rdkit_hydrogen_suppressed(&m) else {
            continue;
        };
        n += 1;
        let a = chematic_smiles::rdkit_canonical_smiles(&m);
        let b = chematic_smiles::rdkit_canonical_smiles(&h);
        if a != b {
            bad += 1;
            println!("{s}\t{a:?}\t{b:?}");
        }
    }
    eprintln!("{bad}/{n} differ");
}
