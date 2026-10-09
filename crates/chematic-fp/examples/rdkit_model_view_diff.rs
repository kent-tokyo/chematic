//! Rows whose RDKit-parity view the RDKit sanitization model corrects.
fn main() {
    let path = std::env::args().nth(1).expect("corpus path");
    let text = std::fs::read_to_string(path).expect("read");
    let (mut n, mut d) = (0, 0);
    for (i, line) in text.lines().enumerate() {
        let Some(s) = line.split_whitespace().next() else {
            continue;
        };
        let Ok(m) = chematic_smiles::parse(s) else {
            continue;
        };
        n += 1;
        if chematic_fp::rdkit_model_view_differs(&m) {
            d += 1;
            println!("{i}\t{s}");
        }
    }
    eprintln!("{d}/{n} rows corrected");
}
