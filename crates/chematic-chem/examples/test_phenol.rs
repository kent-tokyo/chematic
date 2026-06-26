use chematic_chem::descriptors::logp_crippen;
use chematic_smiles::parse;

fn main() {
    let mol = parse("c1ccccc1O").unwrap();
    let logp = logp_crippen(&mol);
    println!("Phenol LogP: {}", logp);
    println!("Expected: 1.3922");
    println!("Diff: {}", logp - 1.3922);
}
