use std::io::{self, Read};
fn main() {
    let mut text = String::new();
    io::stdin().read_to_string(&mut text).unwrap();
    match chematic_mol::parse_gjf(&text) {
        Ok((mol, coords, charge, multiplicity)) => println!(
            "{{\"atom_count\":{},\"coords\":{:?},\"charge\":{},\"multiplicity\":{}}}",
            mol.atom_count(), coords.iter().map(|&(x,y,z)| [x,y,z]).collect::<Vec<_>>(), charge, multiplicity),
        Err(_) => std::process::exit(2),
    }
}
