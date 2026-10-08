//! Reads MOL blocks separated by `$$$$` lines from stdin and prints the
//! on-bits of RDKit's Avalon fingerprint (2048 bits) for each, one line per
//! molecule (`ERROR` where the block is not read).
use std::io::Read;

fn main() {
    let mut input = String::new();
    std::io::stdin().read_to_string(&mut input).unwrap();
    let mut block = String::new();
    for line in input.lines() {
        if line.starts_with("$$$$") {
            match chematic_fp::rdkit_avalon::read_molblock(&block) {
                Some(m) => {
                    let fp = chematic_fp::rdkit_avalon::avalon_fp_bytes(
                        &m,
                        2048,
                        chematic_fp::rdkit_avalon::RDKIT_AVALON_DEFAULT_BIT_FLAGS,
                    );
                    let bits: Vec<String> = chematic_fp::rdkit_avalon::on_bits(&fp)
                        .iter()
                        .map(|b| b.to_string())
                        .collect();
                    println!("{}", bits.join(" "));
                }
                None => println!("ERROR"),
            }
            block.clear();
        } else {
            block.push_str(line);
            block.push('\n');
        }
    }
}
