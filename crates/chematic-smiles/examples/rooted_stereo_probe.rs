//! Original constructed bridgehead cases for the independent RDKit CIP oracle.
#[path = "../tests/support/rooted_bridgehead.rs"]
mod fixture;

fn main() {
    let mut cases = Vec::new();
    for explicit_h in [false, true] {
        for flips in 0..4 {
            for rotation in 0..8 {
                for reverse in [false, true] {
                    let mol = fixture::bridgehead(explicit_h, flips, rotation, reverse);
                    for root in 0..7 {
                        let (text, order) =
                            chematic_smiles::write_rooted_checked(&mol, root).unwrap();
                        cases.push(serde_json::json!({
                            "explicit_h":explicit_h,"flips":flips,"rotation":rotation,"reverse":reverse,
                            "root":root,"smiles":text,"atom_order":order.into_iter().map(|i|i.0).collect::<Vec<_>>()
                        }));
                    }
                }
            }
        }
    }
    println!("{}", serde_json::to_string(&cases).unwrap());
}
