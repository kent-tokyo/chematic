//! `MolToSmarts(MolFromSmarts(s))` / `ReactionToSmarts(ReactionFromSmarts(s))`
//! of the RDKit port for each input line (lines containing `>` are
//! reactions): one JSON object per line.
use serde_json::json;
fn main() {
    let path = std::env::args().nth(1).expect("file");
    for line in std::fs::read_to_string(path).expect("read").lines() {
        if line.is_empty() {
            continue;
        }
        let r = if line.contains('>') {
            chematic_smiles::rdkit_reaction_to_smarts(line)
        } else {
            chematic_smiles::rdkit_smarts_to_smarts(line)
        };
        let v = match r {
            Ok(s) => json!({"ok": s}),
            Err(e) => json!({"err": e.to_string()}),
        };
        println!("{}", json!({"in": line, "out": v}));
    }
}
