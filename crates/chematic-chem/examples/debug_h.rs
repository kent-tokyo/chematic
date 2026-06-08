use chematic_chem::molecular_weight;
use chematic_core::{apply_kekule, implicit_hcount, kekulize};
use chematic_smiles::parse;

fn debug(name: &str, smiles: &str) {
    let mol = parse(smiles).unwrap();
    let kek = match kekulize(&mol) {
        Ok(k) => k,
        Err(e) => {
            println!("{name}: kekulize FAILED: {e}");
            return;
        }
    };
    let kmol = apply_kekule(&mol, &kek);

    let mut total_h = 0;
    let mut issues = vec![];
    for (idx, atom) in kmol.atoms() {
        let h = implicit_hcount(&kmol, idx);
        let bonds: Vec<_> = kmol
            .neighbors(idx)
            .map(|(_, bidx)| format!("{:?}", kmol.bond(bidx).order))
            .collect();
        total_h += h as u32;
        if h > 1 && atom.element.symbol() == "C" && atom.aromatic {
            issues.push(format!(
                "  SUSPECT: atom {:2} ({}) aromatic={} h={} bonds={:?}",
                idx.0,
                atom.element.symbol(),
                atom.aromatic,
                h,
                bonds
            ));
        }
    }
    println!(
        "{name}: MW={:.4} (expected 308.333) total_implicit_H={}",
        molecular_weight(&mol),
        total_h
    );
    for s in issues {
        println!("{}", s);
    }
}

fn main() {
    debug("warfarin", "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O");
}
