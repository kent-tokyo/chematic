//! Applies reaction templates with the RDKit compatibility profile to every
//! molecule of a SMILES file: one JSON line per molecule, per reaction the
//! product sets (canonical SMILES) or the refusal.
//!
//!     cargo run --release -p chematic-rxn --example rdkit_rxn_dump -- REACTIONS.tsv FILE.smi
//!
//! REACTIONS.tsv: `name<TAB>smirks<TAB>partner SMILES separated by spaces`.

use chematic_rxn::{
    RdkitProfileOutcome, ReactionTransformLimits, run_reactants_traced_rdkit_2026_03_6,
};
use serde_json::{Map, Value, json};

fn main() {
    let mut args = std::env::args().skip(1);
    let rxns = std::fs::read_to_string(args.next().expect("reactions")).expect("read");
    let smi = std::fs::read_to_string(args.next().expect("smiles")).expect("read");
    let rxns: Vec<(String, String, Vec<chematic_core::Molecule>)> = rxns
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| {
            let mut f = l.split('\t');
            let name = f.next().unwrap().to_string();
            let smirks = f.next().unwrap().to_string();
            let partners = f
                .next()
                .unwrap_or("")
                .split_whitespace()
                .map(|s| chematic_smiles::parse(s).expect("partner"))
                .collect();
            (name, smirks, partners)
        })
        .collect();
    let limits = ReactionTransformLimits::default();
    for line in smi.lines() {
        let s = line.split_whitespace().next().unwrap_or("");
        if s.is_empty() {
            continue;
        }
        let Ok(mol) = chematic_smiles::parse(s) else {
            println!("{}", json!({"smiles": s, "parse": false}));
            continue;
        };
        let mut out = Map::new();
        for (name, smirks, partners) in &rxns {
            let mut refs = vec![&mol];
            refs.extend(partners.iter());
            let v = match run_reactants_traced_rdkit_2026_03_6(smirks, &refs, &limits) {
                Ok(RdkitProfileOutcome::Report(r)) => json!({"ok": r.products.iter().map(|ps|
                    ps.iter().map(|p| chematic_smiles::canonical_smiles(&p.molecule)).collect::<Vec<_>>()
                ).collect::<Vec<_>>(), "rejected": r.diagnostics.valence_rejected_matches}),
                Ok(RdkitProfileOutcome::Unsupported(u)) => {
                    json!({"err": format!("unsupported: {}", u.reason_code())})
                }
                Err(e) => json!({"err": e.to_string()}),
            };
            out.insert(name.clone(), v);
        }
        println!("{}", json!({"smiles": s, "rxn": Value::Object(out)}));
    }
}
