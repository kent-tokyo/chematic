use chematic_chem::{heavy_atom_count, logp_crippen, molecular_weight, tpsa};
use chematic_smiles::parse;

const PILOT: &[(&str, &str)] = &[
    ("benzene", "c1ccccc1"),
    ("toluene", "Cc1ccccc1"),
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("caffeine", "Cn1cnc2c1c(=O)n(c(=O)n2C)C"),
    ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
    ("paracetamol", "CC(=O)Nc1ccc(O)cc1"),
    ("pyridine", "c1ccncc1"),
    ("aniline", "Nc1ccccc1"),
    ("phenol", "Oc1ccccc1"),
    ("ethanol", "CCO"),
    ("acetic_acid", "CC(=O)O"),
    ("glucose", "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"),
    ("naphthalene", "c1ccc2ccccc2c1"),
    ("indole", "c1ccc2[nH]ccc2c1"),
    ("imidazole", "c1cnc[nH]1"),
    (
        "morphine",
        "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5",
    ),
    ("chlorpromazine", "CN(C)CCCN1c2ccccc2Sc2ccc(Cl)cc21"),
    ("captopril", "CC(CS)C(=O)N1CCCC1C(=O)O"),
    ("methotrexate", "Cn1cnc2nc(N)nc(N)c2c1=O"),
    ("warfarin", "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O"),
];

fn main() {
    println!(
        "{:<18} {:>10} {:>10} {:>10} {:>6}",
        "name", "mw", "logp", "tpsa", "heavy"
    );
    for (name, smi) in PILOT {
        match parse(smi) {
            Ok(mol) => {
                let mw = molecular_weight(&mol);
                let logp = logp_crippen(&mol);
                let tp = tpsa(&mol);
                let heavy = heavy_atom_count(&mol);
                println!(
                    "{:<18} {:>10.4} {:>10.4} {:>10.4} {:>6}",
                    name, mw, logp, tp, heavy
                );
            }
            Err(e) => println!("{:<18} PARSE_ERR: {}", name, e),
        }
    }
}
