//! Example: parse several SMILES strings and print molecular information.

use chematic_smiles::{parse, write};

fn main() {
    let examples = [
        ("methane", "C"),
        ("ethanol", "CCO"),
        ("acetic acid", "CC(=O)O"),
        ("benzene", "c1ccccc1"),
        ("pyridine", "c1ccncc1"),
        ("naphthalene", "c1ccc2ccccc2c1"),
        ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("caffeine", "Cn1cnc2c1c(=O)n(c(=O)n2C)C"),
        ("glucose", "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"),
        ("NaCl", "[Na+].[Cl-]"),
    ];

    for (name, smiles) in &examples {
        match parse(smiles) {
            Ok(mol) => {
                let out = write(&mol);
                println!(
                    "{:<15}  atoms={:>3}  bonds={:>3}  formula={:<12}  smiles_out={}",
                    name,
                    mol.atom_count(),
                    mol.bond_count(),
                    mol.formula(),
                    out,
                );
            }
            Err(e) => {
                eprintln!("{:<15}  ERROR: {}", name, e);
            }
        }
    }
}
