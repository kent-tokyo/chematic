use std::io::{self, BufRead};

fn main() {
    for (index, line) in io::stdin().lock().lines().enumerate() {
        let line = line.expect("stdin");
        let smiles = line.split_whitespace().next().unwrap_or("");
        if smiles.is_empty() {
            println!("{index}\tE");
            continue;
        }
        let result = chematic_smiles::parse(smiles)
            .map(|mol| chematic_chem::add_hydrogens(&mol))
            .map_err(|err| err.to_string())
            .and_then(|mol| {
                let types = chematic_ff::assign_mmff94_numeric_types(&mol)
                    .map_err(|err| err.to_string())?;
                Ok(mol
                    .atoms()
                    .map(|(idx, atom)| {
                        format!("{}:{}", atom.element.atomic_number(), types[idx.0 as usize])
                    })
                    .collect::<Vec<_>>())
            });
        match result {
            Ok(types) => println!("{index}\t{}", types.join(",")),
            Err(_) => println!("{index}\tE"),
        }
    }
}
