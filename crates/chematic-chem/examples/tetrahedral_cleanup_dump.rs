//! Original validation adapter: expose complete source/output graphs and frames.
use chematic_core::{AtomIdx, Chirality, Molecule};
use serde_json::{Value, json};
use std::io::{self, Read};
fn snapshot(m: &Molecule) -> Value {
    json!({"atoms":m.atoms().map(|(i,a)|json!({
        "element":a.element.atomic_number(),"isotope":a.isotope.unwrap_or(0),
        "charge":a.charge,"hydrogens":chematic_core::implicit_hcount(m,i),
        "aromatic":a.aromatic,"map":a.atom_map.unwrap_or(0),"tag":m.atom_tag(i).map(|v|v.get()),
        "chirality":format!("{:?}",a.chirality),"cip":format!("{:?}",a.cip_code),
        "frame":m.stereo_neighbor_order(i)})).collect::<Vec<_>>(),
        "bonds":m.bonds().map(|(_,b)|json!({"begin":b.atom1.0,"end":b.atom2.0,"order":format!("{:?}",b.order)})).collect::<Vec<_>>()})
}
fn main() {
    let mut input = String::new();
    io::stdin().read_to_string(&mut input).unwrap();
    let cases: Vec<Value> = serde_json::from_str(&input).unwrap();
    let mut rows = Vec::new();
    for c in cases {
        let source = chematic_smiles::parse(c["smiles"].as_str().unwrap()).unwrap();
        let perm: Vec<_> = c["order"]
            .as_array()
            .unwrap()
            .iter()
            .map(|i| AtomIdx(i.as_u64().unwrap() as u32))
            .collect();
        let original = source.renumber_atoms(&perm).unwrap();
        let target = AtomIdx(
            perm.iter()
                .position(|a| a.0 == c["atom_index"].as_u64().unwrap() as u32)
                .unwrap() as u32,
        );
        let result =
            chematic_chem::remove_redundant_tetrahedral_stereo(&original, target, 1_000_000)
                .unwrap();
        let mut inverted = original.clone();
        inverted.set_chirality(
            target,
            match original.atom(target).chirality {
                Chirality::Clockwise => Chirality::CounterClockwise,
                Chirality::CounterClockwise => Chirality::Clockwise,
                _ => unreachable!(),
            },
        );
        rows.push(json!({"input":c,"source":snapshot(&original),"inverted":snapshot(&inverted),"output":snapshot(&result.molecule),
            "report":{"removed":result.removed,"atom_index":target.0,"search_visits":result.search_visits,"witness":result.witness.map(|v|v.into_iter().map(|a|a.0).collect::<Vec<_>>())},
            "source_smiles":chematic_smiles::canonical_smiles(&original),"output_smiles":chematic_smiles::canonical_smiles(&result.molecule)}));
    }
    println!("{}", serde_json::to_string(&rows).unwrap());
}
