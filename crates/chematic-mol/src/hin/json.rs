use super::*;
use serde_json::{Value, json};

impl HinRecord {
    /// Versioned inspection of interpreted fields and graph-projection limits.
    pub fn to_json(&self) -> Value {
        json!({
            "schema": "chematic.hin-record.v1",
            "molecule_id": self.molecule_id,
            "title": self.title,
            "source": self.source,
            "atoms": self.atoms.iter().map(|a| json!({
                "source_id": a.source_id, "name": a.name,
                "element": a.element.symbol(), "atom_type": a.atom_type,
                "flags": a.flags, "partial_charge": a.partial_charge,
                "coordinates": a.coordinates,
            })).collect::<Vec<_>>(),
            "bonds": self.bonds.iter().map(|b| json!({
                "atom1": b.atom1, "atom2": b.atom2,
                "bond_type": match b.order {
                    BondOrder::Single => "s", BondOrder::Double => "d",
                    BondOrder::Triple => "t", BondOrder::Aromatic => "a",
                    _ => unreachable!("private HIN record only contains s/d/t/a bonds"),
                },
            })).collect::<Vec<_>>(),
            "graph_projection": {
                "formal_charges": "default_zero_not_inferred_from_partial_charges",
                "atom_types": "not_interpreted", "flags": "not_interpreted",
                "stereo": "not_perceived",
                "implicit_hydrogens": "core_valence_inference",
            },
        })
    }
}
