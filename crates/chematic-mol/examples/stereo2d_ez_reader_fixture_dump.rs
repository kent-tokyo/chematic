//! Differential-validation fixtures for P1-S2 (E/Z direction perception
//! wired into the MOL V2000/V3000 readers via
//! `chematic_perception::apply_ez_directions_from_2d_ex`).
//!
//! Companion to `stereo2d_reader_integration_fixture_dump.rs` (the
//! tetrahedral-parity reader-integration diagnosis) -- same structure, same
//! "generate via `MoleculeBuilder` + the crate's own writers, not hand-typed
//! fixed-width text" convention, this time for double-bond E/Z.
//!
//! Emits one JSON object per fixture line to stdout. Cross-referenced
//! against RDKit by `scripts/stereo2d_ez_reader_diagnosis.py`.
//!
//! Run:
//! ```text
//! cargo run -p chematic-mol --example stereo2d_ez_reader_fixture_dump \
//!     > validation/results/stereo2d_ez_reader_fixture_dump.jsonl
//! ```

use chematic_core::{Atom, BondOrder, Element, Molecule, MoleculeBuilder};
use chematic_mol::mol2000::{MolMetadata, write_mol_with_coords};
use chematic_mol::{read_mol_v3000_with_diagnostics, read_mol_with_diagnostics, write_mol_v3000};
use chematic_perception::EzDirectionDiagnostic;
use serde_json::json;

struct Fixture {
    id: &'static str,
    description: &'static str,
    dialect: &'static str, // "V2000" or "V3000"
    mol_block: String,
}

const Z_COORDS: [(f64, f64); 4] = [(-0.866, 0.5), (0.0, 0.0), (1.5, 0.0), (2.366, 0.5)];
const E_COORDS: [(f64, f64); 4] = [(-0.866, 0.5), (0.0, 0.0), (1.5, 0.0), (2.366, -0.5)];

fn but2ene(coords: [(f64, f64); 4]) -> Molecule {
    let mut b = MoleculeBuilder::new();
    let m0 = b.add_atom(Atom::new(Element::C));
    let m1 = b.add_atom(Atom::new(Element::C));
    let m2 = b.add_atom(Atom::new(Element::C));
    let m3 = b.add_atom(Atom::new(Element::C));
    b.add_bond(m0, m1, BondOrder::Single).unwrap();
    b.add_bond(m1, m2, BondOrder::Double).unwrap();
    b.add_bond(m2, m3, BondOrder::Single).unwrap();
    let _ = coords;
    b.build()
}

// The first few fixtures are pushed straight-line with no intervening
// variable declarations (unlike the later ones, each built inside its own
// `{ ... }` block) -- clippy's vec_init_then_push doesn't fit here since
// the full list mixes both shapes and reads more clearly as one uniform
// push-per-fixture sequence throughout.
#[allow(clippy::vec_init_then_push)]
fn fixtures() -> Vec<Fixture> {
    let mut fx = Vec::new();

    // 1/2. Z- and E-2-butene, V2000.
    fx.push(Fixture {
        id: "z_2butene_v2000",
        description: "(Z)-but-2-ene, V2000",
        dialect: "V2000",
        mol_block: write_mol_with_coords(
            &but2ene(Z_COORDS),
            &MolMetadata::default().with_name("but2ene"),
            &Z_COORDS,
        ),
    });
    fx.push(Fixture {
        id: "e_2butene_v2000",
        description: "(E)-but-2-ene, V2000",
        dialect: "V2000",
        mol_block: write_mol_with_coords(
            &but2ene(E_COORDS),
            &MolMetadata::default().with_name("but2ene"),
            &E_COORDS,
        ),
    });

    // 3. Same 2, V3000.
    fx.push(Fixture {
        id: "z_2butene_v3000",
        description: "(Z)-but-2-ene, V3000",
        dialect: "V3000",
        mol_block: write_mol_v3000(
            &but2ene(Z_COORDS),
            &MolMetadata::default().with_name("but2ene"),
            &Z_COORDS,
        ),
    });
    fx.push(Fixture {
        id: "e_2butene_v3000",
        description: "(E)-but-2-ene, V3000",
        dialect: "V3000",
        mol_block: write_mol_v3000(
            &but2ene(E_COORDS),
            &MolMetadata::default().with_name("but2ene"),
            &E_COORDS,
        ),
    });

    // 4. Trisubstituted alkene.
    {
        let mut b = MoleculeBuilder::new();
        let center = b.add_atom(Atom::new(Element::C));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let ch = b.add_atom(Atom::new(Element::C));
        let me = b.add_atom(Atom::new(Element::C));
        b.add_bond(center, cl, BondOrder::Single).unwrap();
        b.add_bond(center, br, BondOrder::Single).unwrap();
        b.add_bond(center, ch, BondOrder::Double).unwrap();
        b.add_bond(ch, me, BondOrder::Single).unwrap();
        let mol = b.build();
        let coords = vec![
            (0.0, 0.0),
            (-0.866, 0.5),
            (-0.866, -0.5),
            (1.5, 0.0),
            (2.366, 0.5),
        ];
        fx.push(Fixture {
            id: "trisubstituted_alkene_v2000",
            description: "Cl(Br)C=CHCH3, trisubstituted",
            dialect: "V2000",
            mol_block: write_mol_with_coords(
                &mol,
                &MolMetadata::default().with_name("trisub"),
                &coords,
            ),
        });
    }

    // 5. Tetrasubstituted alkene.
    {
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::C));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let c2 = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let i = b.add_atom(Atom::new(Element::I));
        b.add_bond(c1, cl, BondOrder::Single).unwrap();
        b.add_bond(c1, br, BondOrder::Single).unwrap();
        b.add_bond(c1, c2, BondOrder::Double).unwrap();
        b.add_bond(c2, f, BondOrder::Single).unwrap();
        b.add_bond(c2, i, BondOrder::Single).unwrap();
        let mol = b.build();
        let coords = vec![
            (0.0, 0.0),
            (-0.866, 0.5),
            (-0.866, -0.5),
            (1.5, 0.0),
            (2.366, 0.5),
            (2.366, -0.5),
        ];
        fx.push(Fixture {
            id: "tetrasubstituted_alkene_v2000",
            description: "Cl(Br)C=C(F)(I), tetrasubstituted",
            dialect: "V2000",
            mol_block: write_mol_with_coords(
                &mol,
                &MolMetadata::default().with_name("tetrasub"),
                &coords,
            ),
        });
    }

    // 6. Conjugated diene, two independent centers.
    {
        let mut b = MoleculeBuilder::new();
        let me1 = b.add_atom(Atom::new(Element::C));
        let ca = b.add_atom(Atom::new(Element::C));
        let cb = b.add_atom(Atom::new(Element::C));
        let cc = b.add_atom(Atom::new(Element::C));
        let cd = b.add_atom(Atom::new(Element::C));
        let me2 = b.add_atom(Atom::new(Element::C));
        b.add_bond(me1, ca, BondOrder::Single).unwrap();
        b.add_bond(ca, cb, BondOrder::Double).unwrap();
        b.add_bond(cb, cc, BondOrder::Single).unwrap();
        b.add_bond(cc, cd, BondOrder::Double).unwrap();
        b.add_bond(cd, me2, BondOrder::Single).unwrap();
        let mol = b.build();
        let coords = vec![
            (-2.0, 0.5),
            (-1.0, 0.0),
            (0.0, 0.5),
            (1.0, 0.0),
            (2.0, 0.5),
            (3.0, 0.0),
        ];
        fx.push(Fixture {
            id: "conjugated_diene_v2000",
            description: "(2E,4E)-hexa-2,4-diene, all-anti zigzag",
            dialect: "V2000",
            mol_block: write_mol_with_coords(
                &mol,
                &MolMetadata::default().with_name("diene"),
                &coords,
            ),
        });
    }

    // 7. Exocyclic double bond, asymmetric ring.
    {
        let mut b = MoleculeBuilder::new();
        let r1 = b.add_atom(Atom::new(Element::C));
        let r2 = b.add_atom(Atom::new(Element::C));
        let r3 = b.add_atom(Atom::new(Element::C));
        let cl = b.add_atom(Atom::new(Element::CL));
        let exo = b.add_atom(Atom::new(Element::C));
        let br = b.add_atom(Atom::new(Element::BR));
        b.add_bond(r1, r2, BondOrder::Single).unwrap();
        b.add_bond(r2, r3, BondOrder::Single).unwrap();
        b.add_bond(r3, r1, BondOrder::Single).unwrap();
        b.add_bond(r2, cl, BondOrder::Single).unwrap();
        b.add_bond(r1, exo, BondOrder::Double).unwrap();
        b.add_bond(exo, br, BondOrder::Single).unwrap();
        let mol = b.build();
        let coords = vec![
            (0.0, 0.0),
            (1.0, 0.3),
            (0.5, 1.2),
            (2.2, -0.2),
            (-1.2, -0.5),
            (-2.2, -0.2),
        ];
        fx.push(Fixture {
            id: "exocyclic_double_bond_v2000",
            description: "2-chlorocyclopropylidene-CHBr, exocyclic double bond",
            dialect: "V2000",
            mol_block: write_mol_with_coords(
                &mol,
                &MolMetadata::default().with_name("exocyclic"),
                &coords,
            ),
        });
    }

    // 8. Isotopic substituent (V3000, for isotope round-trip).
    {
        let mut b = MoleculeBuilder::new();
        let m0 = b.add_atom(Atom::new(Element::C));
        let m1 = b.add_atom(Atom::new(Element::C));
        let m2 = b.add_atom(Atom::new(Element::C));
        let mut me3 = Atom::new(Element::C);
        me3.isotope = Some(13);
        let m3 = b.add_atom(me3);
        b.add_bond(m0, m1, BondOrder::Single).unwrap();
        b.add_bond(m1, m2, BondOrder::Double).unwrap();
        b.add_bond(m2, m3, BondOrder::Single).unwrap();
        let mol = b.build();
        fx.push(Fixture {
            id: "isotopic_substituent_v3000",
            description: "but-2-ene with a 13C-labeled methyl, V3000",
            dialect: "V3000",
            mol_block: write_mol_v3000(
                &mol,
                &MolMetadata::default().with_name("isotopic"),
                &Z_COORDS,
            ),
        });
    }

    // 9. Wedge + E/Z on different bonds, same molecule.
    {
        let mut b = MoleculeBuilder::new();
        let center = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let ca = b.add_atom(Atom::new(Element::C));
        let cb = b.add_atom(Atom::new(Element::C));
        let me = b.add_atom(Atom::new(Element::C));
        b.add_bond(center, f, BondOrder::Up).unwrap();
        b.add_bond(center, cl, BondOrder::Single).unwrap();
        b.add_bond(center, ca, BondOrder::Single).unwrap();
        b.add_bond(ca, cb, BondOrder::Double).unwrap();
        b.add_bond(cb, me, BondOrder::Single).unwrap();
        let mol = b.build();
        let coords = vec![
            (0.0, 0.0),
            (-1.0, 0.4),
            (-0.5, -1.1),
            (1.2, 0.5),
            (2.4, 0.0),
            (3.4, 0.5),
        ];
        fx.push(Fixture {
            id: "wedge_and_ez_coexist_v2000",
            description: "tetrahedral wedge + independent alkene E/Z in the same molecule",
            dialect: "V2000",
            mol_block: write_mol_with_coords(
                &mol,
                &MolMetadata::default().with_name("coexist"),
                &coords,
            ),
        });
    }

    // 10. Wedge adjacent to (but not shared with) the double bond's carrier.
    {
        let mut b = MoleculeBuilder::new();
        let center = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let ca = b.add_atom(Atom::new(Element::C));
        let cb = b.add_atom(Atom::new(Element::C));
        let me = b.add_atom(Atom::new(Element::C));
        b.add_bond(center, f, BondOrder::Up).unwrap();
        b.add_bond(center, cl, BondOrder::Single).unwrap();
        b.add_bond(center, ca, BondOrder::Single).unwrap();
        b.add_bond(ca, cb, BondOrder::Double).unwrap();
        b.add_bond(cb, me, BondOrder::Single).unwrap();
        let mol = b.build();
        let coords = vec![
            (0.0, 0.0),
            (-1.0, 0.4),
            (-0.5, -1.1),
            (1.2, 0.5),
            (2.4, 0.0),
            (3.4, 0.5),
        ];
        fx.push(Fixture {
            id: "wedge_adjacent_to_double_bond_v2000",
            description: "wedge on stereocenter's OTHER substituent, adjacent to the alkene",
            dialect: "V2000",
            mol_block: write_mol_with_coords(
                &mol,
                &MolMetadata::default().with_name("adjacent"),
                &coords,
            ),
        });
    }

    // 11. Atom renumbering (reverse order).
    {
        let mut b = MoleculeBuilder::new();
        let m3 = b.add_atom(Atom::new(Element::C));
        let m2 = b.add_atom(Atom::new(Element::C));
        let m1 = b.add_atom(Atom::new(Element::C));
        let m0 = b.add_atom(Atom::new(Element::C));
        b.add_bond(m3, m2, BondOrder::Single).unwrap();
        b.add_bond(m2, m1, BondOrder::Double).unwrap();
        b.add_bond(m1, m0, BondOrder::Single).unwrap();
        let mol = b.build();
        let coords = vec![Z_COORDS[3], Z_COORDS[2], Z_COORDS[1], Z_COORDS[0]];
        fx.push(Fixture {
            id: "atom_renumbered_v2000",
            description: "same physical Z-but-2-ene, atoms declared in reverse order",
            dialect: "V2000",
            mol_block: write_mol_with_coords(
                &mol,
                &MolMetadata::default().with_name("renum"),
                &coords,
            ),
        });
    }

    // 14. Rotation.
    {
        let mol = but2ene(Z_COORDS);
        let rotated: Vec<(f64, f64)> = Z_COORDS.iter().map(|&(x, y)| (-y, x)).collect();
        fx.push(Fixture {
            id: "rotated_v2000",
            description: "same physical Z-but-2-ene, rotated 90 degrees",
            dialect: "V2000",
            mol_block: write_mol_with_coords(
                &mol,
                &MolMetadata::default().with_name("rot"),
                &rotated,
            ),
        });
    }

    // 16. Mirror reflection.
    {
        let mol = but2ene(Z_COORDS);
        let mirrored: Vec<(f64, f64)> = Z_COORDS.iter().map(|&(x, y)| (x, -y)).collect();
        fx.push(Fixture {
            id: "mirrored_v2000",
            description: "same physical Z-but-2-ene, reflected (y -> -y)",
            dialect: "V2000",
            mol_block: write_mol_with_coords(
                &mol,
                &MolMetadata::default().with_name("mirror"),
                &mirrored,
            ),
        });
    }

    fx
}

fn dump_report(
    f: &Fixture,
    mol: &Molecule,
    ez_diagnostics: &[EzDirectionDiagnostic],
) -> serde_json::Value {
    let diagnostics_json: Vec<_> = ez_diagnostics
        .iter()
        .map(|d| json!({"bond": d.bond.0, "reason": format!("{:?}", d.reason)}))
        .collect();

    json!({
        "id": f.id,
        "description": f.description,
        "dialect": f.dialect,
        "mol_block": f.mol_block,
        "atom_count": mol.atom_count(),
        "ez_diagnostics": diagnostics_json,
        "write": chematic_smiles::write(mol),
        "canonical_smiles": chematic_smiles::canonical_smiles(mol),
    })
}

fn main() {
    for f in fixtures() {
        let row = if f.dialect == "V3000" {
            match read_mol_v3000_with_diagnostics(&f.mol_block) {
                Err(e) => {
                    json!({"id": f.id, "description": f.description, "mol_block": f.mol_block, "parse_error": format!("{e:?}")})
                }
                Ok(report) => dump_report(&f, &report.mol, &report.ez_diagnostics),
            }
        } else {
            match read_mol_with_diagnostics(&f.mol_block) {
                Err(e) => {
                    json!({"id": f.id, "description": f.description, "mol_block": f.mol_block, "parse_error": format!("{e:?}")})
                }
                Ok(report) => dump_report(&f, &report.mol, &report.ez_diagnostics),
            }
        };
        println!("{row}");
    }
}
