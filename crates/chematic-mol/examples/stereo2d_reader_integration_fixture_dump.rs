//! Differential-validation fixtures for the *reader-integration* work that
//! wires `chematic_perception::apply_local_parity_from_wedges_with_diagnostics`
//! into the V2000/V3000 MOL readers (`chematic_mol::read_mol_with_diagnostics`
//! / `read_mol_v3000_with_diagnostics`).
//!
//! This is a companion to, not a replacement for, `stereo2d_fixture_dump.rs`
//! (the original P1-A0 diagnosis of the *parity math itself*, still frozen
//! and accurate). This file instead validates that wiring that already-
//! calibrated math into the production reader path produces the right
//! results end-to-end: V2000/V3000 agreement, renumbering/reflection/
//! rotation invariance, round-trip losslessness, and the new structured
//! `StereoDiagnostic` API -- categories the original 14-fixture diagnosis
//! did not need to cover, since it predates both the reader wiring and the
//! diagnostics API.
//!
//! Fixtures are built via `MoleculeBuilder` + the crate's own
//! `write_mol_with_coords`/`write_mol_v3000` writers (not hand-typed
//! fixed-width text), so every MOL block is byte-valid by construction.
//!
//! Emits one JSON object per fixture line to stdout. Cross-referenced
//! against RDKit by `scripts/stereo2d_reader_diagnosis.py`.
//!
//! Run:
//! ```text
//! cargo run -p chematic-mol --example stereo2d_reader_integration_fixture_dump \
//!     > validation/results/stereo2d_reader_fixture_dump.jsonl
//! ```

use chematic_chem::{CipMode, assign_cip_with_mode};
use chematic_core::{Atom, BondOrder, Element, Molecule, MoleculeBuilder};
use chematic_mol::mol2000::{MolMetadata, write_mol_with_coords};
use chematic_mol::{read_mol_v3000_with_diagnostics, read_mol_with_diagnostics, write_mol_v3000};
use serde_json::json;

/// Asymmetric, non-degenerate 4-position layout (matches
/// `chematic-perception`'s own `quad_positions()`).
fn quad_positions() -> [(f64, f64); 4] {
    [(-1.0, 0.4), (0.9, 0.7), (-0.5, -1.1), (0.8, -0.6)]
}

struct Fixture {
    id: &'static str,
    description: &'static str,
    dialect: &'static str, // "V2000" or "V3000"
    mol_block: String,
}

fn chfclbr_wedge_v2000() -> String {
    let mut b = MoleculeBuilder::new();
    let c = b.add_atom(Atom::new(Element::C));
    let f = b.add_atom(Atom::new(Element::F));
    let cl = b.add_atom(Atom::new(Element::CL));
    let br = b.add_atom(Atom::new(Element::BR));
    let i = b.add_atom(Atom::new(Element::I));
    b.add_bond(c, f, BondOrder::Up).unwrap();
    b.add_bond(c, cl, BondOrder::Single).unwrap();
    b.add_bond(c, br, BondOrder::Single).unwrap();
    b.add_bond(c, i, BondOrder::Single).unwrap();
    let mol = b.build();
    let quad = quad_positions();
    let coords = vec![(0.0, 0.0), quad[0], quad[1], quad[2], quad[3]];
    write_mol_with_coords(
        &mol,
        &MolMetadata::default().with_name("valid_wedge"),
        &coords,
    )
}

fn fixtures() -> Vec<Fixture> {
    let mut fx = Vec::new();
    let quad = quad_positions();

    // 1. Valid wedge, V2000 -- the baseline case.
    fx.push(Fixture {
        id: "valid_wedge_v2000",
        description: "C(F)(Cl)(Br)(I), solid wedge C->F, V2000",
        dialect: "V2000",
        mol_block: chfclbr_wedge_v2000(),
    });

    // 2. Same physical molecule, V3000.
    {
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let i = b.add_atom(Atom::new(Element::I));
        b.add_bond(c, f, BondOrder::Up).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        b.add_bond(c, br, BondOrder::Single).unwrap();
        b.add_bond(c, i, BondOrder::Single).unwrap();
        let mol = b.build();
        let coords = vec![(0.0, 0.0), quad[0], quad[1], quad[2], quad[3]];
        fx.push(Fixture {
            id: "valid_wedge_v3000",
            description: "same physical molecule as valid_wedge_v2000, V3000 CFG=1",
            dialect: "V3000",
            mol_block: write_mol_v3000(
                &mol,
                &MolMetadata::default().with_name("valid_wedge_v3000"),
                &coords,
            ),
        });
    }

    // 3. Contradictory wedges (two disagreeing wedges from the same center).
    {
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let i = b.add_atom(Atom::new(Element::I));
        b.add_bond(c, f, BondOrder::Up).unwrap();
        b.add_bond(c, cl, BondOrder::Up).unwrap();
        b.add_bond(c, br, BondOrder::Single).unwrap();
        b.add_bond(c, i, BondOrder::Single).unwrap();
        let mol = b.build();
        let coords = vec![(0.0, 0.0), quad[0], quad[1], quad[2], quad[3]];
        fx.push(Fixture {
            id: "contradictory_wedge_v2000",
            description: "C(F)(Cl)(Br)(I), F and Cl both solid wedge (disagreeing), V2000",
            dialect: "V2000",
            mol_block: write_mol_with_coords(
                &mol,
                &MolMetadata::default().with_name("contradictory"),
                &coords,
            ),
        });
    }

    // 4. Atom renumbering: same physical molecule/wedge, atoms added in the
    // opposite order (coords reassigned to match each atom's real position).
    {
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let i = b.add_atom(Atom::new(Element::I));
        let br = b.add_atom(Atom::new(Element::BR));
        let cl = b.add_atom(Atom::new(Element::CL));
        let f = b.add_atom(Atom::new(Element::F));
        b.add_bond(c, i, BondOrder::Single).unwrap();
        b.add_bond(c, br, BondOrder::Single).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        b.add_bond(c, f, BondOrder::Up).unwrap();
        let mol = b.build();
        let coords = vec![(0.0, 0.0), quad[3], quad[2], quad[1], quad[0]];
        fx.push(Fixture {
            id: "atom_renumbered_v2000",
            description: "same physical molecule/wedge as valid_wedge_v2000, atoms added in reverse order",
            dialect: "V2000",
            mol_block: write_mol_with_coords(&mol, &MolMetadata::default().with_name("renumbered"), &coords),
        });
    }

    // 5. Reflection (y -> -y): must flip parity.
    {
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let i = b.add_atom(Atom::new(Element::I));
        b.add_bond(c, f, BondOrder::Up).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        b.add_bond(c, br, BondOrder::Single).unwrap();
        b.add_bond(c, i, BondOrder::Single).unwrap();
        let mol = b.build();
        let coords: Vec<(f64, f64)> = [(0.0, 0.0), quad[0], quad[1], quad[2], quad[3]]
            .iter()
            .map(|&(x, y)| (x, -y))
            .collect();
        fx.push(Fixture {
            id: "reflected_v2000",
            description: "same physical layout as valid_wedge_v2000 but reflected (y -> -y); must flip parity",
            dialect: "V2000",
            mol_block: write_mol_with_coords(&mol, &MolMetadata::default().with_name("reflected"), &coords),
        });
    }

    // 6. Rotation + translation: must preserve parity.
    {
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let i = b.add_atom(Atom::new(Element::I));
        b.add_bond(c, f, BondOrder::Up).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        b.add_bond(c, br, BondOrder::Single).unwrap();
        b.add_bond(c, i, BondOrder::Single).unwrap();
        let mol = b.build();
        // 90-degree rotation + translation: (x,y) -> (-y+5, x+5).
        let coords: Vec<(f64, f64)> = [(0.0, 0.0), quad[0], quad[1], quad[2], quad[3]]
            .iter()
            .map(|&(x, y)| (-y + 5.0, x + 5.0))
            .collect();
        fx.push(Fixture {
            id: "rotated_translated_v2000",
            description: "same physical layout as valid_wedge_v2000, 90-degree rotation + translation; must preserve parity",
            dialect: "V2000",
            mol_block: write_mol_with_coords(&mol, &MolMetadata::default().with_name("rotated"), &coords),
        });
    }

    // 7. Bond declaration order reversed (even permutation for 4 neighbors):
    // must preserve parity.
    {
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let i = b.add_atom(Atom::new(Element::I));
        b.add_bond(c, i, BondOrder::Single).unwrap();
        b.add_bond(c, br, BondOrder::Single).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        b.add_bond(c, f, BondOrder::Up).unwrap();
        let mol = b.build();
        let coords = vec![(0.0, 0.0), quad[0], quad[1], quad[2], quad[3]];
        fx.push(Fixture {
            id: "bond_order_reversed_v2000",
            description: "same physical molecule as valid_wedge_v2000, bond block declared in reverse order",
            dialect: "V2000",
            mol_block: write_mol_with_coords(&mol, &MolMetadata::default().with_name("bond_reversed"), &coords),
        });
    }

    // 8. Charged tetravalent N+ center, 4 distinct heavy substituents (zero
    // implicit H) -- synthetic (not a plausible real molecule) but exercises
    // the 4-explicit-neighbor path on a charged heteroatom, as the task
    // requires ("charged N ... tetrahedral examples where supported").
    {
        let mut n_atom = Atom::new(Element::N);
        n_atom.charge = 1;
        let mut b = MoleculeBuilder::new();
        let n = b.add_atom(n_atom);
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let i = b.add_atom(Atom::new(Element::I));
        b.add_bond(n, f, BondOrder::Up).unwrap();
        b.add_bond(n, cl, BondOrder::Single).unwrap();
        b.add_bond(n, br, BondOrder::Single).unwrap();
        b.add_bond(n, i, BondOrder::Single).unwrap();
        let mol = b.build();
        let coords = vec![(0.0, 0.0), quad[0], quad[1], quad[2], quad[3]];
        fx.push(Fixture {
            id: "charged_n_center_v2000",
            description: "[N+](F)(Cl)(Br)(I), solid wedge N->F -- charged, zero implicit H, 4 explicit neighbors",
            dialect: "V2000",
            mol_block: write_mol_with_coords(&mol, &MolMetadata::default().with_name("charged_n"), &coords),
        });
    }

    // 9. Isotopic substituent: 13C stereocenter (isotope must survive
    // alongside stereo perception; local parity is purely geometric and
    // must be unaffected by isotope).
    {
        let mut c_atom = Atom::new(Element::C);
        c_atom.isotope = Some(13);
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(c_atom);
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let i = b.add_atom(Atom::new(Element::I));
        b.add_bond(c, f, BondOrder::Up).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        b.add_bond(c, br, BondOrder::Single).unwrap();
        b.add_bond(c, i, BondOrder::Single).unwrap();
        let mol = b.build();
        let coords = vec![(0.0, 0.0), quad[0], quad[1], quad[2], quad[3]];
        fx.push(Fixture {
            id: "isotopic_stereocenter_v2000",
            description: "[13C](F)(Cl)(Br)(I), solid wedge C->F -- isotope-13 stereocenter",
            dialect: "V2000",
            mol_block: write_mol_with_coords(
                &mol,
                &MolMetadata::default().with_name("isotopic"),
                &coords,
            ),
        });
    }

    // 10. Ring stereocenter: bromocyclopropane (C1 bonded to 2 ring
    // neighbors + Br + implicit H = 3 heavy neighbors + implicit H, but
    // ring-constrained, unlike any of the 22 existing perception unit tests).
    {
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let c3 = b.add_atom(Atom::new(Element::C));
        let br = b.add_atom(Atom::new(Element::BR));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        b.add_bond(c2, c3, BondOrder::Single).unwrap();
        b.add_bond(c3, c1, BondOrder::Single).unwrap();
        b.add_bond(c1, br, BondOrder::Up).unwrap();
        let mol = b.build();
        // Equilateral-ish triangle for C1/C2/C3, Br off to one side.
        let coords = vec![(0.0, 0.0), (1.0, 0.3), (0.5, 1.2), (-1.0, -0.6)];
        fx.push(Fixture {
            id: "ring_stereocenter_v2000",
            description: "bromocyclopropane: C1(ring)(ring)(Br) + implicit H, solid wedge C1->Br",
            dialect: "V2000",
            mol_block: write_mol_with_coords(
                &mol,
                &MolMetadata::default().with_name("ring"),
                &coords,
            ),
        });
    }

    // 11. Multiple stereocenters: 2,3-dibromobutane, two independent centers.
    {
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let br2 = b.add_atom(Atom::new(Element::BR));
        let c3 = b.add_atom(Atom::new(Element::C));
        let br3 = b.add_atom(Atom::new(Element::BR));
        let c4 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        b.add_bond(c2, br2, BondOrder::Up).unwrap();
        b.add_bond(c2, c3, BondOrder::Single).unwrap();
        b.add_bond(c3, br3, BondOrder::Down).unwrap();
        b.add_bond(c3, c4, BondOrder::Single).unwrap();
        let mol = b.build();
        let coords = vec![
            (-2.0, 0.5),
            (-1.0, 0.0),
            (-1.0, 1.0),
            (0.0, 0.0),
            (0.0, -1.0),
            (1.0, 0.5),
        ];
        fx.push(Fixture {
            id: "multi_stereocenter_v2000",
            description: "2,3-dibromobutane, two independent 3-heavy+implicit-H stereocenters",
            dialect: "V2000",
            mol_block: write_mol_with_coords(
                &mol,
                &MolMetadata::default().with_name("multi"),
                &coords,
            ),
        });
    }

    // 12. Achiral negative control: same skeleton as #1, no wedges at all.
    {
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let i = b.add_atom(Atom::new(Element::I));
        b.add_bond(c, f, BondOrder::Single).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        b.add_bond(c, br, BondOrder::Single).unwrap();
        b.add_bond(c, i, BondOrder::Single).unwrap();
        let mol = b.build();
        let coords = vec![(0.0, 0.0), quad[0], quad[1], quad[2], quad[3]];
        fx.push(Fixture {
            id: "achiral_negative_control_v2000",
            description: "C(F)(Cl)(Br)(I), no wedge/hash bonds at all",
            dialect: "V2000",
            mol_block: write_mol_with_coords(
                &mol,
                &MolMetadata::default().with_name("achiral"),
                &coords,
            ),
        });
    }

    // 13. Wedge on a bond that does not represent tetrahedral stereochemistry
    // (incident to a 2-coordinate atom): dimethyl ether-like fragment,
    // wedge on the C-O bond. O has degree 2 -- not a candidate shape at all,
    // must be silently NotRequested (no diagnostic, no assignment).
    {
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::C));
        let o = b.add_atom(Atom::new(Element::O));
        let c2 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c1, o, BondOrder::Up).unwrap();
        b.add_bond(o, c2, BondOrder::Single).unwrap();
        let mol = b.build();
        let coords = vec![(0.0, 0.0), (1.0, 0.0), (2.0, 0.3)];
        fx.push(Fixture {
            id: "non_tetrahedral_wedge_v2000",
            description: "CH3-O-CH3-like fragment with a wedge on the C-O bond; O has degree 2, not a stereocenter candidate",
            dialect: "V2000",
            mol_block: write_mol_with_coords(&mol, &MolMetadata::default().with_name("nontet"), &coords),
        });
    }

    fx
}

fn dump_report(
    f: &Fixture,
    mol: &Molecule,
    coords: &[(f64, f64)],
    stereo_diagnostics: &[chematic_perception::StereoDiagnostic],
) -> serde_json::Value {
    let chirality_fields: Vec<_> = mol
        .atoms()
        .filter_map(|(idx, a)| {
            if a.chirality == chematic_core::Chirality::None {
                None
            } else {
                Some(json!({
                    "atom": idx.0,
                    "element": a.element.symbol(),
                    "chirality": format!("{:?}", a.chirality),
                    "stereo_neighbor_order": mol.stereo_neighbor_order(idx).map(|v| v.to_vec()),
                }))
            }
        })
        .collect();

    let diagnostics_json: Vec<_> = stereo_diagnostics
        .iter()
        .map(|d| json!({"atom": d.atom.0, "reason": format!("{:?}", d.reason)}))
        .collect();

    let cip = assign_cip_with_mode(mol, CipMode::Accurate);
    let cip_json = cip.ok().map(|assignment| {
        assignment
            .assignments
            .iter()
            .map(|(idx, code)| json!({"atom": idx.0, "cip_code": format!("{code:?}")}))
            .collect::<Vec<_>>()
    });

    json!({
        "id": f.id,
        "description": f.description,
        "dialect": f.dialect,
        "mol_block": f.mol_block,
        "atom_count": mol.atom_count(),
        "coords": coords,
        "chirality": chirality_fields,
        "stereo_diagnostics": diagnostics_json,
        "accurate_cip_labels": cip_json,
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
                Ok(report) => {
                    dump_report(&f, &report.mol, &report.coords, &report.stereo_diagnostics)
                }
            }
        } else {
            match read_mol_with_diagnostics(&f.mol_block) {
                Err(e) => {
                    json!({"id": f.id, "description": f.description, "mol_block": f.mol_block, "parse_error": format!("{e:?}")})
                }
                Ok(report) => {
                    dump_report(&f, &report.mol, &report.coords, &report.stereo_diagnostics)
                }
            }
        };
        println!("{row}");
    }
}
