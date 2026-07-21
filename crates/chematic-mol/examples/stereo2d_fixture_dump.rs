//! P1-A0 diagnostic: frozen fixtures probing chematic's current 2D-stereo
//! pipeline (wedge/dash bonds + 2D coordinates -> R/S / E-Z -> SMILES).
//!
//! Diagnostic only -- calls existing public APIs (`mol2000::parse_mol_with_coords`,
//! `chematic_perception::stereo2d::*`, `chematic_smiles::write`/`canonical_smiles`)
//! exactly as any external caller would; does not modify any reader, writer, or
//! the stereo2d module itself. Each fixture is a hand-built MOL V2000 block (not
//! parsed from an external file) so the exact mechanism under test is explicit
//! and reviewable in this file.
//!
//! Emits one JSON object per fixture line to stdout. Cross-referenced against
//! RDKit by `scripts/stereo2d_diagnosis.py`, which re-parses each fixture's raw
//! MOL block with RDKit and joins the two results into a failure-bucket summary.
//!
//! Run:
//! ```text
//! cargo run -p chematic-mol --example stereo2d_fixture_dump \
//!     > validation/results/stereo2d_fixture_dump.jsonl
//! ```

use chematic_core::Chirality;
use chematic_mol::parse_mol_with_coords;
use chematic_perception::{apply_stereo_from_2d, assign_ez_from_2d, assign_stereo_from_2d};
use serde_json::json;

// ---------------------------------------------------------------------------
// MOL V2000 block builder (spec-accurate fixed-width fields; see
// crates/chematic-mol/src/mol2000.rs's parser for the exact byte offsets this
// must match: x[0..10) y[10..20) z[20..30) sym[31..34) massdiff[34..36) charge[36..39)
// for atom lines; a1[0..3) a2[3..6) type[6..9) stereo[9..12) for bond lines).
// ---------------------------------------------------------------------------

struct AtomSpec {
    x: f64,
    y: f64,
    sym: &'static str,
}

struct BondSpec {
    a1: usize,  // 1-based
    a2: usize,  // 1-based
    order: u8,  // MDL bond type: 1 single, 2 double
    stereo: u8, // MDL bond stereo: 0 none, 1 wedge (up), 6 hash (down)
}

fn atom_line(a: &AtomSpec) -> String {
    format!(
        "{:>10.4}{:>10.4}{:>10.4} {:<3}{:>2}{:>3}  0  0  0  0  0  0  0  0  0  0",
        a.x, a.y, 0.0, a.sym, 0, 0
    )
}

fn bond_line(b: &BondSpec) -> String {
    format!("{:>3}{:>3}{:>3}{:>3}", b.a1, b.a2, b.order, b.stereo)
}

fn mol_block(name: &str, atoms: &[AtomSpec], bonds: &[BondSpec]) -> String {
    let mut s = String::new();
    s.push_str(name);
    s.push('\n');
    s.push_str("  chematic_diag\n");
    s.push('\n');
    s.push_str(&format!(
        "{:>3}{:>3}  0  0  0  0  0  0  0  0  0 V2000\n",
        atoms.len(),
        bonds.len()
    ));
    for a in atoms {
        s.push_str(&atom_line(a));
        s.push('\n');
    }
    for b in bonds {
        s.push_str(&bond_line(b));
        s.push('\n');
    }
    s.push_str("M  END\n");
    s
}

// ---------------------------------------------------------------------------
// Fixture definitions
// ---------------------------------------------------------------------------

/// Reusable "tripod" layout: center at origin, 3 substituents ~120 deg apart.
/// Matches the geometry already proven non-degenerate by
/// `stereo2d.rs::test_r_s_bromochlorofluoromethane`.
fn tripod_atoms(
    center_sym: &'static str,
    a: &'static str,
    b: &'static str,
    c: &'static str,
) -> Vec<AtomSpec> {
    vec![
        AtomSpec {
            x: 0.0,
            y: 0.0,
            sym: center_sym,
        },
        AtomSpec {
            x: -1.0,
            y: -0.5,
            sym: a,
        },
        AtomSpec {
            x: 1.0,
            y: -0.5,
            sym: b,
        },
        AtomSpec {
            x: 0.0,
            y: 1.0,
            sym: c,
        },
    ]
}

struct Fixture {
    id: &'static str,
    description: &'static str,
    mechanism: &'static str,
    mol: String,
}

fn fixtures() -> Vec<Fixture> {
    let mut fx = Vec::new();

    // 1. tetrahedral: 3 heavy neighbors + implicit H.
    fx.push(Fixture {
        id: "tetrahedral_3heavy_implicit_h",
        description: "C(F)(Cl)(Br) + 1 implicit H, solid wedge C->Br",
        mechanism: "tetrahedral_3heavy_implicit_h",
        mol: mol_block(
            "f1_tetrahedral_3heavy_implicit_h",
            &tripod_atoms("C", "F", "Cl", "Br"),
            &[
                BondSpec {
                    a1: 1,
                    a2: 2,
                    order: 1,
                    stereo: 0,
                },
                BondSpec {
                    a1: 1,
                    a2: 3,
                    order: 1,
                    stereo: 0,
                },
                BondSpec {
                    a1: 1,
                    a2: 4,
                    order: 1,
                    stereo: 1,
                },
            ],
        ),
    });

    // 2. tetrahedral: 4 heavy neighbors (explicit H), solid + hash wedge.
    //
    // Note: the two wedged substituents (Br, H) are deliberately placed so
    // neither is collinear-through-center with the other (NOT e.g. Br at
    // (0,1) and H at (0,-1), which RDKit's wedge parser rejects at parse
    // time with "ambiguous stereochemistry - opposing bonds have opposite
    // wedging" -- discovered by actually running this fixture against RDKit
    // during diagnosis; see the design-questions section of the report).
    // Also avoid placing F/Cl at the same y as each other, or Br/H at the
    // same y as each other: fixture #5 (wedge_atom_order_reversed) flips
    // Br's effective z sign, and an earlier, more symmetric coordinate
    // choice here made that flip land exactly on a coplanar (degenerate)
    // 4-point configuration by construction, not because of any real bug --
    // this asymmetric layout keeps the two fixtures' geometry comparable
    // without an accidental coincidence.
    {
        let atoms = vec![
            AtomSpec {
                x: 0.0,
                y: 0.0,
                sym: "C",
            },
            AtomSpec {
                x: -1.0,
                y: 0.4,
                sym: "F",
            },
            AtomSpec {
                x: 0.9,
                y: 0.7,
                sym: "Cl",
            },
            AtomSpec {
                x: -0.5,
                y: -1.1,
                sym: "Br",
            },
            AtomSpec {
                x: 0.8,
                y: -0.6,
                sym: "H",
            },
        ];
        fx.push(Fixture {
            id: "tetrahedral_4heavy_explicit_h",
            description: "C(F)(Cl)(Br)(H explicit), solid wedge C->Br, hash wedge C->H",
            mechanism: "tetrahedral_4heavy",
            mol: mol_block(
                "f2_tetrahedral_4heavy_explicit_h",
                &atoms,
                &[
                    BondSpec {
                        a1: 1,
                        a2: 2,
                        order: 1,
                        stereo: 0,
                    },
                    BondSpec {
                        a1: 1,
                        a2: 3,
                        order: 1,
                        stereo: 0,
                    },
                    BondSpec {
                        a1: 1,
                        a2: 4,
                        order: 1,
                        stereo: 1,
                    },
                    BondSpec {
                        a1: 1,
                        a2: 5,
                        order: 1,
                        stereo: 6,
                    },
                ],
            ),
        });
    }

    // 3. solid wedge (narrow end at center, standard encoding).
    fx.push(Fixture {
        id: "solid_wedge_only",
        description: "C(F)(Cl)(Br) + implicit H, solid wedge C->F, atom1=center (standard)",
        mechanism: "solid_wedge",
        mol: mol_block(
            "f3_solid_wedge_only",
            &tripod_atoms("C", "F", "Cl", "Br"),
            &[
                BondSpec {
                    a1: 1,
                    a2: 2,
                    order: 1,
                    stereo: 1,
                },
                BondSpec {
                    a1: 1,
                    a2: 3,
                    order: 1,
                    stereo: 0,
                },
                BondSpec {
                    a1: 1,
                    a2: 4,
                    order: 1,
                    stereo: 0,
                },
            ],
        ),
    });

    // 4. dashed wedge (same skeleton as #3, hash instead of solid).
    fx.push(Fixture {
        id: "dashed_wedge_only",
        description: "C(F)(Cl)(Br) + implicit H, hash wedge C->F, atom1=center (standard)",
        mechanism: "dashed_wedge",
        mol: mol_block(
            "f4_dashed_wedge_only",
            &tripod_atoms("C", "F", "Cl", "Br"),
            &[
                BondSpec {
                    a1: 1,
                    a2: 2,
                    order: 1,
                    stereo: 6,
                },
                BondSpec {
                    a1: 1,
                    a2: 3,
                    order: 1,
                    stereo: 0,
                },
                BondSpec {
                    a1: 1,
                    a2: 4,
                    order: 1,
                    stereo: 0,
                },
            ],
        ),
    });

    // 5. wedge atom1/atom2 reversed: SAME 4-heavy-neighbor geometry as
    // tetrahedral_4heavy_explicit_h (fixture #2), so assign_rs actually
    // reaches wedge_z() -- a 3-heavy+implicit-H skeleton would bail out on
    // the neighbor-count check before ever exercising the atom-order logic,
    // which is what this fixture exists to test. Only the C->Br wedge bond's
    // atom order is reversed (atom1=Br substituent, atom2=C center) versus
    // fixture #2's atom1=C, atom2=Br -- non-standard per MDL spec (which
    // requires atom1 = narrow end = stereocenter), but real-world files do
    // this. `wedge_z` in stereo2d.rs has an explicit branch for `bond.atom1
    // == center` vs not, so this checks whether that branch actually
    // produces the same 3D interpretation as the standard-order case.
    fx.push(Fixture {
        id: "wedge_atom_order_reversed",
        description: "same graph/geometry as tetrahedral_4heavy_explicit_h, but the C->Br wedge bond is written atom1=Br, atom2=C",
        mechanism: "wedge_atom_order_reversed",
        mol: mol_block(
            "f5_wedge_atom_order_reversed",
            &[
                AtomSpec { x: 0.0, y: 0.0, sym: "C" },
                AtomSpec { x: -1.0, y: 0.4, sym: "F" },
                AtomSpec { x: 0.9, y: 0.7, sym: "Cl" },
                AtomSpec { x: -0.5, y: -1.1, sym: "Br" },
                AtomSpec { x: 0.8, y: -0.6, sym: "H" },
            ],
            &[
                BondSpec { a1: 1, a2: 2, order: 1, stereo: 0 },
                BondSpec { a1: 1, a2: 3, order: 1, stereo: 0 },
                BondSpec { a1: 4, a2: 1, order: 1, stereo: 1 }, // reversed: atom1=Br, atom2=C
                BondSpec { a1: 1, a2: 5, order: 1, stereo: 6 },
            ],
        ),
    });

    // 6. multiple stereocenters: 2,3-dibromobutane, independent wedges.
    fx.push(Fixture {
        id: "multiple_stereocenters",
        description: "CH3-CHBr-CHBr-CH3, two independent 3-heavy+implicit-H stereocenters",
        mechanism: "multiple_stereocenters",
        mol: mol_block(
            "f6_multiple_stereocenters",
            &[
                AtomSpec {
                    x: -2.0,
                    y: 0.5,
                    sym: "C",
                }, // 1: C1 (methyl)
                AtomSpec {
                    x: -1.0,
                    y: 0.0,
                    sym: "C",
                }, // 2: C2 (stereocenter)
                AtomSpec {
                    x: -1.0,
                    y: 1.0,
                    sym: "Br",
                }, // 3: Br on C2
                AtomSpec {
                    x: 0.0,
                    y: 0.0,
                    sym: "C",
                }, // 4: C3 (stereocenter)
                AtomSpec {
                    x: 0.0,
                    y: -1.0,
                    sym: "Br",
                }, // 5: Br on C3
                AtomSpec {
                    x: 1.0,
                    y: 0.5,
                    sym: "C",
                }, // 6: C4 (methyl)
            ],
            &[
                BondSpec {
                    a1: 1,
                    a2: 2,
                    order: 1,
                    stereo: 0,
                },
                BondSpec {
                    a1: 2,
                    a2: 3,
                    order: 1,
                    stereo: 1,
                }, // wedge up at C2
                BondSpec {
                    a1: 2,
                    a2: 4,
                    order: 1,
                    stereo: 0,
                },
                BondSpec {
                    a1: 4,
                    a2: 5,
                    order: 1,
                    stereo: 6,
                }, // wedge down at C3
                BondSpec {
                    a1: 4,
                    a2: 6,
                    order: 1,
                    stereo: 0,
                },
            ],
        ),
    });

    // 7. negative control: same skeleton as #6, no wedges at all.
    fx.push(Fixture {
        id: "no_wedge_negative_control",
        description: "same skeleton as multiple_stereocenters but all bonds plain (stereo=0)",
        mechanism: "negative_control",
        mol: mol_block(
            "f7_no_wedge_negative_control",
            &[
                AtomSpec {
                    x: -2.0,
                    y: 0.5,
                    sym: "C",
                },
                AtomSpec {
                    x: -1.0,
                    y: 0.0,
                    sym: "C",
                },
                AtomSpec {
                    x: -1.0,
                    y: 1.0,
                    sym: "Br",
                },
                AtomSpec {
                    x: 0.0,
                    y: 0.0,
                    sym: "C",
                },
                AtomSpec {
                    x: 0.0,
                    y: -1.0,
                    sym: "Br",
                },
                AtomSpec {
                    x: 1.0,
                    y: 0.5,
                    sym: "C",
                },
            ],
            &[
                BondSpec {
                    a1: 1,
                    a2: 2,
                    order: 1,
                    stereo: 0,
                },
                BondSpec {
                    a1: 2,
                    a2: 3,
                    order: 1,
                    stereo: 0,
                },
                BondSpec {
                    a1: 2,
                    a2: 4,
                    order: 1,
                    stereo: 0,
                },
                BondSpec {
                    a1: 4,
                    a2: 5,
                    order: 1,
                    stereo: 0,
                },
                BondSpec {
                    a1: 4,
                    a2: 6,
                    order: 1,
                    stereo: 0,
                },
            ],
        ),
    });

    // 8. CIP priority tie: two CH3 branches are CIP-identical -> not a real
    // stereocenter despite a wedge being drawn.
    fx.push(Fixture {
        id: "cip_priority_tie",
        description: "C(CH3)(CH3)(F)(Cl), wedge on C->F; the two methyls tie in CIP priority",
        mechanism: "cip_priority_tie",
        mol: mol_block(
            "f8_cip_priority_tie",
            &[
                AtomSpec {
                    x: 0.0,
                    y: 0.0,
                    sym: "C",
                }, // 1: center
                AtomSpec {
                    x: -1.0,
                    y: -0.5,
                    sym: "C",
                }, // 2: methyl A
                AtomSpec {
                    x: 1.0,
                    y: -0.5,
                    sym: "C",
                }, // 3: methyl B
                AtomSpec {
                    x: 0.0,
                    y: 1.0,
                    sym: "F",
                }, // 4: F (wedged)
                AtomSpec {
                    x: 0.0,
                    y: -1.0,
                    sym: "Cl",
                }, // 5: Cl
            ],
            &[
                BondSpec {
                    a1: 1,
                    a2: 2,
                    order: 1,
                    stereo: 0,
                },
                BondSpec {
                    a1: 1,
                    a2: 3,
                    order: 1,
                    stereo: 0,
                },
                BondSpec {
                    a1: 1,
                    a2: 4,
                    order: 1,
                    stereo: 1,
                },
                BondSpec {
                    a1: 1,
                    a2: 5,
                    order: 1,
                    stereo: 0,
                },
            ],
        ),
    });

    // 9. degenerate 2D coordinates: same graph/wedge as #1, but every atom
    // sits at (0,0) -- a "lost 2D layout" file.
    fx.push(Fixture {
        id: "degenerate_2d_coordinates",
        description: "same graph as tetrahedral_3heavy_implicit_h, all atoms at (0,0)",
        mechanism: "degenerate_coordinates",
        mol: mol_block(
            "f9_degenerate_2d_coordinates",
            &[
                AtomSpec {
                    x: 0.0,
                    y: 0.0,
                    sym: "C",
                },
                AtomSpec {
                    x: 0.0,
                    y: 0.0,
                    sym: "F",
                },
                AtomSpec {
                    x: 0.0,
                    y: 0.0,
                    sym: "Cl",
                },
                AtomSpec {
                    x: 0.0,
                    y: 0.0,
                    sym: "Br",
                },
            ],
            &[
                BondSpec {
                    a1: 1,
                    a2: 2,
                    order: 1,
                    stereo: 0,
                },
                BondSpec {
                    a1: 1,
                    a2: 3,
                    order: 1,
                    stereo: 0,
                },
                BondSpec {
                    a1: 1,
                    a2: 4,
                    order: 1,
                    stereo: 1,
                },
            ],
        ),
    });

    // 10. E/Z geometry: (Z)-2-butene-equivalent layout (same coords as the
    // existing stereo2d.rs unit test, this time through the MOL reader).
    fx.push(Fixture {
        id: "ez_geometry_2butene",
        description: "4-carbon chain, central C=C, zigzag 2D layout encoding cis (Z) geometry",
        mechanism: "ez_geometry",
        mol: mol_block(
            "f10_ez_geometry_2butene",
            &[
                AtomSpec {
                    x: -0.866,
                    y: 0.5,
                    sym: "C",
                },
                AtomSpec {
                    x: 0.0,
                    y: 0.0,
                    sym: "C",
                },
                AtomSpec {
                    x: 1.5,
                    y: 0.0,
                    sym: "C",
                },
                AtomSpec {
                    x: 2.366,
                    y: 0.5,
                    sym: "C",
                },
            ],
            &[
                BondSpec {
                    a1: 1,
                    a2: 2,
                    order: 1,
                    stereo: 0,
                },
                BondSpec {
                    a1: 2,
                    a2: 3,
                    order: 2,
                    stereo: 0,
                },
                BondSpec {
                    a1: 3,
                    a2: 4,
                    order: 1,
                    stereo: 0,
                },
            ],
        ),
    });

    // 11. terminal alkene: propene, CH2=CH-CH3 -- no E/Z possible.
    fx.push(Fixture {
        id: "terminal_alkene_propene",
        description: "CH2=CH-CH3, terminal =CH2 has no heavy substituent besides the double bond",
        mechanism: "terminal_alkene",
        mol: mol_block(
            "f11_terminal_alkene_propene",
            &[
                AtomSpec {
                    x: 0.0,
                    y: 1.0,
                    sym: "C",
                }, // terminal =CH2
                AtomSpec {
                    x: 0.0,
                    y: 0.0,
                    sym: "C",
                },
                AtomSpec {
                    x: 1.0,
                    y: 0.0,
                    sym: "C",
                },
            ],
            &[
                BondSpec {
                    a1: 1,
                    a2: 2,
                    order: 2,
                    stereo: 0,
                },
                BondSpec {
                    a1: 2,
                    a2: 3,
                    order: 1,
                    stereo: 0,
                },
            ],
        ),
    });

    // 12. contradictory wedge annotations: two wedges from the same center,
    // both "up" -- over-specified / physically inconsistent.
    fx.push(Fixture {
        id: "contradictory_wedge_annotations",
        description: "C(F)(Cl)(Br) + implicit H, BOTH C->F and C->Cl marked solid wedge",
        mechanism: "contradictory_wedges",
        mol: mol_block(
            "f12_contradictory_wedge_annotations",
            &tripod_atoms("C", "F", "Cl", "Br"),
            &[
                BondSpec {
                    a1: 1,
                    a2: 2,
                    order: 1,
                    stereo: 1,
                },
                BondSpec {
                    a1: 1,
                    a2: 3,
                    order: 1,
                    stereo: 1,
                },
                BondSpec {
                    a1: 1,
                    a2: 4,
                    order: 1,
                    stereo: 0,
                },
            ],
        ),
    });

    fx
}

// ---------------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------------

fn main() {
    for f in fixtures() {
        let parsed = parse_mol_with_coords(&f.mol);
        let row = match parsed {
            Err(e) => json!({
                "id": f.id,
                "mechanism": f.mechanism,
                "description": f.description,
                "mol_block": f.mol,
                "parse_error": format!("{e:?}"),
            }),
            Ok((mol, _meta, coords)) => {
                let rs = assign_stereo_from_2d(&mol, &coords);
                let rs_json: Vec<_> = rs
                    .assignments
                    .iter()
                    .map(|(idx, code)| json!({"atom": idx.0, "cip_code": format!("{code:?}")}))
                    .collect();

                let mut ez_mol = mol.clone();
                assign_ez_from_2d(&mut ez_mol, &coords);
                let ez_json: Vec<_> = ez_mol
                    .atoms()
                    .filter_map(|(idx, a)| {
                        a.cip_code
                            .map(|c| json!({"atom": idx.0, "cip_code": format!("{c:?}")}))
                    })
                    .collect();

                let mut applied = mol.clone();
                apply_stereo_from_2d(&mut applied, &coords);
                let post_apply_fields: Vec<_> = applied
                    .atoms()
                    .map(|(idx, a)| {
                        json!({
                            "atom": idx.0,
                            "element": a.element.symbol(),
                            "cip_code": a.cip_code.map(|c| format!("{c:?}")),
                            "chirality": format!("{:?}", a.chirality),
                            "stereo_neighbor_order": applied.stereo_neighbor_order(idx).map(|v| v.to_vec()),
                        })
                    })
                    .collect();

                // Naive SMILES write immediately after apply_stereo_from_2d --
                // this is exactly what a caller would get today, demonstrating
                // whether the assigned cip_code (if any) reaches the writer.
                let naive_smiles = chematic_smiles::write(&applied);
                let naive_canonical = chematic_smiles::canonical_smiles(&applied);

                json!({
                    "id": f.id,
                    "mechanism": f.mechanism,
                    "description": f.description,
                    "mol_block": f.mol,
                    "atom_count": mol.atom_count(),
                    "coords": coords,
                    "assign_stereo_from_2d_result": rs_json,
                    "assign_ez_from_2d_result": ez_json,
                    "post_apply_stereo_from_2d_fields": post_apply_fields,
                    "naive_smiles_write": naive_smiles,
                    "naive_canonical_smiles": naive_canonical,
                    "chirality_reached_writer": applied.atoms().any(|(_, a)| a.chirality != Chirality::None),
                })
            }
        };
        println!("{row}");
    }

    // 13. coord/atom-count mismatch: not a MOL-file case (a truncated MOL
    // file would just fail to parse -- a different, uninteresting failure
    // mode) but a direct API-level misuse: call assign_stereo_from_2d with
    // fewer coordinates than the molecule has atoms.
    {
        let base = &fixtures()[1]; // tetrahedral_4heavy_explicit_h, 5 atoms
        let (mol, _meta, coords) = parse_mol_with_coords(&base.mol).expect("fixture #2 must parse");
        let short_coords = &coords[..3]; // 3 coords for a 5-atom molecule
        let rs = assign_stereo_from_2d(&mol, short_coords);
        let mut applied = mol.clone();
        apply_stereo_from_2d(&mut applied, short_coords);
        let row = json!({
            "id": "coord_atom_count_mismatch",
            "mechanism": "coord_atom_count_mismatch",
            "description": "assign_stereo_from_2d/apply_stereo_from_2d called with 3 coords for a 5-atom molecule (fixture #2's graph)",
            "mol_block": base.mol,
            "atom_count": mol.atom_count(),
            "coords_provided": short_coords,
            "assign_stereo_from_2d_result": rs.assignments.iter().map(|(idx, code)| json!({"atom": idx.0, "cip_code": format!("{code:?}")})).collect::<Vec<_>>(),
            "panicked": false,
            "note": "no panic, but NOT a safe no-op: the out-of-range neighbors (Br, H) fall back to the CENTER's own (x,y) via unwrap_or(*center_pos) in assign_rs, so their substituted position is wrong (only z from the wedge is still correct) -- yet a CIP code (S) is still returned using this corrupted geometry instead of None or an error. Whether the returned code happens to match the true answer (S, same as fixture #2's full-coords result) is coincidental to this fixture's geometry, not a property of the fallback.",
        });
        println!("{row}");
    }
}
