//! Dumps conformer pairs + chematic's `rmsd_symmetric` result for a handful
//! of molecules chosen to exercise graph-automorphism symmetry (terminal
//! -CF3, neopentane's 4 methyls, benzene's ring, conjugated terminal groups,
//! and isotope-distinct methyls) plus
//! one plain drug-like molecule with no interesting symmetry as a control.
//!
//! One conformer is `dg::generate_coords`'s output; the second is that same
//! conformer under a fixed rigid rotation+translation, with two
//! automorphism-equivalent atoms' positions additionally swapped where the
//! molecule has any (so the comparison genuinely exercises alignment AND
//! symmetry matching together, not a trivial self-comparison).
//!
//! Run:
//! ```text
//! cargo run --release -p chematic-3d --example rmsd_symmetric_oracle_dump \
//!   > /tmp/rmsd_symmetric_oracle_dump.jsonl
//! .venv/bin/python scripts/diagnostics/rmsd_symmetric_oracle_check.py \
//!   /tmp/rmsd_symmetric_oracle_dump.jsonl
//! ```

use chematic_3d::conformer::rmsd_symmetric;
use chematic_3d::coords::{Coords3D, Point3};
use chematic_3d::dg::generate_coords;
use chematic_core::AtomIdx;
use serde_json::json;

/// name, SMILES, (optional) pair of 0-based atom indices to swap in the
/// second conformer (symmetry-equivalent atoms, chosen by inspection of the
/// SMILES atom order).
type Case = (&'static str, &'static str, Option<(usize, usize)>);
const CASES: &[Case] = &[
    ("propane", "CCC", None),
    ("cf3_ethane", "FC(F)(F)C", Some((0, 2))), // two of the three F's
    ("neopentane", "CC(C)(C)C", Some((0, 2))), // two of the four methyl carbons
    ("benzene", "c1ccccc1", Some((0, 3))),     // two para ring carbons
    ("acetate", "CC(=O)[O-]", Some((2, 3))),   // the two formally-different O's
    ("nitro", "C[N+](=O)[O-]", Some((2, 3))),  // terminal oxygens
    ("labeled_isobutane", "[13CH3]C(C)C", Some((0, 2))), // isotope-distinct methyls
    ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", None), // no interesting symmetry, control
];

fn rotate_translate(c: &Coords3D) -> Coords3D {
    // Fixed 37 deg rotation around an arbitrary axis + a translation --
    // deterministic, not random, so this dump is reproducible byte-for-byte.
    let theta = 37f64.to_radians();
    let (cs, sn) = (theta.cos(), theta.sin());
    let n = c.atom_count();
    let mut out = Coords3D::new_zeroed(n);
    for i in 0..n {
        let p = c.get(AtomIdx(i as u32));
        // Rotate around z, then around x, then translate -- enough to mix
        // all three axes without needing a general axis-angle formula.
        let x1 = p.x * cs - p.y * sn;
        let y1 = p.x * sn + p.y * cs;
        let z1 = p.z;
        let y2 = y1 * cs - z1 * sn;
        let z2 = y1 * sn + z1 * cs;
        out.set(AtomIdx(i as u32), Point3::new(x1 + 5.0, y2 - 3.0, z2 + 1.5));
    }
    out
}

fn main() {
    for (name, smiles, swap) in CASES {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let n = mol.atom_count();
        let mut base = generate_coords(&mol);
        if *name == "labeled_isobutane" {
            // The generic DG coordinates can make the two-arm exchange a
            // rigid rotation. Use asymmetric 3D points to test that an
            // isotope-distinct exchange really stays nonzero.
            for (i, [x, y, z]) in [
                [1.0, 1.0, 1.0],
                [0.0, 0.0, 0.0],
                [-1.4, -0.8, 1.2],
                [-1.2, 1.1, -0.9],
            ]
            .into_iter()
            .enumerate()
            {
                base.set(AtomIdx(i as u32), Point3::new(x, y, z));
            }
        }
        let mut moved = rotate_translate(&base);
        if let Some((i, j)) = swap {
            let (ai, aj) = (AtomIdx(*i as u32), AtomIdx(*j as u32));
            let (pi, pj) = (moved.get(ai), moved.get(aj));
            moved.set(ai, pj);
            moved.set(aj, pi);
        }
        let symmetric = rmsd_symmetric(&mol, &base, &moved);
        let base_coords: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                let p = base.get(AtomIdx(i as u32));
                [p.x, p.y, p.z]
            })
            .collect();
        let moved_coords: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                let p = moved.get(AtomIdx(i as u32));
                [p.x, p.y, p.z]
            })
            .collect();
        println!(
            "{}",
            json!({
                "name": name,
                "smiles": smiles,
                "conformer_a": base_coords,
                "conformer_b": moved_coords,
                "chematic_rmsd_symmetric": symmetric,
            })
        );
    }
}
