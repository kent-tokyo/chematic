//! Focused regression measurement for issue #210's legacy-coordinate UFF rescue.
//!
//! Every named molecule is parsed from the issue's frozen SMILES, starts from
//! `generate_coords`, and must either return a finite, stereo-satisfied geometry
//! or be reported as a failure. The process exits unsuccessfully while any named
//! residual remains, making this suitable for the Experimental 3D long-run lane.
//!
//! Run:
//! `cargo run --release -p chematic-3d --example issue210_rescue_measurement`

use chematic_3d::generate_coords;
use chematic_3d::minimize::{EnergyReport, ForceFieldPolicy, MinimizeConfig, minimize_with_policy};
use chematic_3d::stereo_constraints::verify_stereo;
use chematic_smiles::parse;

const CASES: &[(&str, &str)] = &[
    ("ibuprofen_S", "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O"),
    ("naproxen_S", "COc1ccc2cc([C@H](C)C(=O)O)ccc2c1"),
    (
        "testosterone",
        "C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O",
    ),
    (
        "cholesterol",
        "C[C@H](CCCC(C)C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C",
    ),
    (
        "atorvastatin_fragment",
        "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O",
    ),
];

fn main() {
    let config = MinimizeConfig::default();
    let mut failed = 0usize;

    for &(name, smiles) in CASES {
        let mol = parse(smiles).unwrap_or_else(|error| panic!("parse {name}: {error}"));
        let initial = generate_coords(&mol);
        match minimize_with_policy(&mol, initial, ForceFieldPolicy::UffOnly, &config) {
            Ok(result) => {
                let stereo = verify_stereo(&mol, &result.coords);
                let energy_descended = matches!(
                    (result.energy_before, result.energy_after),
                    (
                        EnergyReport::Uff { total: before },
                        EnergyReport::Uff { total: after }
                    ) if after <= before
                );
                let accepted =
                    result.coords.is_finite() && stereo.is_fully_satisfied() && energy_descended;
                println!(
                    "{name}: accepted={accepted} start={:?} energy_descended={energy_descended} stereo={stereo:?}",
                    result.starting_geometry
                );
                if !accepted {
                    failed += 1;
                }
            }
            Err(error) => {
                println!("{name}: accepted=false error={error:?}");
                failed += 1;
            }
        }
    }

    assert_eq!(failed, 0, "issue #210 still has {failed} named residual(s)");
}
