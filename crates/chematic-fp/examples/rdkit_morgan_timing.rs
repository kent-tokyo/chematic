//! Focused timing driver for the RDKit-compatible Morgan public paths.

use std::time::Instant;

use chematic_fp::{
    RdkitMorganConfig, RdkitMorganFpSize, RdkitMorganRadius, rdkit_morgan_bitvec,
    rdkit_morgan_fingerprint,
};

fn carries_stereo(molecule: &chematic_core::Molecule) -> bool {
    molecule
        .atoms()
        .any(|(_, atom)| atom.chirality != chematic_core::Chirality::None)
        || molecule.bonds().any(|(_, bond)| {
            matches!(
                bond.order,
                chematic_core::BondOrder::Up | chematic_core::BondOrder::Down
            )
        })
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(1_000);
    let input = std::fs::read_to_string(&args[1]).unwrap();
    let smiles: Vec<_> = input
        .lines()
        .filter_map(|line| line.split_whitespace().next())
        .take(n)
        .collect();
    for include_chirality in [false, true] {
        let config = RdkitMorganConfig {
            radius: RdkitMorganRadius::R2,
            fp_size: RdkitMorganFpSize::B2048,
            include_chirality,
        };
        let verification_molecules: Vec<_> = smiles
            .iter()
            .filter_map(|smiles| chematic_smiles::parse(smiles).ok())
            .collect();
        let mut mismatches = 0usize;
        for molecule in &verification_molecules {
            let bits = rdkit_morgan_bitvec(molecule, &config).unwrap();
            let detail = rdkit_morgan_fingerprint(molecule, &config).unwrap();
            mismatches += usize::from(bits != detail.fingerprint);
        }
        println!(
            "bit/detail parity chirality={include_chirality}: {}/{}",
            verification_molecules.len() - mismatches,
            verification_molecules.len()
        );
        for (name, detailed) in [("bits", false), ("detail", true)] {
            let molecules: Vec<_> = smiles
                .iter()
                .filter_map(|smiles| chematic_smiles::parse(smiles).ok())
                .collect();
            let start = Instant::now();
            for molecule in &molecules {
                if detailed {
                    std::hint::black_box(rdkit_morgan_fingerprint(molecule, &config).unwrap());
                } else {
                    std::hint::black_box(rdkit_morgan_bitvec(molecule, &config).unwrap());
                }
            }
            println!(
                "{name} chirality={include_chirality}: {:.3} us/mol",
                start.elapsed().as_secs_f64() * 1e6 / molecules.len() as f64
            );
            if include_chirality && !detailed {
                for (subset, expected) in [("plain", false), ("stereo", true)] {
                    let selected: Vec<_> = molecules
                        .iter()
                        .filter(|molecule| carries_stereo(molecule) == expected)
                        .collect();
                    let start = Instant::now();
                    for molecule in &selected {
                        std::hint::black_box(rdkit_morgan_bitvec(molecule, &config).unwrap());
                    }
                    println!(
                        "  {subset}: {} mol, {:.3} us/mol",
                        selected.len(),
                        start.elapsed().as_secs_f64() * 1e6 / selected.len().max(1) as f64
                    );
                    if expected {
                        let start = Instant::now();
                        for molecule in &selected {
                            std::hint::black_box(
                                chematic_smiles::rdkit_legacy_stereo(molecule).unwrap(),
                            );
                        }
                        println!(
                            "  legacy stereo only: {:.3} us/mol",
                            start.elapsed().as_secs_f64() * 1e6 / selected.len() as f64
                        );
                        let start = Instant::now();
                        for molecule in &selected {
                            std::hint::black_box(
                                chematic_cip::assign_cip_accurate_experimental(
                                    molecule,
                                    chematic_cip::CipBudget::default_budget(),
                                )
                                .unwrap(),
                            );
                        }
                        println!(
                            "  accurate CIP only: {:.3} us/mol",
                            start.elapsed().as_secs_f64() * 1e6 / selected.len() as f64
                        );
                    }
                }
            }
        }
    }
}
