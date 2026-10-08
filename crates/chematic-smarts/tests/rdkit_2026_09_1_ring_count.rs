//! `[R<n>]` under `RdkitRingCountModel::RelevantCycles` against RDKit
//! 2026.09.1 (native C++ build, `validation/results/rdkit-native-rebaseline-
//! 2026-03-6-to-2026-09-1-2026-10-04-new.jsonl.gz`). Rows 9 and 23 are two of
//! the six bis-quinolinium macrocycles whose `[R2]`/`[R3]` cells changed from
//! 2026.03.6 and that the 2026.03.6 profile refuses; row 8818 has rings
//! through a copper atom. All 310,000 cells of that grid agree.

use chematic_smarts::{
    RdkitParityConfig, RdkitParityError, RdkitRingCountModel, find_matches_rdkit_parity,
    parse_smarts,
};

fn atoms(
    smiles: &str,
    smarts: &str,
    model: RdkitRingCountModel,
) -> Result<Vec<u32>, RdkitParityError> {
    let mol = chematic_smiles::parse(smiles).unwrap();
    let query = parse_smarts(smarts).unwrap();
    let config = RdkitParityConfig {
        use_rdkit_parity_aromaticity: true,
        ring_count_model: model,
        ..RdkitParityConfig::default()
    };
    let (matches, _) = find_matches_rdkit_parity(&query, &mol, &config)?;
    let mut out: Vec<u32> = matches
        .iter()
        .flat_map(|m| m.values().map(|a| a.0))
        .collect();
    out.sort_unstable();
    out.dedup();
    Ok(out)
}

const CASES: &[(&str, [&[u32]; 3])] = &[
    (
        "c1cc2cc(c1)-c1cccc(c1)C[n+]1ccc(c3ccccc31)NCCCCCCCCCCNc1cc[n+](c3ccccc13)C2",
        [
            &[0, 1, 5, 7, 8, 9, 18, 19, 20, 21, 40, 41, 42, 43],
            &[],
            &[14, 15, 36, 37],
        ],
    ),
    (
        "c1ccc2c(c1)c1cc[n+]2Cc2ccc(cc2)-c2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCN1",
        [&[0, 1, 2, 5, 29, 30, 31, 32], &[], &[]],
    ),
    (
        "[O-][N+]1=CC2=CC=C[O+]2[Cu]13[O+]4C=CC=C4C=[N+]3[O-]",
        [&[1, 2, 4, 5, 6, 10, 11, 12, 14, 15], &[3, 7, 8, 9, 13], &[]],
    ),
];

#[test]
fn relevant_cycle_ring_counts_match_rdkit_2026_09_1() {
    for (smiles, want) in CASES {
        for (k, query) in ["[R1]", "[R2]", "[R3]"].iter().enumerate() {
            assert_eq!(
                atoms(smiles, query, RdkitRingCountModel::RelevantCycles).unwrap(),
                want[k].to_vec(),
                "{smiles} {query}"
            );
        }
    }
}

#[test]
fn the_2026_03_profile_follows_rdkit_ring_order_on_the_macrocycles() {
    // RDKit 2026.03.1 `GetSubstructMatches(MolFromSmarts("[R2]"))` on this
    // macrocycle; the symmetrized-SSSR profile now reproduces RDKit's ring
    // order instead of refusing.
    assert_eq!(
        atoms(CASES[1].0, "[R2]", RdkitRingCountModel::SymmetrizedSssr).unwrap(),
        vec![15, 16, 21, 22, 28, 33]
    );
}
