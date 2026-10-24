use chematic_fp::{rdkit_pattern_fp, rdkit_rdk_fp};
use chematic_smiles::parse;

fn digest(bits: &chematic_fp::BitVec2048) -> u64 {
    bits.to_le_bytes()
        .iter()
        .fold(0xcbf29ce484222325, |hash, byte| {
            (hash ^ u64::from(*byte)).wrapping_mul(0x100000001b3)
        })
}

#[test]
fn aromatic_written_macrocycles_match_rdkit_2026_09_1_fingerprints() {
    // Independent Chem.PatternFingerprint and Chem.RDKFingerprint oracle.
    // RDKit kekulizes neutral O/S bridges in rings of at least nine atoms.
    for (smiles, kekule, pattern_count, pattern_digest, rdk_count, rdk_digest) in [
        (
            "o1cccccccc1",
            "C1=CC=COC=CC=C1",
            54,
            0xcd9a9ae96d04ee61,
            55,
            0x1d349e88442947c6,
        ),
        (
            "s1cccccccc1",
            "C1=CC=CSC=CC=C1",
            55,
            0xe3fd75650350ee55,
            56,
            0xa27210c4717dd2c5,
        ),
        (
            "o1cccccccccccc1",
            "C1=CC=CC=COC=CC=CC=C1",
            74,
            0x23c550b6f153e262,
            57,
            0x3a5f305cdbd73796,
        ),
        (
            "s1cccccccccccc1",
            "C1=CC=CC=CSC=CC=CC=C1",
            75,
            0x1b145b27189279ae,
            58,
            0x9de0c038ea89e215,
        ),
    ] {
        let aromatic = parse(smiles).unwrap();
        let explicit = parse(kekule).unwrap();
        let pattern = rdkit_pattern_fp(&aromatic);
        let rdk = rdkit_rdk_fp(&aromatic);
        assert_eq!(pattern, rdkit_pattern_fp(&explicit), "Pattern: {smiles}");
        assert_eq!(rdk, rdkit_rdk_fp(&explicit), "RDK: {smiles}");
        assert_eq!(
            (pattern.popcount(), digest(&pattern)),
            (pattern_count, pattern_digest),
            "Pattern: {smiles}"
        );
        assert_eq!(
            (rdk.popcount(), digest(&rdk)),
            (rdk_count, rdk_digest),
            "RDK: {smiles}"
        );
    }
}
