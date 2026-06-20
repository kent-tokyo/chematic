//! 3D pharmacophore fingerprints using Euclidean distances.
//!
//! Complements `chematic_fp::pharmacophore_fp_2d` by using 3D Euclidean distances
//! instead of topological distances for feature-pair encoding.

use crate::coords::Coords3D;
use chematic_core::Molecule;
use chematic_fp::BitVec2048;
use chematic_perception::{FeatureType, detect_features};

/// Generate a 3D pharmacophore fingerprint from molecular features and 3D coordinates.
///
/// Uses Euclidean distances in 3D space instead of topological distances.
/// Returns a 2048-bit vector encoding:
/// - Bits 0-5:   Feature type presence (Donor, Acceptor, Aromatic, Hydrophobic, Positive, Negative)
/// - Bits 6-2047: Feature-pair distance bins (6 feature types × 341 distance bins)
pub fn pharmacophore_fp_3d(mol: &Molecule, coords: &Coords3D) -> BitVec2048 {
    let mut fp = BitVec2048::new();
    let features = detect_features(mol);

    if features.is_empty() {
        return fp;
    }

    // Feature type presence bits (0-5)
    let mut types_seen = [false; 6];
    for feat in &features {
        let type_idx = feature_type_to_index(feat.ftype);
        types_seen[type_idx] = true;
    }
    for (i, &seen) in types_seen.iter().enumerate() {
        if seen {
            fp.set(i);
        }
    }

    // Feature-pair 3D distances (6-2047)
    let n = features.len();
    for i in 0..n {
        for j in (i + 1)..n {
            let ti = feature_type_to_index(features[i].ftype);
            let tj = feature_type_to_index(features[j].ftype);

            // Euclidean distance in 3D space
            let pi = coords.get(features[i].atom);
            let pj = coords.get(features[j].atom);
            let dist = pi.distance(&pj);

            let bin = distance_to_bin_3d(dist);

            // Pair-specific bit
            let pair_idx = ti * 6 + tj;
            let bit_pos = 6 + pair_idx * 341 + bin;
            if bit_pos < 2048 {
                fp.set(bit_pos);
            }
        }
    }

    fp
}

/// Tanimoto similarity between two 3D pharmacophore fingerprints.
pub fn tanimoto_pharmacophore_3d(a: &BitVec2048, b: &BitVec2048) -> f64 {
    a.tanimoto(b)
}

// Helper: map FeatureType to index (0-5)
fn feature_type_to_index(ftype: FeatureType) -> usize {
    match ftype {
        FeatureType::Donor => 0,
        FeatureType::Acceptor => 1,
        FeatureType::Aromatic => 2,
        FeatureType::Hydrophobic => 3,
        FeatureType::Positive => 4,
        FeatureType::Negative => 5,
    }
}

// Helper: map Euclidean distance (Angstroms) to bin (0-340)
fn distance_to_bin_3d(dist: f64) -> usize {
    let dist_int = (dist * 10.0) as usize; // Convert to 0.1 Å steps
    if dist_int <= 100 {
        dist_int
    } else {
        // Bin larger distances (> 10 Å)
        100 + ((dist_int - 100) / 10).min(240)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dg::generate_coords;
    use chematic_smiles::parse;

    #[test]
    fn test_pharmacophore_fp_3d_benzene() {
        let mol = parse("c1ccccc1").unwrap();
        let coords = generate_coords(&mol);
        let fp = pharmacophore_fp_3d(&mol, &coords);
        assert!(
            fp.popcount() > 0,
            "benzene should have aromatic features in 3D"
        );
    }

    #[test]
    fn test_pharmacophore_fp_3d_ethanol() {
        let mol = parse("CCO").unwrap();
        let coords = generate_coords(&mol);
        let fp = pharmacophore_fp_3d(&mol, &coords);
        // Ethanol has O atom: should have Donor and Acceptor
        assert!(
            fp.get(0) || fp.get(1),
            "ethanol should have donor/acceptor in 3D"
        );
    }

    #[test]
    fn test_tanimoto_similarity_3d_identical() {
        let mol = parse("c1ccccc1").unwrap();
        let coords = generate_coords(&mol);
        let fp = pharmacophore_fp_3d(&mol, &coords);
        let sim = tanimoto_pharmacophore_3d(&fp, &fp);
        assert!(
            (sim - 1.0).abs() < 1e-6,
            "identical fingerprints should have similarity 1.0"
        );
    }
}
