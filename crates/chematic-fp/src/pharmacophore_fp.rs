//! Pharmacophore fingerprints for virtual screening and similarity search.
//!
//! Converts pharmacophore features (detected by `chematic_perception::detect_features`)
//! into 2D fingerprints suitable for molecular property prediction and lead discovery.
//! For 3D fingerprints, see `chematic-3d::pharmacophore_fp_3d`.
//!
//! **2D Pharmacophore FP**: Feature types + feature-feature distances (2D topological)

use crate::bitvec::BitVec2048;
use chematic_core::Molecule;
use chematic_perception::{FeatureType, detect_features};

/// Generate a 2D pharmacophore fingerprint from molecular features.
///
/// Returns a 2048-bit vector encoding:
/// - Bits 0-5:   Feature type presence (Donor, Acceptor, Aromatic, Hydrophobic, Positive, Negative)
/// - Bits 6-2047: Feature-pair distance bins (6 feature types × 341 distance bins)
///
/// Distance bins for 2D:
/// - Bin 0-10:   Topological distance 1-11
/// - (higher bins reserved for future expansion)
pub fn pharmacophore_fp_2d(mol: &Molecule) -> BitVec2048 {
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

    // Feature-pair distances (6-2047)
    // Each pair (i,j) gets a bin based on topological distance.
    // For 2D, we use a simple heuristic: atom pair distance on molecule graph.
    let n = features.len();
    for i in 0..n {
        for j in (i + 1)..n {
            let ti = feature_type_to_index(features[i].ftype);
            let tj = feature_type_to_index(features[j].ftype);
            // Sort the pair so the bit depends on the pair of feature TYPES, not on
            // which atom happens to appear first in `features` (detection/parse order).
            let (t_lo, t_hi) = if ti <= tj { (ti, tj) } else { (tj, ti) };

            // Simple distance metric: atom graph distance (BFS)
            let dist = topological_distance(mol, features[i].atom, features[j].atom);
            let bin = distance_to_bin_2d(dist);

            // Pair-specific bit: 6 + (t_lo * 6 + t_hi) * 341 + bin
            let pair_idx = t_lo * 6 + t_hi;
            let bit_pos = 6 + pair_idx * 341 + bin;
            if bit_pos < 2048 {
                fp.set(bit_pos);
            }
        }
    }

    fp
}

/// Generate a simple pharmacophore feature count vector (6 elements, one per feature type).
/// Useful for lightweight screening or as a pre-filter before full FP matching.
pub fn pharmacophore_feature_counts(mol: &Molecule) -> [usize; 6] {
    let features = detect_features(mol);
    let mut counts = [0; 6];

    for feat in features {
        let idx = feature_type_to_index(feat.ftype);
        counts[idx] += 1;
    }

    counts
}

/// Tanimoto similarity between two 2D pharmacophore fingerprints.
pub fn tanimoto_pharmacophore_2d(a: &BitVec2048, b: &BitVec2048) -> f64 {
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

// Helper: map topological distance to bin (0-340)
fn distance_to_bin_2d(dist: usize) -> usize {
    if dist == 0 {
        0
    } else if dist <= 11 {
        dist - 1
    } else {
        // Bin larger distances into wider buckets
        11 + ((dist - 12) / 5).min(65)
    }
}

// Helper: compute topological distance between two atoms via BFS
fn topological_distance(
    mol: &Molecule,
    a: chematic_core::AtomIdx,
    b: chematic_core::AtomIdx,
) -> usize {
    if a == b {
        return 0;
    }

    use std::collections::VecDeque;
    let mut visited = vec![false; mol.atom_count()];
    let mut dist = vec![usize::MAX; mol.atom_count()];
    let mut queue = VecDeque::new();

    queue.push_back(a);
    visited[a.0 as usize] = true;
    dist[a.0 as usize] = 0;

    while let Some(curr) = queue.pop_front() {
        if curr == b {
            return dist[b.0 as usize];
        }

        // Iterate through all bonds to find neighbors of curr
        for (_, bond) in mol.bonds() {
            let next = if bond.atom1 == curr {
                Some(bond.atom2)
            } else if bond.atom2 == curr {
                Some(bond.atom1)
            } else {
                None
            };

            if let Some(next_idx) = next
                && !visited[next_idx.0 as usize]
            {
                visited[next_idx.0 as usize] = true;
                dist[next_idx.0 as usize] = dist[curr.0 as usize] + 1;
                queue.push_back(next_idx);
            }
        }
    }

    usize::MAX
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_pharmacophore_fp_2d_pair_bit_order_independent() {
        // Ether O (acceptor, no H) and amine N (donor) are different feature
        // types at a fixed topological distance -- real tie material for the
        // pair_idx computation, since which one is "features[i]" vs
        // "features[j]" flips between these two (same-molecule) spellings.
        let a = pharmacophore_fp_2d(&parse("COCCCCN").unwrap());
        let b = pharmacophore_fp_2d(&parse("NCCCCOC").unwrap());
        assert_eq!(
            a, b,
            "pharmacophore_fp_2d must not depend on which atom is detected first"
        );
    }

    #[test]
    fn test_pharmacophore_fp_2d_benzene() {
        let mol = parse("c1ccccc1").unwrap();
        let fp = pharmacophore_fp_2d(&mol);
        assert!(fp.popcount() > 0, "benzene should have aromatic features");
        assert!(
            fp.get(feature_type_to_index(FeatureType::Aromatic)),
            "aromatic bit should be set"
        );
    }

    #[test]
    fn test_pharmacophore_fp_2d_ethanol() {
        let mol = parse("CCO").unwrap();
        let fp = pharmacophore_fp_2d(&mol);
        // Ethanol has O atom: should have Donor and Acceptor
        assert!(fp.get(0) || fp.get(1), "ethanol should have donor/acceptor");
    }

    #[test]
    fn test_pharmacophore_feature_counts_aspirin() {
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let counts = pharmacophore_feature_counts(&mol);
        // Aspirin: 2 O (acceptors), 1 aromatic ring, 1 C=O (aromatic region)
        assert!(counts[1] > 0, "aspirin should have acceptor features");
        assert!(counts[2] > 0, "aspirin should have aromatic features");
    }

    #[test]
    fn test_tanimoto_similarity_identical() {
        let mol = parse("c1ccccc1").unwrap();
        let fp = pharmacophore_fp_2d(&mol);
        let sim = tanimoto_pharmacophore_2d(&fp, &fp);
        assert!(
            (sim - 1.0).abs() < 1e-6,
            "identical fingerprints should have similarity 1.0"
        );
    }
}
