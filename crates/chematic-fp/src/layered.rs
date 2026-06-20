//! Layered fingerprint: RDKit-like structural layers encoding atom features.

use chematic_core::{AtomIdx, Molecule, implicit_hcount};
use chematic_perception::find_sssr;

use crate::bitvec::BitVec2048;

/// Layered fingerprint — 7 layers of atomic features packed into a 2048-bit vector.
///
/// Each layer uses a subset of the 2048 bits to encode different atom properties:
/// - Layer 0 (0–291): raw atom type (element, H count, charge)
/// - Layer 1 (292–439): heteroatom presence
/// - Layer 2 (440–587): aromaticity
/// - Layer 3 (588–735): ring membership and type
/// - Layer 4 (736–883): charge state and ionization
/// - Layer 5 (884–1031): hydrogen count distribution
/// - Layer 6 (1032–2047): custom features (connectivity, degree, fragment ID)
///
/// The fingerprint is created by hashing (atom_index, layer_bits) pairs.
/// Collisions are expected and provide dimensionality reduction.
pub fn layered_fp(mol: &Molecule) -> BitVec2048 {
    // Compute all 7 layers and merge via bitwise OR
    // (each layer uses disjoint bit ranges, so OR is lossless)
    let layers = layered_fp_by_layer(mol);
    layers
        .iter()
        .fold(BitVec2048::new(), |acc, layer| acc.or(layer))
}

/// Layer-by-layer fingerprint computation with per-layer bitsets.
///
/// Useful for examining which layers contribute to similarity.
/// Returns a vector of 7 BitVec2048 objects, one per layer.
pub fn layered_fp_by_layer(mol: &Molecule) -> [BitVec2048; 7] {
    let mut layers: [BitVec2048; 7] = [
        BitVec2048::new(),
        BitVec2048::new(),
        BitVec2048::new(),
        BitVec2048::new(),
        BitVec2048::new(),
        BitVec2048::new(),
        BitVec2048::new(),
    ];

    // Layer 0: Raw atom type
    for i in 0..mol.atom_count() {
        let atom = mol.atom(AtomIdx(i as u32));
        let elem_code = atom.element.atomic_number() as usize;
        let h_count = implicit_hcount(mol, AtomIdx(i as u32)) as usize;
        let charge = ((atom.charge.clamp(-8, 7) + 8) & 0xF) as usize; // Clamp to [-8, 7] for safe encoding
        let hash = (i * 7919 + elem_code * 199 + h_count * 73 + charge * 11) % 292;
        layers[0].set(hash);
    }

    // Layer 1: Heteroatom presence
    for i in 0..mol.atom_count() {
        let atom = mol.atom(AtomIdx(i as u32));
        let is_heteroatom = matches!(
            atom.element.atomic_number(),
            7 | 8 | 16 | 15 | 9 | 17 | 35 | 53
        );
        if is_heteroatom {
            let elem = atom.element.atomic_number() as usize;
            let hash = (292 + i * 3571 + elem * 43) % 148;
            layers[1].set(292 + hash);
        }
    }

    // Layer 2: Aromaticity
    for i in 0..mol.atom_count() {
        let atom = mol.atom(AtomIdx(i as u32));
        if atom.aromatic {
            let degree = mol.neighbors(AtomIdx(i as u32)).count();
            let hash = (440 + i * 1259 + degree * 29) % 148;
            layers[2].set(440 + hash);
        }
    }

    // Layer 3: Ring membership
    for ring in find_sssr(mol).rings() {
        for atom_idx in ring {
            let i = atom_idx.0 as usize;
            let ring_size = ring.len();
            let hash = (588 + i * 997 + ring_size * 23) % 148;
            layers[3].set(588 + hash);
        }
    }

    // Layer 4: Charge state
    for i in 0..mol.atom_count() {
        let atom = mol.atom(AtomIdx(i as u32));
        if atom.charge != 0 {
            let charge_val = ((atom.charge + 5) & 0xF) as usize;
            let is_positive = atom.charge > 0;
            let charge_type = if is_positive { 1 } else { 0 };
            let hash = (736 + i * 1447 + charge_val * 19 + charge_type * 7) % 148;
            layers[4].set(736 + hash);
        }
    }

    // Layer 5: Hydrogen count
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let h_count = implicit_hcount(mol, idx) as usize;
        if h_count > 0 {
            let hash = (884 + i * 1823 + h_count * 41) % 148;
            layers[5].set(884 + hash);
        }
    }

    // Layer 6: Connectivity and degree
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let degree = mol.neighbors(idx).count();
        let connectivity = degree + (implicit_hcount(mol, idx) as usize);
        let hash = (1032 + i * 5279 + connectivity * 31) % 1016;
        layers[6].set(1032 + hash);
    }

    layers
}

/// Tanimoto similarity between two layered fingerprints.
pub fn tanimoto_layered(a: &BitVec2048, b: &BitVec2048) -> f64 {
    a.tanimoto(b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_layered_fp_methane() {
        let mol = parse("C").unwrap();
        let fp = layered_fp(&mol);
        assert!(fp.popcount() > 0, "methane should have non-zero popcount");
    }

    #[test]
    fn test_layered_fp_benzene_aromatic() {
        let mol = parse("c1ccccc1").unwrap();
        let aromatic_layer = layered_fp_by_layer(&mol)[2].clone();
        assert!(
            aromatic_layer.popcount() > 0,
            "benzene aromatic layer should be non-zero"
        );
    }

    #[test]
    fn test_layered_fp_by_layer_count() {
        let mol = parse("CCO").unwrap();
        let layers = layered_fp_by_layer(&mol);
        assert_eq!(layers.len(), 7, "should have 7 layers");
    }

    #[test]
    fn test_layered_fp_identical_molecules() {
        let mol1 = parse("CC").unwrap();
        let mol2 = parse("CC").unwrap();
        let fp1 = layered_fp(&mol1);
        let fp2 = layered_fp(&mol2);
        let sim = tanimoto_layered(&fp1, &fp2);
        assert_eq!(sim, 1.0, "identical molecules should have tanimoto=1.0");
    }

    #[test]
    fn test_layered_fp_heteroatom_layer() {
        let mol = parse("CC").unwrap(); // no heteroatoms
        let layers = layered_fp_by_layer(&mol);
        assert_eq!(layers[1].popcount(), 0, "ethane has no heteroatoms");

        let mol2 = parse("CCO").unwrap(); // contains O
        let layers2 = layered_fp_by_layer(&mol2);
        assert!(
            layers2[1].popcount() > 0,
            "ethanol heteroatom layer should be non-zero"
        );
    }

    #[test]
    fn test_layered_fp_charged_molecule() {
        let mol = parse("[O-]C(=O)c1ccccc1").unwrap(); // benzoate anion
        let layers = layered_fp_by_layer(&mol);
        assert!(
            layers[4].popcount() > 0,
            "charged molecule should have bits in charge layer"
        );
    }

    #[test]
    fn test_layered_fp_ring_detection() {
        let mol_ring = parse("c1ccccc1").unwrap();
        let mol_linear = parse("CCCCCC").unwrap();
        let layers_ring = layered_fp_by_layer(&mol_ring);
        let layers_linear = layered_fp_by_layer(&mol_linear);
        assert!(
            layers_ring[3].popcount() > 0,
            "benzene should have ring bits"
        );
        assert_eq!(
            layers_linear[3].popcount(),
            0,
            "hexane should have no ring bits"
        );
    }

    #[test]
    fn test_layered_fp_hydrogen_count_layer() {
        let mol = parse("C").unwrap(); // methane: C with 4 H
        let layers = layered_fp_by_layer(&mol);
        assert!(
            layers[5].popcount() > 0,
            "methane hydrogen layer should be non-zero"
        );
    }

    #[test]
    fn test_layered_tanimoto_disjoint() {
        let mol1 = parse("C").unwrap();
        let mol2 = parse("c1ccccc1").unwrap();
        let fp1 = layered_fp(&mol1);
        let fp2 = layered_fp(&mol2);
        let sim = tanimoto_layered(&fp1, &fp2);
        assert!(
            (0.0..=1.0).contains(&sim),
            "tanimoto must be between 0 and 1, got {}",
            sim
        );
    }
}
