//! ECFP (Extended Connectivity Fingerprints) based on the Morgan algorithm.
//!
//! Uses FNV-1a 64-bit hashing for reproducibility and WASM-compatibility.

use chematic_core::{AtomIdx, BondOrder, Molecule, implicit_hcount};
use chematic_perception::find_sssr;

use crate::bitvec::BitVec2048;

// ---------------------------------------------------------------------------
// FNV-1a constants
// ---------------------------------------------------------------------------

const FNV_OFFSET: u64 = 14695981039346656037;
const FNV_PRIME: u64 = 1099511628211;

/// Compute the FNV-1a 64-bit hash of `bytes`.
pub(crate) fn fnv1a(bytes: &[u8]) -> u64 {
    let mut h = FNV_OFFSET;
    for &b in bytes {
        h ^= b as u64;
        h = h.wrapping_mul(FNV_PRIME);
    }
    h
}

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

/// Configuration for ECFP computation.
#[derive(Debug, Clone)]
pub struct EcfpConfig {
    /// Number of iterations (radius). ECFP4 = 2, ECFP6 = 3.
    pub radius: u32,
    /// Output bitvector size (default 2048).
    pub nbits: usize,
}

impl Default for EcfpConfig {
    fn default() -> Self {
        Self { radius: 2, nbits: 2048 }
    }
}

// ---------------------------------------------------------------------------
// Bond type encoding
// ---------------------------------------------------------------------------

/// Map a `BondOrder` to the integer code used in the ECFP hash.
///
/// - Single / Up / Down → 1
/// - Double            → 2
/// - Triple            → 3
/// - Aromatic          → 4
/// - Quadruple         → 5  (not in standard ECFP; assigned a distinct value)
#[inline]
pub(crate) fn bond_type_int(order: BondOrder) -> u8 {
    match order {
        BondOrder::Single | BondOrder::Up | BondOrder::Down => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Aromatic => 4,
        BondOrder::Quadruple => 5,
    }
}

// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------

/// Compute an ECFP fingerprint for `mol` using the given configuration.
///
/// # Algorithm overview
/// 1. Compute initial atom identifiers from atomic properties.
/// 2. Iteratively expand each identifier by incorporating neighbour identifiers
///    (with their bond types) for `config.radius` rounds.
/// 3. After each iteration (including iteration 0), map every identifier to a
///    bit in the output bitvector.
pub fn ecfp(mol: &Molecule, config: &EcfpConfig) -> BitVec2048 {
    let n = mol.atom_count();
    let nbits = config.nbits;
    let mut fp = BitVec2048::new();

    if n == 0 {
        return fp;
    }

    // Pre-compute ring membership for all atoms once.
    let ring_set = find_sssr(mol);

    // --- Step 1: initial atom identifiers (iteration 0) ---
    let mut ids: Vec<u64> = Vec::with_capacity(n);

    for i in 0..n {
        let idx = AtomIdx(i as u32);
        let atom = mol.atom(idx);

        let atomic_number = atom.element.atomic_number();
        let degree = mol.neighbors(idx).count() as u8;
        let h_count = implicit_hcount(mol, idx);
        // Shift formal charge by 8 so that charges in [-8, +7] map to bytes [0, 15].
        // Using i16 arithmetic avoids any overflow on extreme charges.
        let charge_adjusted = (atom.charge as i16 + 8) as u8;
        let is_in_ring: u8 = if ring_set.contains_atom(idx) { 1 } else { 0 };
        let is_aromatic: u8 = if atom.aromatic { 1 } else { 0 };

        let id = fnv1a(&[
            atomic_number,
            degree,
            h_count,
            charge_adjusted,
            is_in_ring,
            is_aromatic,
        ]);

        // Set bit for this atom's identifier.
        fp.set((id % nbits as u64) as usize);
        ids.push(id);
    }

    // --- Step 2: iterative expansion ---
    let mut new_ids: Vec<u64> = vec![0u64; n];

    for r in 1..=config.radius {
        for i in 0..n {
            let idx = AtomIdx(i as u32);

            // Collect (bond_type_int, neighbor_id) pairs.
            let mut neighbor_info: Vec<(u8, u64)> = mol
                .neighbors(idx)
                .map(|(nb_idx, bond_idx)| {
                    let bond = mol.bond(bond_idx);
                    let btype = bond_type_int(bond.order);
                    let nb_id = ids[nb_idx.0 as usize];
                    (btype, nb_id)
                })
                .collect();

            // Sort for canonical ordering.
            neighbor_info.sort_unstable();

            // Build the byte array:
            //   [r as u8, id_bytes(8), (btype(1) ++ nbr_id_bytes(8)) * neighbors]
            let mut bytes: Vec<u8> = Vec::with_capacity(1 + 8 + neighbor_info.len() * 9);
            bytes.push(r as u8);
            bytes.extend_from_slice(&ids[i].to_le_bytes());
            for (btype, nb_id) in &neighbor_info {
                bytes.push(*btype);
                bytes.extend_from_slice(&nb_id.to_le_bytes());
            }

            let new_id = fnv1a(&bytes);
            new_ids[i] = new_id;

            // Set bit for the new identifier.
            fp.set((new_id % nbits as u64) as usize);
        }

        // Swap buffers; new_ids becomes ids for the next iteration.
        core::mem::swap(&mut ids, &mut new_ids);
    }

    fp
}

// ---------------------------------------------------------------------------
// Convenience wrappers
// ---------------------------------------------------------------------------

/// ECFP4 fingerprint (radius = 2, 2048 bits).
pub fn ecfp4(mol: &Molecule) -> BitVec2048 {
    ecfp(mol, &EcfpConfig { radius: 2, nbits: 2048 })
}

/// ECFP6 fingerprint (radius = 3, 2048 bits).
pub fn ecfp6(mol: &Molecule) -> BitVec2048 {
    ecfp(mol, &EcfpConfig { radius: 3, nbits: 2048 })
}

/// Tanimoto similarity between two molecules using ECFP4.
pub fn tanimoto_ecfp4(a: &Molecule, b: &Molecule) -> f64 {
    ecfp4(a).tanimoto(&ecfp4(b))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn benzene() -> Molecule {
        parse("c1ccccc1").unwrap()
    }

    fn ethane() -> Molecule {
        parse("CC").unwrap()
    }

    fn toluene() -> Molecule {
        parse("Cc1ccccc1").unwrap()
    }

    fn aspirin() -> Molecule {
        // Acetylsalicylic acid
        parse("CC(=O)Oc1ccccc1C(=O)O").unwrap()
    }

    fn methane() -> Molecule {
        parse("C").unwrap()
    }

    fn water() -> Molecule {
        parse("O").unwrap()
    }

    #[test]
    fn benzene_ecfp4_nonzero() {
        let fp = ecfp4(&benzene());
        assert!(fp.popcount() > 0, "benzene ECFP4 must be non-zero");
    }

    #[test]
    fn benzene_ecfp4_deterministic() {
        let fp1 = ecfp4(&benzene());
        let fp2 = ecfp4(&benzene());
        assert_eq!(fp1, fp2, "ECFP4 must be deterministic");
    }

    #[test]
    fn ethane_vs_benzene_tanimoto_lt1() {
        let t = tanimoto_ecfp4(&ethane(), &benzene());
        assert!(t < 1.0, "ethane and benzene must differ (tanimoto={t})");
    }

    #[test]
    fn benzene_vs_benzene_tanimoto_eq1() {
        let t = tanimoto_ecfp4(&benzene(), &benzene());
        assert_eq!(t, 1.0, "identical molecules must have tanimoto == 1.0");
    }

    #[test]
    fn benzene_vs_toluene_tanimoto_between() {
        let t = tanimoto_ecfp4(&benzene(), &toluene());
        assert!(t > 0.0, "benzene and toluene share bits (tanimoto={t})");
        assert!(t < 1.0, "benzene and toluene are not identical (tanimoto={t})");
    }

    #[test]
    fn aspirin_ecfp4_many_bits() {
        let fp = ecfp4(&aspirin());
        assert!(
            fp.popcount() > 5,
            "aspirin ECFP4 must have more than 5 bits set (got {})",
            fp.popcount()
        );
    }

    #[test]
    fn ecfp6_vs_ecfp4_benzene_differ() {
        let fp4 = ecfp4(&benzene());
        let fp6 = ecfp6(&benzene());
        // Larger radius explores more environment — the bit counts should differ
        // because radius-3 adds new hash values not present at radius-2.
        assert_ne!(
            fp4.popcount(),
            fp6.popcount(),
            "ECFP6 and ECFP4 should produce different bit counts for benzene"
        );
    }

    #[test]
    fn methane_ecfp4_nonzero() {
        let fp = ecfp4(&methane());
        assert!(fp.popcount() > 0, "methane ECFP4 must be non-zero");
    }

    #[test]
    fn water_ecfp4_nonzero() {
        let fp = ecfp4(&water());
        assert!(fp.popcount() > 0, "water ECFP4 must be non-zero");
    }

    #[test]
    fn tanimoto_ecfp4_benzene_self_is_one() {
        let t = tanimoto_ecfp4(&benzene(), &benzene());
        assert_eq!(t, 1.0, "tanimoto_ecfp4 of identical molecules must be 1.0");
    }

    #[test]
    fn tanimoto_ecfp4_methane_vs_benzene_lt_half() {
        let t = tanimoto_ecfp4(&methane(), &benzene());
        assert!(
            t < 0.5,
            "methane and benzene should be very dissimilar (tanimoto={t})"
        );
    }
}
