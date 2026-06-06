//! ECFP (Extended Connectivity Fingerprints) based on the Morgan algorithm.
//!
//! Uses FNV-1a 64-bit hashing for reproducibility and WASM-compatibility.

use chematic_core::{AtomIdx, BondOrder, Molecule, implicit_hcount};
use chematic_perception::find_sssr;

use crate::bitvec::BitVec2048;

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

/// Configuration for ECFP computation.
#[derive(Debug, Clone)]
pub struct EcfpConfig {
    /// Number of iterations (radius). ECFP4 = 2, ECFP6 = 3.
    pub radius: u32,
    /// Output bitvector size (default 2048).
    pub nbits: usize,
    /// When `true`, include tetrahedral chirality in the initial atom hash so
    /// that R and S enantiomers produce different fingerprints.
    ///
    /// Defaults to `false` (chirality ignored, matching RDKit's `useChirality=False`).
    pub use_chirality: bool,
}

impl Default for EcfpConfig {
    fn default() -> Self {
        Self { radius: 2, nbits: 2048, use_chirality: false }
    }
}

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

    let ring_set = find_sssr(mol);

    // Step 1: initial atom identifiers (iteration 0).
    let mut ids: Vec<u64> = Vec::with_capacity(n);
    for i in 0..n {
        let idx = AtomIdx(i as u32);
        let atom = mol.atom(idx);

        // Shift formal charge by 8 so that charges in [-8, +7] map to bytes [0, 15].
        // i16 arithmetic avoids overflow on extreme charges.
        let charge_adjusted = (atom.charge as i16 + 8) as u8;
        let base_bytes = [
            atom.element.atomic_number(),
            mol.neighbors(idx).count() as u8,
            implicit_hcount(mol, idx),
            charge_adjusted,
            ring_set.contains_atom(idx) as u8,
            atom.aromatic as u8,
        ];
        // When use_chirality is enabled, append the chirality byte so that R
        // and S enantiomers produce different initial identifiers.  The extra
        // byte is only appended when chirality is requested — this preserves
        // bit-compatibility with the default (use_chirality=false) fingerprints.
        let id = if config.use_chirality {
            use chematic_core::Chirality;
            let chirality_byte = match atom.chirality {
                Chirality::None => 0u8,
                Chirality::CounterClockwise => 1u8,
                Chirality::Clockwise => 2u8,
            };
            let mut chiral_bytes = base_bytes.to_vec();
            chiral_bytes.push(chirality_byte);
            fnv1a(&chiral_bytes)
        } else {
            fnv1a(&base_bytes)
        };

        fp.set((id % nbits as u64) as usize);
        ids.push(id);
    }

    // Step 2: iterative expansion.
    let mut new_ids: Vec<u64> = vec![0u64; n];
    for r in 1..=config.radius {
        for i in 0..n {
            let idx = AtomIdx(i as u32);
            let mut neighbor_info: Vec<(u8, u64)> = mol
                .neighbors(idx)
                .map(|(nb_idx, bond_idx)| {
                    let btype = bond_type_int(mol.bond(bond_idx).order);
                    (btype, ids[nb_idx.0 as usize])
                })
                .collect();
            neighbor_info.sort_unstable();

            // Byte layout: [round_u8, self_id(8), (btype(1) ++ nb_id(8))*]
            let mut bytes: Vec<u8> = Vec::with_capacity(1 + 8 + neighbor_info.len() * 9);
            bytes.push(r as u8);
            bytes.extend_from_slice(&ids[i].to_le_bytes());
            for (btype, nb_id) in &neighbor_info {
                bytes.push(*btype);
                bytes.extend_from_slice(&nb_id.to_le_bytes());
            }

            let new_id = fnv1a(&bytes);
            new_ids[i] = new_id;
            fp.set((new_id % nbits as u64) as usize);
        }
        core::mem::swap(&mut ids, &mut new_ids);
    }

    fp
}

/// Count-based Morgan fingerprint: returns a map of `hash → count` for all
/// atom environments up to `radius` iterations.
///
/// Each (atom, iteration) pair contributes its hash to the map.  Unlike the
/// default RDKit behavior, redundant (duplicate) environments are **not**
/// suppressed — every atom contributes at every iteration level.
///
/// This corresponds to `GetMorganFingerprint(mol, radius,
/// useFeatures=False, includeRedundantEnvironments=True)` in RDKit.
pub fn morgan_fp_counts(mol: &Molecule, radius: u32) -> std::collections::HashMap<u64, u32> {
    use std::collections::HashMap;

    let n = mol.atom_count();
    let mut counts: HashMap<u64, u32> = HashMap::new();

    if n == 0 {
        return counts;
    }

    let ring_set = find_sssr(mol);

    // Radius-0: initial atom identifiers.
    let mut ids: Vec<u64> = (0..n)
        .map(|i| {
            let idx = AtomIdx(i as u32);
            let atom = mol.atom(idx);
            let charge_adjusted = (atom.charge as i16 + 8) as u8;
            fnv1a(&[
                atom.element.atomic_number(),
                mol.neighbors(idx).count() as u8,
                implicit_hcount(mol, idx),
                charge_adjusted,
                ring_set.contains_atom(idx) as u8,
                atom.aromatic as u8,
            ])
        })
        .collect();

    for &id in &ids {
        *counts.entry(id).or_insert(0) += 1;
    }

    // Radius 1..=radius: iterative expansion (same hash scheme as ecfp).
    let mut new_ids = vec![0u64; n];
    for r in 1..=radius {
        for i in 0..n {
            let idx = AtomIdx(i as u32);
            let mut neighbor_info: Vec<(u8, u64)> = mol
                .neighbors(idx)
                .map(|(nb_idx, bond_idx)| {
                    (bond_type_int(mol.bond(bond_idx).order), ids[nb_idx.0 as usize])
                })
                .collect();
            neighbor_info.sort_unstable();

            let mut bytes: Vec<u8> = Vec::with_capacity(1 + 8 + neighbor_info.len() * 9);
            bytes.push(r as u8);
            bytes.extend_from_slice(&ids[i].to_le_bytes());
            for (btype, nb_id) in &neighbor_info {
                bytes.push(*btype);
                bytes.extend_from_slice(&nb_id.to_le_bytes());
            }

            let new_id = fnv1a(&bytes);
            new_ids[i] = new_id;
            *counts.entry(new_id).or_insert(0) += 1;
        }
        core::mem::swap(&mut ids, &mut new_ids);
    }

    counts
}

/// ECFP4 fingerprint (radius = 2, 2048 bits).
pub fn ecfp4(mol: &Molecule) -> BitVec2048 {
    ecfp(mol, &EcfpConfig::default())
}

/// ECFP6 fingerprint (radius = 3, 2048 bits).
pub fn ecfp6(mol: &Molecule) -> BitVec2048 {
    ecfp(mol, &EcfpConfig { radius: 3, nbits: 2048, use_chirality: false })
}

/// Tanimoto similarity between two molecules using ECFP4.
pub fn tanimoto_ecfp4(a: &Molecule, b: &Molecule) -> f64 {
    ecfp4(a).tanimoto(&ecfp4(b))
}

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

    // ── morgan_fp_counts ─────────────────────────────────────────────────────

    #[test]
    fn morgan_counts_radius0_atom_count() {
        // At radius 0, one hash per atom.
        let m = benzene();
        let counts = morgan_fp_counts(&m, 0);
        let total: u32 = counts.values().sum();
        assert_eq!(total, m.atom_count() as u32, "radius-0 total count == atom_count");
    }

    #[test]
    fn morgan_counts_radius2_total_grows() {
        // Each additional radius adds one hash per atom → total = n * (radius+1).
        let m = methane();
        let n = m.atom_count() as u32;
        let r = 2u32;
        let counts = morgan_fp_counts(&m, r);
        let total: u32 = counts.values().sum();
        assert_eq!(total, n * (r + 1), "methane total = atom_count * (radius+1)");
    }

    #[test]
    fn morgan_counts_benzene_symmetry() {
        // All 6 benzene C atoms are equivalent → radius-0 yields 1 unique hash.
        let m = benzene();
        let counts = morgan_fp_counts(&m, 0);
        assert_eq!(counts.len(), 1, "benzene has one unique radius-0 environment");
        assert_eq!(*counts.values().next().unwrap(), 6, "that environment appears 6 times");
    }

    #[test]
    fn morgan_counts_empty_mol_is_empty() {
        use chematic_core::MoleculeBuilder;
        let m = MoleculeBuilder::new().build();
        let counts = morgan_fp_counts(&m, 2);
        assert!(counts.is_empty(), "empty molecule yields empty count map");
    }

    #[test]
    fn morgan_counts_deterministic() {
        let m = aspirin();
        let c1 = morgan_fp_counts(&m, 2);
        let c2 = morgan_fp_counts(&m, 2);
        assert_eq!(c1, c2, "morgan_fp_counts must be deterministic");
    }

    #[test]
    fn morgan_counts_consistent_with_ecfp_bits() {
        // Every hash in the count map should be reachable from the ecfp bit set
        // (after folding to 2048 bits).  This checks the same hash scheme.
        let m = toluene();
        let fp = ecfp(&m, &EcfpConfig { radius: 2, nbits: 2048, use_chirality: false });
        let counts = morgan_fp_counts(&m, 2);
        for &hash in counts.keys() {
            let bit = (hash % 2048) as usize;
            assert!(fp.get(bit), "bit {bit} from count map not set in ECFP bitvec");
        }
    }

    // -- Chirality tests ------------------------------------------------------

    #[test]
    fn ecfp4_ignores_chirality_by_default() {
        // L-alanine and D-alanine should produce the same ECFP4 when
        // use_chirality=false (default), since chirality is not in the hash.
        let l_ala = parse("N[C@@H](C)C(=O)O").unwrap();
        let d_ala = parse("N[C@H](C)C(=O)O").unwrap();
        let fp_l = ecfp4(&l_ala);
        let fp_d = ecfp4(&d_ala);
        assert_eq!(fp_l, fp_d, "L/D-alanine ECFP4 should be identical when use_chirality=false");
    }

    #[test]
    fn ecfp4_distinguishes_enantiomers_with_chirality() {
        // With use_chirality=true, L-alanine and D-alanine must have different FPs.
        let l_ala = parse("N[C@@H](C)C(=O)O").unwrap();
        let d_ala = parse("N[C@H](C)C(=O)O").unwrap();
        let config = EcfpConfig { radius: 2, nbits: 2048, use_chirality: true };
        let fp_l = ecfp(&l_ala, &config);
        let fp_d = ecfp(&d_ala, &config);
        assert_ne!(fp_l, fp_d, "L/D-alanine ECFP4 must differ when use_chirality=true");
        // Tanimoto < 1.0 confirms they are not identical.
        assert!(fp_l.tanimoto(&fp_d) < 1.0, "Tanimoto of L/D-alanine must be < 1.0 with use_chirality");
    }
}
