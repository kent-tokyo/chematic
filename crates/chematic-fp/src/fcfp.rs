//! FCFP (Feature-Class Fingerprints) — pharmacophore-based circular fingerprints.
//!
//! Identical to ECFP in the Morgan-expansion step, but the initial atom identifier
//! is derived from pharmacophore *feature classes* rather than atomic properties.
//! This makes FCFP invariant to specific element changes (bioisostere-aware) while
//! remaining sensitive to H-bonding, charge, and hydrophobicity patterns.
//!
//! ## Feature classes
//! Each atom is assigned a bitmask encoding up to 6 binary features:
//!
//! | Bit | Class        | Definition |
//! |-----|-------------|------------|
//! | 0   | Donor       | N or O with ≥ 1 implicit H |
//! | 1   | Acceptor    | N or O (regardless of H) |
//! | 2   | Aromatic    | aromatic atom |
//! | 3   | Hydrophobic | non-aromatic C not bonded to heteroatoms; or halogens |
//! | 4   | PosIonizable| non-amide N with charge ≥ 0 |
//! | 5   | NegIonizable| O with charge ≤ 0 |

use chematic_core::{AtomIdx, BondOrder, Molecule, implicit_hcount};
use chematic_perception::find_sssr;

use crate::bitvec::BitVec2048;
use crate::ecfp::{bond_type_int, fnv1a, EcfpConfig};

// Feature class bit positions
const DONOR: u8 = 1 << 0;
const ACCEPTOR: u8 = 1 << 1;
const AROMATIC: u8 = 1 << 2;
const HYDROPHOBIC: u8 = 1 << 3;
const POS_IONIZABLE: u8 = 1 << 4;
const NEG_IONIZABLE: u8 = 1 << 5;

/// Compute the pharmacophore feature class bitmask for `idx`.
fn feature_class(mol: &Molecule, idx: AtomIdx) -> u8 {
    let atom = mol.atom(idx);
    let an = atom.element.atomic_number();
    let h = implicit_hcount(mol, idx);
    let mut cls = 0u8;

    // Aromatic
    if atom.aromatic {
        cls |= AROMATIC;
    }

    // Donor: N-H or O-H
    if (an == 7 || an == 8) && h > 0 {
        cls |= DONOR;
    }

    // Acceptor: any N or O
    if an == 7 || an == 8 {
        cls |= ACCEPTOR;
    }

    // Hydrophobic: non-aromatic C not bonded to heteroatoms; or halogens (F, Cl, Br, I)
    if an == 6 && !atom.aromatic {
        let bonded_to_hetero = mol.neighbors(idx).any(|(nb, _)| {
            let nban = mol.atom(nb).element.atomic_number();
            nban != 6 && nban != 1
        });
        if !bonded_to_hetero {
            cls |= HYDROPHOBIC;
        }
    }
    if matches!(an, 9 | 17 | 35 | 53) {
        cls |= HYDROPHOBIC;
    }

    // PosIonizable: non-aromatic N with non-negative charge that is not an amide N
    if an == 7 && !atom.aromatic {
        if atom.charge > 0 {
            cls |= POS_IONIZABLE;
        } else {
            // Basic amine: N not attached to a carbonyl carbon
            let is_amide_n = mol.neighbors(idx).any(|(nb, _)| {
                if mol.atom(nb).element.atomic_number() == 6 {
                    mol.neighbors(nb).any(|(nb2, nb2bidx)| {
                        nb2 != idx
                            && mol.atom(nb2).element.atomic_number() == 8
                            && mol.bond(nb2bidx).order == BondOrder::Double
                    })
                } else {
                    false
                }
            });
            if !is_amide_n {
                cls |= POS_IONIZABLE;
            }
        }
    }

    // NegIonizable: O with non-positive charge
    if an == 8 && atom.charge <= 0 {
        cls |= NEG_IONIZABLE;
    }

    cls
}

// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------

/// Compute an FCFP fingerprint for `mol` using the given configuration.
pub fn fcfp(mol: &Molecule, config: &EcfpConfig) -> BitVec2048 {
    let n = mol.atom_count();
    let nbits = config.nbits;
    let mut fp = BitVec2048::new();

    if n == 0 {
        return fp;
    }

    // Pre-compute ring membership (needed for degree calculation).
    let _ring_set = find_sssr(mol);

    // --- Step 1: initial atom identifiers from feature classes ---
    let mut ids: Vec<u64> = Vec::with_capacity(n);

    for i in 0..n {
        let idx = AtomIdx(i as u32);
        let cls = feature_class(mol, idx);
        let id = fnv1a(&[cls]);
        fp.set((id % nbits as u64) as usize);
        ids.push(id);
    }

    // --- Step 2: iterative Morgan expansion ---
    let mut new_ids: Vec<u64> = vec![0u64; n];

    for r in 1..=config.radius {
        for i in 0..n {
            let idx = AtomIdx(i as u32);

            let mut neighbor_info: Vec<(u8, u64)> = mol
                .neighbors(idx)
                .map(|(nb_idx, bond_idx)| {
                    let bond = mol.bond(bond_idx);
                    let btype = bond_type_int(bond.order);
                    let nb_id = ids[nb_idx.0 as usize];
                    (btype, nb_id)
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
            fp.set((new_id % nbits as u64) as usize);
        }

        core::mem::swap(&mut ids, &mut new_ids);
    }

    fp
}

/// FCFP4 fingerprint (radius = 2, 2048 bits).
pub fn fcfp4(mol: &Molecule) -> BitVec2048 {
    fcfp(mol, &EcfpConfig { radius: 2, nbits: 2048 })
}

/// FCFP6 fingerprint (radius = 3, 2048 bits).
pub fn fcfp6(mol: &Molecule) -> BitVec2048 {
    fcfp(mol, &EcfpConfig { radius: 3, nbits: 2048 })
}

/// Tanimoto similarity between two molecules using FCFP4.
pub fn tanimoto_fcfp4(a: &Molecule, b: &Molecule) -> f64 {
    fcfp4(a).tanimoto(&fcfp4(b))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap()
    }

    #[test]
    fn test_fcfp4_same_molecule() {
        let m = mol("c1ccccc1");
        assert!((fcfp4(&m).tanimoto(&fcfp4(&m)) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_fcfp4_different_molecules() {
        let benzene = mol("c1ccccc1");
        let aspirin = mol("CC(=O)Oc1ccccc1C(=O)O");
        let sim = tanimoto_fcfp4(&benzene, &aspirin);
        assert!(sim < 1.0, "FCFP4 of benzene vs aspirin should differ");
        assert!(sim > 0.0, "FCFP4 of benzene vs aspirin should share some bits");
    }

    #[test]
    fn test_fcfp4_bioisostere_closer_than_ecfp4() {
        // Benzene vs pyridine: FCFP4 should be more similar than ECFP4
        // because both ring atoms map to the same Aromatic/Acceptor class.
        use crate::ecfp::tanimoto_ecfp4;
        let benzene = mol("c1ccccc1");
        let pyridine = mol("c1ccncc1");
        let fcfp_sim = tanimoto_fcfp4(&benzene, &pyridine);
        let ecfp_sim = tanimoto_ecfp4(&benzene, &pyridine);
        assert!(fcfp_sim >= ecfp_sim,
            "FCFP4({fcfp_sim:.3}) should be ≥ ECFP4({ecfp_sim:.3}) for bioisosteres");
    }

    #[test]
    fn test_fcfp4_nonempty() {
        let fp = fcfp4(&mol("CCO"));
        assert!(fp.popcount() > 0, "FCFP4 should set at least one bit");
    }
}
