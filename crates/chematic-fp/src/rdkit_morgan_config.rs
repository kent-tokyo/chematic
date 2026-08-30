//! Promoted, config-generalized RDKit-exact Morgan/ECFP fingerprint API.
//!
//! Generalizes [`crate::rdkit_morgan_ecfp4::rdkit_morgan_ecfp4_experimental`]'s single
//! fixed point (radius=2/ECFP4, 2048 bits) to the small, closed set of `(radius,
//! fp_size)` combinations independently re-verified against a live RDKit oracle --
//! see `scripts/gen_ecfp4_rdkit_stable_api_fixtures.py` and
//! `validation/ecfp4_rdkit_stable_api_fixtures.json` (RDKit 2026.03.3, pinned commit
//! `8afba32ec539dcb2369bc84549d802aca3f7eb39`, same pin as
//! [`crate::rdkit_morgan_hash`]).
//!
//! **Deliberately not a refactor of `rdkit_morgan_ecfp4_experimental`** -- that
//! function's body is untouched; its own tests
//! (`hermetic_equivalence_to_diagnostic_default_lifecycle`,
//! `hueckel_fallback_would_be_detectable_if_silently_reintroduced`) are the guard
//! against this module accidentally changing what it computes. This module reuses the
//! same crate-internal primitives ([`crate::rdkit_morgan_hash::expand_one_pass`],
//! [`crate::rdkit_morgan_hash::checked_bond_invariant`]) but is a structurally
//! separate entry point with its own struct/config types, so a bug in the
//! generalized radius/fp_size plumbing here can never regress the frozen r=2/2048-bit
//! path. [`tests::default_config_matches_frozen_rdkit_morgan_ecfp4`] proves the two
//! paths agree at the shared (radius=2, fp_size=2048) point regardless.
//!
//! Radius and fold width are closed enums, not raw integers: every representable
//! value has an oracle-verified cell in the fixture corpus above (`radius_axis`,
//! `fp_size_axis`). There is deliberately no "pass any u32/usize, get a runtime error
//! if unsupported" path -- a caller cannot construct an unverified combination in the
//! first place, so there is nothing to guess or silently coerce.

use chematic_core::{AtomIdx, BondIdx, Molecule};
use rustc_hash::FxHashMap;

use crate::bitvec::BitVecN;
use crate::rdkit_morgan_ecfp4::RdkitMorganError;
use crate::rdkit_morgan_hash::{checked_bond_invariant, expand_one_pass_with_chirality};

/// Morgan/ECFP radius, restricted to the four values independently re-verified
/// against a live RDKit oracle for this API (`validation/ecfp4_rdkit_stable_api_fixtures.json`'s
/// per-fixture `radius_axis`). `R2` is ECFP4 -- the same radius
/// [`crate::rdkit_morgan_ecfp4::rdkit_morgan_ecfp4_experimental`] computes, and the
/// most heavily corpus-tested of the four.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RdkitMorganRadius {
    R0,
    R1,
    /// ECFP4.
    R2,
    R3,
}

impl RdkitMorganRadius {
    fn as_u32(self) -> u32 {
        match self {
            RdkitMorganRadius::R0 => 0,
            RdkitMorganRadius::R1 => 1,
            RdkitMorganRadius::R2 => 2,
            RdkitMorganRadius::R3 => 3,
        }
    }
}

/// Folded fingerprint bit width, restricted to the five values independently
/// re-verified against a live RDKit oracle (`raw_id % fp_size`, RDKit's real
/// Morgan-generator fold convention -- confirmed empirically against
/// `rdFingerprintGenerator.GetMorganGenerator(fpSize=...)`, not assumed). `B2048` is
/// the primary target and matches `rdkit_morgan_ecfp4_experimental`'s fixed size.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RdkitMorganFpSize {
    B128,
    B256,
    B512,
    B1024,
    B2048,
}

impl RdkitMorganFpSize {
    /// Bit width as a plain integer (also the fold modulus).
    pub fn bits(self) -> usize {
        match self {
            RdkitMorganFpSize::B128 => 128,
            RdkitMorganFpSize::B256 => 256,
            RdkitMorganFpSize::B512 => 512,
            RdkitMorganFpSize::B1024 => 1024,
            RdkitMorganFpSize::B2048 => 2048,
        }
    }
}

/// Configuration for [`rdkit_morgan_fingerprint`]. `Default` matches
/// `rdkit_morgan_ecfp4_experimental`'s fixed config (radius=2/ECFP4, 2048 bits).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RdkitMorganConfig {
    pub radius: RdkitMorganRadius,
    pub fp_size: RdkitMorganFpSize,
    /// Include RDKit-compatible tetrahedral chirality contributions.
    /// E/Z bond stereo remains outside this first verified increment.
    pub include_chirality: bool,
}

impl Default for RdkitMorganConfig {
    fn default() -> Self {
        RdkitMorganConfig {
            radius: RdkitMorganRadius::R2,
            fp_size: RdkitMorganFpSize::B2048,
            include_chirality: false,
        }
    }
}

/// Every RDKit-hash-exact Morgan/ECFP view of a molecule at a chosen
/// [`RdkitMorganConfig`] -- the same signals
/// [`crate::rdkit_morgan_ecfp4::RdkitMorganEcfp4`] exposes at its one fixed config,
/// generalized, plus a folded count fingerprint (RDKit's `GetCountFingerprint`
/// shape). `sparse_counts`/`raw_bit_info` depend only on `config.radius` (not
/// `fp_size` -- the fold is applied only when building `fingerprint`/
/// `folded_bit_info`/`folded_counts`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RdkitMorganFingerprint {
    /// Folded fingerprint, `config.fp_size.bits()` wide (`raw_identifier % fp_size`,
    /// OR-combined).
    pub fingerprint: BitVecN,
    /// Raw (unfolded) identifier -> emission count (RDKit's `GetSparseCountFingerprint`
    /// shape).
    pub sparse_counts: FxHashMap<u32, u32>,
    /// Raw identifier -> the `(atom_idx, radius)` environments that produced it.
    pub raw_bit_info: FxHashMap<u32, Vec<(u32, u32)>>,
    /// Folded bit -> the `(atom_idx, radius)` environments that set it.
    pub folded_bit_info: FxHashMap<usize, Vec<(u32, u32)>>,
    /// Folded bit -> emission count (RDKit's `GetCountFingerprint` shape).
    pub folded_counts: FxHashMap<usize, u32>,
}

/// RDKit-bit-exact Morgan/ECFP fingerprint at a caller-chosen [`RdkitMorganConfig`] --
/// generalizes [`crate::rdkit_morgan_ecfp4::rdkit_morgan_ecfp4_experimental`]'s fixed
/// radius=2/2048-bit point to every independently oracle-verified `(radius, fp_size)`
/// cell (see the module docs). Uses the same RDKit-parity aromaticity preprocessing,
/// with the same no-silent-fallback contract: a preprocessing failure is always an
/// `Err`, never a result computed under a different (non-bit-exact) aromaticity
/// engine.
pub fn rdkit_morgan_fingerprint(
    mol: &Molecule,
    config: &RdkitMorganConfig,
) -> Result<RdkitMorganFingerprint, RdkitMorganError> {
    let aromatized = chematic_perception::apply_aromaticity_rdkit_parity_experimental(mol)?;

    let fp_size = config.fp_size.bits();
    let mut result = RdkitMorganFingerprint {
        fingerprint: BitVecN::new(fp_size),
        sparse_counts: FxHashMap::default(),
        raw_bit_info: FxHashMap::default(),
        folded_bit_info: FxHashMap::default(),
        folded_counts: FxHashMap::default(),
    };
    if aromatized.atom_count() == 0 {
        return Ok(result);
    }

    let ring_set = chematic_perception::find_sssr(&aromatized);
    let bond_count = aromatized.bond_count();
    let mut bond_invariants = Vec::with_capacity(bond_count);
    for b in 0..bond_count {
        let bond_idx = BondIdx(b as u32);
        let order = aromatized.bond(bond_idx).order;
        let invariant = checked_bond_invariant(order)
            .ok_or(RdkitMorganError::UnsupportedBondOrder { bond_idx, order })?;
        bond_invariants.push(invariant);
    }

    let cip_codes = if config.include_chirality {
        let assignment = chematic_cip::assign_cip_accurate_experimental(
            &aromatized,
            chematic_cip::CipBudget::default_budget(),
        )
        .map_err(|e| RdkitMorganError::InternalInvariantViolation {
            reason: format!("CIP assignment failed for chiral Morgan fingerprint: {e}"),
        })?;
        Some(
            assignment
                .assignments
                .into_iter()
                .filter(|(_, code)| {
                    matches!(code, chematic_core::CipCode::R | chematic_core::CipCode::S)
                })
                .collect::<rustc_hash::FxHashMap<AtomIdx, chematic_core::CipCode>>(),
        )
    } else {
        None
    };
    let emitted = expand_one_pass_with_chirality(
        &aromatized,
        &ring_set,
        &bond_invariants,
        config.radius.as_u32(),
        true,
        cip_codes.as_ref(),
    );

    for ((atom_idx, radius), raw_id) in emitted {
        let folded = (raw_id as usize) % fp_size;
        result.fingerprint.set(folded);
        *result.sparse_counts.entry(raw_id).or_insert(0) += 1;
        result
            .raw_bit_info
            .entry(raw_id)
            .or_default()
            .push((atom_idx, radius));
        result
            .folded_bit_info
            .entry(folded)
            .or_default()
            .push((atom_idx, radius));
        *result.folded_counts.entry(folded).or_insert(0) += 1;
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::MoleculeBuilder;
    use chematic_smiles::parse;

    /// Proves this module's generalized path agrees with the frozen, untouched
    /// `rdkit_morgan_ecfp4_experimental` at their one shared config point --
    /// promotion didn't fork the math into two subtly-different implementations.
    #[test]
    fn default_config_matches_frozen_rdkit_morgan_ecfp4() {
        for smi in [
            "c1ccncc1",
            "CC(=O)[O-]",
            "c1ccc2ccccc2c1",
            "CC",
            "[13CH4]",
            "C[C@H](N)C(=O)O",
        ] {
            let mol = parse(smi).unwrap();
            let fixed = crate::rdkit_morgan_ecfp4::rdkit_morgan_ecfp4_experimental(&mol).unwrap();
            let general = rdkit_morgan_fingerprint(&mol, &RdkitMorganConfig::default()).unwrap();

            assert_eq!(
                fixed.fingerprint.to_bitvecn(),
                general.fingerprint,
                "fingerprint mismatch for {smi}"
            );
            assert_eq!(
                fixed.sparse_counts, general.sparse_counts,
                "sparse_counts mismatch for {smi}"
            );
            assert_eq!(
                fixed.raw_bit_info, general.raw_bit_info,
                "raw_bit_info mismatch for {smi}"
            );
            assert_eq!(
                fixed.folded_bit_info, general.folded_bit_info,
                "folded_bit_info mismatch for {smi}"
            );
        }
    }

    #[test]
    fn tetrahedral_chirality_matches_rdkit_and_distinguishes_enantiomers() {
        let r = parse("C[C@H](N)C(=O)O").unwrap();
        let s = parse("C[C@@H](N)C(=O)O").unwrap();
        let config = RdkitMorganConfig {
            radius: RdkitMorganRadius::R2,
            fp_size: RdkitMorganFpSize::B2048,
            include_chirality: true,
        };
        let r_result = rdkit_morgan_fingerprint(&r, &config).unwrap();
        let s_result = rdkit_morgan_fingerprint(&s, &config).unwrap();

        let expected_r: FxHashMap<u32, u32> = [
            (101282979, 1),
            (847957139, 1),
            (864662311, 1),
            (864942730, 1),
            (1510328189, 1),
            (1533864325, 1),
            (2245273601, 1),
            (2246699815, 1),
            (2246728737, 1),
            (2599973650, 1),
            (3374146648, 1),
            (3537119515, 1),
            (3855312692, 1),
        ]
        .into_iter()
        .collect();
        let expected_s: FxHashMap<u32, u32> = [
            (803825710, 1),
            (847957139, 1),
            (864662311, 1),
            (864942730, 1),
            (1510328189, 1),
            (1533864325, 1),
            (2245273601, 1),
            (2246699815, 1),
            (2246728737, 1),
            (2599973650, 1),
            (3374146649, 1),
            (3537119515, 1),
            (3855312692, 1),
        ]
        .into_iter()
        .collect();
        assert_eq!(r_result.sparse_counts, expected_r);
        assert_eq!(s_result.sparse_counts, expected_s);
        assert_ne!(r_result.sparse_counts, s_result.sparse_counts);

        let non_chiral = RdkitMorganConfig {
            include_chirality: false,
            ..config
        };
        assert_eq!(
            rdkit_morgan_fingerprint(&r, &non_chiral)
                .unwrap()
                .sparse_counts,
            rdkit_morgan_fingerprint(&s, &non_chiral)
                .unwrap()
                .sparse_counts
        );
    }

    #[test]
    fn radius_zero_has_no_bond_environment_regardless_of_fp_size() {
        // Radius 0 is a pure atom invariant -- every fp_size must fold the same raw
        // ids, just into a narrower/wider bit space.
        let mol = parse("c1ccccc1").unwrap();
        for fp_size in [
            RdkitMorganFpSize::B128,
            RdkitMorganFpSize::B256,
            RdkitMorganFpSize::B2048,
        ] {
            let result = rdkit_morgan_fingerprint(
                &mol,
                &RdkitMorganConfig {
                    radius: RdkitMorganRadius::R0,
                    fp_size,
                    include_chirality: false,
                },
            )
            .unwrap();
            assert_eq!(result.sparse_counts.len(), 1, "benzene has one radius-0 id");
            let (&raw_id, &count) = result.sparse_counts.iter().next().unwrap();
            assert_eq!(count, 6, "all 6 aromatic carbons share one radius-0 id");
            assert!(result.fingerprint.get((raw_id as usize) % fp_size.bits()));
        }
    }

    #[test]
    fn fold_is_plain_modulo_not_xor_fold() {
        // The RDKit real fold convention is `raw_id % fp_size`, NOT BitVec2048::fold's
        // XOR-fold-from-2048 -- verified empirically in the diagnosis RFC (see module
        // docs). Confirm every folded-on bit is directly explained by a raw id modulo
        // fp_size, for a molecule with more than one raw id.
        let mol = parse("CC(=O)[O-]").unwrap();
        let result = rdkit_morgan_fingerprint(
            &mol,
            &RdkitMorganConfig {
                radius: RdkitMorganRadius::R2,
                fp_size: RdkitMorganFpSize::B128,
                include_chirality: false,
            },
        )
        .unwrap();
        for &raw_id in result.sparse_counts.keys() {
            assert!(result.fingerprint.get((raw_id as usize) % 128));
        }
    }

    #[test]
    fn preprocessing_failure_is_err_for_every_config() {
        let smi = "Cc1cn2c(=O)c3ncn(COCCO)c3nc2n1C";
        let mol = parse(smi).unwrap();
        for radius in [
            RdkitMorganRadius::R0,
            RdkitMorganRadius::R1,
            RdkitMorganRadius::R2,
            RdkitMorganRadius::R3,
        ] {
            let result = rdkit_morgan_fingerprint(
                &mol,
                &RdkitMorganConfig {
                    radius,
                    fp_size: RdkitMorganFpSize::B2048,
                    include_chirality: false,
                },
            );
            assert!(
                matches!(result, Err(RdkitMorganError::Aromaticity(_))),
                "expected Aromaticity error at radius {radius:?}, got {result:?}"
            );
        }
    }

    #[test]
    fn empty_molecule_yields_empty_result_not_an_error() {
        let mol = MoleculeBuilder::new().build();
        let result = rdkit_morgan_fingerprint(&mol, &RdkitMorganConfig::default()).unwrap();
        assert_eq!(result.fingerprint.popcount(), 0);
        assert!(result.sparse_counts.is_empty());
    }
}
