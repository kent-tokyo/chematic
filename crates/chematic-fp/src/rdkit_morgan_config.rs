//! Promoted, config-generalized RDKit-exact Morgan/ECFP fingerprint API.
//!
//! Generalizes [`crate::rdkit_morgan_ecfp4::rdkit_morgan_ecfp4_experimental`]'s single
//! fixed point (radius=2/ECFP4, 2048 bits) to the small, closed set of `(radius,
//! fp_size)` combinations independently re-verified against a live RDKit oracle --
//! see `scripts/generators/gen_ecfp4_rdkit_stable_api_fixtures.py` and
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

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule};
use rustc_hash::FxHashMap;

use crate::bitvec::BitVecN;
use crate::rdkit_morgan_ecfp4::{
    RdkitMorganError, reject_known_rdkit_coordination_sanitization_gap,
};
use crate::rdkit_morgan_hash::{
    MorganChirality, checked_bond_invariant, expand_emissions_small_with_chirality,
    expand_one_pass_with_chirality_into,
};

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
    /// RDKit's `includeChirality`: tetrahedral centres (RDKit's legacy
    /// stereo perception and `_CIPCode`) and E/Z double-bond invariants.
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

fn with_rdkit_morgan_input<R>(
    mol: &Molecule,
    config: &RdkitMorganConfig,
    consume: impl FnOnce(&Molecule, &[bool], &[u32], Option<&MorganChirality>) -> R,
) -> Result<R, RdkitMorganError> {
    if let Some(h) = chematic_smiles::rdkit_hydrogen_suppressed(mol) {
        return with_rdkit_morgan_input(&h, config, consume);
    }
    reject_known_rdkit_coordination_sanitization_gap(mol)?;
    // The identity proof preserves every atom and bond field. Borrowing that
    // view avoids an otherwise identical molecule clone and retains its caches.
    chematic_perception::with_rdkit_parity_view(mol, |view| {
        let aromatized = view.map_err(|e| RdkitMorganError::from(e.clone()))?;
        with_rdkit_morgan_aromatized_input(mol, aromatized, config, consume)
    })
}

fn with_rdkit_morgan_aromatized_input<R>(
    mol: &Molecule,
    aromatized: &Molecule,
    config: &RdkitMorganConfig,
    consume: impl FnOnce(&Molecule, &[bool], &[u32], Option<&MorganChirality>) -> R,
) -> Result<R, RdkitMorganError> {
    let ring_atoms = if aromatized.smiles_ring_closure_bond_count() == 1 {
        chematic_perception::ring_atom_flags(aromatized)
    } else {
        let ring_bonds = chematic_perception::ring_bond_flags_shared(aromatized);
        if ring_bonds.iter().any(|&is_ring| is_ring) {
            let mut atoms = vec![false; aromatized.atom_count()];
            for (bond_idx, bond) in aromatized.bonds() {
                if ring_bonds[bond_idx.0 as usize] {
                    atoms[bond.atom1.0 as usize] = true;
                    atoms[bond.atom2.0 as usize] = true;
                }
            }
            atoms
        } else {
            // `connectivity_invariant` treats a missing entry as false, so
            // an empty slice avoids one allocation for acyclic molecules.
            Vec::new()
        }
    };
    let mut bond_invariants = Vec::with_capacity(aromatized.bond_count());
    for b in 0..aromatized.bond_count() {
        let bond_idx = BondIdx(b as u32);
        let order = aromatized.bond(bond_idx).order;
        let invariant = checked_bond_invariant(order)
            .ok_or(RdkitMorganError::UnsupportedBondOrder { bond_idx, order })?;
        bond_invariants.push(invariant);
    }

    // Stereo preparation depends on the molecule, not radius or fold width.
    // Keep it on the original molecule; all graph/stereo edits clear the
    // derived cache. Holding the Arc lets expansion borrow the vectors.
    let stereo = (config.include_chirality && has_stereo_input(mol)).then(|| {
        mol.derived(chematic_core::DerivedSlot::RdkitMorganStereo, || {
            prepare_morgan_stereo(mol, aromatized)
        })
    });
    let chirality = if let Some(stereo) = &stereo {
        let prepared = stereo.as_ref().as_ref().map_err(Clone::clone)?;
        for (b, &st) in prepared.bond_stereo.iter().enumerate() {
            if st != 0 && aromatized.bond(BondIdx(b as u32)).order == BondOrder::Double {
                bond_invariants[b] = 100 + 10 * 2 + u32::from(st);
            }
        }
        Some(&prepared.chirality)
    } else {
        None
    };

    Ok(consume(
        aromatized,
        &ring_atoms,
        &bond_invariants,
        chirality,
    ))
}

struct PreparedMorganStereo {
    chirality: MorganChirality,
    bond_stereo: Vec<u8>,
}

fn prepare_morgan_stereo(
    mol: &Molecule,
    aromatized: &Molecule,
) -> Result<PreparedMorganStereo, RdkitMorganError> {
    // Keep the existing legacy-perception / accurate-CIP fallback boundary.
    match chematic_smiles::rdkit_legacy_stereo(mol) {
        Ok(stereo) if stereo.atom_tagged.len() == aromatized.atom_count() => {
            Ok(PreparedMorganStereo {
                chirality: MorganChirality {
                    tagged: stereo.atom_tagged,
                    code: stereo
                        .atom_cip
                        .iter()
                        .map(|c| match c {
                            Some(b'R') => 3,
                            Some(b'S') => 2,
                            _ => 1,
                        })
                        .collect(),
                },
                bond_stereo: stereo.bond_stereo,
            })
        }
        _ => {
            let assignment = chematic_cip::assign_cip_accurate_experimental(
                aromatized,
                chematic_cip::CipBudget::default_budget(),
            )
            .map_err(|e| RdkitMorganError::InternalInvariantViolation {
                reason: format!("CIP assignment failed for chiral Morgan fingerprint: {e}"),
            })?;
            let codes = assignment
                .assignments
                .into_iter()
                .filter(|(_, code)| {
                    matches!(code, chematic_core::CipCode::R | chematic_core::CipCode::S)
                })
                .collect::<rustc_hash::FxHashMap<AtomIdx, chematic_core::CipCode>>();
            Ok(PreparedMorganStereo {
                chirality: MorganChirality::from_cip(aromatized, &codes),
                bond_stereo: Vec::new(),
            })
        }
    }
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
    let fp_size = config.fp_size.bits();
    let mut result = RdkitMorganFingerprint {
        fingerprint: BitVecN::new(fp_size),
        sparse_counts: FxHashMap::default(),
        raw_bit_info: FxHashMap::default(),
        folded_bit_info: FxHashMap::default(),
        folded_counts: FxHashMap::default(),
    };
    with_rdkit_morgan_input(
        mol,
        config,
        |molecule, ring_atoms, bond_invariants, chirality| {
            expand_one_pass_with_chirality_into(
                molecule,
                ring_atoms,
                bond_invariants,
                config.radius.as_u32(),
                true,
                chirality,
                |atom_idx, radius, raw_id| {
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
                },
            );
        },
    )?;

    Ok(result)
}

/// RDKit-bit-exact folded Morgan bits without materializing sparse counts or
/// bit provenance. The preprocessing and emission lifecycle are shared with
/// [`rdkit_morgan_fingerprint`]; callers needing explanations must use that
/// detailed API.
pub fn rdkit_morgan_bitvec(
    mol: &Molecule,
    config: &RdkitMorganConfig,
) -> Result<BitVecN, RdkitMorganError> {
    let fp_size = config.fp_size.bits();
    let mut fingerprint = BitVecN::new(fp_size);
    with_rdkit_morgan_input(
        mol,
        config,
        |molecule, ring_atoms, bond_invariants, chirality| {
            let mut set_bit = |_: u32, _: u32, raw_id: u32| {
                fingerprint.set((raw_id as usize) % fp_size);
            };
            // Compact environments retain the general path's tetrahedral
            // hash and suppression rules, including tags resolved after R1.
            if expand_emissions_small_with_chirality(
                molecule,
                ring_atoms,
                bond_invariants,
                config.radius.as_u32(),
                chirality,
                &mut set_bit,
            )
            .is_some()
            {
                return;
            }
            expand_one_pass_with_chirality_into(
                molecule,
                ring_atoms,
                bond_invariants,
                config.radius.as_u32(),
                true,
                chirality,
                set_bit,
            );
        },
    )?;
    Ok(fingerprint)
}

/// RDKit's `countSimulation` fingerprint for a Morgan generator with
/// `config` and the default `countBounds` `[1, 2, 4, 8]`: identifiers are
/// folded to `fp_size / 4` buckets, and bucket `b` with count `c` sets bit
/// `4 * b + i` for every bound `bounds[i] <= c`.
pub fn rdkit_morgan_count_simulation(
    mol: &Molecule,
    config: &RdkitMorganConfig,
) -> Result<BitVecN, RdkitMorganError> {
    const BOUNDS: [u32; 4] = [1, 2, 4, 8];
    let result = rdkit_morgan_fingerprint(mol, config)?;
    let fp_size = config.fp_size.bits();
    let effective = (fp_size / BOUNDS.len()) as u32;
    let mut counts: FxHashMap<u32, u32> = FxHashMap::default();
    for (&raw, &count) in &result.sparse_counts {
        *counts.entry(raw % effective).or_insert(0) += count;
    }
    let mut fp = BitVecN::new(fp_size);
    for (&bucket, &count) in &counts {
        for (i, &bound) in BOUNDS.iter().enumerate() {
            if count >= bound {
                fp.set(bucket as usize * BOUNDS.len() + i);
            }
        }
    }
    Ok(fp)
}

/// Whether `mol` carries any stereo specification: a chiral tag or a
/// directional (`/`, `\\`) bond.
fn has_stereo_input(mol: &Molecule) -> bool {
    mol.atoms()
        .any(|(_, a)| a.chirality != chematic_core::Chirality::None)
        || mol
            .bonds()
            .any(|(_, b)| matches!(b.order, BondOrder::Up | BondOrder::Down))
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::MoleculeBuilder;
    use chematic_smiles::parse;

    #[derive(Debug, PartialEq, Eq)]
    struct InputSnapshot {
        atoms: Vec<chematic_core::Atom>,
        bonds: Vec<chematic_core::BondEntry>,
        ring_atoms: Vec<bool>,
        bond_invariants: Vec<u32>,
        chirality: Option<(Vec<bool>, Vec<u32>)>,
    }

    fn snapshot_input(
        mol: &Molecule,
        rings: &[bool],
        bonds: &[u32],
        chirality: Option<&MorganChirality>,
    ) -> InputSnapshot {
        InputSnapshot {
            atoms: mol.atoms().map(|(_, atom)| atom.clone()).collect(),
            bonds: mol.bonds().map(|(_, bond)| bond.clone()).collect(),
            ring_atoms: rings.to_vec(),
            bond_invariants: bonds.to_vec(),
            chirality: chirality.map(|c| (c.tagged.clone(), c.code.clone())),
        }
    }

    // Freeze the former owned-view preparation boundary. The downstream helper
    // is shared so this comparison isolates borrowing and cache ownership.
    fn owned_view_input(
        mol: &Molecule,
        config: &RdkitMorganConfig,
    ) -> Result<InputSnapshot, RdkitMorganError> {
        if let Some(h) = chematic_smiles::rdkit_hydrogen_suppressed(mol) {
            return owned_view_input(&h, config);
        }
        reject_known_rdkit_coordination_sanitization_gap(mol)?;
        let view = chematic_perception::apply_aromaticity_rdkit_parity_shared(mol);
        let aromatized = view
            .as_ref()
            .as_ref()
            .map_err(|e| RdkitMorganError::from(e.clone()))?;
        with_rdkit_morgan_aromatized_input(mol, aromatized, config, snapshot_input)
    }

    fn assert_owned_view_equivalent(mol: &Molecule) {
        for include_chirality in [false, true] {
            let config = RdkitMorganConfig {
                include_chirality,
                ..RdkitMorganConfig::default()
            };
            let reference = mol.clone();
            let expected = owned_view_input(&reference, &config);
            // Repeated calls exercise caches retained on the borrowed input.
            for _ in 0..3 {
                assert_eq!(
                    with_rdkit_morgan_input(mol, &config, snapshot_input),
                    expected,
                );
            }
            assert_eq!(
                with_rdkit_morgan_input(&mol.clone(), &config, snapshot_input),
                expected,
            );
        }
    }

    #[test]
    fn borrowed_parity_input_matches_owned_view() {
        let mut identity = 0;
        let mut nonidentity = 0;
        assert_owned_view_equivalent(&MoleculeBuilder::new().build());
        for smiles in [
            "C",
            "CC(=O)[O-]",
            "[13CH4]",
            "C1CCCCC1",
            "c1ccncc1",
            "c1ccc2ccccc2c1",
            "C1=CC=CC=C1",
            "c1ccccc1.C1=CC=CC=C1",
            "O=C1C=CC(=O)C=C1",
            "C1CC2CCC1C2",
            "C1CC2(CC1)CCCC2",
            "c1ccsc1",
            "C[C@H](N)C(=O)O",
            "C[C@@H](N)C(=O)O",
            "F/C=C\\Cl",
            "C[C@H]1CCCC[C@@H]1O",
            "C1[C@H]2CC[C@@H]1C2",
            "[2H][C@](F)(Cl)Br",
            "[H][C@](F)(Cl)Br",
            "C[C@H](N)C(=O)O.[Na+]",
            "N->[Cu+2]<-N",
            "CN(C)C[C-]12C3=C4C5=C1[Fe++]23456789[C-]%10C6=C7C8=C9%10",
            "c1cccc1",
            "[nH]1cccc1",
            "[O-][N+](=O)c1ccccc1",
        ] {
            let mol = parse(smiles).unwrap();
            if chematic_perception::rdkit_parity_view_is_identity(&mol) {
                identity += 1;
            } else {
                nonidentity += 1;
            }
            assert_owned_view_equivalent(&mol);
        }
        assert!(identity > 0 && nonidentity > 0);
    }

    #[test]
    fn borrowed_parity_input_invalidates_after_graph_and_stereo_edits() {
        let mut mol = parse("C[C@H](O)/C=C/c1ccccc1").unwrap();
        assert_owned_view_equivalent(&mol);
        mol.set_chirality(AtomIdx(1), chematic_core::Chirality::None);
        assert_owned_view_equivalent(&mol);
        mol.set_isotope(AtomIdx(0), Some(13));
        assert_owned_view_equivalent(&mol);
        mol.set_charge(AtomIdx(2), -1);
        mol.set_hydrogen_count(AtomIdx(2), Some(0));
        assert_owned_view_equivalent(&mol);
        mol.set_bond_order(BondIdx(3), BondOrder::Single);
        assert_owned_view_equivalent(&mol);
        mol.remove_bond(BondIdx(5));
        assert_owned_view_equivalent(&mol);
    }

    /// Proves this module's generalized path agrees with the frozen, untouched
    /// `rdkit_morgan_ecfp4_experimental` at their one shared config point --
    /// promotion didn't fork the math into two subtly-different implementations.
    #[test]
    fn default_config_matches_frozen_rdkit_morgan_ecfp4() {
        // `raw_bit_info` is a provenance multimap, not an ordered trace. The
        // generalized path collects through a hash map while the optimized
        // bit-oriented path emits directly, so compare its semantic contents.
        fn sorted_bit_info<K>(map: &FxHashMap<K, Vec<(u32, u32)>>) -> FxHashMap<K, Vec<(u32, u32)>>
        where
            K: Clone + Eq + std::hash::Hash,
        {
            let mut sorted = map.clone();
            for environments in sorted.values_mut() {
                environments.sort_unstable();
            }
            sorted
        }

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
                sorted_bit_info(&fixed.raw_bit_info),
                sorted_bit_info(&general.raw_bit_info),
                "raw_bit_info mismatch for {smi}"
            );
            assert_eq!(
                sorted_bit_info(&fixed.folded_bit_info),
                sorted_bit_info(&general.folded_bit_info),
                "folded_bit_info mismatch for {smi}"
            );
        }
    }

    /// The RDKit-parity view models `cleanUpOrganometallics` (dative
    /// `[C-]`→Fe), so this config path agrees with the fixed ECFP4 one.
    #[test]
    fn feii_coordination_sanitization_is_modeled() {
        let mol = parse("CN(C)C[C-]12C3=C4C5=C1[Fe++]23456789[C-]%10C6=C7C8=C9%10")
            .expect("ferrocene-like SMILES parses");
        let fp = rdkit_morgan_fingerprint(&mol, &RdkitMorganConfig::default())
            .expect("modeled sanitization");
        let ecfp4 = crate::rdkit_morgan_ecfp4::rdkit_morgan_ecfp4_experimental(&mol).unwrap();
        assert_eq!(fp.sparse_counts, ecfp4.sparse_counts);
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
    fn bit_only_path_matches_detailed_path_for_all_supported_configs() {
        for smi in [
            "c1ccncc1",
            "CC(=O)[O-]",
            "C[C@H](N)C(=O)O",
            "F/C=C\\Cl",
            "[13CH3][C@@H](O)Cl",
            "C[C@H](F)[C@@H](Cl)Br",
            "CCC[C@H](F)CCC",
            "C[C@H]1CCCC[C@@H]1O",
            "C1[C@H]2CC[C@@H]1C2",
            "C[C@H](O)/C=C/[C@@H](N)C(=O)O",
            "[2H][C@](F)(Cl)Br",
            "C[C@H](N)C(=O)O.[Na+]",
        ] {
            let mol = parse(smi).unwrap();
            for radius in [
                RdkitMorganRadius::R0,
                RdkitMorganRadius::R1,
                RdkitMorganRadius::R2,
                RdkitMorganRadius::R3,
            ] {
                for fp_size in [
                    RdkitMorganFpSize::B128,
                    RdkitMorganFpSize::B256,
                    RdkitMorganFpSize::B512,
                    RdkitMorganFpSize::B1024,
                    RdkitMorganFpSize::B2048,
                ] {
                    for include_chirality in [false, true] {
                        let config = RdkitMorganConfig {
                            radius,
                            fp_size,
                            include_chirality,
                        };
                        let detailed = rdkit_morgan_fingerprint(&mol, &config).unwrap();
                        let bits = rdkit_morgan_bitvec(&mol, &config).unwrap();
                        assert_eq!(bits, detailed.fingerprint, "{smi} {config:?}");
                    }
                }
            }
        }
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
    fn dummy_bracket_hydrogens_match_rdkit_morgan() {
        // RDKit 2026.03.6 GetMorganGenerator(radius=2, fpSize=2048).
        // A dummy atom has atomic number and mass zero, but its explicit
        // bracket H count participates in total degree and total H.
        for (smiles, raw_id, folded_bit) in [
            ("*", 2_346_609_317, 677),
            ("[*H]", 2_342_069_169, 945),
            ("[*H2]", 2_340_754_421, 1013),
        ] {
            let mol = parse(smiles).unwrap();
            let result = rdkit_morgan_fingerprint(&mol, &RdkitMorganConfig::default()).unwrap();
            assert_eq!(result.sparse_counts, [(raw_id, 1)].into_iter().collect());
            assert_eq!(result.fingerprint.popcount(), 1);
            assert!(result.fingerprint.get(folded_bit));
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
    fn valid_explicit_aromatic_input_is_accepted_for_every_config() {
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
            let result = result.unwrap_or_else(|error| {
                panic!("valid explicit aromatic input rejected at radius {radius:?}: {error:?}")
            });
            assert!(
                !result.sparse_counts.is_empty(),
                "valid explicit aromatic input should produce identifiers at radius {radius:?}"
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
