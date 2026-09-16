//! Production, fallible, RDKit-bit-exact ECFP4 path.
//!
//! Promotes [`crate::rdkit_morgan_hash`]'s source-verified hash port (Milestone M4-A0, PR #124)
//! to a real public API, restricted to the single option envelope M4-A0 actually verified
//! numerically against a live RDKit oracle: radius = 2 (ECFP4), fpSize = 2048,
//! `includeRedundantEnvironments = false` (RDKit's default, suppressed lifecycle),
//! `useChirality = false`, `useBondTypes = true`. **Not** ECFP6/radius = 3 — M4-A0 never
//! compared radius 3 against the oracle, so claiming bit-exactness there would be an
//! unverified extrapolation.
//!
//! Uses [`apply_aromaticity_rdkit_parity_experimental`] internally as a fallible `Result`
//! step. There is no fallback to production Hückel aromaticity anywhere in this module's
//! public path — see the project's own
//! `feedback_fallback_pooling_measurement_error` lesson (M4-A0's original report pooled a
//! Hückel-fallback result into an "RDKit-parity" success count and had to be corrected): this
//! module's whole claim is bit-exactness against real RDKit, so silently substituting a
//! different aromaticity engine on `Err` would silently invalidate that claim on exactly the
//! inputs where a caller most needs to know it doesn't hold.
//!
//! The parity engine may preserve a parser-supplied, internally consistent explicit-aromatic
//! representation when chematic's matching-based Kekulé conversion cannot represent it. This
//! is not a Hückel fallback: no aromaticity is inferred, and self-inconsistent input remains a
//! hard error.

use chematic_core::{AtomIdx, BondIdx, BondOrder, Element, Molecule};
use chematic_perception::AromaticityError;
use rustc_hash::FxHashMap;

use crate::bitvec::BitVec2048;
use crate::rdkit_morgan_hash::{checked_bond_invariant, expand_one_pass};

const ECFP4_RADIUS: u32 = 2;
const ECFP4_FP_SIZE: usize = 2048;
type HypervalentHalogenOxoacidNormalization = (AtomIdx, i8, Vec<(BondIdx, AtomIdx)>);

/// Every RDKit-hash-exact ECFP4 view of a molecule, computed from one shared expansion pass
/// (RDKit's `includeRedundantEnvironments = false` lifecycle) — not independently recomputed
/// per field.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RdkitMorganEcfp4 {
    /// 2048-bit folded fingerprint (`raw_identifier % 2048`, OR-combined).
    pub fingerprint: BitVec2048,
    /// Raw (unfolded) identifier → emission count (RDKit's `GetSparseCountFingerprint`
    /// shape). An identifier's count is how many distinct `(atom, radius)` environments
    /// emitted it, which can exceed 1 on an accidental hash collision.
    pub sparse_counts: FxHashMap<u32, u32>,
    /// Raw identifier → the `(atom_idx, radius)` environments that produced it (RDKit's
    /// `AdditionalOutput.GetBitInfoMap()`, unfolded).
    pub raw_bit_info: FxHashMap<u32, Vec<(u32, u32)>>,
    /// Folded bit (`0..2048`) → the `(atom_idx, radius)` environments that set it (RDKit's
    /// `AdditionalOutput.GetBitInfoMap()` on the folded fingerprint).
    pub folded_bit_info: FxHashMap<usize, Vec<(u32, u32)>>,
}

/// Why [`rdkit_morgan_ecfp4_experimental`] could not produce a result.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RdkitMorganError {
    /// RDKit-parity aromaticity preprocessing failed — see the wrapped
    /// [`AromaticityError`] for the specific reason (kekulization failure or an internal
    /// invariant violation). No fallback to another aromaticity engine is ever attempted; see
    /// the module docs.
    Aromaticity(AromaticityError),
    /// `bond_idx`'s `order` has no real RDKit `Bond::BondType` counterpart (only chematic's
    /// SMARTS-query-only `BondOrder` variants — cannot occur for a SMILES-parsed molecule, but
    /// is checked explicitly rather than assumed unreachable).
    UnsupportedBondOrder { bond_idx: BondIdx, order: BondOrder },
    /// The measured Fe(II) coordination graph whose RDKit sanitization rewrites
    /// covalent input into directed coordination bonds. The compatibility
    /// profile does not model that rewrite, so returning a numeric fingerprint
    /// would falsely imply bit exactness.
    UnsupportedCoordinationSanitization {
        atom_idx: AtomIdx,
        atomic_number: u8,
        degree: usize,
    },
    /// A post-computation sanity check failed that should never happen for chemically valid
    /// input — surfaced as an error rather than a panic or a silently wrong fingerprint.
    InternalInvariantViolation { reason: String },
}

impl std::fmt::Display for RdkitMorganError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RdkitMorganError::Aromaticity(e) => write!(f, "rdkit-exact ecfp4: aromaticity: {e}"),
            RdkitMorganError::UnsupportedBondOrder { bond_idx, order } => write!(
                f,
                "rdkit-exact ecfp4: bond {bond_idx:?} has no RDKit BondType counterpart: {order:?}"
            ),
            RdkitMorganError::UnsupportedCoordinationSanitization {
                atom_idx,
                atomic_number,
                degree,
            } => write!(
                f,
                "rdkit-exact ecfp4: unsupported RDKit coordination sanitization at {atom_idx:?} (Z={atomic_number}, degree={degree})"
            ),
            RdkitMorganError::InternalInvariantViolation { reason } => {
                write!(
                    f,
                    "rdkit-exact ecfp4: internal invariant violation: {reason}"
                )
            }
        }
    }
}

impl std::error::Error for RdkitMorganError {}

impl From<AromaticityError> for RdkitMorganError {
    fn from(e: AromaticityError) -> Self {
        RdkitMorganError::Aromaticity(e)
    }
}

/// RDKit-bit-exact ECFP4 (radius = 2, 2048 bits, `includeRedundantEnvironments = false`,
/// `useChirality = false`, `useBondTypes = true`) — bit-exact against real RDKit for every
/// input where preprocessing succeeds, verified on the full M4-A0 corpus (5,046/5,046
/// `rdkit_parity_success` rows, 100% agreement across raw identifiers, sparse counts, folded
/// bits, and bitInfo). Kekulization/aromaticity failures are reported as `Err`, never silently
/// degraded to a different, non-exact result — see the module docs.
pub fn rdkit_morgan_ecfp4_experimental(
    mol: &Molecule,
) -> Result<RdkitMorganEcfp4, RdkitMorganError> {
    // RDKit sanitizes neutral hypervalent halogen oxoacids while parsing (for
    // example `OCl(=O)(=O)=O` becomes `[O-][Cl+3]([O-])([O-])O`). The core
    // SMILES model intentionally preserves the user spelling, so apply this
    // narrow, profile-specific graph normalization only here: Morgan atom and
    // bond invariants must observe the same graph RDKit fingerprints observe.
    let rdkit_input = normalize_rdkit_hypervalent_halogen_oxoacids(mol);
    reject_known_rdkit_coordination_sanitization_gap(&rdkit_input)?;
    let aromatized =
        chematic_perception::apply_aromaticity_rdkit_parity_experimental(&rdkit_input)?;

    let mut result = RdkitMorganEcfp4::default();
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

    let emitted = expand_one_pass(&aromatized, &ring_set, &bond_invariants, ECFP4_RADIUS, true);

    for ((atom_idx, radius), raw_id) in emitted {
        let folded = (raw_id as usize) % ECFP4_FP_SIZE;
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
    }

    Ok(result)
}

/// Return a typed refusal for the one measured RDKit coordination-sanitization
/// gap: Fe(II), degree 10, with at least two directly bonded anionic carbons.
///
/// This deliberately does *not* reject generic high-coordinate metals (or
/// neutral ferrocene). Those inputs remain within the profile unless a
/// reproducible RDKit mismatch establishes a narrower additional boundary.
pub(crate) fn reject_known_rdkit_coordination_sanitization_gap(
    mol: &Molecule,
) -> Result<(), RdkitMorganError> {
    for (atom_idx, atom) in mol.atoms() {
        let atomic_number = atom.element.atomic_number();
        let degree = mol.degree(atom_idx);
        let anionic_carbon_neighbors = mol
            .neighbors(atom_idx)
            .filter(|(neighbor, _)| {
                let neighbor_atom = mol.atom(*neighbor);
                neighbor_atom.element.atomic_number() == 6 && neighbor_atom.charge < 0
            })
            .count();
        if atomic_number == 26 && atom.charge == 2 && degree == 10 && anionic_carbon_neighbors >= 2
        {
            return Err(RdkitMorganError::UnsupportedCoordinationSanitization {
                atom_idx,
                atomic_number,
                degree,
            });
        }
    }
    Ok(())
}

/// Mirror RDKit's parsed representation for neutral chloric/bromic/iodic
/// oxoacids written with one single-bonded oxygen and one or more terminal
/// neutral double-bonded oxygens. Every converted double O becomes `[O-]` and
/// the central halogen obtains the balancing positive charge. The explicit
/// one-single-O requirement keeps ordinary halogen oxides and unrelated
/// hypervalent graphs out of this compatibility-only conversion.
fn normalize_rdkit_hypervalent_halogen_oxoacids(mol: &Molecule) -> Molecule {
    let mut normalizations: Vec<HypervalentHalogenOxoacidNormalization> = Vec::new();
    for (center, atom) in mol.atoms() {
        if !matches!(atom.element, Element::CL | Element::BR | Element::I) || atom.charge != 0 {
            continue;
        }
        let mut single_oxygens = 0usize;
        let mut double_oxygens = Vec::new();
        let mut eligible = true;
        for (neighbor, bond_idx) in mol.neighbors(center) {
            let neighbor_atom = mol.atom(neighbor);
            if neighbor_atom.element != Element::O
                || neighbor_atom.charge != 0
                || mol.neighbors(neighbor).count() != 1
            {
                eligible = false;
                break;
            }
            match mol.bond(bond_idx).order {
                BondOrder::Single => single_oxygens += 1,
                BondOrder::Double => double_oxygens.push((bond_idx, neighbor)),
                _ => {
                    eligible = false;
                    break;
                }
            }
        }
        if eligible && single_oxygens == 1 && !double_oxygens.is_empty() {
            normalizations.push((center, double_oxygens.len() as i8, double_oxygens));
        }
    }

    if normalizations.is_empty() {
        return mol.clone();
    }
    let mut normalized = mol.clone();
    for (center, charge, oxygens) in normalizations {
        normalized.set_charge(center, charge);
        for (bond_idx, oxygen) in oxygens {
            normalized.set_bond_order(bond_idx, BondOrder::Single);
            normalized.set_charge(oxygen, -1);
        }
    }
    normalized
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::MoleculeBuilder;
    use chematic_smiles::parse;

    #[test]
    fn benzene_matches_rdkit_ground_truth_radius0_and_folds_correctly() {
        let mol = parse("c1ccccc1").unwrap();
        let result = rdkit_morgan_ecfp4_experimental(&mol).unwrap();
        // Radius-0 ground truth pinned in `rdkit_morgan_hash`'s own test module.
        assert!(result.raw_bit_info.contains_key(&3218693969));
        for &raw_id in result.raw_bit_info.keys() {
            let folded = (raw_id as usize) % ECFP4_FP_SIZE;
            assert!(result.fingerprint.get(folded));
        }
    }

    #[test]
    fn hermetic_equivalence_to_diagnostic_default_lifecycle() {
        // The production path's raw_bit_info (inverted) must equal
        // `rdkit_morgan_raw_trace`'s `raw_identifier_default`-Some entries on the same
        // already-aromatized molecule -- proves the promotion didn't drift from the path
        // M4-A0 already validated at 5,046/5,046.
        for smi in ["c1ccncc1", "CC(=O)[O-]", "c1ccc2ccccc2c1", "CC"] {
            let mol = parse(smi).unwrap();
            let aromatized =
                chematic_perception::apply_aromaticity_rdkit_parity_experimental(&mol).unwrap();
            let result = rdkit_morgan_ecfp4_experimental(&mol).unwrap();

            let trace = crate::rdkit_morgan_hash::rdkit_morgan_raw_trace(&aromatized, 2);
            let mut expected: Vec<((u32, u32), u32)> = trace
                .iter()
                .filter_map(|e| {
                    e.raw_identifier_default
                        .map(|rid| ((e.atom_idx, e.radius), rid))
                })
                .collect();
            expected.sort_unstable();

            let mut got: Vec<((u32, u32), u32)> = result
                .raw_bit_info
                .iter()
                .flat_map(|(&raw_id, envs)| envs.iter().map(move |&(a, r)| ((a, r), raw_id)))
                .collect();
            got.sort_unstable();

            assert_eq!(got, expected, "mismatch for {smi}");
        }
    }

    /// Was `kekule_pyridinium_reports_kekulization_failed_not_a_fallback_result`, pinned to
    /// pyridinium's `c1cc[nH+]cc1` as its "kekulize fails" example. `fix/kekulize-charge-aware-k1`
    /// (chematic-core's `atom_must_be_matched`) fixed pyridinium's kekulization -- it now
    /// succeeds and is bit-exact against RDKit (see
    /// `validation/results/ecfp4_bitexact_matrix_summary.json`'s `charged_kekulize_fail` bucket,
    /// 6/6 `verified_bit_exact`) -- so it's no longer a valid "kekulize fails" example. This test
    /// The bridgehead-N purine-like ring remains a useful regression fixture: the
    /// matching-based conversion cannot represent it, but its explicit aromatic
    /// parser representation is valid and must be preserved.
    #[test]
    fn bridgehead_n_purine_preserves_explicit_aromatic_input() {
        let smi = "Cc1cn2c(=O)c3ncn(COCCO)c3nc2n1C";
        let mol = parse(smi).unwrap();
        let result = rdkit_morgan_ecfp4_experimental(&mol).expect("valid explicit aromatic input");
        assert!(!result.sparse_counts.is_empty());
    }

    #[test]
    fn large_explicit_polycyclic_aromatic_matches_pinned_rdkit_bits() {
        // RDKit 2025.09.3 and the official RDKit.js 2026.03.6 package agree
        // on this 59-bit set. Re-perceiving an otherwise self-consistent
        // explicit aromatic graph used to change its partition and produce
        // 18 differing folded bits, so retain the source fixture rather than
        // relying on a browser-only comparison for this regression.
        let smi = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../validation/rdkit_issues/fingerprints/large_polycyclic_aromatic.smi"
        ))
        .trim();
        let mol = parse(smi).expect("large explicit aromatic SMILES parses");
        let result = rdkit_morgan_ecfp4_experimental(&mol).expect("parity preprocessing");
        let actual: Vec<usize> = (0..ECFP4_FP_SIZE)
            .filter(|&bit| result.fingerprint.get(bit))
            .collect();
        let expected = vec![
            13, 54, 58, 61, 80, 149, 182, 236, 264, 290, 383, 434, 451, 486, 567, 602, 617, 653,
            691, 695, 719, 822, 841, 883, 893, 926, 935, 949, 1005, 1019, 1039, 1056, 1057, 1060,
            1063, 1096, 1136, 1145, 1151, 1181, 1337, 1364, 1380, 1438, 1453, 1562, 1576, 1582,
            1594, 1648, 1649, 1693, 1706, 1747, 1950, 1984, 1992, 2022, 2027,
        ];
        assert_eq!(actual, expected);
    }

    #[test]
    fn porphyrin_like_fused_system_matches_pinned_rdkit_bits() {
        // RDKit 2025.09.3 and the official RDKit.js 2026.03.6 package agree.
        // This guards the fused-ring neighbour condition used by RDKit
        // aromaticity: rings sharing two bonds are independent candidates,
        // rather than one larger aromatic subsystem.
        let smi =
            "CC1=C2NC(=C1CCC(O)=O)C=C3N=C(C=C4NC(=CC5=NC(=C2)C(=C5C)C=C)C(=C4C)C=C)C(=C3CCC(O)=O)C";
        let mol = parse(smi).expect("porphyrin-like SMILES parses");
        let result = rdkit_morgan_ecfp4_experimental(&mol).expect("parity preprocessing");
        let actual: Vec<usize> = (0..ECFP4_FP_SIZE)
            .filter(|&bit| result.fingerprint.get(bit))
            .collect();
        assert_eq!(
            actual,
            vec![
                8, 46, 80, 119, 204, 227, 252, 294, 350, 378, 389, 562, 578, 591, 644, 650, 694,
                807, 857, 875, 911, 980, 988, 1057, 1093, 1114, 1115, 1132, 1243, 1257, 1287, 1299,
                1336, 1352, 1366, 1375, 1380, 1470, 1564, 1621, 1626, 1645, 1722, 1734, 1737, 1745,
                1801, 1855, 1873, 1898, 1917, 2034,
            ]
        );
    }

    #[test]
    fn measured_feii_coordination_gap_is_a_typed_refusal() {
        let mol = parse("CN(C)C[C-]12C3=C4C5=C1[Fe++]23456789[C-]%10C6=C7C8=C9%10")
            .expect("ferrocene-like SMILES parses");
        assert!(matches!(
            rdkit_morgan_ecfp4_experimental(&mol),
            Err(RdkitMorganError::UnsupportedCoordinationSanitization {
                atomic_number: 26,
                degree: 10,
                ..
            })
        ));
    }

    #[test]
    fn neutral_high_coordinate_iron_remains_supported() {
        let mol = parse("C12C3=C4C5=C1[Fe]23456789C%10C6=C7C8=C9%10")
            .expect("neutral ferrocene SMILES parses");
        assert!(rdkit_morgan_ecfp4_experimental(&mol).is_ok());
    }

    #[test]
    fn normalizes_neutral_perchloric_acid_for_rdkit_morgan_invariants() {
        let input = parse("OCl(=O)(=O)=O").unwrap();
        let normalized = normalize_rdkit_hypervalent_halogen_oxoacids(&input);
        let chlorine = normalized
            .atoms()
            .find_map(|(idx, atom)| (atom.element == Element::CL).then_some(idx))
            .unwrap();
        assert_eq!(normalized.atom(chlorine).charge, 3);
        let mut negative_oxygens = 0;
        for (neighbor, bond) in normalized.neighbors(chlorine) {
            assert_eq!(normalized.bond(bond).order, BondOrder::Single);
            negative_oxygens += (normalized.atom(neighbor).charge == -1) as usize;
        }
        assert_eq!(negative_oxygens, 3);
        let fingerprint = rdkit_morgan_ecfp4_experimental(&input).unwrap();
        let actual: Vec<usize> = (0..ECFP4_FP_SIZE)
            .filter(|&bit| fingerprint.fingerprint.get(bit))
            .collect();
        // Generated by RDKit 2026.03.6 Morgan radius=2 / 2048 bits for the
        // original neutral input spelling. This is intentionally a profile
        // regression: general SMILES parsing remains lossless.
        assert_eq!(actual, vec![187, 222, 669, 715, 807, 2000]);
    }

    #[test]
    fn degree_zero_atom_emits_only_its_radius0_identifier() {
        // Radius 0 is unconditional for every atom (see `rdkit_morgan_hash`'s own
        // `degree_zero_atom_never_appears_past_radius_zero`) -- only radius >= 1 is suppressed
        // by degree-0 death.
        let mol = parse("[Cl-]").unwrap();
        let result = rdkit_morgan_ecfp4_experimental(&mol).unwrap();
        assert_eq!(result.sparse_counts.len(), 1);
        let envs: Vec<_> = result.raw_bit_info.values().flatten().copied().collect();
        assert_eq!(envs, vec![(0, 0)]);
    }

    #[test]
    fn empty_molecule_yields_empty_result_not_an_error() {
        let mol = MoleculeBuilder::new().build();
        let result = rdkit_morgan_ecfp4_experimental(&mol).unwrap();
        assert_eq!(result, RdkitMorganEcfp4::default());
    }

    /// `BondOrder::Query*` has no real RDKit `Bond::BondType` counterpart (see
    /// [`checked_bond_invariant`]'s doc comment) and cannot arise from `parse()` -- built
    /// programmatically here specifically to prove the explicit-`Err` path, not a guessed
    /// mapping, actually fires.
    #[test]
    fn query_bond_order_is_an_explicit_unsupported_bond_order_error_not_a_guess() {
        use chematic_core::{Atom, Element};

        let mut builder = MoleculeBuilder::new();
        let a = builder.add_atom(Atom::new(Element::C));
        let b = builder.add_atom(Atom::new(Element::C));
        builder.add_bond(a, b, BondOrder::QueryAny).unwrap();
        let mol = builder.build();

        match rdkit_morgan_ecfp4_experimental(&mol) {
            Err(RdkitMorganError::UnsupportedBondOrder { bond_idx, order }) => {
                assert_eq!(bond_idx, BondIdx(0));
                assert_eq!(order, BondOrder::QueryAny);
            }
            other => panic!("expected UnsupportedBondOrder, got {other:?}"),
        }
    }

    /// The explicit-aromatic preservation path must not be replaced by the
    /// ordinary Hückel fallback: the two aromatic partitions are intentionally
    /// observable and the RDKit-compatible path must retain the parser's
    /// representation on this valid bridgehead-N input.
    #[test]
    fn hueckel_fallback_would_be_detectable_if_silently_reintroduced() {
        let smi = "Cc1cn2c(=O)c3ncn(COCCO)c3nc2n1C";
        let mol = parse(smi).unwrap();

        let real = rdkit_morgan_ecfp4_experimental(&mol).expect("explicit aromatic input is valid");

        // The ordinary Hückel result remains a separate, observable model.
        let hueckel_fallback_mol = chematic_perception::apply_aromaticity(&mol);
        let hueckel_aromatic_atoms: Vec<bool> = (0..hueckel_fallback_mol.atom_count())
            .map(|i| {
                hueckel_fallback_mol
                    .atom(chematic_core::AtomIdx(i as u32))
                    .aromatic
            })
            .collect();
        assert!(!real.sparse_counts.is_empty());
        assert_eq!(
            hueckel_aromatic_atoms,
            vec![
                false, true, true, true, true, false, true, true, true, true, false, false, false,
                false, false, true, true, true, true, false
            ]
        );
    }
}
