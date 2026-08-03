//! Rotatable-bond and ring classification (Wave 2 spec §5).
//!
//! For every central (B-C) bond candidate, determines: acyclic vs. ring
//! membership (small ring 3-8 / macrocycle 9+ / fused-or-bridged spanning
//! more than one ring), aromatic-ring-bond, amide-like restricted rotation,
//! double/triple bond, and terminal (non-torsional) bond.
//!
//! Ring-size boundary (`MACROCYCLE_MIN = 9`) is not this crate's own
//! invention -- it is RDKit's own `MIN_MACROCYCLE_SIZE`/`minMacrocycleRingSize`
//! constant (`TorsionPreferences.cpp:35`, `BoundsMatrixBuilder.cpp:36`, both
//! fetched and hashed in the sources manifest), kept identical so this PR's
//! small-ring/macrocycle boundary agrees with the RDKit oracle by
//! construction, not by coincidence.

use std::collections::HashMap;

use chematic_core::{AtomIdx, BondOrder, Molecule};
use chematic_perception::{augmented_ring_set, find_sssr};

/// Small ring: 3-8 members (spec §5). Macrocycle: 9+ (RDKit's own
/// `MIN_MACROCYCLE_SIZE`, see module docs).
pub const SMALL_RING_MIN: usize = 3;
pub const SMALL_RING_MAX: usize = 8;
pub const MACROCYCLE_MIN: usize = 9;

/// Which ring-size bucket (if any) a bond falls into.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RingMembership {
    NotInRing,
    SmallRing(usize),
    Macrocycle(usize),
    /// The bond belongs to more than one ring (fused, bridged, or spiro --
    /// this classifier does not further distinguish those three; see
    /// `containing_ring_sizes` on [`BondClassification`] for the raw data).
    /// `chosen_size` is the **smallest** containing ring, per this PR's
    /// documented judgment call (spec §5: "do not naively use only the
    /// first SSSR ring's size" -- this is the explicit, alternative rule
    /// actually used, recorded here rather than left implicit).
    FusedOrBridged {
        chosen_size: usize,
    },
}

impl RingMembership {
    pub fn is_small_ring_only(&self) -> bool {
        matches!(self, RingMembership::SmallRing(_))
    }
    pub fn is_macrocycle_only(&self) -> bool {
        matches!(self, RingMembership::Macrocycle(_))
    }
    pub fn is_fused_or_bridged(&self) -> bool {
        matches!(self, RingMembership::FusedOrBridged { .. })
    }
    pub fn is_in_any_ring(&self) -> bool {
        !matches!(self, RingMembership::NotInRing)
    }
}

/// Precomputed ring membership for every bond in a molecule, built once per
/// molecule and reused across every candidate bond's classification.
///
/// Uses `chematic_perception::augmented_ring_set` over the SSSR base, not
/// raw `find_sssr` output directly -- per this codebase's own documented
/// gap (`CLAUDE.md`'s aromaticity section, `chematic-perception`'s own
/// aromaticity code): plain SSSR can return one large fundamental cycle
/// (e.g. a fused bicyclic's envelope ring) instead of the smaller component
/// rings that actually determine small-ring/macrocycle bucketing. Using the
/// augmented set recovers those smaller rings before classification runs,
/// so (for example) a fused 6+6 aromatic system's individual 6-rings are
/// seen as 6-membered, not folded into one artificial larger ring.
pub struct RingMembershipIndex {
    by_bond: HashMap<(u32, u32), Vec<usize>>,
    /// Every ring INDEX (position into `rings`) the bond belongs to --
    /// unlike `by_bond` (sizes), this preserves ring IDENTITY, so two
    /// different same-sized rings are distinguishable. Pushed in the same
    /// loop iteration as `by_bond` below, so the two can never desync.
    by_bond_ring_ids: HashMap<(u32, u32), Vec<usize>>,
    /// The augmented ring set itself (atom sequences, not just sizes) --
    /// needed by the basic-chemical-knowledge flat-ring rule, which must
    /// walk a specific ring's atom order to build the A-B-C-D tetrad, not
    /// just know its size.
    rings: Vec<Vec<AtomIdx>>,
}

impl RingMembershipIndex {
    pub fn build(mol: &Molecule) -> Self {
        let sssr = find_sssr(mol);
        let rings = augmented_ring_set(mol, sssr.rings());
        let mut by_bond: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
        let mut by_bond_ring_ids: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
        for (ring_id, ring) in rings.iter().enumerate() {
            let n = ring.len();
            if n == 0 {
                continue;
            }
            for i in 0..n {
                let a = ring[i].0;
                let b = ring[(i + 1) % n].0;
                let key = (a.min(b), a.max(b));
                by_bond.entry(key).or_default().push(n);
                by_bond_ring_ids.entry(key).or_default().push(ring_id);
            }
        }
        Self {
            by_bond,
            by_bond_ring_ids,
            rings,
        }
    }

    /// All rings (atom sequences, in ring order) in this molecule's
    /// augmented ring set.
    pub fn rings(&self) -> &[Vec<AtomIdx>] {
        &self.rings
    }

    /// Every ring size (with duplicates if more than one ring of the same
    /// size contains this bond) the bond `(a,b)` belongs to. Empty if the
    /// bond is not in any ring.
    pub fn ring_sizes_for(&self, a: AtomIdx, b: AtomIdx) -> &[usize] {
        self.by_bond
            .get(&(a.0.min(b.0), a.0.max(b.0)))
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// Every ring index (into [`Self::rings`]) the bond `(a,b)` belongs to.
    /// Use this instead of [`Self::ring_sizes_for`] whenever "does this bond
    /// continue the *same* ring as some other bond" matters -- ring size
    /// alone cannot distinguish two different rings of the same size.
    pub fn ring_ids_for(&self, a: AtomIdx, b: AtomIdx) -> &[usize] {
        self.by_bond_ring_ids
            .get(&(a.0.min(b.0), a.0.max(b.0)))
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }
}

/// Everything spec §5 asks to determine "at minimum" for one candidate
/// central bond, plus the raw ring-size evidence and a human-readable
/// reasoning string for diagnostics.
#[derive(Clone, Debug)]
pub struct BondClassification {
    pub bond: (AtomIdx, AtomIdx),
    pub ring: RingMembership,
    /// Every ring size this bond is part of (see
    /// [`RingMembershipIndex::ring_sizes_for`]) -- kept even when `ring` is
    /// `FusedOrBridged` so a diagnostic can show *why* that bucket was
    /// chosen, not just the final answer.
    pub containing_ring_sizes: Vec<usize>,
    pub aromatic_ring_bond: bool,
    pub amide_like: bool,
    pub double_or_triple: bool,
    pub terminal: bool,
    /// `true` only when the bond is in no ring, is not terminal, is not a
    /// double/triple bond, and is not an aromatic-ring bond. A genuinely
    /// acyclic amide bond (e.g. plain N-methylacetamide's C-N bond) is
    /// still `acyclic_rotatable == true` *and* `amide_like == true` at the
    /// same time -- these are independent facts, not a mutually exclusive
    /// tag (spec §5 asks for at least these facts to all be determined,
    /// not for a single winner among them).
    pub acyclic_rotatable: bool,
    pub reasoning: String,
}

/// Classify one candidate central bond `(a, b)`. Panics if `(a, b)` is not
/// an actual bond in `mol` -- callers are expected to only classify bonds
/// they got from [`candidate_central_bonds`] or a real molecule's own bond
/// list, never arbitrary atom pairs.
pub fn classify_bond(
    mol: &Molecule,
    rings: &RingMembershipIndex,
    a: AtomIdx,
    b: AtomIdx,
) -> BondClassification {
    let (_bond_idx, bond) = mol
        .bond_between(a, b)
        .expect("classify_bond requires an existing bond between a and b");

    let ring_sizes = rings.ring_sizes_for(a, b).to_vec();
    let double_or_triple = matches!(bond.order, BondOrder::Double | BondOrder::Triple);
    let deg_a = mol.neighbors(a).count();
    let deg_b = mol.neighbors(b).count();
    let terminal = deg_a <= 1 || deg_b <= 1;
    let aromatic_ring_bond = bond.order == BondOrder::Aromatic
        || (mol.atom(a).aromatic && mol.atom(b).aromatic && !ring_sizes.is_empty());
    let amide_like = is_amide_like(mol, a, b);

    let ring = if ring_sizes.is_empty() {
        RingMembership::NotInRing
    } else {
        let has_small = ring_sizes
            .iter()
            .any(|&s| (SMALL_RING_MIN..=SMALL_RING_MAX).contains(&s));
        let has_macro = ring_sizes.iter().any(|&s| s >= MACROCYCLE_MIN);
        if (has_small && has_macro) || ring_sizes.len() > 1 {
            // Spans both buckets, or simply belongs to more than one ring
            // (fused/bridged/spiro) even within one bucket: flagged as
            // fused/bridged either way (spec §5's explicit ask), smallest
            // containing ring decides `chosen_size`.
            let smallest = *ring_sizes.iter().min().unwrap();
            RingMembership::FusedOrBridged {
                chosen_size: smallest,
            }
        } else {
            let s = ring_sizes[0];
            if s >= MACROCYCLE_MIN {
                RingMembership::Macrocycle(s)
            } else {
                RingMembership::SmallRing(s)
            }
        }
    };

    let acyclic_rotatable = matches!(ring, RingMembership::NotInRing)
        && !terminal
        && !double_or_triple
        && !aromatic_ring_bond;

    let reasoning = format!(
        "bond (atom{},atom{}): containing_ring_sizes={:?} -> {:?}; aromatic_ring_bond={}; amide_like={}; double_or_triple={}; terminal={}; acyclic_rotatable={}",
        a.0,
        b.0,
        ring_sizes,
        ring,
        aromatic_ring_bond,
        amide_like,
        double_or_triple,
        terminal,
        acyclic_rotatable
    );

    BondClassification {
        bond: (a, b),
        ring,
        containing_ring_sizes: ring_sizes,
        aromatic_ring_bond,
        amide_like,
        double_or_triple,
        terminal,
        acyclic_rotatable,
        reasoning,
    }
}

/// Amide/ester-like restricted-rotation bond: an N or O directly bonded to
/// a carbon that itself carries a real C=O double bond, in either traversal
/// direction. Structural (bond-order-based), not SMARTS-based, so it is
/// cheap to check for every candidate bond before running the (more
/// expensive) SMARTS rule matching in `matcher.rs`.
fn is_amide_like(mol: &Molecule, a: AtomIdx, b: AtomIdx) -> bool {
    let check = |n_or_o: AtomIdx, c: AtomIdx| -> bool {
        let elem = mol.atom(n_or_o).element.atomic_number();
        if elem != 7 && elem != 8 {
            return false;
        }
        if mol.atom(c).element.atomic_number() != 6 {
            return false;
        }
        mol.neighbors(c).any(|(o_idx, _)| {
            mol.atom(o_idx).element.atomic_number() == 8
                && mol
                    .bond_between(c, o_idx)
                    .map(|(_, bd)| bd.order == BondOrder::Double)
                    .unwrap_or(false)
        })
    };
    check(a, b) || check(b, a)
}

/// Every bond in `mol` with both endpoints having at least 2 total
/// neighbors (heavy-atom degree, matching this crate's existing convention
/// of counting only heavy atoms -- implicit H are not represented as
/// separate atoms) -- i.e. every bond for which a real, 4-distinct-atom
/// A-B-C-D dihedral can exist at all. Terminal bonds (either endpoint has
/// degree 1) are excluded here rather than merely flagged, since there is no
/// dihedral to define, let alone score, for a bond with no "A" or no "D".
pub fn candidate_central_bonds(mol: &Molecule) -> Vec<(AtomIdx, AtomIdx)> {
    let mut out = Vec::new();
    for (_, bond) in mol.bonds() {
        let a = bond.atom1;
        let b = bond.atom2;
        if mol.neighbors(a).count() >= 2 && mol.neighbors(b).count() >= 2 {
            out.push((a, b));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn acyclic_ethane_bond_is_acyclic_rotatable() {
        // ethane's own C-C bond has degree-1 endpoints in terms of *other*
        // heavy atoms (each C has only H neighbors besides the other C), so
        // it is actually a *terminal* bond by this classifier's definition
        // (no real A/D exists) -- use butane instead for a genuine acyclic
        // rotatable central bond.
        let mol = parse("CCCC").unwrap();
        let rings = RingMembershipIndex::build(&mol);
        let c = classify_bond(&mol, &rings, AtomIdx(1), AtomIdx(2));
        assert!(c.acyclic_rotatable, "{}", c.reasoning);
        assert_eq!(c.ring, RingMembership::NotInRing);
        assert!(!c.terminal);
        assert!(!c.double_or_triple);
        assert!(!c.aromatic_ring_bond);
        assert!(!c.amide_like);
    }

    #[test]
    fn ethane_central_bond_is_terminal() {
        let mol = parse("CC").unwrap();
        let rings = RingMembershipIndex::build(&mol);
        let c = classify_bond(&mol, &rings, AtomIdx(0), AtomIdx(1));
        assert!(c.terminal, "{}", c.reasoning);
        assert!(!c.acyclic_rotatable);
    }

    #[test]
    fn cyclohexane_ring_bond_is_small_ring_6() {
        let mol = parse("C1CCCCC1").unwrap();
        let rings = RingMembershipIndex::build(&mol);
        let c = classify_bond(&mol, &rings, AtomIdx(0), AtomIdx(1));
        assert_eq!(c.ring, RingMembership::SmallRing(6), "{}", c.reasoning);
        assert!(!c.acyclic_rotatable);
    }

    #[test]
    fn cyclononane_ring_bond_is_macrocycle_9() {
        let mol = parse("C1CCCCCCCC1").unwrap();
        let rings = RingMembershipIndex::build(&mol);
        let c = classify_bond(&mol, &rings, AtomIdx(0), AtomIdx(1));
        assert_eq!(c.ring, RingMembership::Macrocycle(9), "{}", c.reasoning);
    }

    #[test]
    fn cyclododecane_ring_bond_is_macrocycle_12() {
        let mol = parse("C1CCCCCCCCCCC1").unwrap();
        let rings = RingMembershipIndex::build(&mol);
        let c = classify_bond(&mol, &rings, AtomIdx(0), AtomIdx(1));
        assert_eq!(c.ring, RingMembership::Macrocycle(12), "{}", c.reasoning);
    }

    #[test]
    fn benzene_ring_bond_is_aromatic() {
        let mol = parse("c1ccccc1").unwrap();
        let rings = RingMembershipIndex::build(&mol);
        let c = classify_bond(&mol, &rings, AtomIdx(0), AtomIdx(1));
        assert!(c.aromatic_ring_bond, "{}", c.reasoning);
        assert_eq!(c.ring, RingMembership::SmallRing(6));
    }

    #[test]
    fn naphthalene_fusion_bond_is_fused_or_bridged() {
        // Naphthalene's central C-C bond (the ring-fusion bond) belongs to
        // both 6-membered rings -- this is the exact case CLAUDE.md's
        // aromaticity section warns SSSR alone can misrepresent.
        let mol = parse("c1ccc2ccccc2c1").unwrap();
        let rings = RingMembershipIndex::build(&mol);
        // Find the fusion bond: the one bond whose both atoms have 3 ring
        // neighbors (degree 3) instead of 2.
        let fusion = mol
            .bonds()
            .map(|(_, b)| (b.atom1, b.atom2))
            .find(|&(a, b)| mol.neighbors(a).count() == 3 && mol.neighbors(b).count() == 3)
            .expect("naphthalene must have a fusion bond");
        let c = classify_bond(&mol, &rings, fusion.0, fusion.1);
        assert!(c.ring.is_fused_or_bridged(), "{}", c.reasoning);
        assert!(c.containing_ring_sizes.len() >= 2, "{}", c.reasoning);
    }

    #[test]
    fn norbornane_bridgehead_bond_is_fused_or_bridged() {
        // Norbornane (bicyclo[2.2.1]heptane): the two bridgehead carbons'
        // shared bond membership -- bridgeheads sit in more than one ring.
        let mol = parse("C1CC2CCC1C2").unwrap();
        let rings = RingMembershipIndex::build(&mol);
        let any_fused = mol.bonds().any(|(_, b)| {
            classify_bond(&mol, &rings, b.atom1, b.atom2)
                .ring
                .is_fused_or_bridged()
        });
        assert!(
            any_fused,
            "norbornane should have at least one bridged bond"
        );
    }

    #[test]
    fn spiro_shared_atom_neighbors_both_classify_independently() {
        // Spiro[4.5]decane-like: two rings sharing exactly one atom (not a
        // bond) -- neither ring's bonds should be flagged fused/bridged
        // (spiro atoms share an ATOM, not a bond, so no bond belongs to two
        // rings here; this is a negative control for over-eager
        // fused/bridged flagging).
        let mol = parse("C1CCC2(CC1)CCCC2").unwrap();
        let rings = RingMembershipIndex::build(&mol);
        for (_, b) in mol.bonds() {
            let c = classify_bond(&mol, &rings, b.atom1, b.atom2);
            if c.ring.is_in_any_ring() {
                assert!(
                    !c.ring.is_fused_or_bridged(),
                    "spiro ring bonds should not be fused/bridged: {}",
                    c.reasoning
                );
            }
        }
    }

    #[test]
    fn amide_bond_is_amide_like_both_directions() {
        let mol = parse("CC(=O)NC").unwrap(); // N-methylacetamide
        let n_idx = (0..mol.atom_count() as u32)
            .map(AtomIdx)
            .find(|&i| mol.atom(i).element.atomic_number() == 7)
            .unwrap();
        let c_carbonyl = (0..mol.atom_count() as u32)
            .map(AtomIdx)
            .find(|&i| {
                mol.atom(i).element.atomic_number() == 6
                    && mol.neighbors(i).any(|(n, _)| {
                        mol.atom(n).element.atomic_number() == 8
                            && mol
                                .bond_between(i, n)
                                .map(|(_, bd)| bd.order == BondOrder::Double)
                                .unwrap_or(false)
                    })
            })
            .unwrap();
        let rings = RingMembershipIndex::build(&mol);
        let c = classify_bond(&mol, &rings, n_idx, c_carbonyl);
        assert!(c.amide_like, "{}", c.reasoning);
        let c_rev = classify_bond(&mol, &rings, c_carbonyl, n_idx);
        assert!(c_rev.amide_like, "{}", c_rev.reasoning);
    }

    #[test]
    fn double_bond_is_flagged_double_or_triple() {
        let mol = parse("C/C=C/C").unwrap();
        let rings = RingMembershipIndex::build(&mol);
        let c = classify_bond(&mol, &rings, AtomIdx(1), AtomIdx(2));
        assert!(c.double_or_triple, "{}", c.reasoning);
        assert!(!c.acyclic_rotatable);
    }

    #[test]
    fn nitrile_triple_bond_is_flagged_double_or_triple() {
        let mol = parse("CCC#N").unwrap();
        let rings = RingMembershipIndex::build(&mol);
        let c = classify_bond(&mol, &rings, AtomIdx(2), AtomIdx(3));
        assert!(c.double_or_triple, "{}", c.reasoning);
    }

    #[test]
    fn candidate_central_bonds_excludes_terminal_bonds() {
        let mol = parse("CC").unwrap(); // ethane: only one bond, both ends terminal
        let candidates = candidate_central_bonds(&mol);
        assert!(candidates.is_empty(), "{candidates:?}");

        let mol2 = parse("CCCC").unwrap(); // butane: middle bond is a real candidate
        let candidates2 = candidate_central_bonds(&mol2);
        assert_eq!(candidates2.len(), 1, "{candidates2:?}");
    }

    #[test]
    fn ring_membership_index_empty_for_acyclic_molecule() {
        let mol = parse("CCCC").unwrap();
        let rings = RingMembershipIndex::build(&mol);
        assert!(rings.ring_sizes_for(AtomIdx(0), AtomIdx(1)).is_empty());
    }
}
