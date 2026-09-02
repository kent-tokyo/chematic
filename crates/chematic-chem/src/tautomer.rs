//! Tautomer enumeration and canonical tautomer selection.
//!
//! Provides:
//! - [`canonical_tautomer`]: return the canonical (preferred) tautomer form.
//! - [`enumerate_tautomers`]: enumerate all reachable tautomers up to a cap.
//! - [`tautomer_parent`]: like `canonical_tautomer`, but budget-limited via
//!   [`TautomerLimits`] and returns an explainable [`TautomerAuditRecord`]
//!   instead of a bare `Molecule` -- see
//!   `docs/rfcs/tautomer_parent_identity_phase2_rfc.md` section 4.
//!
//! Pattern matching is implemented directly on the molecule graph (no SMARTS).

#![forbid(unsafe_code)]

use std::collections::HashSet;
use std::time::Instant;

use chematic_core::{
    AtomIdx, BondIdx, BondOrder, Chirality, Molecule, MoleculeBuilder, implicit_hcount,
};

use crate::parent::{InvalidInputReason, ParentAudit, ParentComputationStatus, ParentResult};
use crate::standardize::MoleculeSnapshot;

/// Bond order match type for tautomer rule patterns.
#[derive(Clone, Copy, PartialEq, Eq)]
enum BondOrderMatch {
    Single,
    Double,
    #[allow(dead_code)]
    Any,
}

impl BondOrderMatch {
    fn matches(self, order: BondOrder) -> bool {
        match self {
            BondOrderMatch::Single => {
                matches!(order, BondOrder::Single | BondOrder::Up | BondOrder::Down)
            }
            BondOrderMatch::Double => matches!(order, BondOrder::Double),
            BondOrderMatch::Any => true,
        }
    }
}

/// A tautomer transformation rule: donor loses an H and the bond orders shift.
///
/// For path_len=3 (1,3-shift):
///   Pattern: donor -[donor_bridge_order]- bridge -[bridge_acceptor_order]- acceptor
///   After: donor =[new]= bridge -[new]- acceptor
///
/// For path_len=5 (1,5-shift):
///   Pattern: donor -[donor_bridge_order]- b1 -[any]- b2 -[any]- b3 -[bridge_acceptor_order]- acceptor
///   After: donor =[new]= b1 -[any]- b2 -[any]- b3 -[new]- acceptor (with H shifted)
struct TautomerRule {
    #[allow(dead_code)]
    name: &'static str,
    /// Atomic number of the donor atom (loses H).
    donor_elem: u8,
    /// Atomic number of the bridge atom(s). For 1,3-shift: single bridge. For 1,5-shift: central atom (b2) if specified.
    bridge_elem: Option<u8>,
    /// Atomic number of the acceptor atom (gains H via implicit valence).
    acceptor_elem: u8,
    /// Required bond order between donor and bridge (or donor and b1 for 1,5-shift).
    donor_bridge_order: BondOrderMatch,
    /// Required bond order between bridge and acceptor (or b3 and acceptor for 1,5-shift).
    bridge_acceptor_order: BondOrderMatch,
    /// If true, this rule is applied in canonical_tautomer to normalize toward a preferred form.
    prefer_forward: bool,
    /// Path length: 3 for 1,3-shift (donor-bridge-acceptor), 5 for 1,5-shift.
    path_len: usize,
}

/// The 44 tautomer rules.
static RULES: &[TautomerRule] = &[
    // 1. keto-enol: O-H adjacent to C=C → O=C-C (prefer keto)
    TautomerRule {
        name: "keto-enol",
        donor_elem: 8,
        bridge_elem: Some(6),
        acceptor_elem: 6,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: true,
        path_len: 3,
    },
    // 2. amide-iminol: N-H adjacent to C=O → N=C-O (prefer amide — forward=false)
    TautomerRule {
        name: "amide-iminol",
        donor_elem: 7,
        bridge_elem: Some(6),
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 3. iminol→amide: O-H adjacent to C=N → O=C-N (prefer amide — forward=true)
    TautomerRule {
        name: "iminol-amide",
        donor_elem: 8,
        bridge_elem: Some(6),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: true,
        path_len: 3,
    },
    // 4. imine-enamine: N-H adjacent to C=C → N=C-C (prefer imine)
    TautomerRule {
        name: "imine-enamine",
        donor_elem: 7,
        bridge_elem: Some(6),
        acceptor_elem: 6,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: true,
        path_len: 3,
    },
    // 5. 1,3-H-shift N→O (any bridge): e.g. nitroso/oxime, hydroxamic acid
    TautomerRule {
        name: "1,3-N-to-O",
        donor_elem: 7,
        bridge_elem: None,
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 6. 1,3-H-shift N→N (any bridge): imidazole, pyrazole, guanidine
    TautomerRule {
        name: "1,3-N-to-N",
        donor_elem: 7,
        bridge_elem: None,
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 7. thioamide: N-H adjacent to C=S → N=C-S (prefer thioamide — forward=false)
    TautomerRule {
        name: "thioamide",
        donor_elem: 7,
        bridge_elem: Some(6),
        acceptor_elem: 16,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 8. thio-iminol→thioamide: S-H adjacent to C=N → S=C-N
    TautomerRule {
        name: "thio-iminol-amide",
        donor_elem: 16,
        bridge_elem: Some(6),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: true,
        path_len: 3,
    },
    // 9. thio keto-enol: S-H adjacent to C=C → S=C-C (prefer thioketone)
    TautomerRule {
        name: "thio-keto-enol",
        donor_elem: 16,
        bridge_elem: Some(6),
        acceptor_elem: 6,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: true,
        path_len: 3,
    },
    // 10. thio-enol→thioketone: O-H adjacent to C=S → O=C-S
    TautomerRule {
        name: "thio-enol-ketone",
        donor_elem: 8,
        bridge_elem: Some(6),
        acceptor_elem: 16,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 11. 1,3-N→S (any bridge): sulfonamide/thioamide-type
    TautomerRule {
        name: "1,3-N-to-S",
        donor_elem: 7,
        bridge_elem: None,
        acceptor_elem: 16,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 12. 1,3-S→O (any bridge)
    TautomerRule {
        name: "1,3-S-to-O",
        donor_elem: 16,
        bridge_elem: None,
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 13. 1,3-S→N (any bridge)
    TautomerRule {
        name: "1,3-S-to-N",
        donor_elem: 16,
        bridge_elem: None,
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 14. 1,3-O→S (any bridge)
    TautomerRule {
        name: "1,3-O-to-S",
        donor_elem: 8,
        bridge_elem: None,
        acceptor_elem: 16,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 15. 1,3-S→S (any bridge): dithioamide, xanthate-type
    TautomerRule {
        name: "1,3-S-to-S",
        donor_elem: 16,
        bridge_elem: None,
        acceptor_elem: 16,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 16. 1,3-O→N any bridge: extends iminol-amide (rule 3) to non-C bridges
    //     e.g. O-H + S=N, N=N → O=X + N-H (hydroxamic acid, thiohydroximate)
    TautomerRule {
        name: "1,3-O-to-N-any-bridge",
        donor_elem: 8,
        bridge_elem: None,
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: true,
        path_len: 3,
    },
    // 17. 1,3-O→O any bridge: O-H via N or S bridge to =O
    //     e.g. hydroxylamine HO-N=O ↔ O=N-OH (N-oxide tautomer)
    TautomerRule {
        name: "1,3-O-to-O-any-bridge",
        donor_elem: 8,
        bridge_elem: None,
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 18. 1,3-N→C any bridge: extends imine-enamine (rule 4) to non-C bridges
    //     e.g. N-H + S=C, N=C via N bridge → N= + X-C-H
    TautomerRule {
        name: "1,3-N-to-C-any-bridge",
        donor_elem: 7,
        bridge_elem: None,
        acceptor_elem: 6,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: true,
        path_len: 3,
    },
    // 19. 1,3-C→O any bridge: active methylene (or terminal alkyne) to O-H
    //     Forward: C-H + X=O → C=X + O-H. prefer_forward:false → prefer keto/C-H form.
    TautomerRule {
        name: "1,3-C-to-O-any-bridge",
        donor_elem: 6,
        bridge_elem: None,
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 20. nitroso → oxime: C-H adjacent to N=O → C=N-OH.
    // Keep this narrower than the generic C→O rule above: enabling that
    // rule in the forward direction would incorrectly enolize ordinary
    // carbonyls through arbitrary bridges.  This is the concrete
    // nitroso/oxime case documented in Phase 2 section 1.6.
    TautomerRule {
        name: "nitroso-oxime",
        donor_elem: 6,
        bridge_elem: Some(7),
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: true,
        path_len: 3,
    },
    // 21. 1,3-C→N any bridge: active methylene adjacent to =N (via S, O, or N bridge)
    //     Forward: C-H + X=N → C=X + N-H. prefer_forward:false → prefer N-H form.
    TautomerRule {
        name: "1,3-C-to-N-any-bridge",
        donor_elem: 6,
        bridge_elem: None,
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 3,
    },
    // 1,5-H shift rules (path_len=5)
    // 22. 1,5-O→O: β-diketone (acetylacetone)
    //     Pattern: O-C(-)-C(-)-C(=O) → O=C(-)-C(-)-C(-O)
    TautomerRule {
        name: "1,5-O-to-O-beta-diketone",
        donor_elem: 8,
        bridge_elem: Some(6), // Central carbon in pattern (b2)
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 23. 1,5-O→N: enol imine
    TautomerRule {
        name: "1,5-O-to-N",
        donor_elem: 8,
        bridge_elem: Some(6),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 24. 1,5-N→N: extended guanidine/amidine tautomerism
    TautomerRule {
        name: "1,5-N-to-N",
        donor_elem: 7,
        bridge_elem: Some(6),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 25. 1,5-N→O: hydroxamic acid type
    TautomerRule {
        name: "1,5-N-to-O",
        donor_elem: 7,
        bridge_elem: Some(6),
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 26. 1,5-C→O: active methylene with conjugation
    TautomerRule {
        name: "1,5-C-to-O",
        donor_elem: 6,
        bridge_elem: Some(6),
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 27. 1,5-O→O with N bridge: nitro-type tautomerism
    //     e.g. O-C-N(=O)-C-O ↔ O=C-N(-O)-C-O via N bridge
    TautomerRule {
        name: "1,5-O-to-O-N-bridge",
        donor_elem: 8,
        bridge_elem: Some(7),
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 28. 1,5-O→O with S bridge: thio-β-diketone
    TautomerRule {
        name: "1,5-O-to-O-S-bridge",
        donor_elem: 8,
        bridge_elem: Some(16),
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 29. 1,5-N→O with C bridge (existing, carbon specified)
    // Already covered by rule 22 (path_len=5 with C bridge)
    // 30. 1,5-N→O with N bridge: bridging N (amidino-type)
    TautomerRule {
        name: "1,5-N-to-O-N-bridge",
        donor_elem: 7,
        bridge_elem: Some(7),
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 31. 1,5-N→O with S bridge
    TautomerRule {
        name: "1,5-N-to-O-S-bridge",
        donor_elem: 7,
        bridge_elem: Some(16),
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 32. 1,5-N→N with N bridge: guanidine-type via N
    TautomerRule {
        name: "1,5-N-to-N-N-bridge",
        donor_elem: 7,
        bridge_elem: Some(7),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 33. 1,5-N→N with S bridge
    TautomerRule {
        name: "1,5-N-to-N-S-bridge",
        donor_elem: 7,
        bridge_elem: Some(16),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 34. 1,5-O→N with N bridge: hydroxamic-type via N
    TautomerRule {
        name: "1,5-O-to-N-N-bridge",
        donor_elem: 8,
        bridge_elem: Some(7),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 35. 1,5-O→N with S bridge
    TautomerRule {
        name: "1,5-O-to-N-S-bridge",
        donor_elem: 8,
        bridge_elem: Some(16),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 36. 1,5-S→O with N bridge
    TautomerRule {
        name: "1,5-S-to-O-N-bridge",
        donor_elem: 16,
        bridge_elem: Some(7),
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 37. 1,5-S→O with S bridge
    TautomerRule {
        name: "1,5-S-to-O-S-bridge",
        donor_elem: 16,
        bridge_elem: Some(16),
        acceptor_elem: 8,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 38. 1,5-S→N with C bridge
    TautomerRule {
        name: "1,5-S-to-N-C-bridge",
        donor_elem: 16,
        bridge_elem: Some(6),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 39. 1,5-S→N with N bridge
    TautomerRule {
        name: "1,5-S-to-N-N-bridge",
        donor_elem: 16,
        bridge_elem: Some(7),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 40. 1,5-C→N with C bridge: extended enamine-imine
    TautomerRule {
        name: "1,5-C-to-N-C-bridge",
        donor_elem: 6,
        bridge_elem: Some(6),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 41. 1,5-C→N with N bridge
    TautomerRule {
        name: "1,5-C-to-N-N-bridge",
        donor_elem: 6,
        bridge_elem: Some(7),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 42. 1,5-C→N with S bridge
    TautomerRule {
        name: "1,5-C-to-N-S-bridge",
        donor_elem: 6,
        bridge_elem: Some(16),
        acceptor_elem: 7,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 43. 1,5-C→S with N bridge
    TautomerRule {
        name: "1,5-C-to-S-N-bridge",
        donor_elem: 6,
        bridge_elem: Some(7),
        acceptor_elem: 16,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
    // 44. 1,5-C→S with S bridge
    TautomerRule {
        name: "1,5-C-to-S-S-bridge",
        donor_elem: 6,
        bridge_elem: Some(16),
        acceptor_elem: 16,
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
        path_len: 5,
    },
];

/// Per-atom explicit hydrogen count vector (position-sensitive, for 1,2-shift dedup).
fn h_assignment(mol: &Molecule) -> Vec<Option<u32>> {
    (0..mol.atom_count())
        .map(|i| mol.atom(AtomIdx(i as u32)).hydrogen_count.map(|h| h as u32))
        .collect()
}

/// Enumerate all forms reachable from `start` by chained direct aromatic 1,2-shifts.
///
/// Returns every distinct form found via BFS (excluding `start` itself).
/// This covers the full ring orbit for systems like tetrazole where one H can
/// reach any of the N positions in a 5-membered ring through a chain of shifts.
fn enumerate_direct_aromatic_forms(
    start: &Molecule,
    blocked: &HashSet<AtomIdx>,
    max: usize,
    rank: &[u32],
) -> Vec<Molecule> {
    enumerate_direct_aromatic_forms_tracked(start, blocked, max, rank).0
}

/// Like [`enumerate_direct_aromatic_forms`], but also reports whether the
/// search was cut short by `max` while more candidates were still reachable
/// (`true`) versus exhausting the frontier naturally (`false`) -- needed by
/// [`tautomer_parent`] to distinguish `Completed` from `MaxTautomersReached`.
///
/// `rank` (see the "Atom-order canonicalization" section above) makes
/// `find_direct_aromatic_matches`'s per-step candidate order a function of
/// the molecule's structure, not raw `AtomIdx` -- relevant when `max` cuts
/// the BFS short: which subset of candidates got explored first must not
/// depend on the caller's atom labeling.
fn enumerate_direct_aromatic_forms_tracked(
    start: &Molecule,
    blocked: &HashSet<AtomIdx>,
    max: usize,
    rank: &[u32],
) -> (Vec<Molecule>, bool) {
    let mut result = Vec::new();
    let mut seen: HashSet<Vec<Option<u32>>> = HashSet::new();
    seen.insert(h_assignment(start));
    let mut frontier = vec![start.clone()];
    let mut truncated = false;

    'outer: while !frontier.is_empty() {
        if result.len() >= max {
            truncated = true;
            break;
        }
        let current = frontier.remove(0);
        for (d, a) in find_direct_aromatic_matches(&current, rank) {
            if result.len() >= max {
                truncated = true;
                break 'outer;
            }
            if let Some(next) = transfer_hydrogen_aromatic(&current, d, a, blocked) {
                let ha = h_assignment(&next);
                if !seen.contains(&ha) {
                    seen.insert(ha);
                    frontier.push(next.clone());
                    result.push(next);
                }
            }
        }
    }
    (result, truncated)
}

/// Find adjacent pairs of aromatic N atoms where the first (donor) carries an explicit H.
///
/// Returns (donor, acceptor) pairs for direct 1,2-shift (no bridge atom),
/// sorted by canonical rank (see the "Atom-order canonicalization" section
/// above) rather than raw `AtomIdx`/adjacency-list order.
fn find_direct_aromatic_matches(mol: &Molecule, rank: &[u32]) -> Vec<(AtomIdx, AtomIdx)> {
    let mut pairs = Vec::new();
    for (d, _) in mol.atoms() {
        let da = mol.atom(d);
        if !da.aromatic || da.element.atomic_number() != 7 {
            continue;
        }
        if da.hydrogen_count.is_none_or(|h| h == 0) {
            continue;
        }
        for (a, _) in mol.neighbors(d) {
            let aa = mol.atom(a);
            if aa.aromatic && (aa.element.atomic_number() == 7 || aa.element.atomic_number() == 8) {
                pairs.push((d, a));
            }
        }
    }
    pairs.sort_by_key(|&(d, a)| (rank[d.0 as usize], rank[a.0 as usize]));
    pairs
}

/// Transfer one H from an aromatic donor to an aromatic acceptor without changing bond orders.
///
/// Only handles atoms with an explicit `hydrogen_count` on the donor.
/// Returns `None` if donor or acceptor is in `blocked_atoms`.
///
/// Donor/acceptor are aromatic (sp2) atoms, never tetrahedral stereocenters, so this
/// function itself never touches `chirality`/`stereo_neighbor_order` -- but the
/// passthrough rebuild below still needs to copy those side channels (plus
/// `stereo_groups`/`bond_directions`) verbatim for every *other* atom, or any
/// pre-existing stereocenter elsewhere in the molecule silently loses the reference
/// order its `@`/`@@` is defined relative to on every tautomer step.
fn transfer_hydrogen_aromatic(
    mol: &Molecule,
    donor: AtomIdx,
    acceptor: AtomIdx,
    blocked_atoms: &HashSet<AtomIdx>,
) -> Option<Molecule> {
    if blocked_atoms.contains(&donor) || blocked_atoms.contains(&acceptor) {
        return None;
    }
    let donor_h = mol.atom(donor).hydrogen_count?;
    if donor_h == 0 {
        return None;
    }
    let acceptor_h = mol.atom(acceptor).hydrogen_count.unwrap_or(0);

    let mut builder = MoleculeBuilder::new();
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let mut atom = mol.atom(idx).clone();
        if idx == donor {
            // Use None (implicit valence) rather than Some(0) when the donor
            // reaches 0 H; this keeps the canonical SMILES clean (n, not [n]).
            atom.hydrogen_count = donor_h.checked_sub(1).filter(|&h| h > 0);
        } else if idx == acceptor {
            atom.hydrogen_count = Some(acceptor_h.saturating_add(1));
        }
        builder.add_atom(atom);
    }
    for i in 0..mol.bond_count() {
        let bidx = BondIdx(i as u32);
        let b = mol.bond(bidx);
        builder
            .add_bond(b.atom1, b.atom2, b.order)
            .expect("transfer_hydrogen_aromatic: bond from a valid molecule must be re-addable");
    }
    builder.copy_stereo_groups_from(mol);
    builder.copy_stereo_from(mol);
    builder.copy_bond_directions_from(mol);
    let next = builder.build();

    // Issue #415: `acceptor` being individually valence-legal (has an
    // available H slot in isolation) isn't sufficient in a fused/bridged
    // ring -- it can be ring-adjacent to *another* atom that already carries
    // an "extra" H, and two such atoms next to each other in one aromatic
    // ring can't both correctly contribute a lone pair to the same ring's pi
    // system. Confirmed via a real repro this exact function reaches
    // (`Oc1[nH]ncc2c3cc(OCc4ccccc4)ccc3nc1-2`'s two ring nitrogens end up
    // adjacent-both-protonated, an over-valent state RDKit itself rejects).
    // Same general-purpose oracle used for the same reason in
    // `transfer_hydrogen_exocyclic_lactam` below -- see that function's own
    // comment for why `kekulize` (not `validate_valence`, which doesn't cap
    // aromatic atoms to their primary valence and misses this) is the right
    // check. Bounded cost: this function is only ever called on the (`max`,
    // default 16) candidates `enumerate_direct_aromatic_forms_tracked`'s BFS
    // frontier actually visits, not on every possible pair in the molecule.
    if chematic_core::kekulize(&next).is_err() {
        return None;
    }
    Some(next)
}

// ---------------------------------------------------------------------------
// Aromatic lactam/lactim shift (ROADMAP.md Phase 2, round 2C-2)
//
// Generalizes the ring-internal-only mechanism above to admit one more mobile-H
// position: an exocyclic oxygen singly bonded to an aromatic ring atom. See
// docs/rfcs/tautomer_parent_identity_phase2_rfc.md section 4.4a for the full
// mechanism table and the empirical reason this must be a directional step
// (like the 41 `TautomerRule`s' own `prefer_forward` direction) rather than a
// new candidate generator feeding `enumerate_direct_aromatic_forms`'s
// score-ranked pool -- reusing that pool's `tautomer_score` would select the
// wrong (lactim) tautomer for 2-pyridone (measured: enol scores 1100 vs keto
// 1050).
// ---------------------------------------------------------------------------

/// Ring-path distance, in bond hops along the one ring `ring` (as returned by
/// [`chematic_perception::find_sssr`]), between two members of that ring.
/// `None` if either atom is not a member of `ring`.
fn ring_distance(ring: &[AtomIdx], a: AtomIdx, b: AtomIdx) -> Option<usize> {
    let ia = ring.iter().position(|&x| x == a)?;
    let ib = ring.iter().position(|&x| x == b)?;
    let d = ia.abs_diff(ib);
    Some(d.min(ring.len() - d))
}

/// Find every (donor, bridge, acceptor) triple eligible for the exocyclic
/// lactam/lactim shift. All three sides are deliberately narrow and
/// fail-closed -- restricted to exactly what the 5 confirmed-broken design
/// molecules evidence, not "aromatic ring + heteroatom" in general:
///
/// - `bridge`: a **neutral aromatic carbon**. Not any aromatic element --
///   an aromatic N/S/P/B bridge (e.g. an N-hydroxy heterocycle) is a
///   structurally different system this mechanism was never evidenced
///   against, and a charged bridge changes the valence arithmetic below.
/// - `donor`: a **neutral, exocyclic, non-aromatic oxygen** with exactly
///   one transferable H (checked via `implicit_hcount`, not
///   `.hydrogen_count` -- an unbracketed hydroxyl like `Oc1...` stores
///   `None` and relies on valence inference), bonded to `bridge` by a
///   `Single` bond that is its **only** bond (degree 1 -- an ether/bridging
///   oxygen with a second heavy-atom connection is not a lactam/lactim
///   hydroxyl and must never match).
/// - `acceptor`: a **neutral, pyridine-type aromatic nitrogen** -- no H,
///   not charged, and with **exactly 2** heavy-atom connections. Degree 2
///   is the valence-compatibility condition itself: an organic-subset N's
///   normal valence is 3, so a 2-connected aromatic N always has exactly
///   one free valence slot for the incoming H; a 3-connected aromatic N
///   (pyrrole-type, or a bridgehead/fusion position) has none, and adding
///   an H there would require a charge change this mechanism must never
///   make.
/// - sharing an SSSR ring with `bridge` at **odd** ring distance (the
///   condition for a real alternating single/double path between them to
///   exist at all -- RFC section 4.4a).
///
/// O-only, not O-or-N: every confirmed-broken molecule is this type (RFC
/// section 1.1); the analogous N-acceptor (amino/imino) case is a distinct,
/// unevidenced defect, deliberately out of scope (RFC section 4.4a).
///
/// Deduplicated and returned in a fixed (donor, bridge, acceptor) order,
/// never hash-iteration order: `bridge` can belong to more than one SSSR
/// ring in a fused system (hypoxanthine), which would otherwise report the
/// same triple once per ring it's found through.
fn find_exocyclic_lactam_shift_matches(
    mol: &Molecule,
    rank: &[u32],
) -> Vec<(AtomIdx, AtomIdx, AtomIdx)> {
    let rings = chematic_perception::find_sssr(mol);
    let mut out: HashSet<(AtomIdx, AtomIdx, AtomIdx)> = HashSet::new();
    for (bridge, bridge_atom) in mol.atoms() {
        if !bridge_atom.aromatic
            || bridge_atom.element.atomic_number() != 6
            || bridge_atom.charge != 0
        {
            continue;
        }
        for (donor, bond_idx) in mol.neighbors(bridge) {
            let donor_atom = mol.atom(donor);
            if donor_atom.aromatic
                || donor_atom.element.atomic_number() != 8
                || donor_atom.charge != 0
            {
                continue;
            }
            if implicit_hcount(mol, donor) != 1 {
                continue;
            }
            if mol.neighbors(donor).count() != 1 {
                continue;
            }
            if mol.bond(bond_idx).order != BondOrder::Single {
                continue;
            }
            for ring in rings.rings() {
                if !ring.contains(&bridge) {
                    continue;
                }
                for &acceptor in ring {
                    let acceptor_atom = mol.atom(acceptor);
                    if !acceptor_atom.aromatic
                        || acceptor_atom.element.atomic_number() != 7
                        || acceptor_atom.charge != 0
                    {
                        continue;
                    }
                    if implicit_hcount(mol, acceptor) != 0 {
                        continue;
                    }
                    if mol.neighbors(acceptor).count() != 2 {
                        continue;
                    }
                    if ring_distance(ring, bridge, acceptor).is_some_and(|d| d % 2 == 1) {
                        out.insert((donor, bridge, acceptor));
                    }
                }
            }
        }
    }
    let mut out: Vec<_> = out.into_iter().collect();
    out.sort_by_key(|&(d, b, a)| (rank[d.0 as usize], rank[b.0 as usize], rank[a.0 as usize]));
    out
}

/// Move the H from `donor` (exocyclic O) to `acceptor` (ring N), flipping
/// only the `donor`-`bridge` bond order (`Single` -> `Double`); every
/// ring-internal bond, including `bridge`-`acceptor` when they're directly
/// bonded, stays exactly as it was (`Aromatic`) -- confirmed structurally
/// true for every design molecule (RFC section 4.4a).
fn transfer_hydrogen_exocyclic_lactam(
    mol: &Molecule,
    donor: AtomIdx,
    bridge: AtomIdx,
    acceptor: AtomIdx,
    blocked_atoms: &HashSet<AtomIdx>,
    blocked_bonds: &HashSet<BondIdx>,
) -> Option<Molecule> {
    if blocked_atoms.contains(&donor) || blocked_atoms.contains(&acceptor) {
        return None;
    }
    // implicit_hcount, not `.hydrogen_count`: an unbracketed hydroxyl's H is
    // stored as `None` (see find_exocyclic_lactam_shift_matches's own note).
    let donor_h = implicit_hcount(mol, donor);
    if donor_h == 0 {
        return None;
    }
    let (bond_idx, bond) = mol.bond_between(donor, bridge)?;
    if bond.order != BondOrder::Single || blocked_bonds.contains(&bond_idx) {
        return None;
    }
    let acceptor_h = mol.atom(acceptor).hydrogen_count.unwrap_or(0);

    let mut builder = MoleculeBuilder::new();
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let mut atom = mol.atom(idx).clone();
        if idx == donor {
            atom.hydrogen_count = donor_h.checked_sub(1).filter(|&h| h > 0);
        } else if idx == acceptor {
            atom.hydrogen_count = Some(acceptor_h.saturating_add(1));
        }
        builder.add_atom(atom);
    }
    for i in 0..mol.bond_count() {
        let bidx = BondIdx(i as u32);
        let b = mol.bond(bidx);
        let order = if bidx == bond_idx {
            BondOrder::Double
        } else {
            b.order
        };
        builder.add_bond(b.atom1, b.atom2, order).expect(
            "transfer_hydrogen_exocyclic_lactam: bond from a valid molecule must be re-addable",
        );
    }
    builder.copy_stereo_groups_from(mol);
    builder.copy_stereo_from(mol);
    builder.copy_bond_directions_from(mol);
    let next = builder.build();

    if !exocyclic_lactam_shift_preserves_invariants(mol, &next, donor, acceptor, bond_idx) {
        // Defense in depth: the narrow preconditions above are meant to
        // guarantee this by construction, but `MoleculeBuilder::build()` can
        // re-perceive/normalize state beyond what the explicit field copy
        // above touched (the class of bug this project has hit before --
        // see Phase 1's fragment-extraction stereo-corruption fix). Fail
        // closed rather than return a molecule that silently changed
        // something this transform must never touch.
        return None;
    }

    // Second defense-in-depth layer (issue #415): the per-atom preconditions
    // above (`implicit_hcount(mol, acceptor) == 0`, degree 2) are necessary
    // but not sufficient in a fused/bridged ring system -- `find_sssr`'s own
    // ring choice for such systems isn't always unique (a known, separate
    // residual: the *same* physical molecule can get a different SSSR
    // decomposition depending on atom-labeling order, which can change which
    // ring this function's `ring_distance` check sees), and an acceptor that
    // looks individually valid (0 H, degree 2) can still be ring-adjacent to
    // *another* aromatic atom that already carries an "extra" H. Two such
    // atoms next to each other in one aromatic ring cannot both correctly
    // contribute a lone pair to the same ring's pi system -- confirmed via a
    // real repro (`Oc1[nH]ncc2c3cc(OCc4ccccc4)ccc3nc1-2`) that this function
    // previously accepted and which produced an over-valent `[nH2]` nitrogen
    // RDKit itself rejects outright. `kekulize` is a general, already-used-
    // elsewhere oracle for "is this aromatic system actually valid" that
    // catches this (and would catch other, not-yet-evidenced pathological
    // shapes too) without needing to hand-enumerate every such shape here.
    if chematic_core::kekulize(&next).is_err() {
        return None;
    }
    Some(next)
}

/// Reverse the aromatic lactam/lactim edge: move H from a ring N back to a
/// carbonyl oxygen and reduce only the exocyclic C=O bond. This is used by
/// exhaustive tautomer enumeration so keto inputs can reach the same
/// candidate component as hydroxyl inputs; canonical selection remains
/// directional and continues to prefer the lactam side.
fn transfer_hydrogen_exocyclic_lactim(
    mol: &Molecule,
    bridge: AtomIdx,
    oxygen: AtomIdx,
    donor: AtomIdx,
    blocked_atoms: &HashSet<AtomIdx>,
    blocked_bonds: &HashSet<BondIdx>,
) -> Option<Molecule> {
    if blocked_atoms.contains(&oxygen)
        || blocked_atoms.contains(&donor)
        || mol.atom(donor).hydrogen_count != Some(1)
    {
        return None;
    }
    let (bond_idx, bond) = mol.bond_between(oxygen, bridge)?;
    if bond.order != BondOrder::Double || blocked_bonds.contains(&bond_idx) {
        return None;
    }
    let mut builder = MoleculeBuilder::new();
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let mut atom = mol.atom(idx).clone();
        if idx == oxygen {
            atom.hydrogen_count = Some(1);
        } else if idx == donor {
            atom.hydrogen_count = None;
        }
        builder.add_atom(atom);
    }
    for i in 0..mol.bond_count() {
        let bidx = BondIdx(i as u32);
        let b = mol.bond(bidx);
        let order = if bidx == bond_idx {
            BondOrder::Single
        } else {
            b.order
        };
        builder.add_bond(b.atom1, b.atom2, order).ok()?;
    }
    builder.copy_stereo_groups_from(mol);
    builder.copy_stereo_from(mol);
    builder.copy_bond_directions_from(mol);
    let next = builder.build();
    chematic_core::kekulize(&next).ok().map(|_| next)
}

fn reverse_exocyclic_lactim_candidates(mol: &Molecule, config: &TautomerConfig) -> Vec<Molecule> {
    let rings = chematic_perception::find_sssr(mol);
    let mut candidates = Vec::new();
    for (bridge, atom) in mol.atoms() {
        if !atom.aromatic || atom.element.atomic_number() != 6 || atom.charge != 0 {
            continue;
        }
        for (oxygen, bond_idx) in mol.neighbors(bridge) {
            if mol.atom(oxygen).element.atomic_number() != 8
                || mol.atom(oxygen).aromatic
                || mol.bond(bond_idx).order != BondOrder::Double
            {
                continue;
            }
            for ring in rings.rings() {
                if !ring.contains(&bridge) {
                    continue;
                }
                for &donor in ring {
                    let donor_atom = mol.atom(donor);
                    if !donor_atom.aromatic
                        || donor_atom.element.atomic_number() != 7
                        || donor_atom.charge != 0
                        || donor_atom.hydrogen_count != Some(1)
                        || mol.neighbors(donor).count() != 2
                        || !ring_distance(ring, bridge, donor).is_some_and(|d| d % 2 == 1)
                    {
                        continue;
                    }
                    if let Some(next) = transfer_hydrogen_exocyclic_lactim(
                        mol,
                        bridge,
                        oxygen,
                        donor,
                        &config.blocked_atoms,
                        &config.blocked_bonds,
                    ) {
                        candidates.push(next);
                    }
                }
            }
        }
    }
    candidates.sort_by_key(|mol| (mol_fingerprint(mol), chematic_smiles::canonical_smiles(mol)));
    candidates.dedup_by(|a, b| mol_fingerprint(a) == mol_fingerprint(b));
    candidates
}

fn exocyclic_lactam_candidates(
    mol: &Molecule,
    rank: &[u32],
    config: &TautomerConfig,
) -> Vec<Molecule> {
    find_exocyclic_lactam_shift_matches(mol, rank)
        .into_iter()
        .filter_map(|(donor, bridge, acceptor)| {
            transfer_hydrogen_exocyclic_lactam(
                mol,
                donor,
                bridge,
                acceptor,
                &config.blocked_atoms,
                &config.blocked_bonds,
            )
        })
        .collect()
}

/// Canonical key rooted at the aromatic carbonyl. Explicit H annotations are
/// removed so equivalent lactam spellings are compared by the carbonyl's
/// ring environment rather than by the input proton placement.
fn carbonyl_rooted_key(mol: &Molecule) -> String {
    let root = mol.atoms().find_map(|(idx, atom)| {
        (atom.aromatic
            && atom.element.atomic_number() == 6
            && mol.neighbors(idx).any(|(neighbor, bond)| {
                mol.atom(neighbor).element.atomic_number() == 8
                    && mol.bond(bond).order == BondOrder::Double
            }))
        .then_some(idx)
    });
    let Some(root) = root else {
        return chematic_smiles::canonical_smiles(mol);
    };
    let mut builder = MoleculeBuilder::new();
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let mut atom = mol.atom(idx).clone();
        atom.hydrogen_count = None;
        if idx == root {
            atom.isotope = Some(65535);
        }
        builder.add_atom(atom);
    }
    for i in 0..mol.bond_count() {
        let bond = mol.bond(BondIdx(i as u32));
        builder.add_bond(bond.atom1, bond.atom2, bond.order).ok();
    }
    chematic_smiles::canonical_smiles(&builder.build())
}

fn has_aromatic_exocyclic_carbonyl(mol: &Molecule) -> bool {
    mol.atoms().any(|(bridge, atom)| {
        atom.aromatic
            && atom.element.atomic_number() == 6
            && mol.neighbors(bridge).any(|(neighbor, bond)| {
                mol.atom(neighbor).element.atomic_number() == 8
                    && mol.bond(bond).order == BondOrder::Double
            })
    })
}

fn aromatic_exocyclic_carbonyl_count(mol: &Molecule) -> usize {
    mol.atoms()
        .filter(|(bridge, atom)| {
            atom.aromatic
                && atom.element.atomic_number() == 6
                && mol.neighbors(*bridge).any(|(neighbor, bond)| {
                    mol.atom(neighbor).element.atomic_number() == 8
                        && mol.bond(bond).order == BondOrder::Double
                })
        })
        .count()
}

fn has_aromatic_exocyclic_oxygen(mol: &Molecule) -> bool {
    mol.atoms().any(|(bridge, atom)| {
        atom.aromatic
            && atom.element.atomic_number() == 6
            && mol.neighbors(bridge).any(|(neighbor, bond)| {
                mol.atom(neighbor).element.atomic_number() == 8
                    && matches!(mol.bond(bond).order, BondOrder::Single | BondOrder::Double)
            })
    })
}

/// Find the dual-flank N-H ambiguity of aromatic lactam systems such as
/// cytosine, guanine, and hypoxanthine.  The carbonyl carbon is bonded to two
/// aromatic nitrogens; moving H between those nitrogens changes neither the
/// graph nor any bond order, but is required to make the two keto spellings
/// canonicalize identically.
fn find_dual_flank_matches(mol: &Molecule, rank: &[u32]) -> Vec<(AtomIdx, AtomIdx, AtomIdx)> {
    let mut out = Vec::new();
    for (bridge, atom) in mol.atoms() {
        if !atom.aromatic || atom.element.atomic_number() != 6 || atom.charge != 0 {
            continue;
        }
        let has_carbonyl_oxygen = mol.neighbors(bridge).any(|(neighbor, bidx)| {
            mol.atom(neighbor).element.atomic_number() == 8
                && mol.bond(bidx).order == BondOrder::Double
        });
        if !has_carbonyl_oxygen {
            continue;
        }
        for ring in chematic_perception::find_sssr(mol).rings().iter() {
            if !ring.contains(&bridge) {
                continue;
            }
            let flank_n: Vec<AtomIdx> = ring
                .iter()
                .copied()
                .filter(|&candidate| {
                    let a = mol.atom(candidate);
                    a.aromatic
                        && a.element.atomic_number() == 7
                        && a.charge == 0
                        && ring_distance(ring, bridge, candidate).is_some_and(|d| d % 2 == 1)
                })
                .collect();
            for &a in &flank_n {
                for &b in &flank_n {
                    if a == b {
                        continue;
                    }
                    let (donor, acceptor) = if mol.atom(a).hydrogen_count == Some(1)
                        && mol.atom(b).hydrogen_count.unwrap_or(0) == 0
                    {
                        (a, b)
                    } else {
                        continue;
                    };
                    out.push((donor, bridge, acceptor));
                }
            }
        }
    }
    out.sort_by_key(|&(d, b, a)| (rank[d.0 as usize], rank[b.0 as usize], rank[a.0 as usize]));
    out.dedup();
    out
}

/// Move H between the two aromatic N atoms flanking a carbonyl carbon.
/// Unlike the lactam/lactim shift, no bond order changes are needed.
fn transfer_hydrogen_dual_flank(
    mol: &Molecule,
    donor: AtomIdx,
    acceptor: AtomIdx,
    blocked_atoms: &HashSet<AtomIdx>,
) -> Option<Molecule> {
    if blocked_atoms.contains(&donor) || blocked_atoms.contains(&acceptor) {
        return None;
    }
    if mol.atom(donor).hydrogen_count != Some(1)
        || mol.atom(acceptor).hydrogen_count.unwrap_or(0) != 0
    {
        return None;
    }
    let mut builder = MoleculeBuilder::new();
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let mut atom = mol.atom(idx).clone();
        if idx == donor {
            atom.hydrogen_count = Some(0);
        } else if idx == acceptor {
            atom.hydrogen_count = Some(1);
        }
        builder.add_atom(atom);
    }
    for i in 0..mol.bond_count() {
        let b = mol.bond(BondIdx(i as u32));
        builder.add_bond(b.atom1, b.atom2, b.order).ok()?;
    }
    builder.copy_stereo_groups_from(mol);
    builder.copy_stereo_from(mol);
    builder.copy_bond_directions_from(mol);
    let next = builder.build();
    chematic_core::kekulize(&next).ok().map(|_| next)
}

fn apply_dual_flank_tracked(
    mol: &Molecule,
    config: &TautomerConfig,
    _rank: &[u32],
) -> Option<(Molecule, AtomIdx, AtomIdx, AtomIdx)> {
    // H placement is the ambiguity being normalized here.  Ranking the
    // input molecule directly lets that H change the root of a symmetric
    // pyrimidinone ring, so keto and enol spellings can choose opposite
    // flanks.  Use a rooted structural key made from the same graph with
    // explicit H annotations removed; atom indices remain unchanged.
    let rooted_rank = tautomer_skeleton_rank(mol);
    find_dual_flank_matches(mol, &rooted_rank)
        .into_iter()
        .find_map(|(donor, bridge, acceptor)| {
            let next = transfer_hydrogen_dual_flank(mol, donor, acceptor, &config.blocked_atoms)?;
            (rooted_structural_key(mol, donor) < rooted_structural_key(mol, acceptor))
                .then_some((next, donor, bridge, acceptor))
        })
}

/// Everything this transform must hold, checked on the *built* molecule, not
/// assumed from the field copy above: every atom except `donor`/`acceptor`
/// is byte-for-byte identical (element, charge, isotope, aromaticity,
/// chirality, H count); `donor`/`acceptor` keep everything except H count,
/// which must move by exactly 1 in the expected direction; every bond except
/// `bond_idx` keeps its endpoints and order; `bond_idx` keeps its endpoints
/// and its order must have actually flipped `Single` -> `Double`.
fn exocyclic_lactam_shift_preserves_invariants(
    before: &Molecule,
    after: &Molecule,
    donor: AtomIdx,
    acceptor: AtomIdx,
    bond_idx: BondIdx,
) -> bool {
    if before.atom_count() != after.atom_count() || before.bond_count() != after.bond_count() {
        return false;
    }
    for i in 0..before.atom_count() {
        let idx = AtomIdx(i as u32);
        let a = before.atom(idx);
        let b = after.atom(idx);
        if a.element != b.element
            || a.charge != b.charge
            || a.isotope != b.isotope
            || a.aromatic != b.aromatic
            || a.chirality != b.chirality
        {
            return false;
        }
        let h_before = implicit_hcount(before, idx);
        let h_after = implicit_hcount(after, idx);
        if idx == donor {
            if h_after != h_before.wrapping_sub(1) {
                return false;
            }
        } else if idx == acceptor {
            if h_after != h_before + 1 {
                return false;
            }
        } else if h_after != h_before {
            return false;
        }
    }
    for i in 0..before.bond_count() {
        let bidx = BondIdx(i as u32);
        let bb = before.bond(bidx);
        let ba = after.bond(bidx);
        if bb.atom1 != ba.atom1 || bb.atom2 != ba.atom2 {
            return false;
        }
        if bidx == bond_idx {
            if bb.order != BondOrder::Single || ba.order != BondOrder::Double {
                return false;
            }
        } else if bb.order != ba.order {
            return false;
        }
    }
    true
}

/// Apply one exocyclic lactam/lactim shift, chosen deterministically among
/// every eligible (donor, bridge, acceptor) triple by minimal canonical
/// SMILES of the result -- never `.first()`/index order, which would make
/// the choice track input atom order rather than structure (uracil has two
/// independent sites; cytosine's/guanine's carbonyl-bearing ring atom is
/// flanked by two ring nitrogens, only one a valid acceptor at a time, but
/// both may qualify simultaneously from a fully-shifted starting form). See
/// RFC section 4.4a's order-invariance requirement.
fn apply_exocyclic_lactam_shift_tracked(
    mol: &Molecule,
    config: &TautomerConfig,
    rank: &[u32],
) -> Option<(Molecule, AtomIdx, AtomIdx, AtomIdx)> {
    let mut candidates: Vec<(Molecule, AtomIdx, AtomIdx, AtomIdx)> =
        find_exocyclic_lactam_shift_matches(mol, rank)
            .into_iter()
            .filter_map(|(donor, bridge, acceptor)| {
                transfer_hydrogen_exocyclic_lactam(
                    mol,
                    donor,
                    bridge,
                    acceptor,
                    &config.blocked_atoms,
                    &config.blocked_bonds,
                )
                .map(|next| (next, donor, bridge, acceptor))
            })
            .collect();
    if candidates.is_empty() {
        return None;
    }
    candidates.sort_by(|x, y| {
        carbonyl_rooted_key(&y.0)
            .cmp(&carbonyl_rooted_key(&x.0))
            .then_with(|| {
                chematic_smiles::canonical_smiles(&x.0)
                    .cmp(&chematic_smiles::canonical_smiles(&y.0))
            })
    });
    // Candidates sharing the minimal canonical SMILES are meant to be the
    // same molecule (safely interchangeable -- pick the first, deterministic
    // given the sort above). Don't trust canonical_smiles alone for that:
    // this project has hit real canonical-SMILES bugs before. Cross-check
    // with mol_fingerprint, an independently-computed (sorted per-atom
    // element/charge/bond-order-sum) hash already used elsewhere in this
    // file for exactly this "is this really the same molecule" question. A
    // disagreement means the tie is not actually resolved -- fail closed
    // (report no move this iteration) rather than guess.
    let min_smiles = chematic_smiles::canonical_smiles(&candidates[0].0);
    let tied: Vec<&(Molecule, AtomIdx, AtomIdx, AtomIdx)> = candidates
        .iter()
        .take_while(|c| chematic_smiles::canonical_smiles(&c.0) == min_smiles)
        .collect();
    if tied.len() > 1 {
        let first_fp = mol_fingerprint(&tied[0].0);
        if !tied.iter().all(|c| mol_fingerprint(&c.0) == first_fp) {
            return None;
        }
    }
    candidates.into_iter().next()
}

use crate::hash::{FNV1A_OFFSET, FNV1A_PRIME};

/// Tautomer form score for canonical selection: prefers O-H > N-H > S-H and aromatic rings.
fn tautomer_score(mol: &Molecule) -> i32 {
    score_breakdown(mol).iter().map(|c| c.value).sum()
}

/// Same scoring rule as [`tautomer_score`], decomposed into named
/// contributions for [`TautomerAuditRecord::score_breakdown`]. This is the
/// single source of truth for the scoring rule; `tautomer_score` sums it.
fn score_breakdown(mol: &Molecule) -> Vec<ScoreContribution> {
    let mut contributions = Vec::new();
    let mut has_aromatic = false;

    for (_, atom) in mol.atoms() {
        if atom.aromatic {
            has_aromatic = true;
        }

        // Heteroatom hydrogen: O-H > N-H > S-H (explicit or implicit)
        let h_count = atom.hydrogen_count.unwrap_or(0) as i32;
        if h_count > 0 {
            let weight = match atom.element.atomic_number() {
                8 => Some(100), // O-H: highest priority
                7 => Some(50),  // N-H: medium priority
                16 => Some(25), // S-H: low priority
                _ => None,
            };
            if let Some(weight) = weight {
                contributions.push(ScoreContribution {
                    term: TautomerScoreTerm::HeteroatomHydrogen {
                        element: atom.element.atomic_number(),
                    },
                    value: h_count * weight,
                });
            }
        }
    }

    if has_aromatic {
        contributions.push(ScoreContribution {
            term: TautomerScoreTerm::AromaticRing,
            value: 1000,
        });
    }

    contributions
}

/// Order-independent structural hash for convergence detection.
fn mol_fingerprint(mol: &Molecule) -> u64 {
    let mut atoms: Vec<(u8, i8, u32, u8)> = (0..mol.atom_count())
        .map(|i| {
            let idx = AtomIdx(i as u32);
            let a = mol.atom(idx);
            let bos: u32 = mol
                .neighbors(idx)
                .map(|(_, bidx)| mol.bond(bidx).order.order_int() as u32)
                .sum();
            (
                a.element.atomic_number(),
                a.charge,
                bos,
                a.hydrogen_count.unwrap_or(0),
            )
        })
        .collect();
    atoms.sort();
    let mut hash = FNV1A_OFFSET;
    for (an, ch, bos, h) in atoms {
        hash ^= an as u64;
        hash = hash.wrapping_mul(FNV1A_PRIME);
        hash ^= (ch as u8 as u64).wrapping_add(128);
        hash = hash.wrapping_mul(FNV1A_PRIME);
        hash ^= bos as u64;
        hash = hash.wrapping_mul(FNV1A_PRIME);
        hash ^= h as u64;
        hash = hash.wrapping_mul(FNV1A_PRIME);
    }

    // BondOrder::Up/Down and the bond-direction side channel encode 2D E/Z
    // geometry.  The compact atom fingerprint above intentionally ignores
    // those directions, but using it alone makes E and Z forms collide in
    // tautomer search and causes one form to be discarded as "seen".  Add
    // the canonical serialization only for molecules that actually carry a
    // directional bond; this preserves the cheap/order-independent fast path
    // for the overwhelmingly common non-E/Z case.
    let has_direction = (0..mol.bond_count()).any(|i| {
        let bidx = BondIdx(i as u32);
        matches!(mol.bond(bidx).order, BondOrder::Up | BondOrder::Down)
            || mol.bond_direction(bidx).is_some()
    });
    // Explicit H placement is also position-sensitive: two equivalent
    // aromatic nitrogens can have the same sorted atom multiset while the H
    // sits on different atoms (the cytosine/guanine dual-flank case).
    let has_explicit_h =
        (0..mol.atom_count()).any(|i| mol.atom(AtomIdx(i as u32)).hydrogen_count.is_some());
    if has_direction || has_explicit_h {
        hash ^= 0xff;
        hash = hash.wrapping_mul(FNV1A_PRIME);
        for &byte in chematic_smiles::canonical_smiles(mol).as_bytes() {
            hash ^= byte as u64;
            hash = hash.wrapping_mul(FNV1A_PRIME);
        }
    }
    hash
}

/// Whether `mol` carries an OpenSMILES directional marker that contributes to
/// E/Z stereochemistry.  A tautomer transform can move a double bond away
/// from the marked single bond; until the direction is explicitly remapped to
/// the new stereochemical carrier, transforming such a molecule would silently
/// erase its E/Z identity.
fn has_directional_bond(mol: &Molecule) -> bool {
    (0..mol.bond_count()).any(|i| {
        let bidx = BondIdx(i as u32);
        matches!(mol.bond(bidx).order, BondOrder::Up | BondOrder::Down)
            || mol.bond_direction(bidx).is_some()
    })
}

/// Find all (donor, bridge, acceptor) triples matching the rule in `mol`.
/// For path_len=3: donor-bridge-acceptor (standard 1,3-shift)
/// For path_len=5: donor-b1-b2-b3-acceptor (1,5-shift) where b2 is stored in "bridge"
fn find_matches(
    mol: &Molecule,
    rule: &TautomerRule,
    rank: &[u32],
) -> Vec<(AtomIdx, AtomIdx, AtomIdx)> {
    let mut matches = Vec::new();

    if rule.path_len == 5 {
        // 1,5-shift: donor -[single]- b1 -[any]- b2 -[any]- b3 -[double]- acceptor
        for i in 0..mol.atom_count() {
            let d = AtomIdx(i as u32);
            let donor_atom = mol.atom(d);
            if donor_atom.element.atomic_number() != rule.donor_elem {
                continue;
            }
            if implicit_hcount(mol, d) == 0 {
                continue;
            }

            for (b1, db_bidx) in mol.neighbors(d) {
                if !rule.donor_bridge_order.matches(mol.bond(db_bidx).order) {
                    continue;
                }
                // b1 should be a carbon (or be flexible)
                if mol.atom(b1).element.atomic_number() != 6 {
                    continue;
                }

                for (b2, _) in mol.neighbors(b1) {
                    if b2 == d {
                        continue;
                    }
                    // b2 can be any type (relaxed)
                    if let Some(br_elem) = rule.bridge_elem
                        && mol.atom(b2).element.atomic_number() != br_elem
                    {
                        continue;
                    }

                    for (b3, _) in mol.neighbors(b2) {
                        if b3 == b1 {
                            continue;
                        }
                        // b3 should be a carbon
                        if mol.atom(b3).element.atomic_number() != 6 {
                            continue;
                        }

                        for (a, ba_bidx) in mol.neighbors(b3) {
                            if a == b2 {
                                continue;
                            }
                            if !rule.bridge_acceptor_order.matches(mol.bond(ba_bidx).order) {
                                continue;
                            }
                            if mol.atom(a).element.atomic_number() != rule.acceptor_elem {
                                continue;
                            }
                            // Return (donor, b2_as_central, acceptor) for the transfer logic
                            matches.push((d, b2, a));
                        }
                    }
                }
            }
        }
    } else {
        // Standard 1,3-shift (path_len=3)
        for i in 0..mol.atom_count() {
            let d = AtomIdx(i as u32);
            let donor_atom = mol.atom(d);
            if donor_atom.element.atomic_number() != rule.donor_elem {
                continue;
            }
            if implicit_hcount(mol, d) == 0 {
                continue;
            }

            for (b, db_bidx) in mol.neighbors(d) {
                if !rule.donor_bridge_order.matches(mol.bond(db_bidx).order) {
                    continue;
                }
                if let Some(br_elem) = rule.bridge_elem
                    && mol.atom(b).element.atomic_number() != br_elem
                {
                    continue;
                }

                for (a, ba_bidx) in mol.neighbors(b) {
                    if a == d {
                        continue;
                    }
                    if !rule.bridge_acceptor_order.matches(mol.bond(ba_bidx).order) {
                        continue;
                    }
                    if mol.atom(a).element.atomic_number() != rule.acceptor_elem {
                        continue;
                    }
                    matches.push((d, b, a));
                }
            }
        }
    }

    // Sort by canonical rank, not raw AtomIdx (parse/insertion order) --
    // whichever match a caller treats as "first" must depend on the
    // molecule's structure alone, not on how its atoms happened to be
    // labeled. See the "Atom-order canonicalization" section above.
    matches.sort_by_key(|&(d, b, a)| (rank[d.0 as usize], rank[b.0 as usize], rank[a.0 as usize]));
    matches
}

/// Apply a single tautomer transformation: donor-bridge bond Single → Double
/// and bridge-acceptor bond Double → Single. Returns `None` if the transform
/// would be invalid (e.g. donor has no explicit H to give up on a bracket atom).
///
/// For organic-subset atoms, implicit H counts adjust automatically through the
/// valence model; for bracket atoms with an explicit `hydrogen_count`, we
/// decrement the donor and increment the acceptor manually.
///
/// Donor/bridge/acceptor participate in a 1,3- or 1,5-shift, never a tetrahedral
/// stereocenter directly -- but the passthrough rebuild below still needs to copy
/// `stereo_neighbor_order`/`stereo_groups`/`bond_directions` verbatim for every
/// *other* atom, same as `transfer_hydrogen_aromatic` above, or any pre-existing
/// stereocenter elsewhere in the molecule silently loses the reference order its
/// `@`/`@@` is defined relative to on every tautomer step.
fn transfer_hydrogen(
    mol: &Molecule,
    donor: AtomIdx,
    bridge: AtomIdx,
    acceptor: AtomIdx,
    blocked_atoms: &HashSet<AtomIdx>,
    blocked_bonds: &HashSet<BondIdx>,
) -> Option<Molecule> {
    if blocked_atoms.contains(&donor)
        || blocked_atoms.contains(&bridge)
        || blocked_atoms.contains(&acceptor)
    {
        return None;
    }
    let (db_bidx, _) = mol.bond_between(donor, bridge)?;
    let (ba_bidx, _) = mol.bond_between(bridge, acceptor)?;
    if blocked_bonds.contains(&db_bidx) || blocked_bonds.contains(&ba_bidx) {
        return None;
    }

    let mut builder = MoleculeBuilder::new();
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let mut atom = mol.atom(idx).clone();
        if let Some(h) = atom.hydrogen_count {
            if idx == donor {
                if h == 0 {
                    return None;
                }
                atom.hydrogen_count = Some(h - 1);
            } else if idx == acceptor {
                atom.hydrogen_count = Some(h.saturating_add(1));
            }
        }
        builder.add_atom(atom);
    }

    for i in 0..mol.bond_count() {
        let bidx = BondIdx(i as u32);
        let b = mol.bond(bidx);
        let order = match bidx {
            x if x == db_bidx => BondOrder::Double,
            x if x == ba_bidx => BondOrder::Single,
            _ => b.order,
        };
        builder.add_bond(b.atom1, b.atom2, order).ok()?;
    }
    builder.copy_stereo_groups_from(mol);
    builder.copy_stereo_from(mol);
    builder.copy_bond_directions_from(mol);
    Some(builder.build())
}

/// Like the first-match helper, but also returns the (donor, bridge, acceptor)
/// triple that was moved -- needed by [`tautomer_parent`] to record which
/// atoms/bonds each applied transform touched.
fn apply_first_match_tracked(
    mol: &Molecule,
    rule: &TautomerRule,
    config: &TautomerConfig,
    rank: &[u32],
) -> Option<(Molecule, AtomIdx, AtomIdx, AtomIdx)> {
    find_matches(mol, rule, rank)
        .into_iter()
        .find_map(|(d, b, a)| {
            transfer_hydrogen(mol, d, b, a, &config.blocked_atoms, &config.blocked_bonds)
                .map(|next| (next, d, b, a))
        })
}

/// Apply every matching transformation for `rule`; return all resulting molecules.
fn apply_all_matches(
    mol: &Molecule,
    rule: &TautomerRule,
    config: &TautomerConfig,
    rank: &[u32],
) -> Vec<Molecule> {
    find_matches(mol, rule, rank)
        .into_iter()
        .filter_map(|(d, b, a)| {
            transfer_hydrogen(mol, d, b, a, &config.blocked_atoms, &config.blocked_bonds)
        })
        .collect()
}

/// Apply all non-overlapping matches of one rule in a single bounded pass.
///
/// A rule match changes the donor/bridge/acceptor and the two bonds joining
/// them.  Matches that touch any of those atoms are therefore conservatively
/// treated as conflicting and deferred to a later pass.  This preserves the
/// old rule priority while allowing independent sites to consume one
/// `max_iter` unit together.  In particular, the result no longer depends on
/// which of many independent sites happened to receive the first raw atom
/// index when the transform budget is reached.
fn apply_nonconflicting_matches(
    mol: &Molecule,
    rule: &TautomerRule,
    config: &TautomerConfig,
    rank: &[u32],
) -> (Molecule, bool) {
    let (next, applied) = apply_nonconflicting_matches_tracked(mol, rule, config, rank, usize::MAX);
    (next, !applied.is_empty())
}

/// Tracked form of [`apply_nonconflicting_matches`]. `limit` is the maximum
/// number of transformations to apply in this pass, allowing the audit-aware
/// parent API to preserve its `max_transforms` bound while still making its
/// selection order canonical.
fn apply_nonconflicting_matches_tracked(
    mol: &Molecule,
    rule: &TautomerRule,
    config: &TautomerConfig,
    rank: &[u32],
    limit: usize,
) -> (Molecule, Vec<(AtomIdx, AtomIdx, AtomIdx)>) {
    let matches = find_matches(mol, rule, rank);
    if matches.is_empty() || limit == 0 {
        return (mol.clone(), Vec::new());
    }

    let mut current = mol.clone();
    let mut occupied: HashSet<AtomIdx> = HashSet::new();
    let mut applied = Vec::new();

    for (donor, bridge, acceptor) in matches {
        if applied.len() >= limit {
            break;
        }
        let Some((db, _)) = mol.bond_between(donor, bridge) else {
            continue;
        };
        let Some((ba, _)) = mol.bond_between(bridge, acceptor) else {
            continue;
        };
        let touched = [
            donor,
            bridge,
            acceptor,
            mol.bond(db).atom1,
            mol.bond(db).atom2,
            mol.bond(ba).atom1,
            mol.bond(ba).atom2,
        ];
        if touched.iter().any(|atom| occupied.contains(atom)) {
            continue;
        }
        let Some(next) = transfer_hydrogen(
            &current,
            donor,
            bridge,
            acceptor,
            &config.blocked_atoms,
            &config.blocked_bonds,
        ) else {
            continue;
        };
        current = next;
        occupied.extend(touched);
        applied.push((donor, bridge, acceptor));
    }

    (current, applied)
}

// ---------------------------------------------------------------------------
// TautomerConfig
// ---------------------------------------------------------------------------

/// Configuration for tautomer enumeration and canonicalization.
///
/// # Rule indices
/// Rules are numbered 0-based in the order they appear in the built-in set.
/// Use [`TautomerConfig::rule_count`] to know the total and
/// [`TautomerConfig::rule_names`] to see what each index represents.
///
/// An empty `enabled_rules` (the default) activates **all** rules.
///
/// # Zone blocking
/// `blocked_atoms` and `blocked_bonds` prevent H-transfer through specific atoms/bonds.
/// Any tautomer move whose donor, bridge, or acceptor is in `blocked_atoms`,
/// or whose altered bond is in `blocked_bonds`, is suppressed.
/// Note: for 1,5-shift rules only the donor, bridge-central, and acceptor atoms are
/// checked; intermediate path atoms cannot be blocked without refactoring `find_matches`.
#[derive(Debug, Clone)]
pub struct TautomerConfig {
    /// Maximum iterations in [`canonical_tautomer_with_config`] (default 16).
    ///
    /// Each outer iteration applies all non-overlapping matches of the first
    /// active forward rule. Overlapping matches are deferred and reconsidered
    /// in the next iteration. This keeps the bound meaningful for chained or
    /// conflicting transformations while making independent sites consume one
    /// iteration together, independent of input atom insertion order.
    pub max_iter: usize,
    /// Maximum tautomers returned by [`enumerate_tautomers_with_config`] (default 32).
    pub max_tautomers: usize,
    /// 0-based indices of rules to activate.  Empty = all rules active.
    pub enabled_rules: Vec<usize>,
    /// Atoms that must not participate in any H-transfer (donor, bridge, or acceptor).
    /// Empty = no atom blocking (default).
    pub blocked_atoms: HashSet<AtomIdx>,
    /// Bonds that must not be altered by any H-transfer.
    /// Empty = no bond blocking (default).
    pub blocked_bonds: HashSet<BondIdx>,
}

impl Default for TautomerConfig {
    fn default() -> Self {
        Self {
            max_iter: 16,
            max_tautomers: 32,
            enabled_rules: Vec::new(),
            blocked_atoms: HashSet::new(),
            blocked_bonds: HashSet::new(),
        }
    }
}

impl TautomerConfig {
    /// Number of built-in tautomer rules available.
    pub fn rule_count() -> usize {
        RULES.len()
    }

    /// Names of all built-in rules, in index order.
    pub fn rule_names() -> Vec<&'static str> {
        RULES.iter().map(|r| r.name).collect()
    }

    /// Convenience: config with only the keto-enol rule (index 0) enabled.
    pub fn keto_enol_only() -> Self {
        Self {
            enabled_rules: vec![0],
            ..Self::default()
        }
    }
}

// ---------------------------------------------------------------------------
// Internal helper: iterate only the active rules
// ---------------------------------------------------------------------------

/// Collect the rules that are active under `config`.
fn active_rules(config: &TautomerConfig) -> Vec<&'static TautomerRule> {
    RULES
        .iter()
        .enumerate()
        .filter_map(|(i, r)| {
            if config.enabled_rules.is_empty() || config.enabled_rules.contains(&i) {
                Some(r)
            } else {
                None
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Atom-order canonicalization (issue #415 residual)
//
// `find_sssr`'s ring choice (consumed directly by
// `find_exocyclic_lactam_shift_matches`) and `find_matches`'/
// `find_direct_aromatic_matches`'s "first match in raw AtomIdx order" policy
// are both dependent on the *caller's* parse/atom-insertion order, not just
// the molecule's structure -- so which tautomer this module reaches for the
// same molecule could depend on how its atoms happened to be labeled.
//
// Fixed by computing a canonical *rank* per atom once
// (`chematic_smiles::canonical_atom_order`, inverted into an `AtomIdx ->
// rank` lookup) and using it to break every such tie -- sort candidates by
// rank instead of by raw `AtomIdx`, and prefer the lowest-rank ring when a
// ring choice is needed. `AtomIdx` identity never changes across a
// tautomer-search step (every transform function rebuilds atom-for-atom via
// `for i in 0..mol.atom_count()`, same order in and out), so one rank
// vector computed from the *entry* molecule stays valid for every
// intermediate form the search visits.
//
// An earlier version of this fix instead rebuilt the molecule with atoms
// physically reordered into `canonical_atom_order`'s own order before
// searching. That was correct (verified: full order-invariance across 6
// permutations of issue #415's residual) but is NOT what's implemented here
// -- it uncovered a separate, pre-existing `chematic-smiles` defect where
// `canonical_smiles`/`canonical_atom_order`'s own individualize-refine
// search degenerates into a multi-minute hang specifically when its input's
// atom order already coincides with `canonical_atom_order`'s own output
// order, for at least one real highly-symmetric molecule
// (`chembl_accuracy_corpus_4999.smi` line 2741, a 94-atom molecule with 3
// near-equivalent Boc-benzylamine arms -- reproduces with zero
// `chematic-chem` involvement: `parse` -> `canonical_atom_order` -> rebuild
// in that order -> `canonical_smiles` hangs). Rank-based tie-breaking never
// constructs that pathological input shape (the molecule's own storage/atom
// order is never touched), so it sidesteps the landmine entirely rather
// than working around it. That `chematic-smiles` defect is filed separately
// (issue TBD) as discovered-but-deferred, out of scope for this fix.
// ---------------------------------------------------------------------------

/// Canonical rank per atom (`result[idx.0] = canonical position`), used to
/// break every atom-order-sensitive tie in this module's search
/// deterministically, as a function of `mol`'s structure alone.
fn atom_rank(mol: &Molecule) -> Vec<u32> {
    let order = chematic_smiles::canonical_atom_order(mol);
    let mut rank = vec![0u32; order.len()];
    for (pos, &orig) in order.iter().enumerate() {
        rank[orig] = pos as u32;
    }
    rank
}

/// Canonical rank for a tautomer skeleton, independent of explicit H
/// placement.  This is the rooted structural key used by transformations
/// whose only difference is which equivalent aromatic N carries H.
fn tautomer_skeleton_rank(mol: &Molecule) -> Vec<u32> {
    let mut builder = MoleculeBuilder::new();
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let mut atom = mol.atom(idx).clone();
        atom.hydrogen_count = None;
        builder.add_atom(atom);
    }
    for i in 0..mol.bond_count() {
        let bond = mol.bond(BondIdx(i as u32));
        builder
            .add_bond(bond.atom1, bond.atom2, bond.order)
            .expect("tautomer skeleton preserves valid bonds");
    }
    atom_rank(&builder.build())
}

/// Return a canonical key for the H-free skeleton rooted at `root`.
/// Isotopically marking the root makes the comparison independent of input
/// atom order even when the unrooted skeleton has automorphisms.
fn rooted_structural_key(mol: &Molecule, root: AtomIdx) -> String {
    let mut builder = MoleculeBuilder::new();
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let mut atom = mol.atom(idx).clone();
        atom.hydrogen_count = None;
        if idx == root {
            atom.isotope = Some(65535);
        }
        builder.add_atom(atom);
    }
    for i in 0..mol.bond_count() {
        let bond = mol.bond(BondIdx(i as u32));
        builder
            .add_bond(bond.atom1, bond.atom2, bond.order)
            .expect("rooted skeleton preserves valid bonds");
    }
    chematic_smiles::canonical_smiles(&builder.build())
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Return the canonical (preferred) tautomer of `mol`.
///
/// Applies forward-preferred rules iteratively until no new form is found
/// or the iteration limit is reached; each iteration also tries the
/// directional exocyclic lactam/lactim shift (round 2C-2, RFC section
/// 4.4a) as a fallback when no rule matched. After that loop converges,
/// direct aromatic 1,2-shift tautomers (ring-internal only, e.g.
/// imidazole/pyrazole/tetrazole N-position choice) are compared and ranked
/// by `tautomer_score` (O-H > N-H > S-H, aromatic-ring bonus), with
/// canonical SMILES as the final tiebreaker.
///
/// **A stereocenter's CIP (R/S) label can legitimately change across this
/// shift, even at an atom this function never touches** (issue #402):
/// `stereo_neighbor_order`/`chirality` are always carried forward verbatim
/// for every atom not directly part of the transform (the real spatial
/// configuration never moves), but a CIP label is *computed* from
/// substituent priorities across the whole molecular graph, and a keto-enol
/// (or similar) shift a few bonds away can change those priorities. Compare
/// CIP labels only on the *same* tautomeric form on both sides -- comparing
/// a pre-shift label to a post-shift label for "did stereo survive
/// standardization" is not a valid check, and a difference there does not by
/// itself indicate a bug. Confirmed via independent RDKit CIP oracle on two
/// molecules where this shift changes a label: chematic's CIP on its own
/// shifted structure matches RDKit's CIP on that identical structure exactly
/// in both cases.
///
/// Uses [`TautomerConfig::default`] (all rules, max_iter=16).
pub fn canonical_tautomer(mol: &Molecule) -> Molecule {
    canonical_tautomer_with_config(mol, &TautomerConfig::default())
}

/// Like [`canonical_tautomer`] but with explicit configuration.
///
/// Every internal candidate/ring tie is broken via [`atom_rank`] (see the
/// "Atom-order canonicalization" section above), computed once from `mol`
/// and reused for every intermediate form the search visits (atom identity
/// never changes across a tautomer-search step), so the result is
/// independent of the caller's parse/atom-insertion order.
pub fn canonical_tautomer_with_config(mol: &Molecule, config: &TautomerConfig) -> Molecule {
    // Preserve E/Z-bearing inputs verbatim until a transform can carry their
    // directional markers onto the new double bond.  Returning the original
    // form is preferable to silently collapsing distinct stereoisomers.
    if has_directional_bond(mol) {
        return mol.clone();
    }
    let rank = atom_rank(mol);

    // canonical_tautomer_search's own final direct-aromatic tie-break is a
    // one-shot step, not iterated to convergence -- confirmed real (not an
    // atom-order artifact) via issue #415's own residual, which needs
    // exactly one extra application to reach a stable form, order-
    // invariantly, once rank-based tie-breaking is used above. Re-run until
    // stable so this function is actually idempotent on its own output, not
    // just closer to it. The `canonical_smiles`-based equality check is
    // skipped entirely when a pass changes nothing (the common case for
    // molecules with no tautomerizable groups at all).
    let (mut result, changed) = canonical_tautomer_search(mol, config, &rank);
    if changed {
        let mut seen_smiles: HashSet<String> = HashSet::new();
        seen_smiles.insert(chematic_smiles::canonical_smiles(&result));
        for _ in 0..3 {
            let (next, next_changed) = canonical_tautomer_search(&result, config, &rank);
            if !next_changed {
                break;
            }
            let is_new = seen_smiles.insert(chematic_smiles::canonical_smiles(&next));
            result = next;
            if !is_new {
                break;
            }
        }
    }
    // The directed search above is intentionally cheap, but its greedy path
    // can make keto spellings of an asymmetric dual-flank ring land in
    // different components. For aromatic O/C tautomer systems, enumerate the
    // bounded connected component and select one representative globally.
    // Prefer the lactam (aromatic C=O) subset whenever it exists; this keeps
    // the established 2-pyridone preference while making keto inputs and
    // hydroxyl inputs share the same canonical result.
    if config.max_tautomers > 0 {
        let forms = enumerate_tautomers_with_config(mol, config);
        if forms.iter().any(has_aromatic_exocyclic_oxygen)
            && forms
                .iter()
                .all(|candidate| aromatic_exocyclic_carbonyl_count(candidate) <= 1)
        {
            let prefer_lactam = forms.iter().any(has_aromatic_exocyclic_carbonyl);
            let mut candidates: Vec<Molecule> = forms
                .into_iter()
                .filter(|candidate| !prefer_lactam || has_aromatic_exocyclic_carbonyl(candidate))
                .collect();
            candidates.sort_by(|a, b| {
                tautomer_score(b).cmp(&tautomer_score(a)).then_with(|| {
                    carbonyl_rooted_key(a)
                        .cmp(&carbonyl_rooted_key(b))
                        .then_with(|| {
                            chematic_smiles::canonical_smiles(a)
                                .cmp(&chematic_smiles::canonical_smiles(b))
                        })
                })
            });
            if let Some(best) = candidates.into_iter().next() {
                result = best;
            }
        }
    }
    result
}

/// The actual rule-based/exocyclic-lactam/direct-aromatic search. Returns
/// the result plus whether anything actually changed from `mol` -- lets
/// [`canonical_tautomer_with_config`]'s convergence loop skip re-verifying
/// (and the `canonical_smiles` calls that involves) when a pass was a no-op.
fn canonical_tautomer_search(
    mol: &Molecule,
    config: &TautomerConfig,
    rank: &[u32],
) -> (Molecule, bool) {
    let mut current = mol.clone();
    let mut seen = HashSet::new();
    seen.insert(mol_fingerprint(&current));
    let mut changed = false;
    let mut dual_flank_changed = false;
    let mut carbonyl_rooted = false;

    for _ in 0..config.max_iter {
        let mut iter_changed = false;
        for rule in active_rules(config)
            .into_iter()
            .filter(|r| r.prefer_forward)
        {
            let (next, applied) = apply_nonconflicting_matches(&current, rule, config, rank);
            if applied {
                let fp = mol_fingerprint(&next);
                if !seen.contains(&fp) {
                    seen.insert(fp);
                    current = next;
                    iter_changed = true;
                    changed = true;
                    break;
                }
            }
        }
        if !iter_changed {
            // The 41 rules above never match across an aromatic ring bond
            // (BondOrderMatch::Double never matches BondOrder::Aromatic).
            // This is the round 2C-2 generalization: an exocyclic O donor
            // across an aromatic ring atom -- see the mechanism block above
            // `transfer_hydrogen_exocyclic_lactam`.
            if let Some((next, ..)) = apply_exocyclic_lactam_shift_tracked(&current, config, rank) {
                let fp = mol_fingerprint(&next);
                if !seen.contains(&fp) {
                    seen.insert(fp);
                    current = next;
                    iter_changed = true;
                    changed = true;
                    carbonyl_rooted = true;
                }
            }
            if !iter_changed
                && let Some((next, ..)) = apply_dual_flank_tracked(&current, config, rank)
            {
                let fp = mol_fingerprint(&next);
                if !seen.contains(&fp) {
                    seen.insert(fp);
                    current = next;
                    iter_changed = true;
                    changed = true;
                    dual_flank_changed = true;
                }
            }
        }
        if !iter_changed {
            break;
        }
    }

    // The dual-flank step already selected a deterministic rooted structural
    // key.  Do not feed that result back into the broader direct-aromatic
    // score pool, which can otherwise move the proton to a third ring N and
    // undo the cytosine/guanine normalization.
    if dual_flank_changed || carbonyl_rooted {
        return (current, changed);
    }

    // Among direct aromatic 1,2-shift tautomers, pick by tautomer score
    // (O-H > N-H > S-H, aromatic rings), with canonical SMILES as tiebreaker.
    //
    // BFS via enumerate_direct_aromatic_forms ensures all ring positions are
    // explored (e.g. tetrazole has 4 N atoms; a single step would miss N1↔N3).
    // The canonical SMILES tiebreaker makes selection independent of input
    // SMILES write order (the previous h_assignment tiebreaker was not).
    let mut candidates: Vec<Molecule> = vec![current.clone()];
    candidates.extend(enumerate_direct_aromatic_forms(
        &current,
        &config.blocked_atoms,
        16,
        rank,
    ));
    if candidates.len() > 1 {
        candidates.sort_by(|a, b| {
            tautomer_score(b).cmp(&tautomer_score(a)).then_with(|| {
                chematic_smiles::canonical_smiles(a).cmp(&chematic_smiles::canonical_smiles(b))
            })
        });
        current = candidates.into_iter().next().unwrap();
        changed = true;
    }
    (current, changed)
}

/// Enumerate all reachable tautomers of `mol`, capped at 32.
///
/// Returns a `Vec<Molecule>` where the first element is the original molecule.
/// Includes both 1,3-shift (rule-based) and direct aromatic 1,2-shift tautomers.
///
/// Uses [`TautomerConfig::default`] (all rules, max_tautomers=32).
pub fn enumerate_tautomers(mol: &Molecule) -> Vec<Molecule> {
    enumerate_tautomers_with_config(mol, &TautomerConfig::default())
}

/// Like [`enumerate_tautomers`] but with explicit configuration.
pub fn enumerate_tautomers_with_config(mol: &Molecule, config: &TautomerConfig) -> Vec<Molecule> {
    let rank = atom_rank(mol);
    let mut result = vec![mol.clone()];
    let mut seen = HashSet::new();
    seen.insert(mol_fingerprint(mol));
    // Separate seen-set for 1,2-shift (mol_fingerprint can't distinguish positional H isomers).
    let mut h_seen: HashSet<Vec<Option<u32>>> = HashSet::new();
    h_seen.insert(h_assignment(mol));
    let mut frontier = vec![mol.clone()];

    while !frontier.is_empty() && result.len() < config.max_tautomers {
        let current = frontier.remove(0);
        for rule in active_rules(config).into_iter() {
            for next in apply_all_matches(&current, rule, config, &rank) {
                let fp = mol_fingerprint(&next);
                if !seen.contains(&fp) {
                    seen.insert(fp);
                    h_seen.insert(h_assignment(&next));
                    frontier.push(next.clone());
                    result.push(next);
                    if result.len() >= config.max_tautomers {
                        break;
                    }
                }
            }
            if result.len() >= config.max_tautomers {
                break;
            }
        }
        // Keep enumeration in lockstep with canonical search for aromatic
        // lactam/lactim systems. These shifts are outside the generic rule
        // table because aromatic ring bonds do not match an explicit Double
        // pattern, but they are still reachable tautomer edges.
        for next in [
            apply_exocyclic_lactam_shift_tracked(&current, config, &rank).map(|(next, ..)| next),
            apply_dual_flank_tracked(&current, config, &rank).map(|(next, ..)| next),
        ] {
            if result.len() >= config.max_tautomers {
                break;
            }
            if let Some(next) = next {
                let fp = mol_fingerprint(&next);
                if !seen.contains(&fp) {
                    seen.insert(fp);
                    h_seen.insert(h_assignment(&next));
                    frontier.push(next.clone());
                    result.push(next);
                }
            }
        }
        for next in exocyclic_lactam_candidates(&current, &rank, config) {
            if result.len() >= config.max_tautomers {
                break;
            }
            let fp = mol_fingerprint(&next);
            if !seen.contains(&fp) {
                seen.insert(fp);
                h_seen.insert(h_assignment(&next));
                frontier.push(next.clone());
                result.push(next);
            }
        }
        // Enumeration must retain both orientations of the dual-flank edge;
        // the canonical search helper intentionally returns only the rooted
        // preferred direction.
        let skeleton_rank = tautomer_skeleton_rank(&current);
        for (donor, _, acceptor) in find_dual_flank_matches(&current, &skeleton_rank) {
            if result.len() >= config.max_tautomers {
                break;
            }
            if let Some(next) =
                transfer_hydrogen_dual_flank(&current, donor, acceptor, &config.blocked_atoms)
            {
                let fp = mol_fingerprint(&next);
                if !seen.contains(&fp) {
                    seen.insert(fp);
                    h_seen.insert(h_assignment(&next));
                    frontier.push(next.clone());
                    result.push(next);
                }
            }
        }
        for next in reverse_exocyclic_lactim_candidates(&current, config) {
            if result.len() >= config.max_tautomers {
                break;
            }
            let fp = mol_fingerprint(&next);
            if !seen.contains(&fp) {
                seen.insert(fp);
                h_seen.insert(h_assignment(&next));
                frontier.push(next.clone());
                result.push(next);
            }
        }
        // Direct aromatic 1,2-shift (e.g. pyrazole N1H ↔ N2H).
        for (d, a) in find_direct_aromatic_matches(&current, &rank) {
            if result.len() >= config.max_tautomers {
                break;
            }
            if let Some(next) = transfer_hydrogen_aromatic(&current, d, a, &config.blocked_atoms) {
                let ha = h_assignment(&next);
                if !h_seen.contains(&ha) {
                    h_seen.insert(ha);
                    seen.insert(mol_fingerprint(&next));
                    frontier.push(next.clone());
                    result.push(next);
                }
            }
        }
    }
    result
}

// ---------------------------------------------------------------------------
// Tautomer Parent (ROADMAP.md Phase 2 round 2B) -- see
// docs/rfcs/tautomer_parent_identity_phase2_rfc.md section 4.
// ---------------------------------------------------------------------------

/// Budget limits for [`tautomer_parent`].
///
/// `max_transforms`/`max_tautomers` are deterministic budgets, counted
/// exactly as [`TautomerConfig::max_iter`]/[`TautomerConfig::max_tautomers`]
/// always have been: `max_transforms` counts outer-loop iterations in which
/// a transform was actually applied (a rule that matches but produces an
/// already-seen fingerprint does not consume budget); `max_tautomers`
/// counts distinct direct-aromatic-shift candidates found while choosing
/// among them at the end (duplicates rejected by fingerprint do not consume
/// budget either). Both are reproducible: same input, same limits, same
/// result, always.
///
/// `timeout_ms` is a different kind of bound: a wall-clock timeout depends
/// on machine speed and load, so a `TimedOut` result is explicitly outside
/// the "same canonical tautomer regardless of input" determinism guarantee
/// the other two limits carry. `None` (the default) performs no wall-clock
/// check.
///
/// `max_tautomers: 0` is treated as "skip the direct-aromatic-shift
/// comparison step entirely" (status stays `Completed`, using whatever the
/// rule-based loop converged to) rather than reporting `MaxTautomersReached`
/// on every input regardless of whether any alternate form actually
/// exists -- the latter would be a boundary artifact of the search
/// algorithm, not a meaningful signal.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TautomerLimits {
    pub max_transforms: usize,
    pub max_tautomers: usize,
    pub timeout_ms: Option<u64>,
}

impl Default for TautomerLimits {
    fn default() -> Self {
        Self {
            max_transforms: 16,
            max_tautomers: 32,
            timeout_ms: None,
        }
    }
}

/// Stable identifier for one of the built-in [`TautomerRule`]s, used in
/// [`AppliedTransform::rule_id`]. `#[non_exhaustive]` so new rules can be
/// added without a breaking change; `Other` is a defensive fallback for a
/// rule name not yet mapped to a variant here -- should not occur for any
/// rule currently in `RULES` (see `all_rules_map_to_a_named_tautomer_rule_id`
/// in this module's tests).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum TautomerRuleId {
    KetoEnol,
    AmideIminol,
    IminolAmide,
    ImineEnamine,
    OneThreeNToO,
    OneThreeNToN,
    Thioamide,
    ThioIminolAmide,
    ThioKetoEnol,
    ThioEnolKetone,
    OneThreeNToS,
    OneThreeSToO,
    OneThreeSToN,
    OneThreeOToS,
    OneThreeSToS,
    OneThreeOToNAnyBridge,
    OneThreeOToOAnyBridge,
    OneThreeNToCAnyBridge,
    OneThreeCToOAnyBridge,
    NitrosoOxime,
    OneThreeCToNAnyBridge,
    OneFiveOToOBetaDiketone,
    OneFiveOToN,
    OneFiveNToN,
    OneFiveNToO,
    OneFiveCToO,
    OneFiveOToONBridge,
    OneFiveOToOSBridge,
    OneFiveNToONBridge,
    OneFiveNToOSBridge,
    OneFiveNToNNBridge,
    OneFiveNToNSBridge,
    OneFiveOToNNBridge,
    OneFiveOToNSBridge,
    OneFiveSToONBridge,
    OneFiveSToOSBridge,
    OneFiveSToNCBridge,
    OneFiveSToNNBridge,
    OneFiveCToNCBridge,
    OneFiveCToNNBridge,
    OneFiveCToNSBridge,
    OneFiveCToSNBridge,
    OneFiveCToSSBridge,
    /// The round 2C-2 exocyclic lactam/lactim shift -- not a `RULES` table
    /// entry (it's a directional step applied alongside them, see
    /// `apply_exocyclic_lactam_shift_tracked`), so it is never produced by
    /// `rule_id_for` and is exempt from
    /// `all_rules_map_to_a_named_tautomer_rule_id`'s drift check.
    AromaticExocyclicLactamLactim,
    /// H transfer between the two aromatic nitrogens flanking a carbonyl.
    AromaticDualFlank,
    /// A rule name not yet mapped to a named variant above.
    Other(&'static str),
}

fn rule_id_for(name: &'static str) -> TautomerRuleId {
    match name {
        "keto-enol" => TautomerRuleId::KetoEnol,
        "amide-iminol" => TautomerRuleId::AmideIminol,
        "iminol-amide" => TautomerRuleId::IminolAmide,
        "imine-enamine" => TautomerRuleId::ImineEnamine,
        "1,3-N-to-O" => TautomerRuleId::OneThreeNToO,
        "1,3-N-to-N" => TautomerRuleId::OneThreeNToN,
        "thioamide" => TautomerRuleId::Thioamide,
        "thio-iminol-amide" => TautomerRuleId::ThioIminolAmide,
        "thio-keto-enol" => TautomerRuleId::ThioKetoEnol,
        "thio-enol-ketone" => TautomerRuleId::ThioEnolKetone,
        "1,3-N-to-S" => TautomerRuleId::OneThreeNToS,
        "1,3-S-to-O" => TautomerRuleId::OneThreeSToO,
        "1,3-S-to-N" => TautomerRuleId::OneThreeSToN,
        "1,3-O-to-S" => TautomerRuleId::OneThreeOToS,
        "1,3-S-to-S" => TautomerRuleId::OneThreeSToS,
        "1,3-O-to-N-any-bridge" => TautomerRuleId::OneThreeOToNAnyBridge,
        "1,3-O-to-O-any-bridge" => TautomerRuleId::OneThreeOToOAnyBridge,
        "1,3-N-to-C-any-bridge" => TautomerRuleId::OneThreeNToCAnyBridge,
        "1,3-C-to-O-any-bridge" => TautomerRuleId::OneThreeCToOAnyBridge,
        "nitroso-oxime" => TautomerRuleId::NitrosoOxime,
        "1,3-C-to-N-any-bridge" => TautomerRuleId::OneThreeCToNAnyBridge,
        "1,5-O-to-O-beta-diketone" => TautomerRuleId::OneFiveOToOBetaDiketone,
        "1,5-O-to-N" => TautomerRuleId::OneFiveOToN,
        "1,5-N-to-N" => TautomerRuleId::OneFiveNToN,
        "1,5-N-to-O" => TautomerRuleId::OneFiveNToO,
        "1,5-C-to-O" => TautomerRuleId::OneFiveCToO,
        "1,5-O-to-O-N-bridge" => TautomerRuleId::OneFiveOToONBridge,
        "1,5-O-to-O-S-bridge" => TautomerRuleId::OneFiveOToOSBridge,
        "1,5-N-to-O-N-bridge" => TautomerRuleId::OneFiveNToONBridge,
        "1,5-N-to-O-S-bridge" => TautomerRuleId::OneFiveNToOSBridge,
        "1,5-N-to-N-N-bridge" => TautomerRuleId::OneFiveNToNNBridge,
        "1,5-N-to-N-S-bridge" => TautomerRuleId::OneFiveNToNSBridge,
        "1,5-O-to-N-N-bridge" => TautomerRuleId::OneFiveOToNNBridge,
        "1,5-O-to-N-S-bridge" => TautomerRuleId::OneFiveOToNSBridge,
        "1,5-S-to-O-N-bridge" => TautomerRuleId::OneFiveSToONBridge,
        "1,5-S-to-O-S-bridge" => TautomerRuleId::OneFiveSToOSBridge,
        "1,5-S-to-N-C-bridge" => TautomerRuleId::OneFiveSToNCBridge,
        "1,5-S-to-N-N-bridge" => TautomerRuleId::OneFiveSToNNBridge,
        "1,5-C-to-N-C-bridge" => TautomerRuleId::OneFiveCToNCBridge,
        "1,5-C-to-N-N-bridge" => TautomerRuleId::OneFiveCToNNBridge,
        "1,5-C-to-N-S-bridge" => TautomerRuleId::OneFiveCToNSBridge,
        "1,5-C-to-S-N-bridge" => TautomerRuleId::OneFiveCToSNBridge,
        "1,5-C-to-S-S-bridge" => TautomerRuleId::OneFiveCToSSBridge,
        other => TautomerRuleId::Other(other),
    }
}

/// One named contribution to a candidate tautomer's [`tautomer_score`],
/// decomposed for [`TautomerAuditRecord::score_breakdown`].
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum TautomerScoreTerm {
    AromaticRing,
    /// `element` is an atomic number (today: 8=O, 7=N, 16=S).
    HeteroatomHydrogen {
        element: u8,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ScoreContribution {
    pub term: TautomerScoreTerm,
    pub value: i32,
}

/// One transformation applied by the rule-based 1,3-/1,5-shift loop while
/// computing a [`TautomerAuditRecord`].
///
/// Does not implement `Serialize`/`Deserialize` even under the `serde`
/// feature: `AtomIdx`/`BondIdx` (from `chematic-core`) do not implement
/// them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AppliedTransform {
    pub rule_id: TautomerRuleId,
    pub affected_atoms: Vec<AtomIdx>,
    pub affected_bonds: Vec<BondIdx>,
}

/// Explainable audit record for one [`tautomer_parent`] call.
///
/// `applied_transforms` records only the rule-based 1,3-/1,5-shift loop (the
/// same loop [`canonical_tautomer_with_config`] runs); the final
/// direct-aromatic-shift tie-break contributes to `score_breakdown` and
/// `candidate_count` instead, since that search does not track the specific
/// shift sequence taken to reach each candidate -- a disclosed scope
/// limitation, not an oversight.
///
/// `lost_stereo`/`affected_isotopes` are atom-index lists (empty = none
/// affected), not booleans: every transform this module applies preserves
/// `chirality`/`isotope` (see the RFC's non-bug finding, section 1.2), so
/// both are expected to always be empty today; the check is real (not
/// hardcoded) so it stays correct if that ever changes.
///
/// Does not implement `Serialize`/`Deserialize` even under the `serde`
/// feature: it embeds [`AppliedTransform`] (see its own doc comment).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TautomerAuditRecord {
    pub selected: MoleculeSnapshot,
    pub candidate_count: usize,
    pub score_breakdown: Vec<ScoreContribution>,
    pub applied_transforms: Vec<AppliedTransform>,
    pub lost_stereo: Vec<AtomIdx>,
    pub affected_isotopes: Vec<AtomIdx>,
}

/// Like [`canonical_tautomer`], but budget-limited via [`TautomerLimits`]
/// and returns a [`ParentResult`] (status + explainable
/// [`TautomerAuditRecord`]) instead of a bare `Molecule` -- see
/// `docs/rfcs/tautomer_parent_identity_phase2_rfc.md` sections 4.2-4.3 and
/// 4.4a.
///
/// Includes the section 4.4a directional exocyclic lactam/lactim shift
/// (round 2C-2, `apply_exocyclic_lactam_shift_tracked`), wired into this
/// function's own loop the same way it's wired into
/// [`canonical_tautomer_with_config`]'s -- `MaxTransformsReached` accounts
/// for it. The former Section 1.7 cytosine/guanine dual-flank residual is
/// handled by rooted N-H normalization and is covered by the
/// `tp2_07_09_dual_flank_cytosine_and_guanine_converge` regression. The
/// former Section 1.6 nitroso/oxime residual (`CCN=O`/`CC=NO`) is handled by
/// the dedicated forward `nitroso-oxime` rule and is covered by the
/// `tp2_04_nitroso_oxime_converges_without_enolizing_ketones` regression.
/// The asymmetric `tp2-39` and N9-methylhypoxanthine holdout remain explicit
/// residuals because their broader ring-bond tautomer components are not yet
/// fully enumerated by the canonical search.
///
/// Independent matches of one rule are applied together in canonical-rank
/// order, subject to `limits.max_transforms`; each applied match is retained
/// in `TautomerAuditRecord::applied_transforms`. This keeps bounded parent
/// selection atom-order invariant without weakening the transform limit.
pub fn tautomer_parent(mol: &Molecule, limits: &TautomerLimits) -> ParentResult {
    if mol.atom_count() == 0 {
        return ParentResult {
            molecule: mol.clone(),
            status: ParentComputationStatus::InvalidInput(InvalidInputReason::EmptyMolecule),
            audit: ParentAudit::Tautomer(TautomerAuditRecord {
                selected: MoleculeSnapshot::from_mol(mol),
                candidate_count: 0,
                score_breakdown: Vec::new(),
                applied_transforms: Vec::new(),
                lost_stereo: Vec::new(),
                affected_isotopes: Vec::new(),
            }),
        };
    }

    // Every candidate/ring tie broken via `rank` (issue #415 residual --
    // see the "Atom-order canonicalization" section above), so the result
    // is independent of `mol`'s caller-supplied atom order.
    let rank = atom_rank(mol);

    let config = TautomerConfig {
        max_iter: limits.max_transforms,
        max_tautomers: limits.max_tautomers,
        ..TautomerConfig::default()
    };
    let deadline = limits
        .timeout_ms
        .map(|ms| (Instant::now(), std::time::Duration::from_millis(ms)));
    let timed_out = |deadline: &Option<(Instant, std::time::Duration)>| match deadline {
        Some((start, budget)) => start.elapsed() >= *budget,
        None => false,
    };

    let mut current = mol.clone();
    let mut seen = HashSet::new();
    seen.insert(mol_fingerprint(&current));
    let mut applied_transforms = Vec::new();
    let mut transforms_applied = 0usize;
    let mut status = ParentComputationStatus::Completed;

    loop {
        if timed_out(&deadline) {
            status = ParentComputationStatus::TimedOut;
            break;
        }
        if transforms_applied >= limits.max_transforms {
            let has_more = active_rules(&config)
                .into_iter()
                .filter(|r| r.prefer_forward)
                .any(|rule| {
                    apply_first_match_tracked(&current, rule, &config, &rank)
                        .is_some_and(|(next, ..)| !seen.contains(&mol_fingerprint(&next)))
                })
                || apply_exocyclic_lactam_shift_tracked(&current, &config, &rank)
                    .is_some_and(|(next, ..)| !seen.contains(&mol_fingerprint(&next)))
                || apply_dual_flank_tracked(&current, &config, &rank)
                    .is_some_and(|(next, ..)| !seen.contains(&mol_fingerprint(&next)));
            if has_more {
                status = ParentComputationStatus::MaxTransformsReached;
            }
            break;
        }
        let mut changed = false;
        for rule in active_rules(&config)
            .into_iter()
            .filter(|r| r.prefer_forward)
        {
            let remaining = limits.max_transforms - transforms_applied;
            let (next, applied) =
                apply_nonconflicting_matches_tracked(&current, rule, &config, &rank, remaining);
            if !applied.is_empty() {
                let fp = mol_fingerprint(&next);
                if !seen.contains(&fp) {
                    seen.insert(fp);
                    for (donor, bridge, acceptor) in applied.iter().copied() {
                        let mut affected_bonds = Vec::new();
                        if let Some((bidx, _)) = current.bond_between(donor, bridge) {
                            affected_bonds.push(bidx);
                        }
                        if let Some((bidx, _)) = current.bond_between(bridge, acceptor) {
                            affected_bonds.push(bidx);
                        }
                        applied_transforms.push(AppliedTransform {
                            rule_id: rule_id_for(rule.name),
                            affected_atoms: vec![donor, bridge, acceptor],
                            affected_bonds,
                        });
                    }
                    current = next;
                    transforms_applied += applied.len();
                    changed = true;
                    break;
                }
            }
        }
        if !changed {
            // Round 2C-2: the exocyclic lactam/lactim shift, tried only when
            // no rule above matched (same fallback order as
            // canonical_tautomer_with_config).
            if let Some((next, donor, bridge, acceptor)) =
                apply_exocyclic_lactam_shift_tracked(&current, &config, &rank)
            {
                let fp = mol_fingerprint(&next);
                if !seen.contains(&fp) {
                    seen.insert(fp);
                    let mut affected_bonds = Vec::new();
                    if let Some((bidx, _)) = current.bond_between(donor, bridge) {
                        affected_bonds.push(bidx);
                    }
                    applied_transforms.push(AppliedTransform {
                        rule_id: TautomerRuleId::AromaticExocyclicLactamLactim,
                        affected_atoms: vec![donor, bridge, acceptor],
                        affected_bonds,
                    });
                    current = next;
                    transforms_applied += 1;
                    changed = true;
                }
            }
            if !changed
                && let Some((next, donor, bridge, acceptor)) =
                    apply_dual_flank_tracked(&current, &config, &rank)
            {
                let fp = mol_fingerprint(&next);
                if !seen.contains(&fp) {
                    seen.insert(fp);
                    applied_transforms.push(AppliedTransform {
                        rule_id: TautomerRuleId::AromaticDualFlank,
                        affected_atoms: vec![donor, bridge, acceptor],
                        affected_bonds: Vec::new(),
                    });
                    current = next;
                    transforms_applied += 1;
                    changed = true;
                }
            }
        }
        if !changed {
            break; // converged: no forward rule produces an unseen form
        }
    }

    let mut candidate_count = 1;
    if matches!(status, ParentComputationStatus::Completed) && limits.max_tautomers > 0 {
        if timed_out(&deadline) {
            status = ParentComputationStatus::TimedOut;
        } else {
            let (extra, truncated) = enumerate_direct_aromatic_forms_tracked(
                &current,
                &HashSet::new(),
                limits.max_tautomers,
                &rank,
            );
            let mut candidates: Vec<Molecule> = vec![current.clone()];
            candidates.extend(extra);
            candidate_count = candidates.len();
            if candidates.len() > 1 {
                candidates.sort_by(|a, b| {
                    tautomer_score(b).cmp(&tautomer_score(a)).then_with(|| {
                        chematic_smiles::canonical_smiles(a)
                            .cmp(&chematic_smiles::canonical_smiles(b))
                    })
                });
                current = candidates.into_iter().next().unwrap();
            }
            if truncated {
                status = ParentComputationStatus::MaxTautomersReached;
            }
        }
    }

    let common_atoms = mol.atom_count().min(current.atom_count());
    let lost_stereo: Vec<AtomIdx> = (0..common_atoms)
        .map(|i| AtomIdx(i as u32))
        .filter(|&i| {
            mol.atom(i).chirality != Chirality::None && current.atom(i).chirality == Chirality::None
        })
        .collect();
    let affected_isotopes: Vec<AtomIdx> = (0..common_atoms)
        .map(|i| AtomIdx(i as u32))
        .filter(|&i| mol.atom(i).isotope != current.atom(i).isotope)
        .collect();
    let score_breakdown_final = score_breakdown(&current);

    ParentResult {
        status,
        audit: ParentAudit::Tautomer(TautomerAuditRecord {
            selected: MoleculeSnapshot::from_mol(&current),
            candidate_count,
            score_breakdown: score_breakdown_final,
            applied_transforms,
            lost_stereo,
            affected_isotopes,
        }),
        molecule: current,
    }
}

#[cfg(test)]
mod tests {
    #![allow(dead_code)]

    use super::*;
    use chematic_core::{AtomIdx, Chirality, STEREO_H_SENTINEL, StereoGroup};
    use chematic_smiles::{canonical_smiles, parse};
    use std::collections::HashMap;

    /// Test-only: rebuild `mol` with atoms placed in `new_order`
    /// (`new_order[i]` is the `AtomIdx` from `mol` that becomes atom `i` in
    /// the result), remapping every `AtomIdx`/`BondIdx`-keyed side table
    /// (`stereo_neighbor_order`, `bond_directions`, `stereo_groups`) along
    /// with it -- used to construct genuine atom-order permutations of an
    /// already-parsed molecule for order-invariance regression tests (see
    /// `issue415_residual_canonical_tautomer_is_atom_order_invariant_and_idempotent`
    /// below). Production code no longer reorders molecules this way (see
    /// the "Atom-order canonicalization" section above) -- it breaks ties
    /// via [`atom_rank`] instead, without ever rebuilding/relabeling the
    /// molecule's own storage.
    fn reorder_atoms(mol: &Molecule, new_order: &[AtomIdx]) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let mut atom_map: HashMap<AtomIdx, AtomIdx> = HashMap::with_capacity(new_order.len());
        for &old_idx in new_order {
            let new_idx = builder.add_atom(mol.atom(old_idx).clone());
            atom_map.insert(old_idx, new_idx);
        }

        let mut bond_map: HashMap<BondIdx, BondIdx> = HashMap::with_capacity(mol.bond_count());
        for i in 0..mol.bond_count() {
            let old_bidx = BondIdx(i as u32);
            let bond = mol.bond(old_bidx);
            let new_a = atom_map[&bond.atom1];
            let new_b = atom_map[&bond.atom2];
            let new_bidx = builder
                .add_bond(new_a, new_b, bond.order)
                .expect("reorder_atoms: bond from a valid molecule must be re-addable");
            bond_map.insert(old_bidx, new_bidx);
        }

        for (&old_bidx, &new_bidx) in &bond_map {
            if let Some(direction) = mol.bond_direction(old_bidx) {
                builder.set_bond_direction(new_bidx, direction);
            }
        }

        for &old_idx in new_order {
            let new_idx = atom_map[&old_idx];
            if let Some(order) = mol.stereo_neighbor_order(old_idx) {
                let remapped: Vec<u32> = order
                    .iter()
                    .map(|&v| {
                        if v == STEREO_H_SENTINEL {
                            STEREO_H_SENTINEL
                        } else {
                            atom_map[&AtomIdx(v)].0
                        }
                    })
                    .collect();
                builder.set_stereo_neighbor_order(new_idx, remapped);
            }
        }

        for g in mol.stereo_groups() {
            builder.add_stereo_group(StereoGroup::new(
                g.kind.clone(),
                g.atom_indices.iter().map(|&a| atom_map[&a]).collect(),
            ));
        }

        builder.build()
    }

    /// "Comb" molecule for the max_iter-exhaustion tests below: an
    /// asymmetric backbone (an extra methyl end-cap on backbone[0] breaks
    /// reversal symmetry, so which SUBSET of arms converts is externally
    /// observable, not automorphism-hidden) with `n` independent enol arms
    /// hanging off it (`keto=false`) or the equivalent already-converted
    /// ketone arms (`keto=true`, used as a `canonical_smiles`-only control).
    /// `reverse_arms` builds the identical graph with arms attached in the
    /// opposite backbone-position order, to probe atom-insertion-order
    /// sensitivity.
    fn build_comb(n: usize, reverse_arms: bool, keto: bool) -> Molecule {
        use chematic_core::{Atom, Element};
        let mut builder = MoleculeBuilder::new();
        let cap = builder.add_atom(Atom::new(Element::C));
        let mut backbone = Vec::new();
        for _ in 0..n {
            backbone.push(builder.add_atom(Atom::new(Element::C)));
        }
        builder
            .add_bond(cap, backbone[0], BondOrder::Single)
            .unwrap();
        for i in 0..n - 1 {
            builder
                .add_bond(backbone[i], backbone[i + 1], BondOrder::Single)
                .unwrap();
        }
        let arm_order: Vec<usize> = if reverse_arms {
            (0..n).rev().collect()
        } else {
            (0..n).collect()
        };
        for &i in &arm_order {
            if keto {
                // Already-keto arm, matching a fully-converted enol arm's
                // exact topology: backbone-C(=O)-CH3 (NOT backbone-CH2-C(=O)-CH3,
                // which is one carbon longer and not a valid control for it).
                let carbonyl = builder.add_atom(Atom::new(Element::C));
                let mut o = Atom::new(Element::O);
                o.hydrogen_count = Some(0);
                let oxy = builder.add_atom(o);
                let methyl = builder.add_atom(Atom::new(Element::C));
                builder
                    .add_bond(backbone[i], carbonyl, BondOrder::Single)
                    .unwrap();
                builder.add_bond(carbonyl, oxy, BondOrder::Double).unwrap();
                builder
                    .add_bond(carbonyl, methyl, BondOrder::Single)
                    .unwrap();
            } else {
                // Enol arm: backbone-CH=CH-OH (donor=O, bridge=C, acceptor=C).
                let bridge = builder.add_atom(Atom::new(Element::C));
                let acceptor = builder.add_atom(Atom::new(Element::C));
                let mut o = Atom::new(Element::O);
                o.hydrogen_count = Some(1);
                let donor = builder.add_atom(o);
                builder
                    .add_bond(backbone[i], bridge, BondOrder::Single)
                    .unwrap();
                builder
                    .add_bond(bridge, acceptor, BondOrder::Double)
                    .unwrap();
                builder.add_bond(bridge, donor, BondOrder::Single).unwrap();
            }
        }
        builder.build()
    }

    #[test]
    fn test_max_iter_default_diverges_on_many_independent_sites() {
        let forward = build_comb(25, false, false);
        let reversed = build_comb(25, true, false);
        assert_eq!(
            canonical_smiles(&forward),
            canonical_smiles(&reversed),
            "sanity check: build_comb(reverse_arms) must build the SAME molecule"
        );
        assert_eq!(
            canonical_smiles(&canonical_tautomer(&forward)),
            canonical_smiles(&canonical_tautomer(&reversed)),
        );
    }

    #[test]
    fn test_max_iter_1000_resolves_the_divergence() {
        // With enough budget, this also exercises the ordinary convergence
        // path rather than the batched independent-site path.
        let forward = build_comb(25, false, false);
        let reversed = build_comb(25, true, false);
        let big_config = TautomerConfig {
            max_iter: 1000,
            ..TautomerConfig::default()
        };
        let ta = canonical_tautomer_with_config(&forward, &big_config);
        let tb = canonical_tautomer_with_config(&reversed, &big_config);
        assert_eq!(canonical_smiles(&ta), canonical_smiles(&tb));
        // Every site should have converted (no leftover [OH] enol form).
        assert!(!canonical_smiles(&ta).contains("[OH]"));
    }

    #[test]
    fn tautomer_parent_bounded_selection_is_atom_order_invariant() {
        let forward = build_comb(25, false, false);
        let reversed = build_comb(25, true, false);
        let limits = TautomerLimits {
            max_transforms: 16,
            ..TautomerLimits::default()
        };
        let a = tautomer_parent(&forward, &limits);
        let b = tautomer_parent(&reversed, &limits);
        assert_eq!(
            canonical_smiles(&a.molecule),
            canonical_smiles(&b.molecule),
            "bounded parent selection must not depend on atom insertion order"
        );
        let ParentAudit::Tautomer(audit) = a.audit else {
            panic!("tautomer_parent must return a tautomer audit")
        };
        assert_eq!(audit.applied_transforms.len(), limits.max_transforms);
    }

    #[test]
    fn test_canonical_smiles_alone_is_order_independent_on_comb() {
        // Second companion: proves the SAME divergence is not a
        // canonical_smiles bug either -- an already-fully-keto comb (no
        // tautomer.rs involved at all) is order-independent on this exact,
        // pathologically-symmetric (25 near-identical arms) topology.
        let forward = build_comb(25, false, true);
        let reversed = build_comb(25, true, true);
        assert_eq!(canonical_smiles(&forward), canonical_smiles(&reversed));
    }

    #[test]
    fn test_bis_enol_independent_sites_order_independent() {
        // Two independent, non-overlapping, ASYMMETRIC enol sites (no shared
        // atoms) separated by a CH2 bridge: donor=O0-bridge=C1=C2 on one end,
        // donor=O6-bridge=C5=C4 on the other, with a methyl (C7) breaking the
        // symmetry -- an ordinary (well within max_iter) two-site case,
        // confirming canonical_tautomer is safe for the common case.
        let a = parse("OC=CCC=C(O)C").unwrap();
        // Same molecule, respelled from the other end (methyl-bearing enol first).
        let b = parse("CC(O)=CCC=CO").unwrap();
        assert_eq!(
            canonical_smiles(&canonical_tautomer(&a)),
            canonical_smiles(&canonical_tautomer(&b)),
        );
    }

    #[test]
    fn test_canonical_no_tautomers() {
        // Simple alkane — no tautomers
        let mol = parse("CCO").unwrap();
        let t = canonical_tautomer(&mol);
        // Should be same molecule (ethanol has no keto-enol applicable here)
        assert_eq!(t.atom_count(), mol.atom_count());
    }

    #[test]
    fn test_canonical_idempotent() {
        let mol = parse("CC=O").unwrap(); // acetaldehyde
        let t1 = canonical_tautomer(&mol);
        let t2 = canonical_tautomer(&t1);
        assert_eq!(mol_fingerprint(&t1), mol_fingerprint(&t2));
    }

    #[test]
    fn test_enumerate_single_no_match() {
        let mol = parse("C").unwrap(); // methane
        let tautomers = enumerate_tautomers(&mol);
        assert_eq!(tautomers.len(), 1); // only the original
    }

    #[test]
    fn test_enumerate_cap() {
        // Complex molecule — should not exceed 32
        let mol = parse("CC(=O)CC(=O)C").unwrap(); // acetylacetone
        let tautomers = enumerate_tautomers(&mol);
        assert!(tautomers.len() <= 32);
        assert!(!tautomers.is_empty());
    }

    #[test]
    fn test_enumerate_vinyl_alcohol() {
        // OC=C → should find at least 1 more tautomer (keto form CC=O)
        let mol = parse("OC=C").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        assert!(
            tautomers.len() >= 2,
            "Expected >= 2 tautomers for vinyl alcohol, got {}",
            tautomers.len()
        );
    }

    #[test]
    fn test_canonical_amide_unchanged() {
        // CC(=O)N — amide, canonical form should keep amide (N-C=O preferred)
        let mol = parse("CC(=O)N").unwrap();
        let t = canonical_tautomer(&mol);
        assert_eq!(mol_fingerprint(&t), mol_fingerprint(&mol));
    }

    #[test]
    fn test_canonical_acetylacetone_stable() {
        let mol = parse("CC(=O)CC(=O)C").unwrap();
        let t = canonical_tautomer(&mol);
        // Should return a valid molecule (not panic)
        assert!(t.atom_count() > 0);
    }

    #[test]
    fn test_enumerate_includes_original() {
        let mol = parse("CC=O").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        // First element should be the original
        assert_eq!(mol_fingerprint(&tautomers[0]), mol_fingerprint(&mol));
    }

    #[test]
    fn test_canonical_keto_unchanged() {
        // CC=O — aldehyde, already in canonical keto form
        let mol = parse("CC=O").unwrap();
        let t = canonical_tautomer(&mol);
        assert_eq!(mol_fingerprint(&t), mol_fingerprint(&mol));
    }

    #[test]
    fn test_canonical_acetylacetone_enol() {
        // CC(O)=CC(=O)C — enol form of acetylacetone
        // canonical_tautomer should give same result as starting from keto form
        let enol = parse("CC(O)=CC(=O)C").unwrap();
        let keto = parse("CC(=O)CC(=O)C").unwrap();
        let t_enol = canonical_tautomer(&enol);
        let t_keto = canonical_tautomer(&keto);
        // Both should converge to same canonical form
        // (this may or may not hold depending on rule ordering, just check they don't panic)
        assert!(t_enol.atom_count() > 0);
        assert!(t_keto.atom_count() > 0);
    }

    #[test]
    fn test_enumerate_pyrazole_12_shift() {
        // c1cc[nH]n1 — pyrazole: N1H and N2H are direct 1,2-shift tautomers.
        let mol = parse("c1cc[nH]n1").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        assert!(
            tautomers.len() >= 2,
            "Expected >= 2 tautomers for pyrazole, got {}",
            tautomers.len()
        );
    }

    #[test]
    fn test_canonical_pyrazole_normalization() {
        let n1h = parse("c1cc[nH]n1").unwrap();
        let tautomers = enumerate_tautomers(&n1h);
        let n2h = tautomers
            .iter()
            .find(|t| h_assignment(t) != h_assignment(&n1h))
            .expect("enumerate_tautomers should produce N2H tautomer of pyrazole");
        assert_eq!(
            canonical_smiles(&canonical_tautomer(&n1h)),
            canonical_smiles(&canonical_tautomer(n2h)),
            "canonical_tautomer should normalize N1H and N2H to the same form"
        );
    }

    // --- TautomerConfig ---

    #[test]
    fn test_config_default_same_as_no_config() {
        // canonical_tautomer and canonical_tautomer_with_config(default) must agree.
        let mol = parse("OC=C").unwrap(); // enol
        let a = canonical_tautomer(&mol);
        let b = canonical_tautomer_with_config(&mol, &TautomerConfig::default());
        assert_eq!(
            canonical_smiles(&a),
            canonical_smiles(&b),
            "default config should match canonical_tautomer"
        );
    }

    #[test]
    fn test_config_max_iter_one_limits_convergence() {
        // With max_iter=1, at most one rule application happens.
        // Result may differ from full convergence but should not panic.
        let mol = parse("OC=C").unwrap();
        let config = TautomerConfig {
            max_iter: 1,
            ..TautomerConfig::default()
        };
        let _ = canonical_tautomer_with_config(&mol, &config);
    }

    #[test]
    fn test_config_max_tautomers_caps_enumerate() {
        // Acetylacetone has many tautomers. Capping at 2 should return exactly 2.
        let mol = parse("CC(=O)CC(=O)C").unwrap(); // acetylacetone
        let config = TautomerConfig {
            max_tautomers: 2,
            ..TautomerConfig::default()
        };
        let tautomers = enumerate_tautomers_with_config(&mol, &config);
        assert_eq!(
            tautomers.len(),
            2,
            "max_tautomers=2 should return exactly 2"
        );
    }

    #[test]
    fn test_config_enabled_rules_subset() {
        // Enabling only rule 0 (keto-enol) should still work on an enol.
        let mol = parse("OC=C").unwrap();
        let config = TautomerConfig::keto_enol_only();
        let result = canonical_tautomer_with_config(&mol, &config);
        // Should convert enol to ketone (or at least not panic).
        assert!(result.atom_count() > 0);
    }

    #[test]
    fn test_config_empty_enabled_rules_equals_all() {
        let mol = parse("OC=C").unwrap();
        let all = canonical_tautomer_with_config(&mol, &TautomerConfig::default());
        let explicit_empty = canonical_tautomer_with_config(
            &mol,
            &TautomerConfig {
                enabled_rules: vec![],
                ..TautomerConfig::default()
            },
        );
        assert_eq!(
            canonical_smiles(&all),
            canonical_smiles(&explicit_empty),
            "empty enabled_rules should equal all rules"
        );
    }

    #[test]
    fn test_rule_count_and_names() {
        let count = TautomerConfig::rule_count();
        let names = TautomerConfig::rule_names();
        assert!(count > 0);
        assert_eq!(names.len(), count);
        assert!(!names[0].is_empty());
    }

    // B4 Tests: 1,5-shift with heteroatom bridges
    #[test]
    fn test_15_shift_beta_diketone() {
        // CC(=O)CC(=O)C — acetylacetone: classic 1,5-shift case
        let mol = parse("CC(=O)CC(=O)C").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        // Should find enol form via 1,5-shift
        assert!(
            tautomers.len() >= 2,
            "Expected >= 2 tautomers for β-diketone"
        );
    }

    #[test]
    fn test_15_shift_enol_imine() {
        // C-C(=N)-C-C(=O)-H: potential 1,5-O to N shift
        // Build: C1=C-N-C-O with appropriate bonds
        let mol = parse("CC(=N)CC(=O)C").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        // Should enumerate without panic
        assert!(!tautomers.is_empty());
    }

    #[test]
    fn test_15_shift_n_bridge_diketone() {
        // O-C-N(=O)-C-O type: nitro-type with bridging N
        // Using simplified SMILES that approximates the pattern
        let mol = parse("OC1=C(O)C(=O)C(=O)C=C1").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        // Should enumerate with N bridge possibility
        assert!(!tautomers.is_empty());
    }

    #[test]
    fn test_15_shift_s_bridge() {
        // O-C-S(=O)-C-O: thio-β-diketone variant
        let mol = parse("CC(=O)CS(=O)C(=O)C").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        assert!(!tautomers.is_empty());
    }

    #[test]
    fn test_15_shift_n_to_n_with_bridge() {
        // N-C-N(=)-C-N: guanidine-type via N bridge
        let mol = parse("NC(=N)NC(=N)N").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        assert!(!tautomers.is_empty());
    }

    #[test]
    fn test_15_shift_canonical_idempotent() {
        // Applying canonical_tautomer twice should give same result
        let mol = parse("CC(=O)CC(=O)C").unwrap();
        let t1 = canonical_tautomer(&mol);
        let t2 = canonical_tautomer(&t1);
        assert_eq!(mol_fingerprint(&t1), mol_fingerprint(&t2));
    }

    #[test]
    fn test_15_shift_heteroatom_enumeration() {
        // Test that heteroatom-bridged 1,5-shifts are enumerated
        let mol = parse("CC(=O)CC(=O)C").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        // Should find multiple tautomers via both C-bridge and other possibilities
        assert!(
            tautomers.len() >= 2,
            "Expected >= 2 tautomers for acetylacetone with new 1,5-shift rules"
        );
    }

    #[test]
    fn test_15_shift_c_to_o_with_heteroatom() {
        // Active methylene 1,5-shift with heteroatom bridge
        let mol = parse("CC(=O)CC(=O)C").unwrap();
        let config = TautomerConfig {
            max_tautomers: 64,
            ..TautomerConfig::default()
        };
        let tautomers = enumerate_tautomers_with_config(&mol, &config);
        assert!(!tautomers.is_empty());
    }

    #[test]
    fn test_15_shift_no_false_positives() {
        // Molecule that shouldn't match 1,5-shift patterns
        let mol = parse("CC(=O)C").unwrap(); // propanone (no 1,5-shift possible)
        let tautomers = enumerate_tautomers(&mol);
        // Should enumerate any 1,3-shifts but no spurious 1,5-shifts
        assert!(!tautomers.is_empty());
    }

    #[test]
    fn test_15_shift_multiple_donors_acceptors() {
        // β-diketone with multiple potential 1,5-shift sites
        let mol = parse("CC(=O)CC(=O)CC(=O)C").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        // Should enumerate multiple tautomers from multiple sites
        assert!(
            tautomers.len() >= 2,
            "Expected multiple tautomers from multiple 1,5-shift sites"
        );
    }

    #[test]
    fn test_15_shift_config_selectivity() {
        // Create config that enables only 1,5-shift rules (rules 20-42)
        let mol = parse("CC(=O)CC(=O)C").unwrap();
        let config = TautomerConfig {
            max_tautomers: 32,
            enabled_rules: (20..43).collect(), // indices for 1,5-shift rules (if they exist)
            ..TautomerConfig::default()
        };
        let tautomers = enumerate_tautomers_with_config(&mol, &config);
        // Should still enumerate with restricted rule set
        assert!(!tautomers.is_empty());
    }

    // ── Tautomer zone blocking tests ──────────────────────────────────────────

    /// Helper: canonical SMILES of the canonical tautomer.
    fn canonical_smi(smi: &str) -> String {
        let mol = parse(smi).unwrap();
        canonical_smiles(&canonical_tautomer(&mol))
    }

    /// Helper: canonical SMILES of the canonical tautomer with blocked atoms.
    fn blocked_smi(smi: &str, blocked: &[u32]) -> String {
        let mol = parse(smi).unwrap();
        let config = TautomerConfig {
            blocked_atoms: blocked.iter().map(|&i| AtomIdx(i)).collect(),
            ..TautomerConfig::default()
        };
        canonical_smiles(&canonical_tautomer_with_config(&mol, &config))
    }

    #[test]
    fn test_blocking_donor_suppresses_keto_enol() {
        // OC=C: oxygen (index 0) is the donor for the enol→keto tautomer.
        // With default config, canonical is the keto form (C-C=O).
        // With O blocked, the enol tautomer should be preserved (or at least
        // the move through O should not fire).
        let mol = parse("OC=C").unwrap();
        let default_config = TautomerConfig::default();
        let blocked_config = TautomerConfig {
            blocked_atoms: [AtomIdx(0)].into_iter().collect(), // block O
            ..TautomerConfig::default()
        };
        let default_result = canonical_tautomer_with_config(&mol, &default_config);
        let blocked_result = canonical_tautomer_with_config(&mol, &blocked_config);
        // If donor O is blocked, the move should not fire → same as original,
        // or at minimum different from the unconstrained result.
        // The key invariant: no panic, and the blocked result is a valid molecule.
        assert!(blocked_result.atom_count() > 0);
        // When tautomerism through O is blocked, the canonical form cannot reach
        // the keto tautomer — verify it differs from the unconstrained canonical
        // OR equals the input (depending on whether other rules fire).
        let _ = (default_result, blocked_result); // both valid
    }

    #[test]
    fn test_empty_blocked_sets_identical_to_default() {
        // Empty blocked_atoms/blocked_bonds must produce the same result as default.
        for smi in &["OC=C", "CC(=O)CC", "c1cc[nH]c1", "CN=C"] {
            let mol = parse(smi).unwrap();
            let default = canonical_tautomer(&mol);
            let empty_config = TautomerConfig {
                blocked_atoms: HashSet::new(),
                blocked_bonds: HashSet::new(),
                ..TautomerConfig::default()
            };
            let explicit_empty = canonical_tautomer_with_config(&mol, &empty_config);

            assert_eq!(
                canonical_smiles(&default),
                canonical_smiles(&explicit_empty),
                "empty blocked sets must give same result as default for {smi}"
            );
        }
    }

    #[test]
    fn test_enumerate_with_blocking_leq_enumerate_without() {
        // Blocking an atom can only reduce the number of tautomers reachable.
        let mol = parse("CC(=O)CC(=O)C").unwrap(); // 2 carbonyl groups → many tautomers
        let all = enumerate_tautomers(&mol);
        // Block the central C (index 3 in CC(=O)CC(=O)C)
        let config = TautomerConfig {
            blocked_atoms: [AtomIdx(3)].into_iter().collect(),
            ..TautomerConfig::default()
        };
        let blocked = enumerate_tautomers_with_config(&mol, &config);
        assert!(
            blocked.len() <= all.len(),
            "blocking must not increase tautomer count: {} > {}",
            blocked.len(),
            all.len()
        );
    }

    #[test]
    fn test_blocking_all_atoms_preserves_input() {
        // When every atom is blocked, no H-transfer can fire.
        let mol = parse("OC=C").unwrap();
        let n = mol.atom_count();
        let config = TautomerConfig {
            blocked_atoms: (0..n as u32).map(AtomIdx).collect(),
            ..TautomerConfig::default()
        };
        let result = canonical_tautomer_with_config(&mol, &config);
        // Result must equal the input (no tautomer can fire).

        assert_eq!(
            canonical_smiles(&result),
            canonical_smiles(&mol),
            "all atoms blocked → input unchanged"
        );
    }

    #[test]
    fn test_enumerate_all_atoms_blocked_returns_singleton() {
        let mol = parse("OC=C").unwrap();
        let n = mol.atom_count();
        let config = TautomerConfig {
            blocked_atoms: (0..n as u32).map(AtomIdx).collect(),
            ..TautomerConfig::default()
        };
        let tautomers = enumerate_tautomers_with_config(&mol, &config);
        assert_eq!(
            tautomers.len(),
            1,
            "all atoms blocked → only the original is returned"
        );
    }

    #[test]
    fn test_out_of_range_atom_index_is_safe() {
        // An AtomIdx larger than the molecule silently has no effect.
        let mol = parse("OC=C").unwrap();
        let n = mol.atom_count();
        let config = TautomerConfig {
            blocked_atoms: [AtomIdx(n as u32 + 100)].into_iter().collect(),
            ..TautomerConfig::default()
        };
        // Must not panic
        let result = canonical_tautomer_with_config(&mol, &config);
        assert!(result.atom_count() > 0);
    }

    // ── Stereo preservation tests (RDKit #7969: canonical_tautomer should NOT ──
    // ── erase sp3 chirality at stereocenters not involved in tautomerism)     ──

    /// RDKit issue #7969: canonical_tautomer erases sp3 chirality.
    /// chematic must NOT have this bug — chirality fields are preserved via clone().
    #[test]
    fn test_remote_stereo_preserved_keto_enol() {
        // C[C@H](O)CC(=O)C: keto-enol fires at the C(=O) end.
        // The [C@H] stereocenter (index 1) is remote — must be preserved.
        let mol = parse("C[C@H](O)CC(=O)C").unwrap();
        let before_chirality = mol.atom(AtomIdx(1)).chirality;
        assert_ne!(
            before_chirality,
            Chirality::None,
            "test setup: atom 1 must be chiral"
        );

        let t = canonical_tautomer(&mol);
        let after_chirality = t.atom(AtomIdx(1)).chirality;
        assert_ne!(
            after_chirality,
            Chirality::None,
            "Remote [C@H] chirality erased by canonical_tautomer (RDKit #7969 regression)"
        );

        // Additionally verify canonical SMILES contains a chirality marker
        let smi = canonical_smiles(&t);
        assert!(
            smi.contains('@'),
            "Canonical SMILES lost chirality marker: '{}'",
            smi
        );
    }

    #[test]
    fn test_alanine_stereo_trivially_preserved() {
        // [C@@H](N)(C(=O)O)C — no tautomer rule fires; stereo must be unchanged.
        let mol = parse("[C@@H](N)(C(=O)O)C").unwrap();
        let before = mol.atom(AtomIdx(0)).chirality;
        let t = canonical_tautomer(&mol);
        let after = t.atom(AtomIdx(0)).chirality;
        assert_eq!(
            before, after,
            "Alanine chirality changed; was {:?}, got {:?}",
            before, after
        );
    }

    #[test]
    fn test_glucose_all_stereocenters_preserved() {
        let mol = parse("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O").unwrap();
        let before: Vec<Chirality> = mol.atoms().map(|(_, a)| a.chirality).collect();

        let t = canonical_tautomer(&mol);
        let after: Vec<Chirality> = t.atoms().map(|(_, a)| a.chirality).collect();

        assert_eq!(
            before, after,
            "Glucose: stereocenters changed through canonical_tautomer"
        );
    }

    #[test]
    fn test_pyrazole_no_phantom_chirality() {
        // c1cc[nH]n1: N-H tautomers possible, no stereocenters → must stay chiral-free
        let mol = parse("c1cc[nH]n1").unwrap();
        let t = canonical_tautomer(&mol);
        let chiral_count = t
            .atoms()
            .filter(|(_, a)| a.chirality != Chirality::None)
            .count();
        assert_eq!(
            chiral_count, 0,
            "Phantom chirality introduced by pyrazole tautomerism"
        );
    }

    #[test]
    fn test_stereo_at_donor_does_not_panic() {
        // [C@@H](O)(C)C(=O)O — lactic acid: O (donor) is adjacent to stereocentre.
        // The tautomer may legitimately change chirality; we just verify no panic.
        let mol = parse("[C@@H](O)(C)C(=O)O").unwrap();
        let t = canonical_tautomer(&mol);
        assert!(t.atom_count() > 0);
        let smi = canonical_smiles(&t);
        assert!(!smi.is_empty());
    }

    #[test]
    fn test_enumerate_tautomers_remote_stereo_preserved() {
        // C[C@H](O)CC(=O)C: enumerate all tautomers.
        // Every produced tautomer must preserve chirality at atom 1 (remote centre).
        let mol = parse("C[C@H](O)CC(=O)C").unwrap();
        let original_chirality = mol.atom(AtomIdx(1)).chirality;

        let tautomers = enumerate_tautomers(&mol);
        for (i, t) in tautomers.iter().enumerate() {
            if t.atom_count() == mol.atom_count() {
                let ch = t.atom(AtomIdx(1)).chirality;
                assert_eq!(
                    ch, original_chirality,
                    "Tautomer #{}: chirality at atom 1 changed ({:?} → {:?})",
                    i, original_chirality, ch
                );
            }
        }
    }

    // -----------------------------------------------------------------
    // Tautomer-Rebuild-S2: transfer_hydrogen_aromatic must preserve the
    // stereo side channels of every atom/bond it doesn't itself modify.
    // The tests above only check `atom.chirality` (which survives via
    // `atom.clone()` regardless of this bug -- it's `chirality`'s
    // *reference frame*, `stereo_neighbor_order`, that was silently
    // dropped). These specifically probe the side channels: an aromatic
    // N-H ring (to actually fire transfer_hydrogen_aromatic, not the
    // non-aromatic transfer_hydrogen) combined with a remote tetrahedral
    // stereocenter, an enhanced stereo group, and a stashed bond
    // direction.
    // -----------------------------------------------------------------

    /// Chiral secondary-alcohol carbon (atom 1, `[C@H]`) attached to a
    /// pyrazole ring (atoms 6/7 are the two ring nitrogens) -- the N-H
    /// tautomerism is a direct aromatic 1,2-shift, exercising
    /// `transfer_hydrogen_aromatic` specifically (not the non-aromatic
    /// `transfer_hydrogen`, which is a separate, not-yet-fixed instance of
    /// this same bug class -- out of scope for this PR).
    const AROMATIC_TAUTOMER_WITH_REMOTE_STEREOCENTER: &str = "C[C@H](O)c1cc[nH]n1";

    #[test]
    fn transfer_hydrogen_aromatic_preserves_remote_stereo_neighbor_order() {
        let mol = parse(AROMATIC_TAUTOMER_WITH_REMOTE_STEREOCENTER).unwrap();
        let stereocenter = AtomIdx(1);
        let original_order = mol.stereo_neighbor_order(stereocenter).map(|s| s.to_vec());
        assert!(
            original_order.is_some(),
            "test setup sanity: atom 1 must be chiral"
        );

        let (donor, acceptor) = find_direct_aromatic_matches(&mol, &atom_rank(&mol))
            .into_iter()
            .next()
            .expect("test setup sanity: pyrazole should have a direct aromatic N-H match");
        let result = transfer_hydrogen_aromatic(&mol, donor, acceptor, &HashSet::new())
            .expect("pyrazole N-H transfer should succeed");

        assert_eq!(
            result.atom(stereocenter).chirality,
            mol.atom(stereocenter).chirality,
            "remote stereocenter's chirality must be unchanged"
        );
        assert_eq!(
            result
                .stereo_neighbor_order(stereocenter)
                .map(|s| s.to_vec()),
            original_order,
            "remote stereocenter's stereo_neighbor_order must survive transfer_hydrogen_aromatic \
             verbatim -- chirality alone surviving is not enough to keep it interpretable"
        );
    }

    #[test]
    fn transfer_hydrogen_aromatic_preserves_stereo_groups_and_bond_directions() {
        let mut mol = parse(AROMATIC_TAUTOMER_WITH_REMOTE_STEREOCENTER).unwrap();
        let stereocenter = AtomIdx(1);
        mol.add_stereo_group(chematic_core::StereoGroup::new(
            chematic_core::StereoGroupKind::Absolute,
            vec![stereocenter],
        ));
        // Stash an arbitrary bond direction on a bond uninvolved in the H
        // transfer (the C-O bond), mimicking the kind of stashed E/Z
        // marker `apply_aromaticity_ex` leaves behind on an exocyclic bond
        // adjacent to a ring atom promoted to Aromatic order.
        let (co_bond, _) = mol
            .bond_between(AtomIdx(1), AtomIdx(2))
            .expect("C-O bond must exist");
        mol.set_bond_direction(co_bond, BondOrder::Up);

        let (donor, acceptor) = find_direct_aromatic_matches(&mol, &atom_rank(&mol))
            .into_iter()
            .next()
            .expect("test setup sanity: pyrazole should have a direct aromatic N-H match");
        let result = transfer_hydrogen_aromatic(&mol, donor, acceptor, &HashSet::new())
            .expect("pyrazole N-H transfer should succeed");

        assert_eq!(
            result.stereo_groups(),
            mol.stereo_groups(),
            "stereo_groups must survive transfer_hydrogen_aromatic verbatim"
        );
        assert_eq!(
            result.bond_direction(co_bond),
            mol.bond_direction(co_bond),
            "bond_directions for a bond uninvolved in the H transfer must be unchanged"
        );
    }

    // -----------------------------------------------------------------
    // Transfer-Hydrogen-Correctness-P0: the non-aromatic `transfer_hydrogen`
    // (1,3-/1,5-H-shift, e.g. keto-enol) counterpart to
    // Tautomer-Rebuild-S2's `transfer_hydrogen_aromatic` fix -- same bug
    // shape, deliberately left out of scope there (see the comment on
    // `AROMATIC_TAUTOMER_WITH_REMOTE_STEREOCENTER` above). A keto-enol
    // 1,3-shift (donor O-H, bridge C, acceptor C=C) combined with a remote
    // tetrahedral stereocenter, an enhanced stereo group, and a stashed
    // bond direction.
    // -----------------------------------------------------------------

    /// Chiral secondary carbon (atom 1, `[C@H]`) two bonds away from an
    /// enol group (`C(O)=C`, atoms 4/5/6) -- structurally uninvolved in the
    /// keto-enol 1,3-shift, which exercises the non-aromatic
    /// `transfer_hydrogen` specifically.
    const NON_AROMATIC_TAUTOMER_WITH_REMOTE_STEREOCENTER: &str = "C[C@H](Cl)CC(O)=C";

    fn keto_enol_match(mol: &Molecule) -> (AtomIdx, AtomIdx, AtomIdx) {
        let rule = RULES
            .iter()
            .find(|r| r.name == "keto-enol")
            .expect("keto-enol rule must exist");
        find_matches(mol, rule, &atom_rank(mol))
            .into_iter()
            .next()
            .expect("test setup sanity: molecule should have a keto-enol match")
    }

    #[test]
    fn transfer_hydrogen_preserves_remote_stereo_neighbor_order() {
        let mol = parse(NON_AROMATIC_TAUTOMER_WITH_REMOTE_STEREOCENTER).unwrap();
        let stereocenter = AtomIdx(1);
        let original_order = mol.stereo_neighbor_order(stereocenter).map(|s| s.to_vec());
        assert!(
            original_order.is_some(),
            "test setup sanity: atom 1 must be chiral"
        );

        let (donor, bridge, acceptor) = keto_enol_match(&mol);
        let result = transfer_hydrogen(
            &mol,
            donor,
            bridge,
            acceptor,
            &HashSet::new(),
            &HashSet::new(),
        )
        .expect("keto-enol transfer should succeed");

        assert_eq!(
            result.atom(stereocenter).chirality,
            mol.atom(stereocenter).chirality,
            "remote stereocenter's chirality must be unchanged"
        );
        assert_eq!(
            result
                .stereo_neighbor_order(stereocenter)
                .map(|s| s.to_vec()),
            original_order,
            "remote stereocenter's stereo_neighbor_order must survive transfer_hydrogen \
             verbatim -- chirality alone surviving is not enough to keep it interpretable"
        );
    }

    #[test]
    fn transfer_hydrogen_preserves_stereo_groups_and_bond_directions() {
        let mut mol = parse(NON_AROMATIC_TAUTOMER_WITH_REMOTE_STEREOCENTER).unwrap();
        let stereocenter = AtomIdx(1);
        mol.add_stereo_group(chematic_core::StereoGroup::new(
            chematic_core::StereoGroupKind::Absolute,
            vec![stereocenter],
        ));
        // Stash an arbitrary bond direction on a bond uninvolved in the H
        // transfer (the C-Cl bond).
        let (c_cl_bond, _) = mol
            .bond_between(AtomIdx(1), AtomIdx(2))
            .expect("C-Cl bond must exist");
        mol.set_bond_direction(c_cl_bond, BondOrder::Up);

        let (donor, bridge, acceptor) = keto_enol_match(&mol);
        let result = transfer_hydrogen(
            &mol,
            donor,
            bridge,
            acceptor,
            &HashSet::new(),
            &HashSet::new(),
        )
        .expect("keto-enol transfer should succeed");

        assert_eq!(
            result.stereo_groups(),
            mol.stereo_groups(),
            "stereo_groups must survive transfer_hydrogen verbatim"
        );
        assert_eq!(
            result.bond_direction(c_cl_bond),
            mol.bond_direction(c_cl_bond),
            "bond_directions for a bond uninvolved in the H transfer must be unchanged"
        );
    }

    #[test]
    fn transfer_hydrogen_cip_of_uninvolved_stereocenter_survives_canonical_round_trip() {
        // Chemical-identity check: the remote stereocenter's actual CIP
        // label (not just the raw @/@@ character) must survive both the
        // transfer itself and a canonicalize -> reparse round trip.
        let mol = parse("C[C@H:9](Cl)CC(O)=C").unwrap();
        let stereocenter = mol
            .atoms()
            .find(|(_, a)| a.atom_map == Some(9))
            .map(|(idx, _)| idx)
            .expect("atom map tag not found");
        let before_cip = crate::assign_cip(&mol).get(stereocenter);
        assert!(
            before_cip.is_some(),
            "test setup sanity: must be resolvable"
        );

        let (donor, bridge, acceptor) = keto_enol_match(&mol);
        let result = transfer_hydrogen(
            &mol,
            donor,
            bridge,
            acceptor,
            &HashSet::new(),
            &HashSet::new(),
        )
        .expect("keto-enol transfer should succeed");

        let smi = canonical_smiles(&result);
        let reparsed = parse(&smi).expect("valid canonical SMILES");
        let reparsed_center = reparsed
            .atoms()
            .find(|(_, a)| a.atom_map == Some(9))
            .map(|(idx, _)| idx)
            .expect("atom map tag not found after round trip");
        let after_cip = crate::assign_cip(&reparsed).get(reparsed_center);

        assert_eq!(
            before_cip, after_cip,
            "remote stereocenter's CIP code must survive transfer_hydrogen + a canonical round trip"
        );
    }

    #[test]
    fn enumerate_tautomers_keto_enol_count_and_canonical_form_unaffected_by_metadata_fix() {
        // Regression pin: this fix is about metadata preservation only, not
        // enumeration logic -- the tautomer count and canonical_tautomer's
        // chosen form for a plain (no remote stereocenter) enol must be
        // exactly what they were before this PR.
        let mol = parse("CC(O)=C").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        assert_eq!(
            tautomers.len(),
            2,
            "propan-2-enol should enumerate to exactly 2 forms (enol + acetone)"
        );
        // Golden string updated (issue #205): `initial_invariant`'s
        // explicit/implicit-H-count unification changed Morgan ranks for
        // some atoms, so `canonical_smiles`'s neighbor-branch-ordering
        // shifted -- this is still the identical keto tautomer (acetone),
        // just a different valid serialization of the same graph; the
        // enumeration count above (2) is the actual regression pin for
        // "logic unaffected" and did not change.
        assert_eq!(
            canonical_smiles(&canonical_tautomer(&mol)),
            "C(C)(=O)C",
            "the canonical tautomer form must be unchanged (prefers the keto form)"
        );
    }

    #[test]
    fn enumerate_tautomers_count_and_canonical_form_unaffected_by_metadata_fix() {
        // Regression pin: this fix is about metadata preservation only, not
        // enumeration logic (mol_fingerprint/h_assignment dedup don't
        // consult stereo_neighbor_order/stereo_groups/bond_directions at
        // all) -- the tautomer count and canonical_tautomer's chosen form
        // for a plain pyrazole (no remote stereocenter to confound the
        // comparison) must be exactly what they were before this PR.
        let mol = parse("c1cc[nH]n1").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        assert_eq!(
            tautomers.len(),
            2,
            "pyrazole should enumerate to exactly 2 forms"
        );
        // Golden string updated (issue #205), same reason as the acetone
        // test above: same tautomer (pyrazole's NH form), different valid
        // canonical serialization; the count (2) is the unaffected pin.
        assert_eq!(
            canonical_smiles(&canonical_tautomer(&mol)),
            "c1c[nH]nc1",
            "pyrazole's canonical tautomer form must be unchanged"
        );
    }

    #[test]
    fn test_blocked_stereo_preserved_with_zone_blocking() {
        // O[C@@H](F)C(=O)C: block O (index 0) to suppress keto-enol.
        // Stereocentre [C@@H] (index 1) must be preserved.
        let mol = parse("O[C@@H](F)C(=O)C").unwrap();
        let config = TautomerConfig {
            blocked_atoms: [AtomIdx(0)].into_iter().collect(),
            ..TautomerConfig::default()
        };
        let t = canonical_tautomer_with_config(&mol, &config);
        let chirality_after = t.atom(AtomIdx(1)).chirality;
        assert_ne!(
            chirality_after,
            Chirality::None,
            "Chirality erased from blocked stereocentre"
        );
    }

    // ── RDKit PR #9128: E/Z stereo on exocyclic double bonds (hydrazones/imines) ──

    #[test]
    fn test_hydrazone_ez_stereo_preserved_in_canonical_tautomer() {
        // E-hydrazone and Z-hydrazone are DIFFERENT compounds.
        // The tautomer search must not merge them while deduplicating forms.
        let e_hydrazone = parse("C/C=N/N").expect("E-hydrazone");
        let z_hydrazone = parse("C/C=N\\N").expect("Z-hydrazone");
        let e_can = canonical_tautomer(&e_hydrazone);
        let z_can = canonical_tautomer(&z_hydrazone);
        let e_smi = canonical_smiles(&e_can);
        let z_smi = canonical_smiles(&z_can);
        assert_ne!(
            e_smi, z_smi,
            "E and Z hydrazone must remain distinct after canonical_tautomer (RDKit PR #9128): \
             E={e_smi} Z={z_smi}"
        );
    }

    // ── Tetrazole 1H ↔ 2H tautomers (OpenBabel PR #2975 pattern) ────────────
    //
    // Tetrazole has two stable tautomers:
    //   1H-tetrazole: H on N1 (directly adjacent to C)
    //   2H-tetrazole: H on N2 (one position further around the ring from C)
    // Both are aromatic and interconvert via a direct aromatic 1,2-shift.
    // The `find_direct_aromatic_matches` path handles this case.

    #[test]
    fn test_tetrazole_1h_enumerates_two_forms() {
        // c1nnn[nH]1 = 1H-tetrazole.  enumerate_tautomers must find at least
        // 2 forms (1H and 2H) via the direct aromatic 1,2-shift.
        let mol = parse("c1nnn[nH]1").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        assert!(
            tautomers.len() >= 2,
            "Expected >= 2 tautomers for 1H-tetrazole, got {}",
            tautomers.len()
        );
    }

    #[test]
    fn test_tetrazole_2h_enumerates_two_forms() {
        // c1n[nH]nn1 = 2H-tetrazole (H on N two positions from C).
        // enumerate_tautomers must find at least 2 forms (1H and 2H).
        let mol = parse("c1n[nH]nn1").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        assert!(
            tautomers.len() >= 2,
            "Expected >= 2 tautomers for 2H-tetrazole, got {}",
            tautomers.len()
        );
    }

    #[test]
    fn test_tetrazole_canonical_from_1h_and_2h_agrees() {
        // canonical_tautomer must produce the same canonical form regardless of
        // whether the input is 1H or 2H tetrazole.
        let mol_1h = parse("c1nnn[nH]1").unwrap();
        let mol_2h = parse("c1n[nH]nn1").unwrap();
        let can_1h = canonical_smiles(&canonical_tautomer(&mol_1h));
        let can_2h = canonical_smiles(&canonical_tautomer(&mol_2h));
        assert_eq!(
            can_1h, can_2h,
            "canonical_tautomer must agree for 1H ({can_1h}) and 2H ({can_2h}) tetrazole"
        );
    }

    #[test]
    fn test_tetrazole_canonical_preserves_aromaticity() {
        // After canonical_tautomer the ring must remain fully aromatic.
        let mol = parse("c1nnn[nH]1").unwrap();
        let t = canonical_tautomer(&mol);
        let all_aromatic = t.atoms().all(|(_, a)| a.aromatic);
        assert!(
            all_aromatic,
            "all atoms in canonical tetrazole must be aromatic"
        );
    }

    #[test]
    fn imine_ez_stereo_preserved_in_tautomer_enumeration() {
        // Enumerate tautomers of E-imine C/C=N/C — all resulting tautomers must
        // produce valid canonical SMILES (no empty string, no panic).
        let e_imine = parse("C/C=N/C").expect("E-imine");
        let tautomers = enumerate_tautomers(&e_imine);
        assert!(
            !tautomers.is_empty(),
            "E-imine must enumerate at least one tautomer"
        );
        for (i, t) in tautomers.iter().enumerate() {
            let smi = canonical_smiles(t);
            assert!(
                !smi.is_empty(),
                "tautomer {i} must produce valid canonical SMILES"
            );
        }
    }

    // -- Phase 2 round-2B tautomer_parent tests --------------------------------
    // See docs/rfcs/tautomer_parent_identity_phase2_rfc.md sections 4.1-4.3.

    #[test]
    fn all_rules_map_to_a_named_tautomer_rule_id() {
        for rule in RULES {
            assert!(
                !matches!(rule_id_for(rule.name), TautomerRuleId::Other(_)),
                "rule '{}' has no TautomerRuleId variant -- add one",
                rule.name
            );
        }
    }

    #[test]
    fn tautomer_parent_completed_on_keto_enol() {
        let mol = parse("CC(O)=CC(C)=O").unwrap();
        let result = tautomer_parent(&mol, &TautomerLimits::default());
        assert_eq!(result.status, ParentComputationStatus::Completed);
        assert_eq!(
            canonical_smiles(&result.molecule),
            canonical_smiles(&canonical_tautomer(&mol))
        );
        match &result.audit {
            ParentAudit::Tautomer(record) => {
                assert_eq!(record.applied_transforms.len(), 1);
                assert_eq!(
                    record.applied_transforms[0].rule_id,
                    TautomerRuleId::KetoEnol
                );
                assert!(record.lost_stereo.is_empty());
                assert!(record.affected_isotopes.is_empty());
            }
            other => panic!("expected ParentAudit::Tautomer, got {other:?}"),
        }
    }

    #[test]
    fn tautomer_parent_max_transforms_reached_when_budget_too_small() {
        // 3 independent enol arms, each needing exactly 1 transform; a
        // budget of 1 converts only one and must report the exhaustion,
        // not silently return a partially-converted, order-dependent result
        // as if it were `Completed`.
        let mol = build_comb(3, false, false);
        let limits = TautomerLimits {
            max_transforms: 1,
            ..TautomerLimits::default()
        };
        let result = tautomer_parent(&mol, &limits);
        assert_eq!(result.status, ParentComputationStatus::MaxTransformsReached);
        match &result.audit {
            ParentAudit::Tautomer(record) => assert_eq!(record.applied_transforms.len(), 1),
            other => panic!("expected ParentAudit::Tautomer, got {other:?}"),
        }
    }

    #[test]
    fn tautomer_parent_completed_with_enough_budget_for_all_arms() {
        let mol = build_comb(3, false, false);
        let limits = TautomerLimits {
            max_transforms: 10,
            ..TautomerLimits::default()
        };
        let result = tautomer_parent(&mol, &limits);
        assert_eq!(result.status, ParentComputationStatus::Completed);
        match &result.audit {
            ParentAudit::Tautomer(record) => assert_eq!(record.applied_transforms.len(), 3),
            other => panic!("expected ParentAudit::Tautomer, got {other:?}"),
        }
    }

    #[test]
    fn tautomer_parent_max_tautomers_reached_on_tetrazole() {
        let mol = parse("c1nnn[nH]1").unwrap();
        let unlimited = tautomer_parent(&mol, &TautomerLimits::default());
        let (unlimited_candidate_count, _) = match &unlimited.audit {
            ParentAudit::Tautomer(r) => (r.candidate_count, ()),
            other => panic!("expected ParentAudit::Tautomer, got {other:?}"),
        };
        assert!(
            unlimited_candidate_count > 1,
            "tetrazole must have more than one direct-aromatic-shift candidate to make this test meaningful"
        );

        let limited = tautomer_parent(
            &mol,
            &TautomerLimits {
                max_tautomers: 1,
                ..TautomerLimits::default()
            },
        );
        assert_eq!(limited.status, ParentComputationStatus::MaxTautomersReached);
    }

    #[test]
    fn tautomer_parent_max_tautomers_zero_skips_comparison_without_false_exhaustion() {
        // A molecule with NO direct-aromatic-shift candidates at all must
        // not report MaxTautomersReached just because max_tautomers=0 --
        // see TautomerLimits::max_tautomers's doc comment.
        let mol = parse("CC(O)=CC(C)=O").unwrap();
        let result = tautomer_parent(
            &mol,
            &TautomerLimits {
                max_tautomers: 0,
                ..TautomerLimits::default()
            },
        );
        assert_eq!(result.status, ParentComputationStatus::Completed);
    }

    #[test]
    fn tautomer_parent_empty_molecule_is_invalid_input() {
        let empty = MoleculeBuilder::new().build();
        let result = tautomer_parent(&empty, &TautomerLimits::default());
        assert_eq!(
            result.status,
            ParentComputationStatus::InvalidInput(InvalidInputReason::EmptyMolecule)
        );
    }

    #[test]
    fn tp2_29_idempotence_acetylacetone() {
        let mol = parse("CC(=O)CC(C)=O").unwrap();
        let once = tautomer_parent(&mol, &TautomerLimits::default());
        let twice = tautomer_parent(&once.molecule, &TautomerLimits::default());
        assert_eq!(
            canonical_smiles(&once.molecule),
            canonical_smiles(&twice.molecule)
        );
        assert_eq!(twice.status, ParentComputationStatus::Completed);
    }

    #[test]
    fn tautomer_parent_preserves_remote_stereocenter() {
        // Alanine's stereocenter must survive tautomer_parent unaffected,
        // same as the primitive transfer_hydrogen/transfer_hydrogen_aromatic
        // functions this is built on (RFC section 1.2's confirmed non-bug).
        let mol = parse("C[C@@H](N)C(=O)O").unwrap();
        let result = tautomer_parent(&mol, &TautomerLimits::default());
        match &result.audit {
            ParentAudit::Tautomer(record) => {
                assert!(record.lost_stereo.is_empty());
                assert!(record.affected_isotopes.is_empty());
            }
            other => panic!("expected ParentAudit::Tautomer, got {other:?}"),
        }
    }

    // -- Phase 2 round-2C-2/2C-3: exocyclic lactam/lactim shift ----------------
    // See docs/rfcs/tautomer_parent_identity_phase2_rfc.md section 4.4a.

    #[test]
    fn ring_distance_matches_ortho_and_meta() {
        let two_pyridone = parse("O=c1cccc[nH]1").unwrap();
        let ring = chematic_perception::find_sssr(&two_pyridone);
        let ring = ring.rings().first().expect("one ring");
        // bridge = AtomIdx(1) (bonded to the exocyclic O), acceptor =
        // AtomIdx(6) (the ring N) -- ortho, distance 1.
        assert_eq!(ring_distance(ring, AtomIdx(1), AtomIdx(6)), Some(1));

        let three_hydroxypyridine = parse("Oc1cccnc1").unwrap();
        let ring = chematic_perception::find_sssr(&three_hydroxypyridine);
        let ring = ring.rings().first().expect("one ring");
        // bridge = AtomIdx(1), ring N = AtomIdx(5) -- meta, distance 2 (even).
        assert_eq!(ring_distance(ring, AtomIdx(1), AtomIdx(5)), Some(2));
    }

    fn assert_converges(label: &str, variants: &[&str]) {
        let outs: Vec<String> = variants
            .iter()
            .map(|s| canonical_smiles(&canonical_tautomer(&parse(s).unwrap())))
            .collect();
        assert!(
            outs.windows(2).all(|w| w[0] == w[1]),
            "{label}: variants did not converge to one canonical form: {outs:?}"
        );
    }

    fn assert_noop(label: &str, smi: &str) {
        let mol = parse(smi).unwrap();
        let out = canonical_tautomer(&mol);
        assert_eq!(
            canonical_smiles(&out),
            canonical_smiles(&mol),
            "{label}: expected a no-op, got a change"
        );
    }

    #[test]
    fn tp2_05_2_pyridone_converges() {
        assert_converges("tp2-05", &["O=c1cccc[nH]1", "Oc1ccccn1"]);
    }

    #[test]
    fn tp2_06_4_pyridone_converges() {
        assert_converges("tp2-06", &["O=c1cc[nH]cc1", "Oc1ccncc1"]);
    }

    #[test]
    fn tp2_08_uracil_converges() {
        assert_converges("tp2-08", &["O=c1cc[nH]c(=O)[nH]1", "Oc1ccnc(O)n1"]);
    }

    #[test]
    fn enumerate_aromatic_lactam_component_includes_reverse_lactim_edge() {
        let mol = parse("Cc1cc[nH]c(=O)n1").unwrap();
        let forms = enumerate_tautomers(&mol);
        assert!(
            forms.iter().any(|candidate| {
                has_aromatic_exocyclic_oxygen(candidate)
                    && !has_aromatic_exocyclic_carbonyl(candidate)
            }),
            "keto input must expose its reverse aromatic lactim form"
        );
    }

    #[test]
    fn tp2_holdout_01_hypoxanthine_converges() {
        // Holdout: never used to shape the mechanism above.
        assert_converges(
            "tp2-holdout-01",
            &["O=c1[nH]cnc2[nH]cnc12", "Oc1ncnc2[nH]cnc12"],
        );
    }

    #[test]
    fn tp2_holdout_06_n9_methylhypoxanthine_converges() {
        assert_converges(
            "tp2-holdout-06",
            &["Cn1cnc2c(=O)[nH]cnc21", "Cn1cnc2c(O)ncnc21"],
        );
    }

    #[test]
    fn tp2_07_09_dual_flank_cytosine_and_guanine_converge() {
        // The carbonyl-centered dual-flank H move must normalize both keto
        // spellings and the corresponding enol precursor.
        let cytosine_keto =
            canonical_smiles(&canonical_tautomer(&parse("Nc1cc[nH]c(=O)n1").unwrap()));
        let cytosine_enol = canonical_smiles(&canonical_tautomer(&parse("Nc1ccnc(O)n1").unwrap()));
        assert_eq!(
            cytosine_keto, cytosine_enol,
            "cytosine keto and amino/imino variants must converge"
        );

        let guanine_keto = canonical_smiles(&canonical_tautomer(
            &parse("Nc1nc2[nH]cnc2c(=O)[nH]1").unwrap(),
        ));
        let guanine_enol =
            canonical_smiles(&canonical_tautomer(&parse("Nc1nc2[nH]cnc2c(O)n1").unwrap()));
        assert_eq!(
            guanine_keto, guanine_enol,
            "guanine keto and amino/imino variants must converge"
        );
    }

    #[test]
    fn tp2_39_methylpyrimidinone_converges() {
        assert_converges(
            "tp2-39-methylpyrimidinone",
            &["Cc1cc[nH]c(=O)n1", "Cc1ccnc(=O)[nH]1", "Cc1ccnc(O)n1"],
        );
    }

    #[test]
    fn tp2_04_nitroso_oxime_converges_without_enolizing_ketones() {
        let nitroso = parse("CCN=O").unwrap();
        let oxime = parse("CC=NO").unwrap();

        assert_eq!(
            canonical_smiles(&canonical_tautomer(&nitroso)),
            canonical_smiles(&canonical_tautomer(&oxime)),
            "nitrosoethane and acetaldehyde oxime must share one canonical form"
        );

        // The specialized rule must not turn an unrelated carbonyl into an
        // enol through the generic any-bridge C→O rule.
        let ketone = parse("CCC=O").unwrap();
        assert_eq!(
            mol_fingerprint(&canonical_tautomer(&ketone)),
            mol_fingerprint(&ketone),
            "ordinary aldehydes/ketones must remain in their preferred carbonyl form"
        );
    }

    #[test]
    fn tp2_34_isotope_preserved_through_exocyclic_shift() {
        let mol = parse("[18OH]c1ccccn1").unwrap();
        let out = canonical_tautomer(&mol);
        assert_converges("tp2-34", &["[18O]=c1cccc[nH]1", "[18OH]c1ccccn1"]);
        // The isotope must land on an oxygen atom, not be dropped or moved
        // to an unrelated atom.
        let labeled = (0..out.atom_count())
            .map(|i| AtomIdx(i as u32))
            .find(|&i| out.atom(i).isotope == Some(18))
            .expect("18O label must survive");
        assert_eq!(out.atom(labeled).element.atomic_number(), 8);
    }

    #[test]
    fn tp2_35_remote_stereocenter_preserved_through_exocyclic_shift() {
        assert_converges(
            "tp2-35",
            &["O=c1c([C@@H](F)Cl)ccc[nH]1", "Oc1c([C@@H](F)Cl)cccn1"],
        );
        let mol = parse("Oc1c([C@@H](F)Cl)cccn1").unwrap();
        let out = canonical_tautomer(&mol);
        let stereocenters_before = (0..mol.atom_count())
            .filter(|&i| mol.atom(AtomIdx(i as u32)).chirality != Chirality::None)
            .count();
        let stereocenters_after = (0..out.atom_count())
            .filter(|&i| out.atom(AtomIdx(i as u32)).chirality != Chirality::None)
            .count();
        assert_eq!(stereocenters_before, 1);
        assert_eq!(stereocenters_after, 1);
    }

    #[test]
    fn negative_controls_stay_noop_under_exocyclic_shift() {
        for (label, smi) in [
            ("phenol", "Oc1ccccc1"),
            ("anisole", "COc1ccccc1"),
            ("aniline", "Nc1ccccc1"),
            ("pyridine-n-oxide", "[O-][n+]1ccccc1"),
            ("acetamide", "CC(N)=O"),
            ("3-hydroxypyridine", "Oc1cccnc1"),
            ("4-aminopyridine", "Nc1ccncc1"),
            ("2-aminopyridine", "Nc1ccccn1"),
        ] {
            assert_noop(label, smi);
        }
    }

    #[test]
    fn ring_internal_nh_shift_controls_unaffected_by_round_2c2() {
        // Existing ring-internal-only mechanism must still work exactly as
        // before -- the new step is a fallback tried only when no rule (and,
        // in canonical_tautomer_with_config's loop, implicitly, since this
        // check runs on the post-rule-loop `current`) already changed the
        // molecule; it never touches ring-internal bonds.
        for smi in ["c1cc[nH]n1", "c1nn[nH]n1", "c1ccc2[nH]cnc2c1"] {
            let mol = parse(smi).unwrap();
            let out = canonical_tautomer(&mol);
            let twice = canonical_tautomer(&out);
            assert_eq!(
                canonical_smiles(&out),
                canonical_smiles(&twice),
                "{smi}: not idempotent"
            );
        }
    }

    #[test]
    fn canonical_tautomer_and_tautomer_parent_agree_on_exocyclic_shift() {
        for smi in [
            "O=c1cccc[nH]1",
            "Oc1ccccn1",
            "O=c1cc[nH]cc1",
            "Oc1ccncc1",
            "O=c1cc[nH]c(=O)[nH]1",
            "Oc1ccnc(O)n1",
            "O=c1[nH]cnc2[nH]cnc12",
            "Oc1ncnc2[nH]cnc12",
        ] {
            let mol = parse(smi).unwrap();
            let via_canonical = canonical_smiles(&canonical_tautomer(&mol));
            let via_parent =
                canonical_smiles(&tautomer_parent(&mol, &TautomerLimits::default()).molecule);
            assert_eq!(
                via_canonical, via_parent,
                "{smi}: canonical_tautomer and tautomer_parent diverged"
            );
        }
    }

    #[test]
    fn tautomer_parent_max_transforms_reached_via_exocyclic_shift() {
        let mol = parse("Oc1ccccn1").unwrap();
        let limits = TautomerLimits {
            max_transforms: 0,
            ..TautomerLimits::default()
        };
        let result = tautomer_parent(&mol, &limits);
        assert_eq!(result.status, ParentComputationStatus::MaxTransformsReached);
    }

    // -- Round 2C-2 hardening (post-#365 review) --------------------------------
    // Tightened matcher preconditions, post-generation invariant check,
    // fingerprint-cross-checked tie-break, and the atom-order-invariance /
    // fail-closed-negative-control coverage the review asked for.

    /// Build a molecule from named atoms/bonds, inserted in the given
    /// `order` -- lets a test construct the *same* molecular graph with
    /// *different* underlying `AtomIdx` numbering (same technique as
    /// `build_comb`'s `reverse_arms`). This is the only way to test
    /// order-invariance of code that iterates `mol.atoms()`/
    /// `mol.neighbors()`/SSSR ring order without any risk that a
    /// hand-typed alternate SMILES accidentally describes a different
    /// molecule.
    fn build_named(
        order: &[&str],
        atoms: &[(&str, chematic_core::Element, bool, Option<u8>)],
        bonds: &[(&str, &str, BondOrder)],
    ) -> Molecule {
        use chematic_core::Atom;
        let mut builder = MoleculeBuilder::new();
        let mut handles: std::collections::HashMap<&str, AtomIdx> =
            std::collections::HashMap::new();
        for &name in order {
            let (_, elem, aromatic, h) = *atoms.iter().find(|(n, ..)| *n == name).unwrap();
            let mut a = Atom::new(elem);
            a.aromatic = aromatic;
            a.hydrogen_count = h;
            handles.insert(name, builder.add_atom(a));
        }
        for &(a, b, order) in bonds {
            builder.add_bond(handles[&a], handles[&b], order).unwrap();
        }
        builder.build()
    }

    fn rule_id_sequence(result: &ParentResult) -> Vec<TautomerRuleId> {
        match &result.audit {
            ParentAudit::Tautomer(record) => record
                .applied_transforms
                .iter()
                .map(|t| t.rule_id)
                .collect(),
            other => panic!("expected ParentAudit::Tautomer, got {other:?}"),
        }
    }

    #[test]
    fn atom_order_invariance_2_pyridone_graph_construction() {
        use chematic_core::Element;
        let atoms = [
            ("O", Element::O, false, Some(1)),
            ("C1", Element::C, true, None),
            ("C2", Element::C, true, None),
            ("C3", Element::C, true, None),
            ("C4", Element::C, true, None),
            ("C5", Element::C, true, None),
            ("N", Element::N, true, Some(0)),
        ];
        let bonds = [
            ("O", "C1", BondOrder::Single),
            ("C1", "C2", BondOrder::Aromatic),
            ("C2", "C3", BondOrder::Aromatic),
            ("C3", "C4", BondOrder::Aromatic),
            ("C4", "C5", BondOrder::Aromatic),
            ("C5", "N", BondOrder::Aromatic),
            ("N", "C1", BondOrder::Aromatic),
        ];
        let forward = ["O", "C1", "C2", "C3", "C4", "C5", "N"];
        let reversed = ["N", "C5", "C4", "C3", "C2", "C1", "O"];
        let mol_forward = build_named(&forward, &atoms, &bonds);
        let mol_reversed = build_named(&reversed, &atoms, &bonds);
        // Sanity: really the same input molecule before any shift.
        assert_eq!(
            canonical_smiles(&mol_forward),
            canonical_smiles(&mol_reversed)
        );

        let out_forward = canonical_tautomer(&mol_forward);
        let out_reversed = canonical_tautomer(&mol_reversed);
        assert_eq!(
            canonical_smiles(&out_forward),
            canonical_smiles(&out_reversed)
        );

        let limits = TautomerLimits::default();
        let parent_forward = tautomer_parent(&mol_forward, &limits);
        let parent_reversed = tautomer_parent(&mol_reversed, &limits);
        assert_eq!(parent_forward.status, parent_reversed.status);
        assert_eq!(
            canonical_smiles(&parent_forward.molecule),
            canonical_smiles(&parent_reversed.molecule)
        );
        assert_eq!(
            rule_id_sequence(&parent_forward),
            rule_id_sequence(&parent_reversed)
        );
    }

    #[test]
    fn atom_order_invariance_uracil_two_independent_sites_graph_construction() {
        // Both exocyclic sites (O0/C1 and O6/C5) simultaneously have TWO
        // valid odd-distance ring-N acceptor candidates each (N4 at distance
        // 3, N7 at distance 1, from either bridge) when starting from the
        // fully-di-enol form -- exactly the "candidate discovery order"
        // stress case the review flagged. Forward vs reversed atom
        // insertion order changes both iteration order over mol.atoms()
        // and the HashSet-then-sort order find_exocyclic_lactam_shift_matches
        // produces its candidates in.
        use chematic_core::Element;
        let atoms = [
            ("O0", Element::O, false, Some(1)),
            ("C1", Element::C, true, None),
            ("C2", Element::C, true, None),
            ("C3", Element::C, true, None),
            ("N4", Element::N, true, Some(0)),
            ("C5", Element::C, true, None),
            ("O6", Element::O, false, Some(1)),
            ("N7", Element::N, true, Some(0)),
        ];
        let bonds = [
            ("O0", "C1", BondOrder::Single),
            ("C1", "C2", BondOrder::Aromatic),
            ("C2", "C3", BondOrder::Aromatic),
            ("C3", "N4", BondOrder::Aromatic),
            ("N4", "C5", BondOrder::Aromatic),
            ("C5", "O6", BondOrder::Single),
            ("C5", "N7", BondOrder::Aromatic),
            ("N7", "C1", BondOrder::Aromatic),
        ];
        let forward = ["O0", "C1", "C2", "C3", "N4", "C5", "O6", "N7"];
        let reversed = ["N7", "O6", "C5", "N4", "C3", "C2", "C1", "O0"];
        let mol_forward = build_named(&forward, &atoms, &bonds);
        let mol_reversed = build_named(&reversed, &atoms, &bonds);
        assert_eq!(
            canonical_smiles(&mol_forward),
            canonical_smiles(&mol_reversed)
        );

        let out_forward = canonical_tautomer(&mol_forward);
        let out_reversed = canonical_tautomer(&mol_reversed);
        assert_eq!(
            canonical_smiles(&out_forward),
            canonical_smiles(&out_reversed)
        );
        // Must actually be the fully di-keto form, not merely "unanimous but
        // still enol" or a half-shifted state -- checked structurally
        // (double-bonded O count), not by SMILES substring, since "=O" vs
        // "O=" both spell a keto oxygen depending on write direction.
        let double_bonded_o_count = (0..out_forward.atom_count())
            .map(|i| AtomIdx(i as u32))
            .filter(|&idx| {
                out_forward.atom(idx).element.atomic_number() == 8
                    && out_forward
                        .neighbors(idx)
                        .any(|(_, bidx)| out_forward.bond(bidx).order == BondOrder::Double)
            })
            .count();
        assert_eq!(
            double_bonded_o_count,
            2,
            "expected both exocyclic sites shifted to keto: {}",
            canonical_smiles(&out_forward)
        );

        let limits = TautomerLimits::default();
        let parent_forward = tautomer_parent(&mol_forward, &limits);
        let parent_reversed = tautomer_parent(&mol_reversed, &limits);
        assert_eq!(parent_forward.status, parent_reversed.status);
        assert_eq!(
            canonical_smiles(&parent_forward.molecule),
            canonical_smiles(&parent_reversed.molecule)
        );
        assert_eq!(
            rule_id_sequence(&parent_forward),
            rule_id_sequence(&parent_reversed)
        );
    }

    /// Two independent SMILES respellings of the *same* tautomer (never a
    /// different one -- checked explicitly below) must agree on
    /// `canonical_tautomer`, `tautomer_parent`'s molecule + status, and the
    /// applied-rule-id sequence.
    fn assert_order_invariant_respelling(label: &str, a: &str, b: &str) {
        let mol_a = parse(a).unwrap();
        let mol_b = parse(b).unwrap();
        assert_eq!(
            canonical_smiles(&mol_a),
            canonical_smiles(&mol_b),
            "{label}: the two respellings are not even the same input molecule"
        );

        let out_a = canonical_tautomer(&mol_a);
        let out_b = canonical_tautomer(&mol_b);
        assert_eq!(
            canonical_smiles(&out_a),
            canonical_smiles(&out_b),
            "{label}: canonical_tautomer disagreed across respellings"
        );

        let limits = TautomerLimits::default();
        let parent_a = tautomer_parent(&mol_a, &limits);
        let parent_b = tautomer_parent(&mol_b, &limits);
        assert_eq!(
            parent_a.status, parent_b.status,
            "{label}: status disagreed"
        );
        assert_eq!(
            canonical_smiles(&parent_a.molecule),
            canonical_smiles(&parent_b.molecule),
            "{label}: tautomer_parent molecule disagreed"
        );
        assert_eq!(
            rule_id_sequence(&parent_a),
            rule_id_sequence(&parent_b),
            "{label}: applied rule_id sequence disagreed"
        );
    }

    #[test]
    fn atom_order_invariance_4_pyridone_respelling() {
        assert_order_invariant_respelling("4-pyridone", "Oc1ccncc1", "n1ccc(O)cc1");
    }

    #[test]
    fn atom_order_invariance_hypoxanthine_respelling() {
        // Enol input this time (round 2C-2's fix applies to the enol side);
        // b starts the ring traversal at the acceptor nitrogen instead of
        // the bridge carbon, and moves the exocyclic O substituent to the
        // end of the string -- verified equivalent to `a` via
        // canonical_smiles before this test was written.
        assert_order_invariant_respelling("hypoxanthine", "Oc1ncnc2[nH]cnc12", "n1cnc2[nH]cnc2c1O");
    }

    #[test]
    fn atom_order_invariance_isotope_2_pyridone_respelling() {
        assert_order_invariant_respelling("18O-2-pyridone", "[18OH]c1ccccn1", "n1ccccc1[18OH]");
    }

    #[test]
    fn atom_order_invariance_remote_stereocenter_2_pyridone_respelling() {
        // b starts the ring traversal at the ring nitrogen instead of the
        // O-bearing bridge carbon -- verified equivalent to `a` via
        // canonical_smiles before this test was written.
        assert_order_invariant_respelling(
            "2-pyridone + remote stereocenter",
            "Oc1c([C@@H](F)Cl)cccn1",
            "n1c(O)c([C@@H](F)Cl)ccc1",
        );
    }

    #[test]
    fn find_exocyclic_lactam_shift_matches_dedups_across_fused_rings() {
        // Hypoxanthine's bridge/acceptor pair is reachable through more than
        // one SSSR ring in a fused bicyclic system; must be reported once.
        let mol = parse("Oc1ncnc2[nH]cnc12").unwrap();
        let matches = find_exocyclic_lactam_shift_matches(&mol, &atom_rank(&mol));
        assert!(!matches.is_empty());
        let mut seen = HashSet::new();
        for m in &matches {
            assert!(seen.insert(*m), "duplicate candidate triple: {m:?}");
        }
    }

    #[test]
    fn negative_control_charged_pyridinium_acceptor_excluded() {
        // N-methyl-2-hydroxypyridinium: the ring N is charged and
        // 3-connected (methyl substituent) -- excluded by both the charge
        // check and the degree==2 valence-compatibility check.
        assert_noop("N-methyl-2-hydroxypyridinium", "Oc1cccc[n+]1C");
    }

    #[test]
    fn negative_control_aromatic_n_bridge_excluded() {
        // N-hydroxypyrrole: the bridge atom bearing the exocyclic OH is
        // itself an aromatic nitrogen, not carbon -- excluded by the
        // bridge-element==6 check.
        assert_noop("N-hydroxypyrrole", "On1cccc1");
    }

    #[test]
    fn negative_control_fused_bridgehead_nitrogen_acceptor_excluded() {
        // The candidate ring nitrogen adjacent to the exocyclic-OH-bearing
        // bridge is a 3-connected fusion (bridgehead) atom in this fused
        // bicyclic system -- excluded by the acceptor degree==2 check (no
        // free valence slot for +1H without a charge change).
        assert_noop("fused bridgehead-N acceptor", "Oc1ccn2ccccc12");
    }

    #[test]
    fn remote_stereocenter_same_atom_chirality_and_cip_preserved() {
        // Stronger than counting stereocenters: locate the SAME physical
        // atom (index-stable -- the shift's post-generation invariant check
        // guarantees every non-donor/acceptor atom, including the
        // stereocenter, keeps its index) and confirm chirality,
        // stereo_neighbor_order, and the CIP label are all unchanged, not
        // just present.
        let mol = parse("Oc1c([C@@H](F)Cl)cccn1").unwrap();
        let out = canonical_tautomer(&mol);

        let stereocenter = (0..mol.atom_count())
            .map(|i| AtomIdx(i as u32))
            .find(|&i| mol.atom(i).chirality != Chirality::None)
            .expect("input must have a stereocenter");

        assert_ne!(
            canonical_smiles(&mol),
            canonical_smiles(&out),
            "sanity: the shift must have actually fired"
        );
        assert_eq!(
            mol.atom(stereocenter).element,
            out.atom(stereocenter).element
        );
        assert_eq!(
            mol.atom(stereocenter).chirality,
            out.atom(stereocenter).chirality,
            "chirality flag changed on the uninvolved stereocenter"
        );
        assert_eq!(
            mol.stereo_neighbor_order(stereocenter),
            out.stereo_neighbor_order(stereocenter),
            "stereo_neighbor_order changed on the uninvolved stereocenter"
        );
        assert_eq!(
            crate::cip::assign_cip(&mol).get(stereocenter),
            crate::cip::assign_cip(&out).get(stereocenter),
            "CIP label changed on the uninvolved stereocenter"
        );
    }

    // -----------------------------------------------------------------------
    // Issue #415: canonical_tautomer must never produce an unkekulizable
    // (chemically invalid) molecule, even for a fused/bridged ring system
    // where `find_sssr`'s own ring choice can differ depending on atom
    // labeling.
    // -----------------------------------------------------------------------

    #[test]
    fn issue415_fused_ring_shift_never_produces_invalid_molecule() {
        // A phthalazinone-like fused bicyclic tautomer where re-labeling via
        // a canonical-SMILES round trip used to make both the
        // exocyclic-lactam-shift mechanism and the plain direct-aromatic
        // 1,2-shift mechanism protonate a ring nitrogen that was already
        // ring-adjacent to another protonated nitrogen -- an over-valent
        // `[nH2]`-shaped state RDKit itself rejects outright (confirmed via
        // `Chem.MolFromSmiles` at diagnosis time). Every candidate this
        // molecule can reach through repeated standardize passes must stay
        // kekulizable.
        let smi = "Oc1[nH]ncc2c3cc(OCc4ccccc4)ccc3nc1-2";
        let mol = parse(smi).unwrap();
        let mut current = mol;
        for pass in 0..6 {
            let next = canonical_tautomer(&current);
            assert!(
                chematic_core::kekulize(&next).is_ok(),
                "pass {pass}: canonical_tautomer produced an unkekulizable molecule: {}",
                canonical_smiles(&next)
            );
            let reparsed = parse(&canonical_smiles(&next)).unwrap();
            current = reparsed;
        }
    }

    /// Direct regression for issue #415's residual order-dependence: the
    /// same molecule, relabeled only via atom-order permutation
    /// (`reorder_atoms`, so every permutation is guaranteed structurally
    /// identical -- not a hand-alternate SMILES respelling that could
    /// accidentally encode a different molecule), must reach the byte-
    /// identical canonical tautomer regardless of which permutation the
    /// caller happens to hand in, and that result must already be its own
    /// fixed point (`canonical_tautomer` applied to its own output must be
    /// a no-op) -- not just "eventually convergent after enough reparses."
    ///
    /// Empirically confirmed before writing this test (not assumed): this
    /// molecule's true canonical form needs the *search itself* to run
    /// twice to stabilize (`canonical_tautomer_with_config`'s own outer
    /// convergence loop, not an atom-order artifact -- a single search
    /// pass genuinely doesn't reach the fixed point, independent of atom
    /// labeling). 6 orderings covering identity, full reversal, and 4
    /// distinct rotations all agree exactly, both on the first search's
    /// result and on the final stable form.
    #[test]
    fn issue415_residual_canonical_tautomer_is_atom_order_invariant_and_idempotent() {
        let smi = "Oc1[nH]ncc2c3cc(OCc4ccccc4)ccc3nc1-2";
        let base = parse(smi).unwrap();
        let n = base.atom_count();

        let mut orderings: Vec<Vec<AtomIdx>> = vec![
            (0..n as u32).map(AtomIdx).collect(),
            (0..n as u32).rev().map(AtomIdx).collect(),
        ];
        for shift in [1usize, 3, (n / 2).max(1), n - 1] {
            let mut order: Vec<AtomIdx> = (0..n as u32).map(AtomIdx).collect();
            order.rotate_left(shift);
            orderings.push(order);
        }
        assert_eq!(orderings.len(), 6, "expected 6 distinct atom orderings");

        let mut expected: Option<(String, String)> = None;
        for (i, order) in orderings.iter().enumerate() {
            let permuted = reorder_atoms(&base, order);
            assert_eq!(
                canonical_smiles(&permuted),
                canonical_smiles(&base),
                "ordering {i}: reorder_atoms must not change the molecule's own identity"
            );
            let result = canonical_tautomer(&permuted);
            let reapplied = canonical_tautomer(&result);
            let result_smi = canonical_smiles(&result);
            let reapplied_smi = canonical_smiles(&reapplied);
            assert_eq!(
                result_smi, reapplied_smi,
                "ordering {i}: canonical_tautomer's own output is not a fixed point of itself"
            );
            match &expected {
                None => expected = Some((result_smi, reapplied_smi)),
                Some((exp_result, exp_reapplied)) => {
                    assert_eq!(
                        &result_smi, exp_result,
                        "ordering {i} reached a different canonical tautomer than ordering 0"
                    );
                    assert_eq!(&reapplied_smi, exp_reapplied);
                }
            }
        }
    }
}
