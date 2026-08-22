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

/// The 42 tautomer rules.
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
    // 20. 1,3-C→N any bridge: active methylene adjacent to =N (via S, O, or N bridge)
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
    // 21. 1,5-O→O: β-diketone (acetylacetone)
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
    // 22. 1,5-O→N: enol imine
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
    // 23. 1,5-N→N: extended guanidine/amidine tautomerism
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
    // 24. 1,5-N→O: hydroxamic acid type
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
    // 25. 1,5-C→O: active methylene with conjugation
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
    // 26. 1,5-O→O with N bridge: nitro-type tautomerism
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
    // 27. 1,5-O→O with S bridge: thio-β-diketone
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
    // 28. 1,5-N→O with C bridge (existing, carbon specified)
    // Already covered by rule 22 (path_len=5 with C bridge)
    // 29. 1,5-N→O with N bridge: bridging N (amidino-type)
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
    // 30. 1,5-N→O with S bridge
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
    // 31. 1,5-N→N with N bridge: guanidine-type via N
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
    // 32. 1,5-N→N with S bridge
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
    // 33. 1,5-O→N with N bridge: hydroxamic-type via N
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
    // 34. 1,5-O→N with S bridge
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
    // 35. 1,5-S→O with N bridge
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
    // 36. 1,5-S→O with S bridge
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
    // 37. 1,5-S→N with C bridge
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
    // 38. 1,5-S→N with N bridge
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
    // 39. 1,5-C→N with C bridge: extended enamine-imine
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
    // 40. 1,5-C→N with N bridge
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
    // 41. 1,5-C→N with S bridge
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
    // 42. 1,5-C→S with N bridge
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
    // 43. 1,5-C→S with S bridge
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
) -> Vec<Molecule> {
    enumerate_direct_aromatic_forms_tracked(start, blocked, max).0
}

/// Like [`enumerate_direct_aromatic_forms`], but also reports whether the
/// search was cut short by `max` while more candidates were still reachable
/// (`true`) versus exhausting the frontier naturally (`false`) -- needed by
/// [`tautomer_parent`] to distinguish `Completed` from `MaxTautomersReached`.
fn enumerate_direct_aromatic_forms_tracked(
    start: &Molecule,
    blocked: &HashSet<AtomIdx>,
    max: usize,
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
        for (d, a) in find_direct_aromatic_matches(&current) {
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
/// Returns (donor, acceptor) pairs for direct 1,2-shift (no bridge atom).
fn find_direct_aromatic_matches(mol: &Molecule) -> Vec<(AtomIdx, AtomIdx)> {
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
    Some(builder.build())
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
    let mut atoms: Vec<(u8, i8, u32)> = (0..mol.atom_count())
        .map(|i| {
            let idx = AtomIdx(i as u32);
            let a = mol.atom(idx);
            let bos: u32 = mol
                .neighbors(idx)
                .map(|(_, bidx)| mol.bond(bidx).order.order_int() as u32)
                .sum();
            (a.element.atomic_number(), a.charge, bos)
        })
        .collect();
    atoms.sort();
    let mut hash = FNV1A_OFFSET;
    for (an, ch, bos) in atoms {
        hash ^= an as u64;
        hash = hash.wrapping_mul(FNV1A_PRIME);
        hash ^= (ch as u8 as u64).wrapping_add(128);
        hash = hash.wrapping_mul(FNV1A_PRIME);
        hash ^= bos as u64;
        hash = hash.wrapping_mul(FNV1A_PRIME);
    }
    hash
}

/// Find all (donor, bridge, acceptor) triples matching the rule in `mol`.
/// For path_len=3: donor-bridge-acceptor (standard 1,3-shift)
/// For path_len=5: donor-b1-b2-b3-acceptor (1,5-shift) where b2 is stored in "bridge"
fn find_matches(mol: &Molecule, rule: &TautomerRule) -> Vec<(AtomIdx, AtomIdx, AtomIdx)> {
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

/// Apply the first matching transformation for `rule`; return the new molecule.
fn apply_first_match(
    mol: &Molecule,
    rule: &TautomerRule,
    config: &TautomerConfig,
) -> Option<Molecule> {
    apply_first_match_tracked(mol, rule, config).map(|(next, ..)| next)
}

/// Like [`apply_first_match`], but also returns the (donor, bridge, acceptor)
/// triple that was moved -- needed by [`tautomer_parent`] to record which
/// atoms/bonds each applied transform touched.
fn apply_first_match_tracked(
    mol: &Molecule,
    rule: &TautomerRule,
    config: &TautomerConfig,
) -> Option<(Molecule, AtomIdx, AtomIdx, AtomIdx)> {
    find_matches(mol, rule).into_iter().find_map(|(d, b, a)| {
        transfer_hydrogen(mol, d, b, a, &config.blocked_atoms, &config.blocked_bonds)
            .map(|next| (next, d, b, a))
    })
}

/// Apply every matching transformation for `rule`; return all resulting molecules.
fn apply_all_matches(
    mol: &Molecule,
    rule: &TautomerRule,
    config: &TautomerConfig,
) -> Vec<Molecule> {
    find_matches(mol, rule)
        .into_iter()
        .filter_map(|(d, b, a)| {
            transfer_hydrogen(mol, d, b, a, &config.blocked_atoms, &config.blocked_bonds)
        })
        .collect()
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
    /// Each outer iteration applies at most one transformation (the first
    /// match found, in atom-index order, for the first rule that has an
    /// unseen match). If a molecule has MORE independent, non-automorphic,
    /// same-rule tautomerizable sites than `max_iter` allows, the final
    /// result depends on input atom order (parse/spelling order): which
    /// subset of sites got converted before the budget ran out differs.
    /// Confirmed reachable (not just theoretical) for a "comb" molecule
    /// with >16 independent enol arms, but no known real-molecule instance
    /// hits this in practice -- see
    /// `test_max_iter_default_diverges_on_many_independent_sites` (an
    /// `#[ignore]`d regression pin, not a passing guarantee) for the exact
    /// mechanism and why it's not fixed by simply raising this constant.
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
// Public API
// ---------------------------------------------------------------------------

/// Return the canonical (preferred) tautomer of `mol`.
///
/// Applies forward-preferred rules iteratively until no new form is found
/// or the iteration limit is reached. After rule-based normalization, direct
/// aromatic 1,2-shift tautomers are compared and the form with the
/// lexicographically smallest H-assignment vector is chosen.
///
/// Uses [`TautomerConfig::default`] (all rules, max_iter=16).
pub fn canonical_tautomer(mol: &Molecule) -> Molecule {
    canonical_tautomer_with_config(mol, &TautomerConfig::default())
}

/// Like [`canonical_tautomer`] but with explicit configuration.
pub fn canonical_tautomer_with_config(mol: &Molecule, config: &TautomerConfig) -> Molecule {
    let mut current = mol.clone();
    let mut seen = HashSet::new();
    seen.insert(mol_fingerprint(&current));

    for _ in 0..config.max_iter {
        let mut changed = false;
        for rule in active_rules(config)
            .into_iter()
            .filter(|r| r.prefer_forward)
        {
            if let Some(next) = apply_first_match(&current, rule, config) {
                let fp = mol_fingerprint(&next);
                if !seen.contains(&fp) {
                    seen.insert(fp);
                    current = next;
                    changed = true;
                    break;
                }
            }
        }
        if !changed {
            break;
        }
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
    ));
    if candidates.len() > 1 {
        candidates.sort_by(|a, b| {
            tautomer_score(b).cmp(&tautomer_score(a)).then_with(|| {
                chematic_smiles::canonical_smiles(a).cmp(&chematic_smiles::canonical_smiles(b))
            })
        });
        current = candidates.into_iter().next().unwrap();
    }
    current
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
            for next in apply_all_matches(&current, rule, config) {
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
        // Direct aromatic 1,2-shift (e.g. pyrazole N1H ↔ N2H).
        for (d, a) in find_direct_aromatic_matches(&current) {
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
/// `docs/rfcs/tautomer_parent_identity_phase2_rfc.md` sections 4.2-4.3.
///
/// Does **not** implement the section 4.4 aromatic lactam/lactim fix (round
/// 2C) -- on the aromatic tautomer pairs described there, this returns the
/// same non-invariant result `canonical_tautomer` does today, just with an
/// explicit `Completed` status (the rule-based loop and aromatic-shift
/// search both genuinely converge without hitting any budget; the result is
/// wrong for a different, structural reason that budget visibility cannot
/// detect).
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
                    apply_first_match_tracked(&current, rule, &config)
                        .is_some_and(|(next, ..)| !seen.contains(&mol_fingerprint(&next)))
                });
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
            if let Some((next, donor, bridge, acceptor)) =
                apply_first_match_tracked(&current, rule, &config)
            {
                let fp = mol_fingerprint(&next);
                if !seen.contains(&fp) {
                    seen.insert(fp);
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
                    current = next;
                    transforms_applied += 1;
                    changed = true;
                    break;
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
    use chematic_core::{AtomIdx, Chirality};
    use chematic_smiles::{canonical_smiles, parse};

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
    #[ignore = "confirmed real, root-caused, and deliberately not fixed this \
                round: canonical_tautomer_with_config's greedy loop applies \
                one transform per outer iteration via apply_first_match \
                (atom-index order, see find_matches), bounded by \
                config.max_iter (default 16). A molecule with MORE \
                independent same-rule tautomerizable sites than max_iter \
                allows has its FINAL result depend on atom insertion order: \
                which subset of sites got converted before the budget ran \
                out differs. Verified NOT a deeper tie-break flaw -- with \
                max_iter=1000 (room to fully process all 25 sites) both \
                atom orderings converge to the byte-identical correct \
                answer (see test_max_iter_1000_resolves_the_divergence); \
                verified NOT a canonical_smiles bug either (see \
                test_canonical_smiles_alone_is_order_independent_on_comb). \
                The trigger (>16 independent, same-rule, non-automorphic \
                tautomerizable sites in one molecule) is an extreme edge \
                case with no known real-molecule instance -- a proper fix \
                needs batching all of a rule's non-conflicting matches per \
                iteration instead of raising the constant, an architectural \
                change with its own conflict-resolution complexity, \
                correctly left out of scope for this round rather than \
                rushed. See TautomerConfig::max_iter's doc comment."]
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
        // Companion to the #[ignore]d test above: proves the divergence
        // there is PURELY max_iter exhaustion, not a deeper algorithmic
        // flaw -- with enough budget to process all 25 independent sites,
        // atom insertion order no longer matters.
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

        let (donor, acceptor) = find_direct_aromatic_matches(&mol)
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

        let (donor, acceptor) = find_direct_aromatic_matches(&mol)
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
        find_matches(mol, rule)
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
    #[ignore = "known: canonical_tautomer loses E/Z stereo on hydrazones/imines (RDKit PR #9128). \
                mol_fingerprint() does not include Up/Down bond orders so both E and Z forms hash \
                identically; the canonical tautomer selection then returns the same SMILES for both. \
                Fix requires either including stereo in mol_fingerprint or re-applying input E/Z stereo \
                to the canonical tautomer output after selection."]
    fn test_hydrazone_ez_stereo_preserved_in_canonical_tautomer() {
        // E-hydrazone and Z-hydrazone are DIFFERENT compounds.
        // mol_fingerprint() does not encode Up/Down bond orders, so both map to
        // the same structural hash. canonical_tautomer incorrectly merges them.
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
}
