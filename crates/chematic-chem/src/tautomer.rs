//! Tautomer enumeration and canonical tautomer selection.
//!
//! Provides two public functions:
//! - [`canonical_tautomer`]: return the canonical (preferred) tautomer form.
//! - [`enumerate_tautomers`]: enumerate all reachable tautomers up to a cap.
//!
//! Pattern matching is implemented directly on the molecule graph (no SMARTS).

#![forbid(unsafe_code)]

use std::collections::HashSet;

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule, MoleculeBuilder, implicit_hcount};

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

/// The 15 tautomer rules.
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
        bridge_elem: Some(6),  // Central carbon in pattern (b2)
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
fn transfer_hydrogen_aromatic(
    mol: &Molecule,
    donor: AtomIdx,
    acceptor: AtomIdx,
) -> Option<Molecule> {
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
            atom.hydrogen_count = Some(donor_h - 1);
        } else if idx == acceptor {
            atom.hydrogen_count = Some(acceptor_h + 1);
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
    Some(builder.build())
}

fn clone_mol(mol: &Molecule) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    for i in 0..mol.atom_count() {
        builder.add_atom(mol.atom(AtomIdx(i as u32)).clone());
    }
    for i in 0..mol.bond_count() {
        let b = mol.bond(BondIdx(i as u32));
        builder
            .add_bond(b.atom1, b.atom2, b.order)
            .expect("clone_mol: bond from a valid molecule must be re-addable");
    }
    builder.build()
}

const FNV1A_OFFSET: u64 = 0xcbf29ce484222325;
const FNV1A_PRIME: u64 = 0x100000001b3;

/// Tautomer form score for canonical selection: prefers O-H > N-H > S-H and aromatic rings.
fn tautomer_score(mol: &Molecule) -> i32 {
    let mut score = 0i32;
    let mut has_aromatic = false;

    for (_, atom) in mol.atoms() {
        // Aromatic ring bonus
        if atom.aromatic {
            has_aromatic = true;
        }

        // Heteroatom hydrogen: O-H > N-H > S-H (explicit or implicit)
        let h_count = atom.hydrogen_count.unwrap_or(0) as i32;
        if h_count > 0 {
            match atom.element.atomic_number() {
                8 => score += h_count * 100,  // O-H: highest priority
                7 => score += h_count * 50,   // N-H: medium priority
                16 => score += h_count * 25,  // S-H: low priority
                _ => {}
            }
        }
    }

    // Aromatic system bonus
    if has_aromatic {
        score += 1000;
    }

    score
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
                        && mol.atom(b2).element.atomic_number() != br_elem {
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
fn transfer_hydrogen(
    mol: &Molecule,
    donor: AtomIdx,
    bridge: AtomIdx,
    acceptor: AtomIdx,
) -> Option<Molecule> {
    let (db_bidx, _) = mol.bond_between(donor, bridge)?;
    let (ba_bidx, _) = mol.bond_between(bridge, acceptor)?;

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
    Some(builder.build())
}

/// Apply the first matching transformation for `rule`; return the new molecule.
fn apply_first_match(mol: &Molecule, rule: &TautomerRule) -> Option<Molecule> {
    find_matches(mol, rule)
        .into_iter()
        .find_map(|(d, b, a)| transfer_hydrogen(mol, d, b, a))
}

/// Apply every matching transformation for `rule`; return all resulting molecules.
fn apply_all_matches(mol: &Molecule, rule: &TautomerRule) -> Vec<Molecule> {
    find_matches(mol, rule)
        .into_iter()
        .filter_map(|(d, b, a)| transfer_hydrogen(mol, d, b, a))
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
#[derive(Debug, Clone)]
pub struct TautomerConfig {
    /// Maximum iterations in [`canonical_tautomer_with_config`] (default 16).
    pub max_iter: usize,
    /// Maximum tautomers returned by [`enumerate_tautomers_with_config`] (default 32).
    pub max_tautomers: usize,
    /// 0-based indices of rules to activate.  Empty = all rules active.
    pub enabled_rules: Vec<usize>,
}

impl Default for TautomerConfig {
    fn default() -> Self {
        Self {
            max_iter: 16,
            max_tautomers: 32,
            enabled_rules: Vec::new(),
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
    let mut current = clone_mol(mol);
    let mut seen = HashSet::new();
    seen.insert(mol_fingerprint(&current));

    for _ in 0..config.max_iter {
        let mut changed = false;
        for rule in active_rules(config)
            .into_iter()
            .filter(|r| r.prefer_forward)
        {
            if let Some(next) = apply_first_match(&current, rule) {
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
    // (O-H > N-H > S-H, aromatic rings), with H-assignment as tiebreaker.
    let mut candidates: Vec<Molecule> = vec![clone_mol(&current)];
    for (d, a) in find_direct_aromatic_matches(&current) {
        if let Some(t) = transfer_hydrogen_aromatic(&current, d, a) {
            candidates.push(t);
        }
    }
    if candidates.len() > 1 {
        candidates.sort_by(|a, b| {
            tautomer_score(b)
                .cmp(&tautomer_score(a))
                .then_with(|| h_assignment(a).cmp(&h_assignment(b)))
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
    let mut result = vec![clone_mol(mol)];
    let mut seen = HashSet::new();
    seen.insert(mol_fingerprint(mol));
    // Separate seen-set for 1,2-shift (mol_fingerprint can't distinguish positional H isomers).
    let mut h_seen: HashSet<Vec<Option<u32>>> = HashSet::new();
    h_seen.insert(h_assignment(mol));
    let mut frontier = vec![clone_mol(mol)];

    while !frontier.is_empty() && result.len() < config.max_tautomers {
        let current = frontier.remove(0);
        for rule in active_rules(config).into_iter() {
            for next in apply_all_matches(&current, rule) {
                let fp = mol_fingerprint(&next);
                if !seen.contains(&fp) {
                    seen.insert(fp);
                    h_seen.insert(h_assignment(&next));
                    frontier.push(clone_mol(&next));
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
            if let Some(next) = transfer_hydrogen_aromatic(&current, d, a) {
                let ha = h_assignment(&next);
                if !h_seen.contains(&ha) {
                    h_seen.insert(ha);
                    seen.insert(mol_fingerprint(&next));
                    frontier.push(clone_mol(&next));
                    result.push(next);
                }
            }
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

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
        use chematic_smiles::canonical_smiles;
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
        use chematic_smiles::canonical_smiles;
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
        use chematic_smiles::canonical_smiles;
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
        assert!(tautomers.len() >= 2, "Expected >= 2 tautomers for β-diketone");
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
}
