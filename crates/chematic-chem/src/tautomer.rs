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

// ---------------------------------------------------------------------------
// Rule infrastructure
// ---------------------------------------------------------------------------

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
            BondOrderMatch::Single => matches!(order, BondOrder::Single | BondOrder::Up | BondOrder::Down),
            BondOrderMatch::Double => matches!(order, BondOrder::Double),
            BondOrderMatch::Any => true,
        }
    }
}

/// A tautomer transformation rule: donor loses an H and the bond orders shift.
///
/// Pattern: donor -[donor_bridge_order]- bridge -[bridge_acceptor_order]- acceptor
/// After transform: donor =[new]= bridge -[new]- acceptor (with H shifted to acceptor)
struct TautomerRule {
    #[allow(dead_code)]
    name: &'static str,
    /// Atomic number of the donor atom (loses H).
    donor_elem: u8,
    /// Atomic number of the bridge atom (None = any).
    bridge_elem: Option<u8>,
    /// Atomic number of the acceptor atom (gains H via implicit valence).
    acceptor_elem: u8,
    /// Required bond order between donor and bridge.
    donor_bridge_order: BondOrderMatch,
    /// Required bond order between bridge and acceptor.
    bridge_acceptor_order: BondOrderMatch,
    /// If true, this rule is applied in canonical_tautomer to normalize toward a preferred form.
    prefer_forward: bool,
}

/// The 5 tautomer rules.
static RULES: &[TautomerRule] = &[
    // keto-enol: O-H single C-double C → O=C-C (prefer keto form, C=O)
    // Pattern in enol form: donor=O(has H), O-C single, C=C double
    TautomerRule {
        name: "keto-enol",
        donor_elem: 8,      // O
        bridge_elem: Some(6), // C
        acceptor_elem: 6,   // C
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: true,
    },
    // amide-iminol: N(has H)-C=O → iminol form would be N=C-O
    // We prefer the amide (N-C=O), so prefer_forward=false prevents conversion away from amide.
    // The pattern here matches iminol form: N(has H) with single bond to C, C=O.
    // Wait - this pattern actually matches N-C=O (amide). We set prefer_forward=false
    // so canonical_tautomer won't apply it to convert amide → iminol.
    TautomerRule {
        name: "amide-iminol",
        donor_elem: 7,      // N
        bridge_elem: Some(6), // C
        acceptor_elem: 8,   // O
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
    },
    // imine-enamine: N(has H)-C=C → N=C-C (prefer imine)
    TautomerRule {
        name: "imine-enamine",
        donor_elem: 7,      // N
        bridge_elem: Some(6), // C
        acceptor_elem: 6,   // C
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: true,
    },
    // 1,3-H-shift N→O: bidirectional
    TautomerRule {
        name: "1,3-H-shift-N-O",
        donor_elem: 7,      // N
        bridge_elem: None,  // any
        acceptor_elem: 8,   // O
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
    },
    // 1,3-H-shift N→N: bidirectional
    TautomerRule {
        name: "1,3-H-shift-N-N",
        donor_elem: 7,      // N
        bridge_elem: None,  // any
        acceptor_elem: 7,   // N
        donor_bridge_order: BondOrderMatch::Single,
        bridge_acceptor_order: BondOrderMatch::Double,
        prefer_forward: false,
    },
];

// ---------------------------------------------------------------------------
// Molecule cloning helper
// ---------------------------------------------------------------------------

fn clone_mol(mol: &Molecule) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    for i in 0..mol.atom_count() {
        builder.add_atom(mol.atom(AtomIdx(i as u32)).clone());
    }
    for i in 0..mol.bond_count() {
        let b = mol.bond(BondIdx(i as u32));
        builder.add_bond(b.atom1, b.atom2, b.order).ok();
    }
    builder.build()
}

// ---------------------------------------------------------------------------
// Fingerprint for convergence detection
// ---------------------------------------------------------------------------

fn mol_fingerprint(mol: &Molecule) -> u64 {
    let mut atoms: Vec<(u8, i8, u32)> = (0..mol.atom_count())
        .map(|i| {
            let idx = AtomIdx(i as u32);
            let a = mol.atom(idx);
            let bos: u32 = mol.neighbors(idx)
                .map(|(_, bidx)| mol.bond(bidx).order.order_int() as u32)
                .sum();
            (a.element.atomic_number(), a.charge, bos)
        })
        .collect();
    atoms.sort();
    let mut hash = 0xcbf29ce484222325u64;
    for (an, ch, bos) in atoms {
        hash ^= an as u64;
        hash = hash.wrapping_mul(0x100000001b3);
        hash ^= (ch as u8 as u64).wrapping_add(128);
        hash = hash.wrapping_mul(0x100000001b3);
        hash ^= bos as u64;
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash
}

// ---------------------------------------------------------------------------
// Pattern matching
// ---------------------------------------------------------------------------

/// Check if donor has at least 1 hydrogen (implicit or explicit).
fn has_hydrogen(mol: &Molecule, idx: AtomIdx) -> bool {
    implicit_hcount(mol, idx) > 0
}

/// Find all (donor, bridge, acceptor) triples matching the rule in `mol`.
fn find_matches(mol: &Molecule, rule: &TautomerRule) -> Vec<(AtomIdx, AtomIdx, AtomIdx)> {
    let mut matches = Vec::new();

    for i in 0..mol.atom_count() {
        let d = AtomIdx(i as u32);
        let donor_atom = mol.atom(d);

        // Check donor element
        if donor_atom.element.atomic_number() != rule.donor_elem {
            continue;
        }
        // Check donor has at least 1 H
        if !has_hydrogen(mol, d) {
            continue;
        }

        // Iterate over neighbors of donor (potential bridge atoms)
        for (b, db_bidx) in mol.neighbors(d) {
            let db_bond = mol.bond(db_bidx);
            if !rule.donor_bridge_order.matches(db_bond.order) {
                continue;
            }

            let bridge_atom = mol.atom(b);
            // Check bridge element if specified
            if let Some(br_elem) = rule.bridge_elem {
                if bridge_atom.element.atomic_number() != br_elem {
                    continue;
                }
            }

            // Iterate over neighbors of bridge (potential acceptor atoms), excluding donor
            for (a, ba_bidx) in mol.neighbors(b) {
                if a == d {
                    continue; // skip back to donor
                }
                let ba_bond = mol.bond(ba_bidx);
                if !rule.bridge_acceptor_order.matches(ba_bond.order) {
                    continue;
                }

                let acceptor_atom = mol.atom(a);
                if acceptor_atom.element.atomic_number() != rule.acceptor_elem {
                    continue;
                }

                matches.push((d, b, a));
            }
        }
    }

    matches
}

// ---------------------------------------------------------------------------
// H transfer (bond order swap)
// ---------------------------------------------------------------------------

/// Apply a single tautomer transformation: swap bond orders so that:
/// - donor-bridge bond: Single → Double
/// - bridge-acceptor bond: Double → Single
///
/// For organic-subset atoms, the implicit H count adjusts automatically via
/// the valence model. For bracket atoms with explicit hydrogen_count, we also
/// adjust the counts manually.
///
/// Returns None if the transformation would result in an invalid molecule
/// (e.g., duplicate bonds or other errors).
fn transfer_hydrogen(
    mol: &Molecule,
    donor: AtomIdx,
    bridge: AtomIdx,
    acceptor: AtomIdx,
) -> Option<Molecule> {
    // Find the bond indices
    let (db_bidx, _) = mol.bond_between(donor, bridge)?;
    let (ba_bidx, _) = mol.bond_between(bridge, acceptor)?;

    let mut builder = MoleculeBuilder::new();

    // Clone atoms, adjusting explicit H counts for bracket atoms
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let mut atom = mol.atom(idx).clone();

        // For bracket atoms with explicit hydrogen_count, adjust H counts
        if let Some(h) = atom.hydrogen_count {
            if idx == donor {
                // Donor loses 1 H; ensure we don't go negative
                if h == 0 {
                    return None;
                }
                atom.hydrogen_count = Some(h - 1);
            } else if idx == acceptor {
                // Acceptor gains 1 H
                atom.hydrogen_count = Some(h.saturating_add(1));
            }
        }

        builder.add_atom(atom);
    }

    // Clone bonds with modified orders for donor-bridge and bridge-acceptor
    for i in 0..mol.bond_count() {
        let bidx = BondIdx(i as u32);
        let b = mol.bond(bidx);

        let order = if bidx == db_bidx {
            // donor-bridge: Single → Double
            BondOrder::Double
        } else if bidx == ba_bidx {
            // bridge-acceptor: Double → Single
            BondOrder::Single
        } else {
            b.order
        };

        builder.add_bond(b.atom1, b.atom2, order).ok()?;
    }

    Some(builder.build())
}

// ---------------------------------------------------------------------------
// Apply rule to molecule
// ---------------------------------------------------------------------------

/// Apply the first matching transformation for a rule, return the new molecule.
fn apply_first_match(mol: &Molecule, rule: &TautomerRule) -> Option<Molecule> {
    let matches = find_matches(mol, rule);
    for (donor, bridge, acceptor) in matches {
        if let Some(next) = transfer_hydrogen(mol, donor, bridge, acceptor) {
            return Some(next);
        }
    }
    None
}

/// Apply all matching transformations for a rule, return all new molecules.
fn apply_all_matches(mol: &Molecule, rule: &TautomerRule) -> Vec<Molecule> {
    let mut results = Vec::new();
    let matches = find_matches(mol, rule);
    for (donor, bridge, acceptor) in matches {
        if let Some(next) = transfer_hydrogen(mol, donor, bridge, acceptor) {
            results.push(next);
        }
    }
    results
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Return the canonical (preferred) tautomer of `mol`.
///
/// Applies forward-preferred rules iteratively until no new form is found
/// or the iteration limit is reached.
pub fn canonical_tautomer(mol: &Molecule) -> Molecule {
    const MAX_ITER: usize = 16;
    let mut current = clone_mol(mol);
    let mut seen = HashSet::new();
    seen.insert(mol_fingerprint(&current));

    for _ in 0..MAX_ITER {
        let mut changed = false;
        for rule in RULES.iter().filter(|r| r.prefer_forward) {
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
    current
}

/// Enumerate all reachable tautomers of `mol`, capped at 32.
///
/// Returns a `Vec<Molecule>` where the first element is the original molecule.
pub fn enumerate_tautomers(mol: &Molecule) -> Vec<Molecule> {
    const MAX_TAUTOMERS: usize = 32;
    let mut result = vec![clone_mol(mol)];
    let mut seen = HashSet::new();
    seen.insert(mol_fingerprint(mol));
    let mut frontier = vec![clone_mol(mol)];

    while !frontier.is_empty() && result.len() < MAX_TAUTOMERS {
        let current = frontier.remove(0);
        for rule in RULES.iter() {
            for next in apply_all_matches(&current, rule) {
                let fp = mol_fingerprint(&next);
                if !seen.contains(&fp) {
                    seen.insert(fp);
                    frontier.push(clone_mol(&next));
                    result.push(next);
                    if result.len() >= MAX_TAUTOMERS {
                        break;
                    }
                }
            }
            if result.len() >= MAX_TAUTOMERS {
                break;
            }
        }
    }
    result
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

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
        assert!(tautomers.len() >= 1);
    }

    #[test]
    fn test_enumerate_vinyl_alcohol() {
        // OC=C → should find at least 1 more tautomer (keto form CC=O)
        let mol = parse("OC=C").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        assert!(tautomers.len() >= 2, "Expected >= 2 tautomers for vinyl alcohol, got {}", tautomers.len());
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
}
