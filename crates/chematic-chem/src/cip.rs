//! CIP (Cahn–Ingold–Prelog) stereochemistry assignment.
//!
//! Implements R/S (tetrahedral) and E/Z (double bond) assignment for
//! molecules parsed from SMILES with chirality annotations.

use std::collections::{HashMap, HashSet, VecDeque};

use chematic_core::{AtomIdx, BondIdx, BondOrder, CipCode, Chirality, Molecule};

/// The result of a CIP stereochemistry assignment run.
#[derive(Debug)]
pub struct CipAssignment {
    pub assignments: Vec<(AtomIdx, CipCode)>,
}

impl CipAssignment {
    /// Look up the CIP code for a given atom index.
    pub fn get(&self, idx: AtomIdx) -> Option<CipCode> {
        self.assignments
            .iter()
            .find(|(i, _)| *i == idx)
            .map(|(_, c)| *c)
    }
}

/// Run CIP assignment on `mol`.  Returns R/S for chiral tetrahedral centers
/// and E/Z for stereospecified double bonds.
pub fn assign_cip(mol: &Molecule) -> CipAssignment {
    let mut assignments = Vec::new();

    // R/S for tetrahedral centers
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        if let Some(code) = assign_tetrahedral(mol, idx) {
            assignments.push((idx, code));
        }
    }

    // E/Z for double bonds
    for j in 0..mol.bond_count() {
        let bidx = BondIdx(j as u32);
        if let Some((atom_idx, code)) = assign_ez(mol, bidx) {
            assignments.push((atom_idx, code));
        }
    }

    CipAssignment { assignments }
}

/// A single "sphere layer" in a CIP branch expansion: a sorted list of
/// `(atomic_num, isotope)` pairs (sorted descending for lexicographic comparison).
type SphereLayer = Vec<(u8, Option<u16>)>;

/// Get the key `(atomic_num, isotope)` for an atom, used in CIP comparisons.
///
/// For the virtual H sentinel (`AtomIdx(u32::MAX)`), returns `(1, None)`.
fn atom_key(mol: &Molecule, idx: AtomIdx) -> (u8, Option<u16>) {
    if idx.0 == u32::MAX {
        return (1, None);
    }
    let a = mol.atom(idx);
    (a.element.atomic_number(), a.isotope)
}

/// Compare two `(atomic_num, isotope)` keys by CIP priority.
///
/// Higher atomic number wins.  For isotopes: `Some(n)` > `None` (heavier beats
/// unspecified, per CIP rule 2).
fn cmp_key(a: (u8, Option<u16>), b: (u8, Option<u16>)) -> std::cmp::Ordering {
    use std::cmp::Ordering::*;
    match a.0.cmp(&b.0) {
        Equal => match (a.1, b.1) {
            (Some(x), Some(y)) => x.cmp(&y),
            (Some(_), None) => Greater,
            (None, Some(_)) => Less,
            (None, None) => Equal,
        },
        other => other,
    }
}

/// BFS state for one node during sphere expansion.
struct ExpandState {
    node: AtomIdx,
    parent: AtomIdx,
    depth: usize,
    visited: HashSet<AtomIdx>,
}

/// BFS-based CIP sphere expansion for the branch starting at `start`,
/// not going back through `center`.
///
/// At each depth layer the collected `(atomic_num, isotope)` tuples are sorted
/// descending (highest priority first), implementing the hierarchical digraph
/// comparison.
///
/// Phantom atom rules:
/// 1. **Double-bond phantom**: when expanding node B (reached via double bond from A),
///    add a phantom entry for A at the same depth level as B's children.
/// 2. **Ring revisit phantom**: if an already-visited atom is encountered,
///    add a phantom for it but don't expand further.
fn cip_branch_spheres(mol: &Molecule, center: AtomIdx, start: AtomIdx) -> Vec<SphereLayer> {
    let mut layers: HashMap<usize, Vec<(u8, Option<u16>)>> = HashMap::new();
    let max_depth = 8usize;

    // The start atom itself is at depth 1.
    let start_key = atom_key(mol, start);
    layers.entry(1).or_default().push(start_key);

    let mut expand_queue: VecDeque<ExpandState> = VecDeque::new();
    {
        let mut v = HashSet::new();
        v.insert(center);
        v.insert(start);
        expand_queue.push_back(ExpandState {
            node: start,
            parent: center,
            depth: 1,
            visited: v,
        });
    }

    while let Some(state) = expand_queue.pop_front() {
        if state.depth >= max_depth {
            continue;
        }
        let child_depth = state.depth + 1;

        // Phantom of parent: add if the bond used to reach this node was double.
        if let Some((_, bond_to_parent)) = mol.bond_between(state.node, state.parent) {
            if bond_to_parent.order == BondOrder::Double {
                let phantom_key = atom_key(mol, state.parent);
                layers.entry(child_depth).or_default().push(phantom_key);
            }
        }

        for (nb, _) in mol.neighbors(state.node) {
            if nb == state.parent || nb == center {
                continue;
            }
            let child_key = atom_key(mol, nb);
            let layer = layers.entry(child_depth).or_default();

            if state.visited.contains(&nb) {
                // Ring revisit: phantom only, no expansion.
                layer.push(child_key);
            } else {
                layer.push(child_key);
                let mut child_visited = state.visited.clone();
                child_visited.insert(nb);
                expand_queue.push_back(ExpandState {
                    node: nb,
                    parent: state.node,
                    depth: child_depth,
                    visited: child_visited,
                });
            }
        }
    }

    // Sort each layer descending and return as a Vec ordered by depth.
    let max_layer = layers.keys().copied().max().unwrap_or(0);
    let mut result = Vec::new();
    for d in 1..=max_layer {
        let mut layer = layers.remove(&d).unwrap_or_default();
        layer.sort_by(|a, b| cmp_key(*b, *a)); // descending
        result.push(layer);
    }
    result
}

/// Compare two branches from `center` starting at `a` and `b`.
///
/// Returns `Ordering::Greater` if branch `a` has higher CIP priority than `b`.
fn compare_branches(
    mol: &Molecule,
    center: AtomIdx,
    a: AtomIdx,
    b: AtomIdx,
) -> std::cmp::Ordering {
    use std::cmp::Ordering::*;

    // Depth-0 comparison: the substituent atoms themselves.
    let a_key = atom_key(mol, a);
    let b_key = atom_key(mol, b);
    match cmp_key(a_key, b_key) {
        Equal => {}
        other => return other,
    }

    // Sphere-by-sphere comparison.
    let a_spheres = cip_branch_spheres(mol, center, a);
    let b_spheres = cip_branch_spheres(mol, center, b);

    let max_depth = a_spheres.len().max(b_spheres.len());
    for d in 0..max_depth {
        let a_layer = a_spheres.get(d).map(|v| v.as_slice()).unwrap_or(&[]);
        let b_layer = b_spheres.get(d).map(|v| v.as_slice()).unwrap_or(&[]);

        let min_len = a_layer.len().min(b_layer.len());
        for i in 0..min_len {
            match cmp_key(a_layer[i], b_layer[i]) {
                Equal => {}
                other => return other,
            }
        }
        match a_layer.len().cmp(&b_layer.len()) {
            Equal => {}
            other => return other,
        }
    }

    Equal
}

/// Assign CIP priority ranks to `subs` (substituents of `center`).
///
/// Returns `None` if any two substituents have equal priority (tie).
/// Otherwise returns `Vec<u8>` of the same length, where `result[i]` is the
/// rank of `subs[i]` (1 = lowest CIP priority, N = highest).
fn rank_substituents(mol: &Molecule, center: AtomIdx, subs: &[AtomIdx]) -> Option<Vec<u8>> {
    let n = subs.len();
    if n == 0 {
        return Some(vec![]);
    }

    // Sort indices by CIP priority descending.
    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(|&i, &j| compare_branches(mol, center, subs[i], subs[j]).reverse());

    // Check for ties among adjacent elements after sorting.
    for k in 0..n - 1 {
        let i = indices[k];
        let j = indices[k + 1];
        if compare_branches(mol, center, subs[i], subs[j]) == std::cmp::Ordering::Equal {
            return None;
        }
    }

    // Assign ranks: indices[0] gets rank n (highest), indices[n-1] gets rank 1 (lowest).
    let mut ranks = vec![0u8; n];
    for (rank_from_top, &idx) in indices.iter().enumerate() {
        ranks[idx] = (n - rank_from_top) as u8;
    }

    Some(ranks)
}

fn assign_tetrahedral(mol: &Molecule, idx: AtomIdx) -> Option<CipCode> {
    let atom = mol.atom(idx);
    if atom.chirality == Chirality::None {
        return None;
    }

    // Collect neighbors in adjacency order.  In the SMILES parser, bonds are added
    // in the order the SMILES string is processed, so adjacency order matches SMILES
    // encounter order for the non-H substituents.
    let mut neighbors: Vec<AtomIdx> = mol.neighbors(idx).map(|(nb, _)| nb).collect();

    // For bracket atoms with explicit H (e.g. `[C@@H]`), the H occupies a specific
    // position in the SMILES chirality neighbor list:
    //
    //   • If the bracket atom has NO preceding atom in the SMILES chain (it is the
    //     first atom of the fragment, like `[C@@H](F)(Cl)Br`), the H is at position 0
    //     (the "from-viewer" slot) and the non-H neighbors follow.
    //
    //   • If the bracket atom HAS a preceding atom (like `N[C@@H](C)C(=O)O`), the
    //     preceding atom is at position 0, the H is at position 1, and the remaining
    //     non-H neighbors follow.
    //
    // "Preceding atom" = the atom that forms the bond into this atom from the left in
    // the SMILES string.  In the adjacency list, that atom is always added FIRST
    // (before branches and continuations) and therefore has a SMALLER atom index.
    let has_bracket_h = atom.hydrogen_count.map_or(false, |h| h > 0);
    if has_bracket_h {
        // Detect whether a preceding atom is present: the first neighbor, if its index
        // is smaller than `idx`, is the preceding atom.
        let has_preceding = neighbors
            .first()
            .map(|&nb| nb.0 < idx.0)
            .unwrap_or(false);
        let h_insert_pos = if has_preceding { 1 } else { 0 };
        neighbors.insert(h_insert_pos, AtomIdx(u32::MAX));
    }

    if neighbors.len() != 4 {
        return None;
    }

    // ranks[i] = CIP rank of neighbors[i]: 1 = lowest priority, 4 = highest.
    let ranks = rank_substituents(mol, idx, &neighbors)?;

    // --- Parity-based R/S determination -----------------------------------
    //
    // SMILES `@@` means: looking FROM neighbors[0], the sequence
    // neighbors[1]→neighbors[2]→neighbors[3] goes clockwise (CW).
    //
    // CIP R: looking FROM the rank-1 substituent, the sequence
    // rank2→rank3→rank4 (ascending priority) goes CW.
    //
    // Algorithm:
    // 1. Find where rank-1 is in the neighbors list (`lowest_pos`).
    // 2. Moving rank-1 to position 0 takes `lowest_pos` adjacent swaps,
    //    each one flipping CW↔CCW.  So the "effective_cw" (from rank-1's
    //    perspective) = smiles_cw XOR (lowest_pos is odd).
    // 3. After removing rank-1, the remaining three neighbors are in some order.
    //    Count how many swaps are needed to put them in ascending rank order
    //    [rank2, rank3, rank4].  An even number → same orientation; odd → flipped.
    // 4. is_r = effective_cw XOR (remaining_swaps is odd).

    let lowest_pos = ranks.iter().position(|&r| r == 1)?;
    let parity_odd = lowest_pos % 2 == 1;
    let smiles_cw = atom.chirality == Chirality::Clockwise;
    let cw_from_lowest = smiles_cw ^ parity_odd;

    // Remaining ranks in their current positional order (lowest_pos removed).
    let remaining_ranks: Vec<u8> = (0..4usize)
        .filter(|&i| i != lowest_pos)
        .map(|i| ranks[i])
        .collect();

    // Count swaps to reach the ascending-rank target [2, 3, 4].
    let remaining_swaps_odd = {
        let mut r = remaining_ranks.clone();
        let target = [2u8, 3, 4];
        let mut swaps = 0usize;
        for i in 0..3 {
            if r[i] != target[i] {
                if let Some(j_rel) = r[i + 1..].iter().position(|&x| x == target[i]) {
                    r.swap(i, j_rel + i + 1);
                    swaps += 1;
                } else {
                    return None; // invalid ranks (should not happen)
                }
            }
        }
        swaps % 2 == 1
    };

    // R if the effective CW sense matches the ascending-rank arrangement.
    let is_r = cw_from_lowest ^ remaining_swaps_odd;

    Some(if is_r { CipCode::R } else { CipCode::S })
}

/// Determine if a substituent is "up" relative to the alkene end it connects to.
///
/// Returns `Some(true)` = up, `Some(false)` = down, `None` = no stereo bond.
fn substituent_is_up(mol: &Molecule, alkene_end: AtomIdx, sub: AtomIdx) -> Option<bool> {
    let (_, bond) = mol.bond_between(alkene_end, sub)?;
    match bond.order {
        BondOrder::Up => {
            // `/` bond: atom1→atom2 goes "up"
            Some(bond.atom1 == alkene_end)
        }
        BondOrder::Down => {
            // `\` bond: atom1→atom2 goes "down"
            Some(bond.atom1 == sub)
        }
        _ => None,
    }
}

/// Assign E/Z for the double bond at `bond_idx`.
///
/// Returns `Some((atom_idx, E or Z))` using one of the double-bond endpoints
/// as the key atom index.  Returns `None` if the bond isn't double or stereo
/// cannot be determined.
fn assign_ez(mol: &Molecule, bond_idx: BondIdx) -> Option<(AtomIdx, CipCode)> {
    let bond = mol.bond(bond_idx);
    if bond.order != BondOrder::Double {
        return None;
    }

    let a1 = bond.atom1;
    let a2 = bond.atom2;

    // Non-double-bond neighbors for each alkene end (exclude the other alkene atom).
    let subs_a1: Vec<AtomIdx> = mol
        .neighbors(a1)
        .filter(|&(nb, bidx)| nb != a2 && mol.bond(bidx).order != BondOrder::Double)
        .map(|(nb, _)| nb)
        .collect();

    let subs_a2: Vec<AtomIdx> = mol
        .neighbors(a2)
        .filter(|&(nb, bidx)| nb != a1 && mol.bond(bidx).order != BondOrder::Double)
        .map(|(nb, _)| nb)
        .collect();

    if subs_a1.is_empty() || subs_a2.is_empty() {
        return None; // terminal alkene
    }

    // Highest-priority substituent with a stereo (Up/Down) bond at each end.
    let high_sub_a1 = highest_stereo_sub(mol, a1, &subs_a1)?;
    let high_sub_a2 = highest_stereo_sub(mol, a2, &subs_a2)?;

    let up_a1 = substituent_is_up(mol, a1, high_sub_a1)?;
    let up_a2 = substituent_is_up(mol, a2, high_sub_a2)?;

    // Same side → Z (zusammen); opposite → E (entgegen).
    let code = if up_a1 == up_a2 { CipCode::Z } else { CipCode::E };
    Some((a1, code))
}

/// From `subs` at `alkene_end`, return the highest CIP-priority substituent
/// that has an Up/Down bond to `alkene_end`.
fn highest_stereo_sub(mol: &Molecule, alkene_end: AtomIdx, subs: &[AtomIdx]) -> Option<AtomIdx> {
    let mut sorted: Vec<AtomIdx> = subs.to_vec();
    sorted.sort_by(|&a, &b| compare_branches(mol, alkene_end, a, b).reverse());
    sorted
        .into_iter()
        .find(|&sub| substituent_is_up(mol, alkene_end, sub).is_some())
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn cip_at(smiles: &str, atom_idx: usize) -> Option<CipCode> {
        let mol = parse(smiles).unwrap();
        let assignment = assign_cip(&mol);
        assignment.get(AtomIdx(atom_idx as u32))
    }

    // --- Tetrahedral R/S ---

    #[test]
    fn test_l_alanine_s() {
        // N[C@@H](C)C(=O)O — L-alanine, chiral center is atom 1
        assert_eq!(
            cip_at("N[C@@H](C)C(=O)O", 1),
            Some(CipCode::S),
            "L-alanine should be S"
        );
    }

    #[test]
    fn test_d_alanine_r() {
        // N[C@H](C)C(=O)O — D-alanine
        assert_eq!(
            cip_at("N[C@H](C)C(=O)O", 1),
            Some(CipCode::R),
            "D-alanine should be R"
        );
    }

    #[test]
    fn test_chfclbr_r() {
        // [C@@H](F)(Cl)Br — known R configuration
        // CIP priority: Br(35) > Cl(17) > F(9) > H(1)
        // @@H: looking from H, F→Cl→Br is CW → R
        assert_eq!(
            cip_at("[C@@H](F)(Cl)Br", 0),
            Some(CipCode::R),
            "[C@@H](F)(Cl)Br should be R"
        );
    }

    #[test]
    fn test_chfclbr_s() {
        // [C@H](F)(Cl)Br — S
        assert_eq!(
            cip_at("[C@H](F)(Cl)Br", 0),
            Some(CipCode::S),
            "[C@H](F)(Cl)Br should be S"
        );
    }

    #[test]
    fn test_no_chirality() {
        let mol = parse("CC(=O)O").unwrap();
        let assignment = assign_cip(&mol);
        let tetrahedral: Vec<_> = assignment
            .assignments
            .iter()
            .filter(|(_, c)| matches!(c, CipCode::R | CipCode::S))
            .collect();
        assert!(tetrahedral.is_empty(), "acetic acid has no chiral centers");
    }

    #[test]
    fn test_symmetric_center_none() {
        // No @/@@ annotation → no assignment attempted
        let mol = parse("CC(N)(N)CC").unwrap();
        let assignment = assign_cip(&mol);
        assert!(
            assignment.assignments.is_empty(),
            "no stereo annotation → no assignment"
        );
    }

    #[test]
    fn test_assignment_get() {
        let mol = parse("N[C@@H](C)C(=O)O").unwrap();
        let a = assign_cip(&mol);
        assert!(a.get(AtomIdx(1)).is_some(), "atom 1 should have CIP code");
        assert!(a.get(AtomIdx(0)).is_none(), "atom 0 (N) has no chirality");
    }

    #[test]
    fn test_r_lactic_acid_gives_answer() {
        // OC(=O)[C@@H](O)C — lactic acid; should give R or S
        let mol = parse("OC(=O)[C@@H](O)C").unwrap();
        let assignment = assign_cip(&mol);
        let chiral_idx = (0..mol.atom_count())
            .map(|i| AtomIdx(i as u32))
            .find(|&i| mol.atom(i).chirality != Chirality::None)
            .unwrap();
        let code = assignment.get(chiral_idx);
        assert!(
            code == Some(CipCode::R) || code == Some(CipCode::S),
            "should give R or S for lactic acid, got {:?}",
            code
        );
    }

    // --- E/Z double bonds ---

    #[test]
    fn test_trans_2_butene_e() {
        // C/C=C/C — trans-2-butene → E
        let mol = parse("C/C=C/C").unwrap();
        let assignment = assign_cip(&mol);
        let has_e = assignment.assignments.iter().any(|(_, c)| *c == CipCode::E);
        assert!(
            has_e,
            "Expected E for trans-2-butene, got {:?}",
            assignment.assignments
        );
    }

    #[test]
    fn test_cis_2_butene_z() {
        // C/C=C\C — cis-2-butene → Z
        let mol = parse("C/C=C\\C").unwrap();
        let assignment = assign_cip(&mol);
        let has_z = assignment.assignments.iter().any(|(_, c)| *c == CipCode::Z);
        assert!(
            has_z,
            "Expected Z for cis-2-butene, got {:?}",
            assignment.assignments
        );
    }

    #[test]
    fn test_fceccl_e() {
        // F/C=C/Cl → E (F and Cl on opposite sides)
        let mol = parse("F/C=C/Cl").unwrap();
        let assignment = assign_cip(&mol);
        let has_e = assignment.assignments.iter().any(|(_, c)| *c == CipCode::E);
        assert!(
            has_e,
            "Expected E for F/C=C/Cl, got {:?}",
            assignment.assignments
        );
    }

    #[test]
    fn test_fceccl_z() {
        // F/C=C\Cl → Z (F and Cl on same side)
        let mol = parse("F/C=C\\Cl").unwrap();
        let assignment = assign_cip(&mol);
        let has_z = assignment.assignments.iter().any(|(_, c)| *c == CipCode::Z);
        assert!(
            has_z,
            "Expected Z for F/C=C\\Cl, got {:?}",
            assignment.assignments
        );
    }

    #[test]
    fn test_no_ez_no_stereo_bond() {
        // C=C — no Up/Down bonds → no E/Z
        let mol = parse("C=C").unwrap();
        let assignment = assign_cip(&mol);
        let has_ez = assignment
            .assignments
            .iter()
            .any(|(_, c)| matches!(c, CipCode::E | CipCode::Z));
        assert!(!has_ez, "plain C=C should have no E/Z");
    }

    #[test]
    fn test_ez_terminal_no_crash() {
        // /C=C/F — terminal on one side; should not crash
        let mol = parse("/C=C/F").unwrap();
        let _ = assign_cip(&mol);
    }

    #[test]
    fn test_cip_assignment_methane() {
        let mol = parse("C").unwrap();
        let assignment = assign_cip(&mol);
        assert!(assignment.assignments.is_empty());
    }

    #[test]
    fn test_multiple_chiral_centers() {
        // Two chiral centers in a chain
        let mol = parse("N[C@@H](C)[C@H](C)N").unwrap();
        let assignment = assign_cip(&mol);
        let rs_count = assignment
            .assignments
            .iter()
            .filter(|(_, c)| matches!(c, CipCode::R | CipCode::S))
            .count();
        assert_eq!(rs_count, 2, "should assign 2 chiral centers");
    }

    #[test]
    fn test_cip_assignment_struct() {
        let mol = parse("N[C@@H](C)C(=O)O").unwrap();
        let a = assign_cip(&mol);
        assert!(!a.assignments.is_empty());
        let code = a.get(AtomIdx(1));
        assert_eq!(code, Some(CipCode::S));
    }

    #[test]
    fn test_r_s_are_consistent() {
        // The two SMILES should give opposite results
        let r_code = cip_at("[C@@H](F)(Cl)Br", 0);
        let s_code = cip_at("[C@H](F)(Cl)Br", 0);
        assert_ne!(r_code, s_code, "@ and @@ must give opposite results");
    }

    #[test]
    fn test_e_z_are_consistent() {
        // /C=C/C (trans) and /C=C\C (cis) must differ
        let mol_e = parse("C/C=C/C").unwrap();
        let mol_z = parse("C/C=C\\C").unwrap();
        let assign_e = assign_cip(&mol_e);
        let assign_z = assign_cip(&mol_z);
        let has_e = assign_e.assignments.iter().any(|(_, c)| *c == CipCode::E);
        let has_z = assign_z.assignments.iter().any(|(_, c)| *c == CipCode::Z);
        assert!(has_e && has_z, "trans must be E, cis must be Z");
    }
}
