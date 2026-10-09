//! `Canon::canonicalizeFragment` (RDKit 2026.03.1 `Canon.cpp`): the
//! canonical DFS, chirality tags relative to the output order and the
//! `/` `\` directions written around stereo double bonds.

use super::RdkitSmilesError;
use super::mol::{BondDir, BondStereo, BondType, ChiralTag, Mol, insert_implicit_nbors};
use super::stereo::is_atom_potential_tetrahedral_center;

const MAX_NATOMS: i64 = 5000;
const MAX_CYCLES: usize = 1024;
const MAX_BONDTYPE: i64 = 32;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum StackElem {
    Atom(usize),
    /// A bond and the atom written to its left.
    Bond(usize, usize),
    Ring(u32),
    BranchOpen,
    BranchClose,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Color {
    White,
    Grey,
    Black,
}

struct Traversal<'a> {
    mol: &'a Mol,
    ranks: &'a [u32],
    ring_closures: Vec<Vec<usize>>,
    traversal_bond_order: Vec<Vec<usize>>,
    cycles_available: Vec<bool>,
    /// `_TraversalRingClosureBond`.
    ring_closure_idx: Vec<Option<u32>>,
    stack: Vec<StackElem>,
    /// The sorted destinations of every pending frame, back to back.
    possibles: Vec<Possible>,
}

type Possible = (i32, usize, usize);

/// A pending `dfsFindCycles` call: its atom, the bond it came in on, its
/// sorted destinations and the next one to look at.
struct CycleFrame {
    atom: usize,
    /// Its destinations are `possibles[start..end]`.
    start: usize,
    end: usize,
    next: usize,
}

/// A pending `dfsBuildStack` call.
struct StackFrame {
    atom: usize,
    /// Its destinations are `possibles[start..end]`.
    start: usize,
    end: usize,
    next: usize,
    /// The bonds in traversal order, kept for chiral atoms only (the only
    /// ones whose order is read).
    trav: Vec<usize>,
    keep_trav: bool,
    /// A branch was opened for the child being visited.
    close_branch: bool,
}

impl Traversal<'_> {
    fn bond_type_value(&self, b: usize) -> i64 {
        self.mol.bonds[b].bt as i64
    }

    /// Appends the destinations `dfsFindCycles` sorts at `atom` to
    /// `possibles`; returns their range.
    fn cycle_possibles(
        &mut self,
        atom: usize,
        in_bond: Option<usize>,
        colors: &[Color],
    ) -> (usize, usize) {
        let start = self.possibles.len();
        for &b in &self.mol.atom_bonds[atom] {
            if Some(b) == in_bond {
                continue;
            }
            let other = self.mol.bonds[b].other(atom);
            let mut rank = i64::from(self.ranks[other]);
            if colors[other] == Color::Grey {
                rank -= (MAX_BONDTYPE + 1) * MAX_NATOMS * MAX_NATOMS;
                rank += (MAX_BONDTYPE - self.bond_type_value(b)) * MAX_NATOMS;
            } else if self.mol.num_bond_rings(b) != 0 {
                rank += (MAX_BONDTYPE - self.bond_type_value(b)) * MAX_NATOMS * MAX_NATOMS;
            }
            // RDKit computes in `unsigned int` and stores an `int`.
            self.possibles.push((rank as i32, other, b));
        }
        self.possibles[start..].sort_by_key(|p| p.0);
        (start, self.possibles.len())
    }

    /// `dfsFindCycles`, with an explicit stack instead of recursion.
    fn find_cycles(&mut self, start: usize, colors: &mut [Color]) {
        colors[start] = Color::Grey;
        let (s0, e0) = self.cycle_possibles(start, None, colors);
        let mut frames = vec![CycleFrame {
            atom: start,
            start: s0,
            end: e0,
            next: s0,
        }];
        while let Some(frame) = frames.last_mut() {
            if frame.next == frame.end {
                colors[frame.atom] = Color::Black;
                self.possibles.truncate(frame.start);
                frames.pop();
                continue;
            }
            let (_, other, b) = self.possibles[frame.next];
            let atom = frame.atom;
            frame.next += 1;
            match colors[other] {
                Color::White => {
                    colors[other] = Color::Grey;
                    let (s, e) = self.cycle_possibles(other, Some(b), colors);
                    frames.push(CycleFrame {
                        atom: other,
                        start: s,
                        end: e,
                        next: s,
                    });
                }
                Color::Grey => {
                    self.ring_closures[other].push(b);
                    self.ring_closures[atom].push(b);
                }
                Color::Black => {}
            }
        }
    }

    /// The opening part of a `dfsBuildStack` call: the atom, its ring
    /// closures and its sorted destinations.
    fn enter_atom(
        &mut self,
        atom: usize,
        in_bond: Option<usize>,
        colors: &mut [Color],
    ) -> Result<StackFrame, RdkitSmilesError> {
        let mol = self.mol;
        self.stack.push(StackElem::Atom(atom));
        colors[atom] = Color::Grey;
        let keep_trav = mol.atoms[atom].chiral != ChiralTag::Unspecified;
        let mut trav: Vec<usize> = if keep_trav {
            Vec::with_capacity(mol.degree(atom))
        } else {
            Vec::new()
        };
        if let Some(b) = in_bond
            && keep_trav
        {
            trav.push(b);
        }
        let closures = std::mem::take(&mut self.ring_closures[atom]);
        if !closures.is_empty() {
            let mut rings_closed = Vec::new();
            for &b in &closures {
                if keep_trav {
                    trav.push(b);
                }
                if let Some(ring_idx) = self.ring_closure_idx[b] {
                    self.stack.push(StackElem::Bond(b, atom));
                    self.stack.push(StackElem::Ring(ring_idx));
                    rings_closed.push(ring_idx as usize - 1);
                } else {
                    let lowest =
                        self.cycles_available
                            .iter()
                            .position(|&x| x)
                            .ok_or_else(|| {
                                RdkitSmilesError::Unsupported(
                                    "Too many rings open at once. SMILES cannot be generated."
                                        .into(),
                                )
                            })?;
                    self.cycles_available[lowest] = false;
                    let ring_idx = lowest as u32 + 1;
                    self.ring_closure_idx[b] = Some(ring_idx);
                    self.stack.push(StackElem::Ring(ring_idx));
                }
            }
            for r in rings_closed {
                self.cycles_available[r] = true;
            }
        }
        let start = self.possibles.len();
        for &b in &mol.atom_bonds[atom] {
            if Some(b) == in_bond {
                continue;
            }
            let other = mol.bonds[b].other(atom);
            // `seenFromHere`: the atom and its ring-closure partners.
            if colors[other] != Color::White
                || other == atom
                || closures.iter().any(|&c| mol.bonds[c].other(atom) == other)
            {
                continue;
            }
            let mut rank = i64::from(self.ranks[other]);
            if mol.num_bond_rings(b) != 0 {
                rank += (MAX_BONDTYPE - self.bond_type_value(b)) * MAX_NATOMS * MAX_NATOMS;
            }
            self.possibles.push((rank as i32, other, b));
        }
        self.possibles[start..].sort_by_key(|p| p.0);
        self.ring_closures[atom] = closures;
        Ok(StackFrame {
            atom,
            start,
            end: self.possibles.len(),
            next: start,
            trav,
            keep_trav,
            close_branch: false,
        })
    }

    /// `dfsBuildStack`, with an explicit stack instead of recursion.
    fn build_stack(&mut self, start: usize, colors: &mut [Color]) -> Result<(), RdkitSmilesError> {
        let first = self.enter_atom(start, None, colors)?;
        let mut frames = vec![first];
        while let Some(frame) = frames.last_mut() {
            if frame.close_branch {
                self.stack.push(StackElem::BranchClose);
                frame.close_branch = false;
            }
            let np = frame.end;
            let mut child = None;
            while frame.next < np {
                let k = frame.next;
                let (_, other, b) = self.possibles[k];
                frame.next += 1;
                if colors[other] != Color::White {
                    continue;
                }
                if frame.keep_trav {
                    frame.trav.push(b);
                }
                if k + 1 != np {
                    self.stack.push(StackElem::BranchOpen);
                    frame.close_branch = true;
                }
                self.stack.push(StackElem::Bond(b, frame.atom));
                child = Some((other, b));
                break;
            }
            match child {
                Some((other, b)) => {
                    let f = self.enter_atom(other, Some(b), colors)?;
                    frames.push(f);
                }
                None => {
                    let done = frames.pop().expect("frame");
                    self.possibles.truncate(done.start);
                    self.traversal_bond_order[done.atom] = done.trav;
                    colors[done.atom] = Color::Black;
                }
            }
        }
        Ok(())
    }
}

/// Result of [`canonicalize_fragment`].
pub(crate) struct Canonicalized {
    pub stack: Vec<StackElem>,
}

/// `chiralAtomNeedsTagInversion` with the property cache in place.
fn chiral_atom_needs_tag_inversion(
    mol: &Mol,
    a: usize,
    is_first: bool,
    num_closures: usize,
) -> bool {
    let atom = &mol.atoms[a];
    let fourth_valence = atom.num_explicit_hs == 1
        || (!mol.needs_update_property_cache(a) && mol.implicit_valence(a) == 1);
    let unsaturated = mol.atom_bonds[a]
        .iter()
        .any(|&b| mol.bonds[b].bt.as_double() > 1.0);
    mol.degree(a) == 3
        && ((is_first && atom.num_explicit_hs == 1)
            || (!fourth_valence && num_closures == 1 && !unsaturated))
}

/// `Canon::canonicalizeFragment(mol, atomIdx, colors, ranks, molStack, ...,
/// doIsomericSmiles=true, doRandom=false, doChiralInversions=true)` for a
/// whole (connected) molecule. Updates chiral tags and bond directions in
/// `mol` as RDKit does.
pub(crate) fn canonicalize_fragment(
    mol: &mut Mol,
    start: usize,
    ranks: &[u32],
    isomeric: bool,
) -> Result<Canonicalized, RdkitSmilesError> {
    let n = mol.atoms.len();
    let nb = mol.bonds.len();
    let (stack, ring_closures, traversal_bond_order, ring_closure_idx) = {
        let mut t = Traversal {
            mol,
            ranks,
            ring_closures: vec![Vec::new(); n],
            traversal_bond_order: vec![Vec::new(); n],
            cycles_available: vec![true; MAX_CYCLES],
            ring_closure_idx: vec![None; nb],
            stack: Vec::with_capacity(4 * n),
            possibles: Vec::with_capacity(2 * nb + 2),
        };
        let mut tcolors = vec![Color::White; n];
        t.find_cycles(start, &mut tcolors);
        let mut colors = vec![Color::White; n];
        t.build_stack(start, &mut colors)?;
        (
            t.stack,
            t.ring_closures,
            t.traversal_bond_order,
            t.ring_closure_idx,
        )
    };

    // Traversal-order parity of every chiral atom.
    let first_idx = match stack.first() {
        Some(StackElem::Atom(a)) => *a,
        _ => return Err(RdkitSmilesError::Unsupported("empty traversal".into())),
    };
    let mut num_swaps_odd = vec![false; n];
    let mut permutation = vec![0u32; n];
    for a in 0..n {
        if !isomeric || mol.atoms[a].chiral == ChiralTag::Unspecified {
            continue;
        }
        let nontet = mol.atoms[a].chiral.nontet();
        if !is_atom_potential_tetrahedral_center(mol, a) && nontet.is_none() {
            continue;
        }
        let perm = if nontet.is_some() {
            mol.atoms[a].chiral_perm
        } else {
            0
        };
        let first_in_part = a == first_idx;
        let mut order = traversal_bond_order[a].clone();
        if order.len() < mol.degree(a) {
            for &b in &mol.atom_bonds[a] {
                if !order.contains(&b) {
                    order.push(b);
                }
            }
        }
        if perm != 0 {
            // The permutation relative to the output order (with implicit
            // ligands where SMILES puts them).
            let mut probe: Vec<Option<usize>> = order.iter().map(|&b| Some(b)).collect();
            insert_implicit_nbors(&mut probe, nontet.expect("non-tetrahedral"), first_in_part);
            permutation[a] = mol.chiral_permutation(a, &probe, false);
            continue;
        }
        let mut odd = mol.perturbation_is_odd(a, &order);
        if chiral_atom_needs_tag_inversion(mol, a, first_in_part, ring_closures[a].len()) {
            odd = !odd;
        }
        num_swaps_odd[a] = odd;
    }

    let mut atom_visit = vec![0usize; n];
    let mut bond_visit = vec![0usize; nb];
    for (pos, e) in stack.iter().enumerate() {
        match *e {
            StackElem::Atom(a) => atom_visit[a] = pos,
            StackElem::Bond(b, _) => {
                bond_visit[b] = pos;
                mol.bonds[b].dir = BondDir::None;
            }
            _ => {}
        }
    }

    let mut st = DirState {
        bond_dir_counts: vec![0i8; nb],
        atom_dir_counts: vec![0i8; n],
        bond_visit,
        atom_visit,
        ring_closure: ring_closure_idx.iter().map(|x| x.is_some()).collect(),
    };
    canonicalize_double_bonds(mol, &mut st, &stack)?;

    // Chirality relative to the output order (incl. ring stereo).
    let mut ring_adjusted = vec![false; n];
    for e in &stack {
        let StackElem::Atom(a) = *e else { continue };
        if !isomeric || mol.atoms[a].chiral == ChiralTag::Unspecified {
            continue;
        }
        if let Some(rsa) = mol.atoms[a].ring_stereo_atoms.clone() {
            if !ring_adjusted[a] {
                mol.atoms[a].chiral = ChiralTag::Ccw;
                ring_adjusted[a] = true;
            }
            for nbr_v in rsa {
                let nbr = (nbr_v.unsigned_abs() - 1) as usize;
                if !ring_adjusted[nbr] && st.atom_visit[nbr] > st.atom_visit[a] {
                    mol.atoms[nbr].chiral = mol.atoms[a].chiral;
                    if nbr_v < 0 {
                        mol.atoms[nbr].invert_chirality();
                    }
                    if num_swaps_odd[a] {
                        if !num_swaps_odd[nbr] {
                            mol.atoms[nbr].invert_chirality();
                        }
                    } else if num_swaps_odd[nbr] {
                        mol.atoms[nbr].invert_chirality();
                    }
                    ring_adjusted[nbr] = true;
                }
            }
        } else if mol.atoms[a].chiral.is_tetrahedral() {
            if num_swaps_odd[a] {
                mol.atoms[a].invert_chirality();
            }
        } else if permutation[a] != 0 {
            mol.atoms[a].chiral_perm = permutation[a];
        }
    }

    remove_unwanted_bond_dir_specs(mol, &mut st, &stack);
    remove_redundant_bond_dir_specs(mol, &mut st, &stack);
    Ok(Canonicalized { stack })
}

struct DirState {
    bond_dir_counts: Vec<i8>,
    atom_dir_counts: Vec<i8>,
    bond_visit: Vec<usize>,
    atom_visit: Vec<usize>,
    /// Bonds carrying `_TraversalRingClosureBond`.
    ring_closure: Vec<bool>,
}

/// A bond around a double bond with its "flipped" status.
#[derive(Clone, Copy)]
struct Side {
    bond: usize,
    flipped: bool,
}

fn set_direction_from_neighboring_bond(mol: &mut Mol, source: Side, target: Side) {
    let mut dir = mol.bonds[source.bond].dir;
    if source.flipped == target.flipped {
        dir = dir.flipped();
    }
    mol.bonds[target.bond].dir = dir;
}

fn get_reference_direction(
    mol: &Mol,
    dbl: usize,
    ref_atom: usize,
    target_atom: usize,
    ref_ctrl: Side,
    target: Side,
) -> Result<BondDir, RdkitSmilesError> {
    let dbond = &mol.bonds[dbl];
    let mut dir = match dbond.stereo {
        BondStereo::E => mol.bonds[ref_ctrl.bond].dir,
        BondStereo::Z => mol.bonds[ref_ctrl.bond].dir.flipped(),
        _ => BondDir::None,
    };
    if dir == BondDir::None {
        return Err(RdkitSmilesError::Unsupported("stereo not set".into()));
    }
    let stereo_atoms = &dbond.stereo_atoms;
    if mol.degree(ref_atom) == 3
        && !stereo_atoms.contains(&mol.bonds[ref_ctrl.bond].other(ref_atom))
    {
        dir = dir.flipped();
    }
    if mol.degree(target_atom) == 3
        && !stereo_atoms.contains(&mol.bonds[target.bond].other(target_atom))
    {
        dir = dir.flipped();
    }
    if ref_ctrl.flipped != target.flipped {
        dir = dir.flipped();
    }
    Ok(dir)
}

/// `fixConflictAcrossDoubleBond`. RDKit iterates over copies of the two
/// bonds, so its "other bond" is always the first one.
#[allow(clippy::too_many_arguments)]
fn fix_conflict_across_double_bond(
    mol: &Mol,
    st: &mut DirState,
    dbl: usize,
    atom: usize,
    first: Side,
    second: Side,
    ref_atom: usize,
    ref_bond: Side,
) -> Result<bool, RdkitSmilesError> {
    for side in [first, second] {
        let other_bond = first.bond;
        let other_idx = mol.bonds[other_bond].other(atom);
        if st.atom_dir_counts[other_idx] != 2 {
            continue;
        }
        let expected = get_reference_direction(mol, dbl, ref_atom, atom, ref_bond, side)?;
        if expected == mol.bonds[side.bond].dir {
            st.bond_dir_counts[other_bond] = 0;
            st.atom_dir_counts[atom] -= 1;
            st.atom_dir_counts[other_idx] -= 1;
            return Ok(true);
        }
    }
    Ok(false)
}

#[allow(clippy::too_many_arguments)]
fn handle_dir_conflicts_across_double_bond(
    mol: &Mol,
    st: &mut DirState,
    dbl: usize,
    atom1: usize,
    atom1_consistent: bool,
    first1: Side,
    second1: Option<Side>,
    atom2: usize,
    atom2_consistent: bool,
    first2: Side,
    second2: Option<Side>,
) -> Result<bool, RdkitSmilesError> {
    if atom1_consistent && atom2_consistent {
        let expected = get_reference_direction(mol, dbl, atom1, atom2, first1, first2)?;
        Ok(expected == mol.bonds[first2.bond].dir)
    } else if !atom2_consistent && atom1_consistent {
        let second2 = second2.expect("inconsistent side has two bonds");
        fix_conflict_across_double_bond(mol, st, dbl, atom2, first2, second2, atom1, first1)
    } else if !atom1_consistent && atom2_consistent {
        let second1 = second1.expect("inconsistent side has two bonds");
        fix_conflict_across_double_bond(mol, st, dbl, atom1, first1, second1, atom2, first2)
    } else {
        let second1 = second1.expect("inconsistent side has two bonds");
        let second2 = second2.expect("inconsistent side has two bonds");
        for a1b in [first1, second1] {
            for a2b in [first2, second2] {
                let expected = get_reference_direction(mol, dbl, atom1, atom2, a1b, a2b)?;
                if expected == mol.bonds[a2b.bond].dir {
                    // Copies again: the "other" bonds are the first ones.
                    let a1_other = first1.bond;
                    let a1_other_idx = mol.bonds[a1_other].other(atom1);
                    if st.atom_dir_counts[a1_other_idx] != 2 {
                        continue;
                    }
                    let a2_other = first2.bond;
                    let a2_other_idx = mol.bonds[a2_other].other(atom2);
                    if a1_other_idx == a2_other_idx {
                        continue;
                    }
                    if st.atom_dir_counts[a2_other_idx] != 2 {
                        continue;
                    }
                    st.bond_dir_counts[a1_other] = 0;
                    st.atom_dir_counts[atom1] -= 1;
                    st.atom_dir_counts[a1_other_idx] -= 1;
                    st.bond_dir_counts[a2_other] = 0;
                    st.atom_dir_counts[atom2] -= 1;
                    st.atom_dir_counts[a2_other_idx] -= 1;
                    return Ok(true);
                }
            }
        }
        Ok(false)
    }
}

fn same_side_dirs_are_compatible(mol: &Mol, first: Side, second: Side) -> bool {
    let should_match = first.flipped != second.flipped;
    let dirs_match = mol.bonds[first.bond].dir == mol.bonds[second.bond].dir;
    dirs_match == should_match
}

/// `Canon::canonicalizeDoubleBond`.
fn canonicalize_double_bond(
    mol: &mut Mol,
    st: &mut DirState,
    dbl: usize,
) -> Result<(), RdkitSmilesError> {
    let mut atom1 = mol.bonds[dbl].begin;
    let mut atom2 = mol.bonds[dbl].end;
    if !matches!(mol.degree(atom1), 2 | 3) || !matches!(mol.degree(atom2), 2 | 3) {
        return Ok(());
    }
    if st.atom_visit[atom1] >= st.atom_visit[atom2] {
        std::mem::swap(&mut atom1, &mut atom2);
    }
    let num_bonds = mol.bonds.len();
    let find = |mol: &Mol, st: &DirState, atom: usize| {
        let mut first: Option<usize> = None;
        let mut second: Option<usize> = None;
        let mut dir_set = false;
        let mut first_visit = num_bonds + 1;
        for &b in &mol.atom_bonds[atom] {
            if b == dbl || !mol.bonds[b].can_set_double_bond_stereo() {
                continue;
            }
            if st.bond_dir_counts[b] > 0 {
                dir_set = true;
            }
            if first.is_none() || st.bond_visit[b] < first_visit {
                if first.is_some() {
                    second = first;
                }
                first = Some(b);
                first_visit = st.bond_visit[b];
            } else {
                second = Some(b);
            }
        }
        (first, second, dir_set)
    };
    let (f1, s1, dir1_set) = find(mol, st, atom1);
    let (f2, s2, dir2_set) = find(mol, st, atom2);
    let (Some(f1), Some(f2)) = (f1, f2) else {
        return Ok(());
    };
    // "Flipped": the anchor is written after atom1, or before atom2.
    let flipped1 = |b: usize| {
        let anchor = mol.bonds[b].other(atom1);
        (st.atom_visit[atom1] < st.atom_visit[anchor]) != st.ring_closure[b]
    };
    let flipped2 = |b: usize| {
        let anchor = mol.bonds[b].other(atom2);
        (st.atom_visit[anchor] < st.atom_visit[atom2]) != st.ring_closure[b]
    };
    let first1 = Side {
        bond: f1,
        flipped: flipped1(f1),
    };
    let second1 = s1.map(|b| Side {
        bond: b,
        flipped: flipped1(b),
    });
    let first2 = Side {
        bond: f2,
        flipped: flipped2(f2),
    };
    let second2 = s2.map(|b| Side {
        bond: b,
        flipped: flipped2(b),
    });

    if dir1_set && dir2_set {
        let mut atom1_consistent = true;
        if let Some(sec) = second1 {
            if st.bond_dir_counts[first1.bond] == 0 {
                set_direction_from_neighboring_bond(mol, sec, first1);
            } else if st.bond_dir_counts[sec.bond] == 0 {
                set_direction_from_neighboring_bond(mol, first1, sec);
            } else {
                atom1_consistent = same_side_dirs_are_compatible(mol, first1, sec);
            }
            st.bond_dir_counts[sec.bond] += 1;
            st.atom_dir_counts[atom1] += 1;
        }
        st.bond_dir_counts[first1.bond] += 1;
        st.atom_dir_counts[atom1] += 1;

        let mut atom2_consistent = true;
        if let Some(sec) = second2 {
            if st.bond_dir_counts[first2.bond] == 0 {
                set_direction_from_neighboring_bond(mol, sec, first2);
            } else if st.bond_dir_counts[sec.bond] == 0 {
                set_direction_from_neighboring_bond(mol, first2, sec);
            } else {
                atom2_consistent = same_side_dirs_are_compatible(mol, first2, sec);
            }
            st.bond_dir_counts[sec.bond] += 1;
            st.atom_dir_counts[atom2] += 1;
        }
        st.bond_dir_counts[first2.bond] += 1;
        st.atom_dir_counts[atom2] += 1;

        // A failure only logs a warning in RDKit.
        let _ = handle_dir_conflicts_across_double_bond(
            mol,
            st,
            dbl,
            atom1,
            atom1_consistent,
            first1,
            second1,
            atom2,
            atom2_consistent,
            first2,
            second2,
        )?;
        return Ok(());
    }

    let mut set_from_bond1 = true;
    let mut atom1_ctrl = first1;
    let mut atom2_ctrl = first2;
    if !dir1_set && !dir2_set {
        mol.bonds[first1.bond].dir = BondDir::EndUpRight;
        st.bond_dir_counts[first1.bond] += 1;
        st.atom_dir_counts[atom1] += 1;
    } else if !dir2_set {
        if st.bond_dir_counts[first1.bond] > 0 {
            st.bond_dir_counts[first1.bond] += 1;
            st.atom_dir_counts[atom1] += 1;
            if let Some(sec) = second1
                && st.bond_dir_counts[sec.bond] != 0
            {
                // Incompatible directions only log a warning in RDKit.
                st.bond_dir_counts[sec.bond] += 1;
                st.atom_dir_counts[atom1] += 1;
            }
        } else {
            let sec = second1.ok_or_else(|| {
                RdkitSmilesError::Unsupported("inconsistent double bond state".into())
            })?;
            set_direction_from_neighboring_bond(mol, sec, first1);
            st.bond_dir_counts[sec.bond] += 1;
            st.bond_dir_counts[first1.bond] += 1;
            st.atom_dir_counts[atom1] += 2;
            atom1_ctrl = sec;
        }
    } else {
        set_from_bond1 = false;
        if st.bond_dir_counts[first2.bond] > 0 {
            st.bond_dir_counts[first2.bond] += 1;
            st.atom_dir_counts[atom2] += 1;
            if let Some(sec) = second2
                && st.bond_dir_counts[sec.bond] != 0
            {
                st.bond_dir_counts[sec.bond] += 1;
                st.atom_dir_counts[atom2] += 1;
            }
        } else {
            let sec = second2.ok_or_else(|| {
                RdkitSmilesError::Unsupported("inconsistent double bond state".into())
            })?;
            set_direction_from_neighboring_bond(mol, sec, first2);
            st.bond_dir_counts[sec.bond] += 1;
            st.bond_dir_counts[first2.bond] += 1;
            st.atom_dir_counts[atom2] += 2;
            atom2_ctrl = sec;
        }
    }

    if set_from_bond1 {
        let dir = get_reference_direction(mol, dbl, atom1, atom2, atom1_ctrl, first2)?;
        mol.bonds[first2.bond].dir = dir;
        st.bond_dir_counts[first2.bond] += 1;
        st.atom_dir_counts[atom2] += 1;
    } else {
        let dir = get_reference_direction(mol, dbl, atom2, atom1, atom2_ctrl, first1)?;
        mol.bonds[first1.bond].dir = dir;
        st.bond_dir_counts[first1.bond] += 1;
        st.atom_dir_counts[atom1] += 1;
    }

    if mol.degree(atom1) == 3
        && let Some(sec) = second1
        && st.bond_dir_counts[sec.bond] == 0
    {
        set_direction_from_neighboring_bond(mol, first1, sec);
        st.bond_dir_counts[sec.bond] += 1;
        st.atom_dir_counts[atom1] += 1;
    }
    if mol.degree(atom2) == 3
        && let Some(sec) = second2
        && st.bond_dir_counts[sec.bond] == 0
    {
        set_direction_from_neighboring_bond(mol, first2, sec);
        st.bond_dir_counts[sec.bond] += 1;
        st.atom_dir_counts[atom2] += 1;
    }
    Ok(())
}

fn is_stereo_double(mol: &Mol, b: usize) -> bool {
    mol.bonds[b].bt == BondType::Double && mol.bonds[b].stereo > BondStereo::Any
}

/// `Canon::canonicalizeDoubleBonds`.
fn canonicalize_double_bonds(
    mol: &mut Mol,
    st: &mut DirState,
    stack: &[StackElem],
) -> Result<(), RdkitSmilesError> {
    let neighboring_stereo_bond = |mol: &Mol, dbl_atom: usize, nbr_bond: usize| -> Option<usize> {
        let other = mol.bonds[nbr_bond].other(dbl_atom);
        mol.atom_bonds[other]
            .iter()
            .copied()
            .find(|&b| b != nbr_bond && is_stereo_double(mol, b))
    };
    let mut stereo_nbrs: Vec<Vec<usize>> = vec![Vec::new(); mol.bonds.len()];
    let mut queue: Vec<usize> = Vec::new();
    for e in stack {
        let StackElem::Bond(b, _) = *e else { continue };
        mol.bonds[b].dir = BondDir::None;
        if mol.bonds[b].bt != BondType::Double
            || mol.bonds[b].stereo <= BondStereo::Any
            || mol.bonds[b].stereo_atoms.len() < 2
        {
            mol.bonds[b].stereo = BondStereo::None;
            continue;
        }
        let mut current = Vec::new();
        for dbl_atom in [mol.bonds[b].begin, mol.bonds[b].end] {
            for &nbr_bond in &mol.atom_bonds[dbl_atom] {
                if !mol.bonds[nbr_bond].can_have_direction() {
                    continue;
                }
                if let Some(nd) = neighboring_stereo_bond(mol, dbl_atom, nbr_bond) {
                    current.push(nd);
                }
            }
        }
        current.sort_by_key(|&x| st.bond_visit[x]);
        stereo_nbrs[b] = current;
        queue.push(b);
    }
    // std::priority_queue: most stereo neighbours first, then the lowest
    // position in the stack.
    queue.sort_by(|&x, &y| {
        stereo_nbrs[y]
            .len()
            .cmp(&stereo_nbrs[x].len())
            .then(st.bond_visit[x].cmp(&st.bond_visit[y]))
    });
    let mut seen = vec![false; mol.bonds.len()];
    for b in queue {
        if seen[b] {
            continue;
        }
        let mut connected = std::collections::VecDeque::from([b]);
        while let Some(cur) = connected.pop_front() {
            if seen[cur] {
                continue;
            }
            canonicalize_double_bond(mol, st, cur)?;
            seen[cur] = true;
            for &nbr in &stereo_nbrs[cur] {
                if !seen[nbr] {
                    connected.push_back(nbr);
                }
            }
        }
    }
    Ok(())
}

/// `Canon::removeUnwantedBondDirSpecs`.
fn remove_unwanted_bond_dir_specs(mol: &mut Mol, st: &mut DirState, stack: &[StackElem]) {
    for e in stack {
        let StackElem::Bond(b, _) = *e else { continue };
        if mol.bonds[b].bt != BondType::Double || mol.bonds[b].stereo > BondStereo::Any {
            continue;
        }
        let first_atom = mol.bonds[b].begin;
        let second_atom = mol.bonds[b].end;
        if mol.degree(first_atom) == 1 || mol.degree(second_atom) == 1 {
            continue;
        }
        let mut candidates: Vec<usize> = mol.atom_bonds[first_atom]
            .iter()
            .copied()
            .filter(|&x| st.bond_dir_counts[x] != 0)
            .collect();
        if candidates.is_empty() {
            continue;
        }
        if st.atom_dir_counts[first_atom] != 0 {
            candidates.clear();
        }
        let mut on_second = 0;
        for &x in &mol.atom_bonds[second_atom] {
            if st.bond_dir_counts[x] != 0 {
                candidates.push(x);
                on_second += 1;
            }
        }
        if on_second == 0 {
            continue;
        }
        if st.atom_dir_counts[second_atom] != 0 {
            continue;
        }
        candidates.sort_by_key(|&x| st.bond_visit[x]);
        for cand in candidates {
            let cb = &mol.bonds[cand];
            let other = if cb.begin == first_atom || cb.end == first_atom {
                cb.other(first_atom)
            } else {
                cb.other(second_atom)
            };
            if st.atom_dir_counts[other] == 2 {
                st.bond_dir_counts[cand] = 0;
                mol.bonds[cand].dir = BondDir::None;
                st.atom_dir_counts[other] -= 1;
                break;
            }
        }
    }
}

/// `Canon::clearBondDirs`.
fn clear_bond_dirs(mol: &mut Mol, st: &mut DirState, ref_bond: usize, from_atom: usize) {
    let clear_direction = |mol: &mut Mol, st: &mut DirState, bond: usize| {
        st.bond_dir_counts[bond] -= 1;
        if st.bond_dir_counts[bond] == 0 {
            mol.bonds[bond].dir = BondDir::None;
            st.atom_dir_counts[from_atom] -= 1;
            let other = mol.bonds[bond].other(from_atom);
            if st.atom_dir_counts[other] != 0 {
                st.atom_dir_counts[other] -= 1;
            }
        }
    };
    for o_bond in mol.atom_bonds[from_atom].clone() {
        if o_bond != ref_bond && mol.bonds[o_bond].can_have_direction() {
            let (ob, oe) = (mol.bonds[o_bond].begin, mol.bonds[o_bond].end);
            let (rb, re) = (mol.bonds[ref_bond].begin, mol.bonds[ref_bond].end);
            if st.bond_dir_counts[o_bond] >= st.bond_dir_counts[ref_bond]
                && st.atom_dir_counts[ob] != 1
                && st.atom_dir_counts[oe] != 1
            {
                clear_direction(mol, st, o_bond);
            } else if st.atom_dir_counts[rb] != 1 && st.atom_dir_counts[re] != 1 {
                clear_direction(mol, st, ref_bond);
            }
            break;
        }
    }
}

/// `Canon::removeRedundantBondDirSpecs`.
fn remove_redundant_bond_dir_specs(mol: &mut Mol, st: &mut DirState, stack: &[StackElem]) {
    let clear_from_atom = |mol: &mut Mol, st: &mut DirState, t_bond: usize, atom: usize| {
        if st.atom_dir_counts[atom] < 2 {
            return;
        }
        let has_stereo_double = mol.atom_bonds[atom]
            .iter()
            .any(|&b| b != t_bond && is_stereo_double(mol, b));
        if has_stereo_double {
            clear_bond_dirs(mol, st, t_bond, atom);
        }
    };
    for e in stack {
        let StackElem::Bond(t_bond, left) = *e else {
            continue;
        };
        let right = mol.bonds[t_bond].other(left);
        if mol.bonds[t_bond].can_have_direction() && st.bond_dir_counts[t_bond] != 0 {
            clear_from_atom(mol, st, t_bond, left);
            clear_from_atom(mol, st, t_bond, right);
        } else if mol.bonds[t_bond].dir != BondDir::None {
            mol.bonds[t_bond].dir = BondDir::None;
        }
    }
}
