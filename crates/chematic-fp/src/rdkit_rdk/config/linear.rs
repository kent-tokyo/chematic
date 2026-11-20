//! Linear RDK bond sets with the upstream maximum-length closure rule.

use super::{Budget, MAX_CONFIGURED_PATH, RdkitRdkConfig, RdkitRdkConfigError, roots};
use chematic_core::{AtomIdx, Molecule};
use rustc_hash::FxHashSet;
use smallvec::SmallVec;

type BondPath = SmallVec<[u32; MAX_CONFIGURED_PATH]>;

pub(super) fn enumerate(
    mol: &Molecule,
    c: &RdkitRdkConfig,
    excluded: &[bool],
    budget: &mut Budget,
    mut visit: impl FnMut(&[u32], &[u32], u32),
) -> Result<(), RdkitRdkConfigError> {
    let mut adj = vec![Vec::new(); mol.atom_count()];
    for (i, b) in mol.bonds() {
        if !excluded[i.0 as usize] {
            adj[b.atom1.0 as usize].push((b.atom2, i.0));
            adj[b.atom2.0 as usize].push((b.atom1, i.0));
        }
    }
    let mut walker = Walker {
        c,
        adj,
        degrees: vec![0; mol.atom_count()],
        atom_count: 0,
        atoms: SmallVec::new(),
        bonds: BondPath::new(),
        seen: FxHashSet::default(),
        budget,
        visit: &mut visit,
    };
    for root in roots(c) {
        walker.seen.clear();
        if let Some(root) = root {
            walker.walk_from(root)?;
        } else {
            for (a, _) in mol.atoms() {
                walker.walk_from(a)?;
            }
        }
    }
    Ok(())
}

type Visitor<'a> = dyn FnMut(&[u32], &[u32], u32) + 'a;

struct Walker<'a, 'b> {
    c: &'a RdkitRdkConfig,
    adj: Vec<Vec<(AtomIdx, u32)>>,
    degrees: Vec<u32>,
    atom_count: u32,
    atoms: SmallVec<[AtomIdx; MAX_CONFIGURED_PATH + 1]>,
    bonds: BondPath,
    seen: FxHashSet<BondPath>,
    budget: &'b mut Budget,
    visit: &'b mut Visitor<'a>,
}

impl Walker<'_, '_> {
    fn walk_from(&mut self, root: AtomIdx) -> Result<(), RdkitRdkConfigError> {
        self.atoms.push(root);
        self.walk(root)?;
        self.atoms.pop();
        Ok(())
    }

    fn walk(&mut self, end: AtomIdx) -> Result<(), RdkitRdkConfigError> {
        let len = self.bonds.len();
        let final_step = len + 1 == self.c.max_path;
        for pos in 0..self.adj[end.0 as usize].len() {
            self.budget.step()?;
            let (next, bond) = self.adj[end.0 as usize][pos];
            let closes = final_step && self.atoms.len() > 2 && self.atoms[len - 1] != next;
            if self.atoms.contains(&next) && !closes {
                continue;
            }
            self.atoms.push(next);
            self.bonds.push(bond);
            for a in [end, next] {
                let d = &mut self.degrees[a.0 as usize];
                self.atom_count += u32::from(*d == 0);
                *d += 1;
            }
            self.emit()?;
            if !final_step {
                self.walk(next)?;
            }
            for a in [end, next] {
                let d = &mut self.degrees[a.0 as usize];
                *d -= 1;
                self.atom_count -= u32::from(*d == 0);
            }
            self.atoms.pop();
            self.bonds.pop();
        }
        Ok(())
    }

    fn emit(&mut self) -> Result<(), RdkitRdkConfigError> {
        if self.bonds.len() >= self.c.min_path {
            let mut key = self.bonds.clone();
            key.sort_unstable();
            if self.seen.insert(key) {
                self.budget.select()?;
                (self.visit)(&self.bonds, &self.degrees, self.atom_count);
            }
        }
        Ok(())
    }
}
