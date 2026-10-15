//! Bounded RDKit-derived connected subgraphs and rooted linear atom walks.
//! Copyright (C) 2003-2026 Greg Landrum and other RDKit contributors.
//! BSD-3-Clause: see THIRD_PARTY_NOTICES.md.
use super::*;
mod linear;
type Visitor<'a> = dyn FnMut(&[u32], usize) + 'a;
pub(super) struct Budget {
    pub visits: u64,
    pub paths: usize,
    visit_limit: u64,
    path_limit: usize,
}
impl Budget {
    pub fn new(c: &RdkitLayeredConfig) -> Self {
        Self {
            visits: 0,
            paths: 0,
            visit_limit: c.max_visits,
            path_limit: c.max_paths,
        }
    }
    fn step(&mut self) -> Result<(), RdkitLayeredError> {
        if self.visits >= self.visit_limit {
            return Err(error("visit_limit", "path enumeration budget exhausted"));
        }
        self.visits += 1;
        Ok(())
    }
    fn record(&mut self, n: usize) -> Result<(), RdkitLayeredError> {
        if n > self.path_limit - self.paths {
            return Err(error("path_limit", "selected path count exceeds limit"));
        }
        self.paths += n;
        Ok(())
    }
    fn frontier(&self, n: usize) -> Result<(), RdkitLayeredError> {
        if n > self.path_limit {
            return Err(error(
                "frontier_limit",
                "linear path frontier exceeds path limit",
            ));
        }
        Ok(())
    }
}
pub(super) fn enumerate(
    mol: &Molecule,
    c: &RdkitLayeredConfig,
    budget: &mut Budget,
    mut visit: impl FnMut(&[u32], usize),
) -> Result<(), RdkitLayeredError> {
    if !c.branched_paths {
        return linear::enumerate(mol, c, budget, &mut visit);
    }
    let map = build_bond_adjacency(mol);
    let nbrs: Vec<_> = (0..mol.bond_count())
        .map(|i| map.get(&(i as u32)).cloned().unwrap_or_default())
        .collect();
    let mut roots: Vec<_> = map.keys().copied().collect();
    roots.sort_unstable();
    let mut forbidden = vec![false; mol.bond_count()];
    let mut walker = Branched {
        mol,
        c,
        budget,
        visit: &mut visit,
        nbrs,
    };
    for root in roots {
        walker.budget.step()?;
        forbidden[root as usize] = true;
        walker.walk(
            vec![root],
            walker.nbrs[root as usize].clone(),
            forbidden.clone(),
        )?;
    }
    Ok(())
}
struct Branched<'a, 'b> {
    mol: &'a Molecule,
    c: &'a RdkitLayeredConfig,
    budget: &'b mut Budget,
    visit: &'b mut Visitor<'a>,
    nbrs: Vec<Vec<u32>>,
}
impl Branched<'_, '_> {
    fn walk(
        &mut self,
        path: Vec<u32>,
        mut candidates: Vec<u32>,
        mut forbidden: Vec<bool>,
    ) -> Result<(), RdkitLayeredError> {
        self.budget.step()?;
        if path.len() >= self.c.min_path {
            let multiplicity = if self.c.from_atoms.is_empty() {
                1
            } else {
                self.c
                    .from_atoms
                    .iter()
                    .filter(|&&root| {
                        path.iter().any(|&i| {
                            let b = self.mol.bond(BondIdx(i));
                            b.atom1 == root || b.atom2 == root
                        })
                    })
                    .count()
            };
            if multiplicity > 0 {
                self.budget.record(multiplicity)?;
                (self.visit)(&path, multiplicity);
            }
        }
        if path.len() == self.c.max_path {
            return Ok(());
        }
        while let Some(next) = candidates.pop() {
            self.budget.step()?;
            if forbidden[next as usize] {
                continue;
            }
            forbidden[next as usize] = true;
            let mut stack = candidates.clone();
            stack.extend(
                self.nbrs[next as usize]
                    .iter()
                    .filter(|&&i| !forbidden[i as usize])
                    .copied(),
            );
            let mut extended = path.clone();
            extended.push(next);
            self.walk(extended, stack, forbidden.clone())?;
        }
        Ok(())
    }
}
