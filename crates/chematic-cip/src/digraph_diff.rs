//! Structural diff between two [`CipDigraph`]s built for "the same" conceptual
//! stereocenter under different, chemically-equivalent input representations (aromatic
//! notation vs. a Kekulé respelling, an atom-renumbering, a SMILES respelling) --
//! Milestone 3B-0's evidence artifact for *why* the current (pre-M3B) digraph is
//! representation-dependent on mancude ring systems, collected before any comparator
//! change is made.
//!
//! Full per-position tree isomorphism isn't needed for this diagnostic purpose (and
//! would need real atom-identity tracking through arbitrary traversal-order changes,
//! which SMILES respelling in particular doesn't preserve). Instead: at each BFS depth
//! from the root, count how many nodes of each [`CipNodeKind`] tag appear, and report the
//! shallowest depth where the two representations' per-depth counts disagree. That's
//! sufficient to show, concretely, *where* one representation's digraph has (for example)
//! a `MultipleBondDuplicate` the other lacks -- exactly what an aromatic-vs-Kekulé
//! comparison is expected to surface, since `BondOrder::Aromatic.order_int() == 1`
//! contributes none today while an explicit double bond does.

use std::collections::BTreeMap;

use chematic_core::{AtomIdx, Molecule, MoleculeBuilder, STEREO_H_SENTINEL};

use crate::CipError;
use crate::budget::CipBudget;
use crate::digraph::CipDigraph;
use crate::node::CipNodeKind;

/// Rebuild `mol` with atoms inserted in the order given by `perm` (`perm[new_idx] =
/// old_idx`), remapping bond endpoints and `stereo_neighbor_order` accordingly. Returns
/// the permuted molecule plus an `old_idx -> new_idx` map so a caller can locate where a
/// specific atom ended up. Same shape as the `#[cfg(test)]`-private helper in
/// `src/tests.rs`, kept as a separate small `pub` copy here (not shared) since this one
/// must be reachable from `tests/*.rs` integration tests, which only see a crate's public
/// surface, not its internal unit-test module.
pub fn renumber_molecule(mol: &Molecule, perm: &[usize]) -> (Molecule, Vec<u32>) {
    let mut old_to_new = vec![0u32; perm.len()];
    for (new_idx, &old_idx) in perm.iter().enumerate() {
        old_to_new[old_idx] = new_idx as u32;
    }
    let mut builder = MoleculeBuilder::new();
    for &old_idx in perm {
        builder.add_atom(mol.atom(AtomIdx(old_idx as u32)).clone());
    }
    for (_, bond) in mol.bonds() {
        let a = AtomIdx(old_to_new[bond.atom1.0 as usize]);
        let b = AtomIdx(old_to_new[bond.atom2.0 as usize]);
        let _ = builder.add_bond(a, b, bond.order);
    }
    for (old_idx, _) in mol.atoms() {
        if let Some(order) = mol.stereo_neighbor_order(old_idx) {
            let remapped: Vec<u32> = order
                .iter()
                .map(|&v| {
                    if v == STEREO_H_SENTINEL {
                        v
                    } else {
                        old_to_new[v as usize]
                    }
                })
                .collect();
            builder.set_stereo_neighbor_order(AtomIdx(old_to_new[old_idx.0 as usize]), remapped);
        }
    }
    (builder.build(), old_to_new)
}

/// Find the atom in `mol` whose (element, charge, aromatic flag, degree,
/// chirality-tagged-or-not) signature matches `reference_atom` in `reference_mol` --
/// used to locate "the same" stereocenter after a transformation that doesn't preserve
/// atom indices (e.g. a SMILES respelling through a different canonical traversal).
/// Degree (heavy-atom neighbor count) is included because this corpus's own molecules
/// are almost all multi-stereocenter (every `aromatic_mancude` case has 2+ `@`/`@@` tags)
/// -- without it, "element + chirality-tagged" alone is ambiguous on nearly every real
/// case here, since several differently-substituted stereocenters share it. Returns
/// `None` if zero or more than one atom matches even with degree included (an ambiguous
/// match is not used silently, not guessed).
pub fn find_atom_by_signature(
    reference_mol: &Molecule,
    reference_atom: AtomIdx,
    mol: &Molecule,
) -> Option<AtomIdx> {
    let want = reference_mol.atom(reference_atom);
    let want_degree = reference_mol.neighbors(reference_atom).count();
    let mut matches = mol.atoms().filter(|(idx, a)| {
        a.element == want.element
            && a.charge == want.charge
            && a.aromatic == want.aromatic
            && mol.neighbors(*idx).count() == want_degree
            && (a.chirality != chematic_core::Chirality::None)
                == (want.chirality != chematic_core::Chirality::None)
    });
    let (first, _) = matches.next()?;
    if matches.next().is_some() {
        None // ambiguous: more than one candidate, don't guess
    } else {
        Some(first)
    }
}

/// Per-depth histogram of node kind tags reachable from a digraph's root, up to
/// `CipBudget`'s own limits.
fn depth_histogram(
    mol: &Molecule,
    atom: AtomIdx,
    budget: CipBudget,
) -> Result<Vec<BTreeMap<&'static str, usize>>, CipError> {
    let mut graph = CipDigraph::new(mol, atom, budget)?;
    let mut histograms: Vec<BTreeMap<&'static str, usize>> = Vec::new();
    let mut frontier = vec![graph.root()];
    while !frontier.is_empty() {
        let depth = histograms.len();
        let mut level = BTreeMap::new();
        let mut next = Vec::new();
        for node in &frontier {
            let kind = graph.node(*node).kind;
            *level.entry(kind_tag(kind)).or_insert(0) += 1;
            if matches!(kind, CipNodeKind::Atom { .. }) {
                let children = graph.expand_children(*node)?;
                next.extend(children);
            }
        }
        histograms.push(level);
        frontier = next;
        if depth > 64 {
            break; // pathological depth guard; not expected for real corpus molecules
        }
    }
    Ok(histograms)
}

fn kind_tag(kind: CipNodeKind) -> &'static str {
    match kind {
        CipNodeKind::Atom { .. } => "Atom",
        CipNodeKind::MultipleBondDuplicate { .. } => "MultipleBondDuplicate",
        CipNodeKind::RingDuplicate { .. } => "RingDuplicate",
        CipNodeKind::ImplicitHydrogen => "ImplicitHydrogen",
    }
}

/// The shallowest depth where two representations' digraphs disagree on their per-kind
/// node counts, and what each side looked like at that depth.
#[derive(Debug, Clone)]
pub struct DigraphDivergence {
    pub depth: usize,
    pub left: BTreeMap<&'static str, usize>,
    pub right: BTreeMap<&'static str, usize>,
}

/// Build the digraph for `(mol, atom)` and `(other_mol, other_atom)` and report the first
/// depth where their per-kind node-count histograms disagree, or `None` if they agree at
/// every depth both trees reach.
pub fn first_divergence(
    mol: &Molecule,
    atom: AtomIdx,
    other_mol: &Molecule,
    other_atom: AtomIdx,
    budget: CipBudget,
) -> Result<Option<DigraphDivergence>, CipError> {
    let left = depth_histogram(mol, atom, budget)?;
    let right = depth_histogram(other_mol, other_atom, budget)?;
    for depth in 0..left.len().max(right.len()) {
        let l = left.get(depth).cloned().unwrap_or_default();
        let r = right.get(depth).cloned().unwrap_or_default();
        if l != r {
            return Ok(Some(DigraphDivergence {
                depth,
                left: l,
                right: r,
            }));
        }
    }
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::kekulization::{apply_kekule, kekulize};

    #[test]
    fn aromatic_vs_kekule_diverges_on_a_fully_substituted_ipso() {
        // A stereocenter bonded to a phenyl ring: aromatic notation contributes zero
        // MultipleBondDuplicate nodes for the ring bonds (BondOrder::Aromatic.order_int()
        // == 1); the Kekulé-respelled form introduces one where the ring's explicit
        // double bond lands -- this is the root-cause signature M3B-1 needs to fix.
        let mol = chematic_smiles::parse("F[C@](Cl)(Br)c1ccccc1").unwrap();
        let atom = AtomIdx(1);
        let kekule = kekulize(&mol).unwrap();
        let kekule_mol = apply_kekule(&mol, &kekule);

        let divergence =
            first_divergence(&mol, atom, &kekule_mol, atom, CipBudget::default_budget())
                .unwrap()
                .expect("aromatic and Kekulé-respelled forms are expected to diverge today");
        assert!(
            divergence.right.contains_key("MultipleBondDuplicate")
                || divergence.left.contains_key("MultipleBondDuplicate"),
            "divergence should be explained by a MultipleBondDuplicate present on one \
             side and not the other: {divergence:?}"
        );
    }

    #[test]
    fn renumbering_does_not_change_a_non_aromatic_digraphs_kind_histogram() {
        let mol = chematic_smiles::parse("F[C@](Cl)(Br)CC").unwrap();
        let atom = AtomIdx(1);
        let n = mol.atom_count();
        let perm: Vec<usize> = (0..n).rev().collect();
        let (renumbered, old_to_new) = renumber_molecule(&mol, &perm);
        let new_atom = AtomIdx(old_to_new[atom.0 as usize]);

        let divergence = first_divergence(
            &mol,
            atom,
            &renumbered,
            new_atom,
            CipBudget::default_budget(),
        )
        .unwrap();
        assert!(
            divergence.is_none(),
            "renumbering a plain aliphatic molecule must not change its per-depth kind \
             histogram: {divergence:?}"
        );
    }

    #[test]
    fn find_atom_by_signature_locates_the_unique_stereocenter() {
        let mol = chematic_smiles::parse("F[C@](Cl)(Br)CC").unwrap();
        let respelled = chematic_smiles::parse("CC[C@](F)(Cl)Br").unwrap();
        let found = find_atom_by_signature(&mol, AtomIdx(1), &respelled)
            .expect("exactly one chirality-tagged carbon should match");
        assert_eq!(
            respelled.atom(found).chirality,
            mol.atom(AtomIdx(1)).chirality.clone()
        );
    }
}
