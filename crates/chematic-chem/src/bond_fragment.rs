//! Checked bridge-single-bond fragmentation with isotope-labelled dummy atoms.
use std::collections::{HashSet, VecDeque};

use chematic_core::{Atom, AtomIdx, BondIdx, BondOrder, Chirality, Molecule};

/// One dummy attached in place of the original neighbour.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BondAttachment {
    pub source_bond: BondIdx,
    pub anchor: AtomIdx,
    pub dummy: AtomIdx,
    pub isotope: u16,
}

/// Original atom indices are retained; dummy atoms are appended in cut order.
#[derive(Clone)]
pub struct SingleBondFragmentation {
    pub molecule: Molecule,
    pub cuts_per_atom: Vec<usize>,
    pub attachments: Vec<BondAttachment>,
}

/// A checked refusal; the input is never mutated.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BondFragmentError(pub String);

impl std::fmt::Display for BondFragmentError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for BondFragmentError {}

fn refuse<T>(message: impl Into<String>) -> Result<T, BondFragmentError> {
    Err(BondFragmentError(message.into()))
}

fn is_ring_bond(mol: &Molecule, cut: BondIdx) -> bool {
    let bond = mol.bond(cut);
    let mut seen = vec![false; mol.atom_count()];
    let mut queue = VecDeque::from([bond.atom1]);
    seen[bond.atom1.0 as usize] = true;
    while let Some(atom) = queue.pop_front() {
        for (neighbor, edge) in mol.neighbors(atom) {
            if edge == cut || seen[neighbor.0 as usize] {
                continue;
            }
            if neighbor == bond.atom2 {
                return true;
            }
            seen[neighbor.0 as usize] = true;
            queue.push_back(neighbor);
        }
    }
    false
}

pub(crate) fn validate_graph(mol: &Molecule) -> Result<(), BondFragmentError> {
    if mol
        .atoms()
        .any(|(_, a)| a.chirality != Chirality::None || a.cip_code.is_some())
        || mol.has_bond_directions()
        || !mol.stereo_groups().is_empty()
        || mol
            .bonds()
            .any(|(_, b)| matches!(b.order, BondOrder::Up | BondOrder::Down))
    {
        return refuse("bond fragmentation refuses stereo-annotated inputs");
    }
    if mol.atoms().any(|(_, a)| a.wildcard)
        || mol.bonds().any(|(_, b)| {
            !matches!(
                b.order,
                BondOrder::Single | BondOrder::Double | BondOrder::Triple | BondOrder::Aromatic
            )
        })
    {
        return refuse("bond fragmentation requires a concrete standard-bond graph");
    }
    Ok(())
}

fn validate_cuts(
    mol: &Molecule,
    cuts: &[BondIdx],
    labels: Option<&[(u16, u16)]>,
) -> Result<(), BondFragmentError> {
    if mol.atom_count() > 4096 || mol.bond_count() > 8192 || cuts.len() > 1024 {
        return refuse("bond fragmentation exceeds 4096 atoms, 8192 bonds or 1024 cuts");
    }
    if labels.is_some_and(|values| values.len() != cuts.len()) {
        return refuse("dummy labels must have one pair per cut");
    }
    let mut seen = HashSet::new();
    for &cut in cuts {
        if cut.0 as usize >= mol.bond_count() {
            return refuse(format!("bond index {} is out of range", cut.0));
        }
        if !seen.insert(cut) {
            return refuse(format!("duplicate bond index {}", cut.0));
        }
        if mol.bond(cut).order != BondOrder::Single {
            return refuse(format!("bond {} is not a plain single bond", cut.0));
        }
        if is_ring_bond(mol, cut) {
            return refuse(format!("bond {} is in a ring", cut.0));
        }
    }
    if !cuts.is_empty() {
        validate_graph(mol)?;
    }
    Ok(())
}

/// Replace each selected bridge single bond with two singly bonded dummies.
///
/// With `(begin_label, end_label)`, the first dummy is attached to the original
/// end atom and the second to the begin atom, as in RDKit FragmentOnBonds.
/// Default labels are the original begin/end atom indices. Label zero means an
/// unlabelled dummy. All original atom attributes, indices, H counts and bond
/// orders survive. This bounded profile refuses any stereo annotations, query/
/// pre-existing dummy atoms and nonstandard bond orders when cuts are nonempty.
/// Ring cuts, non-single cuts and duplicate indices are refused. No coordinates
/// are generated. Empty cuts clone the graph and report all-zero counters.
pub fn fragment_single_bonds_checked(
    mol: &Molecule,
    cuts: &[BondIdx],
    labels: Option<&[(u16, u16)]>,
) -> Result<SingleBondFragmentation, BondFragmentError> {
    validate_cuts(mol, cuts, labels)?;
    let mut output = mol.clone();
    let mut descending = cuts.to_vec();
    descending.sort_unstable_by_key(|cut| std::cmp::Reverse(cut.0));
    for cut in descending {
        output.remove_bond(cut);
    }
    let mut counts = vec![0; mol.atom_count()];
    let mut attachments = Vec::with_capacity(2 * cuts.len());
    for (i, &cut) in cuts.iter().enumerate() {
        let bond = mol.bond(cut);
        let pair = labels.map_or((bond.atom1.0 as u16, bond.atom2.0 as u16), |p| p[i]);
        for (anchor, isotope) in [(bond.atom2, pair.0), (bond.atom1, pair.1)] {
            let mut atom = Atom::wildcard();
            atom.isotope = (isotope != 0).then_some(isotope);
            let dummy = output.add_atom(atom);
            output
                .add_bond(anchor, dummy, BondOrder::Single)
                .map_err(|error| BondFragmentError(error.to_string()))?;
            counts[anchor.0 as usize] += 1;
            attachments.push(BondAttachment {
                source_bond: cut,
                anchor,
                dummy,
                isotope,
            });
        }
    }
    Ok(SingleBondFragmentation {
        molecule: output,
        cuts_per_atom: counts,
        attachments,
    })
}

#[cfg(test)]
#[path = "bond_fragment_tests.rs"]
mod tests;
