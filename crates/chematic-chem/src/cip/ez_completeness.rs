use super::*;

// ---------------------------------------------------------------------------
// E/Z double-bond stereo completeness
// ---------------------------------------------------------------------------

/// Summary of E/Z-stereogenic double bonds in a molecule.
///
/// The double-bond analog of `chematic_perception::stereo_validation::StereoCompleteness`
/// for tetrahedral centers: `specified` and `unspecified` both count toward `total`
/// equally, so a caller can ask "how many potential E/Z centers exist" independent of
/// whether they were annotated with a `/`/`\` marker in the input.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct EzCompleteness {
    /// Stereogenic double bonds with a `/`/`\` marker that [`assign_ez`] resolves.
    pub specified: usize,
    /// Stereogenic double bonds with no resolvable marker.
    pub unspecified: usize,
    /// `specified + unspecified`.
    pub total: usize,
}

/// Whether `alkene_end`'s two positions (other than the double-bond partner at
/// `other_end`) are distinguishable -- i.e. this end is not a terminal `=CH2`-style
/// end and not symmetrically disubstituted.
///
/// Mirrors the substituent collection [`assign_ez`]/[`highest_stereo_sub`] already use:
/// 0 heavy substituents means both positions are implicit H (terminal, matching
/// `assign_ez`'s own `subs.is_empty()` "terminal alkene" check); 1 heavy substituent is
/// distinguishable from the implicit H filling the other position as long as that
/// substituent isn't *itself* an explicit hydrogen atom (e.g. `[H]C=C[H]` from an SDF
/// with explicit H) -- CIP ranks any other heavy atom above H, so there's no tie to
/// check in that case; 2 heavy substituents are compared via [`compare_branches`], the
/// same CIP machinery `highest_stereo_sub` uses to detect a genuine priority tie (e.g.
/// two identical methyls).
fn stereogenic_end(mol: &Molecule, other_end: AtomIdx, alkene_end: AtomIdx) -> bool {
    let subs: Vec<AtomIdx> = mol
        .neighbors(alkene_end)
        .filter(|&(nb, bidx)| nb != other_end && mol.bond(bidx).order != BondOrder::Double)
        .map(|(nb, _)| nb)
        .collect();

    match subs.len() {
        0 => false, // terminal =CH2-style end
        1 => mol.atom(subs[0]).element.atomic_number() != 1,
        2 => compare_branches(mol, alkene_end, subs[0], subs[1]) != std::cmp::Ordering::Equal,
        _ => false, // not a normal trigonal alkene carbon
    }
}

/// Size of the smallest ring passing through `bidx`, or `None` if `bidx` is acyclic.
///
/// A double bond inside a small ring can't isomerize between cis/trans (the ring
/// constrains the geometry), so it isn't E/Z-stereogenic even when its two ends are
/// otherwise distinguishable. RDKit applies the same cutoff (`minBondRingSize < 8` is
/// excluded, see `Code/GraphMol/Chirality.cpp`'s `isBondPotentialStereoBond`), confirmed
/// directly against a live RDKit oracle here (`Chem.FindPotentialStereo`):
/// cyclohexene/cycloheptene (6- and 7-membered) report no potential stereo bond,
/// cyclooctene (8-membered) does.
///
/// Computed as the shortest path between the bond's two endpoints that doesn't use the
/// bond itself, plus 1 -- deliberately *not* routed through SSSR: SSSR can return a
/// large fundamental cycle instead of a smaller component ring in fused/bridged systems
/// (see this crate's aromaticity handling for the same caveat), which would silently
/// under-flag a genuinely ring-locked bond as large enough to be stereogenic. A BFS
/// shortest cycle is immune to that decomposition choice.
fn bond_min_ring_size(mol: &Molecule, bidx: BondIdx) -> Option<usize> {
    let bond = mol.bond(bidx);
    let (start, goal) = (bond.atom1, bond.atom2);

    let mut dist = vec![usize::MAX; mol.atom_count()];
    dist[start.0 as usize] = 0;
    let mut queue: VecDeque<AtomIdx> = VecDeque::from([start]);
    while let Some(cur) = queue.pop_front() {
        for (nb, b) in mol.neighbors(cur) {
            if b == bidx || dist[nb.0 as usize] != usize::MAX {
                continue;
            }
            dist[nb.0 as usize] = dist[cur.0 as usize] + 1;
            if nb == goal {
                return Some(dist[nb.0 as usize] + 1);
            }
            queue.push_back(nb);
        }
    }
    None
}

/// Rings below this size can't support E/Z isomerism (ring strain locks the geometry).
/// Matches RDKit's `isBondPotentialStereoBond` cutoff -- see [`bond_min_ring_size`].
const MIN_STEREOGENIC_RING_SIZE: usize = 8;

/// Summarise how many E/Z-stereogenic double bonds in `mol` have an explicit
/// `/`/`\` marker vs. are left unannotated.
///
/// A double bond is E/Z-stereogenic when it is non-aromatic (aromatic ring bonds carry
/// `BondOrder::Aromatic`, not `BondOrder::Double`, so they're excluded by construction --
/// the same mechanism [`assign_ez`] relies on), not constrained by a small ring (see
/// [`bond_min_ring_size`]), and both ends have two distinguishable substituents per
/// [`stereogenic_end`]. Non-stereogenic double bonds (terminal, symmetric on either end,
/// or ring-locked) are not counted at all -- matching how
/// `chematic_perception::stereo_validation::stereo_completeness` only counts atoms with
/// 4 distinct neighbors for tetrahedral centers.
pub fn ez_completeness(mol: &Molecule) -> EzCompleteness {
    let mut specified = 0usize;
    let mut unspecified = 0usize;

    for j in 0..mol.bond_count() {
        let bidx = BondIdx(j as u32);
        let bond = mol.bond(bidx);
        if bond.order != BondOrder::Double {
            continue;
        }
        if bond_min_ring_size(mol, bidx).is_some_and(|size| size < MIN_STEREOGENIC_RING_SIZE) {
            continue;
        }
        let a1 = bond.atom1;
        let a2 = bond.atom2;
        if !stereogenic_end(mol, a2, a1) || !stereogenic_end(mol, a1, a2) {
            continue;
        }

        if assign_ez(mol, bidx).is_some() {
            specified += 1;
        } else {
            unspecified += 1;
        }
    }

    EzCompleteness {
        specified,
        unspecified,
        total: specified + unspecified,
    }
}
