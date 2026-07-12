//! Plain acyclic hydrocarbon naming: alkanes, alkenes, alkynes, and their
//! branched forms.

use crate::helpers::{
    alkane_stem, alkane_suffix, alkene_suffix, alkyne_suffix, count_c_chain,
    find_longest_c_chain_candidates, format_substituents, unsaturation_locant,
};
use crate::{IupacError, Namer};
use chematic_core::{AtomIdx, BondOrder};
use std::collections::HashSet;

impl<'a> Namer<'a> {
    // -----------------------------------------------------------------------
    // Acyclic hydrocarbon naming
    // -----------------------------------------------------------------------

    pub(crate) fn name_acyclic_hydrocarbon(
        &self,
        carbons: &[AtomIdx],
    ) -> Result<String, IupacError> {
        let mol = self.mol;
        let n = carbons.len();

        let double_bonds = mol
            .bonds()
            .filter(|(_, b)| b.order == BondOrder::Double)
            .count();
        let triple_bonds = mol
            .bonds()
            .filter(|(_, b)| b.order == BondOrder::Triple)
            .count();
        if double_bonds > 1 || triple_bonds > 1 || (double_bonds > 0 && triple_bonds > 0) {
            return Err(IupacError::NotSupported);
        }

        // Check for branching.
        let c_set: HashSet<AtomIdx> = carbons.iter().copied().collect();
        let is_branched = carbons.iter().any(|&c| {
            mol.neighbors(c)
                .filter(|(nb, _)| c_set.contains(nb))
                .count()
                > 2
        });

        if is_branched {
            // Only saturated branched alkanes supported for now.
            if double_bonds > 0 || triple_bonds > 0 {
                return Err(IupacError::NotSupported);
            }
            return self.name_branched_alkane(carbons);
        }

        if triple_bonds == 1 {
            if n >= 4 {
                let pos = unsaturation_locant(mol, carbons, BondOrder::Triple);
                Ok(format!("{}-{}-yne", alkane_stem(n), pos))
            } else {
                Ok(alkyne_suffix(n))
            }
        } else if double_bonds == 1 {
            if n >= 4 {
                let pos = unsaturation_locant(mol, carbons, BondOrder::Double);
                Ok(format!("{}-{}-ene", alkane_stem(n), pos))
            } else {
                Ok(alkene_suffix(n))
            }
        } else {
            Ok(alkane_suffix(n))
        }
    }

    // -----------------------------------------------------------------------
    // Branched alkane naming (e.g., "2-methylpropane", "2,2-dimethylpropane")
    // -----------------------------------------------------------------------

    pub(crate) fn name_branched_alkane(&self, carbons: &[AtomIdx]) -> Result<String, IupacError> {
        let mol = self.mol;
        let all_c_set: HashSet<AtomIdx> = carbons.iter().copied().collect();

        // For one candidate chain: collect substituents (chain_position_1based,
        // alkyl_length), then apply IUPAC's lowest-locant rule (try forward and
        // reverse numbering, keep whichever gives substituents the smaller
        // first locant). Returns None if a substituent exceeds the supported
        // length (butyl) -- that disqualifies THIS candidate, not the whole
        // molecule, since a differently-chosen tied-length chain might route
        // through what would otherwise be an oversized substituent instead.
        let subs_for = |chain: &[AtomIdx]| -> Option<Vec<(usize, usize)>> {
            let n = chain.len();
            let chain_set: HashSet<AtomIdx> = chain.iter().copied().collect();
            let mut subs: Vec<(usize, usize)> = Vec::new();
            for (pos0, &chain_c) in chain.iter().enumerate() {
                let position = pos0 + 1;
                for (nb, _) in mol.neighbors(chain_c) {
                    if all_c_set.contains(&nb) && !chain_set.contains(&nb) {
                        let sub_len = count_c_chain(mol, nb, chain_c);
                        if sub_len > 4 {
                            return None;
                        }
                        subs.push((position, sub_len));
                    }
                }
            }
            if subs.is_empty() {
                return None;
            }
            let subs_rev: Vec<(usize, usize)> =
                subs.iter().map(|&(pos, len)| (n + 1 - pos, len)).collect();
            let first_fwd = subs.iter().map(|&(p, _)| p).min().unwrap_or(usize::MAX);
            let first_rev = subs_rev.iter().map(|&(p, _)| p).min().unwrap_or(usize::MAX);
            Some(if first_fwd <= first_rev {
                subs
            } else {
                subs_rev
            })
        };

        // Consider every chain tied for maximum length (IUPAC P-44.3: among
        // chains of equal length, prefer the one with MORE substituents) --
        // not just one arbitrary pick. Picking only one, with no substituent-
        // count comparison, is exactly what caused the Round 14 #4 regression
        // (see find_longest_c_chain_candidates's doc comment for the story).
        struct Best {
            chain_len: usize,
            first_locant: usize,
            subs: Vec<(usize, usize)>,
        }

        let candidates = find_longest_c_chain_candidates(mol, carbons);
        let mut best: Option<Best> = None;
        for chain in &candidates {
            let chain_len = chain.len();
            if chain_len < 2 {
                continue;
            }
            let Some(subs) = subs_for(chain) else {
                continue;
            };
            let first_locant = subs.iter().map(|&(p, _)| p).min().unwrap_or(usize::MAX);
            let is_better = match &best {
                None => true,
                Some(b) => {
                    subs.len() > b.subs.len()
                        || (subs.len() == b.subs.len() && first_locant < b.first_locant)
                }
            };
            if is_better {
                best = Some(Best {
                    chain_len,
                    first_locant,
                    subs,
                });
            }
        }

        // Any residual tie (same substituent count AND same first locant across
        // multiple candidates) falls through to whichever candidate
        // find_longest_c_chain_candidates listed first -- deterministic for a
        // fixed input, not spelling-invariant. Documented, out-of-scope-for-now
        // limitation, consistent with this round's E/Z and standardize.rs
        // precedent (see MAX_CHAIN_CANDIDATES's doc comment).
        let Some(Best {
            chain_len: n,
            subs: best_subs,
            ..
        }) = best
        else {
            return Err(IupacError::NotSupported);
        };

        Ok(format!(
            "{}{}",
            format_substituents(&best_subs),
            alkane_suffix(n)
        ))
    }
}
