//! Plain acyclic hydrocarbon naming: alkanes, alkenes, alkynes, and their
//! branched forms.

use crate::helpers::{
    alkane_stem, alkane_suffix, alkene_suffix, alkyne_suffix, count_c_chain, find_longest_c_chain,
    format_substituents, unsaturation_locant,
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

        // Find the principal chain (longest C–C path).
        let chain = find_longest_c_chain(mol, carbons);
        let n = chain.len();
        if n < 2 {
            return Err(IupacError::NotSupported);
        }

        let chain_set: std::collections::HashSet<AtomIdx> = chain.iter().copied().collect();
        let all_c_set: std::collections::HashSet<AtomIdx> = carbons.iter().copied().collect();

        // Collect substituents: (chain_position_1based, alkyl_length).
        let mut subs: Vec<(usize, usize)> = Vec::new();
        for (pos0, &chain_c) in chain.iter().enumerate() {
            let position = pos0 + 1;
            for (nb, _) in mol.neighbors(chain_c) {
                if all_c_set.contains(&nb) && !chain_set.contains(&nb) {
                    // Substituent rooted at `nb`, blocked by chain_c.
                    let sub_len = count_c_chain(mol, nb, chain_c);
                    if sub_len > 4 {
                        return Err(IupacError::NotSupported);
                    }
                    subs.push((position, sub_len));
                }
            }
        }

        if subs.is_empty() {
            return Err(IupacError::NotSupported);
        }

        // Apply IUPAC lowest-locant rule: try forward and reverse numbering.
        let subs_rev: Vec<(usize, usize)> =
            subs.iter().map(|&(pos, len)| (n + 1 - pos, len)).collect();

        let first_fwd = subs.iter().map(|&(p, _)| p).min().unwrap_or(usize::MAX);
        let first_rev = subs_rev.iter().map(|&(p, _)| p).min().unwrap_or(usize::MAX);
        let best_subs = if first_fwd <= first_rev {
            subs
        } else {
            subs_rev
        };

        Ok(format!(
            "{}{}",
            format_substituents(&best_subs),
            alkane_suffix(n)
        ))
    }
}
