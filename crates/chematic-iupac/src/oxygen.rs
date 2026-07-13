//! Oxygen-containing compound naming: alcohols, ethers, aldehydes, ketones,
//! carboxylic acids, esters.

use crate::helpers::{
    alkane_base, alkane_stem, alkane_suffix, anchor_chain_and_substituents, count_c_chain,
    find_longest_c_chain, format_substituents,
};
use crate::{IupacError, Namer};
use chematic_core::{AtomIdx, BondOrder, implicit_hcount};
use std::collections::HashSet;

impl<'a> Namer<'a> {
    // -----------------------------------------------------------------------
    // One-oxygen compound: alcohol / aldehyde / ketone
    // -----------------------------------------------------------------------

    pub(crate) fn name_one_oxygen(
        &self,
        carbons: &[AtomIdx],
        o_idx: AtomIdx,
    ) -> Result<String, IupacError> {
        let mol = self.mol;
        let is_double = mol
            .neighbors(o_idx)
            .any(|(_, bi)| mol.bond(bi).order == BondOrder::Double);

        if !is_double {
            // Ether check: O with 2 C neighbors and no implicit H → R-O-R
            let o_c_nb: Vec<AtomIdx> = mol
                .neighbors(o_idx)
                .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() == 6)
                .map(|(nb, _)| nb)
                .collect();
            if o_c_nb.len() == 2 && implicit_hcount(mol, o_idx) == 0 {
                return self.name_ether(carbons, o_idx, o_c_nb[0], o_c_nb[1]);
            }

            // Find the OH carbon.
            let oh_c = mol
                .neighbors(o_idx)
                .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() == 6)
                .map(|(nb, _)| nb)
                .next()
                .ok_or(IupacError::NotSupported)?;

            // Check for branching.
            let c_set: HashSet<AtomIdx> = carbons.iter().copied().collect();
            let is_branched = carbons.iter().any(|&c| {
                mol.neighbors(c)
                    .filter(|(nb, _)| c_set.contains(nb))
                    .count()
                    > 2
            });
            if is_branched {
                return self.name_branched_alcohol(carbons, oh_c);
            }

            // Straight-chain: determine OH position using longest chain.
            let chain = find_longest_c_chain(mol, carbons);
            let n = chain.len();
            let pos_fwd = chain
                .iter()
                .position(|&c| c == oh_c)
                .map(|p| p + 1)
                .unwrap_or(1);
            let pos = pos_fwd.min(n + 1 - pos_fwd);
            if pos == 1 && n <= 2 {
                // Short common names without locant: methanol, ethanol.
                return Ok(format!("{}anol", alkane_stem(n)));
            }
            return Ok(format!("{}-{}-ol", alkane_base(n), pos));
        }

        // Carbonyl: find the C=O carbon.
        let carbonyl_c = mol
            .neighbors(o_idx)
            .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() == 6)
            .map(|(nb, _)| nb)
            .next()
            .ok_or(IupacError::NotSupported)?;

        if implicit_hcount(mol, carbonyl_c) > 0 {
            // Aldehyde: CHO is position 1; find chain from carbonyl_c, picking
            // among tied-longest candidates by substituent count (IUPAC P-44.3).
            let c_set: HashSet<AtomIdx> = carbons.iter().copied().collect();
            let (chain, subs) = anchor_chain_and_substituents(mol, &c_set, carbonyl_c)
                .ok_or(IupacError::NotSupported)?;
            let n = chain.len();
            let prefix = if subs.is_empty() {
                String::new()
            } else {
                format_substituents(&subs)
            };
            return Ok(format!("{}{}anal", prefix, alkane_stem(n)));
        }

        // Ketone: internal C=O — find principal chain and position.
        let chain = find_longest_c_chain(mol, carbons);
        let n = chain.len();
        if n < 3 {
            return Err(IupacError::NotSupported);
        }
        let chain_set: HashSet<AtomIdx> = chain.iter().copied().collect();
        let all_c_set: HashSet<AtomIdx> = carbons.iter().copied().collect();
        let pos_fwd = chain
            .iter()
            .position(|&c| c == carbonyl_c)
            .map(|p| p + 1)
            .ok_or(IupacError::NotSupported)?;
        let pos = pos_fwd.min(n + 1 - pos_fwd);
        let reversed = pos_fwd > n + 1 - pos_fwd;
        // Collect alkyl substituents on the chain.
        let mut subs: Vec<(usize, usize)> = Vec::new();
        for (idx, &chain_c) in chain.iter().enumerate() {
            let position = idx + 1;
            for (nb, _) in mol.neighbors(chain_c) {
                if all_c_set.contains(&nb) && !chain_set.contains(&nb) {
                    let sub_len = count_c_chain(mol, nb, chain_c);
                    if sub_len > 4 {
                        return Err(IupacError::NotSupported);
                    }
                    let adj_pos = if reversed { n + 1 - position } else { position };
                    subs.push((adj_pos, sub_len));
                }
            }
        }
        let prefix = if subs.is_empty() {
            String::new()
        } else {
            format_substituents(&subs)
        };
        Ok(format!("{}{}-{}-one", prefix, alkane_base(n), pos))
    }

    // -----------------------------------------------------------------------
    // Ether naming (R-O-R → "alkoxyalkane")
    // -----------------------------------------------------------------------

    pub(crate) fn name_ether(
        &self,
        carbons: &[AtomIdx],
        o_idx: AtomIdx,
        side_a: AtomIdx,
        side_b: AtomIdx,
    ) -> Result<String, IupacError> {
        let mol = self.mol;
        // Only unbranched ethers.
        let c_set: HashSet<AtomIdx> = carbons.iter().copied().collect();
        if carbons.iter().any(|&c| {
            mol.neighbors(c)
                .filter(|(nb, _)| c_set.contains(nb))
                .count()
                > 2
        }) {
            return Err(IupacError::NotSupported);
        }
        let len_a = count_c_chain(mol, side_a, o_idx);
        let len_b = count_c_chain(mol, side_b, o_idx);
        let (alkoxy_len, parent_len) = if len_a <= len_b {
            (len_a, len_b)
        } else {
            (len_b, len_a)
        };
        let alkoxy = format!("{}oxy", alkane_stem(alkoxy_len));
        let parent = alkane_suffix(parent_len);
        // Add locant "1-" when parent ≥ 3 C and chains differ (O position is ambiguous).
        if parent_len >= 3 && alkoxy_len != parent_len {
            Ok(format!("1-{alkoxy}{parent}"))
        } else {
            Ok(format!("{alkoxy}{parent}"))
        }
    }

    // -----------------------------------------------------------------------
    // Two-oxygen compound: carboxylic acid or ester
    // -----------------------------------------------------------------------

    pub(crate) fn name_two_oxygens(
        &self,
        carbons: &[AtomIdx],
        o_atoms: &[AtomIdx],
    ) -> Result<String, IupacError> {
        let mol = self.mol;
        let o1 = o_atoms[0];
        let o2 = o_atoms[1];

        let o1_dbl = mol
            .neighbors(o1)
            .any(|(_, bi)| mol.bond(bi).order == BondOrder::Double);
        let o2_dbl = mol
            .neighbors(o2)
            .any(|(_, bi)| mol.bond(bi).order == BondOrder::Double);

        let (carbonyl_o, ester_o) = match (o1_dbl, o2_dbl) {
            (true, false) => (o1, o2),
            (false, true) => (o2, o1),
            _ => return Err(IupacError::NotSupported),
        };

        // Carbonyl C is bonded to the =O oxygen.
        let carbonyl_c = mol
            .neighbors(carbonyl_o)
            .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() == 6)
            .map(|(nb, _)| nb)
            .next()
            .ok_or(IupacError::NotSupported)?;

        // Carbonyl C must also be bonded to the single-bond O.
        if !mol.neighbors(carbonyl_c).any(|(nb, _)| nb == ester_o) {
            return Err(IupacError::NotSupported);
        }

        // Is the single-bond O also bonded to another C (→ ester) or only H (→ acid)?
        let alcohol_c = mol
            .neighbors(ester_o)
            .filter(|(nb, _)| *nb != carbonyl_c && mol.atom(*nb).element.atomic_number() == 6)
            .map(|(nb, _)| nb)
            .next();

        let c_set: HashSet<AtomIdx> = carbons.iter().copied().collect();
        if let Some(alc_c) = alcohol_c {
            // Ester: find acid chain from carbonyl_c (handles branched acid parts),
            // picking among tied-longest candidates by substituent count.
            let (chain_acid, subs) = anchor_chain_and_substituents(mol, &c_set, carbonyl_c)
                .ok_or(IupacError::NotSupported)?;
            let acid_n = chain_acid.len();
            let alcohol_n = count_c_chain(mol, alc_c, ester_o);
            let acid_part = if subs.is_empty() {
                format!("{}anoate", alkane_stem(acid_n))
            } else {
                format!(
                    "{}{}anoate",
                    format_substituents(&subs),
                    alkane_stem(acid_n)
                )
            };
            Ok(format!("{}yl {}", alkane_stem(alcohol_n), acid_part))
        } else {
            // Carboxylic acid — find principal chain from carboxyl C (always
            // position 1), picking among tied-longest candidates by substituent count.
            let (chain, subs) = anchor_chain_and_substituents(mol, &c_set, carbonyl_c)
                .ok_or(IupacError::NotSupported)?;
            let n = chain.len();
            if subs.is_empty() {
                Ok(format!("{}anoic acid", alkane_stem(n)))
            } else {
                Ok(format!(
                    "{}{}anoic acid",
                    format_substituents(&subs),
                    alkane_stem(n)
                ))
            }
        }
    }

    // -----------------------------------------------------------------------
    // Branched alcohol naming (e.g., "propan-2-ol")
    // -----------------------------------------------------------------------

    pub(crate) fn name_branched_alcohol(
        &self,
        carbons: &[AtomIdx],
        oh_c: AtomIdx,
    ) -> Result<String, IupacError> {
        // Find principal chain.
        let chain = find_longest_c_chain(self.mol, carbons);
        let n = chain.len();
        if n < 2 {
            return Err(IupacError::NotSupported);
        }

        let chain_set: HashSet<AtomIdx> = chain.iter().copied().collect();
        let all_c_set: HashSet<AtomIdx> = carbons.iter().copied().collect();

        // The OH carbon must be on the principal chain.
        let pos_on_chain = if chain_set.contains(&oh_c) {
            chain.iter().position(|&c| c == oh_c).map(|p| p + 1)
        } else {
            None
        };

        let pos_fwd = pos_on_chain.ok_or(IupacError::NotSupported)?;
        let pos = pos_fwd.min(n + 1 - pos_fwd);

        // Also collect any alkyl substituents on the chain.
        let mut subs: Vec<(usize, usize)> = Vec::new();
        for (pos0, &chain_c) in chain.iter().enumerate() {
            let position = pos0 + 1;
            for (nb, _) in self.mol.neighbors(chain_c) {
                if all_c_set.contains(&nb) && !chain_set.contains(&nb) {
                    let sub_len = count_c_chain(self.mol, nb, chain_c);
                    if sub_len > 4 {
                        return Err(IupacError::NotSupported);
                    }
                    subs.push((position, sub_len));
                }
            }
        }

        // Re-number subs with the same locant direction as OH position.
        if pos_fwd > n + 1 - pos_fwd {
            // Reverse direction was chosen for OH; re-number subs accordingly.
            subs = subs.iter().map(|&(p, l)| (n + 1 - p, l)).collect();
        }

        let prefix = if subs.is_empty() {
            String::new()
        } else {
            subs.sort_unstable();
            let subs_rev: Vec<(usize, usize)> = subs.iter().map(|&(p, l)| (n + 1 - p, l)).collect();
            let first_fwd = subs.iter().map(|&(p, _)| p).min().unwrap_or(usize::MAX);
            let first_rev = subs_rev.iter().map(|&(p, _)| p).min().unwrap_or(usize::MAX);
            let best = if first_fwd <= first_rev {
                subs.clone()
            } else {
                subs_rev
            };
            format!("{}-", format_substituents(&best))
        };

        Ok(format!("{}{}-{}-ol", prefix, alkane_base(n), pos))
    }
}
