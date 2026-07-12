//! Nitrogen (amide/amine/nitrile), halogen (haloalkane), and sulfur
//! (sulfide/thiol) compound naming.

use crate::helpers::{
    alkane_base, alkane_stem, alkane_suffix, alkyl_prefix, anchor_chain_and_substituents,
    count_c_chain, find_longest_c_chain, format_substituents,
};
use crate::{IupacError, Namer};
use chematic_core::{AtomIdx, BondOrder, implicit_hcount};
use std::collections::HashSet;

impl<'a> Namer<'a> {
    // -----------------------------------------------------------------------
    // Amide: C(=O)–N
    // -----------------------------------------------------------------------

    pub(crate) fn name_amide(
        &self,
        _carbons: &[AtomIdx],
        o_idx: AtomIdx,
        n_idx: AtomIdx,
    ) -> Result<String, IupacError> {
        let mol = self.mol;

        // O must be a carbonyl (C=O).
        if !mol
            .neighbors(o_idx)
            .any(|(_, bi)| mol.bond(bi).order == BondOrder::Double)
        {
            return Err(IupacError::NotSupported);
        }

        let carbonyl_c = mol
            .neighbors(o_idx)
            .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() == 6)
            .map(|(nb, _)| nb)
            .next()
            .ok_or(IupacError::NotSupported)?;

        // Carbonyl C must be bonded to N.
        if !mol.neighbors(carbonyl_c).any(|(nb, _)| nb == n_idx) {
            return Err(IupacError::NotSupported);
        }

        // Only primary/secondary amides (N has ≥ 1 H).
        if implicit_hcount(mol, n_idx) == 0 {
            return Err(IupacError::NotSupported);
        }

        // Amide chain from carbonyl_c (handles branched structures), picking
        // among tied-longest candidates by substituent count (IUPAC P-44.3).
        let c_set: HashSet<AtomIdx> = mol
            .atoms()
            .filter(|(_, a)| a.element.atomic_number() == 6)
            .map(|(i, _)| i)
            .collect();
        let (chain, subs) = anchor_chain_and_substituents(mol, &c_set, carbonyl_c)
            .ok_or(IupacError::NotSupported)?;
        let n = chain.len();
        let prefix = if subs.is_empty() {
            String::new()
        } else {
            format_substituents(&subs)
        };
        Ok(format!("{}{}anamide", prefix, alkane_stem(n)))
    }

    // -----------------------------------------------------------------------
    // Amine naming
    // -----------------------------------------------------------------------

    pub(crate) fn name_amine(
        &self,
        carbons: &[AtomIdx],
        n_idx: AtomIdx,
    ) -> Result<String, IupacError> {
        let mol = self.mol;
        let n_h = implicit_hcount(mol, n_idx);
        let c_sides: Vec<AtomIdx> = mol
            .neighbors(n_idx)
            .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() == 6)
            .map(|(nb, _)| nb)
            .collect();
        let mut chain_lens: Vec<usize> = c_sides
            .iter()
            .map(|&nb| count_c_chain(mol, nb, n_idx))
            .collect();
        chain_lens.sort_unstable_by(|a, b| b.cmp(a)); // descending
        match n_h {
            2 => {
                // Find N-bearing C's position on the principal chain.
                let chain = find_longest_c_chain(mol, carbons);
                let n_chain = chain.len();
                let chain_set: HashSet<AtomIdx> = chain.iter().copied().collect();
                let amine_c = mol
                    .neighbors(n_idx)
                    .filter(|(nb, _)| {
                        mol.atom(*nb).element.atomic_number() == 6 && chain_set.contains(nb)
                    })
                    .map(|(nb, _)| nb)
                    .next()
                    .ok_or(IupacError::NotSupported)?;
                let pos_fwd = chain
                    .iter()
                    .position(|&c| c == amine_c)
                    .map(|p| p + 1)
                    .unwrap_or(1);
                let pos = pos_fwd.min(n_chain + 1 - pos_fwd);
                Ok(format!("{}an-{}-amine", alkane_stem(n_chain), pos))
            }
            1 => {
                if chain_lens.len() != 2 {
                    return Err(IupacError::NotSupported);
                }
                let parent_len = chain_lens[0];
                let sub_len = chain_lens[1];
                Ok(format!(
                    "N-{}yl{}anamine",
                    alkane_stem(sub_len),
                    alkane_stem(parent_len)
                ))
            }
            0 => {
                if chain_lens.len() != 3 {
                    return Err(IupacError::NotSupported);
                }
                let parent_len = chain_lens[0];
                let sub1 = chain_lens[1];
                let sub2 = chain_lens[2];
                if sub1 == sub2 {
                    Ok(format!(
                        "N,N-di{}yl{}anamine",
                        alkane_stem(sub1),
                        alkane_stem(parent_len)
                    ))
                } else {
                    let (lo, hi) = (sub1.min(sub2), sub1.max(sub2));
                    Ok(format!(
                        "N-{}yl-N-{}yl{}anamine",
                        alkane_stem(lo),
                        alkane_stem(hi),
                        alkane_stem(parent_len)
                    ))
                }
            }
            _ => Err(IupacError::NotSupported),
        }
    }

    // -----------------------------------------------------------------------
    // Nitrile naming (R-C≡N → "...nitrile")
    // -----------------------------------------------------------------------

    pub(crate) fn is_nitrile(&self, n_idx: AtomIdx) -> bool {
        self.mol
            .neighbors(n_idx)
            .any(|(_, bi)| self.mol.bond(bi).order == BondOrder::Triple)
    }

    pub(crate) fn name_nitrile(
        &self,
        carbons: &[AtomIdx],
        n_idx: AtomIdx,
    ) -> Result<String, IupacError> {
        let mol = self.mol;
        // Find the C≡N carbon.
        let nitrile_c = mol
            .neighbors(n_idx)
            .filter(|(_, bi)| mol.bond(*bi).order == BondOrder::Triple)
            .map(|(nb, _)| nb)
            .next()
            .ok_or(IupacError::NotSupported)?;
        // Count the total C chain (nitrile C + alkyl chain).
        // count_c_chain gives all C reachable from nitrile_c without crossing N.
        let n_carbons = count_c_chain(mol, nitrile_c, n_idx);
        // n_carbons already includes the nitrile carbon itself.
        if n_carbons == 0 {
            return Err(IupacError::NotSupported);
        }
        // Verify no branching on the C chain
        let c_set: std::collections::HashSet<AtomIdx> = carbons.iter().copied().collect();
        for &c in carbons {
            if mol
                .neighbors(c)
                .filter(|(nb, _)| c_set.contains(nb))
                .count()
                > 2
            {
                return Err(IupacError::NotSupported); // branched nitrile not supported
            }
        }
        Ok(format!("{}enitrile", alkane_base(n_carbons)))
    }

    // -----------------------------------------------------------------------
    // Haloalkane naming
    // -----------------------------------------------------------------------

    pub(crate) fn name_haloalkane(
        &self,
        carbons: &[AtomIdx],
        halogen_atoms: &[AtomIdx],
        prefix: &str,
    ) -> Result<String, IupacError> {
        let mol = self.mol;
        let chain = find_longest_c_chain(mol, carbons);
        let n = chain.len();
        let chain_set: HashSet<AtomIdx> = chain.iter().copied().collect();

        // Find the locant of each halogen on the chain.
        let mut locants: Vec<usize> = Vec::new();
        for &hal in halogen_atoms {
            let hal_c = mol
                .neighbors(hal)
                .filter(|(nb, _)| chain_set.contains(nb))
                .map(|(nb, _)| nb)
                .next()
                .ok_or(IupacError::NotSupported)?;
            let pos = chain
                .iter()
                .position(|&c| c == hal_c)
                .map(|p| p + 1)
                .ok_or(IupacError::NotSupported)?;
            locants.push(pos);
        }

        // Apply lowest-locant rule (compare forward vs reversed numbering).
        let locants_rev: Vec<usize> = locants.iter().map(|&p| n + 1 - p).collect();
        let best = if locants.iter().min() <= locants_rev.iter().min() {
            locants
        } else {
            locants_rev
        };

        let count = halogen_atoms.len();
        let mult = match count {
            1 => prefix.to_string(),
            2 => format!("di{prefix}"),
            3 => format!("tri{prefix}"),
            _ => return Err(IupacError::NotSupported),
        };

        let mut sorted_locs = best;
        sorted_locs.sort_unstable();
        let locant_str = sorted_locs
            .iter()
            .map(|l| l.to_string())
            .collect::<Vec<_>>()
            .join(",");

        // Omit locant for short unambiguous cases (n≤2, single halogen at terminal).
        if n <= 2 && count == 1 {
            Ok(format!("{mult}{}", alkane_suffix(n)))
        } else {
            Ok(format!("{locant_str}-{mult}{}", alkane_suffix(n)))
        }
    }

    // -----------------------------------------------------------------------
    // Sulfide naming (C-S-C, no SH)
    // -----------------------------------------------------------------------

    pub(crate) fn name_sulfide(
        &self,
        carbons: &[AtomIdx],
        s_idx: AtomIdx,
    ) -> Result<String, IupacError> {
        let mol = self.mol;
        // Collect both C neighbors of S.
        let c_neighbors: Vec<AtomIdx> = mol
            .neighbors(s_idx)
            .filter(|(nb, _)| carbons.contains(nb))
            .map(|(nb, _)| nb)
            .collect();
        if c_neighbors.len() != 2 {
            return Err(IupacError::NotSupported);
        }
        // Build the two chains from each side of the S.
        let chain1_len = count_c_chain(mol, c_neighbors[0], s_idx);
        let chain2_len = count_c_chain(mol, c_neighbors[1], s_idx);
        // Use alphabetical order (IUPAC-preferred).
        let mut names = [alkyl_prefix(chain1_len), alkyl_prefix(chain2_len)];
        names.sort();
        Ok(format!("{} {} sulfide", names[0], names[1]))
    }

    // -----------------------------------------------------------------------
    // Thiol naming (R-SH → "...anethiol")
    // -----------------------------------------------------------------------

    pub(crate) fn name_thiol(
        &self,
        carbons: &[AtomIdx],
        s_idx: AtomIdx,
    ) -> Result<String, IupacError> {
        let mol = self.mol;
        if implicit_hcount(mol, s_idx) == 0 {
            return Err(IupacError::NotSupported);
        }
        let chain = find_longest_c_chain(mol, carbons);
        let n = chain.len();
        let chain_set: HashSet<AtomIdx> = chain.iter().copied().collect();
        let thiol_c = mol
            .neighbors(s_idx)
            .filter(|(nb, _)| chain_set.contains(nb))
            .map(|(nb, _)| nb)
            .next()
            .ok_or(IupacError::NotSupported)?;
        let pos_fwd = chain
            .iter()
            .position(|&c| c == thiol_c)
            .map(|p| p + 1)
            .unwrap_or(1);
        let pos = pos_fwd.min(n + 1 - pos_fwd);
        // Terminal SH (pos=1): no locant; internal: add locant.
        if pos == 1 {
            Ok(format!("{}anethiol", alkane_stem(n)))
        } else {
            Ok(format!("{}ane-{}-thiol", alkane_stem(n), pos))
        }
    }
}
