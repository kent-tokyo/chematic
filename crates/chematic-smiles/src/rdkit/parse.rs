//! The molecule RDKit's SMILES parser builds (`toMol`: `CloseMolRings`,
//! `SetUnspecifiedBondTypes`, `AdjustAtomChiralityFlags`), reconstructed
//! from a chematic [`Molecule`] read from SMILES.
//!
//! chematic keeps the parse in its own model; this module converts at the
//! boundary: RDKit's bond numbering ([`Molecule::rdkit_bond_order`]), ring
//! closure begin atoms, `/` `\` directions relative to each bond's begin
//! atom and chiral tags relative to RDKit's bond order.

use chematic_core::{AtomIdx, BondOrder, Chirality, Molecule, STEREO_H_SENTINEL};

use super::RdkitSmilesError;
use super::mol::{Atom, Bond, BondDir, BondType, ChiralTag, Mol, count_swaps};

fn unsupported(what: impl Into<String>) -> RdkitSmilesError {
    RdkitSmilesError::Unsupported(what.into())
}

/// The parser-state molecule (before `removeHs`/sanitization).
pub(crate) fn from_chematic(mol: &Molecule) -> Result<Mol, RdkitSmilesError> {
    let mut out = Mol::default();
    // Stereo stored as CIP labels comes from other formats (MOL files,
    // perception); the SMILES model carries it as tags and directions.
    if mol.atoms().any(|(_, a)| a.cip_code.is_some()) {
        return Err(unsupported(
            "stereo stored as CIP labels (molecule not read from SMILES)",
        ));
    }
    for (_, atom) in mol.atoms() {
        let anum = if atom.wildcard {
            0
        } else {
            u32::from(atom.element.atomic_number())
        };
        let mut a = Atom::new(anum);
        a.isotope = u32::from(atom.isotope.unwrap_or(0));
        a.charge = i32::from(atom.charge);
        a.aromatic = atom.aromatic;
        a.map = atom.atom_map.map(u32::from);
        if !atom.wildcard {
            // Bracket atoms carry their H count and take no implicit Hs.
            if let Some(h) = atom.hydrogen_count {
                a.num_explicit_hs = u32::from(h);
                a.no_implicit = true;
            }
        } else {
            a.num_explicit_hs = u32::from(atom.hydrogen_count.unwrap_or(0));
            a.no_implicit = true;
        }
        out.add_atom(a);
    }

    // RDKit's bond numbering: chain bonds as written, then ring closures.
    let order = mol.rdkit_bond_order();
    let n_closure_ends: usize = (0..mol.atom_count())
        .map(|a| mol.smiles_ring_closure_count(AtomIdx(a as u32)))
        .sum();
    let n_closures = n_closure_ends / 2;
    let first_closure = order.len() - n_closures;
    let mut is_closure = vec![false; mol.bond_count()];
    let mut rd_index = vec![usize::MAX; mol.bond_count()];
    for (k, &bidx) in order.iter().enumerate() {
        rd_index[bidx.0 as usize] = k;
        if k >= first_closure {
            is_closure[bidx.0 as usize] = true;
        }
    }
    for (k, &bidx) in order.iter().enumerate() {
        let bond = mol.bond(bidx);
        let (a1, a2) = (bond.atom1.0 as usize, bond.atom2.0 as usize);
        let (bt, aromatic, dir12) = match bond.order {
            BondOrder::Single => (BondType::Single, false, BondDir::None),
            BondOrder::Up => (BondType::Single, false, BondDir::EndUpRight),
            BondOrder::Down => (BondType::Single, false, BondDir::EndDownRight),
            BondOrder::Double => (BondType::Double, false, BondDir::None),
            BondOrder::Triple => (BondType::Triple, false, BondDir::None),
            BondOrder::Quadruple => (BondType::Quadruple, false, BondDir::None),
            BondOrder::Dative => (BondType::Dative, false, BondDir::None),
            BondOrder::Aromatic => {
                let dir = match mol.bond_direction(bidx) {
                    Some(d) => {
                        let d = if d == BondOrder::Up {
                            BondDir::EndUpRight
                        } else {
                            BondDir::EndDownRight
                        };
                        match mol.bond_direction_anchor(bidx) {
                            Some(anchor) if anchor.0 as usize == a2 => d.flipped(),
                            _ => d,
                        }
                    }
                    None => BondDir::None,
                };
                (BondType::Aromatic, true, dir)
            }
            other => return Err(unsupported(format!("bond order {other:?}"))),
        };
        // `CloseMolRings` keeps the opening partial bond (begin = opening
        // atom) only when it carried an explicit bond order; for an implicit
        // bond or a bare `/` `\` it keeps the closing one (begin = closing
        // atom). chematic stores the opening atom first.
        let swap = is_closure[bidx.0 as usize]
            && match mol.smiles_ring_closure_begins_at_open(bidx) {
                Some(at_open) => !at_open,
                // Not recorded: assume single/aromatic/directional bonds
                // were written bare at the opening digit and double/triple
                // bonds explicitly there.
                None => matches!(
                    bond.order,
                    BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Aromatic
                ),
            };
        let (begin, end, dir) = if swap {
            (a2, a1, dir12.flipped())
        } else {
            (a1, a2, dir12)
        };
        let mut b = Bond::new(begin, end, bt);
        b.aromatic = aromatic;
        b.dir = dir;
        let idx = out.add_bond(b);
        debug_assert_eq!(idx, k);
    }

    // Chiral tags relative to RDKit's bond order (`AdjustAtomChiralityFlags`).
    for (aidx, atom) in mol.atoms() {
        let a = aidx.0 as usize;
        let written = match atom.chirality {
            Chirality::None => continue,
            Chirality::CounterClockwise => ChiralTag::Ccw,
            Chirality::Clockwise => ChiralTag::Cw,
            Chirality::SquarePlanar(_) => {
                return Err(unsupported("non-tetrahedral chirality"));
            }
        };
        let text = mol
            .stereo_neighbor_order(aidx)
            .ok_or_else(|| unsupported("chiral atom without SMILES neighbour order"))?;
        // SMILES bond order (`GetBondOrdering`): the written neighbours
        // without the bracket H.
        let mut smiles_bonds = Vec::with_capacity(text.len());
        for &nb in text {
            if nb == STEREO_H_SENTINEL {
                continue;
            }
            let (bidx, _) = mol
                .bond_between(aidx, AtomIdx(nb))
                .ok_or_else(|| unsupported("stereo neighbour order out of sync"))?;
            smiles_bonds.push(rd_index[bidx.0 as usize]);
        }
        let mut sorted_smiles = smiles_bonds.clone();
        sorted_smiles.sort_unstable();
        let mut sorted_ref = out.atom_bonds[a].clone();
        sorted_ref.sort_unstable();
        if sorted_smiles != sorted_ref {
            return Err(unsupported("stereo neighbour order out of sync"));
        }
        let mut n_swaps = count_swaps(&smiles_bonds, &out.atom_bonds[a]);
        // `chiralAtomNeedsTagInversion` with the parser's state (no property
        // cache yet: only an explicit H is a fourth valence).
        let degree = out.degree(a);
        let is_start = !mol
            .neighbors(aidx)
            .any(|(nb, b)| (nb.0 as usize) < a && !is_closure[b.0 as usize]);
        let n_hs = out.atoms[a].num_explicit_hs;
        let n_closures = mol.smiles_ring_closure_count(aidx);
        let unsaturated = out.atom_bonds[a]
            .iter()
            .any(|&b| out.bonds[b].bt.as_double() > 1.0);
        if degree == 3
            && ((is_start && n_hs == 1) || (n_hs != 1 && n_closures == 1 && !unsaturated))
        {
            n_swaps += 1;
        }
        let mut tag = written;
        if n_swaps % 2 == 1 {
            tag = match tag {
                ChiralTag::Cw => ChiralTag::Ccw,
                _ => ChiralTag::Cw,
            };
        }
        out.atoms[a].chiral = tag;
    }
    Ok(out)
}
