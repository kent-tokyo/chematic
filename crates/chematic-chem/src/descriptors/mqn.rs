// ---------------------------------------------------------------------------
// MQN: Molecular Quantum Numbers (42 integer descriptors)
// ---------------------------------------------------------------------------

/// Molecular Quantum Numbers, RDKit's `CalcMQNs` (42 counts in the
/// Nguyen et al. 2009 order): atom, polarity, charge, topology, bond and
/// ring-size counts over RDKit's ring information (symmetrized SSSR) and
/// aromatic bonds. The former implementation used a different layout.
pub fn mqn(mol: &Molecule) -> Vec<u32> {
    chematic_smiles::with_rdkit_model_molecule(mol, mqn_impl)
}

fn mqn_impl(mol: &Molecule) -> Vec<u32> {
    let mut res = vec![0u32; 42];
    let view = descriptor_aromaticity(mol);
    let view: &Molecule = &view;
    let rings = rdkit_ring_list(view);
    let rings = &rings[..];
    let (atom_rings, bond_rings) = mqn_ring_memberships(view, rings);
    for (idx, at) in view.atoms() {
        accumulate_mqn_atom(&mut res, view, idx, at, atom_rings[idx.0 as usize]);
    }
    res[11] = view
        .atoms()
        .filter(|(_, a)| !a.wildcard && a.element.atomic_number() > 1)
        .count() as u32;
    let n_aromatic = accumulate_mqn_bonds(&mut res, mol, view, &bond_rings);
    res[15] += n_aromatic / 2 + n_aromatic % 2;
    res[16] += n_aromatic / 2;
    res[18] = rotatable_bond_count(mol) as u32;
    accumulate_mqn_ring_sizes(&mut res, rings);
    res
}

fn mqn_ring_memberships(view: &Molecule, rings: &[Vec<AtomIdx>]) -> (Vec<u32>, Vec<u32>) {
    let mut atom_rings = vec![0u32; view.atom_count()];
    let mut bond_rings = vec![0u32; view.bond_count()];
    for ring in rings {
        for (k, &a) in ring.iter().enumerate() {
            atom_rings[a.0 as usize] += 1;
            let b = ring[(k + 1) % ring.len()];
            if let Some((bi, _)) = view.bond_between(a, b) {
                bond_rings[bi.0 as usize] += 1;
            }
        }
    }
    (atom_rings, bond_rings)
}

fn accumulate_mqn_atom(
    res: &mut [u32],
    view: &Molecule,
    idx: AtomIdx,
    at: &chematic_core::Atom,
    n_rings: u32,
) {
        let z = if at.wildcard {
            0
        } else {
            at.element.atomic_number()
        };
        let n_hs = u32::from(implicit_hcount(view, idx));
        let degree = view.degree(idx);
        match z {
            6 => res[0] += 1,
            9 => res[1] += 1,
            17 => res[2] += 1,
            35 => res[3] += 1,
            53 => res[4] += 1,
            16 => res[5] += 1,
            15 => res[6] += 1,
            7 => {
                if n_rings == 0 {
                    res[7] += 1
                } else {
                    res[8] += 1
                }
                if degree != 4 {
                    res[19] += 1;
                    res[20] += 1;
                }
                if n_hs > 0 {
                    res[21] += n_hs;
                    res[22] += 1;
                }
            }
            8 => {
                if n_rings == 0 {
                    res[9] += 1
                } else {
                    res[10] += 1
                }
                res[20] += 1;
                res[19] += if at.charge != -1 { 2 } else { 3 };
                if n_hs > 0 {
                    res[21] += n_hs;
                    res[22] += 1;
                }
            }
            _ => {}
        }
        if at.charge > 0 {
            res[24] += 1;
        } else if at.charge < 0 {
            res[23] += 1;
        }
        if z != 1 {
            match degree {
                1 => res[25] += 1,
                2 => {
                    if n_rings == 0 {
                        res[26] += 1
                    } else {
                        res[29] += 1
                    }
                }
                3 => {
                    if n_rings == 0 {
                        res[27] += 1
                    } else {
                        res[30] += 1
                    }
                }
                4 => {
                    if n_rings == 0 {
                        res[28] += 1
                    } else {
                        res[31] += 1
                    }
                }
                _ => {}
            }
            if n_rings >= 2 {
                res[40] += 1;
            }
        }
}

fn accumulate_mqn_bonds(
    res: &mut [u32],
    mol: &Molecule,
    view: &Molecule,
    bond_rings: &[u32],
) -> u32 {
    let mut n_aromatic = 0u32;
    for (bi, bond) in view.bonds() {
        let n_rings = bond_rings[bi.0 as usize];
        // An aromatic-flagged ring triple bond (`c1ccccc#1`, triple in
        // `mol`, aromatic in the view) counts as aromatic and as triple.
        let order = if bond.order == BondOrder::Aromatic
            && view.bond_count() == mol.bond_count()
            && mol.bond(bi).order == BondOrder::Triple
        {
            n_aromatic += 1;
            BondOrder::Triple
        } else {
            bond.order
        };
        match order {
            BondOrder::Aromatic => n_aromatic += 1,
            BondOrder::Single | BondOrder::Up | BondOrder::Down => {
                if n_rings == 0 {
                    res[12] += 1
                } else {
                    res[15] += 1
                }
            }
            BondOrder::Double => {
                if n_rings == 0 {
                    res[13] += 1
                } else {
                    res[16] += 1
                }
            }
            BondOrder::Triple => {
                if n_rings == 0 {
                    res[14] += 1
                } else {
                    res[17] += 1
                }
            }
            _ => {}
        }
        if n_rings >= 2 {
            res[41] += 1;
        }
    }
    n_aromatic
}

fn accumulate_mqn_ring_sizes(res: &mut [u32], rings: &[Vec<AtomIdx>]) {
    for ring in rings {
        if ring.len() < 10 {
            res[ring.len() + 29] += 1;
        } else {
            res[39] += 1;
        }
    }
}
