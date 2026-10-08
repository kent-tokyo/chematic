//! `MolOps::removeHs` and `MolOps::sanitizeMol` as `MolFromSmiles` runs them
//! (RDKit 2026.03.1 `AddHs.cpp`, `MolOps.cpp`, `ConjugHybrid.cpp`).

use super::RdkitSmilesError;
use super::aromaticity::{count_atom_elec, set_aromaticity};
use super::kekulize::kekulize;
use super::mol::{BondDir, BondStereo, BondType, ChiralTag, Hybridization, Mol};
use super::periodic;
use super::rank::rank_mol_atoms;

/// `MolOps::removeHs(mol, ps, sanitize=true)` with the parameters
/// `MolFromSmiles` uses (defaults plus `updateExplicitCount=true`).
pub(crate) fn remove_hs_and_sanitize(mol: &mut Mol) -> Result<(), RdkitSmilesError> {
    remove_hs(mol, true).map(|_| ())
}

/// `MolOps::removeHs(mol, ps, sanitize=true)` with the default
/// `RemoveHsParameters` but `updateExplicitCount`; returns the original
/// index of every atom kept.
pub(crate) fn remove_hs(
    mol: &mut Mol,
    update_explicit_count: bool,
) -> Result<Vec<usize>, RdkitSmilesError> {
    for a in 0..mol.atoms.len() {
        mol.update_atom_property_cache(a, false)?;
    }
    let to_remove: Vec<bool> = (0..mol.atoms.len())
        .map(|a| should_remove_h(mol, a))
        .collect();
    for idx in (0..mol.atoms.len()).rev() {
        if to_remove[idx] {
            mol_remove_h(mol, idx, update_explicit_count);
        }
    }
    // `atomsToRemove` is never empty for a non-empty molecule.
    if !mol.atoms.is_empty() {
        sanitize_mol(mol)?;
        for a in 0..mol.atoms.len() {
            if !mol.atoms[a].no_implicit
                && mol.atoms[a].chiral != ChiralTag::Unspecified
                && mol.atoms[a].num_explicit_hs > 1
            {
                mol.atoms[a].num_explicit_hs = 0;
                mol.update_atom_property_cache(a, false)?;
            }
        }
    }
    Ok((0..to_remove.len()).filter(|&i| !to_remove[i]).collect())
}

/// `sanitizeMol` on a molecule whose hydrogen atoms all stay graph atoms
/// (RDKit's state after `Chem.AddHs`).
pub(crate) fn sanitize_keeping_hs(mol: &mut Mol) -> Result<(), RdkitSmilesError> {
    for a in 0..mol.atoms.len() {
        mol.update_atom_property_cache(a, false)?;
    }
    if !mol.atoms.is_empty() {
        sanitize_mol(mol)?;
    }
    Ok(())
}

/// `shouldRemoveH` with the default `RemoveHsParameters`.
fn should_remove_h(mol: &Mol, a: usize) -> bool {
    let atom = &mol.atoms[a];
    if atom.anum != 1 {
        return false;
    }
    let degree = mol.degree(a);
    if degree == 0 || degree > 1 {
        return false;
    }
    if atom.isotope != 0 {
        return false;
    }
    let mut only_h_neighbors = true;
    for &b in &mol.atom_bonds[a] {
        let nbr = mol.bonds[b].other(a);
        if mol.atoms[nbr].anum != 1 {
            only_h_neighbors = false;
        }
        // removeDefiningBondStereo=false
        if mol.degree(nbr) == 2 {
            for &nb in &mol.atom_bonds[nbr] {
                let nbond = &mol.bonds[nb];
                if nbond.bt == BondType::Double
                    && (nbond.stereo > BondStereo::Any || mol.bonds[b].dir != BondDir::None)
                {
                    return false;
                }
            }
        }
    }
    !only_h_neighbors
}

/// `may_need_extra_H`: one single and two aromatic bonds, valence three.
fn may_need_extra_h(mol: &Mol, a: usize) -> bool {
    let mut single = 0;
    let mut aromatic = 0;
    for &b in &mol.atom_bonds[a] {
        match mol.bonds[b].bt {
            BondType::Single => single += 1,
            BondType::Aromatic => aromatic += 1,
            _ => return false,
        }
    }
    single == 1 && aromatic == 2 && mol.total_valence(a) == 3
}

/// `molRemoveH(mol, idx, updateExplicitCount)`.
fn mol_remove_h(mol: &mut Mol, idx: usize, update_explicit_count: bool) {
    let bonds = mol.atom_bonds[idx].clone();
    for b in bonds {
        let heavy = mol.bonds[b].other(idx);
        let atom = &mol.atoms[heavy];
        let bump =
            if update_explicit_count || atom.no_implicit || atom.chiral != ChiralTag::Unspecified {
                true
            } else {
                // Issue 228: an H on an aromatic N or P, or on an atom outside
                // its default valence state, must stay counted.
                let anum = atom.anum;
                ((anum == 7 || anum == 15 || may_need_extra_h(mol, heavy))
                    && mol.is_aromatic_atom(heavy))
                    || periodic::valence_list(anum)
                        .iter()
                        .skip(1)
                        .any(|&v| i32::from(v) == mol.total_valence(heavy))
            };
        if bump {
            mol.atoms[heavy].num_explicit_hs += 1;
        }
        if mol.atoms[heavy].chiral != ChiralTag::Unspecified {
            let mut probe: Vec<usize> = mol.atom_bonds[heavy]
                .iter()
                .copied()
                .filter(|&x| x != b)
                .collect();
            probe.push(b);
            if mol.perturbation_is_odd(heavy, &probe) {
                mol.atoms[heavy].invert_chirality();
            }
        }
        if mol.degree(heavy) == 2 {
            for nb in mol.atom_bonds[heavy].clone() {
                if nb != b {
                    if mol.bonds[nb].stereo > BondStereo::Any {
                        mol.bonds[nb].stereo = BondStereo::None;
                        mol.bonds[nb].stereo_atoms.clear();
                    }
                    break;
                }
            }
        }
        let dir = mol.bonds[b].dir;
        if dir != BondDir::None {
            let mut found_a_dir = false;
            let mut o_bond = None;
            for &nb in &mol.atom_bonds[heavy] {
                if nb != b && mol.bonds[nb].bt == BondType::Single {
                    if mol.bonds[nb].dir == BondDir::None {
                        o_bond = Some(nb);
                    } else {
                        found_a_dir = true;
                    }
                }
            }
            if !found_a_dir && let Some(ob) = o_bond {
                let flip = mol.bonds[ob].begin == heavy && mol.bonds[b].begin == heavy;
                mol.bonds[ob].dir = if flip { dir.flipped() } else { dir };
            }
        }
    }
    mol.remove_atom(idx);
}

/// `MolOps::sanitizeMol` (all operations).
pub(crate) fn sanitize_mol(mol: &mut Mol) -> Result<(), RdkitSmilesError> {
    clean_up(mol)?;
    clean_up_organometallics(mol)?;
    mol.update_property_cache(true)?;
    mol.find_rings()?;
    kekulize(mol)?;
    assign_radicals(mol);
    set_aromaticity(mol);
    set_conjugation(mol);
    set_hybridization(mol);
    cleanup_chirality(mol);
    adjust_hs(mol)?;
    mol.update_property_cache(true)?;
    Ok(())
}

/// `cleanUp`: nitro-like nitrogens, P(=O)=C/N, halogen oxides.
fn clean_up(mol: &mut Mol) -> Result<(), RdkitSmilesError> {
    // nitrogensCleanup
    let mut to_consider = vec![false; mol.atoms.len()];
    for a in 0..mol.atoms.len() {
        if mol.atoms[a].anum != 7 || mol.atoms[a].charge != 0 {
            continue;
        }
        if mol.calc_explicit_valence(a, false)? != 5 {
            continue;
        }
        to_consider[a] = true;
        let arom = mol.atoms[a].aromatic;
        mol.atoms[a].aromatic = false;
        let mut update = false;
        for b in mol.atom_bonds[a].clone() {
            let nbr = mol.bonds[b].other(a);
            if mol.atoms[nbr].anum == 8
                && mol.atoms[nbr].charge == 0
                && mol.bonds[b].bt == BondType::Double
            {
                mol.bonds[b].bt = BondType::Single;
                mol.atoms[a].charge = 1;
                mol.atoms[nbr].charge = -1;
                update = true;
                break;
            }
        }
        mol.atoms[a].aromatic = arom;
        if update {
            mol.calc_explicit_valence(a, false)?;
        }
    }
    for a in 0..mol.atoms.len() {
        if !to_consider[a] {
            continue;
        }
        let arom = mol.atoms[a].aromatic;
        mol.atoms[a].aromatic = false;
        let mut update = false;
        for b in mol.atom_bonds[a].clone() {
            let nbr = mol.bonds[b].other(a);
            if mol.atoms[nbr].anum == 7
                && mol.atoms[nbr].charge == 0
                && mol.bonds[b].bt == BondType::Triple
            {
                mol.bonds[b].bt = BondType::Double;
                mol.atoms[a].charge = 1;
                mol.atoms[nbr].charge = -1;
                update = true;
                break;
            }
        }
        mol.atoms[a].aromatic = arom;
        if update {
            mol.calc_explicit_valence(a, false)?;
        }
    }
    for a in 0..mol.atoms.len() {
        match mol.atoms[a].anum {
            15 => phosphorus_cleanup(mol, a)?,
            17 | 35 | 53 => halogen_cleanup(mol, a)?,
            _ => {}
        }
    }
    Ok(())
}

fn phosphorus_cleanup(mol: &mut Mol, a: usize) -> Result<(), RdkitSmilesError> {
    if mol.atoms[a].charge != 0 {
        return Ok(());
    }
    if mol.calc_explicit_valence(a, false)? == 5 && mol.degree(a) == 3 {
        let mut dbl_to_o = None;
        let mut has_double_to_c_or_n = false;
        for &b in &mol.atom_bonds[a] {
            let nbr = mol.bonds[b].other(a);
            let bt = mol.bonds[b].bt;
            if mol.atoms[nbr].anum == 8 && mol.atoms[nbr].charge == 0 && bt == BondType::Double {
                dbl_to_o = Some((b, nbr));
            } else if (mol.atoms[nbr].anum == 6 || mol.atoms[nbr].anum == 7)
                && mol.degree(nbr) >= 2
                && bt == BondType::Double
            {
                has_double_to_c_or_n = true;
            }
        }
        if has_double_to_c_or_n && let Some((b, o)) = dbl_to_o {
            mol.atoms[o].charge = -1;
            mol.bonds[b].bt = BondType::Single;
            mol.atoms[a].charge = 1;
        }
    }
    mol.calc_explicit_valence(a, false)?;
    Ok(())
}

fn halogen_cleanup(mol: &mut Mol, a: usize) -> Result<(), RdkitSmilesError> {
    let ev = mol.calc_explicit_valence(a, false)?;
    if mol.atoms[a].charge == 0 && (ev == 7 || ev == 5 || ev == 3) {
        let all_o = mol.nbrs(a).all(|nb| mol.atoms[nb].anum == 8);
        if all_o {
            let mut formal = 0;
            for b in mol.atom_bonds[a].clone() {
                if mol.bonds[b].bt == BondType::Double {
                    mol.bonds[b].bt = BondType::Single;
                    let other = mol.bonds[b].other(a);
                    formal += 1;
                    mol.atoms[other].charge = -1;
                    mol.calc_explicit_valence(other, false)?;
                }
            }
            mol.atoms[a].charge = formal;
            mol.calc_explicit_valence(a, false)?;
        }
    }
    Ok(())
}

/// `isHypervalentNonMetal`.
fn is_hypervalent_non_metal(mol: &mut Mol, a: usize) -> Result<bool, RdkitSmilesError> {
    if periodic::is_metal(mol.atoms[a].anum) {
        return Ok(false);
    }
    mol.update_atom_property_cache(a, false)?;
    let ev = mol.atoms[a].explicit_valence;
    let eff = mol.atoms[a].anum as i32 - mol.atoms[a].charge;
    if eff <= 0 {
        return Ok(false);
    }
    let vals = periodic::valence_list(eff as u32);
    let max_v = i32::from(*vals.last().unwrap_or(&-1));
    Ok(max_v > 0
        && (ev > max_v || (ev == max_v && mol.atoms[a].aromatic && mol.total_degree(a) == 4)))
}

fn num_dative_bonds(mol: &Mol, a: usize) -> usize {
    mol.atom_bonds[a]
        .iter()
        .filter(|&&b| mol.bonds[b].bt == BondType::Dative)
        .count()
}

fn no_dative(mol: &Mol, a: usize) -> bool {
    matches!(mol.atoms[a].anum, 1 | 2 | 9 | 10)
}

/// `cleanUpOrganometallics`.
fn clean_up_organometallics(mol: &mut Mol) -> Result<(), RdkitSmilesError> {
    // Only a single bond to a metal needs fixing. Without a metal the scan
    // below would only refresh property caches (non-strictly, which cannot
    // fail), and `sanitize_mol` recomputes every one of them next.
    if !mol.atoms.iter().any(|atom| periodic::is_metal(atom.anum)) {
        return Ok(());
    }
    let mut needs_fixing = false;
    for a in 0..mol.atoms.len() {
        if is_hypervalent_non_metal(mol, a)? && !no_dative(mol, a) {
            for &b in &mol.atom_bonds[a] {
                let other = mol.bonds[b].other(a);
                if mol.bonds[b].bt == BondType::Single && periodic::is_metal(mol.atoms[other].anum)
                {
                    needs_fixing = true;
                    break;
                }
            }
        }
        if needs_fixing {
            break;
        }
    }
    if !needs_fixing {
        return Ok(());
    }
    mol.update_property_cache(false)?;
    // `rankMolAtoms` on a molecule without ring information runs
    // `fastFindRings`; the symmetrized SSSR stands in for it.
    let had_rings = mol.rings.is_some();
    if !had_rings {
        mol.find_rings()?;
    }
    let ranks = rank_mol_atoms(mol);
    if !had_rings {
        mol.rings = None;
    }
    let mut atom_ranks: Vec<(usize, u32)> = ranks.iter().copied().enumerate().collect();
    atom_ranks.sort_by_key(|&(_, r)| r);
    for (a, _) in atom_ranks {
        metal_bond_cleanup(mol, a, &ranks)?;
    }
    Ok(())
}

fn metal_bond_cleanup(mol: &mut Mol, a: usize, ranks: &[u32]) -> Result<(), RdkitSmilesError> {
    if is_hypervalent_non_metal(mol, a)? && !no_dative(mol, a) {
        let mut metals: Vec<usize> = Vec::new();
        for &b in &mol.atom_bonds[a] {
            let other = mol.bonds[b].other(a);
            if mol.bonds[b].bt == BondType::Single && periodic::is_metal(mol.atoms[other].anum) {
                metals.push(other);
            }
        }
        if !metals.is_empty() {
            metals.sort_by(|&x, &y| {
                let (dx, dy) = (num_dative_bonds(mol, x), num_dative_bonds(mol, y));
                if dx == dy {
                    ranks[y].cmp(&ranks[x])
                } else {
                    dx.cmp(&dy)
                }
            });
            if let Some(b) = mol.bond_between(a, metals[0]) {
                let bond = &mut mol.bonds[b];
                bond.bt = BondType::Dative;
                bond.begin = a;
                bond.end = metals[0];
            }
        }
    }
    Ok(())
}

/// `assignRadicals`.
fn assign_radicals(mol: &mut Mol) {
    for a in 0..mol.atoms.len() {
        let atom = &mol.atoms[a];
        if !atom.no_implicit || atom.anum == 0 {
            continue;
        }
        let valens = periodic::valence_list(atom.anum);
        let chg = atom.charge;
        let n_outer = periodic::n_outer_elecs(atom.anum);
        if valens.len() != 1 || valens[0] != -1 {
            let mut accum = 0.0f64;
            for &b in &mol.atom_bonds[a] {
                accum += mol.bonds[b].valence_contrib(a);
            }
            accum += f64::from(atom.num_explicit_hs);
            let total_valence = (accum + 0.1) as i32;
            let base_count = if atom.anum == 1 || atom.anum == 2 {
                2
            } else {
                8
            };
            let mut num_radicals = base_count - n_outer - total_valence + chg;
            if num_radicals < 0 {
                num_radicals = 0;
                if valens.len() > 1 {
                    for &v in valens {
                        let v = i32::from(v);
                        if v - total_valence + chg >= 0 {
                            num_radicals = v - total_valence + chg;
                            break;
                        }
                    }
                }
            }
            let num_radicals2 = n_outer - total_valence - chg;
            if num_radicals2 >= 0 {
                num_radicals = num_radicals.min(num_radicals2);
            }
            mol.atoms[a].radicals = num_radicals.max(0) as u32;
        } else if mol.degree(a) > 0 {
            mol.atoms[a].radicals = 0;
        } else {
            let n_valence = (n_outer - chg).max(0);
            mol.atoms[a].radicals = (n_valence % 2) as u32;
        }
    }
}

/// `isAtomConjugCand`.
fn is_atom_conjug_cand(mol: &Mol, a: usize) -> bool {
    let atom = &mol.atoms[a];
    let vals = periodic::valence_list(atom.anum);
    if atom.charge == 0 && vals[0] >= 0 && mol.total_valence(a) > i32::from(vals[0]) {
        return false;
    }
    let nouter = periodic::n_outer_elecs(atom.anum);
    (atom.anum <= 10 || (nouter != 5 && nouter != 6) || (nouter == 6 && mol.total_degree(a) < 2))
        && count_atom_elec(mol, a) > 0
}

/// `setConjugation`.
fn set_conjugation(mol: &mut Mol) {
    for b in &mut mol.bonds {
        b.conjugated = b.aromatic;
    }
    for a in 0..mol.atoms.len() {
        if !is_atom_conjug_cand(mol, a) {
            continue;
        }
        let sbo = mol.degree(a) + mol.total_num_hs(a) as usize;
        if !(2..=3).contains(&sbo) {
            continue;
        }
        let degree = mol.atom_bonds[a].len();
        for i1 in 0..degree {
            let b1 = mol.atom_bonds[a][i1];
            if mol.bonds[b1].valence_contrib(a) < 1.5
                || !is_atom_conjug_cand(mol, mol.bonds[b1].other(a))
            {
                continue;
            }
            for i2 in 0..degree {
                let b2 = mol.atom_bonds[a][i2];
                if b1 == b2 {
                    continue;
                }
                let at2 = mol.bonds[b2].other(a);
                let sbo2 = mol.degree(at2) + mol.total_num_hs(at2) as usize;
                if sbo2 > 3 {
                    continue;
                }
                if is_atom_conjug_cand(mol, at2) {
                    mol.bonds[b1].conjugated = true;
                    mol.bonds[b2].conjugated = true;
                }
            }
        }
    }
}

/// `MolOps::atomHasConjugatedBond`.
pub(crate) fn atom_has_conjugated_bond(mol: &Mol, a: usize) -> bool {
    mol.atom_bonds[a].iter().any(|&b| mol.bonds[b].conjugated)
}

/// `numBondsPlusLonePairs`.
fn num_bonds_plus_lone_pairs(mol: &Mol, a: usize) -> i32 {
    let mut deg = mol.total_degree(a) as i32;
    for &b in &mol.atom_bonds[a] {
        let bond = &mol.bonds[b];
        if bond.bt == BondType::Dative && a != bond.end {
            deg -= 1;
        }
    }
    let atom = &mol.atoms[a];
    if atom.anum <= 1 {
        return deg;
    }
    let nouter = periodic::n_outer_elecs(atom.anum);
    let total_valence = mol.total_valence(a);
    let chg = atom.charge;
    let num_free = nouter - (total_valence + chg);
    if total_valence + nouter - chg < 8 {
        let num_radicals = atom.radicals as i32;
        let num_lone_pairs = (num_free - num_radicals) / 2;
        deg + num_lone_pairs + num_radicals
    } else {
        deg + num_free / 2
    }
}

/// `setHybridization`.
fn set_hybridization(mol: &mut Mol) {
    for a in 0..mol.atoms.len() {
        if mol.atoms[a].anum == 0 {
            mol.atoms[a].hybrid = Hybridization::Unspecified;
            continue;
        }
        if mol.atoms[a].chiral != ChiralTag::Unspecified && mol.total_degree(a) == 4 {
            mol.atoms[a].hybrid = Hybridization::Sp3;
            continue;
        }
        let norbs = if mol.atoms[a].anum < 89 {
            num_bonds_plus_lone_pairs(mol, a)
        } else {
            mol.total_degree(a) as i32
        };
        mol.atoms[a].hybrid = match norbs {
            0 | 1 => Hybridization::S,
            2 => Hybridization::Sp,
            3 => Hybridization::Sp2,
            4 => {
                if mol.total_degree(a) > 3 || !atom_has_conjugated_bond(mol, a) {
                    Hybridization::Sp3
                } else {
                    Hybridization::Sp2
                }
            }
            5 => Hybridization::Sp3d,
            6 => Hybridization::Sp3d2,
            _ => Hybridization::Unspecified,
        };
    }
}

/// `cleanupChirality`.
fn cleanup_chirality(mol: &mut Mol) {
    for atom in &mut mol.atoms {
        if atom.chiral != ChiralTag::Unspecified && atom.hybrid != Hybridization::Sp3 {
            atom.chiral = ChiralTag::Unspecified;
        }
    }
}

/// `adjustHs`.
fn adjust_hs(mol: &mut Mol) -> Result<(), RdkitSmilesError> {
    for a in 0..mol.atoms.len() {
        let orig_implicit = mol.implicit_valence(a);
        mol.calc_explicit_valence(a, false)?;
        let orig_explicit_hs = mol.atoms[a].num_explicit_hs as i32;
        let new_implicit = mol.calc_implicit_valence(a, false)?;
        if new_implicit < orig_implicit {
            mol.atoms[a].num_explicit_hs = (orig_explicit_hs + orig_implicit - new_implicit) as u32;
            mol.calc_explicit_valence(a, false)?;
        }
    }
    Ok(())
}
